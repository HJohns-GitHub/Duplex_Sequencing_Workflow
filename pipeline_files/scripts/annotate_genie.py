#!/usr/bin/env python3
"""Annotate a VEP-annotated VCF with GENIE_Sample_Count.

Matches on (SYMBOL, HGVSp_Short) extracted from the VEP CSQ field against the
lookup table built by collapse_genie.py. Writes GENIE_Sample_Count into INFO
(0 when not found, mirroring how absent COSMIC counts become 0 downstream).
"""
import argparse
import csv
import sys

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}

def hgvsp_to_short(hgvsp):
    """ENSP...:p.Gly386Val -> p.G386V ; returns None if no protein change."""
    if not hgvsp:
        return None
    p = hgvsp.split(":")[-1]  # drop the ENSP00000... prefix
    if not p.startswith("p."):
        return None
    for three, one in AA3TO1.items():
        p = p.replace(three, one)
    return p

def load_lookup(path):
    lut = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            lut[(r["SYMBOL"], r["HGVSp_Short"])] = r["GENIE_Sample_Count"]
    return lut

def get_csq_format(header_lines):
    """Read the CSQ field order from the ##INFO=<ID=CSQ ...> header line."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            fmt = line.split("Format:")[-1].rstrip(">\"\n ")
            return fmt.split("|")
    sys.exit("No CSQ header found in VCF")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("lookup")
    ap.add_argument("out")
    args = ap.parse_args()

    lut = load_lookup(args.lookup)

    with open(args.vcf) as f:
        lines = f.readlines()

    header = [l for l in lines if l.startswith("#")]
    csq_cols = get_csq_format(header)
    sym_i = csq_cols.index("SYMBOL")
    hgvsp_i = csq_cols.index("HGVSp")

    new_info_header = (
        '##INFO=<ID=GENIE_Sample_Count,Number=1,Type=Integer,'
        'Description="Unique GENIE patients with this gene + protein change">\n'
    )

    out = []
    inserted = False
    for line in lines:
        if line.startswith("##"):
            out.append(line)
            continue
        if line.startswith("#CHROM"):
            if not inserted:
                out.append(new_info_header)
                inserted = True
            out.append(line)
            continue

        cols = line.rstrip("\n").split("\t")
        info = cols[7]
        count = 0
        for field in info.split(";"):
            if field.startswith("CSQ="):
                # take the first transcript annotation
                first = field[4:].split(",")[0].split("|")
                symbol = first[sym_i] if len(first) > sym_i else ""
                short = hgvsp_to_short(first[hgvsp_i] if len(first) > hgvsp_i else "")
                if symbol and short:
                    count = lut.get((symbol, short), 0)
                break
        cols[7] = info + f";GENIE_Sample_Count={count}"
        out.append("\t".join(cols) + "\n")

    with open(args.out, "w") as f:
        f.writelines(out)

if __name__ == "__main__":
    main()
