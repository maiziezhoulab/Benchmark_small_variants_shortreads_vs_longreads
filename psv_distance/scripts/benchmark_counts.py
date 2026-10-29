#!/usr/bin/env python3
"""Exact side-specific metric observations from hap.py 0.3.15 GA4GH output.

This implements BlockQuantify::rocEvaluate and addROCValue after qfy has
recomputed BVT and propagated superlocus FILTER annotations. It does NOT
count raw caller VCF rows or infer SNP/INDEL from shared REF/ALT strings.
Official source tag: v0.3.15, commit 84011695b2ff2406c16a335106db6831fb67fdfe.
"""
import argparse
from collections import Counter
import csv
import json


COUNT_COLUMNS = ("TRUTH.TP", "TRUTH.FN", "QUERY.TP", "QUERY.FP", "QUERY.UNK")


def record_regions(record):
    """Return qfy's official trimmed-reference-span region membership."""
    value = record.info.get("Regions", ())
    if isinstance(value, str):
        value = value.split(",")
    return set(value or ())


def passes_filter(record):
    """hap.py regards absent FILTER ('.') and explicit PASS as passing."""
    return not any(x not in ("PASS", "", ".") for x in record.filter.keys())


def observations(record, policy="PASS"):
    """Return each official SNP/INDEL metric observation as a dictionary.

    One shared record can contribute independently to both sides, including
    different BVT classes. Query UNK is included for exclusion accounting but
    never contributes to precision's denominator. Failed-filter truth TP is
    FN under PASS; failed-filter query calls contribute no metric observation.
    Apply this ONLY to qfy-annotated two-sample TRUTH/QUERY benchmark VCFs.
    """
    if policy not in ("PASS", "ALL"):
        raise ValueError("policy must be PASS or ALL")
    if not {"TRUTH", "QUERY"}.issubset(record.samples):
        raise ValueError("Expected benchmark samples TRUTH and QUERY")
    keep_query = policy == "ALL" or passes_filter(record)
    result = []
    for side in ("TRUTH", "QUERY"):
        sample = record.samples[side]
        vtype, original = sample.get("BVT"), sample.get("BD")
        if vtype not in ("SNP", "INDEL"):
            continue
        if side == "TRUTH":
            if original not in ("TP", "FN"):
                continue
            decision = "FN" if original == "TP" and not keep_query else original
        else:
            if original not in ("TP", "FP", "UNK") or not keep_query:
                continue
            decision = original
        result.append({"side": side, "vtype": vtype, "class": vtype,
                       "BD": decision, "decision": decision,
                       "original_decision": original, "count": 1,
                       "policy": policy})
    return result


def load_expected(extended_csv, scopes, policies=("PASS", "ALL")):
    """Read official aggregate rows, including zero expectations for empty cells."""
    expected = Counter()
    seen = set()
    with open(extended_csv, newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["Subset"], row["Filter"], row["Type"])
            if key[0] not in scopes or key[1] not in policies or key[2] not in ("SNP", "INDEL"):
                continue
            if any(row.get(k) != "*" for k in ("Subtype", "Genotype", "QQ")):
                continue
            if key in seen:
                raise ValueError("Duplicate official aggregate row: %r" % (key,))
            seen.add(key)
            for col in COUNT_COLUMNS:
                value = row.get(col, "")
                expected[key + (col,)] = int(float(value)) if value not in ("", None, "nan", "NaN") else 0
    return expected, seen


def compare_counts(observed, extended_csv, scopes, policies=("PASS", "ALL")):
    """Validate a Counter keyed (scope, policy, vtype, SIDE.DECISION)."""
    expected, seen = load_expected(extended_csv, scopes, policies)
    mismatches = []
    for scope in scopes:
        for policy in policies:
            for vtype in ("SNP", "INDEL"):
                for col in COUNT_COLUMNS:
                    key = (scope, policy, vtype, col)
                    if observed[key] != expected[key]:
                        mismatches.append({"key": list(key), "observed": observed[key], "official": expected[key]})
    return {"matched": not mismatches, "mismatches": mismatches,
            "official_rows_present": [list(x) for x in sorted(seen)],
            "unit": "one official BVT/BD observation per benchmark side; VTC is not a numeric weight"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--extended", required=True)
    parser.add_argument("--scope", action="append", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    import pysam
    counts = Counter()
    n_records = 0
    with pysam.VariantFile(args.vcf) as variants:
        for record in variants:
            n_records += 1
            memberships = record_regions(record) & set(args.scope)
            if not memberships:
                continue
            for policy in ("PASS", "ALL"):
                for observation in observations(record, policy):
                    for scope in memberships:
                        counts[(scope, policy, observation["vtype"],
                                observation["side"] + "." + observation["BD"])] += observation["count"]
    report = compare_counts(counts, args.extended, args.scope)
    report["records_scanned"] = n_records
    with open(args.out, "w") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    if not report["matched"]:
        raise SystemExit("Side-specific counts disagree with official hap.py output; see " + args.out)


if __name__ == "__main__":
    main()
