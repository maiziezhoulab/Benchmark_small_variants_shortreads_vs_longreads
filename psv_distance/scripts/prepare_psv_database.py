#!/usr/bin/env python3
"""Prepare a caller-independent DuploMap v0.9 PSV index (run full input in Slurm).

All internal coordinates are zero-based half-open. Graph labels are one-based
inclusive; VCF positions and INFO/paralog positions are one-based. The official
VCF already contains both copies. It must not be expanded a second time.
"""
from __future__ import annotations

import argparse
import bisect
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import struct
import sys
import time


SOURCE_URLS = {
    "database": "https://dl.dropboxusercontent.com/s/bicrz0i7vw2r0e5/hg38.tar.gz",
    "psv_vcf": "https://dl.dropboxusercontent.com/s/fdz38axuyj0bxff/hg38.psvs.vcf.gz",
    "repository": "https://gitlab.com/tprodanov/duplomap",
    "v0.9.0_commit": "d25092450e9f3e3f7092239a3e0921cece303a56",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def effective_span(start0, ref, alt):
    """Changed reference bases, or one left-anchor base for pure insertion.

    Common prefix and suffix are removed without shifting the original aligned
    PSV representation. An insertion at an otherwise empty changed reference
    span uses its immediately preceding reference base. Raw coordinates/alleles
    remain in the index. Del/complex spans therefore detect PSV overlap.
    """
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    r, a = ref[prefix:], alt[prefix:]
    suffix = 0
    while suffix < min(len(r), len(a)) and r[-1-suffix] == a[-1-suffix]:
        suffix += 1
    if suffix:
        r, a = r[:-suffix], a[:-suffix]
    pos = start0 + prefix
    if not r:
        return max(0, pos - 1), max(0, pos - 1) + 1, "INS_ANCHOR"
    kind = "SNP" if len(r) == len(a) == 1 else (
        "MNV" if len(r) == len(a) else "DEL" if not a else "COMPLEX")
    return pos, pos + len(r), kind


def merge_intervals(intervals, abutting=False):
    """Union intervals; copy envelopes deliberately do NOT join abutting nodes."""
    out = []
    for start, end in sorted(set(intervals)):
        if out and (start < out[-1][1] or (abutting and start == out[-1][1])):
            out[-1] = (out[-1][0], max(out[-1][1], end))
        else:
            out.append((start, end))
    return out


def parse_graph(path, family_names):
    """Use original aligned node extents, not fragmented post-divergence BED rows.

    Overlapping nodes in a family/chromosome form a local copy envelope. A
    genomic position cannot be two distinct physical reference copies. The
    coarse family remains a database component, not a biological gene family.
    Edge labels are unused: DuploMap v0.9 writes P for proximity and O for
    paralog alignments, contrary to a tempting interpretation of the letters.
    """
    rx = re.compile(r'^\s*(\d+) \[label="([^:]+):(\d+)-(\d+)\. Component (\d+)"\]')
    grouped = defaultdict(list)
    nodes = []
    with open(path) as fh:
        for line in fh:
            m = rx.match(line)
            if not m:
                continue
            node, chrom, start, end, component = m.groups()
            family = family_names[int(component)]
            start, end = int(start) - 1, int(end)
            nodes.append((int(node), family, chrom, start, end))
            grouped[family, chrom].append((start, end))
    if not nodes:
        raise ValueError("No graph nodes parsed; database schema may have changed")
    copies = []
    for (family, chrom), spans in sorted(grouped.items()):
        for j, (start, end) in enumerate(merge_intervals(spans), 1):
            copies.append((f"{family}|{chrom}|{start}-{end}", family, chrom, start, end))
    return copies, nodes


class IntervalLookup:
    """Small family/chromosome interval lookup supporting overlapping families."""
    def __init__(self, rows):
        self.groups = defaultdict(list)
        for copy_id, family, chrom, start, end in rows:
            self.groups[family, chrom].append((start, end, copy_id))
        self.starts = {}
        for key, group in self.groups.items():
            group.sort()
            self.starts[key] = [x[0] for x in group]

    def containing(self, family, chrom, start, end):
        key = family, chrom
        group = self.groups.get(key, [])
        i = bisect.bisect_right(self.starts.get(key, []), start) - 1
        if i >= 0 and start >= group[i][0] and end <= group[i][1]:
            return group[i][2]
        return None


def read_binary_database(path):
    """Yield typed sections from DuploMap's gzip + big-endian binary format.

    Used only for paired homology windows; preserves reference coordinates and
    orientation and checks serialization end-to-end. No read data are required.
    """
    with gzip.open(path, "rb") as fh:
        def unpack(fmt):
            n = struct.calcsize(fmt)
            data = fh.read(n)
            if len(data) != n:
                raise ValueError(f"Truncated database {path}")
            return struct.unpack(fmt, data)

        nchrom, = unpack(">H")
        chroms = {}
        for _ in range(nchrom):
            cid, namelen = unpack(">HB")
            chroms[cid] = fh.read(namelen).decode()

        def span(directional=False):
            cid, start, end = unpack(">hII" if directional else ">HII")
            strand = "+"
            if directional and cid < 0:
                cid, strand = -cid - 1, "-"
            if end < start:
                raise ValueError(f"Invalid span in {path}")
            return chroms[cid], start, end, strand

        nregions, = unpack(">I")
        for _ in range(nregions):
            yield "region", (span(),)
        nwindows, = unpack(">I")
        for _ in range(nwindows):
            yield "window", (span(True), span(True))
        npsv, = unpack(">I")
        for _ in range(npsv):
            yield "psv", (span(True), span(True))
        data = fh.read(4)
        if data:
            nhigh, = struct.unpack(">I", data)
            for _ in range(nhigh):
                yield "high_cn", (span(),)
        if fh.read(1):
            raise ValueError(f"Unexpected trailing data in {path}")


SCHEMA = """
CREATE TABLE copies(copy_id TEXT PRIMARY KEY, family TEXT, chrom TEXT, start0 INTEGER, end0 INTEGER);
CREATE INDEX copies_chrom ON copies(chrom,start0,end0);
CREATE TABLE coverage(chrom TEXT,start0 INTEGER,end0 INTEGER,family TEXT,copy_id TEXT);
CREATE INDEX coverage_chrom ON coverage(chrom,start0,end0);
CREATE TABLE raw_coverage(chrom TEXT,start0 INTEGER,end0 INTEGER,family TEXT,copy_id TEXT);
CREATE TABLE psv_loci(
 locus_id INTEGER PRIMARY KEY, copy_id TEXT, family TEXT, chrom TEXT,
 start0 INTEGER,end0 INTEGER,variant_class TEXT,psv_hash TEXT,
 raw_start0 INTEGER,ref TEXT,alt TEXT,partner_chrom TEXT,partner_start0 INTEGER,
 partner_strand TEXT,partner_copy_id TEXT,complexity REAL,
 usable INTEGER,status TEXT);
CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
"""


def prepare(args):
    import pysam
    began = time.time()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    dbpath = out / "psv_index.sqlite"
    if dbpath.exists():
        raise FileExistsError(f"Preserving existing output: {dbpath}")
    fasta = pysam.FastaFile(args.reference)
    ref_lengths = dict(zip(fasta.references, fasta.lengths))
    root = Path(args.database)
    family_names = {}
    with open(root / "summary.csv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            family_names[int(row["component"].split("-")[0])] = row["component"]
    all_copies, nodes = parse_graph(root / "duplications.dot", family_names)
    copies = [r for r in all_copies if r[2] in ref_lengths and r[4] <= ref_lengths[r[2]]]
    lookup = IntervalLookup(copies)
    qc = Counter(graph_nodes=len(nodes), graph_copy_envelopes=len(all_copies),
                 reference_compatible_copy_envelopes=len(copies))
    con = sqlite3.connect(dbpath)
    con.executescript(SCHEMA)
    con.executemany("INSERT INTO copies VALUES(?,?,?,?,?)", copies)
    # Keep exact raw coverage; eligible coverage below additionally requires a
    # counterpart in the actual reference and an unambiguous local copy.
    raw_cov = []
    with open(root / "regions.bed") as fh:
        for line in fh:
            chrom, start, end, family = line.rstrip().split("\t")[:4]
            start, end = int(start), int(end)
            raw_cov.append((chrom, start, end, family, lookup.containing(family, chrom, start, end)))
    con.executemany("INSERT INTO raw_coverage VALUES(?,?,?,?,?)", raw_cov)
    coverage_parts = defaultdict(list)
    with gzip.open(out / "homology_windows.tsv.gz", "wt") as winout:
        writer = csv.writer(winout, delimiter="\t", lineterminator="\n")
        writer.writerow(["family", "chrom1", "start1_0", "end1_0", "strand1", "copy1",
                         "chrom2", "start2_0", "end2_0", "strand2", "copy2", "status"])
        for dbfile in sorted(root.glob("*.db")):
            family = dbfile.stem
            for kind, spans in read_binary_database(dbfile):
                qc[f"binary_{kind}_records"] += 1
                if kind != "window":
                    continue
                a, b = spans
                ids = [lookup.containing(family, s[0], s[1], s[2]) for s in spans]
                available = all(s[0] in ref_lengths and s[2] <= ref_lengths[s[0]] for s in spans)
                status = "usable" if available and all(ids) and ids[0] != ids[1] else (
                    "counterpart_not_in_reference" if not available else
                    "same_copy_envelope" if all(ids) and ids[0] == ids[1] else "ambiguous_copy")
                qc[f"homology_window_{status}"] += 1
                writer.writerow([family, *a, ids[0] or "", *b, ids[1] or "", status])
                if status == "usable":
                    for s, copy_id in zip(spans, ids):
                        coverage_parts[family, s[0], copy_id].append((s[1], s[2]))
    coverage = []
    for (family, chrom, copy_id), spans in coverage_parts.items():
        for start, end in merge_intervals(spans, abutting=True):
            coverage.append((chrom, start, end, family, copy_id))
    con.executemany("INSERT INTO coverage VALUES(?,?,?,?,?)", coverage)
    con.commit()

    header_lengths = {}
    batch = []
    mismatch_examples = []
    with gzip.open(args.psv_vcf, "rt") as fh, open(out / "psv_original_header.txt", "w") as header:
        for line in fh:
            if line.startswith("#"):
                header.write(line)
                m = re.match(r"##contig=<ID=([^,]+),length=(\d+)>", line)
                if m:
                    header_lengths[m[1]] = int(m[2])
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 8:
                raise ValueError("Malformed PSV VCF record")
            chrom, pos, _, ref, alt, _, _, info = fields[:8]
            qc["vcf_records"] += 1
            if len(fields) > 8:
                qc["vcf_trailing_field_records"] += 1
            inf = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
            family, psv_hash = inf["db"], inf["psv"]
            partner_chrom, partner_pos, strand = inf["paralog"].rsplit(":", 2)
            raw_start, partner_start = int(pos)-1, int(partner_pos)-1
            start, end, variant_class = effective_span(raw_start, ref, alt)
            copy_id = lookup.containing(family, chrom, raw_start, raw_start + len(ref))
            partner_copy = lookup.containing(family, partner_chrom, partner_start, partner_start + len(alt))
            status = "usable"
            if not set(ref + alt) <= set("ACGT"):
                status = "ambiguous_allele"
            elif chrom not in ref_lengths or partner_chrom not in ref_lengths:
                status = "counterpart_not_in_reference"
            elif header_lengths.get(chrom) != ref_lengths[chrom] or header_lengths.get(partner_chrom) != ref_lengths[partner_chrom]:
                status = "contig_length_mismatch"
            elif copy_id is None or partner_copy is None:
                status = "ambiguous_copy"
            elif copy_id == partner_copy:
                status = "same_copy_envelope"
            else:
                observed_ref = fasta.fetch(chrom, raw_start, raw_start + len(ref)).upper()
                observed_partner = fasta.fetch(partner_chrom, partner_start, partner_start + len(alt)).upper()
                expected_partner = alt if strand == "+" else alt.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                if observed_ref != ref or observed_partner != expected_partner:
                    status = "reference_allele_mismatch"
                    if len(mismatch_examples) < 100:
                        mismatch_examples.append(dict(chrom=chrom,pos=int(pos),ref=ref,observed=observed_ref,
                            partner=inf["paralog"],expected_partner=expected_partner,observed_partner=observed_partner))
                else:
                    qc["reference_validated_vcf_records"] += 1
            qc[f"psv_initial_{status}"] += 1
            batch.append((copy_id, family, chrom, start, end, variant_class, psv_hash,
                          raw_start, ref, alt, partner_chrom, partner_start, strand, partner_copy,
                          float(inf.get("complexity", "nan")), int(status == "usable"), status))
            if len(batch) >= 10000:
                con.executemany("INSERT INTO psv_loci VALUES(NULL," + ",".join("?"*17) + ")", batch)
                con.commit()
                batch.clear()
    if batch:
        con.executemany("INSERT INTO psv_loci VALUES(NULL," + ",".join("?"*17) + ")", batch)
        con.commit()
    con.executescript("""
    CREATE INDEX psv_partner ON psv_loci(family,psv_hash,chrom,raw_start0,partner_chrom,partner_start0);
    CREATE INDEX psv_copy_position ON psv_loci(copy_id,start0,end0);
    CREATE INDEX psv_status ON psv_loci(usable,status);
    """)
    # Pair validation is deliberately performed after all copies have been read.
    con.execute("""UPDATE psv_loci AS a SET usable=0,status='missing_reciprocal_psv'
    WHERE usable=1 AND NOT EXISTS(SELECT 1 FROM psv_loci b WHERE b.family=a.family
    AND b.psv_hash=a.psv_hash AND b.chrom=a.partner_chrom AND b.raw_start0=a.partner_start0
    AND b.partner_chrom=a.chrom AND b.partner_start0=a.raw_start0
    AND b.partner_strand=a.partner_strand AND b.status='usable')""")
    # Keep all edges in the table. Mark exact repeated edge rows so nearest
    # searches and counts can use usable=1 without multiplying the same edge.
    con.execute("""UPDATE psv_loci SET usable=0,status='duplicate_edge_record' WHERE usable=1
    AND locus_id NOT IN (SELECT MIN(locus_id) FROM psv_loci WHERE usable=1
       GROUP BY copy_id,family,chrom,raw_start0,ref,alt,partner_chrom,partner_start0,partner_strand,psv_hash)""")
    con.commit()
    for status, n in con.execute("SELECT status,COUNT(*) FROM psv_loci GROUP BY status"):
        qc[f"psv_final_{status}"] = n
    qc["usable_distinct_reference_spans"] = con.execute(
        "SELECT COUNT(*) FROM (SELECT copy_id,start0,end0 FROM psv_loci WHERE usable=1 GROUP BY copy_id,start0,end0)").fetchone()[0]
    qc["usable_pair_hashes"] = con.execute(
        "SELECT COUNT(*) FROM (SELECT family,psv_hash FROM psv_loci WHERE usable=1 GROUP BY family,psv_hash)").fetchone()[0]
    qc["eligible_coverage_bp_sum"] = sum(r[2] - r[1] for r in coverage)
    qc["raw_coverage_bp_sum"] = sum(r[2] - r[1] for r in raw_cov)
    for name, query in [
            ("copies.tsv", "SELECT * FROM copies ORDER BY chrom,start0"),
            ("coverage.tsv", "SELECT * FROM coverage ORDER BY chrom,start0"),
            ("raw_coverage.tsv", "SELECT * FROM raw_coverage ORDER BY chrom,start0")]:
        cursor = con.execute(query)
        with open(out / name, "w") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow([c[0] for c in cursor.description])
            writer.writerows(cursor)
    with open(out / "coverage.bed", "w") as fh:
        for chrom, start, end, family, copy_id in sorted(coverage):
            fh.write(f"{chrom}\t{start}\t{end}\t{family}\t{copy_id}\n")
    with gzip.open(out / "psv_loci.tsv.gz", "wt") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        cursor = con.execute("SELECT * FROM psv_loci ORDER BY locus_id")
        writer.writerow([c[0] for c in cursor.description])
        writer.writerows(cursor)
    metadata = {
        "sources": SOURCE_URLS, "database_label": "DuploMap hg38 v0.9 All",
        "command": sys.argv, "reference": str(Path(args.reference).resolve()),
        "reference_fai_sha256": sha256(args.reference + ".fai"),
        "psv_vcf_sha256": sha256(args.psv_vcf),
        "database_metadata_sha256": {p.name: sha256(p) for p in
             [root / "regions.bed", root / "duplications.dot", root / "summary.csv", root / "description.txt"]},
        "database_archive_sha256": sha256(args.archive) if args.archive else None,
        "construction_command": (root / "description.txt").read_text().splitlines()[0],
        "coordinate_system": "All index coordinates zero-based half-open; source VCF/paralog 1-based; graph 1-based-inclusive",
        "copy_definition": "Within each database component/chromosome, union strictly overlapping original aligned graph nodes; do not join touching or separated nodes",
        "coverage_definition": "Union of actual binary paired homology windows with both copies present in reference and mapped to distinct unambiguous local copy envelopes",
        "distance_span_definition": "Trim shared REF/ALT prefix and suffix; changed reference bases as half-open interval; pure insertion uses immediately preceding reference base (left anchor)",
        "reliability": "Reference-derived only; QUAL100/PASS do not establish sample-specific PSV reliability",
        "graph_edge_warning": "Do not infer relation from dot P/O labels; v0.9 source emits P for proximity and O for paralog relation",
        "special_statuses": [x[0] for x in con.execute("SELECT DISTINCT status FROM psv_loci")],
        "missing_reference_contigs": sorted(set(header_lengths) - set(ref_lengths)),
        "contig_length_mismatches": {c: [n,ref_lengths[c]] for c,n in header_lengths.items() if c in ref_lengths and n != ref_lengths[c]},
        "mismatch_examples": mismatch_examples, "qc": dict(qc),
        "software": {"python": sys.version, "pysam": pysam.__version__, "sqlite": sqlite3.sqlite_version},
        "elapsed_seconds": time.time() - began,
    }
    con.execute("INSERT INTO metadata VALUES(?,?)", ("provenance", json.dumps(metadata)))
    con.commit()
    con.close()
    fasta.close()
    with open(out / "database_qc.json", "w") as fh:
        json.dump(metadata, fh, indent=2)
    print(json.dumps(metadata["qc"], indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True, help="Extracted hg38 directory with cluster .db files")
    parser.add_argument("--psv-vcf", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--archive", help="Archive path to record SHA256")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    prepare(args)


if __name__ == "__main__":
    main()
