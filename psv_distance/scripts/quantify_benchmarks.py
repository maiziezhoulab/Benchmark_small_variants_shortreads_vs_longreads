#!/usr/bin/env python3
"""Re-quantify a full-context existing xcmp benchmark on a Slurm compute node.

No haplotype rematching, VCF prefiltering, or distance-bin subsetting is done.
PASS and ALL official tables are both retained. Literal GIAB confident BED
boundaries are used; historical automatic insertion padding is not reapplied.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess


ROOT = Path(__file__).resolve().parents[1]
HAPPY_BIN = Path("/path/to/hap.py/bin")


def metadata(path, checksum=False):
    path = Path(path)
    stat = path.stat()
    result = {"path": str(path.resolve()), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if checksum:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        result["sha256"] = digest.hexdigest()
    return result


def corrected_header(benchmark_vcf, reference):
    """Correct stale assembly metadata without altering any VCF record."""
    with open(str(reference) + ".fai") as handle:
        lengths = {row[0]: int(row[1]) for row in (line.split("\t") for line in handle)}
    result, changes = [], []
    saw_reference = False
    with gzip.open(benchmark_vcf, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            if line.startswith("##contig=<"):
                match = re.search(r"(?:<|,)ID=([^,>]+)", line)
                if not match or match.group(1) not in lengths:
                    raise ValueError("Benchmark contig absent from actual reference FAI: " + line.rstrip())
                chrom = match.group(1)
                replacement = "##contig=<ID=%s,length=%d>\n" % (chrom, lengths[chrom])
                if replacement != line:
                    changes.append({"old": line.rstrip(), "new": replacement.rstrip()})
                line = replacement
            elif line.startswith("##reference="):
                saw_reference = True
                replacement = "##reference=" + str(reference) + "\n"
                if replacement != line:
                    changes.append({"old": line.rstrip(), "new": replacement.rstrip()})
                line = replacement
            elif line.startswith("#CHROM"):
                if line.rstrip().split("\t")[-2:] != ["TRUTH", "QUERY"]:
                    raise ValueError("Expected ordered benchmark samples TRUTH and QUERY")
                if not saw_reference:
                    result.append("##reference=" + str(reference) + "\n")
                result.append("##PSVDistanceHeaderCorrection=GRCh38_reference_FAI_metadata_only_no_variant_or_matching_changes\n")
            result.append(line)
    return "".join(result), changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/inputs.json")
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument("--index", type=int, help="Zero-based callset index; suitable for a Slurm array")
    selector.add_argument("--callset", help="Exact callset id in config")
    parser.add_argument("--stratification", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "work/quantified")
    parser.add_argument("--happy-bin", type=Path, default=HAPPY_BIN)
    parser.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    parser.add_argument("--dry-run", action="store_true", help="Validate inputs and print command without execution")
    args = parser.parse_args()
    with args.config.open() as handle:
        config = json.load(handle)
    if args.index is not None:
        if not 0 <= args.index < len(config["callsets"]):
            parser.error("--index outside configured callset list")
        callset = config["callsets"][args.index]
    else:
        candidates = [x for x in config["callsets"] if x["id"] == args.callset]
        if len(candidates) != 1:
            parser.error("--callset must match exactly one configured id")
        callset = candidates[0]
    for path in (callset["benchmark_vcf"], config["reference"], config["highconf_bed"],
                 args.stratification, args.happy_bin / "python", args.happy_bin / "qfy.py",
                 args.happy_bin / "quantify", args.happy_bin / "bcftools", config["reference"] + ".fai"):
        if not Path(path).is_file():
            parser.error("Missing input/executable: " + str(path))
    scopes = []
    with args.stratification.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            label, path = line.rstrip("\n").split("\t", 1)
            if label.startswith("CONF") or label in scopes:
                parser.error("Stratum labels must be unique and must not start CONF")
            bed = Path(path)
            if not bed.is_file():
                bed = args.stratification.resolve().parent / bed
            if not bed.is_file():
                parser.error("Missing stratum BED: " + path)
            scopes.append(label)
    if set(scopes) != {"Segdup", "Segdup_lowmap"}:
        parser.error("Expected exactly the two scopes Segdup and Segdup_lowmap")
    prefix = args.output_dir.resolve() / callset["id"]
    corrected_vcf = str(prefix) + ".input_context.vcf.gz"
    header_path = str(prefix) + ".corrected_header.txt"
    reheader_command = [str(args.happy_bin / "bcftools"), "reheader", "-h", header_path,
                        "-o", corrected_vcf, callset["benchmark_vcf"]]
    index_command = [str(args.happy_bin / "bcftools"), "index", "--tbi", corrected_vcf]
    command = [str(args.happy_bin / "python"), "-B", str(args.happy_bin / "qfy.py"),
               corrected_vcf, "-t", "ga4gh", "-r", config["reference"],
               "-f", config["highconf_bed"], "--stratification", str(args.stratification.resolve()),
               "-V", "--output-vtc", "--preserve-info", "--no-roc", "--threads", str(args.threads),
               "-o", str(prefix)]
    for executable_command in (reheader_command, index_command, command):
        print(shlex.join(executable_command), flush=True)
    if args.dry_run:
        return
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Full VCF quantification must run in a Slurm allocation; use --dry-run on login nodes")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if list(args.output_dir.glob(callset["id"] + ".*")):
        parser.error("Output prefix already exists; choose a new directory to preserve previous outputs")
    env = os.environ.copy()
    env.update({"PATH": str(args.happy_bin) + os.pathsep + env.get("PATH", ""),
                "PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1",
                "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "HGREF": config["reference"]})
    env.pop("PYTHONPATH", None)
    header_text, header_changes = corrected_header(callset["benchmark_vcf"], config["reference"])
    Path(header_path).write_text(header_text)
    provenance = {"callset": callset["id"], "command": command,
                  "preparation_commands": [reheader_command, index_command],
                  "header_correction": {"basis": "actual GRCh38 reference FAI and executed benchmark runinfo",
                                        "changes": header_changes, "variant_records_altered": False},
                  "start_utc": datetime.now(timezone.utc).isoformat(),
                  "slurm_job_id": env.get("SLURM_JOB_ID"), "slurm_array_task_id": env.get("SLURM_ARRAY_TASK_ID"),
                  "host": os.uname().nodename, "matching": "existing full-context xcmp; no rematching",
                  "primary_filter": "official postmatch PASS quantification", "sensitivity_filter": "ALL",
                  "confidence_policy": "literal highconf BED; no --adjust-conf-regions",
                  "inputs": [metadata(args.config, True), metadata(args.stratification, True),
                             metadata(callset["benchmark_vcf"]), metadata(config["reference"]),
                             metadata(config["highconf_bed"], True)],
                  "software": [metadata(args.happy_bin / "qfy.py", True), metadata(args.happy_bin / "quantify", True)],
                  "status": "running"}
    provenance_path = Path(str(prefix) + ".execution.json")
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    with open(str(prefix) + ".stdout.log", "w") as out, open(str(prefix) + ".stderr.log", "w") as err:
        for executable_command in (reheader_command, index_command, command):
            run = subprocess.run(executable_command, env=env, stdout=out, stderr=err)
            if run.returncode:
                break
    provenance.update({"end_utc": datetime.now(timezone.utc).isoformat(), "returncode": run.returncode,
                       "status": "complete" if run.returncode == 0 else "failed"})
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    if run.returncode:
        raise SystemExit(run.returncode)
    for suffix in (".vcf.gz", ".vcf.gz.tbi", ".extended.csv", ".summary.csv"):
        if not Path(str(prefix) + suffix).is_file():
            raise SystemExit("qfy reported success but required output is absent: " + str(prefix) + suffix)


if __name__ == "__main__":
    main()
