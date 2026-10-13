#!/usr/bin/env python3
"""Fetch the GRCh38 host reference into data/reference/ (checked against a pinned SHA-256) and
the PRJEB90335 reads into data/fastq/<batch>/ (checked against data/samplesheet.tsv). Files
already present with the right checksum are kept, so reruns only fetch what is missing.

Uses ENA's submitted files, not fastq_ftp: ENA regenerates the latter, so only the
submitted files are byte-identical to the input of the published results.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import urllib.parse
import urllib.request
from pathlib import Path

FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"

# Host reference used for read filtering: Ensembl release 112, GRCh38 soft-masked primary
# assembly. The SHA-256 was taken from the file used for the published results; Ensembl's own
# CHECKSUMS for this file read "19039 920395" (BSD sum).
HOST_FASTA_URL = (
    "https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/"
    "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
)
HOST_FASTA_SHA256 = "651561b3065aa6083b62290ec96344a90d07388883874070aa368b745fde68fd"
HOST_FASTA = Path("data/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz")


def read_tsv(text: str) -> list[dict[str, str]]:
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def checksum(path: Path, algorithm: str = "md5") -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            h.update(block)
    return h.hexdigest()


def md5sum(path: Path) -> str:
    return checksum(path, "md5")


def submitted_urls(study: str) -> dict[str, str]:
    """File name -> HTTPS URL of every submitted file in the study."""
    query = urllib.parse.urlencode(
        {"accession": study, "result": "read_run", "fields": "submitted_ftp", "format": "tsv"}
    )
    with urllib.request.urlopen(f"{FILEREPORT}?{query}", timeout=120) as resp:
        runs = read_tsv(resp.read().decode())
    urls = [u for r in runs for u in r["submitted_ftp"].split(";") if u]
    return {Path(u).name: f"https://{u}" for u in urls}


def download(url: str, dest: Path, expected: str, algorithm: str = "md5") -> None:
    if dest.exists() and checksum(dest, algorithm) == expected:
        return
    print(f"downloading {dest.name}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    with urllib.request.urlopen(url, timeout=600) as resp, tmp.open("wb") as out:
        while block := resp.read(1 << 20):
            out.write(block)
    if checksum(tmp, algorithm) != expected:
        tmp.unlink()
        sys.exit(f"checksum mismatch for {dest.name}")
    tmp.rename(dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", default="PRJEB90335")
    parser.add_argument("--samplesheet", default=Path("data/samplesheet.tsv"), type=Path)
    parser.add_argument("--outdir", default=Path("data/fastq"), type=Path)
    args = parser.parse_args()

    download(HOST_FASTA_URL, HOST_FASTA, HOST_FASTA_SHA256, "sha256")
    print(f"host reference ready: {HOST_FASTA}")

    sheet = read_tsv(args.samplesheet.read_text())
    files = [
        (args.outdir / row["batch"] / name, md5)
        for row in sheet
        for name, md5 in ((row["forward_file"], row["forward_md5"]), (row["reverse_file"], row["reverse_md5"]))
    ]
    missing = [(dest, md5) for dest, md5 in files if not (dest.exists() and md5sum(dest) == md5)]
    if missing:
        urls = submitted_urls(args.study)
        if not urls:
            sys.exit(f"ENA lists no submitted files for {args.study} (not public yet?)")
        for dest, md5 in missing:
            if dest.name not in urls:
                sys.exit(f"{dest.name}: not among the submitted files of {args.study}")
            download(urls[dest.name], dest, md5)
    print(f"{len(sheet)} samples ready in {args.outdir}")


if __name__ == "__main__":
    main()
