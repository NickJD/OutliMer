"""Profile management commands for OutliMer."""

from __future__ import annotations

import argparse
import csv
import json
import os
from typing import Iterable

from OutliMer import __version__
from OutliMer.OutliMer import load_hash_db, save_hash_db


def _profile_paths(profile_dir: str) -> dict[str, str]:
    return {
        "db_gz": os.path.join(profile_dir, "db.txt.gz"),
        "db": os.path.join(profile_dir, "db.txt"),
        "agg": os.path.join(profile_dir, "agg_counts.csv"),
        "samples": os.path.join(profile_dir, "samples.txt"),
    }


def load_profile(profile_dir: str) -> dict:
    paths = _profile_paths(profile_dir)
    db_path = paths["db_gz"] if os.path.exists(paths["db_gz"]) else paths["db"]
    hashes = load_hash_db(db_path) if os.path.exists(db_path) else set()
    agg_counts: dict[int, int] = {}
    if os.path.exists(paths["agg"]):
        with open(paths["agg"], newline="") as fh:
            for row in csv.reader(fh):
                if len(row) < 2:
                    continue
                try:
                    agg_counts[int(row[0])] = int(row[1])
                except ValueError:
                    continue
    samples: list[str] = []
    if os.path.exists(paths["samples"]):
        with open(paths["samples"]) as fh:
            samples = [line.strip() for line in fh if line.strip()]
    return {
        "path": profile_dir,
        "hashes": hashes,
        "agg_counts": agg_counts,
        "samples": samples,
    }


def save_profile(
    profile_dir: str,
    hashes: set[int],
    samples: Iterable[str] = (),
    agg_counts: dict[int, int] | None = None,
) -> None:
    os.makedirs(profile_dir, exist_ok=True)
    paths = _profile_paths(profile_dir)
    save_hash_db(paths["db_gz"], hashes)
    with open(paths["samples"], "w") as fh:
        for sample in sorted(set(samples)):
            fh.write(sample + "\n")
    with open(paths["agg"], "w", newline="") as fh:
        writer = csv.writer(fh)
        for h, count in sorted((agg_counts or {}).items()):
            writer.writerow([h, count])


def describe_profile(profile_dir: str) -> dict:
    profile = load_profile(profile_dir)
    return {
        "profile_dir": profile_dir,
        "n_hashes": len(profile["hashes"]),
        "n_agg_hashes": len(profile["agg_counts"]),
        "n_samples": len(profile["samples"]),
        "samples": profile["samples"],
    }


def compare_profiles(profile_a: str, profile_b: str) -> dict:
    a = load_profile(profile_a)
    b = load_profile(profile_b)
    hashes_a = a["hashes"]
    hashes_b = b["hashes"]
    shared = hashes_a & hashes_b
    union = hashes_a | hashes_b
    return {
        "profile_a": profile_a,
        "profile_b": profile_b,
        "n_hashes_a": len(hashes_a),
        "n_hashes_b": len(hashes_b),
        "n_shared_hashes": len(shared),
        "n_a_only": len(hashes_a - hashes_b),
        "n_b_only": len(hashes_b - hashes_a),
        "jaccard": len(shared) / len(union) if union else 0.0,
    }


def _write_json_or_print(payload: dict, out_path: str | None) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True)
    if out_path:
        with open(out_path, "w") as fh:
            fh.write(text + "\n")
    else:
        print(text)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Manage OutliMer hash profiles.")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    describe = sub.add_parser("describe", help="Describe a profile directory")
    describe.add_argument("profile_dir")
    describe.add_argument("--json-out")

    compare = sub.add_parser("compare", help="Compare two profile directories")
    compare.add_argument("profile_a")
    compare.add_argument("profile_b")
    compare.add_argument("--json-out")

    build = sub.add_parser("build", help="Build a profile from a hash DB")
    build.add_argument("profile_dir")
    build.add_argument("--hash-db", required=True)
    build.add_argument("--sample", action="append", default=[])

    update = sub.add_parser("update", help="Merge hashes into an existing profile")
    update.add_argument("profile_dir")
    update.add_argument("--hash-db", required=True)
    update.add_argument("--sample", action="append", default=[])

    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.command == "describe":
        _write_json_or_print(describe_profile(args.profile_dir), args.json_out)
        return 0

    if args.command == "compare":
        _write_json_or_print(
            compare_profiles(args.profile_a, args.profile_b),
            args.json_out,
        )
        return 0

    if args.command == "build":
        hashes = load_hash_db(args.hash_db)
        save_profile(args.profile_dir, hashes, samples=args.sample)
        print(f"Wrote profile with {len(hashes)} hashes to {args.profile_dir}")
        return 0

    if args.command == "update":
        profile = load_profile(args.profile_dir)
        new_hashes = load_hash_db(args.hash_db)
        samples = set(profile["samples"]) | set(args.sample)
        hashes = set(profile["hashes"]) | set(new_hashes)
        save_profile(
            args.profile_dir,
            hashes,
            samples=samples,
            agg_counts=profile["agg_counts"],
        )
        print(f"Updated profile with {len(hashes)} hashes at {args.profile_dir}")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
