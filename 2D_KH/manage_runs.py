#!/usr/bin/env python3
"""
Check whether every parameter file has a complete dataset and submit missing/incomplete runs.

Usage examples:
  ./manage_runs.py                    # print status table for every params_cases/*.nml
  ./manage_runs.py --submit           # resubmit missing/incomplete runs via 2D_KH_mpi.sh
  ./manage_runs.py --dry-run --submit # show qsub commands without running them
  ./manage_runs.py --include-unknown --submit # also re-run entries lacking tend/dtout info
  ./manage_runs.py --submit --no-chain # submit all jobs at once (default chains them)
  ./manage_runs.py --summary-json summary.json --summary-csv summary.csv # persist status tables
"""
from __future__ import annotations

import argparse
import csv
import json
import shlex
import subprocess
import sys
from dataclasses import dataclass
from decimal import ROUND_FLOOR, Decimal, InvalidOperation, localcontext
from pathlib import Path
import re
import time
from typing import Sequence

FIELD_RX = re.compile(r"field-(\d+)\.dat$")


@dataclass
class RunInfo:
    path: Path
    field_count: int
    max_index: int | None
    zero_sized: bool

    def is_complete(self, expected: int | None) -> bool:
        return (
            expected is not None
            and self.field_count >= expected
            and self.max_index is not None
            and self.max_index >= expected - 1
            and not self.zero_sized
        )


@dataclass
class ParamStatus:
    param_file: Path
    base_name: str
    expected_fields: int | None
    actual_fields: int | None
    status: str
    detail: str
    run_dir: Path | None


def parse_namelist(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line or line.startswith(("&", "/")) or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip().lower()] = value.strip()
    return values


def strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] in ("'", '"') and value[-1] == value[0]:
        return value[1:-1]
    return value


def compute_expected_snapshots(values: dict[str, str]) -> int | None:
    if "tend" not in values or "dtout" not in values:
        return None
    try:
        tend = Decimal(strip_quotes(values["tend"]))
        dtout = Decimal(strip_quotes(values["dtout"]))
    except (InvalidOperation, ValueError):
        return None
    if dtout <= 0 or tend < 0:
        return None
    with localcontext() as ctx:
        ctx.prec = 28
        ratio = (tend / dtout).to_integral_value(rounding=ROUND_FLOOR)
    expected = int(ratio) + 1
    return expected if expected > 0 else 1


def scan_run_dir(run_dir: Path) -> RunInfo:
    indices: list[int] = []
    zero = False
    for field_path in sorted(run_dir.glob("field-*.dat")):
        match = FIELD_RX.match(field_path.name)
        if not match:
            continue
        indices.append(int(match.group(1)))
        try:
            if field_path.stat().st_size == 0:
                zero = True
        except OSError:
            zero = True
    return RunInfo(
        path=run_dir,
        field_count=len(indices),
        max_index=max(indices) if indices else None,
        zero_sized=zero,
    )


def evaluate_param(param_file: Path, data_root: Path) -> ParamStatus:
    values = parse_namelist(param_file)
    expected = compute_expected_snapshots(values)
    base_name = param_file.stem
    run_dirs = sorted(data_root.glob(f"{base_name}_*"))
    if not run_dirs:
        detail = f"no {base_name}_* directory found under {data_root}"
        return ParamStatus(param_file, base_name, expected, None, "missing", detail, None)

    latest_record: tuple[Path, RunInfo] | None = None
    for run_dir in reversed(run_dirs):
        info = scan_run_dir(run_dir)
        if info.is_complete(expected):
            detail = f"{run_dir.name}: {info.field_count}/{expected} field files"
            return ParamStatus(param_file, base_name, expected, info.field_count, "ok", detail, run_dir)
        if latest_record is None:
            latest_record = (run_dir, info)

    run_dir, info = latest_record  # type: ignore[misc]
    if expected is None:
        detail = f"{run_dir.name}: cannot compute expected outputs (check tend/dtout)"
        return ParamStatus(param_file, base_name, expected, info.field_count, "unknown", detail, run_dir)

    detail = f"{run_dir.name}: {info.field_count}/{expected} files"
    if info.max_index is not None:
        detail += f", max index {info.max_index:05d}"
    if info.zero_sized:
        detail += ", zero-sized snapshot detected"
    return ParamStatus(param_file, base_name, expected, info.field_count, "incomplete", detail, run_dir)


def summarize(statuses: Sequence[ParamStatus]) -> None:
    header = f"{'parameter file':30} {'status':10} details"
    print(header)
    print("-" * len(header))
    for status in statuses:
        name = status.param_file.name
        print(f"{name:30} {status.status:10} {status.detail}")


def status_to_record(status: ParamStatus) -> dict[str, str | int | None]:
    return {
        "param_file": str(status.param_file),
        "base_name": status.base_name,
        "status": status.status,
        "expected_fields": status.expected_fields,
        "actual_fields": status.actual_fields,
        "detail": status.detail,
        "run_dir": str(status.run_dir) if status.run_dir else None,
    }


def write_summary_csv(path: Path, statuses: Sequence[ParamStatus]) -> None:
    records = [status_to_record(s) for s in statuses]
    fieldnames = [
        "param_file",
        "base_name",
        "status",
        "expected_fields",
        "actual_fields",
        "detail",
        "run_dir",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    print(f"[INFO] Wrote CSV summary to {path}")


def write_summary_json(path: Path, statuses: Sequence[ParamStatus]) -> None:
    records = [status_to_record(s) for s in statuses]
    with path.open("w") as fh:
        json.dump(records, fh, indent=2)
    print(f"[INFO] Wrote JSON summary to {path}")


def wait_for_job(job_id: str, interval: int) -> None:
    """Block until qstat no longer sees job_id."""
    if interval <= 0:
        interval = 30
    print(f"[INFO] Waiting for job {job_id} to finish...")
    while True:
        result = subprocess.run(
            ["qstat", job_id],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if result.returncode == 0:
            time.sleep(interval)
        else:
            break


def submit_needed(
    statuses: Sequence[ParamStatus],
    job_script: Path,
    extra_args: Sequence[str],
    include_unknown: bool,
    dry_run: bool,
    chain: bool,
    wait_interval: int,
) -> None:
    eligible = {"missing", "incomplete"}
    if include_unknown:
        eligible.add("unknown")
    targets = [s for s in statuses if s.status in eligible]
    if not targets:
        print("[INFO] Nothing to submit.")
        return

    for status in targets:
        env_arg = f"PARAM_LIST={status.param_file.resolve()}"
        cmd = ["qsub", *extra_args, "-v", env_arg, str(job_script)]
        printable = " ".join(shlex.quote(part) for part in cmd)
        if dry_run:
            print(f"[DRY] {printable}")
            continue
        print(f"[INFO] Submitting {status.param_file.name} ...")
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        response = result.stdout.strip()
        print(f"[INFO]  -> {response}")
        if chain and response:
            job_id = response.split()[0]
            if job_id:
                wait_for_job(job_id, wait_interval)
            else:
                print("[WARN] Could not parse job ID from qsub output; skipping wait.")


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=script_dir, help="Base directory that contains params_cases, data, and the job script (default: folder containing this file).")
    parser.add_argument("--param-dir", default="params_cases", help="Relative directory with parameter files.")
    parser.add_argument("--data-dir", default="data", help="Relative directory with simulation outputs.")
    parser.add_argument("--job-script", default="2D_KH_mpi.sh", help="PBS script used to run ap.out.")
    parser.add_argument("--qsub-extra", action="append", default=[], help="Extra tokens appended to the qsub command (repeat if needed).")
    parser.add_argument("--submit", action="store_true", help="Submit qsub jobs for every missing/incomplete parameter file.")
    parser.add_argument("--include-unknown", action="store_true", help="Allow --submit to re-run entries whose completeness could not be determined.")
    parser.add_argument("--dry-run", action="store_true", help="Print qsub commands instead of executing them.")
    parser.add_argument("--no-chain", action="store_true", help="Submit all qsub jobs immediately without waiting for the previous one to finish.")
    parser.add_argument("--wait-interval", type=int, default=60, help="Seconds between qstat polls while chaining jobs (default: 60).")
    parser.add_argument("--summary-json", type=Path, help="Write the status table to this JSON file.")
    parser.add_argument("--summary-csv", type=Path, help="Write the status table to this CSV file.")
    args = parser.parse_args()

    root = args.root.resolve()
    param_dir = (root / args.param_dir).resolve()
    data_dir = (root / args.data_dir).resolve()
    job_script = (root / args.job_script).resolve()

    if not param_dir.is_dir():
        sys.exit(f"[ERROR] Parameter directory {param_dir} not found.")
    if not job_script.is_file():
        sys.exit(f"[ERROR] Job script {job_script} not found.")

    param_files = sorted(param_dir.glob("*.nml"))
    if not param_files:
        sys.exit(f"[ERROR] No *.nml files found under {param_dir}.")

    statuses = [evaluate_param(p, data_dir) for p in param_files]
    summarize(statuses)

    should_submit = args.submit or args.dry_run
    final_statuses = statuses
    if should_submit:
        extra_tokens: list[str] = []
        for item in args.qsub_extra:
            extra_tokens.extend(shlex.split(item))
        submit_needed(
            statuses=statuses,
            job_script=job_script,
            extra_args=extra_tokens,
            include_unknown=args.include_unknown,
            dry_run=args.dry_run,
            chain=not args.no_chain,
            wait_interval=args.wait_interval,
        )
        final_statuses = [evaluate_param(p, data_dir) for p in param_files]

    if args.summary_json:
        write_summary_json(args.summary_json, final_statuses)
    if args.summary_csv:
        write_summary_csv(args.summary_csv, final_statuses)


if __name__ == "__main__":
    main()
