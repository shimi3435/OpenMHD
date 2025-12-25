#!/usr/bin/env python3
"""
Create multiple parameter files by perturbing values from a template.

Example:
  ./generate_params.py --template params.nml --count 5 \
      --float-range vy_amplitude:0.6:1.2 --float-range by_wavenumber:1.0:3.0 \
      --int-range nx:128:384:32 --output-dir params_cases --prefix basic_case_
"""
from __future__ import annotations

import argparse
import random
import re
from pathlib import Path


def parse_float_range(spec: str) -> tuple[str, float, float]:
    key, min_v, max_v = spec.split(":")
    return key.strip(), float(min_v), float(max_v)


def parse_int_range(spec: str) -> tuple[str, int, int, int]:
    parts = spec.split(":")
    if len(parts) not in (3, 4):
        raise ValueError(f"Invalid int-range spec '{spec}' (expected key:min:max[:step])")
    key = parts[0].strip()
    min_v = int(parts[1])
    max_v = int(parts[2])
    step = int(parts[3]) if len(parts) == 4 else 1
    if step <= 0:
        raise ValueError("step must be positive")
    return key, min_v, max_v, step


def parse_enum(spec: str) -> tuple[str, list[str]]:
    key, values = spec.split(":")
    return key.strip(), [v.strip() for v in values.split("|") if v.strip()]


def read_template(path: Path) -> list[str]:
    return path.read_text().splitlines(True)


def write_param_file(path: Path, lines: list[str]) -> None:
    path.write_text("".join(lines))


def update_value(lines: list[str], key: str, value: str) -> None:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=", re.IGNORECASE)
    fmt_value = value
    for idx, line in enumerate(lines):
        if pattern.match(line):
            prefix = line.split("=", 1)[0]
            lines[idx] = f"{prefix}= {fmt_value}\n"
            return
    raise KeyError(f"Key '{key}' not found in template.")


def wrap_string(value: str) -> str:
    if (value.startswith("'") and value.endswith("'")) or (
        value.startswith('"') and value.endswith('"')
    ):
        return value
    return f"'{value}'"


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", default="params.nml", help="Template params file")
    parser.add_argument("--output-dir", default="params_cases", help="Directory for new files")
    parser.add_argument("--count", type=int, default=5, help="Number of parameter files to generate")
    parser.add_argument("--prefix", default="auto_case_", help="Prefix for file base names")
    parser.add_argument("--start-index", type=int, default=1, help="Starting index for naming")
    parser.add_argument("--seed", type=int, help="Random seed (optional)")
    parser.add_argument("--runlist", help="Optional file to list generated parameter paths")
    parser.add_argument(
        "--float-range",
        action="append",
        default=[],
        metavar="KEY:MIN:MAX",
        help="Random float range for a key (can be repeated)",
    )
    parser.add_argument(
        "--int-range",
        action="append",
        default=[],
        metavar="KEY:MIN:MAX[:STEP]",
        help="Random integer range for a key (can be repeated)",
    )
    parser.add_argument(
        "--enum",
        action="append",
        default=[],
        metavar="KEY:VALA|VALB|...",
        help="Choose randomly from fixed values for a key (can be repeated)",
    )
    parser.add_argument(
        "--fixed",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Set a fixed value for every file (quoted automatically if needed)",
    )
    return parser


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    template_path = Path(args.template)
    if not template_path.exists():
        raise SystemExit(f"Template {template_path} not found.")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.seed is not None:
        random.seed(args.seed)

    float_specs = [parse_float_range(spec) for spec in args.float_range]
    int_specs = [parse_int_range(spec) for spec in args.int_range]
    enum_specs = [parse_enum(spec) for spec in args.enum]

    fixed_pairs = []
    for item in args.fixed:
        if "=" not in item:
            raise ValueError(f"Invalid fixed spec '{item}', expected KEY=VALUE.")
        key, value = item.split("=", 1)
        fixed_pairs.append((key.strip(), value.strip()))

    runlist_lines = []
    template_lines = read_template(template_path)
    idx = args.start_index
    for n in range(args.count):
        lines = template_lines.copy()
        case_name = f"{args.prefix}{idx+n:04d}"
        case_dir = Path("data") / case_name
        update_value(lines, "output_dir", wrap_string(str(case_dir)))

        for key, lo, hi in float_specs:
            value = random.uniform(lo, hi)
            update_value(lines, key, f"{value:.8f}")

        for key, lo, hi, step in int_specs:
            choices = list(range(lo, hi + 1, step))
            value = random.choice(choices)
            update_value(lines, key, str(value))

        for key, choices in enum_specs:
            if not choices:
                continue
            value = random.choice(choices)
            if not re.match(r"^[\d\.\-+eE]+$", value):
                value = wrap_string(value)
            update_value(lines, key, value)

        for key, value in fixed_pairs:
            fmt = value
            if not re.match(r"^[\d\.\-+eE]+$", value):
                fmt = wrap_string(value)
            update_value(lines, key, fmt)

        params_path = out_dir / f"{case_name}.nml"
        write_param_file(params_path, lines)
        runlist_lines.append(str(params_path))
        print(f"[INFO] wrote {params_path}")

    if args.runlist:
        Path(args.runlist).write_text("\n".join(runlist_lines) + "\n")
        print(f"[INFO] saved list of files to {args.runlist}")


if __name__ == "__main__":
    main()
