from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openmhd

# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

VAR_MAP = {
    "vx": vx,
    "vy": vy,
    "vz": vz,
    "pr": pr,
    "ro": ro,
    "bx": bx,
    "by": by,
    "bz": bz,
    "ps": ps,
}


def find_cases(base_dir: Path) -> list[str]:
    cases = sorted([p.name for p in base_dir.iterdir() if p.is_dir()])
    if any(base_dir.glob("field-*.dat")):
        cases.insert(0, ".")
    return cases


def find_fields(case_dir: Path) -> list[Path]:
    return sorted(case_dir.glob("field-*.dat"))


def choose_case_and_field(base_dir: Path, case: str|None, index: int|None) -> Path:
    cases = find_cases(base_dir)
    if not cases:
        raise SystemExit(f"No case directories found under {base_dir}")

    if case is None:
        print("Available cases:")
        for idx, name in enumerate(cases):
            label = base_dir if name == "." else base_dir / name
            print(f"  [{idx}] {label}")
        selection = input("Enter case index: ").strip()
        if not selection.isdigit():
            raise SystemExit("Invalid input.")
        case_idx = int(selection)
        if case_idx < 0 or case_idx >= len(cases):
            raise SystemExit("Index out of range.")
        case = cases[case_idx]
    elif case not in cases:
        raise SystemExit(f"Case '{case}' not found under {base_dir}")

    case_dir = base_dir if case == "." else (base_dir / case)
    fields = find_fields(case_dir)
    if not fields:
        raise SystemExit(f"No field files found in {case_dir}")

    if index is None:
        print("Available field files:")
        for idx, path in enumerate(fields):
            print(f"  [{idx}] {path.name}")
        selection = input("Enter field index: ").strip()
        if not selection.isdigit():
            raise SystemExit("Invalid input.")
        idx = int(selection)
    else:
        idx = index
    if idx < 0 or idx >= len(fields):
        raise SystemExit("Field index out of range.")
    return fields[idx]


def compute_vector_potential(x: np.ndarray, y: np.ndarray, bx_data: np.ndarray, by_data: np.ndarray) -> np.ndarray:
    az = np.zeros((x.size, y.size), dtype=np.double)
    fx = 0.5 * (x[1] - x[0])
    fy = 0.5 * (y[1] - y[0])
    az[0,0] = fy * bx_data[0,0] - fx * by_data[0,0]
    for j in range(1, y.size):
        az[0, j] = az[0, j-1] + fy * (bx_data[0, j-1] + bx_data[0, j])
    for i in range(1, x.size):
        az[i, :] = az[i-1, :] - fx * (by_data[i-1, :] + by_data[i, :])
    return az


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot OpenMHD data from arbitrary cases.")
    parser.add_argument("--base", default="data", help="Root directory containing case folders")
    parser.add_argument("--case", help="Case directory to use (defaults to interactive selection)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--field-index", type=int, help="Index in sorted field list")
    group.add_argument("--field-file", help="Specific field file name/path")
    parser.add_argument("--quantity", default="pr", choices=sorted(VAR_MAP.keys()), help="Variable to visualize")
    parser.add_argument("--cmap", default="jet", help="Matplotlib colormap name")
    parser.add_argument("--no-contours", action="store_true", help="Disable B-field line contours")
    return parser.parse_args()


def main():
    args = parse_args()
    base_dir = Path(args.base).expanduser().resolve()
    if not base_dir.exists():
        raise SystemExit(f"Base directory '{base_dir}' not found.")
    if args.field_file:
        field_path = Path(args.field_file).expanduser()
        if not field_path.is_file():
            raise SystemExit(f"Field file '{field_path}' not found.")
    else:
        field_path = choose_case_and_field(base_dir, args.case, args.field_index)

    print(f"[INFO] Loading {field_path}")
    x,y,t,data = openmhd.data_read(str(field_path))

    idx = VAR_MAP[args.quantity]
    plt.clf()
    extent=[x[0],x[-1],y[0],y[-1]]

    tmp = np.asarray(data[:,:,idx])
    mymax = max(tmp.max(), -tmp.min()) if tmp.max() > 0.0 else 0.0
    mymin = min(tmp.min(), -tmp.max()) if tmp.min() < 0.0 else 0.0
    myimg = plt.imshow(tmp.T,origin='lower',vmin=mymin,vmax=mymax,cmap=args.cmap,extent=extent)

    plt.xlabel("X",size=16)
    plt.ylabel("Y",size=16)
    plt.title(f"{args.quantity.upper()} (t = {t:6.2f})", size=20)
    plt.colorbar()

    if not args.no_contours:
        az = compute_vector_potential(x, y, data[:,:,bx], data[:,:,by])
        plt.contour(az.T,extent=extent,colors='w',linestyles='solid',linewidths=0.8)

    plt.show()


if __name__ == "__main__":
    main()
