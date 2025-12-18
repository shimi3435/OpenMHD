import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import openmhd

# dummy index
vx=0;vy=1;vz=2;pr=3;ro=4;bx=5;by=6;bz=7;ps=8

def find_cases(base_dir: Path) -> list[str]:
    return sorted([p.name for p in base_dir.iterdir() if p.is_dir()])


def find_fields(case_dir: Path) -> list[Path]:
    return sorted(case_dir.glob("field-*.dat"))


def choose_case_and_field(base_dir: Path, case: str|None, index: int|None) -> Path:
    cases = find_cases(base_dir)
    if not cases:
        raise SystemExit(f"No case directories found under {base_dir}")

    if case is None:
        print("Available cases:")
        for idx, name in enumerate(cases):
            print(f"  [{idx}] {name}")
        selection = input("Enter case index: ").strip()
        if not selection.isdigit():
            raise SystemExit("Invalid input.")
        case_idx = int(selection)
        if case_idx < 0 or case_idx >= len(cases):
            raise SystemExit("Index out of range.")
        case = cases[case_idx]
    elif case not in cases:
        raise SystemExit(f"Case '{case}' not found under {base_dir}")

    case_dir = base_dir / case
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot OpenMHD data from arbitrary cases.")
    parser.add_argument("--base", default="data", help="Root directory containing case folders")
    parser.add_argument("--case", help="Case directory to use (defaults to interactive selection)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--field-index", type=int, help="Index in sorted field list")
    group.add_argument("--field-file", help="Specific field file name/path")
    parser.add_argument("--quiet", action="store_true", help="Disable interactive prompts")
    return parser.parse_args()


def main():
    args = parse_args()
    base_dir = Path(args.base).expanduser().resolve()
    if args.field_file:
        field_path = Path(args.field_file)
        if not field_path.is_file():
            raise SystemExit(f"Field file '{field_path}' not found.")
    else:
        field_path = choose_case_and_field(base_dir, args.case, args.field_index)

    print(f"[INFO] Loading {field_path}")
    x,y,t,data = openmhd.data_read(str(field_path))
    # reading the data (partial domain: [ix1,ix2] x [jx1,jx2])
    # x,y,t,data = openmhd.data_read("data/field-00015.dat",ix2=101,jx1=50,jx2=101)

    # clearing the current figure, if any
    plt.clf()
    # extent: [left, right, bottom, top]
    extent=[x[0],x[-1],y[0],y[-1]]
    # 2D plot (vmin/mymin: minimum value, vmax/mymax: max value)
    # Note: ().T is necessary for 2-D plot routines (imshow/pcolormesh...)
    tmp = np.ndarray((x.size,y.size),np.double)
    tmp[:,:] = data[:,:,ro]
    mymax = max(tmp.max(), -tmp.min()) if( tmp.max() > 0.0 ) else 0.0
    mymin = min(tmp.min(), -tmp.max()) if( tmp.min() < 0.0 ) else 0.0
    myimg = plt.imshow(tmp.T,origin='lower',vmin=mymin,vmax=mymax,cmap='bwr',extent=extent)

# image operations (e.g. colormaps)
# myimg.set_cmap('jet')
# myimg.set_cmap('RdBu_r')  # colortable(70,/reverse) in IDL
# myimg.set_cmap('seismic')
# myimg.set_cmap('bwr')
# myimg.set_cmap('gnuplot2')

# useful options
# plt.grid()
    plt.xlabel("X",size=16)
    plt.ylabel("Y",size=16)
    plt.title('Density (t = %6.1f)' % t, size=20)

# color bar (next to myimg)
    plt.colorbar()

# flow vectors
# Note: ().T is necessary for a 2-D array
    myxsub = 40; myysub = 20; myxsub0 = int(myxsub/2); myysub0 = int(myysub/2)
    myvec = plt.quiver( x[myxsub0::myxsub],
                        y[myysub0::myysub],
                        data[myxsub0::myxsub,myysub0::myysub,vx].T,
                        data[myxsub0::myxsub,myysub0::myysub,vy].T,
                        alpha=0.7,angles='xy')

    # plot
    plt.show()

    # image file
    # plt.savefig('output.png', dpi=144)

if __name__ == "__main__":
    main()
