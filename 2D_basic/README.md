# 2D_basic: Orszag-Tang 渦テスト

## 概要
- 2D_basic は Orszag-Tang 渦を模した 2 次元 MHD テストケースです。正弦波状の速度・磁場を周期境界で発達させ、ショック形成や磁気島の生成など基本挙動を確認できます。
- 実装は Kelvin-Helmholtz 版と同等の OpenMHD コア (HLL/HLLD 流束、TVD Runge-Kutta 時間積分、GLM 発散制御) を共有し、`model.f90` で初期条件のみを切り替えています。
- Fortran + MPI + OpenMP で実装されており、MPI 並列は `mpi_nums` (x, y 分割数) で制御します。結果は `data/` 配下に step ごとの `field-*.dat` が生成されます。

## 主なファイル
- `mainp.f90` / `modelp.f90` / `mpibc.f90`: 並列時間積分と境界通信。
- `model.f90`: 背景密度・圧力 (`density0`, `pressure0`) と正弦波パターン (`vx_*`, `bx_*`, `by_*`) を定義。
- `params.nml` と `params_cases/`: Namelist 入力。複数ケースを置くとスクリプトが自動検出。
- `2D_basic_mpi.sh`: PBS 用の投入スクリプト。
- `plot.py`, `plot_movie.py`, `plot.ipynb`: Python 可視化サンプル。

## 依存関係とビルド
1. MPI 対応 Fortran コンパイラ (例: `mpif90 -fopenmp`) を用意します。
2. `Makefile` で使用するコンパイラや最適化オプションを確認します (デフォルトは `mpif90 -fopenmp -O2`)。
3. ビルド手順:
   ```bash
   cd 2D_basic
   make runp   # ap.out (MPI+OpenMP)
   # または make run でシリアル版 a.out
   ```
4. 不要な生成物は `make clean` で削除できます。

## パラメータの要点 (`params.nml`)
- `nx`, `ny`, `domain_x`, `domain_y_min`: 格子と領域 (デフォルトは 0〜2π の正方領域)。
- `density0`, `pressure0`: 一様な初期プラズマ状態。
- `vx_amplitude`, `vy_amplitude`, `bx_amplitude`, `by_amplitude`: 速度・磁場の振幅。`*_wavenumber` で波数を調整できます。
- `cfl`, `flux_type`, `lm_type`, `time_type`: 数値安定性やリーマンソルバー選択。
- `bz0`, `ps0`: z 成分や GLM 追加変数の初期値。
- `bc_periodicity = .true., .true.` なので両方向周期境界。必要に応じて変更可能です。
- 出力制御 (`dtout`, `output_dir`, `io_type`) や再起動 (`n_start`) も同ファイルで設定します。

## 実行方法
### 直接実行
```bash
mpirun -np 4 ./ap.out params.nml data/OT_ref
```
- 実行時引数は `<param_file> <出力ディレクトリ>`。存在しない場合は自動作成され、実行時点の `params.nml` がコピーされます。
- `PARAM_LIST` を指定すると複数の Namelist を一度に実行できます。

### PBS 投入 (`2D_basic_mpi.sh`)
- `#PBS -l select=4:ncpus=24:ompthreads=24` など NIFS のバッチ設定を含みます。別環境ではキュー名やノード数を調整してください。
- `module load openmpi/5.0.7/rocm6.3.3` で MPI を読み込み、`OMP_NUM_THREADS=24` を設定したのち、`mpirun --map-by NUMA -x UCX_MAX_RNDV_RAILS=4 ./ap.out ...` を実行します。
- `PARAM_LIST` 未設定時は `params_cases/*.nml` を列挙し、無ければ `params.nml` を利用します。スクリプト実行後、`data/<case>_<timestamp>` が生成され、ログは `2D_basic_mpi_PS.o*` にまとめられます。

## 後処理
- `plot.py`: 任意の `data/<run>` を読み込み、密度・圧力・磁場などを matplotlib で可視化。
- `plot_movie.py`: `ffmpeg` と組み合わせて時間発展アニメーションを作成。
- `summary.csv`/`summary.json`: `manage_runs.py` が作るランサマリ。パラメータ掃引の整理に便利です。

## Tips
- パラメータ掃引を行う場合は `param_runs.list` と `run_param_sweep.sh` を利用すると、リスト化したケースを順次実行できます。
- 高解像度で実行する際は `select` や `OMP_NUM_THREADS`、`mpi_nums` を揃えて領域が割り切れるようにしてください。
- 途中再開は `n_start` と出力ディレクトリを既存 run に合わせるだけで可能です。
