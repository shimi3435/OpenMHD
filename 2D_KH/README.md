# 2D_KH: Kelvin-Helmholtz 不安定性テスト

## 概要
- 圧縮性 MHD 方程式を 2 次元で解く Fortran コードです。MPI 並列 (領域分割) と OpenMP に対応しており、`mainp.f90`/`modelp.f90` を中心に Kelvin-Helmholtz (KH) せん断層の成長を追跡します。
- せん断速度と密度勾配は `model.f90` で定義され、GLM 法による発散制御 (`glm_ss2.f90`) や HLLD など複数のリーマンソルバー (`flux_solver.f90`) を切り替えてベンチマークできます。
- 代表的なパラメータは `params.nml` もしくは `params_cases/*.nml` に Namelist 形式で記述します。`output_dir` 以下に `field-*.dat` や `summary.*` が保存され、`plot.py`/`plot.ipynb` や Web UI の Plotly 表示で可視化できます。

## 依存関係とビルド
1. MPI 対応 Fortran コンパイラ (例: `mpif90`) と OpenMP を利用できる環境を用意します。
2. 必要に応じて `Makefile` の `F90`/`FFLAGS` を編集します。
3. 並列版をビルド:
   ```bash
   cd 2D_KH
   make runp   # ap.out を生成 (run を使えば単一ノード向け a.out)
   ```
4. ビルド済みバイナリは `ap.out` (MPI+OpenMP) / `a.out` (シリアル) です。`clean` で一時ファイルを削除できます。

## パラメータ設定
- `params.nml` の主な項目:
  - `nx`, `ny`: 内部格子数。MPI 分割 (`mpi_nums`) と整合させます。
  - `domain_x`, `domain_y_min`: 計算領域。`domain_y_min` から `dx` 間隔で上端が決まります。
  - `alpha`, `vx_shear`, `vy_amplitude`: KH せん断プロファイルと摂動の強さ。
  - `bx0`, `by0`, `bz0`: 一様背景磁場。
  - `cfl`, `flux_type`, `lm_type`, `time_type`: 数値スキームの設定。
  - `dtout`, `output_every_step`, `output_dir`, `io_type`: 出力タイミングと I/O 方法 (MPI-IO か単一ファイルか)。
  - `bc_periodicity`: 各方向の境界条件。標準では x 周期 / y 非周期。
- 複数ケースを回す場合は `params_cases/*.nml` を作成し、`PARAM_LIST="params_cases/case1.nml ..."` のように渡します。`generate_params.py` や `manage_runs.py` を使うと一括生成・実行が容易です。
- PBS 経由で `RUN_DIR` を指定する場合は、`PARAM_LIST` を 1 つに絞ってください（出力先を固定するため）。

## 実行方法
### ローカル / 対話実行
```bash
mpirun -np 4 ./ap.out params.nml data/KH_run1
# もしくは PARAM_LIST を利用
PARAM_LIST="params_cases/kx1.nml params_cases/kx2.nml" ./2D_KH_mpi_PS.sh
```
- 引数は `./ap.out <param_file> <出力先ディレクトリ>` です。スクリプトでは `data/<ケース名>_<タイムスタンプ>` が自動生成され、`params.nml` が各 run_dir にコピーされます。
- `run_param_sweep.sh` は `param_runs.list` を読み込んで複数ケースを連続実行します。

### PBS ジョブ (`2D_KH_mpi_PS.sh` / `2D_KH_serial.sh`)
- NIFS B 系列を想定した qsub スクリプトです。`PARAM_LIST` が未指定の場合は `params_cases/*.nml` を列挙し、存在しなければ `params.nml` を実行します。
- `module load openmpi/5.0.7/rocm6.3.3`、`OMP_NUM_THREADS=24` を設定し、`mpirun --map-by NUMA` で 4 ノード×24 スレッド構成を採用しています。必要に応じて `#PBS` リソース行や `module load` を自分の環境に合わせて変更してください。
- 標準出力/エラーはジョブ名に紐づく `.oXXXXX` ファイルにまとめて出力されます。
- `2D_KH_serial.sh` は `a.out` を逐次実行するシリアル版です。`#PBS -l` と `OMP_NUM_THREADS` は環境に合わせて調整してください。
- Web アプリ経由の実行は、PBS チェックを有効にした場合は `params_cases` のケースを選択して投入し、無効にした場合はローカルでサブプロセス実行します。いずれも `OpenMHD/2D_KH/data/` に出力します。

## 解析・可視化
- `plot.py`, `plot.ipynb`: Python (matplotlib) を使った 2D スライスや時間発展の確認用サンプル。`PYTHONPATH` にリポジトリルートを含めるか、`%run -i plot.py` として実行します。
- `summary.csv`/`summary.json`: `manage_runs.py` が run ディレクトリのメタデータを収集した結果。後処理の整理に利用できます。

## 再実行/再開のヒント
- 計算途中で停止した場合は `data/<run>/checkpoint` などに保存したステップ番号を `params.nml` の `n_start` へ指定し、同じ `output_dir` を渡すことで再開できます。
- 解析済みデータを残したまま追加実行したい場合は `params_cases` に別名を付け、スクリプトでタイムスタンプ付きディレクトリを生成する運用が安全です。
