# mcal: 有機半導体結晶の移動度テンソル計算プログラム
[![Python](https://img.shields.io/badge/python-3.9%20or%20newer-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![docs](https://img.shields.io/badge/docs-here-11419572)](https://matsui-lab-yamagata.github.io/mcal/)


# 概要
`mcal`は有機半導体の移動度テンソルを計算するツールです。結晶構造から移動積分と再配列エネルギーを計算し、異方性と経路の連続性を考慮して移動度テンソルを算出します。

# 必要環境
* Python 3.9以降
* NumPy
* Pandas
* Matplotlib
* yu-tcal==4.0.2

## 量子化学計算ツール
以下のうちいずれか一つが必要です：
* Gaussian 09または16
* PySCF (macOS / Linux / WSL2(Windows Subsystem for Linux))
* GPU4PySCF (macOS / Linux / WSL2(Windows Subsystem for Linux))

# 注意事項
* Gaussianのパスが設定されている必要があります。
* PySCFはmacOS / Linux / WSL2(Windows Subsystem for Linux)でのみサポートされています。

# インストール
## Gaussian 09または16を使用する場合（PySCFなし）
```
pip install yu-mcal
```

## PySCFを使用する場合（CPUのみ、macOS / Linux / WSL2）
```
pip install "yu-mcal[pyscf]"
```

## PySCFでGPUアクセラレーションを使用する場合（macOS / Linux / WSL2）
### 1. インストールされているCUDA Toolkitのバージョンを確認する
```
nvcc --version
```

### 2. GPUアクセラレーション対応のmcalをインストールする
CUDA Toolkitのバージョンが12.xの場合:
```
pip install "yu-mcal[gpu4pyscf-cuda12]"
```
CUDA Toolkitのバージョンが11.xの場合:
```
pip install "yu-mcal[gpu4pyscf-cuda11]"
```


## インストールの確認

インストール後、以下のコマンドで確認できます：

```bash
mcal --help
```

# mcal 使用マニュアル

## 基本的な使用方法

```bash
mcal <cif_filename or pkl_filename> <osc_type> [オプション]
```

### 必須引数

- `cif_filename`: CIFファイルのパス
- `pkl_filename`: pickleファイルのパス
- `osc_type`: 何型の有機半導体として計算するか
  - `p`: p型半導体 (移動積分にHOMOを使用)
  - `n`: n型半導体 (移動積分にLUMOを使用)

### 基本例

```bash
# p型半導体として計算
mcal xxx.cif p

# n型半導体として計算
mcal xxx.cif n
```

## オプション

### 計算設定

#### `-M, --method <method>`
計算手法を指定します（GaussianおよびPySCF共通）。
- **デフォルト**: `B3LYP/6-31G(d,p)`
- **例**: `mcal xxx.cif p -M "B3LYP/6-31G(d)"`

#### `-c, --cpu <number>`
使用するCPU数を指定します。
- **デフォルト**: `4`
- **例**: `mcal xxx.cif p -c 8`

#### `-m, --mem <memory>`
メモリ量をGB単位で指定します。
- **デフォルト**: `10`
- **例**: `mcal xxx.cif p -m 16`

#### `-g, --g09`
Gaussian 09を使用します（デフォルトはGaussian 16）。
- **例**: `mcal xxx.cif p -g`

### PySCF設定

#### `--pyscf`
Gaussianの代わりにPySCFを使用して計算します。`yu-mcal[pyscf]`のインストールが必要です。
- **例**: `mcal xxx.cif p --pyscf`

#### `--gpu4pyscf`
gpu4pyscfによるGPUアクセラレーションを使用します。`--pyscf`を指定しなくても自動的にPySCFモードが有効になります。`yu-mcal[gpu4pyscf-cuda11]`または`yu-mcal[gpu4pyscf-cuda12]`のインストールが必要です。
- **例**: `mcal xxx.cif p --gpu4pyscf`

#### `--cart`
球面調和関数の代わりにデカルト基底関数を使用します（PySCFのみ）。
- **例**: `mcal xxx.cif p --pyscf --cart`

### 計算制御

#### `-r, --read`
計算を実行せずに既存のファイルから結果を読み取ります。Gaussianの場合はlogファイル、PySCFの場合はチェックポイント（`.chk`）ファイルを読み取ります。
- **例**: `mcal xxx.cif p -r`

#### `-rp, --read_pickle`
計算を実行せずに既存のpickleファイルから結果を読み取ります。
- **例**: `mcal xxx_result.pkl p -rp`

#### `--resume`
既存の結果を使用して計算を再開します。Gaussianの場合はlogファイルの正常終了を確認し、PySCFの場合はチェックポイント（`.chk`）ファイルの存在を確認します。
- **例**: `mcal xxx.cif p --resume`

#### `--fullcal`
慣性モーメントと重心間距離を使用した高速化処理を無効にし、すべてのペアに対して移動積分を計算します。
- **例**: `mcal xxx.cif p --fullcal`

#### `--cellsize <number>`
移動積分計算のために中央単位格子の周りに各方向に拡張する単位格子数を指定します。
- **デフォルト**: `2`（5×5×5のスーパーセルを作成）
- **例**: 
  - `mcal xxx.cif p --cellsize 1`（3×3×3のスーパーセルを作成）
  - `mcal xxx.cif p --cellsize 3`（7×7×7のスーパーセルを作成）

### 出力設定

#### `-p, --pickle`
計算結果をpickleファイルに保存します。
- **例**: `mcal xxx.cif p -p`

#### `--plot-plane <plane>`
指定した結晶学的面上で移動度テンソルを2D極座標プロットで描画します。
- **利用可能な面**: `ab`, `ac`, `ba`, `bc`, `ca`, `cb`
- **デフォルト**: None（プロットは生成されません）
- **例**: 
  - `mcal xxx.cif p --plot-plane ab`（ab面にプロット）
  - `mcal xxx.cif p --plot-plane bc`（bc面にプロット）

## 使用例

### 基本的な計算
```bash
# p型xxxの移動度を計算
mcal xxx.cif p

# 8CPUと16GBメモリを使用
mcal xxx.cif p -c 8 -m 16
```

### 高精度計算
```bash
# すべてのペアに対して移動積分を計算（高精度、時間がかかる）
mcal xxx.cif p --fullcal

# より大きなスーパーセルを使用して移動積分計算範囲を拡大
mcal xxx.cif p --cellsize 3

# 異なる基底関数セットを使用
mcal xxx.cif p -M "B3LYP/6-311G(d,p)"
```

### PySCFを使用した計算
```bash
# PySCFで計算（CPU）
mcal xxx.cif p --pyscf

# PySCFでGPUアクセラレーションを使用（--pyscf不要）
mcal xxx.cif p --gpu4pyscf

# 8CPUと16GBメモリを使用してPySCFで計算
mcal xxx.cif p --pyscf -c 8 -m 16

# 中断されたPySCF計算を再開
mcal xxx.cif p --pyscf --resume

# 既存のPySCFチェックポイントファイルから読み取り
mcal xxx.cif p --pyscf -r
```

### 結果の再利用
```bash
# 既存の計算結果から読み取り
mcal xxx.cif p -r

# 既存のpickleファイルから読み取り
mcal xxx_result.pkl p -rp

# 中断された計算を再開
mcal xxx.cif p --resume

# 結果をpickleファイルに保存
mcal xxx.cif p -p
```

## 出力

### 標準出力
- 再配列エネルギー
- 各ペアの移動積分
- 拡散係数テンソル
- 移動度テンソル
- 移動度の固有値と固有ベクトル

### 生成ファイル

#### 再配列エネルギーのファイル

再配列エネルギー計算では以下のファイルが生成されます（`c` = p型のカチオン、`a` = n型のアニオン）：

##### Gaussian
- `xxx_opt_n.gjf` / `xxx_opt_n.log` — 中性分子の構造最適化
- `xxx_c.gjf` / `xxx_c.log`（または `xxx_a`）— 中性構造でのイオンのSPエネルギー計算
- `xxx_opt_c.gjf` / `xxx_opt_c.log`（または `xxx_opt_a`）— イオンの構造最適化
- `xxx_n.gjf` / `xxx_n.log` — イオン構造での中性分子のSPエネルギー計算

##### PySCF
- `xxx_opt_n.xyz` / `xxx_opt_n.chk` — 中性分子の構造最適化
- `xxx_c.chk`（または `xxx_a.chk`）— 中性構造でのイオンのSPエネルギー計算
- `xxx_opt_c.xyz` / `xxx_opt_c.chk`（または `xxx_opt_a`）— イオンの構造最適化
- `xxx_n.chk` — イオン構造での中性分子のSPエネルギー計算

#### 移動積分のファイル

mcal は `(s_t_i_j_k)` という記法でファイルを生成します：

| 記号 | 意味 |
|------|------|
| `s` | 参照ユニットセル (0,0,0) 内の分子インデックス |
| `t` | 隣接ユニットセル内の分子インデックス |
| `i` | **a** 軸方向の並進インデックス |
| `j` | **b** 軸方向の並進インデックス |
| `k` | **c** 軸方向の並進インデックス |

**例：** `xxx-(0_0_1_0_0)` は (0,0,0) セルの 0 番目の分子と (1,0,0) セルの 0 番目の分子間の移動積分を表します。

##### Gaussian
- `xxx-(s_t_i_j_k).gjf` / `xxx-(s_t_i_j_k).log` — ダイマー
- `xxx-(s_t_i_j_k)_m1.gjf` / `xxx-(s_t_i_j_k)_m1.log` — モノマー1
- `xxx-(s_t_i_j_k)_m2.gjf` / `xxx-(s_t_i_j_k)_m2.log` — モノマー2

##### PySCF
- `xxx-(s_t_i_j_k).xyz` / `xxx-(s_t_i_j_k).chk` — ダイマー
- `xxx-(s_t_i_j_k)_m1.chk` — モノマー1
- `xxx-(s_t_i_j_k)_m2.chk` — モノマー2

## 注意事項

1. **計算時間**: 分子数とセルサイズによって計算時間は大きく変わります
2. **メモリ使用量**: 大きなシステムでは十分なメモリを確保してください
3. **Gaussianのインストール**: Gaussian 09またはGaussian 16が必要です
4. **依存関係**: 必要なPythonライブラリがすべてインストールされていることを確認してください

## トラブルシューティング

### 計算が途中で停止した場合
```bash
# --resumeオプションで再開
mcal xxx.cif p --resume
```

### メモリ不足エラーの場合
```bash
# メモリ量を増やす
mcal xxx.cif p -m 32
```

### 計算時間を短縮するには
```bash
# 高速化処理を有効にする（デフォルト）
mcal xxx.cif p

# より小さなスーパーセルを使用して高速計算
mcal xxx.cif p --cellsize 1

# CPU数を増やす
mcal xxx.cif p -c 16
```

# 著者
[山形大学 有機エレクトロニクス研究センター (ROEL) 松井研究室](https://matsui-lab.yz.yamagata-u.ac.jp/index.html)  
松井 弘之、尾沢 昂輝  
Email: h-matsui[at]yz.yamagata-u.ac.jp  
[at]を@に置き換えてください

# 謝辞
本研究はJSPS特別研究員奨励費、JP25KJ0647の助成を受けたものです。  

# 参考文献
[1] Qiming Sun et al., Recent developments in the PySCF program package, *J. Chem. Phys.* **2020**, *153*, 024109.  
[2] Lee-Ping Wang, Chenchen Song, Geometry optimization made simple with translation and rotation coordinates, *J. Chem. Phys.* **2016**, *144*, 214108.  
