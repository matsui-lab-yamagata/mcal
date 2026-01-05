# mcal: 有機半導体結晶の移動度テンソル計算プログラム
[![Python](https://img.shields.io/badge/python-3.9%20or%20newer-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![docs](https://img.shields.io/badge/docs-here-11419572)](https://matsui-lab-yamagata.github.io/mcal/)


# 概要
`mcal.py`は有機半導体の移動度テンソルを計算するツールです。結晶構造から移動積分と再配列エネルギーを計算し、異方性と経路の連続性を考慮して移動度テンソルを算出します。

# 必要環境
* Python 3.9以降
* NumPy
* Pandas
* Gaussian 09または16

# 注意事項
* Gaussianのパスが設定されている必要があります。

# インストール
```bash
pip install yu-mcal
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
Gaussianで使用する計算手法を指定します。
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

### 計算制御

#### `-r, --read`
Gaussianを実行せずに既存のlogファイルから結果を読み取ります。
- **例**: `mcal xxx.cif p -r`

#### `-rp, --read_pickle`
計算を実行せずに既存のpickleファイルから結果を読み取ります。
- **例**: `mcal xxx_result.pkl p -rp`

#### `--resume`
ログファイルが正常に終了している場合、既存の結果を使用して計算を再開します。
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

### 拡散係数計算手法

#### `--mc`
モンテカルロ法を使用して拡散係数テンソルを計算します。(テスト用)
- **例**: `mcal xxx.cif p --mc`

#### `--ode`
常微分方程式法を使用して拡散係数テンソルを計算します。(テスト用)
- **例**: `mcal xxx.cif p --ode`

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

### 拡散係数の比較
```bash
# 通常計算 + モンテカルロ法 + ODEで比較
mcal xxx.cif p --mc --ode
```

## 出力

### 標準出力
- 再配列エネルギー
- 各ペアの移動積分
- 拡散係数テンソル
- 移動度テンソル
- 移動度の固有値と固有ベクトル

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
