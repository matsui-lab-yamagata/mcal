# mcal: 有機半導体結晶の移動度テンソル計算プログラム
[![Python](https://img.shields.io/badge/python-3.7%20or%20newer-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# 概要
`mcal.py`は有機半導体の移動度テンソルを計算するツールです。結晶構造から移動積分と再組織化エネルギーを計算し、異方性と経路の連続性を考慮して移動度テンソルを決定します。

# 必要環境
* Python 3.7以降
* NumPy
* Pandas
* Gaussian 09または16

# 重要な注意事項
* Gaussianのパスが設定されている必要があります。

# mcal 使用マニュアル

## 基本的な使用方法

```bash
python mcal.py <cif_filename> <osc_type> [オプション]
```

### 必須引数

- `cif_filename`: CIFファイルのパス
- `osc_type`: 有機半導体の種類
  - `p`: p型半導体（HOMOレベルを使用）
  - `n`: n型半導体（LUMOレベルを使用）

### 基本例

```bash
# p型半導体として計算
python mcal.py xxx.cif p

# n型半導体として計算
python mcal.py xxx.cif n
```

## オプション

### 計算設定

#### `-M, --method <method>`
Gaussian計算で使用する計算手法を指定します。
- **デフォルト**: `B3LYP/6-31G(d,p)`
- **例**: `python mcal.py xxx.cif p -M "B3LYP/6-31G(d)"`

#### `-c, --cpu <number>`
使用するCPU数を指定します。
- **デフォルト**: `4`
- **例**: `python mcal.py xxx.cif p -c 8`

#### `-m, --mem <memory>`
メモリ量をGB単位で指定します。
- **デフォルト**: `10`
- **例**: `python mcal.py xxx.cif p -m 16`

#### `-g, --g09`
Gaussian 09を使用します（デフォルトはGaussian 16）。
- **例**: `python mcal.py xxx.cif p -g`

### 計算制御

#### `-r, --read`
Gaussianを実行せずに既存のログファイルから結果を読み取ります。
- **例**: `python mcal.py xxx.cif p -r`

#### `--resume`
ログファイルが正常に終了している場合、既存の結果を使用して計算を再開します。
- **例**: `python mcal.py xxx.cif p --resume`

#### `--fullcal`
慣性モーメントと重心間距離を使用した高速化処理を無効にし、すべてのペアに対して移動積分を計算します。
- **例**: `python mcal.py xxx.cif p --fullcal`

#### `--cellsize <number>`
移動積分計算のために中央単位格子の周りに各方向に拡張する単位格子数を指定します。
- **デフォルト**: `2`（5×5×5のスーパーセルを作成）
- **例**: 
  - `python mcal.py xxx.cif p --cellsize 1`（3×3×3のスーパーセルを作成）
  - `python mcal.py xxx.cif p --cellsize 3`（7×7×7のスーパーセルを作成）

### 出力設定

#### `-p, --pickle`
計算結果をpickleファイルに保存します。
- **例**: `python mcal.py xxx.cif p -p`

### 拡散係数計算手法

#### `--mc`
モンテカルロ法を使用して拡散係数テンソルを計算します。
- **例**: `python mcal.py xxx.cif p --mc`

#### `--pde`
偏微分方程式法を使用して拡散係数テンソルを計算します。
- **例**: `python mcal.py xxx.cif p --pde`

## 実用的な使用例

### 基本的な計算
```bash
# p型xxxの移動度を計算
python mcal.py xxx.cif p

# 8CPUと16GBメモリを使用
python mcal.py xxx.cif p -c 8 -m 16
```

### 高精度計算
```bash
# すべてのペアに対して移動積分を計算（高精度、時間がかかる）
python mcal.py xxx.cif p --fullcal

# より小さなスーパーセルを使用して高速計算
python mcal.py xxx.cif p --cellsize 1

# より大きなスーパーセルを使用して移動積分計算範囲を拡大
python mcal.py xxx.cif p --cellsize 3

# 異なる基底関数セットを使用
python mcal.py xxx.cif p -M "B3LYP/6-311G(d,p)"
```

### 結果の再利用
```bash
# 既存の計算結果から読み取り
python mcal.py xxx.cif p -r

# 中断された計算を再開
python mcal.py xxx.cif p --resume

# 結果をpickleファイルに保存
python mcal.py xxx.cif p -p
```

### 拡散係数の比較
```bash
# 通常計算 + モンテカルロ + PDE法で比較
python mcal.py xxx.cif p --mc --pde
```

## 出力

### 標準出力
- 再組織化エネルギー
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
python mcal.py xxx.cif p --resume
```

### メモリ不足エラーの場合
```bash
# メモリ量を増やす
python mcal.py xxx.cif p -m 32
```

### 計算時間を短縮するには
```bash
# 高速化処理を有効にする（デフォルト）
python mcal.py xxx.cif p

# CPU数を増やす
python mcal.py xxx.cif p -c 16
```
