# neutrinoAnalysis

📖 **[English](#english) · [日本語](#日本語)**

---

<a id="english"></a>
## English

Neutrino-flux optimization pipeline. The flux is reconstructed by χ² minimization
from the observed rates (`Ratebin`) and the response matrix (`CRmat`), and
Monte-Carlo (Neyman construction) confidence intervals are obtained for each
flux parameter.

Two interchangeable solver backends (`scipy` / `osqp`) are provided.

### Directory layout

```
README.md
1eV/                              # working directory for the 1 eV threshold
  neutrino_analysis_band.py       # main module (per-threshold copy; differs)
  CRmat/originalUnit/             # response matrices CRmat<intervals>_originalUnit.csv
  Ratebin/                        # Ratebin2 / Ratebin7
  Danny’s files/                  # theoretical flux curves
  confidence_band_<scen>_T<T>.ipynb
  comparison*.ipynb
  T3/  T20/  T300/                # per-T output roots (auto-created)
    scenario_bkg_<scen>/
      flux_*.pdf
      bands/band_bkg<scen>_idx<NNN>.json

5eV/                              # 5 eV threshold, same layout
  neutrino_analysis_fast.py       # 5eV also keeps the fast module alongside band
  neutrino_analysis_band.py
  ...
```

**Important:** run notebooks with the cwd set to `1eV/` or `5eV/`. The code reads
`CRmat/originalUnit/` etc. by relative path from there. Outputs are written
under `T<T>/scenario_bkg_<scen>/`, so threshold × T × background-scenario
combinations never overwrite each other.

### Dependencies

```bash
pip install numpy scipy matplotlib pandas numba
pip install cvxpy osqp        # required when using solver='osqp'
pip install joblib            # optional; enables scipy-backend MC parallelism
```

`cvxpy` brings CLARABEL / SCS / HiGHS, used for solver fallback and vertex
selection. Tested on Python 3.12.

### Basic usage

Run with cwd `1eV/` or `5eV/` (the data paths anchor here).

```python
# At 5eV both fast and band can be imported; at 1eV only band.
from neutrino_analysis_band import NeutrinoAnalysis

# T is a constructor argument (default 3). Output goes to T<T>/scenario_bkg_<scen>/.
a = NeutrinoAnalysis(background_scenario='flat', intervals='180',
                     GeV=0.32e16, solver='scipy', T=3)

res = a.optimize(a.data_vector)      # χ² minimization
print(res.fun / a.c)                 # χ² / c
a.plot_flux_comparison(save=True)    # optimized flux vs. theoretical models
print(a.scenario_dir, a.bands_dir)   # check the output paths
```

#### Switching solvers

```python
a.set_solver('osqp')   # fast backend (cvxpy + OSQP/CLARABEL/SCS)
```

- `scipy` — `trust-constr` with analytic Jacobian and constant Hessian. Stable.
- `osqp` — formulates χ² as a QP. The free solve is fast; fixed-parameter
  solves auto-route to CLARABEL. For underdetermined problems (params > bins),
  a HiGHS-simplex vertex selection returns a piecewise-constant (staircase)
  flux (`_OSQPBackend.vertex_select`, default `True`).

Both backends compute the same χ². For underdetermined problems the flux
*shape* is not unique, so backends may return different equally-optimal
solutions (see Caveats below).

### Δχ² scan (manual grid)

```python
a.optimize(a.data_vector)            # compute the best fit first
scan = a.scan_fixed_parameter(
    fixed_index=39, scan_range=0.95, num_points=21,
    num_pseudo_data=100, seed=42, n_jobs=1,
)
import pandas as pd
pd.DataFrame(scan)                   # included / cutoff / Δχ² per fixed value
```

`scan` is a list of dicts; `pd.DataFrame(scan)` turns it into a table directly.

### Confidence band by root finding (`neutrino_analysis_band.py`)

Instead of tuning grid endpoints by hand, locate the band edges directly with
bracketing + geometric bisection.

```python
from neutrino_analysis_band import NeutrinoAnalysis, load_band

a = NeutrinoAnalysis(background_scenario='flat', intervals='360',
                     GeV=0.32e16, solver='osqp', T=3)
a.optimize(a.data_vector)

# locate the 1σ / 90% / 2σ edges for one index in a single call
band = a.find_confidence_band(
    fixed_index=0,
    levels=(0.678, 0.90, 0.954),
    num_pseudo_data=50,   # bracketing stage (coarse)
    n_pseudo_edge=500,    # bisection at the edges (larger → more stable, slower)
    step=1.5, rel_tol=0.03, seed=42, verbose=True,
)
```

The outer reach is bracketed with the widest level (2σ); all narrower-level
edges nest inside that interval. Use larger `n_pseudo_edge` near the edges to
suppress Monte-Carlo noise in the cutoff. Fix `seed` so each evaluation is
reproducible (this stops the bisection from jittering).

#### Save per index and overlay on the optimize result

```python
for idx in [0, 5, 10, 20, 40]:
    a.find_and_save_band(idx,
                         num_pseudo_data=50, n_pseudo_edge=500,
                         step=1.5, rel_tol=0.03, seed=42)
# default outdir = a.bands_dir = 'T3/scenario_bkg_flat/bands'

# overlay the saved bands on the optimized flux
a.plot_flux_with_bands(f'{a.bands_dir}/band_*.json',
                       levels=(0.678, 0.90, 0.954), ylim=(0, 3e13))
```

- `save_band` / `load_band` — JSON, one file per index
- `find_and_save_band` — search one index and save immediately
- `plot_flux_with_bands` — overlay per-index asymmetric error bars on the
  scatter (centres = saved best fits, so they line up with the optimize result)

#### Multi-scenario comparison

```python
a_flat = NeutrinoAnalysis(background_scenario='flat', intervals='180',
                          GeV=0.32e16, solver='osqp', T=3); a_flat.optimize(a_flat.data_vector)
a_a    = NeutrinoAnalysis(background_scenario='a',    intervals='180',
                          GeV=0.32e16, solver='osqp', T=3); a_a.optimize(a_a.data_vector)

a_flat.plot_band_comparison(
    {'flat': f'{a_flat.bands_dir}/band_*.json',
     'a':    f'{a_a.bands_dir}/band_*.json'},
    level=0.954,
    optimized={'flat': a_flat, 'a': a_a},   # overlay each scenario's optimized flux
    ylim=(0, 3e13), save=True,
)
```

- Each scenario uses its own colour. Bands are in physical units
  (cm⁻² s⁻¹) and are comparable regardless of the GeV unit choice.
- `optimized` takes `{label: NeutrinoAnalysis or flux array}` and scatters each
  scenario's optimized flux in the matching colour (the labels should match
  those in `groups`).

### Band-search parameter tuning

| Argument | Role | Recommended |
|---|---|---|
| `num_pseudo_data` | pseudo-data sample size during bracketing | enough to resolve the requested percentile; **≥50 for 2σ** |
| `n_pseudo_edge` | pseudo-data sample size during edge bisection (cutoff precision) | **≥500** (seed-to-seed noise converges) |
| `step` | bracketing expansion factor (reach = `v0 · step^max_bracket`) | 1.5. Too small (e.g. 1.05) with default `max_bracket=25` cannot reach a far upper edge and returns `inf` |
| `rel_tol` | relative tolerance of the edge (its resolution) | 0.03. **Sub-`rel_tol` differences carry no physical meaning** |
| `seed` | RNG seed (bisection reproducibility) | fix it |

Notes:
- Edge values carry MC noise from the cutoff. Small `n_pseudo_edge` makes the
  2σ width fluctuate by a few percent between seeds (90% is more central and
  stable; **2σ is an extreme percentile and therefore noisier**).
- If a scenario-to-scenario difference is smaller than `rel_tol`, its sign is
  not physical. To claim a difference, report **edge ± uncertainty from a
  seed ensemble** and show that the difference exceeds the uncertainty.
- When the upper edge is weakly constrained (shallow profile), report a
  **one-sided limit** rather than chasing precision.

### Backgrounds

`background_scenario` ∈ `{'a', 'b', 'b2', 'c', 'flat', 'none'}`. `set_background()`
switches it on the fly (the OSQP cache is cleared automatically).

### Caveats

- **Underdeterminacy.** With `intervals='180'` there are 180 flux parameters
  versus 29 observed bins (rank 29). The χ² minimizer is a 151-dimensional
  face, so the flux itself is not uniquely determined. `scipy` and `osqp` may
  return different points on that face. **The Δχ² profile (and the confidence
  interval) depends only on χ² and is unaffected**; physical conclusions are
  stable.
- **OSQP scaling.** The natural-unit constant `c` is ~1e70, so the OSQP
  backend column-scales internally. Without that, OSQP/CLARABEL can return a
  wrong point still labelled `optimal`.
- **Band search assumption.** The band is taken to be connected (one crossing
  on each side).

### Backend differences

| | scipy | osqp |
|---|---|---|
| Solver | trust-constr | OSQP → CLARABEL → SCS (CLARABEL first for fixed) |
| Free solve | stable | fast |
| Fixed-parameter solve | stable | CLARABEL (OSQP does not converge at large `n`) |
| MC parallelism (`n_jobs`) | works | serial |
| Flux shape | algorithm-dependent vertex | explicit piecewise-constant vertex via HiGHS |

---

<a id="日本語"></a>
## 日本語

ニュートリノフラックスの最適化パイプライン。χ² 最小化で観測レート（`Ratebin`）と
応答行列（`CRmat`）からフラックスを再構成し、モンテカルロ（Neyman 構成）で
各フラックスパラメータの信頼区間を求める。

2 つの交換可能なソルバーバックエンド（`scipy` / `osqp`）を持つ。

### ディレクトリ構成

```
README.md
1eV/                              # 検出器しきい値 1 eV の作業ディレクトリ
  neutrino_analysis_band.py       # 本体（しきい値別に独立コピー、差分あり）
  CRmat/originalUnit/             # 応答行列 CRmat<intervals>_originalUnit.csv
  Ratebin/                        # Ratebin2 / Ratebin7
  Danny’s files/                  # 理論フラックス曲線
  confidence_band_<scen>_T<T>.ipynb
  comparison*.ipynb
  T3/  T20/  T300/                # T 別の出力ルート (実行時に自動生成)
    scenario_bkg_<scen>/
      flux_*.pdf
      bands/band_bkg<scen>_idx<NNN>.json

5eV/                              # しきい値 5 eV、同じ構成
  neutrino_analysis_fast.py       # 5eV 側のみ fast 版も同居
  neutrino_analysis_band.py
  ...（1eV と同じレイアウト）
```

**重要:** notebook の実行は **`1eV/` または `5eV/` をカレントディレクトリに**して行う。
コードはそこから `CRmat/originalUnit/` などを相対パスで読み込む。
出力は `T<T>/scenario_bkg_<scen>/` 配下に自動で振り分けられるため、しきい値・T・背景シナリオ
の組み合わせが互いに上書きすることはない。

### 依存パッケージ

```bash
pip install numpy scipy matplotlib pandas numba
pip install cvxpy osqp        # solver='osqp' を使う場合に必須
pip install joblib            # scipy バックエンドのモンテカルロ並列化 (任意)
```

`cvxpy` は CLARABEL / SCS / HiGHS も同梱する（バックエンドのフォールバックと頂点選択に使用）。
Python 3.12 で動作確認。

### 基本的な使い方

cwd を `1eV/` または `5eV/` にして実行する（データの相対パスがそこに紐づく）。

```python
# 5eV では fast / band どちらも import 可能。1eV では band のみ。
from neutrino_analysis_band import NeutrinoAnalysis

# T を引数で渡す（既定 T=3）。出力先は自動で T<T>/scenario_bkg_<scen>/...
a = NeutrinoAnalysis(background_scenario='flat', intervals='180',
                     GeV=0.32e16, solver='scipy', T=3)

res = a.optimize(a.data_vector)      # χ² 最小化
print(res.fun / a.c)                 # χ²/c
a.plot_flux_comparison(save=True)    # 最適フラックス vs 理論曲線
print(a.scenario_dir, a.bands_dir)   # 出力先の確認
```

#### ソルバーの切り替え

```python
a.set_solver('osqp')   # 高速バックエンド (cvxpy + OSQP/CLARABEL/SCS)
```

- `scipy` — `trust-constr` を解析的ヤコビアン・定数ヘッシアンで解く。安定。
- `osqp` — χ² を二次計画として解く。フリー解は高速。固定パラメータ解は内部で
  CLARABEL に切り替わる。劣決定（パラメータ数 > ビン数）の場合は HiGHS シンプレックスで
  区分定数（階段状）の頂点解を選ぶ（`_OSQPBackend.vertex_select`、既定 True）。

両バックエンドが返す χ² は同一定義。劣決定の場合フラックス *形状* は一意でないため、
バックエンドにより異なる等価最適解になりうる（詳細は下記「注意点」）。

### Δχ² スキャン（手動グリッド）

```python
a.optimize(a.data_vector)            # 先にベストフィットを求める
scan = a.scan_fixed_parameter(
    fixed_index=39, scan_range=0.95, num_points=21,
    num_pseudo_data=100, seed=42, n_jobs=1,
)
import pandas as pd
pd.DataFrame(scan)                   # 各固定値の included / cutoff / Δχ²
```

`scan` は辞書のリスト。`pd.DataFrame(scan)` でそのまま表になる。

### 信頼バンドの根探索（`neutrino_analysis_band.py`）

グリッドの上下限を手で詰める代わりに、ブラケット＋幾何二分法でバンドの端を直接探す。

```python
from neutrino_analysis_band import NeutrinoAnalysis, load_band

a = NeutrinoAnalysis(background_scenario='flat', intervals='360',
                     GeV=0.32e16, solver='osqp', T=3)
a.optimize(a.data_vector)

# 1 インデックスで 1σ/90%/2σ を同時に求める
band = a.find_confidence_band(
    fixed_index=0,
    levels=(0.678, 0.90, 0.954),
    num_pseudo_data=50,   # ブラケット段階 (粗い)
    n_pseudo_edge=500,    # 端の二分法 (大きいほど安定・遅い)
    step=1.5, rel_tol=0.03, seed=42, verbose=True,
)
```

最も広い 2σ で外側をブラケットし、その区間に全レベルの端が入れ子で収まる。
端付近のみ `n_pseudo_edge` を増やして cutoff のモンテカルロノイズを抑える。
`seed` 固定で各点を再現可能にし、二分法のジッターを防ぐ。

#### 複数インデックスを保存して重ねる

```python
for idx in [0, 5, 10, 20, 40]:
    a.find_and_save_band(idx,
                         num_pseudo_data=50, n_pseudo_edge=500,
                         step=1.5, rel_tol=0.03, seed=42)
# 既定の保存先 = a.bands_dir = 'T3/scenario_bkg_flat/bands'

# 保存したバンドを optimize 結果の上に重ねる
a.plot_flux_with_bands(f'{a.bands_dir}/band_*.json',
                       levels=(0.678, 0.90, 0.954), ylim=(0, 3e13))
```

- `save_band` / `load_band` … バンドを 1 インデックス 1 ファイルの JSON で保存・読込
- `find_and_save_band` … 1 インデックスを探索して即保存
- `plot_flux_with_bands` … 保存済みバンドを各インデックスの誤差棒として散布図に重ねる
  （バンド中心は保存時の `self.result.x[index]` なので最適フラックスと一致）

#### 複数シナリオのバンドを重ねて比較

```python
a_flat = NeutrinoAnalysis(background_scenario='flat', intervals='180',
                          GeV=0.32e16, solver='osqp', T=3); a_flat.optimize(a_flat.data_vector)
a_a    = NeutrinoAnalysis(background_scenario='a',    intervals='180',
                          GeV=0.32e16, solver='osqp', T=3); a_a.optimize(a_a.data_vector)

a_flat.plot_band_comparison(
    {'flat': f'{a_flat.bands_dir}/band_*.json',
     'a':    f'{a_a.bands_dir}/band_*.json'},
    level=0.954,
    optimized={'flat': a_flat, 'a': a_a},   # 各シナリオの optimize 結果も重ねる
    ylim=(0, 3e13), save=True,
)
```

- 各シナリオを別色の誤差棒で重ねる。バンドは物理単位なので GeV 単位の違いに依らず比較可能。
- `optimized` に `{ラベル: NeutrinoAnalysis または flux 配列}` を渡すと、最適フラックス散布図を
  バンドと同色で重ねる（ラベルを `groups` と一致させると同色になる）。

### バンド探索パラメータの指針

| 引数 | 役割 | 推奨 |
|---|---|---|
| `num_pseudo_data` | ブラケット段階の擬似データ数 | 要求パーセンタイルを表現できる数。2σ なら **≥50** |
| `n_pseudo_edge` | 端の二分法での擬似データ数（cutoff 精度） | **≥500**（seed 依存ノイズが収束） |
| `step` | ブラケットの拡大率（到達範囲 `v0·step^max_bracket`） | **1.5**。小さすぎ（例 1.05）+ 既定 `max_bracket=25` だと遠い上端に届かず端が `inf` になる |
| `rel_tol` | 端の相対許容幅（端の分解能） | 0.03。**これ以下の桁の差は意味を持たない** |
| `seed` | 乱数固定（二分法の再現性） | 固定する |

注意:
- 端の値は cutoff のモンテカルロ揺らぎを持つ。小さい `n_pseudo_edge` では seed を変えるだけで
  2σ 幅が数 % ブレる（90% は中心寄りで安定、**2σ は極端な分位点なので特にブレやすい**）。
- シナリオ間の差が `rel_tol`（端の分解能）より小さい場合、その大小に物理的意味はない。
  差を主張するなら、**seed アンサンブルで端 ± 誤差を出し、差が誤差を超えることを示す**。
- 上端が弱くしか制約されない（プロファイルが浅い）場合は、精度を上げるより**片側極限として報告**する方が筋が良い。

### 背景シナリオ

`background_scenario` は `'a' / 'b' / 'b2' / 'c' / 'flat' / 'none'` から選ぶ。
`set_background()` で実行中に切り替え可能（OSQP のキャッシュは自動でクリアされる）。

### 注意点

- **劣決定性**: 例えば `intervals='180'` はフラックス 180 パラメータに対し観測は 29 ビン
  （rank 29）。χ² 最小解は 151 次元の面となり、フラックスは一意に定まらない。
  scipy / osqp は等価最適解の中の異なる点を返しうる。**Δχ² プロファイル（信頼区間）は
  χ² のみに依存し、この非一意性に左右されない**ため、物理的結論は安定。
- **OSQP のスケーリング**: 自然単位定数 `c` が ~1e70 と巨大なため、内部で列スケーリング
  してから解く。これがないと OSQP/CLARABEL は誤った点を `optimal` として返す。
- **バンド探索の前提**: バンドが連結（上下に交差点 1 つずつ）であることを仮定する。

### バックエンドの違い（まとめ）

| | scipy | osqp |
|---|---|---|
| ソルバー | trust-constr | OSQP→CLARABEL→SCS（固定解は CLARABEL 優先） |
| フリー解 | 安定 | 高速 |
| 固定パラメータ解 | 安定 | CLARABEL（OSQP は大規模で未収束のため） |
| モンテカルロ並列 (`n_jobs`) | 有効 | 無効（逐次） |
| フラックス形状 | アルゴリズム依存の頂点 | HiGHS で区分定数の頂点を明示選択 |
