# 脂滴尺寸离散度 / 分布形状分析

> 2026-06-02 的均值/中位数分析显示 NTC 与 CIDEC KD 无差异。
> 本次（2026-06-08）改从**分布的离散度与形状**入手，复用已有的 114 万脂滴测量，
> 不重新分割。

## 动机

fusion ↔ fission 本质上改变的是**分布的离散度**，而非中心位置：

| 表型 | 预期分布变化 |
|------|-------------|
| **Fusion**（CIDEC 正常）| 少数脂滴变大 → 重右尾 → 方差↑、CV↑、偏度↑、峰度↑、大脂滴比例↑ |
| **Fission**（CIDEC KD）| 脂滴趋于均匀 → 方差↓、CV↓、偏度↓、大脂滴比例↓ |

因此若 CIDEC KD 真的导致 fission（失去 fusion），预期：
**CV↓、skew↓、kurtosis↓、frac_large↓、frac_small↑**。

均值/中位数会把这些差异抹平，所以专门计算离散度指标。

## 方法

- 数据源：`output/2026-06-02/lipid_droplet_measurements.csv`（1,138,839 脂滴）
- 脚本：`dispersion_analysis.py`
- 每个 field / well 计算：std, CV, QCD（稳健 CV）, IQR, 偏度, 峰度, p90, p99, 小/大/超大脂滴比例
- 检验：
  - 逐 field (n=54 vs 54) 与逐 well (n=6 vs 6, 正确生物学重复单位) 的 Welch t-test
  - Pooled 全脂滴：Levene（方差齐性）、KS（分布形状）、Mann-Whitney（位移）

阈值：small <20 μm²，large >100 μm²，xlarge >200 μm²

## 结果

### 逐 field 比较 (n = 54 vs 54)

| 参数 | NTC | CIDEC KD | Δ% | p (Welch) | 符合 fission? |
|------|-----|----------|-----|-----------|:---:|
| **std** | 92.99 | 90.10 | −3.1% | 0.038 * | ✅ |
| **CV** | 1.621 | 1.565 | −3.5% | **0.0031 \*\*** | ✅ |
| **QCD** | 0.482 | 0.475 | −1.5% | 0.0089 ** | ✅ |
| skew | 5.74 | 5.86 | +2.1% | 0.045 * | ❌ |
| kurtosis | 43.06 | 45.14 | +4.8% | 0.061 ns | ❌ |
| p90 | 106.3 | 106.2 | −0.06% | 0.94 ns | ✅ |
| p99 | 519.0 | 499.4 | −3.8% | 0.063 ns | ✅ |
| frac_small | 0.245 | 0.224 | −8.5% | 4.4e-5 *** | ❌ |
| frac_large | 0.109 | 0.110 | +1.0% | 0.38 ns | ❌ |
| frac_xlarge | 0.0433 | 0.0408 | −5.8% | 0.018 * | ✅ |

### 逐 well 比较 (n = 6 vs 6，正确统计单位)

**所有指标均不显著**（CV p=0.064，QCD p=0.057，std p=0.069，其余更高）。
field 水平的显著性是**伪重复（pseudo-replication）**造成的。

### Pooled 全脂滴检验 (114 万)

| 检验 | 统计量 | p | 效应量 |
|------|--------|---|--------|
| Levene（方差齐性） | W=14.17 | 1.7e-4 *** | variance ratio KD/NTC = **0.938**（仅差 6%）|
| KS（分布形状） | D=0.0257 | ~1e-164 | **D 极小** |
| Mann-Whitney | — | ~1e-163 | — |

> ECDF 与右尾 survival 曲线几乎完全重合（见 `dispersion_analysis.png`）。
> 极小的 p 值纯粹来自 114 万的样本量，并非实质效应。

## 综合解读

**离散度分析比均值分析多挖出一丝苗头，但不足以作为可靠表型信号。**

- ✅ **方向自洽的部分指向 fission**：CIDEC KD 组脂滴尺寸更均匀
  （std / CV / QCD 一致下降 ~3%，CV 在 field 水平 p=0.003），
  符合"CIDEC KD → 失去 fusion → 更均匀"的预期。
- ⚠️ **三点致命折扣**：
  1. **样本量伪影** — pooled p≈1e-164 仅因 N 巨大；真实效应量极小（KS D=0.026，方差只差 6%）。
  2. **正确统计单位下消失** — 以 well 为重复 (n=6) 时所有指标不显著。
  3. **内部不自洽** — std/CV↓ 指向 fission，但 skew↑、frac_small↓ 指向相反；20 个比较 12 个符合预测 ≈ 随机水平。

**结论**：进一步印证 2026-06-02 的判断——瓶颈是明场成像本身无法特异识别脂滴，
而非统计方法不够细。真正的 next step 仍是**荧光脂滴染色（BODIPY/Nile Red）**，
届时这套离散度指标会给出干净得多的读数。

## 输出文件

| 文件 | 说明 |
|------|------|
| `field_dispersion.csv` | 逐 field 离散度指标 (108 rows) |
| `well_dispersion.csv` | 逐 well 离散度指标 (12 rows) |
| `dispersion_tests.csv` | 所有比较的检验汇总（field + well 水平）|
| `dispersion_analysis.png` | ECDF + 右尾 survival + CV/skew/kurtosis/frac_large box plots |

## 数据路径

```
脚本:  /home/QYJI/das/lypolysis/dispersion_analysis.py
输入:  /home/QYJI/das/lypolysis/output/2026-06-02/lipid_droplet_measurements.csv
输出:  /home/QYJI/das/lypolysis/output/2026-06-08/
```
