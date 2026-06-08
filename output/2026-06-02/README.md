# KDE Plate Lipid Droplet Analysis

## 实验信息

| 项目 | 内容 |
|------|------|
| 数据来源 | Yang Huan (Sr Scientist I, Target Discovery) |
| 日期 | 2026-05-27 拍摄, 2026-06-02 分析 |
| 仪器 | PerkinElmer Sonata |
| 板型 | 96-well Corning 353072 |
| 物镜 | 10X Air, NA 0.3 |
| 通道 | Brightfield (单通道) |
| Binning | 2×2 |
| 像素尺寸 | 1.196 μm/pixel |
| 图像大小 | 1080×1080 px |
| 拍摄范围 | Row 2-7, Col 2-11 (60 wells) |
| 每孔视野 | 9 fields |
| Z 层数 | 5 planes |
| 总图片数 | 2700 tiff |

## 实验设计

- **Plate ID**: EEBKT1&2Batch10Plate2KDE-Lonza-NN0006G6JZ-10XAi
- **NTC (Negative Control)**: Col 3, Rows 2-7
- **CIDEC KD (Positive)**: Col 10, Rows 2-7
- **目的**: CIDEC 促进脂滴 fusion → KD 后预期脂滴变小变多 (fission phenotype)
- **背景**: 探索性活动，为未来 adipocyte-based imaging assay 开发做准备

## 分析方法

### Pipeline (v3)

1. **Best focal plane selection** — Brenner gradient 在 5 个 Z-plane 中选择最佳聚焦面
2. **Lipid droplet segmentation** — Local background subtraction + watershed:
   - Gaussian smooth (σ=1) → signal
   - Gaussian smooth (σ=15) → local background
   - Signal - Background > 0.015 → binary mask
   - Remove small objects (<6 px) + fill holes
   - **Distance-transform watershed** 拆分聚集的脂滴 cluster (min_distance=3 px)
3. **Quantification** — `regionprops_table` 提取 area, equivalent_diameter
4. **Size filtering** — 保留 6-800 px (8.6-1145 μm²) 的物体
5. **Statistical comparison** — Per-field aggregation + t-test

### 参数

```python
SEG_SIGMA_SMOOTH = 1.0      # 噪声平滑
SEG_SIGMA_BG = 15.0         # 背景估计尺度
SEG_THRESHOLD = 0.015       # 检测阈值
SEG_MIN_AREA = 6 px         # ~8.6 μm²
SEG_MAX_AREA = 800 px       # ~1145 μm² (d~34μm)
SEG_WATERSHED_DIST = 3 px   # watershed seed 最小间距
```

### 方法迭代记录

| 版本 | 方法 | 检出数/field | 速度 |
|------|------|-------------|------|
| v1 | CLAHE + local threshold | ~1000 | 慢 |
| v2 | Local background subtraction | ~9600 | 0.18s/img |
| **v3** | **+ watershed splitting** | **~10500** | 1.1s/img |

Watershed 将紧挨的小脂滴 cluster 拆分，物体检出数 +13.5%。

## 结果

### 1. 个体脂滴统计 (v3 watershed)

| 指标 | NTC (mean ± std) | CIDEC KD (mean ± std) | p-value | 方向 |
|------|-------------------|------------------------|---------|------|
| Droplets/field | 10523 ± 482 | 10567 ± 441 | 0.62 ns | 无差异 |
| Mean area (μm²) | 57.3 ± 2.1 | 57.6 ± 1.9 | 0.46 ns | 无差异 |
| Median area (μm²) | 32.1 ± 1.9 | 33.6 ± 2.1 | 0.0002 *** | KD ↑ 5% |
| Mean diameter (μm) | 7.46 ± 0.11 | 7.55 ± 0.14 | 0.0012 ** | KD ↑ 1% |

- 总测量脂滴数: 1,138,839
- NTC 54 fields, CIDEC KD 54 fields
- 效应量极小（中位面积差仅 ~1.5%）

### 2. 前景/背景面积分析

直接二值分割（不分割个体），比较总脂滴覆盖面积占比：

| 聚合层级 | NTC | CIDEC KD | p-value |
|---------|-----|----------|---------|
| Per-field (n=54) | 41.39% ± 1.04% | 40.99% ± 1.02% | 0.047 * |
| Per-well (n=6) | 41.39% ± 0.50% | 40.99% ± 0.77% | 0.35 ns |

以 well 为统计单位（更合理）时**无显著差异**。这排除了"小脂滴分割遗漏导致漏检差异"的可能——总前景面积已包含所有亮信号区域，两组仍一致。

### 3. 阈值/分割方法探索

针对"背景纹理被误判为前景"的问题测试了多种方案：

| 方法 | 前景占比 | 评价 |
|------|---------|------|
| DoG th=0.015 (当前) | 40.9% | 含部分背景纹理 |
| DoG th=0.030 | 35.3% | 更严格 |
| DoG th=0.020 + opening | 36.8% | 平衡 |
| Otsu adaptive | 11.6% | 太严格，丢暗脂滴 |
| Local Gaussian/Mean | 50-51% | **更差**，背景噪声被放大 |
| Sauvola / Niblack | 50-66% | **更差**，文档二值化假设不适用 |

**结论**：局部自适应阈值（Sauvola/Niblack/local mean）假设每个窗口内都有前景+背景两类，但明场图的真实背景是平坦无信号区，局部阈值会把噪声强行二分，反而误检更多。明场 10X 下脂滴与细胞质无干净强度分界，调阈值收益有限，此方向暂停。

### 综合解读

**三个独立量化维度（个体计数/大小、总前景面积、分割形态）一致显示 NTC 与 CIDEC KD 无显著差异。** 与经典 "CIDEC 促 fusion → KD 导致 fission" 模型不符。

可能原因：
1. **KD 效率不足** — 未有效降低 CIDEC 功能（需 qPCR/WB 确认）
2. **Brightfield 局限性** — 10X 明场无法捕捉 CIDEC 介导的脂滴形态变化；细胞质 texture 与脂滴难以区分
3. **缺乏荧光标记** — 无法特异性识别 LD

### 建议

- **改用荧光脂滴染色**（BODIPY/Nile Red）获得干净的前景背景分离，是推进此 assay 的关键
- 加入 Hoechst 做 cell count 归一化
- 确认 KD 效率（per-well qPCR 或 western）

### 后续分析：离散度 / 分布形状（2026-06-08）

考虑到 fusion/fission 改变的是分布**离散度**而非均值，进一步计算了 std / CV / 偏度 / 峰度 / 尾部比例（复用已有 114 万脂滴测量，未重新分割）。

- **苗头**：CIDEC KD 组脂滴尺寸略更均匀（std/CV/QCD 一致 ↓~3%，CV field 水平 p=0.003），方向符合 fission 预期。
- **但不可靠**：① pooled p≈1e-164 仅为 114 万样本量伪影（KS D=0.026，方差只差 6%）；② **以 well 为统计单位时全部不显著**；③ skew/frac_small 方向相反，内部不自洽。
- **结论**：再次印证瓶颈是明场成像本身，非统计方法。详见 [`../2026-06-08/README.md`](../2026-06-08/README.md)。

## 输出文件

| 文件 | 说明 |
|------|------|
| `lipid_droplet_measurements.csv` | 所有单个脂滴测量 (1.14M rows) |
| `field_summary.csv` | 每个视野的统计汇总 (108 rows) |
| `foreground_area_analysis.csv` | 前景面积占比逐视野数据 |
| `lipid_droplet_comparison.png` | 个体脂滴统计比较图 (6 panels) |
| `foreground_area_comparison.png` | 前景面积占比比较图 |
| `foreground_mask_visualization.png` | 前景/背景掩膜全视野可视化 |
| `segmentation_v3_zoomed.png` | v3 watershed 分割放大 QC |
| `watershed_comparison.png` | watershed 前后对比 |
| `threshold_comparison.png` | 全局阈值方法对比 |
| `adaptive_threshold_comparison.png` | 局部自适应阈值对比 |
| `segmentation_examples.png` | 分割质量示例 (4 wells) |
| `segmentation_v2_qc.png` | 分割参数调优 QC |
| `segmentation_debug.png` | 方法比较 (tophat/local thresh) |
| `blob_detection_test.png` | LoG blob detection 测试 |
| `preview/` | Huan 提供的代表性 PNG |

## 数据路径

```
源数据: /TDE_TV/shared_folder/OFGM/20260527-KDE plate/
脚本:  /home/QYJI/das/lypolysis/kde_lipid_droplet_analysis.py
输出:  /home/QYJI/das/lypolysis/output/2026-06-02/
```
