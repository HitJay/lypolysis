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

### Pipeline (v2)

1. **Best focal plane selection** — Brenner gradient 在 5 个 Z-plane 中选择最佳聚焦面
2. **Lipid droplet segmentation** — Local background subtraction:
   - Gaussian smooth (σ=1) → signal
   - Gaussian smooth (σ=15) → local background
   - Signal - Background > 0.015 → binary mask
   - Remove small objects (<6 px)
   - Fill holes
   - Connected component labeling
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
```

## 结果

### 统计摘要

| 指标 | NTC (mean ± std) | CIDEC KD (mean ± std) | p-value | 方向 |
|------|-------------------|------------------------|---------|------|
| Droplets/field | 9270 ± 467 | 8845 ± 724 | 0.0005 *** | KD ↓ 5% |
| Mean area (μm²) | 51.0 ± 3.8 | 54.8 ± 4.5 | <0.0001 *** | KD ↑ 7% |
| Median area (μm²) | 30.8 ± 2.6 | 32.4 ± 2.3 | 0.0014 ** | KD ↑ 5% |
| Mean diameter (μm) | 7.17 ± 0.27 | 7.40 ± 0.28 | <0.0001 *** | KD ↑ 3% |

- 总测量脂滴数: 978,200
- NTC 54 fields, CIDEC KD 54 fields

### 解读

**CIDEC KD 后脂滴略少(-5%)且略大(+7%)**，与经典 "CIDEC 促 fusion → KD 导致 fission" 模型方向不一致。

可能原因：
1. **Effect size 很小** — 5-7% 差异，虽统计显著但生物学意义待确认
2. **Well-to-well 变异** — CIDEC KD 组 R02-R04 和 R05-R07 表现不同（可能 KD 效率差异）
3. **Brightfield 局限性** — 10X 明场下细胞间隙/texture 可能被误识为小脂滴，掩盖真实差异
4. **缺乏荧光标记** — 无法区分 LD 和其他亮的细胞结构

### 建议

- 后续可用 BODIPY/Nile Red 荧光标记脂滴，提高检测特异性
- 考虑加入 Hoechst 做 cell count 归一化
- 如需确认 KD 效率，可加入 per-well qPCR 或 western 数据关联

## 输出文件

| 文件 | 说明 |
|------|------|
| `lipid_droplet_measurements.csv` | 所有单个脂滴测量 (978k rows) |
| `field_summary.csv` | 每个视野的统计汇总 (108 rows) |
| `lipid_droplet_comparison.png` | 统计比较图 (6 panels) |
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
