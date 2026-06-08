"""
KDE Plate Lipid Droplet Analysis (v2)
======================================
CIDEC KD vs NTC: Assess lipid droplet fusion vs fission phenotype.
CIDEC promotes LD fusion → KD should yield smaller, more numerous droplets (fission).

Data: 20260527-KDE plate, Brightfield 10X, 96-well (rows 2-7, cols 2-11)
- NTC: Col 3 (negative control)
- CIDEC KD: Col 10 (positive/target)

Pipeline:
1. Best focal plane selection (max Brenner gradient across Z-stack)
2. Lipid droplet segmentation (local background subtraction + watershed)
3. Morphometric quantification
4. Statistical comparison NTC vs CIDEC KD
"""

import os
import numpy as np
from PIL import Image
from pathlib import Path
from scipy import ndimage
from skimage import filters, morphology, measure, exposure, segmentation
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

# === Paths ===
IMG_DIR = Path("/TDE_TV/shared_folder/OFGM/20260527-KDE plate/"
               "20260527-EEBKT1&2Batch10Plate2KDE-Lonza-NN0006G6JZ-10XAi__"
               "2026-05-27T14_51_09-Measurement 1/Images")
OUT_DIR = Path("/home/QYJI/das/lypolysis/output/2026-06-02")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# === Plate layout ===
# From tU/ folder: NTC = col3, CIDEC KD = col10
# Rows with data: 2-7
NTC_WELLS = [(r, 3) for r in range(2, 8)]
CIDEC_WELLS = [(r, 10) for r in range(2, 8)]
N_FIELDS = 9
N_PLANES = 5
PIXEL_SIZE_UM = 1.196  # μm/pixel

# Segmentation parameters (tuned for 10X brightfield adipocyte LDs)
SEG_SIGMA_SMOOTH = 1.0      # Gaussian smoothing sigma
SEG_SIGMA_BG = 15.0         # Background estimation sigma
SEG_THRESHOLD = 0.015       # Intensity above local background
SEG_MIN_AREA = 6            # Min droplet area (px) ~ 8.6 μm²
SEG_MAX_AREA = 800          # Max droplet area (px) ~ 1145 μm² (d~34μm)
SEG_WATERSHED_DIST = 3      # Min distance between watershed seeds


def load_image(row, col, field, plane):
    """Load a single TIFF image."""
    fname = f"r{row:02d}c{col:02d}f{field:02d}p{plane:02d}-ch1sk1fk1fl1.tiff"
    fpath = IMG_DIR / fname
    if not fpath.exists():
        return None
    img = Image.open(fpath)
    return np.array(img).astype(np.float32)


def brenner_focus(img):
    """Brenner gradient focus metric - higher = more in focus."""
    dx = img[2:, :] - img[:-2, :]
    return np.mean(dx ** 2)


def select_best_plane(row, col, field):
    """Select the best focal plane by Brenner gradient."""
    scores = []
    for p in range(1, N_PLANES + 1):
        img = load_image(row, col, field, p)
        if img is None:
            scores.append(0)
        else:
            scores.append(brenner_focus(img))
    best_p = np.argmax(scores) + 1
    return best_p, scores


def segment_lipid_droplets(img):
    """
    Segment lipid droplets from brightfield image.
    
    Strategy: Local background subtraction + watershed splitting.
    In brightfield 10X, LDs appear as bright refractile objects above
    local background. We subtract a large-sigma Gaussian (background)
    from a small-sigma Gaussian (signal) to isolate small bright objects.
    Then use distance-transform watershed to split merged clusters.
    """
    # Normalize to 0-1
    img_norm = (img - img.min()) / (img.max() - img.min() + 1e-8)
    
    # Smooth to reduce noise
    img_smooth = filters.gaussian(img_norm, sigma=SEG_SIGMA_SMOOTH)
    
    # Estimate local background with large Gaussian
    local_bg = filters.gaussian(img_norm, sigma=SEG_SIGMA_BG)
    
    # Bright objects = above local background
    diff = img_smooth - local_bg
    binary = diff > SEG_THRESHOLD
    
    # Remove small noise
    binary = morphology.remove_small_objects(binary, min_size=SEG_MIN_AREA)
    
    # Fill holes in objects
    binary = ndimage.binary_fill_holes(binary)
    
    # Watershed splitting to separate touching droplets
    # Distance transform finds center of each blob
    distance = ndimage.distance_transform_edt(binary)
    # Smooth distance map slightly to avoid over-fragmentation
    distance = filters.gaussian(distance, sigma=0.8)
    # Find local maxima as watershed seeds
    from skimage.feature import peak_local_max
    coords = peak_local_max(distance, min_distance=SEG_WATERSHED_DIST,
                            labels=binary.astype(int))
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers = measure.label(mask)
    # Apply watershed on inverted distance
    labels = segmentation.watershed(-distance, markers, mask=binary)
    
    return labels


def quantify_droplets(labels):
    """Extract morphometric features from segmented droplets (vectorized)."""
    if labels.max() == 0:
        return pd.DataFrame()
    
    # Use regionprops_table for vectorized extraction
    # Only compute cheap properties for speed (~10k objects/image)
    table = measure.regionprops_table(labels, properties=(
        'label', 'area', 'equivalent_diameter'
    ))
    df = pd.DataFrame(table)
    
    # Filter by area range
    df = df[(df['area'] >= SEG_MIN_AREA) & (df['area'] <= SEG_MAX_AREA)].copy()
    
    if len(df) == 0:
        return pd.DataFrame()
    
    # Convert units
    df['area_px'] = df['area']
    df['area_um2'] = df['area'] * PIXEL_SIZE_UM ** 2
    df['equivalent_diameter_um'] = df['equivalent_diameter'] * PIXEL_SIZE_UM
    
    df = df[['area_px', 'area_um2', 'equivalent_diameter_um']].reset_index(drop=True)
    return df


def process_well(row, col, condition):
    """Process all fields in a well."""
    well_results = []
    
    for field in range(1, N_FIELDS + 1):
        # Select best focal plane
        best_p, _ = select_best_plane(row, col, field)
        img = load_image(row, col, field, best_p)
        
        if img is None:
            continue
        
        # Segment lipid droplets
        labels = segment_lipid_droplets(img)
        
        # Quantify
        df = quantify_droplets(labels)
        if len(df) > 0:
            df['row'] = row
            df['col'] = col
            df['field'] = field
            df['best_plane'] = best_p
            df['condition'] = condition
            df['well'] = f"R{row:02d}C{col:02d}"
            df['n_droplets'] = len(df)
            well_results.append(df)
    
    if well_results:
        return pd.concat(well_results, ignore_index=True)
    return pd.DataFrame()


def main():
    print("=" * 60)
    print("KDE Plate Lipid Droplet Analysis")
    print("CIDEC KD vs NTC — Fusion or Fission?")
    print("=" * 60)
    
    # === Step 1: Process all wells ===
    all_results = []
    
    print("\n[1/4] Processing NTC wells (Col 3)...")
    for row, col in NTC_WELLS:
        print(f"  Well R{row:02d}C{col:02d}...", end=" ")
        df = process_well(row, col, "NTC")
        if len(df) > 0:
            all_results.append(df)
            print(f"{df['n_droplets'].iloc[0]} droplets/field avg")
        else:
            print("no data")
    
    print("\n[2/4] Processing CIDEC KD wells (Col 10)...")
    for row, col in CIDEC_WELLS:
        print(f"  Well R{row:02d}C{col:02d}...", end=" ")
        df = process_well(row, col, "CIDEC_KD")
        if len(df) > 0:
            all_results.append(df)
            print(f"{df['n_droplets'].iloc[0]} droplets/field avg")
        else:
            print("no data")
    
    if not all_results:
        print("ERROR: No data processed!")
        return
    
    results = pd.concat(all_results, ignore_index=True)
    results.to_csv(OUT_DIR / "lipid_droplet_measurements.csv", index=False)
    print(f"\nTotal droplets measured: {len(results)}")
    
    # === Step 2: Per-field summary ===
    print("\n[3/4] Computing per-field statistics...")
    field_summary = results.groupby(['condition', 'well', 'field']).agg(
        n_droplets=('area_um2', 'count'),
        mean_area_um2=('area_um2', 'mean'),
        median_area_um2=('area_um2', 'median'),
        total_area_um2=('area_um2', 'sum'),
        mean_diameter_um=('equivalent_diameter_um', 'mean'),
        median_diameter_um=('equivalent_diameter_um', 'median'),
    ).reset_index()
    
    field_summary.to_csv(OUT_DIR / "field_summary.csv", index=False)
    
    # === Step 3: Statistical comparison ===
    print("\n[4/4] Statistical comparison...")
    ntc_fields = field_summary[field_summary['condition'] == 'NTC']
    cidec_fields = field_summary[field_summary['condition'] == 'CIDEC_KD']
    
    from scipy import stats
    
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    
    metrics = ['n_droplets', 'mean_area_um2', 'median_area_um2', 'mean_diameter_um']
    metric_labels = ['# Droplets/field', 'Mean area (μm²)', 'Median area (μm²)', 'Mean diameter (μm)']
    
    for metric, label in zip(metrics, metric_labels):
        ntc_vals = ntc_fields[metric].values
        cidec_vals = cidec_fields[metric].values
        t_stat, p_val = stats.ttest_ind(ntc_vals, cidec_vals)
        
        print(f"\n{label}:")
        print(f"  NTC:      {np.mean(ntc_vals):.2f} ± {np.std(ntc_vals):.2f}")
        print(f"  CIDEC KD: {np.mean(cidec_vals):.2f} ± {np.std(cidec_vals):.2f}")
        print(f"  p-value:  {p_val:.4f} {'***' if p_val<0.001 else '**' if p_val<0.01 else '*' if p_val<0.05 else 'ns'}")
    
    # Interpretation
    ntc_n = ntc_fields['n_droplets'].mean()
    cidec_n = cidec_fields['n_droplets'].mean()
    ntc_size = ntc_fields['mean_area_um2'].mean()
    cidec_size = cidec_fields['mean_area_um2'].mean()
    
    print("\n" + "=" * 60)
    print("INTERPRETATION:")
    if cidec_n > ntc_n and cidec_size < ntc_size:
        print("  → CIDEC KD: MORE droplets, SMALLER size")
        print("  → Phenotype: FISSION (consistent with CIDEC fusion role)")
    elif cidec_n < ntc_n and cidec_size > ntc_size:
        print("  → CIDEC KD: FEWER droplets, LARGER size")
        print("  → Phenotype: FUSION (unexpected for CIDEC KD)")
    else:
        print(f"  → CIDEC KD: n={cidec_n:.1f} vs NTC n={ntc_n:.1f}")
        print(f"  → CIDEC KD: size={cidec_size:.1f} vs NTC size={ntc_size:.1f}")
        print("  → Mixed phenotype, needs further investigation")
    print("=" * 60)
    
    # === Step 4: Generate figures ===
    generate_figures(results, field_summary)
    
    print(f"\nAll outputs saved to: {OUT_DIR}")


def generate_figures(results, field_summary):
    """Generate summary visualization figures."""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("CIDEC KD vs NTC: Lipid Droplet Analysis\n(KDE Plate 20260527, Brightfield 10X)", 
                 fontsize=14, fontweight='bold')
    
    colors = {'NTC': '#4477AA', 'CIDEC_KD': '#CC6677'}
    
    # 1. Droplet count per field
    ax = axes[0, 0]
    for cond in ['NTC', 'CIDEC_KD']:
        data = field_summary[field_summary['condition'] == cond]['n_droplets']
        ax.hist(data, bins=15, alpha=0.6, label=cond, color=colors[cond])
    ax.set_xlabel('# Droplets per field')
    ax.set_ylabel('Count')
    ax.set_title('Droplet Count Distribution')
    ax.legend()
    
    # 2. Droplet area distribution
    ax = axes[0, 1]
    for cond in ['NTC', 'CIDEC_KD']:
        data = results[results['condition'] == cond]['area_um2']
        ax.hist(data, bins=50, alpha=0.6, label=cond, color=colors[cond], density=True)
    ax.set_xlabel('Area (μm²)')
    ax.set_ylabel('Density')
    ax.set_title('Individual Droplet Area')
    ax.legend()
    ax.set_xlim(0, np.percentile(results['area_um2'], 95))
    
    # 3. Droplet diameter distribution
    ax = axes[0, 2]
    for cond in ['NTC', 'CIDEC_KD']:
        data = results[results['condition'] == cond]['equivalent_diameter_um']
        ax.hist(data, bins=50, alpha=0.6, label=cond, color=colors[cond], density=True)
    ax.set_xlabel('Equivalent Diameter (μm)')
    ax.set_ylabel('Density')
    ax.set_title('Individual Droplet Diameter')
    ax.legend()
    ax.set_xlim(0, np.percentile(results['equivalent_diameter_um'], 95))
    
    # 4. Box plot - droplets per field
    ax = axes[1, 0]
    bp_data = [field_summary[field_summary['condition'] == c]['n_droplets'].values 
               for c in ['NTC', 'CIDEC_KD']]
    bp = ax.boxplot(bp_data, labels=['NTC', 'CIDEC KD'], patch_artist=True)
    for patch, color in zip(bp['boxes'], [colors['NTC'], colors['CIDEC_KD']]):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_ylabel('# Droplets per field')
    ax.set_title('Droplet Count Comparison')
    
    # 5. Box plot - mean area per field
    ax = axes[1, 1]
    bp_data = [field_summary[field_summary['condition'] == c]['mean_area_um2'].values 
               for c in ['NTC', 'CIDEC_KD']]
    bp = ax.boxplot(bp_data, labels=['NTC', 'CIDEC KD'], patch_artist=True)
    for patch, color in zip(bp['boxes'], [colors['NTC'], colors['CIDEC_KD']]):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_ylabel('Mean area (μm²)')
    ax.set_title('Mean Droplet Area Comparison')
    
    # 6. Scatter: number vs size (each dot = one field)
    ax = axes[1, 2]
    for cond in ['NTC', 'CIDEC_KD']:
        df = field_summary[field_summary['condition'] == cond]
        ax.scatter(df['n_droplets'], df['mean_area_um2'], 
                  alpha=0.6, label=cond, color=colors[cond], s=40)
    ax.set_xlabel('# Droplets per field')
    ax.set_ylabel('Mean area (μm²)')
    ax.set_title('Number vs Size (per field)')
    ax.legend()
    
    plt.tight_layout()
    fig.savefig(OUT_DIR / "lipid_droplet_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  → Saved: lipid_droplet_comparison.png")
    
    # === Segmentation example figure ===
    generate_example_segmentation()


def generate_example_segmentation():
    """Show example segmentation for QC."""
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig.suptitle("Segmentation Examples (Best focal plane)", fontsize=12)
    
    examples = [
        (3, 3, 5, "NTC R03C03 F5"),   # NTC example
        (4, 3, 5, "NTC R04C03 F5"),
        (3, 10, 5, "CIDEC R03C10 F5"),  # CIDEC example
        (4, 10, 5, "CIDEC R04C10 F5"),
    ]
    
    for i, (row, col, field, title) in enumerate(examples):
        best_p, _ = select_best_plane(row, col, field)
        img = load_image(row, col, field, best_p)
        if img is None:
            continue
        
        labels = segment_lipid_droplets(img)
        
        # Raw image
        ax = axes[0, i]
        ax.imshow(img, cmap='gray')
        ax.set_title(f"{title}\nPlane {best_p}")
        ax.axis('off')
        
        # Overlay
        ax = axes[1, i]
        img_norm = (img - img.min()) / (img.max() - img.min())
        ax.imshow(img_norm, cmap='gray')
        # Overlay contours
        contours = measure.find_contours(labels > 0, 0.5)
        for contour in contours:
            ax.plot(contour[:, 1], contour[:, 0], 'r-', linewidth=0.5)
        n_obj = labels.max()
        ax.set_title(f"Segmented: {n_obj} objects")
        ax.axis('off')
    
    plt.tight_layout()
    fig.savefig(OUT_DIR / "segmentation_examples.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("  → Saved: segmentation_examples.png")


if __name__ == "__main__":
    main()
