"""
Lipid Droplet Dispersion / Distribution-Shape Analysis
=======================================================
CIDEC KD vs NTC — look for signal in the *shape* of the size distribution
rather than its center.

Rationale
---------
CIDEC promotes LD fusion. The classic readouts (mean / median size, count)
showed no biologically meaningful difference. But fusion vs fission is
fundamentally a change in *dispersion*, not central tendency:

  - Fusion  -> a few droplets grow large  -> heavy right tail
              -> higher variance, higher CV, positive skew, fat tail (kurtosis),
                 larger fraction of "big" droplets.
  - Fission -> droplets become uniform     -> lower variance / CV / skew.

So CIDEC KD (loss of fusion) should, if any signal exists, show:
  CV down, skew down, kurtosis down, frac_large down, frac_small up.

This script re-uses the existing per-droplet measurements (no re-segmentation)
and computes spread / shape statistics at field and well level, plus pooled
two-sample tests that are explicitly sensitive to dispersion differences
(Levene for variance, KS / Mann-Whitney for distribution shift).
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
import matplotlib.pyplot as plt

# === Paths ===
SRC = Path("/home/QYJI/das/lypolysis/output/2026-06-02/lipid_droplet_measurements.csv")
OUT_DIR = Path("/home/QYJI/das/lypolysis/output/2026-06-08")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Size thresholds for tail fractions (um^2)
SMALL_TH = 20.0    # ~ d < 5 um
LARGE_TH = 100.0   # ~ d > 11 um
XLARGE_TH = 200.0  # ~ d > 16 um

COLORS = {"NTC": "#4477AA", "CIDEC_KD": "#CC6677"}


def dispersion_stats(a):
    """Spread / shape statistics for one group of droplet areas."""
    a = np.asarray(a, dtype=float)
    mean = a.mean()
    std = a.std(ddof=1)
    return pd.Series({
        "n": a.size,
        "mean": mean,
        "std": std,
        "cv": std / mean,                       # coefficient of variation
        "iqr": np.subtract(*np.percentile(a, [75, 25])),
        "qcd": (np.percentile(a, 75) - np.percentile(a, 25)) /
               (np.percentile(a, 75) + np.percentile(a, 25)),  # robust CV
        "skew": stats.skew(a),
        "kurtosis": stats.kurtosis(a),          # excess kurtosis (Fisher)
        "p90": np.percentile(a, 90),
        "p99": np.percentile(a, 99),
        "frac_small": np.mean(a < SMALL_TH),
        "frac_large": np.mean(a > LARGE_TH),
        "frac_xlarge": np.mean(a > XLARGE_TH),
    })


def fmt_p(p):
    return f"{p:.4g} " + ("***" if p < 0.001 else "**" if p < 0.01 else
                          "*" if p < 0.05 else "ns")


def main():
    print("=" * 70)
    print("Lipid Droplet Dispersion / Shape Analysis  (NTC vs CIDEC KD)")
    print("=" * 70)

    df = pd.read_csv(SRC, usecols=["area_um2", "condition", "well", "field"])
    print(f"Loaded {len(df):,} droplets "
          f"({df['condition'].value_counts().to_dict()})")

    # ---- Per-field dispersion (n = 54 vs 54) ----
    field_disp = (df.groupby(["condition", "well", "field"])["area_um2"]
                    .apply(dispersion_stats).unstack().reset_index())
    field_disp.to_csv(OUT_DIR / "field_dispersion.csv", index=False)

    # ---- Per-well dispersion (n = 6 vs 6, pools all droplets in a well) ----
    well_disp = (df.groupby(["condition", "well"])["area_um2"]
                   .apply(dispersion_stats).unstack().reset_index())
    well_disp.to_csv(OUT_DIR / "well_dispersion.csv", index=False)

    metrics = ["std", "cv", "qcd", "skew", "kurtosis", "p90", "p99",
               "frac_small", "frac_large", "frac_xlarge"]
    expected = {  # predicted direction if CIDEC KD -> fission (less fusion)
        "std": "KD<NTC", "cv": "KD<NTC", "qcd": "KD<NTC", "skew": "KD<NTC",
        "kurtosis": "KD<NTC", "p90": "KD<NTC", "p99": "KD<NTC",
        "frac_small": "KD>NTC", "frac_large": "KD<NTC", "frac_xlarge": "KD<NTC",
    }

    def compare(table, label, n_note):
        print("\n" + "-" * 70)
        print(f"{label}  ({n_note})")
        print("-" * 70)
        print(f"{'metric':<12}{'NTC':>12}{'CIDEC_KD':>12}{'delta%':>9}"
              f"{'p(t-test)':>14}  expect")
        rows = []
        ntc = table[table.condition == "NTC"]
        kd = table[table.condition == "CIDEC_KD"]
        for m in metrics:
            a, b = ntc[m].values, kd[m].values
            t, p = stats.ttest_ind(a, b, equal_var=False)
            delta = (b.mean() - a.mean()) / a.mean() * 100
            direction = "KD>NTC" if b.mean() > a.mean() else "KD<NTC"
            match = "OK" if direction == expected[m] else "opp"
            print(f"{m:<12}{a.mean():>12.4g}{b.mean():>12.4g}{delta:>+8.2f}%"
                  f"{fmt_p(p):>14}  {expected[m]} [{match}]")
            rows.append({"metric": m, "NTC_mean": a.mean(), "KD_mean": b.mean(),
                         "delta_pct": delta, "p_value": p,
                         "observed_dir": direction, "expected_dir": expected[m],
                         "match_expected": direction == expected[m]})
        return pd.DataFrame(rows)

    field_test = compare(field_disp, "PER-FIELD comparison (Welch t-test)",
                         "n = 54 vs 54")
    field_test.insert(0, "level", "per_field")
    well_test = compare(well_disp, "PER-WELL comparison (Welch t-test)",
                        "n = 6 vs 6  (proper biological replicate unit)")
    well_test.insert(0, "level", "per_well")

    pd.concat([field_test, well_test], ignore_index=True).to_csv(
        OUT_DIR / "dispersion_tests.csv", index=False)

    # ---- Pooled droplet-level tests (dispersion-sensitive) ----
    a = df[df.condition == "NTC"]["area_um2"].values
    b = df[df.condition == "CIDEC_KD"]["area_um2"].values
    print("\n" + "-" * 70)
    print("POOLED droplet-level tests (all droplets)")
    print("-" * 70)
    # Levene: tests equality of variance (the core dispersion question)
    lev_W, lev_p = stats.levene(a, b, center="median")
    print(f"  Levene (variance equality):  W={lev_W:.2f}  p={fmt_p(lev_p)}")
    # KS: any distributional difference
    ks_D, ks_p = stats.ks_2samp(a, b)
    print(f"  KS (distribution shape):     D={ks_D:.4f}  p={fmt_p(ks_p)}")
    # Mann-Whitney: stochastic dominance / shift
    mw_U, mw_p = stats.mannwhitneyu(a, b, alternative="two-sided")
    print(f"  Mann-Whitney U:              p={fmt_p(mw_p)}")
    print(f"  Variance:  NTC={a.var(ddof=1):.2f}  KD={b.var(ddof=1):.2f}  "
          f"ratio KD/NTC={b.var(ddof=1)/a.var(ddof=1):.4f}")

    # ---- Figures ----
    make_figures(df, field_disp)

    print("\n" + "=" * 70)
    n_match = field_test["match_expected"].sum() + well_test["match_expected"].sum()
    print(f"Metrics matching fission prediction: "
          f"{n_match}/{2*len(metrics)}")
    print(f"Outputs -> {OUT_DIR}")
    print("=" * 70)


def make_figures(df, field_disp):
    # ECDF of pooled areas (the most direct view of distribution shape)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle("Lipid Droplet Size DISPERSION: NTC vs CIDEC KD",
                 fontsize=14, fontweight="bold")

    ax = axes[0, 0]
    for cond in ["NTC", "CIDEC_KD"]:
        a = np.sort(df[df.condition == cond]["area_um2"].values)
        y = np.arange(1, a.size + 1) / a.size
        ax.plot(a, y, label=cond, color=COLORS[cond], lw=1.5)
    ax.set_xscale("log")
    ax.set_xlabel("Area (um^2, log)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Pooled ECDF (full distribution)")
    ax.legend()

    # Zoom on the large-droplet tail (where fusion signal lives)
    ax = axes[0, 1]
    for cond in ["NTC", "CIDEC_KD"]:
        a = np.sort(df[df.condition == cond]["area_um2"].values)
        y = 1 - np.arange(1, a.size + 1) / a.size  # survival
        ax.plot(a, y, label=cond, color=COLORS[cond], lw=1.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Area (um^2, log)")
    ax.set_ylabel("P(Area > x)  (log)")
    ax.set_title("Right-tail survival (fusion signal)")
    ax.legend()

    # Per-field dispersion box plots
    box_metrics = [("cv", "Coefficient of variation"),
                   ("skew", "Skewness"),
                   ("kurtosis", "Excess kurtosis"),
                   ("frac_large", "Fraction area > 100 um^2")]
    positions = [(0, 2), (1, 0), (1, 1), (1, 2)]
    for (m, title), pos in zip(box_metrics, positions):
        ax = axes[pos]
        data = [field_disp[field_disp.condition == c][m].values
                for c in ["NTC", "CIDEC_KD"]]
        bp = ax.boxplot(data, tick_labels=["NTC", "CIDEC KD"],
                        patch_artist=True)
        for patch, c in zip(bp["boxes"], [COLORS["NTC"], COLORS["CIDEC_KD"]]):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        t, p = stats.ttest_ind(data[0], data[1], equal_var=False)
        ax.set_title(f"{title}\nWelch p={p:.3g}")
        ax.set_ylabel(m)

    plt.tight_layout()
    fig.savefig(OUT_DIR / "dispersion_analysis.png", dpi=150,
                bbox_inches="tight")
    plt.close()
    print("  -> Saved: dispersion_analysis.png")


if __name__ == "__main__":
    main()
