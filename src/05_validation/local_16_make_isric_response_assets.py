#!/usr/bin/env python
"""Create manuscript/response-ready tables and figures for ISRIC validation."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\Codex\260724 小论文")
IN_DIR = ROOT / "revise" / "isric_salinity_validation"
OUT_DIR = IN_DIR


def fmt(x: float) -> str:
    if pd.isna(x):
        return ""
    return f"{x:.3f}"


def main() -> None:
    corr = pd.read_csv(IN_DIR / "isric_spearman_correlations.csv")
    binary = pd.read_csv(IN_DIR / "isric_saline_vs_nonsaline_harsei_minus_rsei.csv")
    gradient = pd.read_csv(IN_DIR / "isric_harsei_vs_rsei_salinity_gradient_test.csv")

    annual_corr = corr[corr["year"].astype(str) != "pooled"].copy()
    annual_corr["year"] = annual_corr["year"].astype(int)
    binary["year"] = binary["year"].astype(int)
    gradient["year"] = gradient["year"].astype(int)
    corr_wide = annual_corr.pivot(index="year", columns="variable", values="spearman_rho_with_salinity")
    corr_n = annual_corr.groupby("year")["n"].max()

    delta = binary[binary["group"] == "saline_minus_non_saline"].copy().set_index("year")
    gradient = gradient.set_index("year")

    rows = []
    for year in sorted(corr_wide.index.astype(int)):
        rows.append(
            {
                "year": year,
                "n_pixels": int(corr_n.loc[year]),
                "rho_SRSI_salinity": corr_wide.loc[year, "SRSI"],
                "rho_RSEI_salinity": corr_wide.loc[year, "RSEI"],
                "rho_HARSEI_salinity": corr_wide.loc[year, "HARSEI"],
                "abs_rho_HARSEI_minus_abs_rho_RSEI": abs(corr_wide.loc[year, "HARSEI"]) - abs(corr_wide.loc[year, "RSEI"]),
                "delta_median_SRSI_saline_minus_non_saline": delta.loc[year, "SRSI_median"],
                "delta_median_RSEI_saline_minus_non_saline": delta.loc[year, "RSEI_median"],
                "delta_median_HARSEI_saline_minus_non_saline": delta.loc[year, "HARSEI_median"],
                "delta_median_HARSEI_minus_RSEI_saline_minus_non_saline": delta.loc[year, "HARSEI_minus_RSEI_median"],
                "HARSEI_clearer_decline_than_RSEI": gradient.loc[year, "HARSEI_clearer_decline_than_RSEI"],
            }
        )

    summary = pd.DataFrame(rows)
    csv_path = OUT_DIR / "table_isric_validation_summary_for_response.csv"
    summary.to_csv(csv_path, index=False, encoding="utf-8-sig")

    md = summary.copy()
    for col in md.columns:
        if col not in {"year", "n_pixels", "HARSEI_clearer_decline_than_RSEI"}:
            md[col] = md[col].map(fmt)
    md_path = OUT_DIR / "table_isric_validation_summary_for_response.md"
    md_path.write_text(md.to_markdown(index=False), encoding="utf-8")

    try:
        import matplotlib.pyplot as plt

        years = summary["year"].astype(str)
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), dpi=300)

        axes[0].axhline(0, color="0.25", lw=0.8)
        axes[0].plot(years, summary["rho_SRSI_salinity"], marker="o", label="SRSI")
        axes[0].plot(years, summary["rho_RSEI_salinity"], marker="s", label="RSEI")
        axes[0].plot(years, summary["rho_HARSEI_salinity"], marker="^", label="HARSEI")
        axes[0].set_ylabel("Spearman rho with ISRIC salinity")
        axes[0].set_xlabel("Year")
        axes[0].legend(frameon=False)
        axes[0].grid(True, color="0.88", lw=0.8)

        width = 0.25
        x = range(len(years))
        axes[1].axhline(0, color="0.25", lw=0.8)
        axes[1].bar([i - width for i in x], summary["delta_median_SRSI_saline_minus_non_saline"], width=width, label="SRSI")
        axes[1].bar(list(x), summary["delta_median_RSEI_saline_minus_non_saline"], width=width, label="RSEI")
        axes[1].bar([i + width for i in x], summary["delta_median_HARSEI_saline_minus_non_saline"], width=width, label="HARSEI")
        axes[1].set_xticks(list(x))
        axes[1].set_xticklabels(years)
        axes[1].set_ylabel("Median difference: saline - non-saline")
        axes[1].set_xlabel("Year")
        axes[1].legend(frameon=False)
        axes[1].grid(True, axis="y", color="0.88", lw=0.8)

        fig.tight_layout()
        fig.savefig(OUT_DIR / "fig_isric_validation_summary_for_response.png")
        plt.close(fig)
    except Exception as exc:  # pragma: no cover
        print(f"Figure generation skipped: {exc}")

    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
