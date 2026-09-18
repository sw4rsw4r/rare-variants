"""Build overview figures from deposited summary tables."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "overview"
OUT.mkdir(parents=True, exist_ok=True)

COLS = {
    "Unadjusted logistic": "#D55E00",
    "PC-adjusted logistic": "#0072B2",
    "CLR 1:6": "#E69F00",
    "CLR 1:4": "#CC79A7",
    "CLR 1:1": "#009E73",
}
P_MAP = {
    "Unadjusted logistic": "P_PLAIN",
    "PC-adjusted logistic": "P_PC_ADJUSTED",
    "CLR 1:6": "P_1TO6",
    "CLR 1:4": "P_1TO4",
    "CLR 1:1": "P_1TO1",
}


def style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def qq_from_p(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    p = np.sort(p[np.isfinite(p) & (p > 0) & (p <= 1)])
    n = len(p)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p)
    return expected, observed


def plot_qq(path: Path, title: str, out: Path) -> None:
    d = pd.read_csv(path, sep="\t")
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    ax.plot([0, 8], [0, 8], ls="--", color="0.45", lw=1)
    for method, col in P_MAP.items():
        exp, obs = qq_from_p(d[col].to_numpy(dtype=float))
        ax.plot(exp, obs, color=COLS[method], lw=1.2, label=method)
        ax.scatter(exp, obs, s=12, color=COLS[method], alpha=0.7, linewidths=0)
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 8)
    ax.set_xlabel(r"Expected $-\log_{10}(P)$")
    ax.set_ylabel(r"Observed $-\log_{10}(P)$")
    ax.set_title(title)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_simulation(path: Path, out: Path) -> None:
    d = pd.read_csv(path)
    keep = d["method"].isin(["Oracle GWAS", "PC-adjusted GWAS", "Matching + CLR"])
    d = d.loc[keep].copy()
    d["ratio_label"] = pd.Categorical(d["ratio_label"], ["1:6", "1:4", "1:1"], ordered=True)
    method_col = {
        "Oracle GWAS": "#6A3D9A",
        "PC-adjusted GWAS": "#0072B2",
        "Matching + CLR": "#009E73",
    }
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), sharex=True)
    panels = [
        ("Negative control", "Type I error", 0.05, axes[0]),
        ("Positive control", "Power", None, axes[1]),
    ]
    for scenario, ylab, href, ax in panels:
        sub = d[d["scenario"] == scenario]
        for method, color in method_col.items():
            m = sub[sub["method"] == method].sort_values("ratio_label")
            ax.errorbar(
                m["ratio_label"].astype(str),
                m["rejection_rate"],
                yerr=[
                    m["rejection_rate"] - m["ci_lower"],
                    m["ci_upper"] - m["rejection_rate"],
                ],
                color=color,
                marker="o",
                lw=1.4,
                capsize=3,
                label=method,
            )
        if href is not None:
            ax.axhline(href, ls="--", color="0.45", lw=1)
        ax.set_ylim(0, 1)
        ax.set_ylabel(ylab)
        ax.set_xlabel("Matching ratio")
        ax.set_title(scenario)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_pc(path: Path, out: Path) -> None:
    d = pd.read_csv(path, sep="\t")
    d["pc_num"] = d["VARIABLE"].str.replace("PC", "", regex=False).astype(int)
    ratios = ["1:6", "1:4", "1:1"]
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 4.2), sharey=True)
    for ax, ratio in zip(axes, ratios):
        sub = d[d["Ratio"] == ratio].sort_values("pc_num", ascending=False)
        y = np.arange(len(sub))
        ax.axvspan(-0.1, 0.1, color="#EAF4F1", zorder=0)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.scatter(sub["PRE_SMD"], y, s=28, color="#D55E00", label="Before", zorder=2)
        ax.scatter(sub["POST_SMD"], y, s=28, color="#009E73", label="After", zorder=3)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["VARIABLE"])
        ax.set_title(ratio)
        ax.set_xlabel("SMD")
        ax.set_xlim(-0.18, 0.18)
    axes[0].legend(frameon=False, loc="lower right")
    fig.suptitle("PC balance before and after matching", y=1.02)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_centres(path: Path, out: Path) -> None:
    d = pd.read_csv(path, sep="\t")
    d = d[d["Ratio"] == "1:6"].copy()
    d["abs_pre"] = d["PRE_SMD"].abs()
    d = d.sort_values("abs_pre", ascending=True)
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(6.4, 6.6))
    ax.axvspan(-0.1, 0.1, color="#EAF4F1", zorder=0)
    ax.axvline(0, color="0.5", lw=0.8)
    ax.scatter(d["PRE_SMD"], y, s=22, color="#D55E00", label="Before matching")
    ax.scatter(d["POST_SMD"], y, s=22, color="#009E73", label="After matching")
    ax.set_yticks(y)
    ax.set_yticklabels(d["CENTRE"])
    ax.set_xlabel("SMD")
    ax.set_title("Assessment centres at 1:6")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_rs4988235(path: Path, out: Path) -> None:
    d = pd.read_csv(path, sep="\t")
    order = [
        "Plain logistic",
        "PC-adjusted logistic",
        "CLR 1:6",
        "CLR 1:4",
        "CLR 1:1",
    ]
    d["Method"] = pd.Categorical(d["Method"], order, ordered=True)
    d = d.sort_values("Method")
    y = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.axvline(1, ls="--", color="0.45", lw=1)
    ax.errorbar(
        d["OR"],
        y,
        xerr=[d["OR"] - d["OR_lower"], d["OR_upper"] - d["OR"]],
        fmt="o",
        color="#009E73",
        capsize=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(d["Method"])
    ax.set_xlabel("Odds ratio per A allele")
    ax.set_title("rs4988235-A and lactose intolerance")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def main() -> None:
    style()
    real = ROOT / "summary" / "real"
    plot_qq(real / "lct_tv_comparison.tsv", "LCT region and TV watching", OUT / "qq_lct_tv.png")
    plot_qq(real / "chr2_tv_comparison.tsv", "Chromosome 2 rare variants and TV watching", OUT / "qq_chr2_tv.png")
    plot_simulation(ROOT / "summary" / "simulation" / "summary_metrics.csv", OUT / "simulation.png")
    plot_pc(real / "pc1_pc10_balance.tsv", OUT / "pc_smd.png")
    plot_centres(real / "recruitment_centre_balance_detailed.tsv", OUT / "centres.png")
    plot_rs4988235(real / "rs4988235_positive_control.tsv", OUT / "rs4988235.png")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
