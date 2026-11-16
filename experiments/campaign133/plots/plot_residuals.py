#!/usr/bin/env python3
"""Plot raw history columns; run from any directory with Matplotlib installed."""
import csv
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SERIES = (
    ("S01-warm-unmodified64", "Warm unmodified (64 evaluations)"),
    ("S01-warm-FD64", "Warm FD (64 evaluations)"),
)


def history(name):
    path = ROOT / "evidence" / name / "history.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    expected = 64
    if len(rows) != expected:
        raise ValueError(f"{path}: expected {expected} rows, got {len(rows)}")
    for row in rows:
        for column in ("r_K", "r_Σ"):
            value = float(row[column])
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{path}: {column} must be finite and positive")
    return rows


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: plot_residuals.py OUTPUT.png")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True, sharey=True)
    for name, label in SERIES:
        rows = history(name)
        for ax, column in zip(axes, ("r_K", "r_Σ")):
            ax.semilogy(range(1, len(rows) + 1), [float(row[column]) for row in rows], label=label)
    for ax, title in zip(axes, (r"$r_K$", r"$r_\Sigma$")):
        ax.axhline(1e-8, color="black", linestyle="--", label="strict 1e-8")
        ax.set(xlim=(1, 64), ylim=(5e-9, 2), xlabel="Map evaluation within each run", ylabel="Raw residual", title=title)
        ax.grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize="small")
    fig.text(0.5, 0.01, "Both unresolved | Fixed U, no Poisson | Map counts reset\n"
             "Unmodified: evaluation 1 repeats source boundary; FD repair at initialization only",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    fig.savefig(sys.argv[1], dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
