#!/usr/bin/env python3
"""
plot_contig_lengths.py
Reads a FASTA assembly file and creates a contig length distribution plot.
Usage: python plot_contig_lengths.py <assembly.fasta> <output.pdf>
"""

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO


def main():
    if len(sys.argv) != 3:
        print("Usage: python plot_contig_lengths.py <assembly.fasta> <output.pdf>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = sys.argv[2]

    # ---- Read contig lengths from the FASTA ----
    lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths.append(len(record.seq))

    lengths.sort(reverse=True)

    # ---- Compute basic stats ----
    total = sum(lengths)
    n_contigs = len(lengths)
    largest = max(lengths)

    # Calculate N50
    running_sum = 0
    n50 = 0
    for l in lengths:
        running_sum += l
        if running_sum >= total / 2:
            n50 = l
            break

    # ---- Plot ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: histogram of contig lengths
    axes[0].hist(lengths, bins=50, color="steelblue", edgecolor="black")
    axes[0].set_xlabel("Contig Length (bp)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Contig Length Distribution")
    axes[0].axvline(n50, color="red", linestyle="--", label=f"N50 = {n50:,} bp")
    axes[0].legend()

    # Right: cumulative length (Nx) plot
    cumulative = []
    running = 0
    for l in lengths:
        running += l
        cumulative.append(running)

    axes[1].plot(range(1, n_contigs + 1), cumulative, color="steelblue")
    axes[1].set_xlabel("Contig Index (sorted by length)")
    axes[1].set_ylabel("Cumulative Length (bp)")
    axes[1].set_title("Cumulative Assembly Length")
    axes[1].axhline(total, color="gray", linestyle=":", label=f"Total = {total:,} bp")
    axes[1].legend()

    # Add summary text
    summary = (
        f"Contigs: {n_contigs}\n"
        f"Total: {total:,} bp\n"
        f"Largest: {largest:,} bp\n"
        f"N50: {n50:,} bp"
    )
    fig.text(0.02, 0.02, summary, fontsize=9, family="monospace",
             bbox=dict(boxstyle="round", facecolor="lightyellow"))

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"Plot saved to {output_file}")


if __name__ == "__main__":
    main()
