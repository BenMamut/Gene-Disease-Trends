import math
import matplotlib.pyplot as plt
from pathlib import Path
import logging
import config

logger = logging.getLogger(__name__)

MIN_TOTAL_MENTIONS = config.MIN_TOTAL_MENTIONS
RECENT_YEARS = config.RECENT_YEARS
UNDERRESEARCH_TOP_N = config.UNDERRESEARCH_TOP_N


def _ensure_int_years(counts):
    """Convert year keys to integers (in case they’re loaded as strings)."""
    return {int(k): v for k, v in counts.items()}


def stacked_trend_top5(df, outdir, term):
    """
    Plot stacked bar chart of top-5 gene–disease pairs in the most recent year,
    with an 'Other' category showing all remaining pairs.
    """
    df = df.copy()
    df["annual_counts"] = df["annual_counts"].apply(_ensure_int_years)

    df = df[df["total_mentions"] >= MIN_TOTAL_MENTIONS]
    if df.empty:
        print(f"[viz] no pairs exceed {MIN_TOTAL_MENTIONS} total mentions — nothing to plot.")
        return

    all_years = sorted({y for d in df["annual_counts"] for y in d})
    last_year = all_years[-1]

    df["count_last_year"] = df["annual_counts"].apply(lambda d: d.get(last_year, 0))
    top5 = df.nlargest(5, "count_last_year")
    labels = [f"{r['gene_name']}–{r['disease_name']}" for _, r in top5.iterrows()]

    data = {label: [] for label in labels}
    data["Other"] = []

    for year in all_years:
        total_all = sum(d.get(year, 0) for d in df["annual_counts"])
        top5_counts = [r["annual_counts"].get(year, 0) for _, r in top5.iterrows()]
        other = total_all - sum(top5_counts)
        for label, count in zip(labels, top5_counts):
            data[label].append(count)
        data["Other"].append(max(other, 0))

    fig, ax = plt.subplots()
    bottom = [0] * len(all_years)
    for label in labels + ["Other"]:
        ax.bar(all_years, data[label], bottom=bottom, label=label)
        bottom = [b + c for b, c in zip(bottom, data[label])]

    ax.set_xlabel("Year")
    ax.set_ylabel("Number of co-mentions")
    ax.set_title(f"Top 5 gene–disease trends for '{term}' (through {last_year})")
    ax.legend(ncol=2, fontsize="small", loc="upper left")
    plt.xticks(all_years[:: max(1, len(all_years) // 10)])
    plt.tight_layout()

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"{term.replace(' ', '_')}_top5_trends.png"
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"[viz] stacked_trend_top5 → {outpath}")


def stacked_trend_top5_percent(df, outdir, term):
    """
    Plot 100% stacked bar chart showing the relative share of top-5 gene–disease pairs over time.
    """
    df = df.copy()
    df["annual_counts"] = df["annual_counts"].apply(_ensure_int_years)

    df = df[df["total_mentions"] >= MIN_TOTAL_MENTIONS]
    if df.empty:
        print(f"[viz] no pairs exceed {MIN_TOTAL_MENTIONS} total mentions — nothing to plot.")
        return

    all_years = sorted({y for d in df["annual_counts"] for y in d})
    last_year = all_years[-1]

    df["count_last_year"] = df["annual_counts"].apply(lambda d: d.get(last_year, 0))
    top5 = df.nlargest(5, "count_last_year")
    labels = [f"{r['gene_name']}–{r['disease_name']}" for _, r in top5.iterrows()]

    pct_data = {label: [] for label in labels}
    pct_data["Other"] = []

    for year in all_years:
        total = sum(d.get(year, 0) for d in df["annual_counts"])
        top_counts = [r["annual_counts"].get(year, 0) for _, r in top5.iterrows()]
        for label, val in zip(labels, top_counts):
            pct_data[label].append((val / total * 100) if total else 0)
        other_pct = ((total - sum(top_counts)) / total * 100) if total else 0
        pct_data["Other"].append(max(other_pct, 0))

    fig, ax = plt.subplots()
    bottom = [0] * len(all_years)
    for label in labels + ["Other"]:
        ax.bar(all_years, pct_data[label], bottom=bottom, label=label)
        bottom = [b + c for b, c in zip(bottom, pct_data[label])]

    ax.set_xlabel("Year")
    ax.set_ylabel("Percentage of co-mentions (%)")
    ax.set_title(f"Top 5 gene–disease relative share for '{term}' (through {last_year})")
    ax.legend(ncol=2, fontsize="small", loc="upper left")
    plt.xticks(all_years[:: max(1, len(all_years) // 10)])
    plt.tight_layout()

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"{term.replace(' ', '_')}_top5_trends_percent.png"
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"[viz] stacked_trend_top5_percent → {outpath}")


def underresearched_pairs(df, outdir, term):
    """
    Identify and visualize top under-researched gene–disease pairs.

    Scoring combines:
      - Low recent activity (low ratio of recent mentions)
      - High burstiness (most mentions clustered in one year)
      - Moderate volume (log-scaled total_mentions)

    score = (1 - recent_ratio) * burstiness * log(total_mentions)
    """
    df2 = df.copy()
    df2["annual_counts"] = df2["annual_counts"].apply(_ensure_int_years)

    all_years = {y for d in df2["annual_counts"] for y in d}
    current_year = max(all_years)

    def recent_count(counts):
        return sum(v for y, v in counts.items() if y >= current_year - RECENT_YEARS + 1)

    df2["recent_count"] = df2["annual_counts"].apply(recent_count)
    df2["recent_ratio"] = df2.apply(
        lambda r: r["recent_count"] / r["total_mentions"] if r["total_mentions"] > 0 else 0,
        axis=1
    )
    df2["last_year"] = df2["annual_counts"].apply(lambda d: max(d.keys()))
    df2["time_since_last"] = df2["last_year"].apply(lambda y: current_year - y)

    df2 = df2[df2["total_mentions"] >= MIN_TOTAL_MENTIONS]
    if df2.empty:
        print(f"[viz] no pairs exceed {MIN_TOTAL_MENTIONS} total mentions — nothing to plot.")
        return

    df2["under_score"] = df2.apply(
        lambda r: (1 - r["recent_ratio"]) * r["burstiness"] * math.log(r["total_mentions"]),
        axis=1
    )

    top = df2.nlargest(UNDERRESEARCH_TOP_N, "under_score")

    labels = [f"{r['gene_name']}–{r['disease_name']}" for _, r in top.iterrows()]
    scores = top["under_score"].tolist()

    fig, ax = plt.subplots(figsize=(10, max(5, len(labels) * 0.6)))
    ax.barh(labels[::-1], scores[::-1])
    ax.set_xlabel("Under-researched score")

    # Use suptitle to avoid it being cut off
    fig.suptitle(f"Under-researched pairs for '{term}'", fontsize=14, ha='center')

    # Leave room for the title
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / f"{term.replace(' ', '_')}_underresearched.png"
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    print(f"[viz] underresearched_pairs → {outpath}")
