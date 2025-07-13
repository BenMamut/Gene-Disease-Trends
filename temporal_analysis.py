import pickle
import logging
import json
from collections import Counter

import pandas as pd

logger = logging.getLogger(__name__)


def load_graph(filepath: str):
    with open(filepath, "rb") as f:
        G = pickle.load(f)
    logger.info("Loaded graph with %d nodes and %d edges",
                G.number_of_nodes(), G.number_of_edges())
    return G


def analyze_temporal_edges(G) -> pd.DataFrame:
    """
    Analyze the temporal behavior of gene–disease links.
    Computes first mention, peak year, burstiness, and annual counts.
    """
    records = []
    skipped = 0

    for u, v, attrs in G.edges(data=True):
        tu = G.nodes[u].get("type")
        tv = G.nodes[v].get("type")

        if tu == "gene" and tv == "disease":
            gene, disease = u, v
        elif tv == "gene" and tu == "disease":
            gene, disease = v, u
        else:
            logger.warning("Edge %s-%s has bad types %s-%s", u, v, tu, tv)
            skipped += 1
            continue

        years = [int(y) for y in attrs.get("years", []) if y.isdigit()]
        if not years:
            continue

        cnt = Counter(years)
        first = min(cnt)
        peak, peak_count = max(cnt.items(), key=lambda x: x[1])
        total = sum(cnt.values())
        burst = peak_count / total if total else 0.0

        records.append({
            "gene": gene,
            "disease": disease,
            "first_mention": first,
            "peak_year": peak,
            "time_to_peak": peak - first,
            "total_mentions": total,
            "peak_count": peak_count,
            "burstiness": burst,
            "annual_counts": dict(sorted(cnt.items()))
        })

    logger.info("Analyzed %d edges (skipped %d)", len(records), skipped)
    df = pd.DataFrame(records)
    return df.sort_values("total_mentions", ascending=False).reset_index(drop=True)


def save_report(df: pd.DataFrame, csv_path: str):
    """
    Save the DataFrame to CSV, serializing the annual_counts column.
    """
    df_out = df.copy()
    df_out["annual_counts"] = df_out["annual_counts"].apply(json.dumps)
    df_out.to_csv(csv_path, index=False)
    logger.info("Report saved (with annual_counts) to %s", csv_path)
