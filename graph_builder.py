import networkx as nx
import pickle
import logging
from typing import Optional, List, Dict, Any

from ner import extract_gene_disease_pairs

logger = logging.getLogger(__name__)


def is_valid_gene_id(gene_id: str) -> bool:
    return gene_id.isdigit()


def is_valid_disease_id(disease_id: str) -> bool:
    return (
        disease_id.startswith(("MESH:", "OMIM:")) or
        (disease_id.startswith("D") and disease_id[1:].isdigit())
    )


def build_temporal_graph(
    term: str,
    batch_size: int = 100,
    records: Optional[List[Dict[str, Any]]] = None
) -> nx.Graph:
    """
    Build the gene-disease co-mention graph.
    Uses pre-fetched `records` if provided, otherwise retrieves data.
    """
    pairs = extract_gene_disease_pairs(term, batch_size=batch_size, records=records)

    G = nx.Graph()
    skipped = 0

    for pmid, gene, disease, year in pairs:
        if not is_valid_gene_id(gene) or not is_valid_disease_id(disease):
            logger.warning("Skipping invalid pair %s %s %s", pmid, gene, disease)
            skipped += 1
            continue

        if not G.has_node(gene):
            G.add_node(gene, type="gene")
        if not G.has_node(disease):
            G.add_node(disease, type="disease")

        if G.has_edge(gene, disease):
            G[gene][disease]["pmids"].append(pmid)
            G[gene][disease]["years"].append(year)
        else:
            G.add_edge(gene, disease, pmids=[pmid], years=[year])

    logger.info(
        "Built graph for '%s': %d nodes, %d edges (skipped %d invalid)",
        term, G.number_of_nodes(), G.number_of_edges(), skipped
    )
    return G


def save_graph(G: nx.Graph, filepath: str):
    with open(filepath, "wb") as f:
        pickle.dump(G, f)
    logger.info("Graph saved to %s", filepath)
