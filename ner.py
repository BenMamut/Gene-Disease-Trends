import time
import logging
import requests
import json
from typing import List, Tuple, Dict, Optional

from data_retrieval import retrieve_full_data

# PubTator BioC JSON export endpoint
PUBTATOR_URL = (
    "https://www.ncbi.nlm.nih.gov/research/pubtator-api/"
    "publications/export/biocjson"
)

logger = logging.getLogger(__name__)


def _split_ids(identifier: str) -> List[str]:
    """Split an identifier string on semicolons or commas."""
    for sep in (';', ','):
        if sep in identifier:
            return [p.strip() for p in identifier.split(sep) if p.strip()]
    return [identifier.strip()]


def annotate_pmids(pmids: List[str], batch_size: int = 100) -> Dict[str, Dict]:
    annotations: Dict[str, Dict] = {}

    for i in range(0, len(pmids), batch_size):
        chunk = pmids[i:i + batch_size]
        resp = requests.get(PUBTATOR_URL, params={
            "pmids": ",".join(chunk),
            "concepts": "Gene,Disease"
        })

        try:
            data = resp.json()
        except json.JSONDecodeError:
            logger.warning("Non‑JSON response for PMIDs %s, skipping batch.", chunk)
            time.sleep(0.3)
            continue

        docs = data.get("documents") or data.get("PubTator3") or []

        for doc in docs:
            pmid = doc.get("id")
            ann = {"genes": [], "diseases": []}

            for passage in doc.get("passages", []):
                text = passage.get("text", "")
                denots = (
                    passage.get("annotations")
                    or passage.get("denotations")
                    or passage.get("annotation")
                    or []
                )

                for inf in denots:
                    infons = inf.get("infons", {})
                    t = infons.get("type", "").lower()
                    raw_id = infons.get("identifier", "").strip()

                    if not raw_id:
                        continue

                    locs = inf.get("locations") or inf.get("location") or []
                    if not locs:
                        continue

                    offset = locs[0]["offset"]
                    length = locs[0]["length"]
                    snippet = text[offset:offset + length]

                    for ident in _split_ids(raw_id):
                        if t == "gene":
                            ann["genes"].append((offset, offset + length, snippet, ident))
                        elif t == "disease":
                            ann["diseases"].append((offset, offset + length, snippet, ident))

            annotations[pmid] = ann

        time.sleep(0.3)

    return annotations


def extract_gene_disease_pairs(
    term: str,
    batch_size: int = 100,
    records: Optional[List[Dict]] = None
) -> List[Tuple[str, str, str, str]]:
    """
    Full pipeline: retrieve → annotate → extract co‑mentions.
    If `records` is provided, retrieval is skipped.
    Returns list of (pmid, gene_id, disease_id, year).
    """
    if records is None:
        records = retrieve_full_data(term)

    pmids = [r["pmid"] for r in records]
    year_map = {r["pmid"]: r["year"] for r in records}
    ann = annotate_pmids(pmids, batch_size=batch_size)

    results: List[Tuple[str, str, str, str]] = []

    for pmid, tags in ann.items():
        genes = {gid for *_, gid in tags["genes"] if gid.isdigit()}

        diseases = set()
        for *_, did in tags["diseases"]:
            if ':' in did:
                prefix, core = did.split(':', 1)
                if core.isalnum():
                    diseases.add(did)
            elif did.startswith("D") and did[1:].isdigit():
                diseases.add(did)

        bad_genes = {gid for *_, gid in tags["genes"] if not gid.isdigit()}
        bad_diseases = {did for *_, did in tags["diseases"] if did and did not in diseases}

        if bad_genes or bad_diseases:
            logger.warning(
                "PMID %s: filtered bad genes=%s bad diseases=%s",
                pmid, bad_genes, bad_diseases
            )

        for g in genes:
            for d in diseases:
                results.append((pmid, g, d, year_map.get(pmid, "")))

    logger.info("Extracted %d gene–disease pairs for '%s'", len(results), term)
    return results
