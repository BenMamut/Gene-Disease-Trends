import os
import json
import time
import logging
import requests
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Iterator
from datetime import datetime

import config

# PubMed E-utilities
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Ensure data directory exists
DATA_DIR = config.DATA_DIR
DATA_DIR.mkdir(exist_ok=True)

logger = logging.getLogger(__name__)


def _term_file(term: str) -> Path:
    safe = term.replace(" ", "_")
    return DATA_DIR / f"{safe}.json"


def _load_data(term: str) -> Dict:
    path = _term_file(term)
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {"completed_years": [], "records": []}


def _save_data(term: str, data: Dict):
    path = _term_file(term)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def _chunked(iterable: List, size: int) -> Iterator[List]:
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]


def _search_year(term: str, year: int, max_retries: int = 3) -> List[str]:
    """
    Search PubMed for PMIDs in a given year, with retries on HTTP 429.
    """
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": 100000,
        "retmode": "json",
        "datetype": "pdat",
        "mindate": str(year),
        "maxdate": str(year),
        **({"api_key": os.getenv("NCBI_API_KEY")} if os.getenv("NCBI_API_KEY") else {})
    }

    for attempt in range(1, max_retries + 1):
        resp = requests.get(ESEARCH_URL, params=params)
        if resp.status_code == 429:
            wait = 2 ** attempt
            logger.warning("HTTP 429 on year %d (attempt %d/%d), retrying in %ds", year, attempt, max_retries, wait)
            time.sleep(wait)
            continue
        try:
            resp.raise_for_status()
        except requests.HTTPError:
            logger.error("HTTP error for year %d: %s", year, resp.status_code)
            return []
        return resp.json().get("esearchresult", {}).get("idlist", [])

    logger.error("Failed to fetch PMIDs for %d after %d retries", year, max_retries)
    return []


def _fetch_pmids(pmids: List[str], batch_size: int, max_retries: int = 3) -> List[Dict]:
    records: List[Dict] = []
    for chunk in _chunked(pmids, batch_size):
        for attempt in range(1, max_retries + 1):
            params = {
                "db": "pubmed",
                "id": ",".join(chunk),
                "retmode": "xml",
                **({"api_key": os.getenv("NCBI_API_KEY")} if os.getenv("NCBI_API_KEY") else {})
            }
            resp = requests.get(EFETCH_URL, params=params)
            if resp.status_code == 429:
                wait = 2 ** attempt
                logger.warning("HTTP 429 for efetch (attempt %d/%d), retrying in %ds", attempt, max_retries, wait)
                time.sleep(wait)
                continue
            try:
                resp.raise_for_status()
            except requests.HTTPError:
                logger.error("Failed to fetch PMIDs: %s", resp.status_code)
                break

            root = ET.fromstring(resp.text)
            for art in root.findall(".//PubmedArticle"):
                pmid = art.findtext(".//PMID", default="")
                title = art.findtext(".//ArticleTitle", default="")
                abstract = " ".join(p.text or "" for p in art.findall(".//AbstractText"))
                year = (
                    art.findtext(".//PubDate/Year")
                    or art.findtext(".//PubDate/MedlineDate", default="").split(" ")[0]
                )
                records.append({
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "year": year
                })
            time.sleep(0.34)
            break  # break retry loop on success
    return records


def retrieve_full_data(term: str,
                       start_year: int = config.START_YEAR,
                       end_year: int = None) -> List[Dict]:
    """
    Retrieve all PubMed records for a term across a year range, using caching and retries.
    """
    end_year = end_year or datetime.now().year
    data = _load_data(term)
    done = set(int(y) for y in data.get("completed_years", []))
    existing = {r["pmid"] for r in data.get("records", [])}
    updated = False

    for year in range(start_year, end_year + 1):
        is_current = (year == datetime.now().year)

        if year in done and not is_current:
            logger.info("Skipping year %d (already done)", year)
            continue

        logger.info("Year %d: searching PubMed for '%s'...", year, term)
        pmids = _search_year(term, year)

        if not pmids:
            logger.info(" → no PMIDs found for %d", year)
            if not is_current:
                data["completed_years"].append(year)
                updated = True
            continue

        new_pmids = [p for p in pmids if p not in existing]
        if not new_pmids:
            logger.info(" → no new PMIDs in %d", year)
            if not is_current:
                data["completed_years"].append(year)
                updated = True
            continue

        logger.info(" → fetching %d abstracts", len(new_pmids))
        recs = _fetch_pmids(new_pmids, config.BATCH_SIZE)
        data["records"].extend(recs)
        existing.update(new_pmids)

        if not is_current:
            data["completed_years"].append(year)

        updated = True

    if updated:
        _save_data(term, data)
        logger.info("Data file updated for '%s'", term)
    else:
        logger.info("Data file for '%s' already up to date", term)

    return data["records"]
