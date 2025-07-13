import os
import re
import time
import json
import logging
import requests
from pathlib import Path
from typing import List, Dict, Optional

# Entrez and MeSH endpoints
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
MESH_LOOKUP_URL = "https://id.nlm.nih.gov/mesh/lookup/descriptor"

# Cache directory for resolved names
CACHE_DIR = Path.home() / ".temporal_gdg_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

logger = logging.getLogger(__name__)

# Valid MeSH IDs: D###### or C######
MESH_UI_PATTERN = re.compile(r"^[DC]\d{6}$")


def _cache_path(prefix: str, id_: str) -> Path:
    return CACHE_DIR / f"{prefix}_{id_}.json"


def _load_cached_name(prefix: str, id_: str) -> Optional[str]:
    path = _cache_path(prefix, id_)
    if path.exists():
        return json.loads(path.read_text()).get("name")
    return None


def _save_cached_name(prefix: str, id_: str, name: str):
    _cache_path(prefix, id_).write_text(json.dumps({"name": name}, indent=2))


def _batch_fetch(db: str, ids: List[str], batch_size: int = 200) -> Dict[str, Optional[dict]]:
    out = {}
    for i in range(0, len(ids), batch_size):
        chunk = ids[i:i + batch_size]
        params = {
            "db": db,
            "id": ",".join(chunk),
            "retmode": "json",
            **({"api_key": os.getenv("NCBI_API_KEY")} if os.getenv("NCBI_API_KEY") else {}),
        }
        resp = requests.get(ESUMMARY_URL, params=params)
        if resp.ok:
            data = resp.json().get("result", {})
            for uid in data.get("uids", []):
                out[uid] = data.get(uid)
        else:
            logger.warning("ESummary failed for %s: HTTP %s", chunk, resp.status_code)
            for uid in chunk:
                out[uid] = None
        time.sleep(0.34)
    return out


def get_gene_names(gene_ids: List[str]) -> Dict[str, str]:
    """
    Resolve Entrez gene IDs to names via ESummary.
    Uses disk cache to avoid repeated lookups.
    """
    unique = [gid for gid in set(gene_ids) if gid.isdigit()]
    res, to_fetch = {}, []

    for gid in unique:
        cached = _load_cached_name("gene", gid)
        if cached:
            res[gid] = cached
        else:
            to_fetch.append(gid)

    if to_fetch:
        summaries = _batch_fetch("gene", to_fetch)
        for gid, summary in summaries.items():
            name = summary.get("name") if summary else gid
            res[gid] = name
            _save_cached_name("gene", gid, name)

    return res


def _batch_fetch_mesh_names(cores: List[str], batch_size: int = 50) -> Dict[str, str]:
    """
    Batch lookup for MeSH UI codes using the descriptor API.
    Returns mapping from code → name.
    """
    out = {}
    valid = [c for c in cores if MESH_UI_PATTERN.match(c)]

    for i in range(0, len(valid), batch_size):
        chunk = valid[i:i + batch_size]
        resp = requests.get(MESH_LOOKUP_URL, params={"descriptor": ",".join(chunk)})

        if resp.status_code == 400:
            logger.warning("MeSH lookup bad request for %s; skipping.", chunk)
            continue

        if resp.ok:
            for entry in resp.json():
                ui = entry.get("descriptor")
                label = entry.get("label")
                if ui and label:
                    out[ui] = label
        else:
            logger.warning("MeSH lookup failed for %s: HTTP %s", chunk, resp.status_code)

        time.sleep(0.34)

    return out


def get_disease_names(mesh_ids: List[str]) -> Dict[str, str]:
    """
    Resolve MeSH/OMIM IDs to readable labels.
    Valid MeSH codes (D###### or C######) are resolved via batch API.
    Others are returned as-is.
    """
    cores = []
    res: Dict[str, str] = {}

    for mid in set(mesh_ids):
        prefix, core = mid.split(":", 1) if ":" in mid else (None, mid)

        if MESH_UI_PATTERN.match(core):
            cached = _load_cached_name("mesh", core)
            if cached and cached != core:
                res[mid] = cached
            else:
                cores.append(core)
        else:
            res[mid] = mid  # nonstandard, keep raw ID

    if cores:
        labels = _batch_fetch_mesh_names(cores)
        for core in cores:
            name = labels.get(core) or core
            if core == name:
                logger.warning("No MeSH label for %s; using code.", core)
            _save_cached_name("mesh", core, name)

            for mid in mesh_ids:
                if mid.split(":", 1)[-1] == core:
                    res[mid] = name

    return res
