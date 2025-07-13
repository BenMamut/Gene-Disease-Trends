# Gene-Disease-Trends

**Gene-Disease-Trends** is a Python pipeline for analyzing temporal trends in gene–disease co-mentions across biomedical literature. It extracts and aggregates co-mentions from PubMed abstracts and highlights how interest in specific gene–disease relationships has evolved over time. The project also identifies potentially under-researched gene–disease links using a heuristic scoring method.

## Overview

This project performs:
- Automated year-by-year data retrieval from PubTator and PubMed
- Named Entity Recognition (NER) to extract gene–disease pairs
- Graph-based aggregation of co-mention frequency over time
- Temporal analysis of research interest and burstiness
- Heuristic-based detection of potentially under-researched links
- Visualizations of trends and research gaps

## Why It Matters

Analyzing temporal co-mention trends can:
- Surface research areas that may be overlooked
- Help researchers identify candidates for further study
- Track the evolution of biomedical focus areas
- Provide a lightweight, scalable way to explore literature trends

## Visual Outputs

The pipeline generates:
- Top 5 gene–disease co-mention trends (raw count per year)
- Top 5 co-mention percentage share over time
- Top potentially under-researched gene–disease pairs

All plots and reports are saved to the `output/` directory.

## How to Run

Step 1: Clone the repository  
→ `git clone https://github.com/BenMamut/Gene-Disease-Trends.git`  
→ `cd Gene-Disease-Trends`

Step 2: Install required dependencies  

Step 3: Configure search terms, dates, and pipeline activity
→ Edit `config.py` and modify the `TERMS` list with queries of interest.  
Example:
```python
TERMS = [
    "Sickle Cell Disease",
    "BRCA1 Cancer",
    "Alzheimer's Disease",
]
```

Step 4: Run the pipeline  
→ `python pipeline.py`

## Output

- CSV report of temporal co-mention metrics per gene–disease pair
- PNG visualizations of yearly trends and under-researched candidates
- All files saved under output/ (configurable)

## Files

- `pipeline.py`: Main orchestration of the pipeline
- `data_retrieval.py`: Pulls abstracts via PubMed
- `ner.py`: Analyzes abstracts via PubTator
- `graph_builder.py`: Builds the temporal co-mention graph
- `temporal_analysis.py`: Extracts trends, burstiness, and metrics
- `visualizations.py`: Generates output graphs
- `gene_disease_utils.py`: Maps gene/disease IDs to readable names
- `config.py`: Central configuration file

## Limitations

- Uses abstract-level annotations; full-text data is not included
- Heuristic scoring may yield false positives or biologically irrelevant pairs
- Focuses on temporal patterns, not mechanistic relationships
