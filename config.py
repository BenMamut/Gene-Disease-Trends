from pathlib import Path

START_YEAR = 1990        # earliest year to fetch
END_YEAR = None          # latest year to fetch; if None, defaults to current calendar year
BATCH_SIZE = 100         # how many IDs to fetch from PubTator per request

TERMS = [
    "Sickle cell disease",
]

OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

PIPELINE_STEPS = {
    "retrieve":  True,   # pull/update the data file
    "annotate":  True,   # run NER
    "graph":     True,   # build graph
    "analyze":   True,   # temporal analysis
    "report":    True,   # save CSV
    "visualize": True,    # stub only
}

VISUALIZATIONS = [
    "stacked_trend_top5",
    "stacked_trend_top5_percent",
    "underresearched_pairs",
]

MIN_TOTAL_MENTIONS = 10           # minimum mentions to include in plots
RECENT_YEARS = 5                  # how many recent years to consider for underresearched
UNDERRESEARCH_TOP_N = 10          # number of underresearched pairs to show
