import logging
import pandas as pd

import config
from data_retrieval import retrieve_full_data
from graph_builder import build_temporal_graph, save_graph
from temporal_analysis import load_graph, analyze_temporal_edges, save_report
from gene_disease_utils import get_gene_names, get_disease_names
import visualizations

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    for term in config.TERMS:
        safe_name = term.replace(" ", "_")
        graph_path = config.OUTPUT_DIR / f"{safe_name}_graph.pkl"
        report_path = config.OUTPUT_DIR / f"{safe_name}_temporal_report.csv"

        # Step 1: Retrieve
        if config.PIPELINE_STEPS.get("retrieve"):
            records = retrieve_full_data(
                term,
                start_year=config.START_YEAR,
                end_year=config.END_YEAR
            )
        else:
            logger.info("Skipping data retrieval.")
            records = None

        # Step 2: Build graph
        if config.PIPELINE_STEPS.get("graph"):
            G = build_temporal_graph(
                term,
                batch_size=config.BATCH_SIZE,
                records=records
            )
            save_graph(G, str(graph_path))
        else:
            logger.info("Skipping graph build.")

        # Step 3: Analyze
        if config.PIPELINE_STEPS.get("analyze"):
            G = load_graph(str(graph_path))
            df = analyze_temporal_edges(G)
        else:
            logger.info("Skipping analysis; loading existing CSV.")
            df = pd.read_csv(report_path)

        # Step 4: Report
        if config.PIPELINE_STEPS.get("report"):
            gene_map = get_gene_names(df["gene"].tolist())
            disease_map = get_disease_names(df["disease"].tolist())
            df["gene_name"] = df["gene"].map(gene_map)
            df["disease_name"] = df["disease"].map(disease_map)
            save_report(df, str(report_path))
        else:
            logger.info("Skipping report save.")

        # Step 5: Visualizations
        if config.PIPELINE_STEPS.get("visualize"):
            if not config.PIPELINE_STEPS.get("analyze"):
                if "annual_counts" not in df.columns:
                    logger.error("'annual_counts' missing in %s", report_path)
                    continue
                df["annual_counts"] = df["annual_counts"].apply(eval)

            for viz_name in config.VISUALIZATIONS:
                viz_fn = getattr(visualizations, viz_name, None)
                if not viz_fn:
                    logger.warning("Visualization '%s' not found.", viz_name)
                    continue
                try:
                    viz_fn(df, config.OUTPUT_DIR, safe_name)
                except Exception as e:
                    logger.error("Error running visualization '%s': %s", viz_name, e)
        else:
            logger.info("Skipping visualizations.")

        print(f"→ Pipeline complete for '{term}'. Outputs in {config.OUTPUT_DIR}")


if __name__ == "__main__":
    main()
