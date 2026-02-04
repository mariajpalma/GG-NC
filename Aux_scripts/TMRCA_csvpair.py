#!/usr/bin/env python3
import tskit
import numpy as np
import pandas as pd
import itertools
import json
import multiprocessing
from tqdm import tqdm
import argparse
import logging
import os
import csv

def parse_args():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description="Calculate pairwise TMRCAs between individuals in a tree sequence"
    )
    parser.add_argument("tree_sequence", help="Input .trees file")
    parser.add_argument("output_csv", help="Output CSV file for pairwise TMRCAs")
    parser.add_argument(
        "--pairs-file",
        help="CSV file with specific pairs to process (columns: ID1,ID2)"
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Number of CPU processes to use"
    )
    parser.add_argument(
        "--verbosity",
        "-v",
        action="count",
        default=0,
        help="Increase output verbosity"
    )
    return parser.parse_args()

def extract_individual_id(ind):
    """Your ID extraction method from metadata"""
    metadata = json.loads(ind.metadata) if ind.metadata else {}
    return (
        metadata.get("sample") or 
        metadata.get("individual_id") or 
        metadata.get("sample_id") or 
        str(ind.id)  # Fallback to internal ID
    )

def load_filtered_individuals(ts):
    """Versión optimizada para árboles grandes"""
    individuals = {}
    id_to_nodes = {}
    node_to_ind = {}  # Mapeo auxiliar: nodo → ID individual
    
    # Primera pasada: construir node_to_ind
    for node in ts.nodes():
        if node.individual != -1:
            node_to_ind[node.id] = node.individual
    
    # Segunda pasada: procesar individuos
    for individual in ts.individuals():
        ext_id = extract_individual_id(individual)
        nodes = [n for n, ind in node_to_ind.items() if ind == individual.id]
        
        individuals[individual.id] = {
            'external_id': ext_id,
            'nodes': nodes,
            'metadata': individual.metadata
        }
        id_to_nodes[ext_id] = nodes
    
    return individuals, id_to_nodes


def process_individual_pair(args):
    """Calculate TMRCAs for one pair"""
    pair, id_to_nodes, ts_path = args
    ts = tskit.load(ts_path)
    
    try:
        nodes1 = id_to_nodes[pair[0]]
        nodes2 = id_to_nodes[pair[1]]
    except KeyError:
        return (pair[0], pair[1], np.nan, np.nan, np.nan)
    
    total_tmrca = 0.0
    total_span = 0.0
    min_tmrca = float('inf')
    max_tmrca = -float('inf')
    
    for tree in ts.trees():
        span = tree.span
        for n1 in nodes1:
            for n2 in nodes2:
                mrca_node = tree.mrca(n1, n2)
                tmrca = tree.time(mrca_node)
                total_tmrca += tmrca * span
                total_span += span
                
                if tmrca < min_tmrca:
                    min_tmrca = tmrca
                if tmrca > max_tmrca:
                    max_tmrca = tmrca
    
    mean_tmrca = total_tmrca / total_span if total_span > 0 else np.nan
    return (pair[0], pair[1], mean_tmrca, min_tmrca, max_tmrca)

def read_pairs_file(pairs_file):
    """Read pairs from CSV file"""
    with open(pairs_file) as f:
        reader = csv.reader(f)
        # Skip header if present
        try:
            if 'ID1' in next(reader):
                pass
        except StopIteration:
            pass
        return [tuple(row[:2]) for row in reader]

def calculate_pairwise_tmrcas(ts_path, id_to_nodes, pairs, num_processes):
    """Parallel processing of specified pairs"""
    args = [(pair, id_to_nodes, ts_path) for pair in pairs]
    
    results = []
    with multiprocessing.Pool(processes=num_processes) as pool:
        for result in tqdm(
            pool.imap(process_individual_pair, args),
            total=len(pairs),
            desc="Processing pairs"
        ):
            results.append(result)
    return results

def main():
    args = parse_args()
    
    # Configure logging
    log_level = logging.WARNING if args.verbosity == 0 else logging.INFO
    logging.basicConfig(level=log_level, format="%(message)s")
    
    logging.info(f"Loading tree sequence from {args.tree_sequence}")
    ts = tskit.load(args.tree_sequence)
    
    logging.info("Identifying individuals...")
    individuals, id_to_nodes = load_filtered_individuals(ts)
    
    if args.pairs_file:
        logging.info(f"Reading pairs from {args.pairs_file}")
        pairs = read_pairs_file(args.pairs_file)
    else:
        logging.info("Generating all possible pairs")
        pairs = list(itertools.combinations(id_to_nodes.keys(), 2))
    
    logging.info(f"Processing {len(pairs)} pairs")
    results = calculate_pairwise_tmrcas(
        args.tree_sequence,
        id_to_nodes,
        pairs,
        args.processes
    )
    
    logging.info(f"Saving results to {args.output_csv}")
    df = pd.DataFrame(results, columns=['ID1', 'ID2', 'tmrca_mean', 'tmrca_min', 'tmrca_max'])
    df.to_csv(args.output_csv, index=False)
    logging.info("Done!")

if __name__ == "__main__":
    main()
