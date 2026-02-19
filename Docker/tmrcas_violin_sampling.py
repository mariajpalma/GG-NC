import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import time
import argparse
from scipy import stats
import os

# ====================== FUNCIONES OPTIMIZADAS ======================
def read_community_matrix(matrix_file):
    """Read the community assignment matrix."""
    return np.loadtxt(matrix_file, dtype=np.int32, delimiter=',')

def read_id_file(id_file):
    """Read the individual ID file (one ID per line)."""
    with open(id_file) as f:
        return [line.strip() for line in f if line.strip()]

def read_tmrca(tmrca_file):
    """Read the pairwise TMRCA file with selected columns."""
    return pd.read_csv(
        tmrca_file,
        sep=",", 
        header=0,
        usecols=['ID1', 'ID2', 'tmrca_mean'],
        dtype={'ID1': 'category', 'ID2': 'category', 'tmrca_mean': 'float32'}
    )

def process_communities(community_matrix, ids, exclude_label=1):
    """Construct per-resolution community membership dictionaries."""
    return [
        {
            comm: [ids[idx] for idx in np.where(row == comm)[0] if comm != exclude_label]
            for comm in np.unique(row) if comm != exclude_label
        }
        for row in community_matrix
    ]

def sample_community(args):
    comm_id, members, tmrca_df, n_samples, sample_size = args
    
    if len(members) < sample_size:
        return None
    
    results = []
    members_arr = np.array(members)
    
    for _ in range(n_samples):
        sampled = np.random.choice(members_arr, size=sample_size, replace=True)
        pairs = {(min(m1, m2), max(m1, m2)) 
                for m1, m2 in combinations(sampled, 2) if m1 != m2}
        
        if not pairs:  # Si no hay pares válidos
            continue
            
        try:
            pairs_df = pd.DataFrame(list(pairs), columns=['ID1', 'ID2'])
            merged = pd.merge(
                pairs_df,
                tmrca_df[['ID1', 'ID2', 'tmrca_mean']],  # Selección explícita
                on=['ID1', 'ID2'],
                how='inner'
            )
            
            if not merged.empty:
                results.append({
                    'Community': f"C{comm_id}",
                    'Mean_TMRCA': merged['tmrca_mean'].mean(),
                    'Size': len(members)
                })
        except KeyError as e:
            print(f"Error en merge: {e}")
            print("Pairs sample:", pairs)
            print("TMRCA columns:", tmrca_df.columns)
            continue
    
    return results if results else None


def process_resolution(args):
    """Process all communities within a single resolution."""
    resolution, communities, tmrca_df, n_samples, sample_size = args
    resolution_results = []
    
    for comm_id, members in communities.items():
        samples = sample_community((comm_id, members, tmrca_df, n_samples, sample_size))
        if samples:
            for sample in samples:
                sample['Resolution'] = resolution
                sample['Resolution_Label'] = f"Resolution {resolution}"
            resolution_results.extend(samples)
    
    return (resolution, resolution_results)

def init_worker():
    """Initialize each worker process with a unique random seed."""
    np.random.seed(os.getpid())



def calculate_sampled_tmrcas(community_results, tmrca_df, n_samples=100, sample_size=6, requested_cores=None, max_cores=8):
    """Compute sampled TMRCA values across all resolutions."""
    tasks = [
        (res, communities, tmrca_df, n_samples, sample_size)
        for res, communities in enumerate(community_results)
    ]
    
    # Ajuste dinámico de cores
    available_cores = cpu_count()
    cap = max_cores if max_cores is not None else available_cores
    req = requested_cores if requested_cores is not None else cap
    optimal_cores = max(1, min(req, cap, available_cores, len(tasks)))
    
    print(f"Using {optimal_cores} worker processes for {len(tasks)} resolutions "
      f"(requested={requested_cores}, cap={cap}, available={available_cores})")
    
    results_by_resolution = defaultdict(list)
    with Pool(optimal_cores, initializer=init_worker) as pool:
        for res, res_results in pool.imap(process_resolution, tasks):
            results_by_resolution[res] = res_results
    
    # Reconstrucción ordenada
    final_data = []
    for res in sorted(results_by_resolution.keys()):
        final_data.extend(results_by_resolution[res])
    
    return pd.DataFrame(final_data)



# ====================== NUEVA FUNCIÓN DE PLOTEO ======================
def plot_sampled_tmrca(df, output_file, upper_lim, lower_lim, n_resolutions):
    """Plot median sampled TMRCA per resolution."""
    lower_lim = float(lower_lim)
    upper_lim = float(upper_lim)
    integer_values = np.arange(lower_lim, upper_lim + 1)

    # Calcular posiciones de las etiquetas
    positions = [
        int(round((val - lower_lim) / (upper_lim - lower_lim) * (n_resolutions - 1)))
        for val in integer_values
    ]

    
    # Calcular estadísticas por resolución
    stats_df = df.groupby('Resolution')['Mean_TMRCA'].agg(
        median='median',
        mean='mean',
        q025=lambda x: np.quantile(x, 0.025),
        q975=lambda x: np.quantile(x, 0.975),
        count='count'
    ).reset_index()


    with PdfPages(output_file) as pdf:
        plt.figure(figsize=(14, 8))
        
        # Línea de medianas
        plt.plot(stats_df['Resolution'], stats_df['median'], 
                'b-', linewidth=2, label='Median')
        
        # Intervalos de error estándar
        plt.fill_between(
            x=stats_df['Resolution'],
            y1=stats_df['q025'],
            y2=stats_df['q975'],
            color= "gray",
            alpha=0.3,
            label='95% confidence interval'
        )
        
        # Configurar ejes
        plt.xticks(
            ticks=positions,
            labels=[str(val) for val in integer_values],
            rotation=45
        )
        
        for pos in positions:
            plt.axvline(x=pos, color="gray", linestyle=":", alpha=0.3)
        
        plt.title("TMRCA per resolution (sampling with replacement)", fontsize=14)
        plt.xlabel("Resolution", fontsize=12)
        plt.ylabel("TMRCA", fontsize=12)
        plt.legend()
        plt.grid(True, axis='y', linestyle='--', alpha=0.6)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

# ====================== MAIN ======================
def main(matrix_file, id_file, tmrca_file, output_name, upper_lim, lower_lim, cores=None, max_cores=8):
    print("=== ANALYSIS STARTED ===")
    start_time = time.time()
    
    # 1. Lectura de datos
    print("Reading input files...")
    t0 = time.time()
    community_matrix = read_community_matrix(matrix_file)
    ids = read_id_file(id_file)
    tmrca_df = read_tmrca(tmrca_file)
    print(f"✓ Input files loaded in {time.time()-t0:.2f} seconds")
    
    # 2. Procesamiento de comunidades
    print("Processing communities...")
    t0 = time.time()
    community_results = process_communities(community_matrix, ids)
    print(f"✓ Communities processed in {time.time()-t0:.2f} seconds")
    
    # 3. Cálculo de TMRCAs con muestreo (paralelo)
    print("Computing sampled TMRCAs (parallel execution)...")
    t0 = time.time()
    tmrca_samples = calculate_sampled_tmrcas(community_results, tmrca_df, requested_cores=cores)
    print(f"✓ Sampled TMRCAs computed in {time.time()-t0:.2f} seconds")
    
    # 4. Guardar resultados
    tmrca_samples.to_csv(output_name + "_samples.csv", index=False)
    
    # 5. Generar gráficos
    print("Generating plots...")
    t0 = time.time()
    n_resolutions = community_matrix.shape[0]
    plot_sampled_tmrca(tmrca_samples, output_name + ".pdf", upper_lim, lower_lim, n_resolutions)
    print(f"✓ Plots generated in {time.time()-t0:.2f} seconds")
    elapsed = time.time() - start_time
    print(f"\n=== ANALYSIS COMPLETED IN {elapsed:.2f} SECONDS ===")
    print("Outputs written to:")
    print(f"- {output_name}_samples.csv")
    print(f"- {output_name}.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Análisis de comunidades basado en TMRCA con muestreo')
    parser.add_argument('matrix_file', help='CSV file containing the community assignment matrix (resolutions × individuals).')
    parser.add_argument('id_file', help='Text file containing individual IDs (one per line, ordered to match the matrix columns)')
    parser.add_argument('tmrca_file', help='CSV file with pairwise TMRCA estimates (columns: ID1, ID2, tmrca_mean).')
    parser.add_argument('output_name', help='Base name (prefix) for output files.')
    parser.add_argument('upper_lim', help="Upper bound of the explored resolution space (for x-axis labeling).")
    parser.add_argument('lower_lim', help="Lower bound of the explored resolution space (for x-axis labeling).")
    parser.add_argument("cores", type=int, default=None, help="Requested number of worker processes (will be capped).")
    
    args = parser.parse_args()
    
    main(args.matrix_file, args.id_file, args.tmrca_file, args.output_name, args.upper_lim, args.lower_lim, args.cores)
