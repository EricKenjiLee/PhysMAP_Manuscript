"""
PhysMAP Utility Functions
Core processing functions for multi-modal neural data analysis.
Extracted from processJianing.ipynb for use in Streamlit app.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from anndata import AnnData
from muon import MuData
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Tuple, Optional, List
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# Data Processing Functions
# =============================================================================

def clr_normalize(X: np.ndarray) -> np.ndarray:
    """
    Centered Log-Ratio (CLR) normalization per cell (row).
    Equivalent to Seurat's NormalizeData(method='CLR', margin=4).

    Args:
        X: Data matrix (cells x features)

    Returns:
        CLR normalized matrix
    """
    X_pseudo = X + 1
    log_X = np.log(X_pseudo)
    geometric_mean = np.mean(log_X, axis=1, keepdims=True)
    X_clr = log_X - geometric_mean
    return X_clr


def load_csv_data(
    wf_file,
    isi_file,
    psth_file,
    metadata_file
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """
    Load CSV files for each modality and metadata.

    Args:
        wf_file: File-like object or path to waveform CSV
        isi_file: File-like object or path to ISI CSV
        psth_file: File-like object or path to PSTH CSV
        metadata_file: File-like object or path to metadata CSV

    Returns:
        Tuple of (modality_dict, metadata_df)
    """
    # Load each modality
    wf_df = pd.read_csv(wf_file, index_col=0)
    isi_df = pd.read_csv(isi_file, index_col=0)
    psth_df = pd.read_csv(psth_file, index_col=0)
    metadata_df = pd.read_csv(metadata_file, index_col=0)

    # Ensure consistent cell ordering
    common_cells = list(set(wf_df.index) & set(isi_df.index) &
                        set(psth_df.index) & set(metadata_df.index))
    common_cells.sort()

    modality_dict = {
        'WF': wf_df.loc[common_cells],
        'ISI': isi_df.loc[common_cells],
        'PSTH': psth_df.loc[common_cells]
    }

    metadata_df = metadata_df.loc[common_cells]

    return modality_dict, metadata_df


def create_anndata_from_df(
    data_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    modality_name: str
) -> AnnData:
    """
    Create an AnnData object from a pandas DataFrame.

    Args:
        data_df: Feature matrix (cells x features)
        metadata_df: Cell metadata
        modality_name: Name of the modality

    Returns:
        AnnData object
    """
    adata = AnnData(
        X=data_df.values.astype(np.float32),
        obs=metadata_df.copy()
    )
    adata.var_names = [f"{modality_name}_{i}" for i in range(data_df.shape[1])]
    adata.obs_names = data_df.index.tolist()

    return adata


def process_modality(
    adata: AnnData,
    modality_name: str,
    n_pcs: int = 30,
    normalize: bool = True,
    n_neighbors: int = 20,
    metric: str = 'cosine',
    min_dist: float = 0.2,
    resolution: float = 2.0,
    random_state: int = 42
) -> AnnData:
    """
    Process a single modality through the full pipeline.

    Steps:
    1. CLR normalization (optional)
    2. Scale data
    3. PCA
    4. Find neighbors
    5. UMAP (2D)
    6. Leiden clustering

    Args:
        adata: AnnData object
        modality_name: Name of the modality
        n_pcs: Number of principal components
        normalize: Whether to apply CLR normalization
        n_neighbors: Number of neighbors for kNN graph
        metric: Distance metric for neighbors
        min_dist: UMAP minimum distance
        resolution: Leiden clustering resolution
        random_state: Random seed

    Returns:
        Processed AnnData object
    """
    adata = adata.copy()

    # Store raw data
    adata.layers['raw'] = adata.X.copy()

    # CLR normalization
    if normalize:
        adata.X = clr_normalize(adata.X)
        adata.layers['normalized'] = adata.X.copy()

    # Use all features as variable
    adata.var['highly_variable'] = True

    # Scale data
    sc.pp.scale(adata, max_value=None)

    # PCA
    n_pcs_actual = min(n_pcs, min(adata.n_obs, adata.n_vars) - 1)
    sc.tl.pca(adata, n_comps=n_pcs_actual, svd_solver='arpack', random_state=random_state)

    # Find neighbors
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs_actual,
                    metric=metric, random_state=random_state)

    # UMAP
    sc.tl.umap(adata, random_state=random_state, min_dist=min_dist)

    # Leiden clustering
    sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
    adata.obs[f'{modality_name}_leiden'] = adata.obs['leiden']

    return adata


def run_wnn_integration(
    adata_dict: Dict[str, AnnData],
    n_neighbors: int = 20,
    n_bandwidth_neighbors: int = 20,
    min_dist: float = 0.2,
    resolution: float = 2.0,
    random_state: int = 42
) -> Tuple[MuData, AnnData]:
    """
    Run weighted nearest neighbors integration on multiple modalities.

    Args:
        adata_dict: Dictionary mapping modality names to AnnData objects
        n_neighbors: Number of neighbors for WNN
        n_bandwidth_neighbors: Bandwidth neighbors for WNN
        min_dist: UMAP minimum distance
        resolution: Leiden clustering resolution
        random_state: Random seed

    Returns:
        Tuple of (MuData object, AnnData for clustering)
    """
    # Create MuData
    mdata = MuData(adata_dict)

    # Run Muon's WNN
    mu.pp.neighbors(
        mdata,
        key_added='wnn',
        n_neighbors=n_neighbors,
        n_bandwidth_neighbors=n_bandwidth_neighbors,
        n_multineighbors=200,
        metric='euclidean',
        low_memory=False
    )

    # Compute UMAP
    mu.tl.umap(mdata, neighbors_key='wnn', random_state=random_state, min_dist=min_dist)

    # Create temporary AnnData for clustering
    import anndata
    adata_cluster = anndata.AnnData(
        X=np.zeros((mdata.n_obs, 1)),
        obs=mdata.obs.copy()
    )

    # Copy WNN graph
    connectivities_key = 'wnn_connectivities' if 'wnn_connectivities' in mdata.obsp else 'connectivities'
    distances_key = 'wnn_distances' if 'wnn_distances' in mdata.obsp else 'distances'

    adata_cluster.obsp['connectivities'] = mdata.obsp[connectivities_key]
    if distances_key in mdata.obsp:
        adata_cluster.obsp['distances'] = mdata.obsp[distances_key]

    adata_cluster.uns['neighbors'] = {
        'connectivities_key': 'connectivities',
        'distances_key': 'distances',
        'params': {'n_neighbors': n_neighbors, 'method': 'wnn'}
    }

    # Leiden clustering
    sc.tl.leiden(adata_cluster, resolution=resolution, random_state=random_state)
    mdata.obs['wnn_leiden'] = adata_cluster.obs['leiden'].values

    return mdata, adata_cluster


# =============================================================================
# Visualization Functions
# =============================================================================

def plot_umap_categorical(
    embedding: np.ndarray,
    labels: np.ndarray,
    title: str = "UMAP",
    figsize: Tuple[int, int] = (8, 6),
    point_size: int = 30,
    alpha: float = 0.7
) -> plt.Figure:
    """
    Plot UMAP colored by categorical labels.

    Args:
        embedding: UMAP coordinates (n_cells, 2)
        labels: Category labels for each cell
        title: Plot title
        figsize: Figure size
        point_size: Size of scatter points
        alpha: Point transparency

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    unique_labels = np.unique(labels)
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_labels)))

    for i, label in enumerate(unique_labels):
        mask = labels == label
        ax.scatter(embedding[mask, 0], embedding[mask, 1],
                   c=[colors[i]], label=str(label), s=point_size, alpha=alpha)

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax.set_aspect('equal')

    plt.tight_layout()
    return fig


def plot_umap_continuous(
    embedding: np.ndarray,
    values: np.ndarray,
    title: str = "UMAP",
    cmap: str = 'viridis',
    figsize: Tuple[int, int] = (8, 6),
    point_size: int = 30,
    alpha: float = 0.7,
    colorbar_label: str = ""
) -> plt.Figure:
    """
    Plot UMAP colored by continuous values.

    Args:
        embedding: UMAP coordinates (n_cells, 2)
        values: Continuous values for coloring
        title: Plot title
        cmap: Colormap name
        figsize: Figure size
        point_size: Size of scatter points
        alpha: Point transparency
        colorbar_label: Label for colorbar

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(embedding[:, 0], embedding[:, 1],
                         c=values, cmap=cmap, s=point_size, alpha=alpha)
    plt.colorbar(scatter, ax=ax, label=colorbar_label)

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_title(title)
    ax.set_aspect('equal')

    plt.tight_layout()
    return fig


def plot_modality_weights_histogram(
    weights_df: pd.DataFrame,
    figsize: Tuple[int, int] = (10, 4)
) -> plt.Figure:
    """
    Plot histogram of modality weights.

    Args:
        weights_df: DataFrame with columns for each modality's weight
        figsize: Figure size

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    colors = {'WF': '#F39922', 'ISI': '#12A84B', 'PSTH': '#0B78BE'}

    for col in weights_df.columns:
        mod_name = col.replace('_weight', '').replace(':mod_weight', '')
        color = colors.get(mod_name, '#666666')
        ax.hist(weights_df[col], bins=20, alpha=0.6, color=color, label=mod_name)

    ax.set_xlabel('Weight')
    ax.set_ylabel('Count')
    ax.set_title('Modality Weight Distributions')
    ax.legend()

    plt.tight_layout()
    return fig


def plot_weights_by_celltype(
    weights_df: pd.DataFrame,
    cell_types: np.ndarray,
    figsize: Tuple[int, int] = (10, 5)
) -> plt.Figure:
    """
    Plot stacked bar chart of average weights by cell type.

    Args:
        weights_df: DataFrame with weight columns
        cell_types: Array of cell type labels
        figsize: Figure size

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Combine weights and cell types
    df = weights_df.copy()
    df['CellType'] = cell_types

    # Calculate mean weights per cell type
    weight_cols = [c for c in df.columns if 'weight' in c.lower()]
    mean_weights = df.groupby('CellType')[weight_cols].mean()

    # Rename columns for plotting
    rename_map = {}
    for col in weight_cols:
        for mod in ['WF', 'ISI', 'PSTH']:
            if mod in col:
                rename_map[col] = mod
                break
    mean_weights = mean_weights.rename(columns=rename_map)

    colors = ['#F39922', '#12A84B', '#0B78BE']
    mean_weights.plot(kind='bar', stacked=True, ax=ax, color=colors)

    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Mean Weight')
    ax.set_title('Average Modality Weights by Cell Type')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    plt.tight_layout()
    return fig


def plot_dominant_modality(
    embedding: np.ndarray,
    weights_df: pd.DataFrame,
    figsize: Tuple[int, int] = (8, 6),
    point_size: int = 30
) -> plt.Figure:
    """
    Plot UMAP colored by dominant modality for each cell.

    Args:
        embedding: UMAP coordinates
        weights_df: DataFrame with weight columns
        figsize: Figure size
        point_size: Size of scatter points

    Returns:
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    colors = {'WF': '#F39922', 'ISI': '#12A84B', 'PSTH': '#0B78BE'}

    # Find dominant modality
    weight_cols = list(weights_df.columns)
    dominant = weights_df.idxmax(axis=1)

    for col in weight_cols:
        # Extract modality name
        for mod in ['WF', 'ISI', 'PSTH']:
            if mod in col:
                mask = dominant == col
                ax.scatter(embedding[mask, 0], embedding[mask, 1],
                          c=colors[mod], label=mod, s=point_size, alpha=0.7)
                break

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_title('Dominant Modality per Cell')
    ax.legend()
    ax.set_aspect('equal')

    plt.tight_layout()
    return fig


def get_data_summary(
    modality_dict: Dict[str, pd.DataFrame],
    metadata_df: pd.DataFrame
) -> Dict:
    """
    Generate summary statistics for the loaded data.

    Args:
        modality_dict: Dictionary of modality DataFrames
        metadata_df: Metadata DataFrame

    Returns:
        Dictionary of summary statistics
    """
    summary = {
        'n_cells': len(metadata_df),
        'modalities': {},
        'cell_types': {},
        'metadata_columns': list(metadata_df.columns)
    }

    for name, df in modality_dict.items():
        summary['modalities'][name] = {
            'n_features': df.shape[1],
            'shape': df.shape
        }

    if 'CellType' in metadata_df.columns:
        summary['cell_types'] = metadata_df['CellType'].value_counts().to_dict()

    if 'Layer' in metadata_df.columns:
        summary['layers'] = metadata_df['Layer'].value_counts().to_dict()

    return summary
