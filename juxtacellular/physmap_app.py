"""
PhysMAP: Multi-Modal Neural Data Explorer
Interactive Streamlit application for analyzing multi-modal electrophysiology data.
"""

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO

# Import utility functions
from physmap_utils import (
    load_csv_data,
    create_anndata_from_df,
    process_modality,
    run_wnn_integration,
    plot_umap_categorical,
    plot_umap_continuous,
    plot_modality_weights_histogram,
    plot_weights_by_celltype,
    plot_dominant_modality,
    get_data_summary
)

# =============================================================================
# Page Configuration
# =============================================================================

st.set_page_config(
    page_title="PhysMAP",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("PhysMAP: Multi-Modal Neural Data Explorer")
st.markdown("*Analyze multi-modal electrophysiology data with unimodal and multimodal integration*")

# =============================================================================
# Session State Initialization
# =============================================================================

if 'processed' not in st.session_state:
    st.session_state.processed = False
    st.session_state.modality_dict = None
    st.session_state.metadata_df = None
    st.session_state.adata_wf = None
    st.session_state.adata_isi = None
    st.session_state.adata_psth = None
    st.session_state.mdata_wnn = None
    st.session_state.weights_df = None

# =============================================================================
# Sidebar - Data Upload and Parameters
# =============================================================================

with st.sidebar:
    st.header("Data Upload")

    wf_file = st.file_uploader("Waveforms (WF)", type=['csv', 'xlsx'],
                                help="CSV with cells as rows, WF features as columns")
    isi_file = st.file_uploader("ISI Data", type=['csv', 'xlsx'],
                                 help="CSV with cells as rows, ISI features as columns")
    psth_file = st.file_uploader("PSTH Data", type=['csv', 'xlsx'],
                                  help="CSV with cells as rows, PSTH features as columns")
    meta_file = st.file_uploader("Metadata", type=['csv', 'xlsx'],
                                  help="CSV with CellType, Layer, and other metadata")

    st.divider()
    st.header("Processing Parameters")

    n_pcs = st.slider("Number of PCs", 10, 50, 30,
                      help="Number of principal components for dimensionality reduction")
    n_neighbors = st.slider("N Neighbors", 5, 50, 20,
                            help="Number of neighbors for kNN graph construction")
    min_dist = st.slider("UMAP Min Distance", 0.0, 1.0, 0.2, 0.05,
                         help="Minimum distance between points in UMAP")
    metric = st.selectbox("Distance Metric", ["cosine", "euclidean", "manhattan"],
                          help="Metric for neighbor calculation")
    resolution = st.slider("Clustering Resolution", 0.1, 3.0, 2.0, 0.1,
                           help="Resolution parameter for Leiden clustering")
    random_seed = st.number_input("Random Seed", value=42, min_value=0,
                                   help="Seed for reproducibility")

    st.divider()
    st.header("Normalization")

    normalize_wf = st.checkbox("Normalize WF (CLR)", value=False,
                               help="Apply CLR normalization to waveforms")
    normalize_isi = st.checkbox("Normalize ISI (CLR)", value=True,
                                help="Apply CLR normalization to ISI")
    normalize_psth = st.checkbox("Normalize PSTH (CLR)", value=True,
                                 help="Apply CLR normalization to PSTH")

    st.divider()
    st.header("Visualization")

    point_size = st.slider("Point Size", 10, 100, 30)

    st.divider()

    # Run button
    run_button = st.button("Run Analysis", type="primary", use_container_width=True)

# =============================================================================
# Main Content
# =============================================================================

# Check if all files are uploaded
all_files_uploaded = all([wf_file, isi_file, psth_file, meta_file])

if not all_files_uploaded:
    st.info("Please upload all four data files (WF, ISI, PSTH, Metadata) to begin analysis.")

    # Show expected format
    with st.expander("Expected CSV Format"):
        st.markdown("""
        **waveforms.csv** (cells x features):
        ```
        cell_id,WF_0,WF_1,WF_2,...
        cell_0,0.123,0.456,...
        cell_1,0.234,0.567,...
        ```

        **isi.csv** (cells x features):
        ```
        cell_id,ISI_0,ISI_1,...
        cell_0,0.1,0.2,...
        ```

        **psth.csv** (cells x features):
        ```
        cell_id,PSTH_0,PSTH_1,...
        cell_0,0.5,0.3,...
        ```

        **metadata.csv** (cells x metadata):
        ```
        cell_id,CellType,Layer,Depth,latency,width,ratio_p2t
        cell_0,E,4,500,10.5,0.3,1.2
        ```
        """)

# Run analysis when button is clicked
if run_button and all_files_uploaded:
    with st.spinner("Loading data..."):
        try:
            modality_dict, metadata_df = load_csv_data(wf_file, isi_file, psth_file, meta_file)
            st.session_state.modality_dict = modality_dict
            st.session_state.metadata_df = metadata_df
        except Exception as e:
            st.error(f"Error loading data: {e}")
            st.stop()

    # Create layerCellType if both columns exist
    if 'CellType' in metadata_df.columns and 'Layer' in metadata_df.columns:
        metadata_df['layerCellType'] = metadata_df['CellType'].astype(str) + '-' + metadata_df['Layer'].astype(str)

    progress_bar = st.progress(0)
    status_text = st.empty()

    # Process WF
    status_text.text("Processing WF modality...")
    adata_wf = create_anndata_from_df(modality_dict['WF'], metadata_df, 'WF')
    adata_wf = process_modality(
        adata_wf, 'WF', n_pcs=n_pcs, normalize=normalize_wf,
        n_neighbors=n_neighbors, metric=metric, min_dist=min_dist,
        resolution=resolution, random_state=random_seed
    )
    st.session_state.adata_wf = adata_wf
    progress_bar.progress(25)

    # Process ISI
    status_text.text("Processing ISI modality...")
    adata_isi = create_anndata_from_df(modality_dict['ISI'], metadata_df, 'ISI')
    adata_isi = process_modality(
        adata_isi, 'ISI', n_pcs=n_pcs, normalize=normalize_isi,
        n_neighbors=n_neighbors, metric=metric, min_dist=min_dist,
        resolution=resolution, random_state=random_seed
    )
    st.session_state.adata_isi = adata_isi
    progress_bar.progress(50)

    # Process PSTH
    status_text.text("Processing PSTH modality...")
    adata_psth = create_anndata_from_df(modality_dict['PSTH'], metadata_df, 'PSTH')
    adata_psth = process_modality(
        adata_psth, 'PSTH', n_pcs=n_pcs, normalize=normalize_psth,
        n_neighbors=n_neighbors, metric=metric, min_dist=min_dist,
        resolution=resolution, random_state=random_seed
    )
    st.session_state.adata_psth = adata_psth
    progress_bar.progress(75)

    # Run WNN integration
    status_text.text("Running multi-modal integration (WNN)...")
    adata_dict = {'WF': adata_wf, 'ISI': adata_isi, 'PSTH': adata_psth}
    mdata_wnn, _ = run_wnn_integration(
        adata_dict, n_neighbors=n_neighbors,
        min_dist=min_dist, resolution=resolution, random_state=random_seed
    )
    st.session_state.mdata_wnn = mdata_wnn

    # Extract modality weights
    weight_cols = [c for c in mdata_wnn.obs.columns if 'weight' in c.lower()]
    if weight_cols:
        st.session_state.weights_df = mdata_wnn.obs[weight_cols].copy()

    progress_bar.progress(100)
    status_text.text("Analysis complete!")

    st.session_state.processed = True
    st.rerun()

# =============================================================================
# Display Results
# =============================================================================

if st.session_state.processed:
    # Retrieve from session state
    adata_wf = st.session_state.adata_wf
    adata_isi = st.session_state.adata_isi
    adata_psth = st.session_state.adata_psth
    mdata_wnn = st.session_state.mdata_wnn
    metadata_df = st.session_state.metadata_df
    weights_df = st.session_state.weights_df

    # Get available color options
    color_options = ['leiden']
    if 'CellType' in metadata_df.columns:
        color_options.append('CellType')
    if 'Layer' in metadata_df.columns:
        color_options.append('Layer')
    if 'layerCellType' in metadata_df.columns:
        color_options.append('layerCellType')

    # Create tabs
    tab1, tab2, tab3 = st.tabs([
        "Individual Modalities",
        "Multi-Modal Integration",
        "Data Summary"
    ])

    # =========================================================================
    # Tab 1: Individual Modalities
    # =========================================================================
    with tab1:
        st.header("Individual Modality UMAPs")

        color_by = st.selectbox("Color by:", color_options, key='tab1_color')

        col1, col2, col3 = st.columns(3)

        with col1:
            st.subheader("WF")
            if color_by == 'leiden':
                labels = adata_wf.obs['WF_leiden'].values
            else:
                labels = adata_wf.obs[color_by].values
            fig_wf = plot_umap_categorical(
                adata_wf.obsm['X_umap'], labels,
                title=f"WF - {color_by}", point_size=point_size
            )
            st.pyplot(fig_wf)
            plt.close()

        with col2:
            st.subheader("ISI")
            if color_by == 'leiden':
                labels = adata_isi.obs['ISI_leiden'].values
            else:
                labels = adata_isi.obs[color_by].values
            fig_isi = plot_umap_categorical(
                adata_isi.obsm['X_umap'], labels,
                title=f"ISI - {color_by}", point_size=point_size
            )
            st.pyplot(fig_isi)
            plt.close()

        with col3:
            st.subheader("PSTH")
            if color_by == 'leiden':
                labels = adata_psth.obs['PSTH_leiden'].values
            else:
                labels = adata_psth.obs[color_by].values
            fig_psth = plot_umap_categorical(
                adata_psth.obsm['X_umap'], labels,
                title=f"PSTH - {color_by}", point_size=point_size
            )
            st.pyplot(fig_psth)
            plt.close()

        # Cluster counts
        st.subheader("Cluster Counts")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.write("**WF Clusters:**")
            st.write(adata_wf.obs['WF_leiden'].value_counts())
        with col2:
            st.write("**ISI Clusters:**")
            st.write(adata_isi.obs['ISI_leiden'].value_counts())
        with col3:
            st.write("**PSTH Clusters:**")
            st.write(adata_psth.obs['PSTH_leiden'].value_counts())

    # =========================================================================
    # Tab 2: Multi-Modal Integration
    # =========================================================================
    with tab2:
        st.header("Weighted Nearest Neighbors (WNN) Integration")

        # WNN UMAP
        col1, col2 = st.columns([2, 1])

        with col1:
            color_by_wnn = st.selectbox("Color by:", color_options + ['wnn_leiden'], key='tab2_color')

            if color_by_wnn == 'wnn_leiden':
                labels = mdata_wnn.obs['wnn_leiden'].values
            elif color_by_wnn == 'leiden':
                labels = mdata_wnn.obs.get('wnn_leiden', mdata_wnn.mod['WF'].obs['WF_leiden']).values
            else:
                labels = mdata_wnn.mod['WF'].obs[color_by_wnn].values

            fig_wnn = plot_umap_categorical(
                mdata_wnn.obsm['X_umap'], labels,
                title=f"WNN UMAP - {color_by_wnn}",
                figsize=(10, 8), point_size=point_size
            )
            st.pyplot(fig_wnn)
            plt.close()

        with col2:
            st.subheader("WNN Clusters")
            st.write(mdata_wnn.obs['wnn_leiden'].value_counts())

        # Modality weights
        if weights_df is not None and len(weights_df.columns) > 0:
            st.subheader("Modality Weights")

            col1, col2 = st.columns(2)

            with col1:
                fig_hist = plot_modality_weights_histogram(weights_df)
                st.pyplot(fig_hist)
                plt.close()

                # Weight statistics
                st.write("**Weight Statistics:**")
                stats_df = weights_df.describe().T[['mean', 'std', 'min', 'max']]
                stats_df.index = [idx.replace(':mod_weight', '').replace('_weight', '')
                                  for idx in stats_df.index]
                st.dataframe(stats_df.round(3))

            with col2:
                fig_dominant = plot_dominant_modality(
                    mdata_wnn.obsm['X_umap'], weights_df, point_size=point_size
                )
                st.pyplot(fig_dominant)
                plt.close()

            # Weights by cell type
            if 'CellType' in metadata_df.columns:
                st.subheader("Weights by Cell Type")
                fig_by_type = plot_weights_by_celltype(
                    weights_df, mdata_wnn.mod['WF'].obs['CellType'].values
                )
                st.pyplot(fig_by_type)
                plt.close()

    # =========================================================================
    # Tab 3: Data Summary
    # =========================================================================
    with tab3:
        st.header("Data Summary")

        modality_dict = st.session_state.modality_dict
        summary = get_data_summary(modality_dict, metadata_df)

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Dataset Overview")
            st.metric("Total Cells", summary['n_cells'])

            st.write("**Modality Dimensions:**")
            for mod_name, mod_info in summary['modalities'].items():
                st.write(f"- {mod_name}: {mod_info['n_features']} features")

            st.write("**Metadata Columns:**")
            st.write(", ".join(summary['metadata_columns']))

        with col2:
            if 'cell_types' in summary and summary['cell_types']:
                st.subheader("Cell Type Distribution")
                ct_df = pd.DataFrame.from_dict(summary['cell_types'], orient='index', columns=['Count'])
                st.bar_chart(ct_df)

            if 'layers' in summary and summary['layers']:
                st.subheader("Layer Distribution")
                layer_df = pd.DataFrame.from_dict(summary['layers'], orient='index', columns=['Count'])
                st.bar_chart(layer_df)

        # Download buttons
        st.subheader("Download Results")

        col1, col2, col3 = st.columns(3)

        with col1:
            # Download WNN UMAP coordinates
            umap_df = pd.DataFrame(
                mdata_wnn.obsm['X_umap'],
                columns=['UMAP1', 'UMAP2'],
                index=mdata_wnn.obs_names
            )
            umap_df['wnn_leiden'] = mdata_wnn.obs['wnn_leiden'].values
            csv = umap_df.to_csv()
            st.download_button(
                "Download WNN UMAP",
                csv,
                "wnn_umap_coordinates.csv",
                "text/csv"
            )

        with col2:
            # Download modality weights
            if weights_df is not None:
                csv = weights_df.to_csv()
                st.download_button(
                    "Download Modality Weights",
                    csv,
                    "modality_weights.csv",
                    "text/csv"
                )

        with col3:
            # Download cluster assignments
            clusters_df = pd.DataFrame({
                'WF_cluster': adata_wf.obs['WF_leiden'].values,
                'ISI_cluster': adata_isi.obs['ISI_leiden'].values,
                'PSTH_cluster': adata_psth.obs['PSTH_leiden'].values,
                'WNN_cluster': mdata_wnn.obs['wnn_leiden'].values
            }, index=metadata_df.index)
            csv = clusters_df.to_csv()
            st.download_button(
                "Download Cluster Assignments",
                csv,
                "cluster_assignments.csv",
                "text/csv"
            )

        # Processing parameters used
        st.subheader("Processing Parameters")
        params_df = pd.DataFrame({
            'Parameter': ['N PCs', 'N Neighbors', 'Min Distance', 'Metric', 'Resolution', 'Random Seed'],
            'Value': [n_pcs, n_neighbors, min_dist, metric, resolution, random_seed]
        })
        st.dataframe(params_df, hide_index=True)

# =============================================================================
# Footer
# =============================================================================

st.divider()
st.markdown(
    "<div style='text-align: center; color: gray;'>"
    "PhysMAP | Multi-Modal Neural Data Analysis"
    "</div>",
    unsafe_allow_html=True
)
