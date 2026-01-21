"""
Generate example CSV files from the Jianing MATLAB data.
Run this script to create example_data/ CSVs for testing the Streamlit app.
"""

import numpy as np
import pandas as pd
from scipy.io import loadmat
import os

def extract_matlab_strings(arr):
    """Extract strings from MATLAB cell arrays."""
    arr = np.asarray(arr)
    if arr.ndim > 1:
        arr = arr.flatten()

    result = []
    for item in arr:
        if isinstance(item, np.ndarray):
            if item.size > 0:
                val = item.flat[0]
                result.append(str(val).strip())
            else:
                result.append('')
        elif isinstance(item, str):
            result.append(item.strip())
        else:
            result.append(str(item).strip())

    return np.array(result)

def main():
    # Paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    mat_file = os.path.join(script_dir, 'JianingData', 'MergedData.mat')
    output_dir = os.path.join(script_dir, 'example_data')

    # Load MATLAB data
    print(f"Loading {mat_file}...")
    mat_data = loadmat(mat_file, squeeze_me=False)

    # Get dimensions
    depth_arr = np.asarray(mat_data['depth']).flatten()
    n_cells = len(depth_arr)
    cell_ids = [f"cell_{i}" for i in range(n_cells)]

    # Extract matrices
    ISI = np.asarray(mat_data['ISI']).astype(np.float64)
    WF = np.asarray(mat_data['waveformV']).astype(np.float64)
    PSTH = np.asarray(mat_data['AllResp']).astype(np.float64)

    # Ensure correct orientation (cells x features)
    if ISI.shape[0] != n_cells:
        ISI = ISI.T
    if WF.shape[0] != n_cells:
        WF = WF.T
    if PSTH.shape[0] != n_cells:
        PSTH = PSTH.T

    # Extract metadata
    depth = depth_arr.astype(np.float64)
    latency = np.asarray(mat_data['latency']).flatten().astype(np.float64)
    latency = np.nan_to_num(latency, nan=0.0)
    onset_rate = np.asarray(mat_data['onset_rate']).flatten().astype(np.float64)
    onset_rate = np.nan_to_num(onset_rate, nan=0.0)

    layer = extract_matlab_strings(mat_data['layer'])
    cell_type = extract_matlab_strings(mat_data['true_cell_type'])

    # Try to load width and ratio_p2t from .Rda files
    try:
        import pyreadr
        width_path = os.path.join(script_dir, 'JianingData', 'width.Rda')
        ratio_path = os.path.join(script_dir, 'JianingData', 'ratio_p2t.Rda')

        width_result = pyreadr.read_r(width_path)
        width = list(width_result.values())[0].values.flatten()

        ratio_result = pyreadr.read_r(ratio_path)
        ratio_p2t = list(ratio_result.values())[0].values.flatten()
    except:
        print("Warning: Could not load .Rda files, using zeros for width/ratio_p2t")
        width = np.zeros(n_cells)
        ratio_p2t = np.zeros(n_cells)

    # Create DataFrames
    print("Creating CSV files...")

    # Waveforms
    wf_df = pd.DataFrame(WF, index=cell_ids,
                         columns=[f"WF_{i}" for i in range(WF.shape[1])])
    wf_df.index.name = 'cell_id'
    wf_df.to_csv(os.path.join(output_dir, 'waveforms.csv'))
    print(f"  waveforms.csv: {wf_df.shape}")

    # ISI
    isi_df = pd.DataFrame(ISI, index=cell_ids,
                          columns=[f"ISI_{i}" for i in range(ISI.shape[1])])
    isi_df.index.name = 'cell_id'
    isi_df.to_csv(os.path.join(output_dir, 'isi.csv'))
    print(f"  isi.csv: {isi_df.shape}")

    # PSTH
    psth_df = pd.DataFrame(PSTH, index=cell_ids,
                           columns=[f"PSTH_{i}" for i in range(PSTH.shape[1])])
    psth_df.index.name = 'cell_id'
    psth_df.to_csv(os.path.join(output_dir, 'psth.csv'))
    print(f"  psth.csv: {psth_df.shape}")

    # Metadata
    meta_df = pd.DataFrame({
        'CellType': cell_type,
        'Layer': layer,
        'Depth': depth,
        'latency': latency,
        'onset': onset_rate,
        'width': width,
        'ratio_p2t': ratio_p2t
    }, index=cell_ids)
    meta_df.index.name = 'cell_id'
    meta_df.to_csv(os.path.join(output_dir, 'metadata.csv'))
    print(f"  metadata.csv: {meta_df.shape}")

    print(f"\nExample data saved to {output_dir}/")
    print(f"Total cells: {n_cells}")
    print(f"Cell types: {np.unique(cell_type)}")
    print(f"Layers: {np.unique(layer)}")

if __name__ == '__main__':
    main()
