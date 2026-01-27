from pathlib import Path
import pandas as pd
from Bio import SeqIO
from protFeat.feature_extracter import extract_protein_feature
import sys

INPUT_FOLDER = Path(r"C:\Users\91830\Desktop\classfier_pfeature\protfeat\MAPH")
OUTPUT_FOLDER = Path(r"C:\Users\91830\Desktop\classfier_pfeature\protfeat\feature_extraction_output")
FASTA_NAME = "MAPH_permissive_noPfam_250AA"
FEATURES = ["CTriad"]
OUTPUT_FILE = "MAPH_CTriad_feat.csv"

FEATURE_DIM_LIB = {
    "AAC": 20, "CKSAAP": 2400, "DPC": 400, "GAAC": 5, "CKSAAGP": 150,
    "GDPC": 25, "Moran": 240, "Geary": 240, "NMBroto": 240,
    "CTDC": 39, "CTDT": 39, "CTDD": 195, "CTriad": 343,
    "SOCNumber": 60, "QSOrder": 100, "PAAC": 50, "APAAC": 80
}

def ensure_fasta_extension(fasta_name):
    return fasta_name if fasta_name.endswith(".fasta") else fasta_name + ".fasta"

def extract_features():
    for feature in FEATURES:
        output_txt = OUTPUT_FOLDER / f"{FASTA_NAME.replace('.fasta','')}_{feature}.txt"
        if output_txt.exists():
            continue
        extract_protein_feature(feature, 0, str(INPUT_FOLDER), FASTA_NAME.replace(".fasta", ""))

def load_protein_ids(fasta_path):
    if not fasta_path.exists():
        sys.exit(1)
    ids = [record.id for record in SeqIO.parse(str(fasta_path), "fasta")]
    print(f"FASTA sequences: {len(ids)}")
    return ids

def get_clean_feature_df(feature, protein_ids):
    txt_path = OUTPUT_FOLDER / f"{FASTA_NAME.replace('.fasta','')}_{feature}.txt"
    df = pd.read_csv(txt_path, sep="\t", header=None)
    if df[0].dtype == 'object' or isinstance(df.iloc[0, 0], str):
        try:
            pd.to_numeric(df[0])
            if df.shape[1] == FEATURE_DIM_LIB.get(feature, 0) + 1:
                df = df.drop(columns=[0])
        except ValueError:
            df = df.drop(columns=[0])
    df = df.dropna(axis=1, how='all')
    expected_dim = FEATURE_DIM_LIB.get(feature)
    if expected_dim:
        df = df.iloc[:, :expected_dim]
    df.columns = [f"{feature}_{i+1}" for i in range(df.shape[1])]
    df.insert(0, "IDs", protein_ids[:len(df)])
    print(f"{feature}: {df.shape[0]} rows x {df.shape[1]} columns")
    return df

def main():
    global FASTA_NAME
    FASTA_NAME = ensure_fasta_extension(FASTA_NAME)
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
    extract_features()
    fasta_path = INPUT_FOLDER / FASTA_NAME
    protein_ids = load_protein_ids(fasta_path)
    final_df = pd.DataFrame({"IDs": protein_ids})
    for feature in FEATURES:
        df = get_clean_feature_df(feature, protein_ids)
        final_df = pd.merge(final_df, df, on="IDs", how="left")
    output_path = OUTPUT_FOLDER / OUTPUT_FILE
    final_df.to_csv(output_path, index=False, header=True)  # HEADER INCLUDED
    print(f"Final dataset: {final_df.shape[0]} rows x {final_df.shape[1]} columns")
    print(f"Saved to: {output_path}")

if __name__ == "__main__":
    main()
