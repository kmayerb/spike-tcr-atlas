"""
v1_Fig_S1.py  (formerly 01/ID159_catalog_figures_00_prelim_graphs.py)

Supplementary Figure S1: TCR-beta similarity networks for HLA-associated
TCRs, one panel for HLA class I (A/B/C) and one for HLA class II (DRB1).

  Nodes : unique (V-family + CDR3b, HLA allele) features
  Edges : approximate TCRdist(beta) <= MAX_DIST  AND  same HLA allele
  Color : HLA allele ('binary' column)
  Size  : number of OXF samples in which the TCR was HLA-detected (x10)

Outputs (in OUT_DIR):
  v1_Fig_S1_network_dr.pdf / .color_map.csv     class II network
  v1_Fig_S1_network_abc.pdf / .color_map.csv    class I network
  *_counts_*.csv                                per-allele node counts (for R)
"""

import os
import sys

import networkx as nx
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
TCR2HLA_PATH = "/fh/fast/gilbert_p/kmayerbl/TCR2HLA/"
TCRDISTGPU_PATH = "/fh/fast/gilbert_p/kmayerbl/TCRdist_GPU/"
DATA_DIR = "/fh/fast/gilbert_p/kmayerbl/uox_lambde_sap_2025/checkpoints_data"
BREADTH_CSV = os.path.join(DATA_DIR, "2025112_fold_change_query100_hla_breadth_raw.csv")
FOLD_CHANGE_CSV = os.path.join(DATA_DIR, "2025112_fold_change_query100.csv")

OUT_DIR = "Manuscript/figures"          # relative to the working directory
OUT_STEM_CLASS_II = "v1_Fig_S1_network_dr"
OUT_STEM_CLASS_I = "v1_Fig_S1_network_abc"

MAX_DIST = 25              # TCRdist(beta) edge threshold
MIN_COMPONENT_SIZE = 1     # 1 keeps everything (incl. singletons); 2 drops singletons
PLOT_SIZE_IN = 12          # width and height passed to draw_categories
NODE_SIZE_SCALE = 10       # node size = hla_detect * NODE_SIZE_SCALE
NODE_SIZE_CLIP = (1, 250)

CLASS_I_PREFIXES = ("A_", "B_", "C_")
CLASS_II_PREFIXES = ("DR",)

# Manual colour tweaks applied on top of draw_categories' automatic palette.
CLASS_I_COLOR_OVERRIDES = {"B_5101": "black", "C_0202": "#808080"}
CLASS_II_COLOR_OVERRIDES = {}   # e.g. fill from a fixed DRB1 palette if wanted

sys.path.append(TCR2HLA_PATH)
sys.path.append(TCRDISTGPU_PATH)
from tcrdistgpu.distance import TCRgpu          # freely available: kmayerb/tcrdistgpu  v1.0-manuscript,b54d004
from tcrdistgpu.draw import draw_categories     # freely available: kmayerb/tcr2hla 0b437d3b1ed9ecbca4be2655e821004d45dc811b

os.makedirs(OUT_DIR, exist_ok=True)

# V-family code (e.g. "V05") -> representative TRBV allele used for TCRdist.
# The exact gene within a family is unknown, so the first functional-looking
# member is used as a stand-in. "# P" marks families that are mostly pseudogenes.
TRBV_FAMILY_TO_GENE = {
    "V01": "TRBV1*01",     "V02": "TRBV2*01",     "V03": "TRBV3-1*01",
    "V04": "TRBV4-1*01",   "V05": "TRBV5-1*01",   "V06": "TRBV6-1*01",
    "V07": "TRBV7-2*01",   "V08": "TRBV8-1*01",   "V09": "TRBV9*01",     # P
    "V10": "TRBV10-1*01",  "V11": "TRBV11-1*01",  "V12": "TRBV12-3*01",
    "V13": "TRBV13*01",    "V14": "TRBV14*01",    "V15": "TRBV15*01",
    "V16": "TRBV16*01",    "V17": "TRBV17*01",    "V18": "TRBV18*01",
    "V19": "TRBV19*01",    "V20": "TRBV20-1*01",  "V21": "TRBV21-1*01",
    "V22": "TRBV22*01",    "V23": "TRBV23-1*01",  "V24": "TRBV24-1*01",  # V22: P
    "V25": "TRBV25-1*01",  "V26": "TRBV26*01",    "V27": "TRBV27*01",
    "V28": "TRBV28*01",    "V29": "TRBV29-1*01",  "V30": "TRBV30*01",
}

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def to_tribble(df):
    """Return R `tribble(...)` source code reproducing a small DataFrame."""
    def fmt(val):
        if pd.isna(val):
            return "NA"
        if isinstance(val, str):
            return '"' + val.replace('"', '\\"') + '"'
        if isinstance(val, (bool, np.bool_)):
            return str(val).upper()
        return str(val)

    lines = ["tribble(", "  " + ", ".join(f"~{c}" for c in df.columns) + ","]
    rows = ["  " + ", ".join(fmt(v) for v in row)
            for row in df.itertuples(index=False)]
    lines.append(",\n".join(rows))
    lines.append(")")
    return "\n".join(lines)


def is_valid_cdr3(cdr3, min_len=5):
    """True if cdr3 is a string of >= min_len standard amino acids."""
    return isinstance(cdr3, str) and len(cdr3) >= min_len and set(cdr3) <= STANDARD_AA


def compute_tcrdist_csr(df, chain="b", max_dist=14):
    """
    Compute a sparse approximate-TCRdist matrix among the rows of `df`.

    Rows with an unrecognised V gene or invalid CDR3 are dropped first.

    Returns
    -------
    csr : scipy.sparse matrix, shape (n_valid, n_valid)
    clean : DataFrame of the n_valid rows that were kept, re-indexed 0..n_valid-1.
        Row i of `clean` is node i of the matrix -- always use `clean`
        (not the input df) for anything downstream of the matrix.
    """
    if chain not in ("a", "b", "ab"):
        raise ValueError(f"Invalid chain option: {chain!r}")

    tg = TCRgpu(mode="cpu",          # set to "cuda" on a GPU node
                chunk_size=100,      # rows per chunk; limits RAM use
                cdr3a_col="cdr3_a_aa", cdr3b_col="cdr3_b_aa",
                va_col="v_a_gene", vb_col="v_b_gene")

    keep = pd.Series(True, index=df.index)
    if "a" in chain:
        keep &= df["v_a_gene"].isin(tg.params_vec.keys())
        keep &= df["cdr3_a_aa"].apply(is_valid_cdr3)
    if "b" in chain:
        keep &= df["v_b_gene"].isin(tg.params_vec.keys())
        keep &= df["cdr3_b_aa"].apply(is_valid_cdr3)

    clean = df[keep].reset_index(drop=True)
    print(f"\tTCRdist({max_dist}, {chain}): kept {len(clean)} / {len(df)} rows "
          f"({len(df) - len(clean)} dropped as invalid)")

    encode = {"a": tg.encode_tcrs_a, "b": tg.encode_tcrs_b, "ab": tg.encode_tcrs}[chain]
    enc = encode(clean)
    csr = tg.compute_csr_experimental(enc, enc, max_dist=max_dist)
    return csr, clean


def build_same_allele_graph(csr, data, group_col="binary", min_component_size=1):
    """
    Graph whose edges are TCRdist neighbours (from `csr`) that also share the
    same value of `group_col`. Node i corresponds to row i of `data`.
    """
    if csr.shape[0] != len(data):
        raise ValueError(f"csr has {csr.shape[0]} nodes but data has {len(data)} rows")

    G_dist = nx.from_scipy_sparse_array(csr)
    group = data[group_col].to_numpy()

    G = nx.Graph()
    G.add_nodes_from(range(len(data)))
    G.add_edges_from((u, v) for u, v in G_dist.edges() if group[u] == group[v])
    print(f"\tDistance graph: {G_dist}")
    print(f"\tSame-{group_col} graph: {G}")

    if min_component_size > 1:
        small = [n for comp in nx.connected_components(G)
                 if len(comp) < min_component_size for n in comp]
        G.remove_nodes_from(small)
        print(f"\tAfter removing components < {min_component_size}: {G}")
    return G


def plot_network(csr, data, savefig, color_overrides=None):
    """Build the same-allele network and draw it. Returns {allele: colour}."""
    G = build_same_allele_graph(csr, data, min_component_size=MIN_COMPONENT_SIZE)
    data = data.assign(node_size=data["hla_detect"] * NODE_SIZE_SCALE)
    key_to_color = draw_categories(G=G, data=data,
                                   color_col="binary",
                                   size_col="node_size",
                                   size_clip=NODE_SIZE_CLIP,
                                   key_to_color_override=color_overrides,
                                   legend=True,
                                   width=PLOT_SIZE_IN, height=PLOT_SIZE_IN,
                                   savefig=savefig)
    print(f"\tWrote {savefig}")
    return key_to_color


def save_color_map(key_to_color, path):
    pd.Series(key_to_color, name="color").rename_axis("binary").to_csv(path)


def export_allele_counts(data, stem):
    """Per-allele node counts at three hla_detect thresholds -> CSV + R tribble."""
    thresholds = {"all": data["hla_detect"] >= 0,
                  "detect_gt0": data["hla_detect"] > 0,
                  "detect_gt1": data["hla_detect"] > 1}
    for label, mask in thresholds.items():
        counts = (data.loc[mask, "binary"].value_counts()
                  .rename_axis("binary").reset_index(name="n"))
        counts.to_csv(os.path.join(OUT_DIR, f"{stem}_counts_{label}.csv"), index=False)
        print(f"# {stem}: {label}")
        print(to_tribble(counts))


# --------------------------------------------------------------------------
# 1. Load data and summarise per TCR feature across OXF samples
# --------------------------------------------------------------------------
# dfb = pd.read_csv('/fh/fast/gilbert_p/kmayerbl/uox_lambde_sap_2025/checkpoints_data/2025112_fold_change_query100_hla_breadth_raw.csv')
# df1 = pd.read_csv('/fh/fast/gilbert_p/kmayerbl/uox_lambde_sap_2025/checkpoints_data/2025112_fold_change_query100.csv')
# # MERGE QUERY TO 
# df1x = df1.merge(dfb[['name_1','index','vfamcdr3','binary','detect','hla_detect']], how = "left", on =['name_1','index','vfamcdr3','binary'] )
# df1x['detect'] = df1x['detect'].fillna(0)
# df1x['hla_detect'] = df1x['hla_detect'].fillna(0)
# # only OXF samples
# dfox = df1x[df1x.name_1.str.startswith('OXF')].reset_index(drop = True)
# # Summary across featurs how often to we get an HLA detect
# dfox.groupby(['index','vfamcdr3','binary'])[['detect','hla_detect']].mean()
# #df = df1.query('name_1 == "OXFMV251_TCRB"').reset_index(drop = True)
# df = dfox.groupby(['index','vfamcdr3','binary','name_1']).head(1).groupby(['index','vfamcdr3','binary'])[['detect','hla_detect']].sum().reset_index(drop = False)
# 



KEYS = ["index", "vfamcdr3", "binary"]

breadth = pd.read_csv(BREADTH_CSV)
fold_change = pd.read_csv(FOLD_CHANGE_CSV)

merged = fold_change.merge(
    breadth[["name_1"] + KEYS + ["detect", "hla_detect"]],
    how="left", on=["name_1"] + KEYS,
)
merged[["detect", "hla_detect"]] = merged[["detect", "hla_detect"]].fillna(0)

oxf = merged[merged["name_1"].str.startswith("OXF")]

# One row per (feature, sample), then sum across samples:
#   detect / hla_detect = number of OXF samples in which the feature was detected
features = (oxf.drop_duplicates(subset=KEYS + ["name_1"])
               .groupby(KEYS)[["detect", "hla_detect"]].sum()
               .reset_index())

features["v_b_gene"] = features["vfamcdr3"].str[:3].map(TRBV_FAMILY_TO_GENE)
features["cdr3_b_aa"] = features["vfamcdr3"].str[3:]

data_i = features[features["binary"].str.startswith(CLASS_I_PREFIXES)].reset_index(drop=True)
data_ii = features[features["binary"].str.startswith(CLASS_II_PREFIXES)].reset_index(drop=True)

# --------------------------------------------------------------------------
# 2. Class II (DRB1) network
# --------------------------------------------------------------------------
print("## Class II")
csr_ii, data_ii = compute_tcrdist_csr(data_ii, chain="b", max_dist=MAX_DIST)
colors_ii = plot_network(csr_ii, data_ii,
                         savefig=os.path.join(OUT_DIR, f"{OUT_STEM_CLASS_II}.pdf"),
                         color_overrides=CLASS_II_COLOR_OVERRIDES or None)
save_color_map(colors_ii, os.path.join(OUT_DIR, f"{OUT_STEM_CLASS_II}.color_map.csv"))
export_allele_counts(data_ii, OUT_STEM_CLASS_II)

# --------------------------------------------------------------------------
# 3. Class I (A/B/C) network
# --------------------------------------------------------------------------
# draw_categories assigns colours automatically; to tweak a few, we take its
# palette from a first pass and redraw with the overrides applied.
print("## Class I")
csr_i, data_i = compute_tcrdist_csr(data_i, chain="b", max_dist=MAX_DIST)
draft_path = os.path.join(OUT_DIR, "v1_Fig_S1_network_class_i_draft.pdf")
colors_i = plot_network(csr_i, data_i, savefig=draft_path)
colors_i.update(CLASS_I_COLOR_OVERRIDES)
colors_i = plot_network(csr_i, data_i,
                        savefig=os.path.join(OUT_DIR, f"{OUT_STEM_CLASS_I}.pdf"),
                        color_overrides=colors_i)
save_color_map(colors_i, os.path.join(OUT_DIR, f"{OUT_STEM_CLASS_I}.color_map.csv"))
export_allele_counts(data_i, OUT_STEM_CLASS_I)
