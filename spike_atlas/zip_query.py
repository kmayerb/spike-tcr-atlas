import os
import zipfile
import functools
from typing import Callable, Optional
import pandas as pd
import parmap
import numpy as np
import pandas as pd
DEFAULT_RANDOM_STATE = 103083

def get_valid_names(zip_path: str) -> list[str]:
    """Return data-file member names, stripping Mac metadata entries."""
    with zipfile.ZipFile(zip_path) as zf:
        return [
            f for f in zf.namelist()
            if f.endswith(('.tsv', '.csv'))
            and not f.startswith('__MACOSX')
            and not os.path.basename(f).startswith('._')
        ]

def open_df_from_zip(zip_path: str, filename: str) -> pd.DataFrame:
    """
    Open a single CSV or TSV member from a zip archive without extracting it.
    Separator is inferred from the file extension.
    """
    ext = os.path.splitext(filename)[-1].lower()
    if ext == '.tsv':
        sep = '\t'
    elif ext == '.csv':
        sep = ','
    else:
        raise ValueError(
            f"Unsupported file extension {ext!r} for {filename!r}. "
            "Expected '.csv' or '.tsv'."
        )

    with zipfile.ZipFile(zip_path) as zf:
        with zf.open(filename) as fh:
            return pd.read_csv(fh, sep=sep, low_memory=False)

# ─────────────────────────────────────────────────────────────────────────────
#  PARSE & QUERY STAGES  (each independently swappable)
# ─────────────────────────────────────────────────────────────────────────────
def apply_parse(
    parse_func: Callable[[pd.DataFrame, str], pd.DataFrame],
    df: pd.DataFrame,
    filename: str,
    **kwargs,
) -> pd.DataFrame:
    """
    Apply an arbitrary parse function to a raw DataFrame.

    The parse_func signature must be:
        parse_func(df: pd.DataFrame, filename: str, **kwargs) -> pd.DataFrame
    """
    return parse_func(df=df, f=filename, **kwargs)

def apply_query(
    query_func: Callable[..., pd.DataFrame],
    df: pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
    """
    Apply an arbitrary query/filter/merge function to a parsed DataFrame.

    The query_func signature must be:
        query_func(df: pd.DataFrame, **kwargs) -> pd.DataFrame
    """
    return query_func(df=df, **kwargs)


#
# Deafult Parser
#
def parse_adapative_v2_gt1(df,f, gt1= False, out_cols = ['v','vfamcdr3','templates','frequency','productive_frequency','v_family','v_gene','j_gene','v_resolved','amino_acid','rearrangement','total_rearrangements','sample_id']):
    cols = df.columns
    # THIS IS FOR "V2" ADAPTIVE COLUMNS
    if 'aminoAcid' in cols:
        df['amino_acid'] = df['aminoAcid']
    if 'count (templates/reads)' in cols:
        df['templates'] = df['count (templates/reads)']
    if 'nucleotide' in cols:
        df['rearrangement'] = df['nucleotide']
    if 'frequencyCount (%)' in cols:
        df['frequency'] = df['frequencyCount (%)']
    if 'sequenceStatus' in cols:
        df['frame_type'] = df['sequenceStatus']
    if 'vFamilyName' in cols:
        df['v_family'] = df['vFamilyName']
    if 'vGeneName' in cols:
        df['v_gene'] = df['vGeneName']
    if 'jGeneName' in cols:
        df['j_gene'] = df['jGeneName']
    if 'vMaxResolved' in cols:
        df['v_resolved'] = df['vMaxResolved']
    if 'jMaxResolved' in cols:
        df['j_resolved'] = df['jMaxResolved']


    if 'v_resolved' not in cols:
        df['v_resolved'] = df['v_gene'].replace('unresolved', np.nan)
        df['v_resolved'] = df['v_resolved'].fillna(df['v_family'])
    
    if 'frequency' not in cols:
        df['frequency'] = np.nan
    if 'rearrangement' not in cols:
        df['rearrangement'] = np.nan

    for x in ['amino_acid', 'templates',  'v_gene', 'j_gene', 'v_family', 'v_resolved','frequency']:
       assert x in df.columns, f"{x} is missing"

    df['v'] = df['v_resolved'].apply(lambda x: x.split("*")[0].split("-")[0].replace("TCRB","") if isinstance(x,str) else np.nan)
    #df['v'] = df['v_family'].str.replace("TCRB","")
    if 'frame_type' in cols:
        df = df.query('frame_type == "In"').reset_index()
    
    df['productive_frequency'] = df['templates'] / df['templates'].sum()
    
    ix1 = df['amino_acid'].str.len() >= 6
    ix2 = df['amino_acid'].notna()
    df= df[ix1&ix2].reset_index(drop = True)
    ix3 = df['amino_acid'].apply(lambda x: x.find("*") == -1 )
    df= df[ix3].reset_index(drop = True)


    df = df[df['v'] != "unresolved"].reset_index(drop = True)
    df = df[df['v'] != "unknown"].reset_index(drop = True)
    df = df[df['v'] != np.nan].reset_index(drop = True)
    df['V'] = df['v']

    df['total_rearrangements'] = df.shape[0]

    if gt1:
        df= df[df['templates'] > 1].reset_index(drop = True)

    df['sample_id'] = os.path.basename(f).replace(".tsv","").replace('.csv','')
    df['vfamcdr3'] = df['v']+df['amino_acid']
    #print(os.path.basename(f).replace(".tsv","").replace('.csv',''))

    dfout = df[ out_cols].sort_values('templates', ascending = False)
    return(dfout)

def parse_v1x_gt1(df,f, gt1= False, out_cols = ['v','vfamcdr3','templates','frequency','productive_frequency','v_family','v_gene','j_gene','v_resolved','amino_acid','rearrangement','total_rearrangements','sample_id']):
    cols = df.columns
    # THIS IS FOR "V2" ADAPTIVE COLUMNS
    if 'aminoAcid' in cols:
        df['amino_acid'] = df['aminoAcid']
    if 'count (templates/reads)' in cols:
        df['templates'] = df['count (templates/reads)']
    if 'nucleotide' in cols:
        df['rearrangement'] = df['nucleotide']
    else:
        df['nucleotide'] = df['rearrangement']
    if 'frequencyCount (%)' in cols:
        df['frequency'] = df['frequencyCount (%)']
    if 'sequenceStatus' in cols:
        df['frame_type'] = df['sequenceStatus']
    if 'vFamilyName' in cols:
        df['v_family'] = df['vFamilyName']
    if 'vGeneName' in cols:
        df['v_gene'] = df['vGeneName']
    if 'jGeneName' in cols:
        df['j_gene'] = df['jGeneName']
    if 'vMaxResolved' in cols:
        df['v_resolved'] = df['vMaxResolved']
    if 'jMaxResolved' in cols:
        df['j_resolved'] = df['jMaxResolved']

    # DEAL WITH MISSING TEMPLATES IN PART OF EMERSON
    if df.templates.isna().mean() == 1.0:
       df['templates'] = df['seq_reads']


    if 'v_resolved' not in cols:
        df['v_resolved'] = df['v_gene'].replace('unresolved', np.nan)
        df['v_resolved'] = df['v_resolved'].fillna(df['v_family'])

    if 'frequency' not in cols:
        df['frequency'] = np.nan
    if 'rearrangement' not in cols:
        df['rearrangement'] = np.nan

    for x in ['amino_acid', 'templates',  'v_gene', 'j_gene', 'v_family', 'v_resolved','frequency']:
       assert x in df.columns, f"{x} is missing"

    df['v'] = df['v_resolved'].apply(lambda x: x.split("*")[0].split("-")[0].replace("TCRB","") if isinstance(x,str) else np.nan)
    #df['v'] = df['v_family'].str.replace("TCRB","")
    if 'frame_type' in cols:
        df = df.query('frame_type == "In"').reset_index()

    df['productive_frequency'] = df['templates'] / df['templates'].sum()

    ix1 = df['amino_acid'].str.len() >= 6
    ix2 = df['amino_acid'].notna()
    df= df[ix1&ix2].reset_index(drop = True)
    ix3 = df['amino_acid'].apply(lambda x: x.find("*") == -1 )
    df= df[ix3].reset_index(drop = True)


    df = df[df['v'] != "unresolved"].reset_index(drop = True)
    df = df[df['v'] != "unknown"].reset_index(drop = True)
    df = df[df['v'] != np.nan].reset_index(drop = True)
    df['V'] = df['v']

    df['total_rearrangements'] = df.shape[0]

    if gt1:
        df= df[df['templates'] > 1].reset_index(drop = True)

    df['sample_id'] = os.path.basename(f).replace(".tsv","").replace('.csv','')
    df['vfamcdr3'] = df['v']+df['amino_acid']
    #print(os.path.basename(f).replace(".tsv","").replace('.csv',''))

    dfout = df[ out_cols].sort_values('templates', ascending = False)
    return(dfout)


from functools import partial
parse_v2x_gt1_vfamcdr3_dna =  partial(parse_adapative_v2_gt1, 
     gt1= False, 
     out_cols = ['vfamcdr3','templates','productive_frequency','V','nucleotide','v_resolved','j_resolved'])
parse_v1x_gt1_vfamcdr3_dna =  partial(parse_v1x_gt1, 
     gt1= False, 
     out_cols = ['vfamcdr3','templates','productive_frequency','V','nucleotide','v_resolved','j_resolved'])






from typing import Optional
def fix_depth(d, target = 100000, rs =DEFAULT_RANDOM_STATE ):
    """
    The idea here is to fix depth
    """
    d['i'] = d.index
    total_starting = d['templates'].sum()
    #if total_starting < target:
        #print(f"{total_starting} < {target}")
    try:
        tmp= d.sample(n = target, 
        weights = d['productive_frequency'],
        random_state = rs, 
        replace=True)[['i','vfamcdr3', 'nucleotide','v_resolved', 'j_resolved']]
        #print("SUCCESSFUL DOWNSAMPLE")
    except ValueError:
        print("VALUE ERROR")
        print(d.head())
        print(d['productive_frequency'].sort_values().head())
        print(d['productive_frequency'].sort_values().tail())
        print("ADDING MIN")
        d['productive_frequency'] = d['productive_frequency'].apply(lambda x : max(1E-6,x))
        tmp= d.sample(n = target, 
        weights = d['productive_frequency'],
        random_state = rs, 
        replace=True)[['i','vfamcdr3', 'nucleotide','v_resolved', 'j_resolved']]
        #print("SUCCESSFUL DOWNSAMPLE AFTER FIXING WEIGHTS")
    tmp['templates'] =1
    ds = tmp.groupby(['vfamcdr3', 'nucleotide','v_resolved', 'j_resolved']).\
        sum().\
        reset_index().\
        sort_values('templates', ascending = False)
    ds['productive_frequency'] = ds['templates'] / ds['templates'].sum()
        
    return ds

# ─────────────────────────────────────────────────────────────────────────────
#  MERGE / QUERY IMPLEMENTATION  (the default query_func)
# ─────────────────────────────────────────────────────────────────────────────
def merge_query(
    df: pd.DataFrame,
    filename: str,                          # passed through from pipeline
    query: Optional[pd.DataFrame] = None,
    depthsample: Optional[int] = None,
    random_state: int = 103083,
    clone_col: str = "nucleotide",
    group_col: str = "vfamcdr3",
    freq_col: str = "productive_frequency",
    templates_col: str = "templates",
    min_frequency: Optional[float] = None,
    aggregate: bool = False,
) -> pd.DataFrame:
    """
    Filter / normalise / merge a repertoire DataFrame against an optional query.

    Parameters
    ----------
    df          : Raw repertoire DataFrame (already loaded from zip).
    filename    : Original member filename — used to derive sample_name.
    query       : Optional clone list to left-join on group_col.
    (remaining parameters documented in merge_query_to_repertoire)
    """
    source_label = filename

    if depthsample is not None:
        #print(f"  [{os.path.basename(filename)}] normalising to {depthsample:,} reads …")
        df = fix_depth(d=df, target=depthsample, rs=random_state)

    for col in (clone_col, group_col, freq_col):
        if col not in df.columns:
            raise ValueError(f"Expected column {col!r} not found in {source_label!r}")

    df = df[df[clone_col].notna()].reset_index(drop=True)
    df = df[df[group_col].notna()].reset_index(drop=True)
    df["clone_ref"] = df.index

    if min_frequency is not None:
        df = df[df[freq_col] > min_frequency].reset_index(drop=True)

    if aggregate:
        df = df.groupby([clone_col, group_col]).sum().reset_index()

    df = df.sort_values(templates_col, ascending=False).reset_index(drop=True)
    df["rank"] = df.index

    sample_name = (
        os.path.basename(source_label)
        .removesuffix(".tsv")
        .removesuffix(".csv")
    )

    result = query.merge(df, how="left", on=group_col) if query is not None else df
    result["name_1"]            = sample_name
    result["total_templates_1"] = df[templates_col].sum()
    result["total_unique_1"]    = len(df)

    return result


# ─────────────────────────────────────────────────────────────────────────────
#  4SINGLE-FILE PIPELINE WORKER  (this is what parmap calls per filename)
# ─────────────────────────────────────────────────────────────────────────────

def _pipeline_worker(
    filename: str,
    *,
    zip_path: str,
    parse_func: Optional[Callable] = None,
    parse_kwargs: dict,
    query_func: Callable,
    query_kwargs: dict,
) -> pd.DataFrame:
    """
    Full pipeline for one zip member:
        open → (optional) parse → query
    Safe to run in a subprocess: opens the zip file independently.
    """
    df = open_df_from_zip(zip_path, filename)

    if parse_func is not None:
        df = apply_parse(parse_func, df, filename, **parse_kwargs)

    query_df = apply_query(query_func, df=df, filename=filename, **query_kwargs)
    return query_df

# ─────────────────────────────────────────────────────────────────────────────
#  5. PUBLIC ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def parmap_zip_pipeline(
    zip_path: str,
    query_func: Callable = merge_query,
    query_kwargs: Optional[dict] = None,
    parse_func: Optional[Callable] = None,
    parse_kwargs: Optional[dict] = None,
    cpus: int = 4,
    pm_pbar: bool = True,
) -> list[pd.DataFrame]:
    """
    Run a parse + query pipeline over every data file in a zip archive
    using parmap for parallelism.  No intermediate files are written.

    Parameters
    ----------
    zip_path      : Path to the .zip archive.
    query_func    : Callable applied after optional parse step.
                    Signature: query_func(df, filename, **query_kwargs) -> DataFrame
                    Defaults to merge_query.
    query_kwargs  : Keyword arguments forwarded to query_func.
    parse_func    : Optional callable applied to the raw DataFrame first.
                    Signature: parse_func(df, filename, **parse_kwargs) -> DataFrame
                    Pass None to skip.
    parse_kwargs  : Keyword arguments forwarded to parse_func.
    cpus          : Number of parallel worker processes.
    pm_pbar       : Show parmap progress bar.

    Returns
    -------
    list[pd.DataFrame]  — one DataFrame per archive member, in namelist order.
    """
    filenames = get_valid_names(zip_path)
    print(f"Found {len(filenames)} files in {os.path.basename(zip_path)}")

    worker = functools.partial(
        _pipeline_worker,
        zip_path=zip_path,
        parse_func=parse_func,
        parse_kwargs=parse_kwargs or {},
        query_func=query_func,
        query_kwargs=query_kwargs or {},
    )

    return parmap.map(worker, filenames, pm_pbar=pm_pbar, pm_processes=cpus)


def add_AH(fold, A, hf_long):
    # NEXT STEP IS JOINING fold to A
    foldA = fold.merge(A[['index','locus','binary','vfamcdr3','logit_odds', 'yule_q_odds','odds_exact']], how = "left", on = ['binary','vfamcdr3','index'])
    assert fold.shape[0] == foldA.shape[0]
    # NEXT STEP IS APPLYING PARTICIPANT HLA 
    hf_long = hf_long
    foldAH = foldA.merge(hf_long, how = "left", 
                        left_on =['name_1','binary'], 
                        right_on = ['sample_id','binary'],
                        )
    assert fold.shape[0] == foldAH.shape[0]
    # ADD A DETECTION
    foldAH['detect'] = (foldAH['productive_frequency'] > 1E-6).astype('int64')
    foldAH['hla_detect'] = ((foldAH['detect'] ==1) & (foldAH['hla_value'] ==1) ).astype('int64')
    foldAH['x1w1'] = foldAH['detect'] * 1
    foldAH['yule_q_odds'] = foldAH['yule_q_odds'].apply(lambda x: max(0, x))
    foldAH['logit_odds'] = foldAH['logit_odds'].apply(lambda x: max(0, x))
    foldAH['x1w2'] = foldAH['detect'] * foldAH['yule_q_odds']
    foldAH['x1w3'] = foldAH['detect'] * foldAH['logit_odds']
    return foldAH




"""
# THERE IS AN EXPLANATION OF WHAT IS HAPPENENING
# parmap iterates filenames
#         ▼
# _pipeline_worker(filename)
#     ├── open_df_from_zip()      # reads zip member into DataFrame, no disk I/O
#     ├── apply_parse()           # optional, swap in any parse_func
#     └── apply_query()           # swap in any query_func (default: merge_query)
### USAGES EXAMPLE
DEFAULT_RANDOM_STATE = 103083
CPUS = 24
from functools import partial
from tcrtest.ui import VfamCDR3, parse_v2x_gt1
# LOAD A PARSE FUNCTION, MUST ACCEPT df, f arguments
parse_v2x_gt1_vfamcdr3_dna =  partial(parse_adapative_v2_gt1, 
     gt1= False, 
     out_cols = ['vfamcdr3','templates','productive_frequency','V','nucleotide','v_resolved','j_resolved'])
     
data_path = Path("/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools") / "data" / "scov_s_hla_tcrs_v1.csv"
q =  pd.read_csv(data_path)
zip_path = '/fh/fast/gilbert_p/kmayerbl/large_cohorts/2026_L1/hsTCRb_v4b/hstcrb_v4b_ctrl_sampleExport.2026-05-12_02-09-16.zip'

results = parmap_zip_pipeline(
    zip_path     = zip_path,
    parse_func   = parse_v2x_gt1_vfamcdr3_dna,
    parse_kwargs = {},
    query_func   = merge_query,
    query_kwargs = dict(query=q,random_state = DEFAULT_RANDOM_STATE, depthsample = 100_000 ),
    cpus         = CPUS,
)
combined = pd.concat(results, ignore_index=True)
"""









