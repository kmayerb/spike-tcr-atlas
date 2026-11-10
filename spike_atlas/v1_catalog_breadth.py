"""
V1 Manuscript Run : May 27, 2026 (20260527)
catalog_breadth.py

"""
import pandas as pd
import sys
import numpy as np
import parmap
from zip_query import parmap_zip_pipeline, merge_query, parse_v2x_gt1_vfamcdr3_dna, parse_v1x_gt1_vfamcdr3_dna, add_AH
# Note: zip_query.py has tools for parsing files and querying a catalog directly against main files at once 
# using zipfile with parmap, as run path to the assembled data files:
# folder_path = '/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/SPIKE_ATLAS_DATA'
# sys.path.append(folder_path)

DEFAULT_RANDOM_STATE = 103083
CPUS = 24

def get_catalog_detects(
    name,
    zip_path,
    catalog_path,
    parse_function,
    merge_function,
    random_state=DEFAULT_RANDOM_STATE,
    cpus=CPUS,
    depth_sample=100_000):
    """
    Match a TCR catalog to each repertoire in a zip file and return combined hits.

    Parameters
    ----------
    name : str
        Label assigned to the run, stored in the 'query_run_name' column.
    zip_path : str
        Path to the zip file containing raw repertoire files.
    catalog_path : str
        Path to CSV file containing the TCR catalog to query against.
    parse_function : callable
        Parses a raw repertoire file and assigns 'vfamcdr3' from TRBV gene and CDR3.
    merge_function : callable
        Downsamples a repertoire and queries it against the catalog via left join on 'vfamcdr3'.
    random_state : int, optional
        Seed for reproducible downsampling. Default is DEFAULT_RANDOM_STATE.
    cpus : int, optional
        Number of parallel workers. Default is CPUS.
    depth_sample : int, optional
        Number of reads to downsample each repertoire to. Default is 100_000.

    Returns
    -------
    pd.DataFrame
        Combined catalog-hit rows across all repertoires, with 'query_run_name' column
        set to `name`. Column 'name_1' identifies the source repertoire file.
    """
    print(f"Running <get_catalog_detects> with {cpus} cpus with readsampling to {depth_sample} reads.")
    print(f"Loading query from <catalog_path>: {catalog_path}")
    q = pd.read_csv(catalog_path)
    merge_kwargs = dict(query=q, random_state=random_state, depthsample=depth_sample)

    results = parmap_zip_pipeline(
        zip_path     = zip_path,
        parse_func   = parse_function,
        parse_kwargs = {},
        query_func   = merge_function,
        query_kwargs = merge_kwargs,
        cpus         = cpus)
    combined = pd.concat(results)
    combined['query_run_name'] = name
    return combined


def join_detect_to_hla(catalog_hits_df, association_file, hla_file):
    """
    Join catalog hits with HLA association statistics and per-sample HLA predictions.

    Parameters
    ----------
    catalog_hits_df : pd.DataFrame
        Output of `get_catalog_detects`. Must contain columns 'index', 'name_1', and 'binary'.
    association_file : str
        Path to CSV containing HLA association statistics. Must contain column 'i'
        (used as join key) and 'binary' allele labels.
    hla_file : str
        Path to long-format CSV of per-sample HLA predictions. Must contain columns
        'sample_id', 'binary', and 'pred'.

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with HLA association and prediction columns appended,
        as produced by `add_AH`.

    Raises
    ------
    ValueError
        If any 'index' values in `catalog_hits_df` are absent from the association file,
        or if any alleles in `catalog_hits_df` are absent from the HLA file.
    """
    print(f"Loading association file {association_file}")
    A = pd.read_csv(association_file)
    A['index'] = A['i']

    print("Checking that all index in catalog_hits_df found in association file.")
    missing = set(catalog_hits_df['index']) - set(A['index'])
    if missing:
        raise ValueError(f"Indices not found in association file: {missing}")

    print(f"Loading HLA file {hla_file}")
    H = pd.read_csv(hla_file)
    H['hla_value'] = H['pred'].astype('float')

    assert 'name_1' in catalog_hits_df.columns, "'name_1' column missing from catalog_hits_df."
    assert 'sample_id' in H.columns, "'sample_id' column missing from HLA file."

    print("Checking that all alleles are present in HLA file.")
    all_alleles = catalog_hits_df['binary'].unique()
    all_H_allele = H.binary.unique()
    for allele in all_alleles:
        assert allele in all_H_allele, f"{allele} not found in HLA file."

    print("Joining HLA file to catalog matched hits.")
    return add_AH(fold=catalog_hits_df, A=A, hf_long=H)


def summarize_hla_breadth(df):
    """
    Summarize HLA-matched TCR breadth per sample.

    Breadth is defined as the sum of weighted clone counts ('x1w2') for HLA-matched
    detections, divided by total unique clones in the repertoire.

    Parameters
    ----------
    df : pd.DataFrame
        Output of `join_detect_to_hla`. Must contain columns 'hla_detect',
        'name_1', 'total_unique_1', and 'x1w2'.

    Returns
    -------
    pd.DataFrame
        One row per sample with columns: 'name_1', 'total_unique_1', 'hla_detect',
        'x1w2', 'breadth', and 'log2_breadth'.
    """
    df2 = df.query('hla_detect == 1').copy()
    df2_sum = df2.groupby(['name_1', 'total_unique_1', 'hla_detect'])['x1w2'].sum().reset_index()
    df2_sum['breadth'] = df2_sum['x1w2'] / df2_sum['total_unique_1']
    df2_sum['log2_breadth'] = np.log2(df2_sum['breadth'])
    return df2_sum


def summarize_xhla_breadth(df):
    """
    Summarize TCR breadth for detections that lack HLA restriction (cross-HLA hits).

    Selects clones that are catalog-detected ('detect == 1') but not HLA-matched
    ('hla_detect == 0'), then computes breadth as for `summarize_hla_breadth`.

    Parameters
    ----------
    df : pd.DataFrame
        Output of `join_detect_to_hla`. Must contain columns 'hla_detect', 'detect',
        'name_1', 'total_unique_1', and 'x1w2'.

    Returns
    -------
    pd.DataFrame
        One row per sample with columns: 'name_1', 'total_unique_1', 'hla_detect',
        'x1w2', 'breadth', and 'log2_breadth'.
    """
    df2 = df.query('hla_detect == 0 & detect == 1').copy()
    df2_sum = df2.groupby(['name_1', 'total_unique_1', 'hla_detect'])['x1w2'].sum().reset_index()
    df2_sum['breadth'] = df2_sum['x1w2'] / df2_sum['total_unique_1']
    df2_sum['log2_breadth'] = np.log2(df2_sum['breadth'])
    return df2_sum


def hla_breadth_pipeline(
    name,
    zip_path,
    hla_path,
    catalog_path     = 'data/scov_s_hla_tcrs_v1.csv', # as run '/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/scov_s_hla_tcrs_v1.csv', # catalog of features 
    association_path = 'data/scov_s_hla_tcrs_v1_associations.csv',# as run '/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/20250904_7135_HLA_S_feature_set_odds_yule_logistic_post_vac_all_studies_vs_unexposed.csv', # feature frequency in unexposed donors
    parse_function   = parse_v2x_gt1_vfamcdr3_dna,
    merge_function   = merge_query,
    random_state     = DEFAULT_RANDOM_STATE,
    cpus             = CPUS,
    depth_sample     = 100_000):
    """
    End-to-end pipeline: catalog detection, HLA joining, and breadth summarization.

    Runs catalog matching across all repertoires in a zip archive, joins HLA
    association and prediction data, and returns HLA-restricted and cross-HLA
    breadth summaries. Samples with zero HLA-matched detections are retained
    with breadth set to 0 and flagged via 'zero_breadth'.

    Parameters
    ----------
    name : str
        Label for this pipeline run, passed through to `get_catalog_detects`.
    zip_path : str
        Path to zip archive of repertoire files.
    hla_path : str
        Path to long-format CSV of per-sample HLA predictions.
    catalog_path : str, optional
        Path to the SARS-CoV-2 spike HLA-TCR catalog CSV.
    association_path : str, optional
        Path to CSV of HLA association statistics (odds ratios, etc.).
    parse_function : callable, optional
        Repertoire parse function; default parses v2 export with DNA sequences.
    merge_function : callable, optional
        Catalog query function; default performs a downsampled left-join merge.
    cpus : int, optional
        Number of parallel workers. Default is CPUS.
    depth_sample : int, optional
        Reads to downsample per repertoire. Default is 100_000.

    Returns
    -------
    dict with keys:
        catalog_hits_hla : pd.DataFrame
            Full joined DataFrame of all catalog hits with HLA columns.
        catalog_hits_hla_detects : pd.DataFrame
            Subset of catalog_hits_hla where 'detect == 1'.
        catalog_hits_hla_summary : pd.DataFrame
            Per-sample HLA-restricted breadth summary, including zero-breadth
            samples flagged with 'zero_breadth == True'.
        catalog_hits_xhla_summary : pd.DataFrame
            Per-sample cross-HLA breadth summary.
        all_samples : pd.Series
            Sorted list of all sample names found in the repertoire zip.

    Notes
    -----
    Zero-breadth samples (no HLA-matched detections) are appended to
    'catalog_hits_hla_summary' with breadth=0 and zero_breadth=True rather than
    dropped, so downstream analyses retain the full sample set. Use caution: some
    zero-breadth samples may reflect very low read depth rather than true absence
    of response.
    """
    catalog_hits = get_catalog_detects(
        name           = name,
        zip_path       = zip_path,
        catalog_path   = catalog_path,
        parse_function = parse_function,
        merge_function = merge_function,
        random_state   = random_state,
        cpus           = cpus,
        depth_sample   = depth_sample)

    all_samples = pd.Series(sorted(list(catalog_hits.name_1.unique())))
    all_hla_samples = pd.Series(sorted(list(pd.read_csv(hla_path)['sample_id'].unique())))

    missing_hla = all_samples[~all_samples.isin(all_hla_samples)]

    print(f"Number of samples in repertoire zip:    {len(all_samples)}")
    print(f"Number of samples with HLA data:        {all_samples.isin(all_hla_samples).sum()}")
    print(f"Number of samples missing HLA data:     {len(missing_hla)}")
    if not missing_hla.empty:
        print("Samples missing HLA data:")
        for s in missing_hla:
            print(f"  {s}")

    catalog_hits_hla = join_detect_to_hla(
        catalog_hits_df  = catalog_hits,
        association_file = association_path,
        hla_file         = hla_path)

    catalog_hits_hla_detects  = catalog_hits_hla.query('detect == 1')
    catalog_hits_hla_summary  = summarize_hla_breadth(catalog_hits_hla)
    catalog_hits_xhla_summary = summarize_xhla_breadth(catalog_hits_hla)

    results = {
        'catalog_hits_hla'         : catalog_hits_hla,
        'catalog_hits_hla_detects' : catalog_hits_hla_detects,
        'catalog_hits_hla_summary' : catalog_hits_hla_summary,
        'catalog_hits_xhla_summary': catalog_hits_xhla_summary,
        'all_samples'              : all_samples,
    }

    # Retain zero-breadth samples so they are not silently dropped downstream.
    ix = results['all_samples'].isin(results['catalog_hits_hla_summary']['name_1'])
    if ix.mean() != 1.0:
        zero_breadth_samples = results['all_samples'][~ix]
        zb = (results['catalog_hits_hla']
              [results['catalog_hits_hla']['name_1'].isin(zero_breadth_samples)]
              [['name_1', 'total_unique_1']]
              .groupby('name_1').head(1)
              .reset_index(drop=True))
        zb['hla_detect']   = 1
        zb['x1w2']         = 0
        zb['breadth']      = 0
        zb['zero_breadth'] = True
        # itialize zero breadth as false
        results['catalog_hits_hla_summary']['zero_breadth'] = False
       
        
        results['catalog_hits_hla_summary'] = pd.concat(
            [results['catalog_hits_hla_summary'], zb]).reset_index(drop=True)
        
         # indicate samples that lacked any hla information and will have zero breadth for that reason
        results['catalog_hits_hla_summary']['missing_hla'] = results['catalog_hits_hla_summary']['name_1'].isin(missing_hla)

    return results

"""
#### Repertoire Files
Repertoire files are massive and each is downloadable from the primary source.
The script processes them as individual `.zip` folders by source project.

Datasets generated for this study -- 10.5281/zenodo.23965835
* 2.1G cov_boost_day0_day14.zip (Oxford - This Study)
* 833M cov2.zip (Oxford - Prior Trial Data)

Previously published datasets:
* 9.1G cov3006.zip (https://doi.org/10.6084/m9.figshare.31807288)
* 489M ad26_alter.zip (Ad26.COV2.S; paper DOI 10.1038/s41586-021-03681-2)
* 2.9G delmonte.zip (https://doi.org/10.21417/OMD2023JACI)
* 3.4G elyanow.zip (clients.adaptivebiotech.com/pub/elyanow-2022-jci/10.21417/RE2022JCI)
* 3.9G emerson_natgen_v2.zip (clients.adaptivebiotech.com/pub/emerson-2017-natgen/ 10.21417/B7001Z)
* 705M v4b_ctrl.zip (https://doi.org/10.21417/ADPT2020V4CD)


#### Inferred HLA genotype prediction files

We used [TCR2HLA](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013767) to make high-confidence predictions for common HLA alleles.
These predictions are in long format and are joined by the pipeline to the
query outputs to designate HLA-matched hits.

* cov_boost_hla_predictions.csv
* cov2_hla_predictions.csv
* cov3006_hla_predictions.csv
* ad26_alter_hla_predictions.csv
* delmonte_hla_predictions.csv
* elyanow_hla_predictions.csv
* emerson_natgen_v2_hla_predictions.csv
* v4b_ctrl_hla_predictions.csv
"""


if __name__ == "__main__":
    import os
    #archive_name = 'cov_boost'
    DEFAULT_RANDOM_STATE = 103083
    CPUS = 24
    SPIKE_ATLAS_DATA = '/fh/fast/gilbert_p/kmayerbl/spike_atlas_data'
    project = {
        'ad26_alter': {
            'name'     : 'ad26_alter',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/ad26_alter.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/ad26_alter_hla_predictions.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
        },
        'elyanow': {
            'name'     : 'elyanow',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/elyanow.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/elyanow_hla_predictions.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
        },
        'cov2': {
            'name'     : 'cov2',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/cov2.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/cov2_hla_predictions.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
        },
        'cov3006': {
            'name'     : 'cov3006',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/cov3006.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/cov3006_hla_predictions.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
        },
        'delmonte': {
            'name'     : 'delmonte',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/delmonte.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/delmonte_hla_predictions.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
        },
        'v4b_ctrl': {
            'name'     : 'v4b_ctrl',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/v4b_ctrl.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/v4b_ctrl_hla_predictions.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
        },
        'emerson_natgen_v2': {
            'name' :'emerson_natgen_v2',
            'zip_path'   : f'{SPIKE_ATLAS_DATA}/emerson_natgen_v2.zip',
            'hla_path'   : f'{SPIKE_ATLAS_DATA}/emerson_natgen_v2_hla_predictions.csv',
            'outdir'     : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'}, 
        'emerson_natgen_v2_genotyped': {
            'name' :'emerson_natgen_v2_genotyped',
            'zip_path'   : f'{SPIKE_ATLAS_DATA}/emerson_natgen_v2.zip',
            'hla_path'   : f'{SPIKE_ATLAS_DATA}/emerson_natgen_v2_hla_genotypes.csv',
            'outdir'     : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'},

        'cov2_genotyped': {
            'name'     : 'cov2_genotyped',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/cov2.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/cov2_hla_genotypes.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
            },
        'cov_boost_day0_day14_genotyped': {
            'name'     : 'cov_boost_day0_day14_genotyped',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/cov_boost_day0_day14.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/cov_boost_hla_genotypes.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
        },
        'cov_boost_day0_day14': {
            'name'     : 'cov_boost_day0_day14',
            'zip_path' : f'{SPIKE_ATLAS_DATA}/cov_boost_day0_day14.zip',
            'hla_path' : f'{SPIKE_ATLAS_DATA}/cov_boost_hla_predictions.csv',
            'outdir'   : f'{SPIKE_ATLAS_DATA}/estimate_20260527/'
        }
        }

    import os
    print("Checking project file paths...\n")
    all_present = True
    for proj_name, cfg in project.items():
        zip_ok = os.path.exists(cfg['zip_path'])
        hla_ok = os.path.exists(cfg['hla_path'])
        if not zip_ok or not hla_ok:
            all_present = False
            print(f"[MISSING] {proj_name}")
            if not zip_ok:
                print(f"          zip_path : {cfg['zip_path']}")
            if not hla_ok:
                print(f"          hla_path : {cfg['hla_path']}")
    if all_present:
        print("All zip and HLA paths present for all projects.")
    else:
        print("\nResolve missing files before running the pipeline.")

    # Peform pipline with imputed HLA for all emerson data
    # Datasets generated for this study -- 10.5281/zenodo.23965835 ('cov_boost_day0_day14_genotyped','cov2')
    # Other datasets are from public data 
    for run_name in ['cov_boost_day0_day14_genotyped','cov2']: # 'cov_boost_day0_day14', 'ad26_alter', 'elyanow', 'cov_boost', 'cov2', 'cov3006', 'delmonte', 'v4b_ctrl', 'emerson_natgen_v2', 'emerson_natgen_v2_genotyped']:
        print(f"Computing SARS-S-CoV-2 HLA Spike BREADTH for {run_name}")
        run = project[run_name]
        x = hla_breadth_pipeline(
            name     = run['name'],
            zip_path = run['zip_path'],
            hla_path = run['hla_path'],
            catalog_path     = 'data/scov_s_hla_tcrs_v1.csv', # as actually run  on HPC '/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/scov_s_hla_tcrs_v1.csv', # catalog of features 
            association_path = 'data/scov_s_hla_tcrs_v1_associations.csv',# as actually run on HPC '/fh/fast/gilbert_p/kmayerbl/prepost_mrna_vac_code_only/projects/00_ID159_tools/data/20250904_7135_HLA_S_feature_set_odds_yule_logistic_post_vac_all_studies_vs_unexposed.csv', # feature frequency in post vaccine vs unexposed donors
            parse_function   = parse_v2x_gt1_vfamcdr3_dna,
            merge_function   = merge_query,
            random_state     = DEFAULT_RANDOM_STATE,
            cpus             = CPUS,
            depth_sample     = 100_000)
        outfile = f"{run['name']}_spike_breadth_estimate_x1w2.csv"
        outfile_path = os.path.join(run['outdir'], outfile )
        x['catalog_hits_hla_summary'].to_csv(outfile_path , index = False)
        outfile2 = f"{run['name']}_spike_breadth_all_hits.csv"
        outfile_path2 = os.path.join(run['outdir'], outfile2 )
        x['catalog_hits_hla'].query('detect == 1').to_csv(outfile_path2 , index = False)
        print(outfile_path)
        #!head {outfile_path}

