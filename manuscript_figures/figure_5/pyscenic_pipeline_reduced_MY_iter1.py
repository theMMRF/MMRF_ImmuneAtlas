from tkinter.ttk import Progressbar
import pandas as pd
import numpy as np
import os, glob
import pickle

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

#from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.binarization import binarize

import scanpy as sc
import seaborn as sns

import matplotlib.pyplot as plt
from dask.distributed import Client, LocalCluster

import loompy as lp
os.chdir('/opt/megaseq-data/sarthak/projects/MMRF')

DATA_FOLDER="./scenic_MY_reduced" ## Change here
RESOURCES_FOLDER="./data"
DATABASE_FOLDER = "./data"

DATABASES_GLOB = os.path.join("/opt/pmedseq-data/sarthak/resources/scenic/resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/*feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join("/opt/pmedseq-data/sarthak/resources/scenic/resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl")

HG_TFS_FNAME = os.path.join("outs/allTFs_hg38.txt")
#SC_EXP_FNAME = os.path.join(RESOURCES_FOLDER, "GSE60361_C1-3005-Expression.txt")
ADATA_FNAME = os.path.join("./objects/anndata_4000_42266Cells_MY_REDUCE_counts.h5ad") ## chanfe here
MIN_GENES = 20

def coexpression(ex_matrix, tf_names, adjacencies_fname, module_fname):

    ### Phase 1 ###

    adjancencies = grnboost2(expression_data=ex_matrix, tf_names=tf_names, verbose=False)
    adjancencies.to_csv(adjacencies_fname, index=False, sep='\t')

    modules = list(modules_from_adjacencies(adjancencies, ex_matrix, min_genes = MIN_GENES)) #### Min 10 genes

    with open(module_fname, 'wb') as f:
        pickle.dump(modules, f)

    return modules

def pruning(dbs, modules, motifs_fname, regulons_fname):

    ### Phase 2 ###

    #MOTIF_ANNOTATIONS_FNAME = os.path.join('data', "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, num_workers=6)
    df.to_csv(motifs_fname)

    regulons = df2regulons(df)

    with open(regulons_fname, 'wb') as f:
        pickle.dump(regulons, f)
    return regulons

if __name__ == "__main__":

    print("[*] Loading tf")
    tf_names = load_tf_names(HG_TFS_FNAME)

    db_fnames = glob.glob(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]

    print("[*] Loading Ranking DB..")
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

    for i in range(0,25):


        group = 'subcluster_V03072023'
        adata = sc.read_h5ad(ADATA_FNAME) ## Change
        ex_matrix=adata.to_df()

        DATASET_ID = 'NB_run'+str(i) # change
        adjacencies_fname = os.path.join(DATA_FOLDER, "{}_adjacencies.tsv".format(DATASET_ID))
        module_fname = os.path.join(DATA_FOLDER, "{}_modules.p".format(DATASET_ID))
        motifs_fname = os.path.join(DATA_FOLDER, "{}_motifs.csv".format(DATASET_ID))
        regulons_fname = os.path.join(DATA_FOLDER, "{}_regulons.p".format(DATASET_ID))
        AUCELL_FNAME = os.path.join(DATA_FOLDER, "{}_aucell.csv".format(DATASET_ID))

        BIN_MTX_FNAME = os.path.join(DATA_FOLDER, '{}_bin.csv'.format(DATASET_ID))
        #THR_FNAME = os.path.join(DATA_FOLDER, '{}_thresholds.csv'.format(DATASET_ID))
        DF_RESULTS_FNAME = os.path.join(DATA_FOLDER, '{}_df_results.csv'.format(DATASET_ID))

        print("[*] Adj..")
        modules = coexpression(ex_matrix, tf_names, adjacencies_fname, module_fname)

        print("[*] Pruning..")
        regulons = pruning(dbs, modules, motifs_fname, regulons_fname)

        with open(module_fname, 'rb') as f:
            modules = pickle.load(f)

        ### Get pruned modules for targets with cis regulatory footprints 
        with open(regulons_fname, 'rb') as f:
            regulons = pickle.load(f)


    #     auc_mtx = aucell(ex_matrix, regulons,  num_workers=4)
    #     auc_mtx.to_csv(AUCELL_FNAME)

    #     add_scenic_metadata(adata, auc_mtx, regulons)

    #     #bin_mtx, thresholds = binarize(auc_mtx)
    #    # bin_mtx.to_csv(BIN_MTX_FNAME)
    #     #thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(THR_FNAME)

    #     ##regulons = [r.rename(r.name.replace('(',' (')) for r in regulons]


    #     df_obs = adata.obs
    #     signature_column_names = list(df_obs.select_dtypes('number').columns)
    #     signature_column_names = list(filter(lambda s: s.startswith('Regulon('), signature_column_names))
    #     df_scores = df_obs[signature_column_names + [group]]
    #     df_results = ((df_scores.groupby(by=group).mean() - df_obs[signature_column_names].mean())/ df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})
    #     df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
    #     df_results[(df_results.Z >= 0)].sort_values('Z', ascending=False)
        
    #     df_results.to_csv(DF_RESULTS_FNAME)