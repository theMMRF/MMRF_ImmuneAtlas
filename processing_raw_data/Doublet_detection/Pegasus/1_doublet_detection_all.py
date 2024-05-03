# -*- coding: utf-8 -*-
import pegasus as pg
import os

def main():
  with open('/diskmnt/Projects/MMRF_analysis/QC/Pegasus/input_files.txt') as file:
      lines = file.readlines()
  i = 0
  fails = []
  for line in lines:
    input = line[:-1]
    output = '/diskmnt/Projects/MMRF_analysis/QC/Pegasus/data/' + line.split('/')[5] + '/output/' + line.split('/')[5] + '_raw_feature_bc_matrix.h5ad'
    i += 1
    print(i)
    print(input)
    print(output)
    try :
        if not os.path.exists(output):
            doublet_detection(input,output)
    except Exception:
        fails.append(input)
  with open('fail_files.txt',"w") as f:
    for file in fails:
        file = file + '\n'
        f.write(file)

def doublet_detection(input,output):
  data = pg.read_input(input)
  #data
  #data.X

  #data.obs.head()

  #data.obs['Channel'].value_counts()

  #data.var.head()

  #data.uns['genome']

  #data.uns['modality']

  pg.qc_metrics(data, min_genes=500, max_genes=6000, mito_prefix='MT-', percent_mito=10)


  df_qc = pg.get_filter_stats(data)
  #df_qc


  pg.qcviolin(data, plot_type='gene', dpi=100)



  pg.qcviolin(data, plot_type='count', dpi=100)


  pg.qcviolin(data, plot_type='mito', dpi=100)


  pg.filter_data(data)

  pg.identify_robust_genes(data)

  #data.obs['Channel'].value_counts()


  pg.log_norm(data)

  data_trial = data.copy()

  pg.highly_variable_features(data_trial)

  data_trial.var.loc[data_trial.var['highly_variable_features']].sort_values(by='hvf_rank')

  pg.hvfplot(data_trial, dpi=200)

  pg.pca(data_trial)

  coord_pc1 = data_trial.uns['PCs'][:, 0]
  #coord_pc1

  data_trial.var.loc[data_trial.var['highly_variable_features']].index.values

  #data_trial.obsm['X_pca'].shape

  pg.neighbors(data_trial)

  print(f"Get {data_trial.obsm['pca_knn_indices'].shape[1]} nearest neighbors (excluding itself) for each cell.")
  #data_trial.obsm['pca_knn_indices']

  #data_trial.obsm['pca_knn_distances']

  pg.louvain(data_trial)

  #data_trial.obs['louvain_labels'].value_counts()

  pg.compo_plot(data_trial, 'louvain_labels', 'Channel')

  pg.tsne(data_trial)

  pg.scatter(data_trial, attrs=['louvain_labels', 'Channel'], basis='tsne')

  pg.highly_variable_features(data, batch='Channel') 
  pg.pca(data)
  pca_key = pg.run_harmony(data)

  #data.obsm['X_pca_harmony'].shape

  pg.neighbors(data, rep=pca_key)
  pg.louvain(data, rep=pca_key)

  pg.compo_plot(data, 'louvain_labels', 'Channel')

  pg.tsne(data, rep=pca_key)

  pg.scatter(data, attrs=['louvain_labels', 'Channel'], basis='tsne')

  pg.umap(data, rep=pca_key)

  pg.scatter(data, attrs=['louvain_labels', 'Channel'], basis='umap')

  pg.de_analysis(data, cluster='louvain_labels')

  marker_dict = pg.markers(data)

  marker_dict['1']['up'].sort_values(by='log2FC', ascending=False)

  pg.volcano(data, cluster_id = '1', dpi=200)

  #pg.write_results_to_excel(marker_dict, "MantonBM_subset.de.xlsx")

  celltype_dict = pg.infer_cell_types(data, markers = 'human_immune')
  cluster_names = pg.infer_cluster_names(celltype_dict)

  celltype_dict['1']


  pg.annotate(data, name='anno', based_on='louvain_labels', anno_dict=cluster_names)
  #data.obs['anno'].value_counts()

  pg.scatter(data, attrs='anno', basis='tsne', dpi=100)

  pg.scatter(data, attrs='anno', basis='umap', legend_loc='on data', dpi=150)

  #data

  data.select_matrix('raw.X')
  #data

  data.select_matrix('X')

  pg.diffmap(data, rep=pca_key)

  # data.obsm['X_diffmap'].shape

  pg.fle(data)

  pg.scatter(data, attrs='anno', basis='fle')

  #pg.write_output(data, "MantonBM_result.zarr.zip")

  # -----------------------------------------------------------------------------------------------------------------
  # doublet detection code

  pg.infer_doublets(data, channel_attr = 'Channel', clust_attr = 'anno',plot_hist = None)

  pg.scatter(data,attrs=['anno','doublet_score'], basis='umap', wspace=1.2)

  # data.uns['pred_dbl_cluster']

  pg.mark_doublets(data)
  # data.obs['demux_type'].value_counts()
  pg.scatter(data, attrs = ['anno', 'demux_type'], legend_loc = ['on data', 'right margin'], 
            wspace = 0.1,alpha = [1.0, 0.8], palettes = 'demux_type:gainsboro,red')

  pg.qc_metrics(data, select_singlets=True)
  pg.filter_data(data)

  pg.highly_variable_features(data, batch='Channel')
  pg.pca(data)
  pca_key = pg.run_harmony(data)
  pg.neighbors(data,rep=pca_key)
  pg.louvain(data,rep=pca_key)
  pg.umap(data,rep=pca_key)

  pg.de_analysis(data, cluster='louvain_labels')
  celltype_dict = pg.infer_cell_types(data, markers = 'human_immune',de_test='mwu')
  cluster_names = pg.infer_cluster_names(celltype_dict)
  pg.annotate(data, name='anno', based_on='louvain_labels', anno_dict=cluster_names)

  pg.scatter(data,attrs='anno',legend_loc='on data',basis='umap')

  pg.write_output(data, output)
  

if __name__ == '__main__':
    main()
