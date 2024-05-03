#Lijun Yao 
#lijunyao@wustl.edu
#script was adapted from Reyka Jayasinghe

#https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb

import scrublet as scr
#If error in identifying above rerun `pip install scrublet` should fix all issues
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import gzip
import pandas as pd

import sys, getopt

def main(argv):
   sample = ''
   cellranger = ''
   try:
      opts, args = getopt.getopt(argv,"hs:c:",["sname=","celloutput="])
   except getopt.GetoptError:
      print ('2_Multiplet.py -s <sample> -c <cellranger>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('2_Multiplet.py -s <sample> -c <cellranger>')
         print('python 2_Multiplet.py -s RM024R1-XBn1_1 -c /diskmnt/Datasets/mmy_scratch/MetNet/Brain/cellranger/RM024R1-XBn1_1/RM024R1-XBn1_1/')
         sys.exit()
      elif opt in ("-s", "--sname"):
         sample = arg
      elif opt in ("-c", "--celloutput"):
         cellranger = arg
   print ('Sample is "', sample)
   print ('Cellranger Directory is "', cellranger)

if __name__ == "__main__":
   main(sys.argv[1:])

#Define input arguments
cellranger=sys.argv[4]
sample=sys.argv[2]

#input_dir = '/diskmnt/Datasets/mmy_scratch/MetNet/Brain/cellranger/RM024R1-XBn1_1/RM024R1-XBn1_1/outs/filtered_feature_bc_matrix/'
#Only works on filtered_feature matrix - raw matrix takes too much memory
input_dir=cellranger+"outs/raw_feature_bc_matrix/"
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
barcodes=pd.read_csv(input_dir+'barcodes.tsv.gz',compression='gzip',header=None)

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
scrub.call_doublets(threshold=0.20)


###Barcode and Doublet Prediction File
#https://github.com/swolock/scrublet/issues/5
outputfile=sample+"_scrublet_output_table.csv"
df = pd.DataFrame({
    'doublet_score': scrub.doublet_scores_obs_,
    'predicted_doublet': scrub.predicted_doublets_
})
barcodes.columns = ['Barcodes']
dffinal=pd.concat([barcodes,df],axis=1)
dffinal.to_csv(outputfile, index=False)

###DOUBLET HISTOGRAM and PREDICTED DOUBLET RATE
plt.figure(figsize=(3, 3))
scrub.plot_histogram()
plt.savefig(sample+"_doublets_hist.pdf")
plt.close()

print('Running UMAP...')
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
print('Done.')

plt.figure(figsize=(3, 3))
scrub.plot_embedding('UMAP', order_points=True);
plt.savefig(sample+"_doublets_umap.pdf")
plt.close()
