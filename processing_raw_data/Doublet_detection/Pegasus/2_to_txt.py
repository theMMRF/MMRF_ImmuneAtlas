# -*- coding: utf-8 -*-
import scanpy as sc
import os

def main():
  with open('/diskmnt/Projects/MMRF_analysis/QC/Pegasus/input_files.txt') as file:
      lines = file.readlines()
  i = 0
  fail = []
  for line in lines:
    input = '/diskmnt/Projects/MMRF_analysis/QC/Pegasus/data/' + line.split('/')[5] + '/output/' + line.split('/')[5] + '_raw_feature_bc_matrix.h5ad'
    output = '/diskmnt/Projects/MMRF_analysis/QC/Pegasus/data/' + line.split('/')[5] + '/output/' + line.split('/')[5] + '_doublet_detection_result.txt'
    i += 1
    print(i)
    print(input)
    print(output)
    if os.path.exists(input):
        to_txt(input,output)
    else:
        fail.append(input)
  print(fail)

def to_txt(input,output):
  data = sc.read_h5ad(input)
  data = data.obs
  data.to_csv(output, header=True, index=True, sep='\t', mode='a')



if __name__ == '__main__':
    main()
