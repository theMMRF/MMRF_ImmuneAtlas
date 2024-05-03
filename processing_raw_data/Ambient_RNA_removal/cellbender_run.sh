
docker run -v /path/to/input/folder/${sample_id}/outs:/inputs -v /path/to/input/folder/cellbender:/outputs -t us.gcr.io/broad-dsde-methods/cellbender:latest cellbender remove-background --input /inputs/raw_feature_bc_matrix.h5 --output /outputs/raw_feature_bc_matrix_cb_${f}.h5 --expected-cells 5000 --total-droplets-included 20000 --fpr 0.01 --epochs 150"


