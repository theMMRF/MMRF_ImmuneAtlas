while read LINE  
do  
  sample_id=$(echo $LINE | cut -d'/' -f6)

  echo $sample_id
  mkdir -p /diskmnt/Projects/MMRF_analysis/QC/Pegasus/data/$sample_id/output/
    
done <input_files.txt
