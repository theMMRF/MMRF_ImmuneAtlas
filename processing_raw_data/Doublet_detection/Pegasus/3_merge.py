# this is a python script to merge all the output txt file to a single one.In adition, add the sample_id after the barcode, like ATTTTTCGGTA-1_5479_56789_SM.
import os


w = open('merged_results.txt', 'a')

# write the head of the file
f = open('data/' + os.listdir('data/')[0] + '/output/' + os.listdir('data/')[0] + '_doublet_detection_result.txt', "r")
head = f.readlines()[0]
w.write(head)
f.close()

# write all the other records expect head of each file
for files in os.listdir('data/'):
	path = 'data/' + files + '/output/' + files + '_doublet_detection_result.txt'
	if not os.path.exists(path):
		continue
	print('path=', path)
	print('sampleid=', files)

	f = open(path, "r")
	lines = f.readlines()[1:]
	for line in lines:
		line = line.split('\t')
		line[0] = line[0]+'-'+files
		line = line[:4]+line[5:]
		line = "\t".join(line)
		w.write(line)

	f.close()

w.close()

