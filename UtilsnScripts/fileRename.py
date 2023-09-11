import os
import re
import numpy as np

path = './rename/' # <--- You add path here
files = os.listdir(path)
ordered_files = sorted(files, key=lambda x: (int(re.sub('\D','',x)),x))
iter = 0


for file in ordered_files:
	name = 'snapshot-'
	if iter < 10:
	    name += '000' + str(iter)
	elif iter >= 10 and iter < 100:
	    name += '00' + str(iter)
	elif iter >= 100 and iter < 1000:
	    name += '0' + str(iter)

	os.rename(os.path.join(path, file), os.path.join(path, name + '.tiff'))
	iter = iter+1
