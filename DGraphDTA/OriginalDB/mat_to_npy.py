import os 

os.getcwd()
# os.listdir()

# os.chdir('/mnt/ailon/pdb_test/non-redundant-PDB')
# os.getcwd()
os.chdir('/mnt/ailon/pdb_test')

mat_dir = mat_dir = os.path.join('non-redundant-PDB')


# convert .mat files to numpy files
import os
import numpy as np
import scipy.io
import random

mat_dir = os.path.join('non-redundant-PDB')
output_dir = os.path.join('npy_PDB')
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
file_list = os.listdir(mat_dir)
random.shuffle(file_list)
#inputs = []
#outputs = []
for file in file_list:
    input_file = os.path.join(mat_dir, file)
    input_file = scipy.io.loadmat(input_file)
    mat_value = input_file['mat']
    output_file = os.path.join(output_dir, file.split('.mat')[0]+'.npy')
    np.save(output_file, mat_value
