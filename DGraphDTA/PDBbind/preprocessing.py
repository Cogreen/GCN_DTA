#
# Data: 20201106  
# Writer: Gimlet

# creating data when data is changed as an incorrect dataset 


## Load sample data
# load data
import pandas as pd
import os 
os.getcwd()
os.chdir('/mnt/ailon/pdb_test/pdb_sample_work/make_kd_candidate')
pdb_sample = pd.read_csv('pdb_kd_004_sample.csv')
pdb_sample

# search mol2 file 
mol2_path = '/mnt/ailon/pdb_test/mol2/'
mol_list = []
for pdb_code in pdb_sample['PDB_code']:
    # mole(ligand)
    mol_list.append(os.path.join(mol2_path + pdb_code + '_ligand.mol2'))

mol_list[0]

pdb_sample



## Test for making Ligand Smiles 
import rdkit
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

a = rdkit.Chem.rdmolfiles.MolFromMol2File(mol_list[1])
mol_list[0]
a = Chem.MolFromMol2File(mol_list[1])
a
print(a)

b = Chem.MolToSmiles(a)
b = Chem.MolToSmiles(a,isomericSmiles=False)
b


# mol file from mol2 file: via rdkit and openbabel # Use openbable if rdkit doesn't work
real_list = []
drugs = []

for i in range(len(mol_list)):
    try:
        mol2 = Chem.MolFromMol2File(mol_list[i])
        if mol2 != None:
            real_list.append(mol_list[i])
            drugs.append(Chem.MolFromMol2File(mol_list[i]))
        else:
            pass
            #os.system('obabel -imol2 ' + mol_list[i] + ' -omol -O ' + mol_list[i].split('.')[0] + '.mol')
            # obabel -i바꿀파일확장자 바꿀파일이름 -o저장할파일확장자 -O 저장할파일이름.mol
            #mol = Chem.MolFromMolFile(mol_list[i].split('.')[0] + '.mol')
            #if mol != None:
            #    real_list.append(mol_list[i].split('.')[0] + '.mol')
            #    drugs.append(Chem.MolFromMolFile(mol_list[i].split('.')[0] + '.mol'))
            #else:
            #    pass
    except:
        pass
        
len(drugs)
       
    
    
    
## Make ligand Dictionary: {pdb_code: smiles}
real_list_pdb_code = [] 

for i in real_list:
    real_list_pdb_code.append(i.split('/')[5][0:4])
len(real_list_pdb_code)

drug_smiles = []
for i in drugs:
    drug_smiles.append(Chem.MolToSmiles(i,isomericSmiles=False))
len(drug_smiles)


# Make dictionary form
ligand_dic = {}
for real_list_pdb_code, drug_smiles in zip(real_list_pdb_code, drug_smiles):
    ligand_dic[real_list_pdb_code] = drug_smiles


# save ligand_dic
import json

os.chdir('/mnt/ailon/pdb_test/')
with open('ligand_pdb_dic.txt', 'w') as f:
    f.write(json.dumps(ligand_dic))


# load dataset
from collections import OrderedDict
import json, pickle

current_path = os.getcwd() + '/'
ligands = json.load(open(current_path + 'ligand_pdb_dic.txt'), object_pairs_hook=OrderedDict)

ligands



#####################################################################
# Make affinity Matrix : 77 x 77 (available dataset)                #
#####################################################################

# Make matrix 77 x 77 with value 0
import numpy as np

matrix_0 = [[0 for x in range(77)] for x in range(77)]
matrix_array = np.array(matrix_0, dtype = np.float64)
matrix_array.shape


# Available affinity from available lignad
import json, pickle
ligand_dir = '/mnt/ailon/pdb_test/pdb_sample_work/make_kd_candidate/'
ligands = json.load(open(ligand_dir + 'ligand_sample.json'))
#     protein_dir = '/mnt/ailon/pdb_test/pdb_sample_work/make_kd_candidate/'
#     proteins = json.load(open(protein_dir + 'protein_sample.json'))  


# Load the protein 
# Load sample data
import pandas as pd
import os 

os.getcwd()
os.chdir('/mnt/ailon/pdb_test/pdb_sample_work/make_kd_candidate')

pdb_sample = pd.read_csv('pdb_kd_004_sample.csv')
pdb_sample


# Make available affinity
affinityAvailable = []
for i in ligands.keys():
    for j in range(len(pdb_sample['PDB_code'])):
        if i == pdb_sample['PDB_code'][j]:
            affinityAvailable.append(pdb_sample['affinity_nM'][j])
            
            
# change the matrix value based on affinity value(unit: nM)
for i in range(len(affinityAvailable)):
    matrix_array[i][i] = affinityAvailable[i]
matrix_array


len(matrix_array)


# Save affinityAvailable data
import pickle
with open('Y_sampleAvailable1106.txt', 'wb') as f:
    pickle.dump(matrix_array, f)

# Check the saved path
os.getcwd()


import pandas as pd

pdb = pd.read_csv('pdb_kd_samp_matching.csv')
pdb

pdb_sample = pd.read_csv('/mnt/ailon/pdb_test/pdb_sample_work/pdb_sample_1104_test.csv')
len(pdb_sample)



########################## remake smaple ###################33

tmp = []
for i in range(len(pdb_sample)):
    for j in range(len(pdb)):
        if pdb_sample['target_key'][i] == pdb['PDB_code'][j]:
            pdb_sample['affinity'][i] = pdb['affinity_nM'][j]
#            tmp.append(pdb['affinity_nM'][j])
#            print(pdb_sample['affinity'][i])
        else:
            pass
            
            
import numpy as np
pdb_sample['affinity'] = [-np.log10(y / 1e9) for y in pdb_sample['affinity']]
            
            
import os
os.getcwd()

pdb_sample = pd.DataFrame(pdb_sample)


pdb_sample.to_csv('pdb_sample_1106.csv', header = True, index = False)


pdb_sample = pd.read_csv('pdb_sample_1106.csv')


pdb_sample
            



## affinity matrix

# Make matrix 77 x 77 with value 0
import numpy as np
 
matrix_0 = [[0 for x in range(77)] for x in range(77)]
matrix_array = np.array(matrix_0, dtype = np.float64)
matrix_array.shape


# change the matrix value based on affinity value(unit: nM)
for i in range(len(pdb_sample['affinity'])):
    matrix_array[i][i] = pdb_sample['affinity'][i]    
    
matrix_array
    
    
# Save affinityAvailable data
import pickle
sample_path = '/mnt/ailon/pdb_test/pdb_sample_work/make_kd_candidate/'
with open(os.path.join(sample_path + 'Y_1106.txt'), 'wb') as f:
    pickle.dump(matrix_array, f)
    
# Check the saved path
os.getcwd()


dataset_path = '/mnt/ailon/pdb_test/pdb_sample_work/make_kd_candidate/'
affinity = pickle.load(open(dataset_path + 'Y_1106.txt', 'rb'), encoding='latin1')



affinity
