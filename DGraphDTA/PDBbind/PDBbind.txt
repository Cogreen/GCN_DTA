1. PDBbind having kd affinity (sample test: 100->77; samples are reduced due to smiles format) 
2. only Kinase among PDBbind (sample test: 100)


# PDB_bind_DgraphDTA

This page is bulit to test PDB_bind sample dataset via DGraphDTA model.

## Sampling
1. Sampling(100): affinity distribution.py
2. Ligand smiles: make_SmilesfromMol.py
3. PDB sequence file: 
4. Affinity Matrix(for sampling data): make_affinity_matrix.py
5. Making Fasta files: make_fasta_file.py


## Make alignment files: *.aln
* Execute: bahs pdb_fasta_test.sh (on terminal)
- pdb_fasta_test.sh
- pdb_fasta_list_txt
- hhblits.sh
- aln_gen_v2.py
- unicluster


# Results for Making Sampling data 
> protein 
- protein_sample.json
- protein_sampleAvailable_2.json
- pdb_kd_004_sample.csv

> ligand 
- ligand_sample.json
- ligand_sampleAvailable.json

> affinity
- Y_sample.txt
- Y_sampleAvailable.txt
- Y_sampleAvailable3.txt


# Check residue_table 
residue_table.py
> output: 

# Refer
mat_to_npy.py
