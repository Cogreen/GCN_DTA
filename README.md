# GCN_DTA
## Graph Convolution Network for Drug Target Affinity

## Dependedncies
numpy == 1.17.4
kreas == 2.3.1
Pconsc4 == 0.4
pytorch == 1.3.0
PyG (torch-geometric) == 1.3.2
hhsuite (https://github.com/soedinglab/hh-suite)
rdkit == 2019.03.4.0
ccmpred (https://github.com/soedinglab/CCMpred)


### ■ Apply DGraphDTA model to the original dataset
- Davis dataset 
- KIBA dataset 

### ■ Apply DGraphDTA model to PBDbind dataset and Kinase data


### ■ Correct the model considering on Data distribution 
- re-training: applying only dataset having available affinity
- re-training2: affinity one-hot encoding depeding on whether affinity exists or not
