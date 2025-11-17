import MarsGT 
from MarsGT.conv import *
from MarsGT.egrn import *
from MarsGT.marsgt_model import *
from MarsGT.utils import *
import anndata as ad
from collections import Counter
import copy
import dill
from functools import partial
import json
import math
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
from operator import itemgetter
import random
import scipy.sparse as sp
from scipy.io import mmread
from scipy.sparse import hstack, vstack, coo_matrix
import seaborn as sb
from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.decomposition import IncrementalPCA
from sklearn.decomposition import SparsePCA
from sklearn.metrics import accuracy_score
from sklearn.metrics.cluster import normalized_mutual_info_score
import time
import torch
import torch.cuda as cuda
from torch import nn
from torch.autograd import Variable
import torch.distributions as D
import torch.nn.functional as F
import torch_geometric.data as Data
from torch_geometric.nn import GCNConv, GATConv
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.nn.inits import glorot, uniform
from torch_geometric.utils import softmax as Softmax
from torchmetrics.functional import pairwise_cosine_similarity
import warnings
from warnings import filterwarnings
import xlwt
import argparse
from tqdm import tqdm
import scanpy as sc
from scipy import sparse

from os.path import exists


import pkg_resources
marsgt_version = pkg_resources.get_distribution("MarsGT").version
print('Please verify with the official website whether you are using the latest version.')
print('You are using version :',marsgt_version)


filterwarnings("ignore")
seed = 0
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
os.environ['PYTHONHASHSEED'] = str(seed)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

parser = argparse.ArgumentParser(description='Training GNN on gene cell graph')
parser.add_argument('--fi', type=int, default=0) # This parameter is used for the benchmark to specify the starting sequence number of the created files
parser.add_argument('--labsm', type=float, default=0.1) # The rate of LabelSmoothing
parser.add_argument('--wd', type=float, default=0.1) # The 'weight_decay' parameter is used to specify the strength of L2 regularization
parser.add_argument('--lr', type=float, default=0.0005) # learning rate
parser.add_argument('--n_hid', type=int, default=104) # The number of layers should be a multiple of 'n_head' in order to make any modifications
parser.add_argument('--nheads', type=int, default=8) # The 'heads' parameter represents the number of attention heads in the attention mechanism
parser.add_argument('--nlayers', type=int, default=3) # The number of layers in network
parser.add_argument('--cell_size', type=int, default=30) # The number of cells per subgraph (batch)
parser.add_argument('--neighbor', type=int, default=20) # The number of neighboring nodes to be selected for each cell in the subgraph
parser.add_argument('--egrn', type=bool, default=False) # Whether to output the Enhancer-Gene regulatory network
parser.add_argument('--epochs', type=int, default=3) # The epoch number of NodeDimensionReduction
parser.add_argument('--num_epochs', type=int, default=3) # The epoch number of MarsGT-Model
parser.add_argument('--output_file', type=str, default='Tutorial_example\output') # Please choose an output path to replace this path on your own.
args = parser.parse_args([])

output_file = args.output_file
fi=args.fi
labsm = args.labsm
lr = args.lr
wd = args.wd
n_hid = args.n_hid
nheads = args.nheads
nlayers = args.nlayers
cell_size = args.cell_size
neighbor = args.neighbor
egrn = args.egrn
epochs = args.epochs
num_epochs = args.num_epochs

path = '/data/Home/fabotao/Projects/Rare_cell/data/'

experiments = ["Mouse_colon_RNA_ATAC","Mouse_gdT_RNA_ATAC","Human_retina_rpe_choroid_RNA_ATAC","Human_Gray_matter_RNA_ATAC"]
labsm_all = [0, 0.1, 0.3]
wd_all = [0, 0.1, 0.3]
lr_all = [0.001, 0.0005]


for experiment in experiments:
  print(experiment)
  gene_cell = ad.read_h5ad(path + experiment + '/Gene_Cell.h5ad')
  peak_cell = ad.read_h5ad(path + experiment + '/Peak_Cell.h5ad')
  gene_peak = ad.read_h5ad(path + experiment + '/Gene_Peak.h5ad')
    
  gene_names = pd.read_csv(path + experiment + '/Gene_names.tsv', sep='\t', header=None)
  peak_names = pd.read_csv(path + experiment + '/Peak_names.tsv', sep='\t', header=None)

  RNA_matrix = gene_cell.X
  ATAC_matrix = peak_cell.X
  RP_matrix = gene_peak.X
  Gene_Peak = gene_peak.X
    
  Gene_Peak.obs_names = gene_peak.obs_names
  Gene_Peak.var_names = gene_peak.var_names
  ATAC_matrix.obs_names = peak_cell.obs_names
  ATAC_matrix.var_names = peak_cell.var_names
  RNA_matrix.obs_names = gene_cell.obs_names
  RNA_matrix.var_names = gene_cell.var_names

  cell_num = RNA_matrix.shape[1]
  gene_num = RNA_matrix.shape[0]
  peak_num = ATAC_matrix.shape[0]
  
  device = torch.device("cuda" if cuda.is_available() else "cpu")
  print('You will use : ',device)
  # clustering result by scanpy
  initial_pre = initial_clustering(RNA_matrix) 
  # number of every cluster
  cluster_ini_num = len(set(initial_pre)) 
  ini_p1 = [int(i) for i in initial_pre] 
  # partite the data into batches
  indices, Node_Ids, dic = batch_select_whole(RNA_matrix, ATAC_matrix, neighbor = [neighbor], cell_size=cell_size)
  n_batch = len(indices)

  for labsm in labsm_all:
    print(labsm)
    for lr in lr_all:
      print(lr)
      for wd in wd_all:
            print(wd)
            postfix = '_labsm_'+ str(labsm) + '_lr_' + str(lr) + '_wd_' + str(wd)
            file_path = path + experiment + '/pred_GPU' + postfix +'.npy'
            if(exists(file_path)):
              continue
    
            if __name__ == "__main__":

                # Reduce the dimensionality of features for cell, gene, and peak data.
                node_model = NodeDimensionReduction(RNA_matrix, ATAC_matrix, indices, ini_p1, n_hid=n_hid, n_heads=nheads,n_layers=nlayers,labsm=labsm, lr=lr, wd=wd, device=device, num_types=3, num_relations=2, epochs=10)
                gnn,cell_emb,gene_emb,peak_emb,h = node_model.train_model(n_batch=n_batch)
    
                # Instantiate the MarsGT_model
                MarsGT_model = MarsGT(gnn=gnn, h=h, labsm=labsm, n_hid=n_hid, n_batch=n_batch, device=device,lr=lr,wd=wd, num_epochs=10)
                # Train the model
                MarsGT_gnn = MarsGT_model.train_model(indices=indices,RNA_matrix=RNA_matrix, ATAC_matrix=ATAC_matrix, Gene_Peak=Gene_Peak, ini_p1=ini_p1)
                # The result of MarsGT
                MarsGT_result = MarsGT_pred(RNA_matrix, ATAC_matrix, RP_matrix, egrn=egrn, MarsGT_gnn=MarsGT_gnn, indices=indices, nodes_id=Node_Ids, cell_size=cell_size, device=device, gene_names=gene_names, peak_names=peak_names)
        
                # Save numpy arrays to files
                np.save(path + experiment + "/Node_Ids_GPU" + postfix +".npy", Node_Ids)
                np.save(path + experiment + "/pred_GPU"+ postfix +".npy", MarsGT_result['pred_label'])
                np.save(path + experiment + "/cell_embedding_GPU"+ postfix +".npy", MarsGT_result['cell_embedding'])
        
                del MarsGT_result, MarsGT_gnn, MarsGT_model, gnn
  del gene_cell, peak_cell, gene_peak, RNA_matrix, ATAC_matrix, RP_matrix, Gene_Peak, Node_Ids

 
   



## Output npy to txt

path = '/data/Home/fabotao/Projects/Rare_cell/data/'

experiments = ["Mouse_colon_RNA_ATAC","Mouse_gdT_RNA_ATAC","Human_retina_rpe_choroid_RNA_ATAC","Human_Gray_matter_RNA_ATAC"]
labsm_all = [0, 0.1, 0.3]
wd_all = [0, 0.1, 0.3]
lr_all = [0.001, 0.0005]


for experiment in experiments:
  print(experiment)
  for labsm in labsm_all:
    print(labsm)
    for lr in lr_all:
      print(lr)
      for wd in wd_all:
            print(wd)
            postfix = '_labsm_'+ str(labsm) + '_lr_' + str(lr) + '_wd_' + str(wd)
            pred_file_path = path + experiment + '/pred_GPU' + postfix +'.npy'
            Id_file_path = path + experiment + '/Node_Ids_GPU' + postfix +'.npy'
            pred = np.load(pred_file_path)
            pred_df = pd.DataFrame(pred)
            pred_df.to_csv(path + experiment + '/pred_GPU' + postfix +'.csv', index=False, header=False)
            id = np.load(Id_file_path)
            id_df = pd.DataFrame(id)
            id_df.to_csv(path + experiment + '/Node_Ids_GPU' + postfix +'.csv', index=False, header=False)







