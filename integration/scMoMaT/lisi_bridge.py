import numpy as np
import pandas as pd
import scanpy as sc
from datetime import datetime
from scipy.sparse.csgraph import connected_components
from sklearn.metrics.cluster import silhouette_samples, silhouette_score
from scib.metrics import lisi
import torch
import math
from sklearn.neighbors import NearestNeighbors
import scipy
import gc

def print_results(results, exclude_keys=['df_lisi']):
    for k, v in results.items():
        if k not in exclude_keys:
            print(f'{k}={v:.5f}')

def batch_gpu_pairdist(emb1, emb2, batch_size=1024):
    n_batch = math.ceil(emb2.shape[0] / batch_size)
    emb2_gpu = torch.FloatTensor(emb2).cuda()
    emb2_gpu = emb2_gpu / torch.linalg.norm(emb2_gpu, ord=2, dim=1, keepdim=True)
    
    st = 0
    dist = []
    for i in range(n_batch):
        bsz = min(batch_size, emb1.shape[0] - i*batch_size)
        emb1_batch_gpu = torch.FloatTensor(emb1[st:st+bsz]).cuda()
        emb1_batch_gpu /= torch.linalg.norm(emb1_batch_gpu, ord=2, dim=1, keepdim=True)
        
        _ = -emb1_batch_gpu @ emb2_gpu.T  # 0-similarity => dist
        dist.append(_.cpu().numpy())
        st = st+bsz
        
        del emb1_batch_gpu
        torch.cuda.empty_cache()
        gc.collect()
    
    del emb2_gpu
    torch.cuda.empty_cache()
    gc.collect()
    
    dist = np.vstack(dist)
    return dist

def eval_FOSCTTM(adata1, adata2, use_rep='X_emb', return_dist=False):
    dist = batch_gpu_pairdist(adata1.obsm[use_rep], adata2.obsm[use_rep], batch_size=2048)
    foscttm_x = (dist < dist.diagonal().reshape(-1, 1)).mean(axis=1)
    foscttm_y = (dist < dist.diagonal()).mean(axis=0)
    foscttm = (foscttm_x+foscttm_y).mean()/2

    return foscttm, dist

def snn_scores(
        x, y, k=1
):
    '''
        return: matching score matrix
    '''

    # print(f'{k} neighbors to consider during matching')

    x = x / np.linalg.norm(x, axis=1, keepdims=True)
    y = y / np.linalg.norm(y, axis=1, keepdims=True)

    ky = k or min(round(0.01 * y.shape[0]), 1000)   
    nny = NearestNeighbors(n_neighbors=ky).fit(y)
    x2y = nny.kneighbors_graph(x)
    y2y = nny.kneighbors_graph(y)

    kx = k or min(round(0.01 * x.shape[0]), 1000)
    nnx = NearestNeighbors(n_neighbors=kx).fit(x)
    y2x = nnx.kneighbors_graph(y)
    x2x = nnx.kneighbors_graph(x)

    x2y_intersection = x2y @ y2y.T
    y2x_intersection = y2x @ x2x.T
    jaccard = x2y_intersection + y2x_intersection.T
    jaccard.data = jaccard.data / (2 * kx + 2 * ky - jaccard.data)
    matching_matrix = jaccard.multiply(1 / jaccard.sum(axis=1)).tocsr()
    return matching_matrix

def eval_matching_score(
        mod1, mod2, split_by='batch', k=1, use_rep='X'
):  
    '''
        return: scipy.sparse.csr_matrix
    '''

    mod1_splits = set(mod1.obs[split_by])
    mod2_splits = set(mod2.obs[split_by])
    splits = mod1_splits | mod2_splits
    
    matching_matrices, mod1_obs_names, mod2_obs_names = [], [], []
    for split in splits:
        mod1_split = mod1[mod1.obs[split_by] == split]
        mod2_split = mod2[mod2.obs[split_by] == split]
        mod1_obs_names.append(mod1_split.obs_names)
        mod2_obs_names.append(mod2_split.obs_names)
        
        matching_matrices.append(
            snn_scores(mod1_split.X, mod2_split.X, k)
            if use_rep=='X' else
            snn_scores(mod1_split.obsm[use_rep], mod2_split.obsm[use_rep], k)
        )
        
    mod1_obs_names = pd.Index(np.concatenate(mod1_obs_names))
    mod2_obs_names = pd.Index(np.concatenate(mod2_obs_names))
    combined_matrix = scipy.sparse.block_diag(matching_matrices, format="csr")
    score_matrix = combined_matrix[
        mod1_obs_names.get_indexer(mod1.obs_names), :
    ][
        :, mod2_obs_names.get_indexer(mod2.obs_names)
    ]

    score = (score_matrix.diagonal() / score_matrix.sum(axis=1).A1).mean()
    return score

def eval_lisi(
        adata,
        batch_keys=['domain', 'batch'],
        use_rep='X_emb', use_neighbors=False,
    ):
    res = {}
    for key in batch_keys:
        adata.obs[key] = adata.obs[key].astype('category')

        _lisi = lisi.ilisi_graph(
            adata,
            key,
            'embed' if not use_neighbors else 'knn',
            use_rep=use_rep,
            k0=90,
            subsample=None,
            scale=True,
            n_cores=1,
            verbose=False,
        )
        res[key] = _lisi
    df = pd.DataFrame.from_dict(res, orient='index').T
    df.columns = [_+'_LISI' for _ in df.columns]
    return df

def eval_bridge(
        adata1, adata2,
        label_key='cell_type',
        batch_key='batch',
        use_rep='X_emb',
        use_fosc=True, use_acc=False, use_score=True,
    ):
    # foscttm
    fosc, _dist = eval_FOSCTTM(adata1, adata2, use_rep=use_rep, return_dist=True)
    results = {'FOSCTTM': fosc}

    if use_score:
        score = eval_matching_score(adata1, adata2, split_by=batch_key, k=1, use_rep=use_rep)
        results.update({'Match_score':score})

    print_results(results)
    return results