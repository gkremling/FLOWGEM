#   THIS IS AN EXAMPLE!
# Make sure all your python dependencies are installed and imported in the file.
# To install python package run, for instance: py_install("pandas").
#
# import all the dependencies below:
import torch
import torch.nn as nn
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
import torch.optim as optim
import numpy as np
from joblib import Parallel, delayed
from sklearn.model_selection import KFold
import os
import torch.nn.functional as F

# define your python imputation function



#================ nn ================#

class NPnet(nn.Module):
    def __init__(self, n, m):
        super().__init__()
        self.Wb = nn.Parameter(torch.zeros(n, m+1))

    def forward(self, x):
        return self.Wb

#================ utils ================#

def kernel_comp(x, y, sigma):
    x = x.view(x.shape[0], -1)
    y = y.view(y.shape[0], -1)
    # compute a kernel matrix
    t1 = torch.tile(torch.sum(x**2, dim=1, keepdim=True), (1, y.shape[0]))
    t2 = -2*torch.matmul(x, y.T)
    t3 = torch.tile(torch.sum(y**2, dim=1, keepdim=True).T, (x.shape[0], 1))
    return torch.exp(- (t1 + t2 + t3)/2/(sigma**2))

def comp_dist(x,y):
    # in case x and y are not vectors
    x = x.view(x.shape[0], -1)
    y = y.view(y.shape[0], -1)
    t1 = torch.tile(torch.sum(x**2, dim=1, keepdim=True), (1, y.shape[0]))
    t2 = -2*torch.matmul(x, y.T)
    t3 = torch.tile(torch.sum(y**2, dim=1, keepdim=True).T, (x.shape[0], 1))
    return t1 + t2 + t3

#================ wgf ================#

# psicon = lambda d: torch.exp(d - 1)
psicon = lambda d: 0.5*d*d + d
rbfkernel = lambda x, y, sigma: kernel_comp(x, y, sigma)


### solve systems blockwise to ensure memory usage linear in sample size (otherwise Jupyter crashes for large n)
### and regularize systems of linear equations to avoid (near-)singular system matrices A
def opt_wb_batched_lin_reg(xp, xq, x, sigma, block_size=1000, reg=True):
    """
    Memory-linear version of opt_wb_batched.
    Memory scales O(block_size * max(np, nq))
    """

    n, d = x.shape
    np = xp.shape[0]
    nq = xq.shape[0]

    device = x.device

    w = torch.empty(n, d, device=device, dtype=torch.float64)
    b = torch.empty(n, 1, device=device, dtype=torch.float64)

    # Precompute norms of xp and xq
    xp_norm = (xp**2).sum(dim=1, keepdim=True).T  # (1,np)
    xq_norm = (xq**2).sum(dim=1, keepdim=True).T  # (1,nq)

    for start in range(0, n, block_size):
        end = min(start + block_size, n)
        xb = x[start:end]                          # (bs,d)
        bs = xb.shape[0]

        xb_norm = (xb**2).sum(dim=1, keepdim=True)

        # ---------- xp kernel block ----------
        dist2_p = xb_norm - 2 * xb @ xp.T + xp_norm
        kpx = torch.exp(-dist2_p / (2 * sigma**2))      # (bs,np)

        kpx_mean = kpx.mean(dim=1, keepdim=True)
        Xpkpx = (kpx @ xp) / np

        # ---------- xq kernel block ----------
        dist2_q = xb_norm - 2 * xb @ xq.T + xq_norm
        kqx = torch.exp(-dist2_q / (2 * sigma**2))      # (bs,nq)

        kqx_mean = kqx.mean(dim=1, keepdim=True)
        Xqkqx = (kqx @ xq) / nq

        # ---------- build c ----------
        c1 = Xpkpx - Xqkqx
        c2 = kpx_mean - kqx_mean
        c = torch.cat([c1, c2], dim=1)                  # (bs,d+1)

        # ---------- build A ----------
        # A11 = 1/nq sum k(x_i,xq_j) xq_j xq_j^T
        A11 = torch.einsum('bk,kd,ke->bde', kqx, xq, xq) / nq

        A12 = Xqkqx                                     # (bs,d)
        A22 = kqx_mean.squeeze(1)                       # (bs,)

        A = torch.zeros(bs, d+1, d+1, dtype=torch.float64, device=device)
        A[:, :d, :d] = A11
        A[:, :d, d] = A12
        A[:, d, :d] = A12
        A[:, d, d] = A22

        # -------- regularize -------
        #if reg:
        #  print("reg")
        #    eps = 1e-6
        #    A += eps * torch.eye(d+1, device=A.device)

        # ---------- solve ----------
        eps = 1e-5
        A += eps * torch.eye(d+1, device=A.device)
        res = torch.linalg.solve(A, c)
        # res = torch.linalg.solve(A, c)

        w[start:end] = res[:, :d]
        b[start:end] = res[:, d:]

    return w, b


def sigma_heuristic_pairs(X, num_pairs=1000000):
    n = X.shape[0]

    i = torch.randint(0, n, (num_pairs,))
    j = torch.randint(0, n, (num_pairs,))

    dvals = torch.norm(X[i] - X[j], dim=1)

    return (0.5 * dvals.median()).sqrt().item()
    

def evaluate_sigma(xp, xq, sigma, splits_p, splits_q):
    obj = 0
    for (train_p, val_p), (train_q, val_q) in zip(splits_p, splits_q):
        w, b = opt_wb_batched_lin_reg(xp[train_p, :], xq[train_q, :], torch.cat([xp[val_p, :], xq[val_q, :]], dim=0), sigma)
        obj += torch.mean((w[:len(val_p), :] * xp[val_p, :]).sum(dim=1) + b[:len(val_p)].T).item()
        obj -= torch.mean(psicon((w[len(val_p):, :] * xq[val_q, :]).sum(dim=1) + b[len(val_p):].T)).item()

    return obj


def sample_wgf_new(X0, X_obs, M, T=500, eta=0.01, sigma_vals=None):
    Xt = X0.clone().to(dtype = torch.float64, device = device)
    M = M.clone().to(dtype = torch.float64, device = device)

    # determine all possible missingness patterns M
    unique_M, inverse = torch.unique(M, dim=0, return_inverse=True)

    # precompute observed indices for each m (improves computation efficiency)
    mask_idx = {tuple(m.tolist()): (m == 1).nonzero(as_tuple=True)[0] for m in unique_M}

    # samples drawn from p, i.e. conditional distribution X^(m)|M=m for different values of m
    # (dict with M as key and corresponding X^(m)'s as values)
    zp = {m: X_obs[inverse == i][:, mask_idx[m]].clone().to(dtype = torch.float64, device = device) for i,m in enumerate(mask_idx)}

    # weird to call this "net", it is actually just a class to store (w,b) but derived from nn.Module to enable easy SGD
    # OPTIONAL: defining it beforehand enables warm start
    # nets = {m: NPnet(Xt.shape[0], len(mask_idx[m])).to(device) for m in mask_idx}
        
    # X1_quant = torch.quantile(Xt[:, 0], 0.1).item()
    # print('t:{:3d} Quantile: {:.6f}'.format(0, X1_quant))

    # prepare cross validation for sigma
    if sigma_vals is not None:
        n_jobs = min(len(sigma_vals), os.cpu_count()-1)
        kf = KFold(n_splits=3, shuffle=True, random_state=42)
        splits_p = {m: list(kf.split(zp[m])) for m in mask_idx}
        splits_q = list(kf.split(Xt))
    
    grads = {}
    
    Xhats = []
    for t in range(T):
        print(f"Starting t = {t+1}")
        # alternatively use heuristic (doesn't yield very good results)
        # sigma_heur = (.5*comp_dist(Xt, Xt).flatten().median()).sqrt().item()
        sigma_heur = sigma_heuristic_pairs(Xt)
        #print("heuristic sigma:", sigma_heur)
        if sigma_vals is None: #and t % 10 == 0:
            sigma = sigma_heur
        else:
            print("heuristic sigma:", sigma_heur)
        for m in mask_idx:
            # choose sigma by cross validation (see App. I.1)
            if sigma_vals is not None:
                #print(f"m = {m}")
                sigma_scores = Parallel(n_jobs=n_jobs)(
                    delayed(evaluate_sigma)(zp[m], Xt[:, mask_idx[m]], sig, splits_p[m], splits_q)
                    for sig in sigma_vals
                )
                sigma = sigma_vals[np.argmax(sigma_scores)]
                #print(f"Chosen sigma: {sigma}")
            # approximately optimizes (w,b) via SGD according to the objective in equation (11)
            grads[m] = opt_wb_batched_lin_reg(zp[m], Xt[:, mask_idx[m]], Xt[:, mask_idx[m]], sigma)[0]
        
        for m in mask_idx:
            # enlarge grad to R^d
            grad_full = torch.zeros_like(Xt)
            grad_full.index_copy_(1, mask_idx[m], grads[m])
            
            # update all particles according to equation (1) 
            Xt = Xt + eta * (zp[m].shape[0]/X_obs.shape[0]) * grad_full
        
        # X1_quant = torch.quantile(Xt[:, 0], 0.1).item()
        # print('t:{:3d} Quantile: {:.6f}'.format(t+1, X1_quant))
        if sigma_vals is not None:
            print()

        Xhats.append(Xt.clone().cpu().numpy())

    return Xhats


def impute_bootstrap_per_col(X_obs, M):
    X = X_obs.clone()
    n, d = X.shape
    for j in range(d):
        observed_mask = M[:, j] == 1
        missing_mask = M[:, j] == 0

        observed_values = X[observed_mask, j]
        if observed_values.numel() == 0:
            raise ValueError(f"Value of column {j} is always missing. Need to be observed at least once.")

        num_missing = missing_mask.sum()
        if num_missing > 0:
            rand_idx = torch.randint(
                0, observed_values.shape[0],
                (num_missing,),
                device=X.device
            )
            X[missing_mask, j] = observed_values[rand_idx]
    return X


def impute_wgf_python(X_obs):
    X = torch.tensor(X_obs.values)
    M = torch.where(X.isnan(), 0, 1)

    # (1) Standardize each column using observed values only
    col_std = torch.nanmean((X - torch.nanmean(X, dim=0)) ** 2, dim=0) ** 0.5  # std per column
    col_mean = torch.nanmean(X, dim=0)

    X0 = (X - col_mean) / col_std

    X0 = impute_bootstrap_per_col(X0, M)
    Xhats = sample_wgf_new(X0, X0, M, T=150, eta=0.01, sigma_vals=None)
    
    # (2) Reverse transformation
    X_reconstructed = torch.tensor(Xhats[-1]) * col_std + col_mean

    return X_reconstructed
