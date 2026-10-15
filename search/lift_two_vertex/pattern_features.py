"""Pattern-level features of shape-103 lifts over Z_N x|_r Z_M, computed on the 2M-vertex bi-circulant quotient Q
(the lift is the Z_N-cover of Q): Q has (0,y) ~ (0, y +- i), (1, y + l_k) and (1,y) ~ (1, y +- j), (0, y - l_k) for the
loop levels i, j and link levels l_k.  features() returns the quotient girth (shortest closed non-backtracking walk),
the distance spheres of Q from both roots, and the score ratio (sum over the 2M cells of N e^(-W_cell/N), W_cell = the
number of words of length <= D from the roots landing in the cell, divided by the uniform value 2 n e^(-Moore/n))."""
import math
from collections import Counter

import numpy as np


def moore(d, D):
    return 1 + sum(d * (d - 1) ** (l - 1) for l in range(1, D + 1))


def quotient_adj(M, i, j, L):
    """adjacency (2M x 5) of Q; dart columns: loop+, loop-, link1, link2, link3"""
    y = np.arange(M)
    f0 = np.stack([(y + i) % M, (y - i) % M] + [M + (y + l) % M for l in L], axis=1)
    f1 = np.stack([M + (y + j) % M, M + (y - j) % M] + [(y - l) % M for l in L], axis=1)
    return np.concatenate([f0, f1], axis=0)


def nb_walks_quotient(adj, root, D):
    """non-backtracking walk counts per quotient vertex and length, by darts; backtracking = the inverse dart
    (loop+ <-> loop-, link k <-> link k back).  Returns W (D+1, 2M)."""
    n2, d = adj.shape
    M = n2 // 2
    # inverse dart of (u, k): loop+ <-> loop- (k 0 <-> 1), link k (k>=2) <-> link k of the other fibre
    W = np.zeros((D + 1, n2), dtype=np.int64)
    W[0, root] = 1
    cur = np.zeros((n2, d), dtype=np.int64)
    cur[root, :] = 1
    W[1] = np.bincount(adj[root], minlength=n2)
    inv_col = np.array([1, 0, 2, 3, 4])
    for l in range(2, D + 1):
        atv = np.bincount(adj.ravel(), weights=cur.ravel(), minlength=n2).astype(np.int64)
        # new[v, k] = atv[v] - (walks that arrived at v by the dart whose reverse is (v, k)) = cur[u, inv_col[k]] where adj[v,k] = u
        new = np.empty_like(cur)
        for k in range(d):
            u = adj[:, k]
            new[:, k] = atv - cur[u, inv_col[k]]
        cur = new
        W[l] = np.bincount(adj.ravel(), weights=cur.ravel(), minlength=n2).astype(np.int64)
    return W


def bfs_q(adj, root):
    n2 = adj.shape[0]
    dist = np.full(n2, -1)
    dist[root] = 0
    fr = [root]
    k = 0
    while fr:
        k += 1
        nxt = []
        for u in fr:
            for v in adj[u]:
                if dist[v] < 0:
                    dist[v] = k
                    nxt.append(int(v))
        fr = nxt
    return dist


def features(N, M, D, i, j, L):
    adj = quotient_adj(M, i, j, L)
    n = 2 * N * M
    lam = moore(5, D) / n
    uni = n * math.exp(-lam)
    out = dict(i=i, j=j, L=tuple(L))
    scores = []
    sph = []
    cgirth = 99
    minmass = []
    for root in (0, M):
        W = nb_walks_quotient(adj, root, D)
        masses = W.sum(axis=0)  # words of length <= D per quotient vertex = per (fibre, layer) cell
        scores.append(float(np.sum(N * np.exp(-masses / N))))
        minmass.append(float(masses.min() / masses.mean()))
        dist = bfs_q(adj, root)
        cnt = Counter(dist.tolist())
        sph.append((cnt.get(1, 0), cnt.get(2, 0), cnt.get(3, 0), int(dist.max())))
        for l in range(1, D + 1):
            if W[l, root] > 0:
                cgirth = min(cgirth, l)
                break
    out["score"] = sum(scores)
    out["score_ratio"] = sum(scores) / (2 * uni)
    out["sph"] = sph
    out["full2"] = all(s[0] == 5 and s[1] == 20 for s in sph)
    out["full1"] = all(s[0] == 5 for s in sph)
    out["qdiam"] = max(s[3] for s in sph)
    out["qgirth"] = cgirth
    out["minmass"] = min(minmass)
    return out
