#!/usr/bin/env python3
"""aut_k.py: automorphisms of K commuting with alpha for a gf.py spec, block by block, written as a AUTK file for table_sweep.c
(magic 'AUTK', int32 nB, int32 nK, uint16 B[nB * nK]: each row a permutation of K in gf.py's mixed-radix indexing).

Per block the centraliser of alpha_block inside a known subgroup of Aut(block):
  cyc  Z_q          : all units u (Aut(Z_q) is abelian)
  vec  F_p^e, A     : the centraliser of A in GL(e, p) (nullspace of X -> XA - AX, invertible elements) when it has at most
                      CAP elements, else the invertible elements of F_p[A] (polynomials in A), else scalars x <A>
  heis p^{1+2}_+, A : (X, l) -> ((v, z) -> (X v, det(X) z + l.v)) with X in C(A) (as for vec) and l A = det(A) l
  frob F_p^2 x| <Bf>: (N, j, w) -> ((v, 0) -> (N v, 0), c -> (w, c^j)) over the normaliser list of <Bf> with N A = A N,
                      every w in F_p^2, kept when the permutation commutes with the block's alpha (numeric check, so twisted
                      blocks are handled too)
  ab2  Z_p^a x Z_p^b: brute force over Aut, numeric commutation check
B = the direct product of the block groups, largest blocks first while the product stays <= CAP (the other blocks
contribute the identity only), so B is a group.  Every element is checked numerically to be an automorphism of K
commuting with alpha (table_sweep.c re-checks at load).

usage: aut_k.py SPEC_JSON|@file OUT.aut [--cap 20000]
"""
import itertools
import json
import sys

import numpy as np

import gf  # noqa: E402
import enum_family as ef  # noqa: E402


def nullspace_mod_p(rows, p, nvar):
    """basis of {x : rows x = 0} over F_p (rows: list of coefficient lists)"""
    Mx = [r[:] for r in rows]
    piv_cols = []
    r = 0
    for c in range(nvar):
        piv = next((i for i in range(r, len(Mx)) if Mx[i][c] % p), None)
        if piv is None:
            continue
        Mx[r], Mx[piv] = Mx[piv], Mx[r]
        iv = pow(Mx[r][c], -1, p)
        Mx[r] = [v * iv % p for v in Mx[r]]
        for i in range(len(Mx)):
            if i != r and Mx[i][c] % p:
                f = Mx[i][c]
                Mx[i] = [(a - f * b) % p for a, b in zip(Mx[i], Mx[r])]
        piv_cols.append(c)
        r += 1
    free = [c for c in range(nvar) if c not in piv_cols]
    basis = []
    for fc in free:
        x = [0] * nvar
        x[fc] = 1
        for i, pc in enumerate(piv_cols):
            x[pc] = (-Mx[i][fc]) % p
        basis.append(x)
    return basis


def det_mod_p(X, p):
    e = len(X)
    Mx = [r[:] for r in X]
    d = 1
    for c in range(e):
        piv = next((i for i in range(c, e) if Mx[i][c] % p), None)
        if piv is None:
            return 0
        if piv != c:
            Mx[c], Mx[piv] = Mx[piv], Mx[c]
            d = -d
        d = d * Mx[c][c] % p
        iv = pow(Mx[c][c], -1, p)
        for i in range(c + 1, e):
            f = Mx[i][c] * iv % p
            Mx[i] = [(a - f * b) % p for a, b in zip(Mx[i], Mx[c])]
    return d % p


def centralizer_matrices(A, p, cap):
    """invertible X with XA = AX: full centraliser if small, else F_p[A]^*, else scalars x <A>"""
    e = len(A)
    rows = []
    for i in range(e):
        for j in range(e):
            coef = [0] * (e * e)
            # (XA)_{ij} = sum_k X_{ik} A_{kj};  (AX)_{ij} = sum_k A_{ik} X_{kj}
            for k in range(e):
                coef[i * e + k] = (coef[i * e + k] + A[k][j]) % p
                coef[k * e + j] = (coef[k * e + j] - A[i][k]) % p
            rows.append(coef)
    basis = nullspace_mod_p(rows, p, e * e)
    d = len(basis)
    mats = []
    if p ** d <= cap:
        for coeffs in itertools.product(range(p), repeat=d):
            X = [[0] * e for _ in range(e)]
            for cf, bv in zip(coeffs, basis):
                if cf:
                    for i in range(e):
                        for j in range(e):
                            X[i][j] = (X[i][j] + cf * bv[i * e + j]) % p
            if det_mod_p(X, p):
                mats.append(X)
        return mats, "full"
    # polynomials in A
    pows = [[[int(i == j) for j in range(e)] for i in range(e)]]
    while True:
        nxt = gf.matmul(pows[-1], A, p)
        # linear independence of nxt from pows: rank test
        vecs = [[v for row in P for v in row] for P in pows + [nxt]]
        if ef.rank_mod_p(vecs, p) < len(vecs):
            break
        pows.append(nxt)
    m = len(pows)
    if p ** m <= cap:
        for coeffs in itertools.product(range(p), repeat=m):
            X = [[0] * e for _ in range(e)]
            for cf, P in zip(coeffs, pows):
                if cf:
                    for i in range(e):
                        for j in range(e):
                            X[i][j] = (X[i][j] + cf * P[i][j]) % p
            if det_mod_p(X, p):
                mats.append(X)
        return mats, "poly"
    o = gf.matorder(A, p)
    seen = set()
    for s in range(1, p):
        P = [[int(i == j) for j in range(e)] for i in range(e)]
        for _ in range(o):
            X = [[s * v % p for v in row] for row in P]
            key = tuple(v for row in X for v in row)
            if key not in seen:
                seen.add(key); mats.append(X)
            P = gf.matmul(P, A, p)
    return mats, "scalar_x_powers"


def perm_vec(p, e, X):
    n = p ** e
    idx = np.arange(n)
    coords = np.stack([(idx // p ** (e - 1 - i)) % p for i in range(e)], axis=1)
    img = (coords @ np.array(X, dtype=np.int64).T) % p
    out = np.zeros(n, dtype=np.int64)
    for i in range(e):
        out = out * p + img[:, i]
    return out


def block_autos(bl, cap):
    """list of permutations (np.int64 arrays) of the block commuting with its alpha, as a group; plus a label"""
    t = bl["t"]
    nb, T, alpha, _ = gf.build_block(bl)
    alpha = alpha.astype(np.int64)
    perms = []
    label = ""
    if t == "cyc":
        q = bl["q"]
        for u in range(1, q):
            if np.gcd(u, q) == 1:
                perms.append((u * np.arange(q)) % q)
        label = f"cyc{q}:units{len(perms)}"
    elif t == "vec":
        p, A = bl["p"], bl["A"]
        mats, how = centralizer_matrices(A, p, cap)
        perms = [perm_vec(p, len(A), X) for X in mats]
        label = f"vec{p}^{len(A)}:C(A)={len(perms)}({how})"
    elif t == "heis":
        p, A = bl["p"], bl["A"]
        mats, how = centralizer_matrices(A, p, cap)
        dA = gf.det2(A, p)
        # l with l A = dA l  (row vector): solve
        rows = [[(A[0][0] - dA) % p, A[1][0] % p], [A[0][1] % p, (A[1][1] - dA) % p]]   # (lA - dA l)_j = sum_i l_i A_ij - dA l_j
        ls_basis = nullspace_mod_p(rows, p, 2)
        ls = []
        for coeffs in itertools.product(range(p), repeat=len(ls_basis)):
            l = [0, 0]
            for cf, bv in zip(coeffs, ls_basis):
                l[0] = (l[0] + cf * bv[0]) % p; l[1] = (l[1] + cf * bv[1]) % p
            ls.append(tuple(l))
        ls = sorted(set(ls))
        n3 = p ** 3
        idx = np.arange(n3); z_ = idx % p; tt = idx // p; b_ = tt % p; a_ = tt // p
        for X in mats:
            dX = gf.det2(X, p)
            ia = (X[0][0] * a_ + X[0][1] * b_) % p
            ib = (X[1][0] * a_ + X[1][1] * b_) % p
            for l in ls:
                zz = (dX * z_ + l[0] * a_ + l[1] * b_) % p
                perms.append((ia * p + ib) * p + zz)
        label = f"heis{p}:C(A)={len(mats)}({how})xl{len(ls)}"
    elif t == "frob":
        p, Bf, A, e = bl["p"], bl["B"], bl["A"], bl.get("e", 1)
        k = gf.matorder(Bf, p)
        Nrm = ef.frob_normalizer_list(p, Bf, k)
        Bp = [np.array(gf.matpow(Bf, i, p), dtype=np.int64) for i in range(k)]
        n2 = p * p * k
        idx = np.arange(n2); i_ = idx % k; tt = idx // k; b_ = tt % p; a_ = tt // p
        V = np.stack([a_, b_], axis=1)
        cnt_nj = 0
        for N, j in Nrm:
            if gf.matmul(N, A, p) != gf.matmul(A, N, p):
                continue
            cnt_nj += 1
            Nv = (V @ np.array(N, dtype=np.int64).T) % p          # N v
            # phi(v, i) = (N v, 0) (w, j)^i = (N v + S_i(w), j i),  S_i(w) = sum_{s<i} B^{js} w
            # precompute the matrices S_i = sum_{s<i} B^{js}
            S = []
            acc = np.zeros((2, 2), dtype=np.int64)
            for i in range(k):
                S.append(acc.copy())
                acc = (acc + Bp[(j * i) % k]) % p
            Smat = np.stack(S)                                      # (k, 2, 2)
            for w in itertools.product(range(p), repeat=2):
                wv = np.array(w, dtype=np.int64)
                Sw = (Smat @ wv) % p                                # (k, 2)
                va = (Nv[:, 0] + Sw[i_, 0]) % p
                vb = (Nv[:, 1] + Sw[i_, 1]) % p
                perm = (va * p + vb) * k + (j * i_) % k
                if np.array_equal(perm[alpha], alpha[perm]):
                    perms.append(perm)
        label = f"frob{p}^2x|{k}:NA=AN:{cnt_nj},commuting{len(perms)}"
    elif t == "ab2":
        p, a, b = bl["p"], bl["a"], bl["b"]
        qa, qb = p ** a, p ** b
        n_ = qa * qb
        idx = np.arange(n_); y_ = idx % qb; x_ = idx // qb
        for x1 in range(qa):
            for y1 in range(qb):
                for x2 in range(qa):
                    for y2 in range(qb):
                        if (qa * x1) % qa or (qa * y1) % qb or (qb * x2) % qa or (qb * y2) % qb:
                            continue
                        img = ((x1 * x_ + x2 * y_) % qa) * qb + (y1 * x_ + y2 * y_) % qb
                        if len(set(img.tolist())) != n_:
                            continue
                        if np.array_equal(img[alpha], alpha[img]):
                            perms.append(img)
        label = f"ab2{p}^{a},{b}:commuting{len(perms)}"
    else:
        raise ValueError(t)
    # verify: bijection, homomorphism, commuting with alpha (exact, on the block table)
    T64 = T.astype(np.int64)
    good = []
    for pm in perms:
        pm = np.asarray(pm, dtype=np.int64)
        if len(np.unique(pm)) != nb:
            continue
        if not np.array_equal(pm[alpha], alpha[pm]):
            continue
        if not np.array_equal(pm[T64], T64[pm][:, pm]):
            continue
        good.append(pm)
    assert len(good) == len(perms), f"{label}: {len(perms) - len(good)} candidates failed the automorphism check"
    return good, label


def main():
    spec = gf.load_spec(sys.argv[1])
    outfn = sys.argv[2]
    cap = int(sys.argv[sys.argv.index("--cap") + 1]) if "--cap" in sys.argv else 20000
    blocks = spec["blocks"]
    built = [gf.build_block(bl) for bl in blocks]
    sizes = [b[0] for b in built]
    nK = int(np.prod(sizes))
    groups = []
    labels = []
    for bl in blocks:
        g, lab = block_autos(bl, cap)
        groups.append(g); labels.append(lab)
    # keep whole block groups, largest first, while the product stays <= cap
    order = sorted(range(len(blocks)), key=lambda i: -len(groups[i]))
    use = [False] * len(blocks)
    prod = 1
    for i in order:
        if prod * len(groups[i]) <= cap:
            use[i] = True; prod *= len(groups[i])
    chosen = [groups[i] if use[i] else [np.arange(sizes[i])] for i in range(len(blocks))]
    # mixed radix: first block most significant
    idx = np.arange(nK)
    coords = []
    rem = idx.copy()
    for s in reversed(sizes):
        coords.append(rem % s); rem //= s
    coords = coords[::-1]
    out = np.zeros((prod, nK), dtype=np.uint16)
    r = 0
    for combo in itertools.product(*chosen):
        img = np.zeros(nK, dtype=np.int64)
        for pm, c, s in zip(combo, coords, sizes):
            img = img * s + pm[c]
        out[r] = img; r += 1
    assert r == prod
    with open(outfn, "wb") as f:
        f.write(b"AUTK")
        f.write(np.array([prod, nK], dtype=np.int32).tobytes())
        f.write(out.tobytes())
    print(json.dumps({"nK": nK, "nB": prod, "blocks": labels, "used": use}))


if __name__ == "__main__":
    main()
