// Chart co-design SA in C: alphabet a, a^-1, b, b^-1, s over F_3^4 (2 inverse pairs + 1 involution), pool of
// anyword-OK controllers (ctrl5 --dump lines of ONE base).  Objective = min over controllers of
//   #needed pairs (controller distance 1/3/4) with no controllable NB-4 word
// + w0 * (C0 uncovered points) + w2 * (C2 missing pencil points).
// Usage: cdsa pool.dump steps seed [--w0 W] [--w2 W] [--T0 T] [--T1 T] [--out chart.txt] [--init chart.txt]
//        [--top M] [--time S]
// Output chart format (text): 5 blocks "A (16 entries row-major) d (4 entries)", plus the best controller line.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define P 3
#define M 4
#define K 5
#define NW 320
#define MAXPAIRS 700
#define MAXCTRL 20000

typedef struct { int a[M][M]; } Mat;   // a[i][j] = row i, col j
typedef struct { int v[M]; } Vec;

static int symInv[K] = {1, 0, 3, 2, 4};
static int miss2[K * (K - 1)];  // C2 missing pencil points per reduced 2-word
static Mat A[K]; static Vec d[K];
static int words[NW][4];
static uint64_t rng_state = 88172645463325252ULL;
static uint64_t xr(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }
static int rnd(int n) { return (int)(xr() % (uint64_t)n); }
static double rndu(void) { return (xr() >> 11) * (1.0 / 9007199254740992.0); }

static int md(int x) { x %= P; return x < 0 ? x + P : x; }
static Vec matvec(const Mat *m, const Vec *v) { Vec r; for (int i = 0; i < M; i++) { int s = 0; for (int j = 0; j < M; j++) s += m->a[i][j] * v->v[j]; r.v[i] = md(s); } return r; }
static Mat matmul(const Mat *x, const Mat *y) { Mat r; for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { int s = 0; for (int k = 0; k < M; k++) s += x->a[i][k] * y->a[k][j]; r.a[i][j] = md(s); } return r; }
static int modinv(int a) { a = md(a); return a == 1 ? 1 : (a == 2 ? 2 : 0); }
static int mat_inv(const Mat *m, Mat *inv) {
    int W[M][2 * M];
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { W[i][j] = md(m->a[i][j]); W[i][M + j] = (i == j); }
    for (int c = 0; c < M; c++) {
        int piv = -1; for (int i = c; i < M; i++) if (W[i][c]) { piv = i; break; }
        if (piv < 0) return 0;
        if (piv != c) for (int j = 0; j < 2 * M; j++) { int t = W[c][j]; W[c][j] = W[piv][j]; W[piv][j] = t; }
        int iv = modinv(W[c][c]); for (int j = 0; j < 2 * M; j++) W[c][j] = md(W[c][j] * iv);
        for (int i = 0; i < M; i++) if (i != c && W[i][c]) { int f = W[i][c]; for (int j = 0; j < 2 * M; j++) W[i][j] = md(W[i][j] - f * W[c][j]); }
    }
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) inv->a[i][j] = W[i][M + j];
    return 1;
}
// rank of a set of column vectors
static int rank_vecs(Vec *vs, int n) {
    int W[8][M]; for (int i = 0; i < n; i++) for (int j = 0; j < M; j++) W[i][j] = vs[i].v[j];
    int rk = 0;
    for (int c = 0; c < M && rk < n; c++) {
        int piv = -1; for (int i = rk; i < n; i++) if (W[i][c]) { piv = i; break; }
        if (piv < 0) continue;
        for (int j = 0; j < M; j++) { int t = W[rk][j]; W[rk][j] = W[piv][j]; W[piv][j] = t; }
        int iv = modinv(W[rk][c]); for (int j = 0; j < M; j++) W[rk][j] = md(W[rk][j] * iv);
        for (int i = 0; i < n; i++) if (i != rk && W[i][c]) { int f = W[i][c]; for (int j = 0; j < M; j++) W[i][j] = md(W[i][j] - f * W[rk][j]); }
        rk++;
    }
    return rk;
}
static int vec_index(const Vec *v) { return v->v[0] + 3 * v->v[1] + 9 * v->v[2] + 27 * v->v[3]; }
static int is_zero(const Vec *v) { return !(v->v[0] | v->v[1] | v->v[2] | v->v[3]); }

static void build_words(void) {
    int n = 0;
    for (int w0 = 0; w0 < K; w0++) for (int w1 = 0; w1 < K; w1++) for (int w2 = 0; w2 < K; w2++) for (int w3 = 0; w3 < K; w3++) {
        if (w1 == symInv[w0] || w2 == symInv[w1] || w3 == symInv[w2]) continue;
        words[n][0] = w0; words[n][1] = w1; words[n][2] = w2; words[n][3] = w3; n++;
    }
    if (n != NW) { fprintf(stderr, "word count %d\n", n); exit(1); }
}

// ---- chart maintenance ----
static void derive_pair(int t) {   // t even: set A[t+1] = A[t]^-1, d[t+1] = A[t]^-1 d[t]
    Mat inv; if (!mat_inv(&A[t], &inv)) { fprintf(stderr, "singular pair\n"); exit(1); }
    A[t + 1] = inv; d[t + 1] = matvec(&inv, &d[t]);
}
static Vec random_nonzero(void) { Vec v; do { for (int i = 0; i < M; i++) v.v[i] = rnd(P); } while (is_zero(&v)); return v; }
static Mat random_gl(void) { Mat m, inv; do { for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) m.a[i][j] = rnd(P); } while (!mat_inv(&m, &inv)); return m; }
// random invertible P with first column = dv (so P e_0 = dv); A = P D P^-1 with D = diag(1, s1, s2, s3) is an involution fixing dv
static int inv_signs[3] = {2, 1, 1};
static double targeted_frac = 0.4;
static int freeze_pairs = 0;
static Mat Pinv_mat;  // stabiliser element used
static void random_involution(void) {
    Mat Pm, Pi;
    do {
        for (int i = 0; i < M; i++) Pm.a[i][0] = d[4].v[i];
        for (int j = 1; j < M; j++) for (int i = 0; i < M; i++) Pm.a[i][j] = rnd(P);
    } while (!mat_inv(&Pm, &Pi));
    Mat D; memset(&D, 0, sizeof D); D.a[0][0] = 1; for (int j = 1; j < M; j++) D.a[j][j] = inv_signs[j - 1];
    Mat t = matmul(&Pm, &D); A[4] = matmul(&t, &Pi);
}
static void check_chart(void) {
    Mat I; memset(&I, 0, sizeof I); for (int i = 0; i < M; i++) I.a[i][i] = 1;
    for (int t = 0; t < 4; t += 2) { Mat pr = matmul(&A[t], &A[t + 1]); if (memcmp(&pr, &I, sizeof I)) { fprintf(stderr, "pair check failed\n"); exit(1); } }
    Mat sq = matmul(&A[4], &A[4]); if (memcmp(&sq, &I, sizeof I)) { fprintf(stderr, "involution check failed\n"); exit(1); }
    Vec ad = matvec(&A[4], &d[4]); if (memcmp(&ad, &d[4], sizeof ad)) { fprintf(stderr, "involution d check failed\n"); exit(1); }
}
static void random_chart(void) {
    for (int t = 0; t < 4; t += 2) { A[t] = random_gl(); d[t] = random_nonzero(); derive_pair(t); }
    d[4] = random_nonzero(); random_involution();
}

// ---- objective pieces ----
static uint64_t badmask[5]; static int nbad;
static void compute_bad(void) {
    memset(badmask, 0, sizeof badmask); nbad = 0;
    for (int w = 0; w < NW; w++) {
        Vec cols[4];
        for (int i = 0; i < 4; i++) { Vec v = d[words[w][i]]; for (int j = i + 1; j < 4; j++) v = matvec(&A[words[w][j]], &v); cols[i] = v; }
        if (rank_vecs(cols, 4) < M) { badmask[w >> 6] |= 1ULL << (w & 63); nbad++; }
    }
}
static int c0_uncovered(void) {
    static unsigned char cov[81]; memset(cov, 0, sizeof cov);
    for (int t = 0; t < K; t++) {
        Vec Lt = d[symInv[t]]; Mat *Ati = &A[symInv[t]];
        for (int u = 0; u < K; u++) {
            Vec Lu = d[symInv[u]];
            for (int pass = 0; pass < 2; pass++) {
                Vec w;
                if (pass == 0) w = Lu; else { if (u == symInv[t]) continue; w = matvec(Ati, &Lu); }
                for (int l1 = 0; l1 < P; l1++) for (int l2 = 0; l2 < P; l2++) {
                    Vec x; for (int i = 0; i < M; i++) x.v[i] = md(l1 * Lt.v[i] + l2 * w.v[i]);
                    cov[vec_index(&x)] = 1;
                }
            }
        }
    }
    int c = 0; for (int i = 0; i < 81; i++) c += !cov[i];
    return c;
}
// annihilator functionals of a 2-dim W: two independent linear functionals vanishing on W (rows of a 2x4 matrix)
static int annihilator(Vec *W, int f[2][M]) {
    // solve f . w1 = 0, f . w2 = 0: nullspace of the 2x4 matrix with rows w1, w2
    int R[2][M]; for (int i = 0; i < 2; i++) for (int j = 0; j < M; j++) R[i][j] = W[i].v[j];
    int pivcol[2] = {-1, -1}, rk = 0;
    for (int c = 0; c < M && rk < 2; c++) {
        int piv = -1; for (int i = rk; i < 2; i++) if (R[i][c]) { piv = i; break; }
        if (piv < 0) continue;
        for (int j = 0; j < M; j++) { int t = R[rk][j]; R[rk][j] = R[piv][j]; R[piv][j] = t; }
        int iv = modinv(R[rk][c]); for (int j = 0; j < M; j++) R[rk][j] = md(R[rk][j] * iv);
        for (int i = 0; i < 2; i++) if (i != rk && R[i][c]) { int fct = R[i][c]; for (int j = 0; j < M; j++) R[i][j] = md(R[i][j] - fct * R[rk][j]); }
        pivcol[rk] = c; rk++;
    }
    if (rk != 2) return 0;
    // free columns -> basis of nullspace
    int nf = 0;
    for (int c = 0; c < M; c++) {
        if (c == pivcol[0] || c == pivcol[1]) continue;
        for (int j = 0; j < M; j++) f[nf][j] = 0;
        f[nf][c] = 1;
        for (int i = 0; i < 2; i++) f[nf][pivcol[i]] = md(-R[i][c]);
        nf++;
    }
    return nf == 2;
}
static int c2_missing(void) {
    int total = 0, idx = 0;
    for (int t = 0; t < K; t++) for (int u = 0; u < K; u++) {
        if (u == symInv[t]) continue;
        int myidx = idx++;
        Vec W[2]; W[0] = matvec(&A[u], &d[t]); W[1] = d[u];
        int f[2][M];
        if (rank_vecs(W, 2) < 2 || !annihilator(W, f)) { total += 13; miss2[myidx] = 13; continue; }  // degenerate geodesic plane: heavy penalty
        int hit[4] = {0, 0, 0, 0};
        Mat AuAt = matmul(&A[u], &A[t]);
        for (int r = 0; r < K; r++) {
            Vec Lr = d[symInv[r]];
            Vec dirs[3]; dirs[0] = Lr; dirs[1] = matvec(&A[u], &Lr); dirs[2] = matvec(&AuAt, &Lr);
            for (int i = 0; i < 3; i++) {
                int q0 = 0, q1 = 0; for (int j = 0; j < M; j++) { q0 += f[0][j] * dirs[i].v[j]; q1 += f[1][j] * dirs[i].v[j]; }
                q0 = md(q0); q1 = md(q1);
                if (!q0 && !q1) continue;
                // projective point of F_3^2: normalise first nonzero coordinate to 1
                if (q0) { int iv = modinv(q0); q1 = md(q1 * iv); hit[q1] = 1; /* (1, q1): indices 0,1,2 */ }
                else hit[3] = 1;   /* (0, 1) */
            }
        }
        int mm = 0; for (int i = 0; i < 4; i++) mm += !hit[i];
        miss2[myidx] = mm; total += mm;
    }
    return total;
}

// ---- pool ----
static int nctrl = 0;
static uint64_t (*pm)[5];      // pair masks, concatenated
static int *pw2;               // geodesic 2-word index for distance-2 pairs, -1 otherwise
static int *pstart, *pcount;   // per controller
static char **plines;
static int load_pool(const char *fn, int top) {
    FILE *f = fopen(fn, "r"); if (!f) { perror("pool"); exit(2); }
    size_t cap = 1 << 16; char *buf = malloc(cap);
    pm = malloc(sizeof(*pm) * MAXCTRL * (size_t)MAXPAIRS); pw2 = malloc(sizeof(int) * MAXCTRL * (size_t)MAXPAIRS); pstart = malloc(MAXCTRL * sizeof(int)); pcount = malloc(MAXCTRL * sizeof(int)); plines = malloc(MAXCTRL * sizeof(char *));
    long tot = 0;
    while (fgets(buf, cap, f)) {
        size_t L = strlen(buf);
        while (L == cap - 1 && buf[L - 1] != '\n') { cap *= 2; buf = realloc(buf, cap); if (!fgets(buf + L, cap - L, f)) break; L = strlen(buf); }
        char *bar = strchr(buf, '|'); if (!bar) continue;
        *bar = 0;
        plines[nctrl] = strdup(buf);
        pstart[nctrl] = (int)tot; int cnt = 0;
        char *p = bar + 1;
        while (1) {
            int c, v, dist, w2; unsigned long long m0, m1, m2, m3, m4; int used;
            if (sscanf(p, " %d %d %d %d %llx %llx %llx %llx %llx%n", &c, &v, &dist, &w2, &m0, &m1, &m2, &m3, &m4, &used) != 9) break;
            p += used;
            if (dist == 1 || dist == 3 || dist == 4 || (dist == 2 && w2 >= 0)) { pm[tot][0] = m0; pm[tot][1] = m1; pm[tot][2] = m2; pm[tot][3] = m3; pm[tot][4] = m4; pw2[tot] = dist == 2 ? w2 : -1; tot++; cnt++; }
        }
        pcount[nctrl] = cnt; nctrl++;
        if (nctrl >= MAXCTRL) break;
    }
    fclose(f); free(buf);
    if (top > 0 && top < nctrl) {
        // keep the `top` controllers with the fewest singleton pairs (stable selection by index order)
        int *sing = malloc(nctrl * sizeof(int)); int *idx = malloc(nctrl * sizeof(int));
        for (int c = 0; c < nctrl; c++) { int s = 0; for (int i = 0; i < pcount[c]; i++) { if (pw2[pstart[c] + i] >= 0) continue; int bits = 0; for (int j = 0; j < 5; j++) bits += __builtin_popcountll(pm[pstart[c] + i][j]); s += (bits == 1); } sing[c] = s; idx[c] = c; }
        for (int i = 1; i < nctrl; i++) { int t = idx[i], j = i - 1; while (j >= 0 && sing[idx[j]] > sing[t]) { idx[j + 1] = idx[j]; j--; } idx[j + 1] = t; }
        int *ns = malloc(top * sizeof(int)), *nc = malloc(top * sizeof(int)); char **nl = malloc(top * sizeof(char *));
        for (int i = 0; i < top; i++) { ns[i] = pstart[idx[i]]; nc[i] = pcount[idx[i]]; nl[i] = plines[idx[i]]; }
        memcpy(pstart, ns, top * sizeof(int)); memcpy(pcount, nc, top * sizeof(int)); memcpy(plines, nl, top * sizeof(char *));
        fprintf(stderr, "kept %d of %d controllers, singleton pairs %d .. %d\n", top, nctrl, sing[idx[0]], sing[idx[top - 1]]);
        nctrl = top;
    }
    return nctrl;
}
static int best_ctrl_idx = 0;
static int unc_list[MAXPAIRS], n_unc = 0;   // uncovered pair indices (global row index) of the controller found by pool_min
static int pool_min(int limit) {
    int best = 1 << 30; uint64_t nb[5]; for (int j = 0; j < 5; j++) nb[j] = ~badmask[j];
    for (int c = 0; c < nctrl; c++) {
        int u = 0; uint64_t (*rows)[5] = pm + pstart[c];
        int *w2 = pw2 + pstart[c];
        for (int i = 0; i < pcount[c] && u < best; i++) {
            if (w2[i] >= 0 && miss2[w2[i]] == 0) continue;
            uint64_t any = (rows[i][0] & nb[0]) | (rows[i][1] & nb[1]) | (rows[i][2] & nb[2]) | (rows[i][3] & nb[3]) | (rows[i][4] & nb[4]);
            u += (any == 0);
        }
        if (u < best) { best = u; best_ctrl_idx = c; if (best <= limit) break; }
    }
    // after the early exit, still look for an exact zero among the remaining controllers (cheap: each scan stops at
    // its first uncovered pair), so that a hit is never missed in the hot phase where `limit` is loose
    if (best > 0) {
        for (int c = best_ctrl_idx + 1; c < nctrl; c++) {
            uint64_t (*rows)[5] = pm + pstart[c]; int *w2 = pw2 + pstart[c]; int zero = 1;
            for (int i = 0; i < pcount[c]; i++) {
                if (w2[i] >= 0 && miss2[w2[i]] == 0) continue;
                uint64_t any = (rows[i][0] & nb[0]) | (rows[i][1] & nb[1]) | (rows[i][2] & nb[2]) | (rows[i][3] & nb[3]) | (rows[i][4] & nb[4]);
                if (any == 0) { zero = 0; break; }
            }
            if (zero) { best = 0; best_ctrl_idx = c; break; }
        }
    }
    // uncovered list for the chosen controller
    n_unc = 0;
    { int c = best_ctrl_idx; uint64_t (*rows)[5] = pm + pstart[c]; int *w2 = pw2 + pstart[c];
      for (int i = 0; i < pcount[c]; i++) {
          if (w2[i] >= 0 && miss2[w2[i]] == 0) continue;
          uint64_t any = (rows[i][0] & nb[0]) | (rows[i][1] & nb[1]) | (rows[i][2] & nb[2]) | (rows[i][3] & nb[3]) | (rows[i][4] & nb[4]);
          if (any == 0 && n_unc < MAXPAIRS) unc_list[n_unc++] = pstart[c] + i;
      } }
    return best;
}
// targeted move: pick an uncovered pair of the current best controller, one of its words, and move d of its first or
// last symbol out of the span of the other three controllability columns
static int targeted_move(void) {
    if (n_unc == 0) return 0;
    int row = unc_list[rnd(n_unc)];
    int ws[NW], nw = 0;
    for (int w = 0; w < NW; w++) if (pm[row][w >> 6] >> (w & 63) & 1) ws[nw++] = w;
    if (nw == 0) return 0;
    int w = ws[rnd(nw)];
    Vec cols[4];
    for (int i = 0; i < 4; i++) { Vec v = d[words[w][i]]; for (int j = i + 1; j < 4; j++) v = matvec(&A[words[w][j]], &v); cols[i] = v; }
    int pos = rnd(2) ? 3 : 0;
    int t = words[w][pos];
    Vec others[3]; int no = 0; for (int i = 0; i < 4; i++) if (i != pos) others[no++] = cols[i];
    if (rank_vecs(others, 3) < 3) return 0;
    // new column value c_new outside span(others); translate back to d_t
    Vec cnew; int tries = 0;
    do { cnew = random_nonzero(); Vec tmp[4] = {others[0], others[1], others[2], cnew}; if (rank_vecs(tmp, 4) == 4) break; } while (++tries < 100);
    if (tries >= 100) return 0;
    Vec dnew = cnew;
    if (pos == 0) { // column = A_{t4} A_{t3} A_{t2} d_{t1}: invert the product
        for (int j = 3; j >= 1; j--) { Mat inv; if (!mat_inv(&A[words[w][j]], &inv)) return 0; dnew = matvec(&inv, &dnew); }
    }
    if (t == 4) {
        // involution: conjugate A_s by Q with Q d_old = d_new
        Mat Q, Qi; Vec dold = d[4];
        int tr = 0;
        do { Q = random_gl(); // force Q d_old = d_new by adjusting: build Q with first basis vector d_old -> d_new
             Mat Bm, Bi; Vec basis[M]; int nb_ = 1; basis[0] = dold;
             for (int j = 0; j < M && nb_ < M; j++) { Vec e; memset(&e, 0, sizeof e); e.v[j] = 1; Vec tmp[M]; for (int i = 0; i < nb_; i++) tmp[i] = basis[i]; tmp[nb_] = e; if (rank_vecs(tmp, nb_ + 1) == nb_ + 1) basis[nb_++] = e; }
             for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) Bm.a[i][j] = basis[j].v[i];
             if (!mat_inv(&Bm, &Bi)) return 0;
             Mat Im; for (int i = 0; i < M; i++) Im.a[i][0] = dnew.v[i];
             for (int j = 1; j < M; j++) for (int i = 0; i < M; i++) Im.a[i][j] = rnd(P);
             Q = matmul(&Im, &Bi);
        } while (!mat_inv(&Q, &Qi) && ++tr < 50);
        if (tr >= 50) return 0;
        Mat t1 = matmul(&Q, &A[4]); A[4] = matmul(&t1, &Qi); d[4] = dnew;
    } else {
        int tp = t & ~1;   // the pair's even symbol
        if (t == tp) { d[tp] = dnew; }
        else { d[tp] = matvec(&A[tp], &dnew); }   // d_{t^-1} = A_t^-1 d_t  =>  d_t = A_t d_{t^-1}
        derive_pair(tp);
    }
    return 1;
}

static void save_chart(const char *fn, int sc[3], int ctrl_idx) {
    FILE *f = fopen(fn, "w");
    fprintf(f, "# unc=%d c0=%d c2=%d nbad=%d\n", sc[0], sc[1], sc[2], nbad);
    for (int t = 0; t < K; t++) {
        fprintf(f, "A");
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) fprintf(f, " %d", A[t].a[i][j]);
        fprintf(f, " d"); for (int i = 0; i < M; i++) fprintf(f, " %d", d[t].v[i]);
        fprintf(f, "\n");
    }
    fprintf(f, "controller %s\n", plines[ctrl_idx]);
    fclose(f);
}
static void load_chart(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror("init"); exit(2); }
    char line[4096];
    int t = 0;
    while (fgets(line, sizeof line, f) && t < K) {
        if (line[0] != 'A') continue;
        char *p = line + 1;
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { A[t].a[i][j] = (int)strtol(p, &p, 10); }
        while (*p && *p != 'd') p++; p++;
        for (int i = 0; i < M; i++) d[t].v[i] = (int)strtol(p, &p, 10);
        t++;
    }
    fclose(f);
    if (t != K) { fprintf(stderr, "init chart incomplete\n"); exit(2); }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 2; }
    const char *pool = argv[1]; long steps = atol(argv[2]); uint64_t seed = atoll(argv[3]);
    double w0 = 1.0, w2 = 1.0, T0 = 2.0, T1 = 0.1, tlimit = 1e9; const char *out = "cdsa_out.txt", *init = NULL; int top = 0;
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "--w0")) w0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--w2")) w2 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T0")) T0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T1")) T1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--out")) out = argv[++i];
        else if (!strcmp(argv[i], "--init")) init = argv[++i];
        else if (!strcmp(argv[i], "--top")) top = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--time")) tlimit = atof(argv[++i]);
        else if (!strcmp(argv[i], "--targeted")) targeted_frac = atof(argv[++i]);
        else if (!strcmp(argv[i], "--freeze-pairs")) freeze_pairs = 1;
    }
    rng_state ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    build_words();
    load_pool(pool, top);
    fprintf(stderr, "pool: %d controllers\n", nctrl);
    if (init) load_chart(init); else random_chart();
    check_chart();
    int sc[3]; compute_bad(); sc[1] = c0_uncovered(); sc[2] = c2_missing(); sc[0] = pool_min(0);
    double curv = sc[0] + w0 * sc[1] + w2 * sc[2];
    int best[3] = {sc[0], sc[1], sc[2]}; double bestv = curv; int bestc = best_ctrl_idx;
    Mat bestA[K]; Vec bestd[K]; memcpy(bestA, A, sizeof A); memcpy(bestd, d, sizeof d);
    fprintf(stderr, "init unc=%d c0=%d c2=%d nbad=%d\n", sc[0], sc[1], sc[2], nbad);
    time_t t0 = time(NULL);
    Mat saveA[K]; Vec saved[K];
    long accepted = 0;
    for (long st = 0; st < steps; st++) {
        if ((st & 1023) == 0 && difftime(time(NULL), t0) > tlimit) break;
        double Temp = T0 * pow(T1 / T0, (double)st / (steps > 1 ? steps - 1 : 1));
        memcpy(saveA, A, sizeof A); memcpy(saved, d, sizeof d);
        int sym = freeze_pairs ? 2 : rnd(3);   // 0: pair a, 1: pair b, 2: involution
        if (!freeze_pairs && targeted_frac > 0 && rndu() < targeted_frac && targeted_move()) { sym = -1; }
        if (sym < 0) { }
        else if (sym < 2) {
            int t = 2 * sym; int mv = rnd(10);
            if (mv < 4) { Mat inv; int j = rnd(M); for (int tries = 0; tries < 50; tries++) { for (int i = 0; i < M; i++) A[t].a[i][j] = rnd(P); if (mat_inv(&A[t], &inv)) break; if (tries == 49) A[t] = saveA[t]; } }
            else if (mv < 7) { int i0 = rnd(M), j0 = rnd(M), l = 1 + rnd(P - 1); if (i0 == j0) j0 = (j0 + 1) % M; for (int i = 0; i < M; i++) A[t].a[i][j0] = md(A[t].a[i][j0] + l * A[t].a[i][i0]); }
            else d[t] = random_nonzero();
            derive_pair(t);
        } else {
            int mv = rnd(10);
            if (mv < 3) { d[4] = random_nonzero(); random_involution(); }
            else if (mv < 5) { inv_signs[rnd(3)] = (rnd(2) ? 1 : 2); if (inv_signs[0] == 1 && inv_signs[1] == 1 && inv_signs[2] == 1) inv_signs[rnd(3)] = 2; random_involution(); }
            else {
                // conjugate by a random element fixing d: Q = [d | small change of a basis]
                Mat Q, Qi;
                // basis through d: d, e_j for the first independent standard vectors
                Vec basis[M]; int nb_ = 1; basis[0] = d[4];
                for (int j = 0; j < M && nb_ < M; j++) { Vec e; memset(&e, 0, sizeof e); e.v[j] = 1; Vec tmp[M]; for (int i = 0; i < nb_; i++) tmp[i] = basis[i]; tmp[nb_] = e; if (rank_vecs(tmp, nb_ + 1) == nb_ + 1) basis[nb_++] = e; }
                Mat Bm; for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) Bm.a[i][j] = basis[j].v[i];
                Mat Bi; mat_inv(&Bm, &Bi);
                Mat Im = Bm; int j = 1 + rnd(M - 1);
                do { for (int i = 0; i < M; i++) Im.a[i][j] = rnd(P); } while (!mat_inv(&Im, &Qi));
                Q = matmul(&Im, &Bi);   // Q basis_i = img_i, Q d = d
                if (!mat_inv(&Q, &Qi)) { A[4] = saveA[4]; }
                else { Mat t1 = matmul(&Q, &A[4]); A[4] = matmul(&t1, &Qi); }
            }
        }
        int nsc[3]; compute_bad(); nsc[1] = c0_uncovered(); nsc[2] = c2_missing();
        double partial = w0 * nsc[1] + w2 * nsc[2];
        // pool min with early exit relative to the acceptance threshold
        int limit = (int)floor(curv - partial + Temp * 5);
        nsc[0] = pool_min(limit < 0 ? 0 : limit);
        double v = nsc[0] + partial;
        if (v <= curv || rndu() < exp((curv - v) / Temp)) {
            curv = v; accepted++;
            if (v < bestv) {
                bestv = v; best[0] = nsc[0]; best[1] = nsc[1]; best[2] = nsc[2]; bestc = best_ctrl_idx;
                memcpy(bestA, A, sizeof A); memcpy(bestd, d, sizeof d);
                fprintf(stderr, "step %ld T=%.3f best unc=%d c0=%d c2=%d nbad=%d ctrl=%d (%.0fs)\n", st, Temp, best[0], best[1], best[2], nbad, bestc, difftime(time(NULL), t0));
                { Mat tA[K]; Vec td[K]; memcpy(tA, A, sizeof A); memcpy(td, d, sizeof d); memcpy(A, bestA, sizeof A); memcpy(d, bestd, sizeof d); save_chart(out, best, bestc); memcpy(A, tA, sizeof A); memcpy(d, td, sizeof d); }
                if (bestv == 0) break;
            }
        } else { memcpy(A, saveA, sizeof A); memcpy(d, saved, sizeof d); }
        if ((st % 20000) == 0) fprintf(stderr, "step %ld T=%.3f cur=%.2f best=%d/%d/%d acc=%ld (%.0fs)\n", st, Temp, curv, best[0], best[1], best[2], accepted, difftime(time(NULL), t0));
    }
    memcpy(A, bestA, sizeof A); memcpy(d, bestd, sizeof d); check_chart();
    compute_bad(); save_chart(out, best, bestc);
    printf("RESULT unc=%d c0=%d c2=%d nbad=%d ctrl=%d saved=%s\n", best[0], best[1], best[2], nbad, bestc, out);
    return 0;
}
