// exlift.c: exact route-chart lift evaluator and chart SA for SMALL controllers with general alphabets.
//
// Controller file (text):
//   N K
//   inv[0..K-1]          symbol inverses: pairs are consecutive (t, t+1) with inv[t] = t+1; self-inverse: inv[t] = t
//   dim[0..K-1]          translation dimension of the symbol: 1 = line symbol (p neighbours), 0 = plain symbol (1 neighbour)
//   N lines: nbr[v][0..K-1]   target of symbol t at v (must satisfy nbr[nbr[v][t]][inv[t]] = v)
//   R r_1 ... r_R        root controller vertices (deck-group orbit representatives); R = 0 means all N
// Fibre F_p^m (p = 3 by default; vectors x = sum x_i p^i).  Symbol t: (v, x) ~ (nbr[v][t], M_t x + T_t) with
//   T_t = {lambda d_t : lambda in F_p} (dim 1) or {c_t} (dim 0).
// Chart constraints (so that the lift is an undirected graph):
//   pair (t, t+1):   M_{t+1} = M_t^-1,  d_{t+1} = M_t^-1 d_t   (dim 1)   /   c_{t+1} = -M_t^-1 c_t   (dim 0)
//   self-inverse t:  M_t^2 = I,  M_t d_t = d_t                 (dim 1)   /   M_t c_t = -c_t           (dim 0)
// Objective: number of ordered pairs ((r, x), u) with lift distance > D, roots = the R base vertices x all x in F_p^m
//   (bit-parallel BFS over all roots at once).  0 = the lift has diameter <= D from every root; with the deck group
//   acting transitively on the fibres over each orbit this is the diameter of the whole lift.
// Usage: exlift ctrl.txt [--p P] [--m M] [--D D] [--steps S] [--time SEC] [--seed X] [--T0 t] [--T1 t]
//                        [--init chart.txt] [--out chart.txt] [--edges FILE] [--report] [--eval]
//   --eval: evaluate the --init chart only (prints objective, degree checks, distance distribution from root 0)
//   --edges FILE: write the sorted edge list of the final (best) lift when its objective is 0 (or always with --eval)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MAXM 8
#define MAXK 8
#define MAXF 729
#define MAXN 512
// dim 2 (pair symbols only): T_t = span(d_t, d2_t), p^2 neighbours per symbol (URSA's q = 4 regime over F_2^8)

static int P = 3, Mdim = 4, F = 81, D = 4, pneg = 0;   // pneg: plain self-inverse symbols frozen at M = -I (c free)
static int N, K, inv_[MAXK], dim_[MAXK], nbr[MAXN][MAXK], R, roots[MAXN];
static int sym_free[MAXK];  // 1 if the symbol is the representative of its class (even index of a pair, or self-inverse)
typedef struct { int a[MAXM][MAXM]; } Mat;   // a[i][j] row i col j, action y_i = sum_j a[i][j] x_j
typedef struct { int v[MAXM]; } Vec;
static Mat Mt[MAXK]; static Vec dt[MAXK], dt2[MAXK];      // chart: matrix and direction (dim 1; dim 2: dt, dt2) / constant (dim 0)
static uint64_t rng_state = 88172645463325252ULL;
static uint64_t xr(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }
static int rnd(int n) { return (int)(xr() % (uint64_t)n); }
static double rndu(void) { return (xr() >> 11) * (1.0 / 9007199254740992.0); }
static int md(int x) { x %= P; return x < 0 ? x + P : x; }
static int modinv(int a) { a = md(a); for (int i = 1; i < P; i++) if ((a * i) % P == 1) return i; return 0; }
static Vec matvec(const Mat *m, const Vec *v) { Vec r; for (int i = 0; i < Mdim; i++) { int s = 0; for (int j = 0; j < Mdim; j++) s += m->a[i][j] * v->v[j]; r.v[i] = md(s); } return r; }
static Mat matmul(const Mat *x, const Mat *y) { Mat r; memset(&r, 0, sizeof r); for (int i = 0; i < Mdim; i++) for (int j = 0; j < Mdim; j++) { int s = 0; for (int k = 0; k < Mdim; k++) s += x->a[i][k] * y->a[k][j]; r.a[i][j] = md(s); } return r; }
static int mat_inv(const Mat *m, Mat *inv) {
    int W[MAXM][2 * MAXM];
    for (int i = 0; i < Mdim; i++) for (int j = 0; j < Mdim; j++) { W[i][j] = md(m->a[i][j]); W[i][Mdim + j] = (i == j); }
    for (int c = 0; c < Mdim; c++) {
        int piv = -1; for (int i = c; i < Mdim; i++) if (W[i][c]) { piv = i; break; }
        if (piv < 0) return 0;
        if (piv != c) for (int j = 0; j < 2 * Mdim; j++) { int t = W[c][j]; W[c][j] = W[piv][j]; W[piv][j] = t; }
        int iv = modinv(W[c][c]); for (int j = 0; j < 2 * Mdim; j++) W[c][j] = md(W[c][j] * iv);
        for (int i = 0; i < Mdim; i++) if (i != c && W[i][c]) { int f = W[i][c]; for (int j = 0; j < 2 * Mdim; j++) W[i][j] = md(W[i][j] - f * W[c][j]); }
    }
    for (int i = 0; i < Mdim; i++) for (int j = 0; j < Mdim; j++) inv->a[i][j] = W[i][Mdim + j];
    return 1;
}
static int is_zero(const Vec *v) { for (int i = 0; i < Mdim; i++) if (v->v[i]) return 0; return 1; }
static Vec random_vec(int nonzero) { Vec v; do { for (int i = 0; i < Mdim; i++) v.v[i] = rnd(P); } while (nonzero && is_zero(&v)); return v; }
static Mat random_gl(void) { Mat m, inv; do { for (int i = 0; i < Mdim; i++) for (int j = 0; j < Mdim; j++) m.a[i][j] = rnd(P); } while (!mat_inv(&m, &inv)); return m; }
static Mat identity(void) { Mat I; memset(&I, 0, sizeof I); for (int i = 0; i < Mdim; i++) I.a[i][i] = 1; return I; }
// random invertible Pm whose first column is d (if d != 0), returns Pm and its inverse
static void random_frame(const Vec *d, Mat *Pm, Mat *Pi) {
    do {
        for (int j = 0; j < Mdim; j++) for (int i = 0; i < Mdim; i++) Pm->a[i][j] = rnd(P);
        if (d && !is_zero(d)) for (int i = 0; i < Mdim; i++) Pm->a[i][0] = d->v[i];
    } while (!mat_inv(Pm, Pi));
}
// random exact involution: M = Pm diag(1, s_1..s_{m-1}) Pm^-1 with Pm e_0 = d (so M d = d); signs random, not all +1 unless allow_id
static Mat random_involution(const Vec *d, int fix_first, int allow_trivial) {
    Mat Pm, Pi, Dg; random_frame(fix_first ? d : NULL, &Pm, &Pi);
    memset(&Dg, 0, sizeof Dg);
    int all1 = 1, allm = 1;
    for (int j = 0; j < Mdim; j++) { int s = (j == 0 && fix_first) ? 1 : (rnd(2) ? 1 : P - 1); Dg.a[j][j] = s; if (s != 1) all1 = 0; if (s != P - 1) allm = 0; }
    if (!allow_trivial && (all1 || allm)) { int j = fix_first ? 1 + rnd(Mdim - 1) : rnd(Mdim); Dg.a[j][j] = all1 ? P - 1 : 1; }
    Mat t = matmul(&Pm, &Dg); return matmul(&t, &Pi);
}
// conjugate M by a random element Q with Q d = d (Q = frame through d with one other image changed)
static void conjugate_fixing(Mat *M, const Vec *d) {
    Mat Bm, Bi, Im, Qi, Q;
    // basis through d: d then standard vectors
    Vec basis[MAXM]; int nb = 0;
    if (d && !is_zero(d)) basis[nb++] = *d;
    for (int j = 0; j < Mdim && nb < Mdim; j++) {
        Vec e; memset(&e, 0, sizeof e); e.v[j] = 1;
        for (int i = 0; i < Mdim; i++) for (int k = 0; k < nb; k++) Bm.a[i][k] = basis[k].v[i];
        for (int i = 0; i < Mdim; i++) Bm.a[i][nb] = e.v[i];
        // test independence: rank of first nb+1 columns -> use a quick elimination on a copy
        int W[MAXM][MAXM]; for (int i = 0; i < Mdim; i++) for (int k = 0; k <= nb; k++) W[i][k] = Bm.a[i][k];
        int rk = 0; for (int c = 0; c <= nb && rk < Mdim; c++) { int piv = -1; for (int i = rk; i < Mdim; i++) if (W[i][c]) { piv = i; break; } if (piv < 0) continue; for (int k = 0; k <= nb; k++) { int tt = W[rk][k]; W[rk][k] = W[piv][k]; W[piv][k] = tt; } int iv = modinv(W[rk][c]); for (int k = 0; k <= nb; k++) W[rk][k] = md(W[rk][k] * iv); for (int i = 0; i < Mdim; i++) if (i != rk && W[i][c]) { int f = W[i][c]; for (int k = 0; k <= nb; k++) W[i][k] = md(W[i][k] - f * W[rk][k]); } rk++; }
        if (rk == nb + 1) basis[nb++] = e;
    }
    for (int i = 0; i < Mdim; i++) for (int k = 0; k < Mdim; k++) Bm.a[i][k] = basis[k].v[i];
    if (!mat_inv(&Bm, &Bi)) return;
    int first = (d && !is_zero(d)) ? 1 : 0;
    int tries = 0;
    do { Im = Bm; int j = first + rnd(Mdim - first); for (int i = 0; i < Mdim; i++) Im.a[i][j] = rnd(P); } while (!mat_inv(&Im, &Qi) && ++tries < 50);
    if (tries >= 50) return;
    Q = matmul(&Im, &Bi); if (!mat_inv(&Q, &Qi)) return;
    Mat t1 = matmul(&Q, M); *M = matmul(&t1, &Qi);
}
static void derive_pair(int t) {
    Mat inv; if (!mat_inv(&Mt[t], &inv)) { fprintf(stderr, "singular pair matrix\n"); exit(1); }
    Mt[t + 1] = inv; Vec w = matvec(&inv, &dt[t]);
    if (dim_[t] >= 1) dt[t + 1] = w; else { for (int i = 0; i < Mdim; i++) dt[t + 1].v[i] = md(-w.v[i]); }
    if (dim_[t] == 2) dt2[t + 1] = matvec(&inv, &dt2[t]);
}
static int indep2(const Vec *a, const Vec *b) {   // a, b linearly independent?
    if (is_zero(a) || is_zero(b)) return 0;
    for (int l = 1; l < P; l++) { int eq = 1; for (int i = 0; i < Mdim; i++) if (md(a->v[i] - l * b->v[i])) { eq = 0; break; } if (eq) return 0; }
    return 1;
}
// char-2 involution fixing d (optional): I + N with N = sum u_i f_i^T, f_i(u_j) = 0, f_i(d) = 0, rank 1..3
static Mat random_involution2(const Vec *d) {
    for (int tries = 0; tries < 10000; tries++) {
        int r = 1 + rnd(3); Vec u[3], f[3]; int ok = 1;
        for (int i = 0; i < r; i++) { u[i] = random_vec(1); f[i] = random_vec(1); }
        for (int i = 0; i < r && ok; i++) {
            if (d) { int s = 0; for (int q = 0; q < Mdim; q++) s += f[i].v[q] * d->v[q]; if (s & 1) ok = 0; }
            for (int j = 0; j < r && ok; j++) { int s = 0; for (int q = 0; q < Mdim; q++) s += f[i].v[q] * u[j].v[q]; if (s & 1) ok = 0; }
        }
        if (!ok) continue;
        Mat Nm; memset(&Nm, 0, sizeof Nm); int nz = 0;
        for (int i = 0; i < r; i++) for (int a = 0; a < Mdim; a++) for (int b = 0; b < Mdim; b++) { Nm.a[a][b] ^= (u[i].v[a] & f[i].v[b]); }
        for (int a = 0; a < Mdim; a++) for (int b = 0; b < Mdim; b++) nz |= Nm.a[a][b];
        if (!nz) continue;
        Mat sq = matmul(&Nm, &Nm); int z = 1; for (int a = 0; a < Mdim; a++) for (int b = 0; b < Mdim; b++) if (sq.a[a][b]) z = 0;
        if (!z) continue;
        Mat Mm = Nm; for (int a = 0; a < Mdim; a++) Mm.a[a][a] ^= 1;
        return Mm;
    }
    fprintf(stderr, "random_involution2 failed\n"); exit(1);
}
static Mat rand_invol(const Vec *d, int fix_first, int allow_trivial) { return P == 2 ? random_involution2(fix_first ? d : NULL) : random_involution(d, fix_first, allow_trivial); }
static void random_symbol(int t) {
    if (inv_[t] != t) { Mt[t] = random_gl(); dt[t] = random_vec(dim_[t] >= 1); if (dim_[t] == 2) { do dt2[t] = random_vec(1); while (!indep2(&dt[t], &dt2[t])); } derive_pair(t); }
    else if (dim_[t] == 1) { dt[t] = random_vec(1); Mt[t] = rand_invol(&dt[t], 1, 0); }
    else { Mt[t] = rand_invol(NULL, 0, 1); Vec w = random_vec(0); Vec Mw = matvec(&Mt[t], &w); for (int i = 0; i < Mdim; i++) dt[t].v[i] = md(w.v[i] - Mw.v[i]); }
}
static void check_chart(void) {
    Mat I = identity();
    for (int t = 0; t < K; t++) {
        Mat inv; if (!mat_inv(&Mt[t], &inv)) { fprintf(stderr, "chart: singular M_%d\n", t); exit(1); }
        if (dim_[t] >= 1 && is_zero(&dt[t])) { fprintf(stderr, "chart: zero d_%d\n", t); exit(1); }
        if (dim_[t] == 2 && (inv_[t] == t || !indep2(&dt[t], &dt2[t]))) { fprintf(stderr, "chart: dim-2 symbol %d must be a pair with independent d, d2\n", t); exit(1); }
        if (inv_[t] == t) {
            Mat sq = matmul(&Mt[t], &Mt[t]); if (memcmp(&sq, &I, sizeof I)) { fprintf(stderr, "chart: M_%d not an involution\n", t); exit(1); }
            Vec w = matvec(&Mt[t], &dt[t]);
            for (int i = 0; i < Mdim; i++) if (md(w.v[i] - (dim_[t] == 1 ? dt[t].v[i] : -dt[t].v[i]))) { fprintf(stderr, "chart: M_%d d_%d condition violated\n", t, t); exit(1); }
        } else if (inv_[t] == t + 1) {
            Mat pr = matmul(&Mt[t], &Mt[t + 1]); if (memcmp(&pr, &I, sizeof I)) { fprintf(stderr, "chart: pair %d not inverse\n", t); exit(1); }
            Vec w = matvec(&Mt[t + 1], &dt[t]);
            for (int i = 0; i < Mdim; i++) if (md(w.v[i] - (dim_[t] >= 1 ? dt[t + 1].v[i] : -dt[t + 1].v[i]))) { fprintf(stderr, "chart: pair %d d condition violated\n", t); exit(1); }
            if (dim_[t] == 2) { Vec w2 = matvec(&Mt[t + 1], &dt2[t]); if (memcmp(&w2, &dt2[t + 1], sizeof w2)) { fprintf(stderr, "chart: pair %d d2 condition violated\n", t); exit(1); } }
        }
    }
}
// ---- lift evaluation ----
static int addT[MAXF][MAXF];
static int lin[MAXK][MAXF], tsz[MAXK], tset[MAXK][16];
static int vec_index(const Vec *v) { int x = 0, pw = 1; for (int i = 0; i < Mdim; i++) { x += v->v[i] * pw; pw *= P; } return x; }
static Vec index_vec(int x) { Vec v; for (int i = 0; i < Mdim; i++) { v.v[i] = x % P; x /= P; } return v; }
static int n, deg, W;              // lift order, degree, root words
static int *adj; static uint64_t *reach, *reach2;
static int nroots;
static void build_tables(void) {
    for (int t = 0; t < K; t++) {
        for (int x = 0; x < F; x++) { Vec v = index_vec(x); Vec w = matvec(&Mt[t], &v); lin[t][x] = vec_index(&w); }
        if (dim_[t] == 1) { tsz[t] = P; for (int l = 0; l < P; l++) { Vec w; for (int i = 0; i < Mdim; i++) w.v[i] = md(l * dt[t].v[i]); tset[t][l] = vec_index(&w); } }
        else if (dim_[t] == 2) { tsz[t] = P * P; for (int l = 0; l < P; l++) for (int l2 = 0; l2 < P; l2++) { Vec w; for (int i = 0; i < Mdim; i++) w.v[i] = md(l * dt[t].v[i] + l2 * dt2[t].v[i]); tset[t][l * P + l2] = vec_index(&w); } }
        else { tsz[t] = 1; tset[t][0] = vec_index(&dt[t]); }
    }
}
static void build_adj(void) {
    for (int v = 0; v < N; v++) for (int x = 0; x < F; x++) {
        int *row = adj + (long)(v * F + x) * deg; int c = 0;
        for (int t = 0; t < K; t++) { int w = nbr[v][t]; int lx = lin[t][x]; for (int j = 0; j < tsz[t]; j++) row[c++] = w * F + addT[lx][tset[t][j]]; }
    }
}
static long evaluate(int *simple_out) {
    build_tables(); build_adj();
    memset(reach, 0, (size_t)n * W * sizeof(uint64_t));
    for (int r = 0; r < R; r++) for (int x = 0; x < F; x++) { int u = roots[r] * F + x; int b = r * F + x; reach[(long)u * W + (b >> 6)] |= 1ULL << (b & 63); }
    for (int round = 0; round < D; round++) {
        for (int u = 0; u < n; u++) {
            uint64_t *out = reach2 + (long)u * W; const uint64_t *self = reach + (long)u * W;
            for (int w = 0; w < W; w++) out[w] = self[w];
            const int *row = adj + (long)u * deg;
            for (int i = 0; i < deg; i++) { const uint64_t *nb = reach + (long)row[i] * W; for (int w = 0; w < W; w++) out[w] |= nb[w]; }
        }
        uint64_t *tmp = reach; reach = reach2; reach2 = tmp;
    }
    long covered = 0;
    for (long i = 0; i < (long)n * W; i++) covered += __builtin_popcountll(reach[i]);
    if (simple_out) {   // loops / duplicate neighbours at the roots' fibres (G-invariant) -> degree defects
        int bad = 0;
        for (int r = 0; r < R && !bad; r++) for (int x = 0; x < F && !bad; x++) {
            int u = roots[r] * F + x; const int *row = adj + (long)u * deg;
            for (int i = 0; i < deg && !bad; i++) { if (row[i] == u) bad = 1; for (int j = i + 1; j < deg; j++) if (row[j] == row[i]) { bad = 1; break; } }
        }
        *simple_out = !bad;
    }
    return (long)nroots * n - covered;
}
static void save_chart(const char *fn, long obj, const char *tag) {
    FILE *f = fopen(fn, "w"); if (!f) return;
    fprintf(f, "# obj=%ld p=%d m=%d D=%d %s\n", obj, P, Mdim, D, tag ? tag : "");
    for (int t = 0; t < K; t++) {
        fprintf(f, "A");
        for (int i = 0; i < Mdim; i++) for (int j = 0; j < Mdim; j++) fprintf(f, " %d", Mt[t].a[i][j]);
        fprintf(f, dim_[t] >= 1 ? " d" : " c"); for (int i = 0; i < Mdim; i++) fprintf(f, " %d", dt[t].v[i]);
        if (dim_[t] == 2) { fprintf(f, " e"); for (int i = 0; i < Mdim; i++) fprintf(f, " %d", dt2[t].v[i]); }
        fprintf(f, "\n");
    }
    fclose(f);
}
static void load_chart(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror("init chart"); exit(2); }
    char line[8192]; int t = 0;
    while (fgets(line, sizeof line, f) && t < K) {
        if (line[0] != 'A') continue;
        char *p = line + 1;
        for (int i = 0; i < Mdim; i++) for (int j = 0; j < Mdim; j++) Mt[t].a[i][j] = md((int)strtol(p, &p, 10));
        while (*p && *p != 'd' && *p != 'c') p++; if (*p) p++;
        for (int i = 0; i < Mdim; i++) dt[t].v[i] = md((int)strtol(p, &p, 10));
        if (dim_[t] == 2) { while (*p && *p != 'e') p++; if (*p) p++; for (int i = 0; i < Mdim; i++) dt2[t].v[i] = md((int)strtol(p, &p, 10)); }
        t++;
    }
    fclose(f);
    if (t != K) { fprintf(stderr, "init chart has %d symbols, need %d\n", t, K); exit(2); }
}
static void write_edges(const char *fn) {
    // sorted edge list u < v
    long ne = (long)n * deg / 2; long *E = malloc(sizeof(long) * (size_t)n * deg); long c = 0;
    for (int u = 0; u < n; u++) { const int *row = adj + (long)u * deg; for (int i = 0; i < deg; i++) if (row[i] > u) E[c++] = (long)u * n + row[i]; }
    int cmp(const void *a, const void *b) { long x = *(const long *)a, y = *(const long *)b; return x < y ? -1 : x > y; }
    qsort(E, c, sizeof(long), cmp);
    FILE *f = fopen(fn, "w"); long dups = 0;
    for (long i = 0; i < c; i++) { if (i && E[i] == E[i - 1]) { dups++; continue; } fprintf(f, "%ld %ld\n", E[i] / n, E[i] % n); }
    fclose(f); free(E);
    if (dups || c != ne) fprintf(stderr, "edge list: %ld entries, %ld duplicates (not simple)\n", c, dups);
}
static void report(void) {
    // per (root base vertex, controller target) uncovered counts for the current reach (after evaluate)
    for (int r = 0; r < R; r++) {
        printf("root %d:", roots[r]);
        for (int w = 0; w < N; w++) {
            long unc = 0;
            for (int x = 0; x < F; x++) { int b = r * F + x; for (int y = 0; y < F; y++) { int u = w * F + y; if (!(reach[(long)u * W + (b >> 6)] >> (b & 63) & 1)) unc++; } }
            if (unc) printf(" %d:%ld", w, unc);
        }
        printf("\n");
    }
}
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: exlift ctrl.txt [options]\n"); return 2; }
    long steps = 200000; double tlimit = 1e9, T0 = 2.0, T1 = 0.05; uint64_t seed = 1; const char *init = NULL, *out = "exlift_best.txt", *efile = NULL; int do_report = 0, eval_only = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--pneg")) pneg = 1;
        else if (!strcmp(argv[i], "--p")) P = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--m")) Mdim = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--D")) D = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--steps")) steps = atol(argv[++i]);
        else if (!strcmp(argv[i], "--time")) tlimit = atof(argv[++i]);
        else if (!strcmp(argv[i], "--seed")) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--T0")) T0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T1")) T1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--init")) init = argv[++i];
        else if (!strcmp(argv[i], "--out")) out = argv[++i];
        else if (!strcmp(argv[i], "--edges")) efile = argv[++i];
        else if (!strcmp(argv[i], "--report")) do_report = 1;
        else if (!strcmp(argv[i], "--eval")) eval_only = 1;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    if (Mdim > MAXM) { fprintf(stderr, "m too large\n"); return 2; }
    F = 1; for (int i = 0; i < Mdim; i++) F *= P;
    if (F > MAXF) { fprintf(stderr, "F too large\n"); return 2; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror("ctrl"); return 2; }
    if (fscanf(f, "%d %d", &N, &K) != 2 || N > MAXN || K > MAXK) { fprintf(stderr, "bad controller header\n"); return 2; }
    for (int t = 0; t < K; t++) if (fscanf(f, "%d", &inv_[t]) != 1) return 2;
    for (int t = 0; t < K; t++) if (fscanf(f, "%d", &dim_[t]) != 1) return 2;
    for (int v = 0; v < N; v++) for (int t = 0; t < K; t++) if (fscanf(f, "%d", &nbr[v][t]) != 1) return 2;
    if (fscanf(f, "%d", &R) != 1) return 2;
    if (R == 0) { R = N; for (int r = 0; r < N; r++) roots[r] = r; } else for (int r = 0; r < R; r++) if (fscanf(f, "%d", &roots[r]) != 1) return 2;
    fclose(f);
    for (int t = 0; t < K; t++) { if (inv_[t] == t) sym_free[t] = 1; else if (inv_[t] == t + 1 && inv_[t + 1] == t) sym_free[t] = 1; else if (inv_[t] == t - 1 && inv_[t - 1] == t) sym_free[t] = 0; else { fprintf(stderr, "bad inv table\n"); return 2; } }
    for (int v = 0; v < N; v++) for (int t = 0; t < K; t++) if (nbr[nbr[v][t]][inv_[t]] != v) { fprintf(stderr, "controller: symbol %d at %d not invertible by %d\n", t, v, inv_[t]); return 2; }
    deg = 0; for (int t = 0; t < K; t++) deg += dim_[t] == 2 ? P * P : (dim_[t] == 1 ? P : 1);
    n = N * F; nroots = R * F; W = (nroots + 63) / 64;
    for (int a = 0; a < F; a++) for (int b = 0; b < F; b++) { Vec va = index_vec(a), vb = index_vec(b), s; for (int i = 0; i < Mdim; i++) s.v[i] = md(va.v[i] + vb.v[i]); addT[a][b] = vec_index(&s); }
    adj = malloc(sizeof(int) * (size_t)n * deg); reach = calloc((size_t)n * W, sizeof(uint64_t)); reach2 = calloc((size_t)n * W, sizeof(uint64_t));
    rng_state ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 20; i++) xr();
    fprintf(stderr, "controller N=%d K=%d deg=%d fibre F_%d^%d lift n=%d roots %d (%d words)\n", N, K, deg, P, Mdim, n, nroots, W);
    if (init) load_chart(init); else for (int t = 0; t < K; t++) if (sym_free[t]) random_symbol(t);
    if (pneg) for (int t = 0; t < K; t++) if (inv_[t] == t && dim_[t] == 0) { Mt[t] = identity(); for (int i = 0; i < Mdim; i++) Mt[t].a[i][i] = P - 1; if (!init) dt[t] = random_vec(0); }
    check_chart();
    int simple; long cur = evaluate(&simple);
    if (eval_only) {
        printf("objective %ld simple_at_roots %d n %d deg %d\n", cur, simple, n, deg);
        // distance distribution from root 0 by a plain BFS
        int *dist = malloc(n * sizeof(int)); int *q = malloc(n * sizeof(int)); for (int i = 0; i < n; i++) dist[i] = -1;
        int s = roots[0] * F, qh = 0, qt = 0; dist[s] = 0; q[qt++] = s;
        while (qh < qt) { int x = q[qh++]; const int *row = adj + (long)x * deg; for (int i = 0; i < deg; i++) if (dist[row[i]] < 0) { dist[row[i]] = dist[x] + 1; q[qt++] = row[i]; } }
        long hist[64] = {0}; int maxd = 0; for (int i = 0; i < n; i++) { if (dist[i] >= 0) { hist[dist[i]]++; if (dist[i] > maxd) maxd = dist[i]; } else hist[63]++; }
        printf("distance distribution from (root %d, 0):", roots[0]); for (int i = 0; i <= maxd; i++) printf(" %ld", hist[i]); if (hist[63]) printf(" unreachable %ld", hist[63]); printf("\n");
        if (do_report) report();
        if (efile) write_edges(efile);
        return 0;
    }
    long best = cur; Mat bestM[MAXK]; Vec bestd[MAXK], bestd2[MAXK]; memcpy(bestM, Mt, sizeof Mt); memcpy(bestd, dt, sizeof dt); memcpy(bestd2, dt2, sizeof dt2);
    fprintf(stderr, "init objective %ld (simple %d)\n", cur, simple);
    save_chart(out, best, "init");
    time_t t0 = time(NULL); long accepted = 0; Mat saveM[MAXK]; Vec saved[MAXK], saved2[MAXK];
    int nfree = 0, freelist[MAXK]; for (int t = 0; t < K; t++) if (sym_free[t]) freelist[nfree++] = t;
    long st;
    for (st = 0; st < steps; st++) {
        if ((st & 255) == 0 && difftime(time(NULL), t0) > tlimit) break;
        double Temp = T0 * pow(T1 / T0, (double)st / (double)(steps > 1 ? steps - 1 : 1));
        memcpy(saveM, Mt, sizeof Mt); memcpy(saved, dt, sizeof dt); memcpy(saved2, dt2, sizeof dt2);
        int t = freelist[rnd(nfree)]; int mv = rnd(10);
        if (inv_[t] != t) {   // pair
            if (mv < 3) { Mat inv; int j = rnd(Mdim); for (int tries = 0; tries < 50; tries++) { for (int i = 0; i < Mdim; i++) Mt[t].a[i][j] = rnd(P); if (mat_inv(&Mt[t], &inv)) break; if (tries == 49) Mt[t] = saveM[t]; } }
            else if (mv < 6) { int i0 = rnd(Mdim), j0 = rnd(Mdim), l = 1 + rnd(P - 1); if (i0 == j0) j0 = (j0 + 1) % Mdim; for (int i = 0; i < Mdim; i++) Mt[t].a[i][j0] = md(Mt[t].a[i][j0] + l * Mt[t].a[i][i0]); }
            else if (mv < 8) {
                if (dim_[t] == 2 && rnd(2)) { int j = rnd(Mdim); dt2[t].v[j] = md(dt2[t].v[j] + 1 + rnd(P - 1)); if (!indep2(&dt[t], &dt2[t])) dt2[t] = saved2[t]; }
                else { int j = rnd(Mdim); dt[t].v[j] = md(dt[t].v[j] + 1 + rnd(P - 1)); if ((dim_[t] == 1 && is_zero(&dt[t])) || (dim_[t] == 2 && !indep2(&dt[t], &dt2[t]))) dt[t] = saved[t]; }
            }
            else { dt[t] = random_vec(dim_[t] >= 1); if (dim_[t] == 2) { do dt2[t] = random_vec(1); while (!indep2(&dt[t], &dt2[t])); } }
            derive_pair(t);
        } else if (dim_[t] == 1) {   // self-inverse line symbol
            if (mv < 3) { dt[t] = random_vec(1); Mt[t] = rand_invol(&dt[t], 1, 0); }
            else if (mv < 5) { Mt[t] = rand_invol(&dt[t], 1, 0); }
            else if (mv < 8) conjugate_fixing(&Mt[t], &dt[t]);
            else { // move d within the +1 eigenspace: d' = d + (w + M w) for random w (M(w + Mw) = Mw + w)
                Vec w = random_vec(0), Mw = matvec(&Mt[t], &w); Vec nd; for (int i = 0; i < Mdim; i++) nd.v[i] = md(dt[t].v[i] + w.v[i] + Mw.v[i]);
                if (!is_zero(&nd)) dt[t] = nd;
            }
        } else if (pneg) {   // plain symbol frozen at M = -I: only c moves (any vector is in the (-1)-eigenspace)
            if (mv < 5) { int j = rnd(Mdim); dt[t].v[j] = md(dt[t].v[j] + 1 + rnd(P - 1)); } else dt[t] = random_vec(0);
        } else {   // self-inverse plain symbol: M involution, c in the (-1)-eigenspace
            if (mv < 3) { Mt[t] = rand_invol(NULL, 0, 1); Vec w = random_vec(0), Mw = matvec(&Mt[t], &w); for (int i = 0; i < Mdim; i++) dt[t].v[i] = md(w.v[i] - Mw.v[i]); }
            else if (mv < 6) { conjugate_fixing(&Mt[t], NULL); Vec w = random_vec(0), Mw = matvec(&Mt[t], &w); Vec c = dt[t]; Vec Mc = matvec(&Mt[t], &c); int ok = 1; for (int i = 0; i < Mdim; i++) if (md(Mc.v[i] + c.v[i])) ok = 0; if (!ok) for (int i = 0; i < Mdim; i++) dt[t].v[i] = md(w.v[i] - Mw.v[i]); }
            else { Vec w = random_vec(0), Mw = matvec(&Mt[t], &w); for (int i = 0; i < Mdim; i++) dt[t].v[i] = md(w.v[i] - Mw.v[i]); }
        }
        long v = evaluate(NULL);
        if (v <= cur || rndu() < exp((cur - v) / Temp)) {
            cur = v; accepted++;
            if (v < best) {
                best = v; memcpy(bestM, Mt, sizeof Mt); memcpy(bestd, dt, sizeof dt); memcpy(bestd2, dt2, sizeof dt2);
                fprintf(stderr, "step %ld T=%.3f best %ld (%.0fs)\n", st, Temp, best, difftime(time(NULL), t0));
                save_chart(out, best, "best");
                if (best == 0) break;
            }
        } else { memcpy(Mt, saveM, sizeof Mt); memcpy(dt, saved, sizeof dt); memcpy(dt2, saved2, sizeof dt2); }
        if (st % 20000 == 0) fprintf(stderr, "step %ld T=%.3f cur=%ld best=%ld acc=%ld (%.0fs)\n", st, Temp, cur, best, accepted, difftime(time(NULL), t0));
    }
    memcpy(Mt, bestM, sizeof Mt); memcpy(dt, bestd, sizeof dt); memcpy(dt2, bestd2, sizeof dt2); check_chart();
    long fin = evaluate(&simple);
    save_chart(out, fin, "final");
    printf("RESULT obj=%ld simple=%d n=%d deg=%d steps=%ld accepted=%ld seconds=%.0f chart=%s\n", fin, simple, n, deg, st, accepted, difftime(time(NULL), t0), out);
    if (do_report) report();
    if (fin == 0 && simple && efile) { write_edges(efile); printf("LIFT_OK n=%d deg=%d D=%d edges=%s\n", n, deg, D, efile); }
    return 0;
}
