// mixlift: exact lift evaluator / SA for route-chart lifts with MIXED alphabets over F_3^4
//   line symbols 0 = a, 1 = a^-1, 2 = b, 3 = b^-1: (v, x) ~ (t(v), A_t x + lambda d_t), lambda in F_3  (3 neighbours each)
//   plain symbols 4 .. 3 + P:                     (v, x) ~ (q(v), M_q x + c_q)                     (1 neighbour each)
// Controller: regular lift over a group table of a k-vertex base (alpha0, beta0; voltages va, vb; mixctrl conventions,
// vertex (c, g) = c * nG + g) with the plain symbols given by CENTRAL group elements (mixctrl's --plain items):
//   inv z     a self-inverse symbol s: (c, g) -> (c, g z); chart M_s = -I (scalar), c_s free (M^2 = I, (M + I) c = 0)
//   pair z    p, p^-1: (c, g) -> (c, g z^{+-1}); chart of p: M = +I or -I (--Mp s1[,s2] per pair, default -1), c free;
//             chart of p^-1: (M^-1, -M^-1 c)
//   e.g. "inv 42" (degree 13), "pair 1" or the older spelling "z3 1" (14), "inv 13 13" (14), "pair 5 inv 13" (15), "pair 3 pair 9" (16)
// Objective (exact): number of (root, vertex) pairs at distance > 4, roots (c, e, x) for the k base vertices and
// all 81 x (the deck group G acts on the lift), bit-parallel BFS.  Total ordered far pairs = objective x |G|.
// Usage: mixlift group.tbl k a0 b0 "va=.. vb=.." --plain ITEMS steps seed [--T0 T] [--T1 T] [--time S] [--out F]
//        [--init chart] [--edges F] [--tag T] [--eval] [--fibres] [--polish S] [--Mp s1[,s2]]
// Chart file: 4 lines "A <16 entries row-major> d <4>" (a, a^-1, b, b^-1), P lines "P <16 entries> c <4>" (plain symbols
//   in order), "controller ...".  Six-symbol chart files load unchanged.
// --perdart (k >= 2): one chart per base DART instead of per symbol: the line dart (c, a) carries (A, d), its reverse
//   (alpha0(c), a^-1) the derived (A^-1, A^-1 d), same for b; the plain loops at c carry their own constants.  k x (2 line
//   charts + P plain constants) free parameters; the hard / soft screen is unchanged (the plain letters are loops, so the four
//   insertions of a plain letter into a 3-word apply the same dart charts).  Chart file lines "D c t <16> d <4>" / "Q c t <16> c <4>".
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define P3 3
#define M 4
#define F 81
#define KL 4
#define MAXP 4
#define MAXSYM (KL + MAXP)
#define MAXN 1024
#define MAXK 8
#define RW ((MAXK * 81 + 63) / 64)

typedef struct { int a[M][M]; } Mat;
typedef struct { int v[M]; } Vec;

static Mat A[MAXSYM]; static Vec dv[MAXSYM];   // line: (A, d); plain: (M, c)   (per symbol)
static Mat Ad[MAXK][MAXSYM]; static Vec dd[MAXK][MAXSYM];   // per base dart (c, t): what build_adj uses
static int perdart = 0;
static int k = 1, nG, N, K = 4, P = 0, deg = 12;
static int *T, *ginv; static int a0[MAXK], b0[MAXK], va[MAXK], vb[MAXK];
// plain groups: ngrp items; grp_kind[g] 0 = inv (one symbol), 1 = pair (two symbols); grp_sym[g] = first symbol index
static int ngrp = 0, grp_kind[MAXP], grp_sym[MAXP], grp_sign[MAXP], pz[MAXP];
static char plainspec[64] = "";
static int fibres = 0, polish = 0;
static int nbr[MAXN][MAXSYM];
static uint64_t rng_state = 88172645463325252ULL;
static uint64_t xr(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }
static int rnd(int n) { return (int)(xr() % (uint64_t)n); }
static double rndu(void) { return (xr() >> 11) * (1.0 / 9007199254740992.0); }
static int md(int x) { x %= P3; return x < 0 ? x + P3 : x; }
static Vec matvec(const Mat *m, const Vec *v) { Vec r; for (int i = 0; i < M; i++) { int s = 0; for (int j = 0; j < M; j++) s += m->a[i][j] * v->v[j]; r.v[i] = md(s); } return r; }
static Mat matmul(const Mat *x, const Mat *y) { Mat r; for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { int s = 0; for (int kk = 0; kk < M; kk++) s += x->a[i][kk] * y->a[kk][j]; r.a[i][j] = md(s); } return r; }
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
static int vec_index(const Vec *v) { return v->v[0] + 3 * v->v[1] + 9 * v->v[2] + 27 * v->v[3]; }
static Vec index_vec(int x) { Vec v; for (int i = 0; i < M; i++) { v.v[i] = x % 3; x /= 3; } return v; }
static int is_zero(const Vec *v) { return !(v->v[0] | v->v[1] | v->v[2] | v->v[3]); }
static Vec random_nonzero(void) { Vec v; do { for (int i = 0; i < M; i++) v.v[i] = rnd(P3); } while (is_zero(&v)); return v; }
static Vec random_vec(void) { Vec v; for (int i = 0; i < M; i++) v.v[i] = rnd(P3); return v; }
static Mat random_gl(void) { Mat m, inv; do { for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) m.a[i][j] = rnd(P3); } while (!mat_inv(&m, &inv)); return m; }
static Vec neg(const Vec *v) { Vec r; for (int i = 0; i < M; i++) r.v[i] = md(-v->v[i]); return r; }
static Mat scalar(int s) { Mat m; memset(&m, 0, sizeof m); for (int i = 0; i < M; i++) m.a[i][i] = md(s); return m; }

// ---- chart maintenance ----
static void derive_pair(int t) {   // line pair, t even
    Mat inv; if (!mat_inv(&A[t], &inv)) { fprintf(stderr, "singular pair\n"); exit(1); }
    A[t + 1] = inv; dv[t + 1] = matvec(&inv, &dv[t]);
}
static void derive_plain(int g) {  // the second symbol of a pair from the first
    if (grp_kind[g] != 1) return;
    int t = grp_sym[g];
    Mat inv; if (!mat_inv(&A[t], &inv)) { fprintf(stderr, "singular M_p\n"); exit(1); }
    A[t + 1] = inv; Vec w = matvec(&inv, &dv[t]); dv[t + 1] = neg(&w);
}
// two involution symbols on the same element must carry different constants (simple lift)
static int plain_clash(void) {
    for (int i = 0; i < P; i++) for (int j = i + 1; j < P; j++) if (pz[i] == pz[j] && !memcmp(&dv[KL + i], &dv[KL + j], sizeof(Vec)) && !memcmp(&A[KL + i], &A[KL + j], sizeof(Mat))) return 1;
    return 0;
}
static void random_plain(void) {
    for (int g = 0; g < ngrp; g++) { int t = grp_sym[g]; A[t] = scalar(grp_kind[g] == 0 ? -1 : grp_sign[g]); dv[t] = random_vec(); derive_plain(g); }
    for (int tries = 0; plain_clash() && tries < 100; tries++) { int g = rnd(ngrp); dv[grp_sym[g]] = random_vec(); derive_plain(g); }
}
static void check_chart(void) {
    Mat I = scalar(1);
    for (int t = 0; t < KL; t += 2) {
        Mat pr = matmul(&A[t], &A[t + 1]); if (memcmp(&pr, &I, sizeof I)) { fprintf(stderr, "pair check failed\n"); exit(1); }
        Vec w = matvec(&A[t + 1], &dv[t]); if (memcmp(&w, &dv[t + 1], sizeof w)) { fprintf(stderr, "pair d check failed\n"); exit(1); }
        if (is_zero(&dv[t])) { fprintf(stderr, "zero d\n"); exit(1); }
    }
    for (int g = 0; g < ngrp; g++) {
        int t = grp_sym[g];
        if (grp_kind[g] == 0) {
            Mat sq = matmul(&A[t], &A[t]); if (memcmp(&sq, &I, sizeof I)) { fprintf(stderr, "plain involution check failed\n"); exit(1); }
            Vec w = matvec(&A[t], &dv[t]); Vec nw = neg(&w); if (memcmp(&nw, &dv[t], sizeof w)) { fprintf(stderr, "plain (M+I)c = 0 check failed\n"); exit(1); }
        } else {
            Mat pr = matmul(&A[t], &A[t + 1]); if (memcmp(&pr, &I, sizeof I)) { fprintf(stderr, "plain pair check failed\n"); exit(1); }
            Vec w = matvec(&A[t + 1], &dv[t]); Vec s; for (int i = 0; i < M; i++) s.v[i] = md(w.v[i] + dv[t + 1].v[i]); if (!is_zero(&s)) { fprintf(stderr, "plain pair c check failed\n"); exit(1); }
        }
    }
    if (plain_clash()) { fprintf(stderr, "double plain edge with equal charts (lift not simple)\n"); exit(1); }
}
static void random_chart(void) {
    for (int t = 0; t < KL; t += 2) { A[t] = random_gl(); dv[t] = random_nonzero(); derive_pair(t); }
    random_plain();
}
// ---- per-dart charts ----
static void sync_charts(void) { if (perdart) return; for (int c = 0; c < k; c++) for (int t = 0; t < K; t++) { Ad[c][t] = A[t]; dd[c][t] = dv[t]; } }
static int partner(int c, int t) { return t == 0 ? a0[c] : b0[c]; }   // base vertex at the head of the dart (c, t), t in {0, 2}
static void derive_pair_d(int c, int t) {   // dart (c, t) -> reverse dart (partner, t + 1)
    Mat inv; if (!mat_inv(&Ad[c][t], &inv)) { fprintf(stderr, "singular pair\n"); exit(1); }
    int h = partner(c, t); Ad[h][t + 1] = inv; dd[h][t + 1] = matvec(&inv, &dd[c][t]);
}
static void derive_plain_d(int c, int g) {
    if (grp_kind[g] != 1) return;
    int t = grp_sym[g];
    Mat inv; if (!mat_inv(&Ad[c][t], &inv)) { fprintf(stderr, "singular M_p\n"); exit(1); }
    Ad[c][t + 1] = inv; Vec w = matvec(&inv, &dd[c][t]); dd[c][t + 1] = neg(&w);
}
static int plain_clash_d(void) {
    for (int c = 0; c < k; c++) for (int i = 0; i < P; i++) for (int j = i + 1; j < P; j++) if (pz[i] == pz[j] && !memcmp(&dd[c][KL + i], &dd[c][KL + j], sizeof(Vec)) && !memcmp(&Ad[c][KL + i], &Ad[c][KL + j], sizeof(Mat))) return 1;
    return 0;
}
static void random_chart_d(void) {
    for (int c = 0; c < k; c++) for (int t = 0; t < KL; t += 2) { Ad[c][t] = random_gl(); dd[c][t] = random_nonzero(); derive_pair_d(c, t); }
    for (int c = 0; c < k; c++) for (int g = 0; g < ngrp; g++) { int t = grp_sym[g]; Ad[c][t] = scalar(grp_kind[g] == 0 ? -1 : grp_sign[g]); dd[c][t] = random_vec(); derive_plain_d(c, g); }
    for (int tries = 0; plain_clash_d() && tries < 100; tries++) { int c = rnd(k), g = rnd(ngrp); dd[c][grp_sym[g]] = random_vec(); derive_plain_d(c, g); }
}
static void check_chart_d(void) {
    Mat I = scalar(1);
    for (int c = 0; c < k; c++) {
        for (int t = 0; t < KL; t += 2) {
            int h = partner(c, t);
            Mat pr = matmul(&Ad[c][t], &Ad[h][t + 1]); if (memcmp(&pr, &I, sizeof I)) { fprintf(stderr, "dart pair check failed\n"); exit(1); }
            Vec w = matvec(&Ad[h][t + 1], &dd[c][t]); if (memcmp(&w, &dd[h][t + 1], sizeof w)) { fprintf(stderr, "dart pair d check failed\n"); exit(1); }
            if (is_zero(&dd[c][t])) { fprintf(stderr, "zero d\n"); exit(1); }
        }
        for (int g = 0; g < ngrp; g++) {
            int t = grp_sym[g];
            if (grp_kind[g] == 0) { Mat sq = matmul(&Ad[c][t], &Ad[c][t]); Vec w = matvec(&Ad[c][t], &dd[c][t]); Vec nw = neg(&w);
                if (memcmp(&sq, &I, sizeof I) || memcmp(&nw, &dd[c][t], sizeof w)) { fprintf(stderr, "dart plain involution check failed\n"); exit(1); } }
            else { Mat pr = matmul(&Ad[c][t], &Ad[c][t + 1]); Vec w = matvec(&Ad[c][t + 1], &dd[c][t]); Vec s; for (int i = 0; i < M; i++) s.v[i] = md(w.v[i] + dd[c][t + 1].v[i]);
                if (memcmp(&pr, &I, sizeof I) || !is_zero(&s)) { fprintf(stderr, "dart plain pair check failed\n"); exit(1); } }
        }
    }
    if (plain_clash_d()) { fprintf(stderr, "double plain edge with equal charts (lift not simple)\n"); exit(1); }
}

// ---- lift + BFS ----
static int *adj;
static uint64_t (*reach)[RW], (*front)[RW], (*nxt)[RW];
static long n;
static void build_adj(void) {
    for (int v = 0; v < N; v++) for (int x = 0; x < F; x++) {
        int *row = adj + ((long)v * F + x) * deg; int c = 0;
        Vec xv = index_vec(x);
        int bc = v / nG;
        for (int t = 0; t < K; t++) {
            int w = nbr[v][t];
            Vec y = matvec(&Ad[bc][t], &xv);
            if (t < KL) {
                for (int l = 0; l < P3; l++) { Vec z; for (int i = 0; i < M; i++) z.v[i] = md(y.v[i] + l * dd[bc][t].v[i]); row[c++] = w * F + vec_index(&z); }
            } else {
                Vec z; for (int i = 0; i < M; i++) z.v[i] = md(y.v[i] + dd[bc][t].v[i]); row[c++] = w * F + vec_index(&z);
            }
        }
    }
}
static long check_symmetry(long *loops, long *dups) {
    long asym = 0; *loops = 0; *dups = 0;
    for (long u = 0; u < n; u++) {
        int *row = adj + u * deg;
        for (int i = 0; i < deg; i++) {
            long w = row[i]; if (w == u) (*loops)++;
            for (int j = i + 1; j < deg; j++) if (row[j] == w) (*dups)++;
            int *rw = adj + w * deg, ok = 0; for (int j = 0; j < deg; j++) if (rw[j] == u) { ok = 1; break; }
            asym += !ok;
        }
    }
    return asym;
}
static int nrw;
static long evaluate(void) {
    memset(reach, 0, n * sizeof(*reach)); memset(front, 0, n * sizeof(*front));
    for (int c = 0; c < k; c++) for (int x = 0; x < F; x++) {
        long r = c * F + x; long v = ((long)c * nG) * F + x;
        reach[v][r >> 6] |= 1ULL << (r & 63); front[v][r >> 6] |= 1ULL << (r & 63);
    }
    for (int layer = 0; layer < 4; layer++) {
        memset(nxt, 0, n * sizeof(*nxt));
        for (long v = 0; v < n; v++) {
            uint64_t any = 0; for (int w = 0; w < nrw; w++) any |= front[v][w];
            if (!any) continue;
            int *row = adj + v * deg;
            for (int i = 0; i < deg; i++) { uint64_t *nx = nxt[row[i]]; for (int w = 0; w < nrw; w++) nx[w] |= front[v][w]; }
        }
        for (long v = 0; v < n; v++) for (int w = 0; w < nrw; w++) { uint64_t nw = nxt[v][w] & ~reach[v][w]; front[v][w] = nw; reach[v][w] |= nw; }
    }
    long unc = 0; long nr = (long)k * F;
    for (long v = 0; v < n; v++) for (int w = 0; w < nrw; w++) {
        uint64_t mask = (nr - 64L * w >= 64) ? ~0ULL : ((1ULL << (nr - 64L * w)) - 1);
        unc += __builtin_popcountll(~reach[v][w] & mask);
    }
    return unc;
}
static int build_controller(void) {
    static int alpha[MAXN], beta[MAXN], alphai[MAXN], betai[MAXN];
    for (int c = 0; c < k; c++) for (int g = 0; g < nG; g++) {
        int v = c * nG + g;
        alpha[v] = a0[c] * nG + T[g * nG + va[c]];
        beta[v] = b0[c] * nG + T[g * nG + vb[c]];
    }
    for (int v = 0; v < N; v++) { alphai[alpha[v]] = v; betai[beta[v]] = v; }
    for (int c = 0; c < k; c++) for (int g = 0; g < nG; g++) {
        int v = c * nG + g;
        nbr[v][0] = alpha[v]; nbr[v][1] = alphai[v]; nbr[v][2] = beta[v]; nbr[v][3] = betai[v];
        for (int j = 0; j < P; j++) nbr[v][KL + j] = c * nG + T[g * nG + pz[j]];
    }
    for (int v = 0; v < N; v++) {
        for (int i = 0; i < K; i++) { if (nbr[v][i] == v) return 0; for (int j = i + 1; j < K; j++) if (nbr[v][i] == nbr[v][j] && (i < KL || j < KL || pz[i - KL] != pz[j - KL])) return 0; }
    }
    return 1;
}
static void save_chart(const char *fn, long best, const char *ctrl_desc) {
    FILE *f = fopen(fn, "w"); if (!f) return;
    fprintf(f, "# mixlift unc=%ld deg=%d N=%d n=%ld plain=%s%s\n", best, deg, N, n, plainspec, perdart ? " perdart" : "");
    if (perdart) {
        for (int c = 0; c < k; c++) for (int t = 0; t < K; t++) {
            fprintf(f, "%s %d %d", t < KL ? "D" : "Q", c, t);
            for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) fprintf(f, " %d", Ad[c][t].a[i][j]);
            fprintf(f, t < KL ? " d" : " c"); for (int i = 0; i < M; i++) fprintf(f, " %d", dd[c][t].v[i]);
            fprintf(f, "\n");
        }
    } else
    for (int t = 0; t < K; t++) {
        fprintf(f, t < KL ? "A" : "P");
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) fprintf(f, " %d", A[t].a[i][j]);
        fprintf(f, t < KL ? " d" : " c"); for (int i = 0; i < M; i++) fprintf(f, " %d", dv[t].v[i]);
        fprintf(f, "\n");
    }
    fprintf(f, "controller %s\n", ctrl_desc);
    fclose(f);
}
static void load_chart(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror("init"); exit(2); }
    char line[4096]; int t = 0;
    int nd = 0;
    while (fgets(line, sizeof line, f) && t < K) {
        if (perdart) {
            if (line[0] != 'D' && line[0] != 'Q') continue;
            char *p = line + 1; int c = (int)strtol(p, &p, 10), tt = (int)strtol(p, &p, 10);
            if (c < 0 || c >= k || tt < 0 || tt >= K) { fprintf(stderr, "bad dart line\n"); exit(2); }
            for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) Ad[c][tt].a[i][j] = (int)strtol(p, &p, 10);
            while (*p && *p != 'd' && *p != 'c') p++;
            p++;
            for (int i = 0; i < M; i++) dd[c][tt].v[i] = (int)strtol(p, &p, 10);
            nd++; continue;
        }
        if (line[0] != 'A' && line[0] != 'P') continue;
        char *p = line + 1;
        for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) A[t].a[i][j] = (int)strtol(p, &p, 10);
        while (*p && *p != 'd' && *p != 'c') p++;
        p++;
        for (int i = 0; i < M; i++) dv[t].v[i] = (int)strtol(p, &p, 10);
        t++;
    }
    fclose(f);
    if (perdart) { if (nd != k * K) { fprintf(stderr, "init chart incomplete (%d of %d dart lines)\n", nd, k * K); exit(2); } return; }
    if (t != K) { fprintf(stderr, "init chart incomplete (%d of %d symbol lines)\n", t, K); exit(2); }
}
static int parse_perm(const char *s, int *out) {
    int cnt = 0; const char *p = s;
    while (*p && *p != ' ' && cnt < MAXK) { out[cnt++] = atoi(p); while (*p && *p != ',' && *p != ' ') p++; if (*p == ',') p++; }
    return cnt;
}
static void parse_volts(const char *s, int *va, int *vb) {
    const char *p;
    if ((p = strstr(s, "va="))) parse_perm(p + 3, va);
    if ((p = strstr(s, "vb="))) parse_perm(p + 3, vb);
}
static void ctrl_string(char *buf, size_t sz, char **argv) {
    int p = snprintf(buf, sz, "base=%s/%s va=", argv[3], argv[4]);
    for (int c = 0; c < k; c++) p += snprintf(buf + p, sz - p, "%d%s", va[c], c + 1 < k ? "," : "");
    p += snprintf(buf + p, sz - p, " vb=");
    for (int c = 0; c < k; c++) p += snprintf(buf + p, sz - p, "%d%s", vb[c], c + 1 < k ? "," : "");
    p += snprintf(buf + p, sz - p, " plain=%s%s Mp=", plainspec, perdart ? " perdart" : "");
    for (int g = 0; g < ngrp; g++) p += snprintf(buf + p, sz - p, "%d%s", grp_kind[g] ? grp_sign[g] : -1, g + 1 < ngrp ? "," : "");
    snprintf(buf + p, sz - p, " group=%s", argv[1]);
}
static int is_int(const char *s) { if (*s == '-') s++; if (!*s) return 0; while (*s) { if (*s < '0' || *s > '9') return 0; s++; } return 1; }
static int gorder(int z) { int o = 1, x = z; while (x != 0) { x = T[x * nG + z]; o++; } return o; }
static int add_plain(const char *kind, int z) {
    if (z <= 0 || z >= nG) { fprintf(stderr, "plain element %d out of range\n", z); return 0; }
    for (int g = 0; g < nG; g++) if (T[g * nG + z] != T[z * nG + g]) { fprintf(stderr, "plain element %d not central\n", z); return 0; }
    int o = gorder(z), g = ngrp;
    if (!strcmp(kind, "inv")) {
        if (o != 2) { fprintf(stderr, "plain element %d not an involution\n", z); return 0; }
        if (P + 1 > MAXP) { fprintf(stderr, "too many plain symbols\n"); return 0; }
        grp_kind[g] = 0; grp_sym[g] = KL + P; pz[P++] = z;
        snprintf(plainspec + strlen(plainspec), sizeof plainspec - strlen(plainspec), "%sinv:%d", g ? "," : "", z);
    } else {
        if (o < 3) { fprintf(stderr, "plain element %d is an involution: use inv\n", z); return 0; }
        if (P + 2 > MAXP) { fprintf(stderr, "too many plain symbols\n"); return 0; }
        grp_kind[g] = 1; grp_sym[g] = KL + P; pz[P++] = z; pz[P++] = ginv[z];
        snprintf(plainspec + strlen(plainspec), sizeof plainspec - strlen(plainspec), "%spair:%d", g ? "," : "", z);
    }
    grp_sign[g] = -1; ngrp++;
    return 1;
}
static void write_edges(const char *efile) {
    FILE *f = fopen(efile, "w");
    for (long u = 0; u < n; u++) { int *row = adj + u * deg; for (int i = 0; i < deg; i++) if (row[i] > u) fprintf(f, "%ld %d\n", u, row[i]); }
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc < 9) { fprintf(stderr, "usage: mixlift group.tbl k a0 b0 \"va=.. vb=..\" --plain ITEMS steps seed [opts]\n"); return 2; }
    k = atoi(argv[2]);
    if (parse_perm(argv[3], a0) != k || parse_perm(argv[4], b0) != k) { fprintf(stderr, "perm length\n"); return 2; }
    parse_volts(argv[5], va, vb);
    FILE *gf = fopen(argv[1], "r"); if (!gf) { perror("group"); return 2; }
    if (fscanf(gf, "%d", &nG) != 1) return 2;
    T = malloc(nG * nG * sizeof(int)); ginv = malloc(nG * sizeof(int));
    for (int i = 0; i < nG * nG; i++) if (fscanf(gf, "%d", &T[i]) != 1) return 2;
    fclose(gf);
    for (int g = 0; g < nG; g++) for (int h = 0; h < nG; h++) if (T[g * nG + h] == 0) ginv[g] = h;
    int ai = 6;
    if (strcmp(argv[ai], "--plain")) { fprintf(stderr, "--plain expected\n"); return 2; }
    ai++;
    while (ai < argc && !is_int(argv[ai])) {
        const char *kind = argv[ai++];
        if (strcmp(kind, "inv") && strcmp(kind, "pair") && strcmp(kind, "z3")) { fprintf(stderr, "plain item kind %s\n", kind); return 2; }
        if (!strcmp(kind, "z3")) kind = "pair";
        if (ai >= argc || !is_int(argv[ai])) { fprintf(stderr, "plain item needs an element\n"); return 2; }
        if (!add_plain(kind, atoi(argv[ai++]))) return 2;
        // "inv z1 z2" (two involutions): a second integer right after an inv item is another involution
        if (!strcmp(kind, "inv") && ai + 2 < argc && is_int(argv[ai]) && is_int(argv[ai + 1]) && is_int(argv[ai + 2])) { if (!add_plain("inv", atoi(argv[ai++]))) return 2; }
    }
    if (!P) { fprintf(stderr, "--plain needs items\n"); return 2; }
    K = KL + P; deg = 3 * KL + P;
    if (argc < ai + 2) { fprintf(stderr, "steps seed missing\n"); return 2; }
    long steps = atol(argv[ai]); uint64_t seed = atoll(argv[ai + 1]);
    double T0 = 2.0, T1 = 0.1, tlimit = 1e9; const char *out = "mixlift_out.txt", *init = NULL, *efile = NULL, *tag = "x"; int eval_only = 0;
    for (int i = ai + 2; i < argc; i++) {
        if (!strcmp(argv[i], "--T0")) T0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T1")) T1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--time")) tlimit = atof(argv[++i]);
        else if (!strcmp(argv[i], "--out")) out = argv[++i];
        else if (!strcmp(argv[i], "--init")) init = argv[++i];
        else if (!strcmp(argv[i], "--edges")) efile = argv[++i];
        else if (!strcmp(argv[i], "--tag")) tag = argv[++i];
        else if (!strcmp(argv[i], "--eval")) eval_only = 1;
        else if (!strcmp(argv[i], "--Mp")) { const char *p = argv[++i]; int g = 0; while (*p && g < ngrp) { int s = (int)strtol(p, (char **)&p, 10); if (grp_kind[g] == 1) grp_sign[g] = s < 0 ? -1 : 1; g++; if (*p == ',') p++; } }
        else if (!strcmp(argv[i], "--fibres")) fibres = 1;
        else if (!strcmp(argv[i], "--polish")) polish = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--perdart")) perdart = 1;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    rng_state ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    N = k * nG; if (N > MAXN) { fprintf(stderr, "N too big\n"); return 2; }
    nrw = (k * F + 63) / 64;
    if (!build_controller()) { fprintf(stderr, "controller not simple\n"); return 2; }
    n = (long)N * F;
    adj = malloc(n * deg * sizeof(int));
    reach = malloc(n * sizeof(*reach)); front = malloc(n * sizeof(*front)); nxt = malloc(n * sizeof(*nxt));
    if (perdart && k < 2) { fprintf(stderr, "--perdart needs k >= 2\n"); return 2; }
    if (perdart && polish) { fprintf(stderr, "--polish is not implemented with --perdart\n"); return 2; }
    if (init) load_chart(init); else if (perdart) random_chart_d(); else random_chart();
    if (perdart) check_chart_d(); else { check_chart(); sync_charts(); }
    build_adj();
    long loops, dups; long asym = check_symmetry(&loops, &dups);
    if (asym || loops || dups) { fprintf(stderr, "lift not a simple symmetric graph: asym %ld loops %ld dups %ld\n", asym, loops, dups); return 3; }
    char ctrl_desc[1024]; ctrl_string(ctrl_desc, sizeof ctrl_desc, argv);
    long cur = evaluate(), best = cur; double curv = cur;
    fprintf(stderr, "plain %s deg %d N %d n %ld roots %ld: init uncovered %ld (of %ld)\n", plainspec, deg, N, n, (long)k * F, cur, (long)k * F * n);
    if (eval_only) {
        save_chart(out, cur, ctrl_desc);
        if (fibres) {
            // per (root base vertex c, target controller vertex v): number of (x, y) with (c, e, x) !~ (v, y) within 4
            long tot = 0, nf = 0;
            for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) {
                long miss = 0;
                for (int y = 0; y < F; y++) { long u = (long)v * F + y; for (int x = 0; x < F; x++) { long r = c * F + x; if (!(reach[u][r >> 6] >> (r & 63) & 1)) miss++; } }
                if (miss) { printf("FIBRE c=%d v=%d missing=%ld of %d\n", c, v, miss, F * F); tot += miss; nf++; }
            }
            printf("FIBRES incomplete=%ld missing_total=%ld\n", nf, tot);
        }
        printf("EVAL tag=%s deg=%d N=%d n=%ld unc=%ld far_ordered=%ld far_unordered=%ld ctrl=\"%s\"\n", tag, deg, N, n, cur, cur * nG, cur * nG / 2, ctrl_desc);
        if (efile) write_edges(efile);
        return 0;
    }
    Mat bestA[MAXSYM]; Vec bestd[MAXSYM]; memcpy(bestA, A, sizeof A); memcpy(bestd, dv, sizeof dv);
    static Mat bestAd[MAXK][MAXSYM], saveAd[MAXK][MAXSYM]; static Vec bestdd[MAXK][MAXSYM], savedd[MAXK][MAXSYM];
    memcpy(bestAd, Ad, sizeof Ad); memcpy(bestdd, dd, sizeof dd);
    time_t t0 = time(NULL); long accepted = 0;
    Mat saveA[MAXSYM]; Vec saved[MAXSYM];
    for (long st = 0; st < steps && best > 0; st++) {
        if ((st & 63) == 0 && difftime(time(NULL), t0) > tlimit) break;
        double Temp = T0 * pow(T1 / T0, (double)st / (steps > 1 ? steps - 1 : 1));
        memcpy(saveA, A, sizeof A); memcpy(saved, dv, sizeof dv); memcpy(saveAd, Ad, sizeof Ad); memcpy(savedd, dd, sizeof dd);
        int sym = rnd(2 + ngrp);   // 0: pair a, 1: pair b, 2..: a plain group's constant
        int c = perdart ? rnd(k) : 0;
        // in per-dart mode the move acts on the dart (c, .) through Ad / dd; otherwise on the symbol chart A / dv (all darts)
        Mat *At = perdart ? &Ad[c][0] : &A[0]; Vec *dt = perdart ? &dd[c][0] : &dv[0];
        if (sym < 2) {
            int t = 2 * sym; int mv = rnd(10);
            if (mv < 4) { Mat inv; int j = rnd(M); for (int tries = 0; tries < 50; tries++) { for (int i = 0; i < M; i++) At[t].a[i][j] = rnd(P3); if (mat_inv(&At[t], &inv)) break; if (tries == 49) At[t] = perdart ? saveAd[c][t] : saveA[t]; } }
            else if (mv < 7) { int i0 = rnd(M), j0 = rnd(M), l = 1 + rnd(P3 - 1); if (i0 == j0) j0 = (j0 + 1) % M; for (int i = 0; i < M; i++) At[t].a[i][j0] = md(At[t].a[i][j0] + l * At[t].a[i][i0]); }
            else dt[t] = random_nonzero();
            if (perdart) derive_pair_d(c, t); else derive_pair(t);
        } else {
            int g = sym - 2, t = grp_sym[g], mv = rnd(10);
            if (mv < 6) dt[t] = random_vec(); else { int i = rnd(M); dt[t].v[i] = md(dt[t].v[i] + 1 + rnd(2)); }
            if (perdart) derive_plain_d(c, g); else derive_plain(g);
            if (perdart ? plain_clash_d() : plain_clash()) { memcpy(A, saveA, sizeof A); memcpy(dv, saved, sizeof dv); memcpy(Ad, saveAd, sizeof Ad); memcpy(dd, savedd, sizeof dd); continue; }
        }
        sync_charts();
        build_adj();
        long v = evaluate();
        if (v <= cur || rndu() < exp((curv - v) / Temp)) {
            cur = v; curv = v; accepted++;
            if (v < best) {
                best = v; memcpy(bestA, A, sizeof A); memcpy(bestd, dv, sizeof dv); memcpy(bestAd, Ad, sizeof Ad); memcpy(bestdd, dd, sizeof dd);
                fprintf(stderr, "step %ld T=%.3f best=%ld (%.0fs)\n", st, Temp, best, difftime(time(NULL), t0));
                save_chart(out, best, ctrl_desc);
            }
        } else { memcpy(A, saveA, sizeof A); memcpy(dv, saved, sizeof dv); memcpy(Ad, saveAd, sizeof Ad); memcpy(dd, savedd, sizeof dd); }
        if ((st % 5000) == 0) fprintf(stderr, "step %ld T=%.3f cur=%ld best=%ld acc=%ld (%.0fs)\n", st, Temp, cur, best, accepted, difftime(time(NULL), t0));
    }
    memcpy(A, bestA, sizeof A); memcpy(dv, bestd, sizeof dv); memcpy(Ad, bestAd, sizeof Ad); memcpy(dd, bestdd, sizeof dd);
    if (perdart) check_chart_d(); else { check_chart(); sync_charts(); }
    build_adj();
    long fin = evaluate();
    for (int sweep = 0; sweep < polish && fin > 0; sweep++) {
        long before = fin;
        // plain constants: every vector per plain group
        for (int g = 0; g < ngrp; g++) {
            int t = grp_sym[g]; Vec save = dv[t];
            for (int x = 0; x < F; x++) {
                dv[t] = index_vec(x); derive_plain(g);
                if (plain_clash()) continue;
                sync_charts(); build_adj(); long v = evaluate();
                if (v < fin) { fin = v; save = dv[t]; } else dv[t] = save;
            }
            dv[t] = save; derive_plain(g);
        }
        // line directions d_a, d_b: every nonzero vector
        for (int t = 0; t < KL; t += 2) {
            Vec save = dv[t];
            for (int x = 1; x < F; x++) {
                dv[t] = index_vec(x); derive_pair(t); sync_charts(); build_adj(); long v = evaluate();
                if (v < fin) { fin = v; save = dv[t]; } else { dv[t] = save; derive_pair(t); }
            }
        }
        // line matrices: every column replacement keeping invertibility
        for (int t = 0; t < KL; t += 2) for (int j = 0; j < M; j++) {
            Mat save = A[t];
            for (int x = 1; x < F; x++) {
                Vec col = index_vec(x); Mat inv;
                for (int i = 0; i < M; i++) A[t].a[i][j] = col.v[i];
                if (!mat_inv(&A[t], &inv)) { A[t] = save; continue; }
                derive_pair(t); sync_charts(); build_adj(); long v = evaluate();
                if (v < fin) { fin = v; save = A[t]; } else { A[t] = save; derive_pair(t); }
            }
        }
        sync_charts(); build_adj(); fin = evaluate();
        fprintf(stderr, "polish sweep %d: %ld -> %ld\n", sweep, before, fin);
        if (fin == before) break;
    }
    save_chart(out, fin, ctrl_desc);
    printf("RESULT tag=%s deg=%d N=%d n=%ld best=%ld accepted=%ld secs=%.0f saved=%s ctrl=\"%s\"\n", tag, deg, N, n, fin, accepted, difftime(time(NULL), t0), out, ctrl_desc);
    if (fin == 0 && efile) {
        write_edges(efile);
        printf("LIFT0 n=%ld deg=%d file=%s\n", n, deg, efile);
    }
    return 0;
}
