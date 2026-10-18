// univ.c: chart-only SA for UNIVERSAL route charts over F_p^4 with K line symbols (pairs a, a^-1, b, b^-1, ...): minimise the
// number of reduced 4-words t1 t2 t3 t4 (t_{i+1} != t_i^-1) whose controllability columns [A4 A3 A2 d1 | A4 A3 d2 | A4 d3 | d4]
// are dependent (108 words for two pairs).  Universality (0 bad words) means every geodesic of a 4-regular controller of
// diameter 4 covers its target fibre, which with one pause symbol leaves only the chart-only conditions at distance
// <= 3.  usage: univ npairs seed steps restarts [--p P] [--out chart.txt] [--T0 t] [--T1 t]
// Output: per restart the best bad-word count; the best chart is written in exlift format (pairs, dims 1) plus one pause
// symbol line (M = I, random d) so that it can be used as --init for a PP..Z controller.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#define M 4
#define MAXK 8
static int P = 3, K;
typedef struct { int a[M][M]; } Mat; typedef struct { int v[M]; } Vec;
static Mat A[MAXK]; static Vec d[MAXK]; static int inv_[MAXK];
static uint64_t rs = 88172645463325252ULL;
static uint64_t xr(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static int rnd(int n) { return (int)(xr() % (uint64_t)n); }
static int md(int x) { x %= P; return x < 0 ? x + P : x; }
static int modinv(int a) { a = md(a); for (int i = 1; i < P; i++) if ((a * i) % P == 1) return i; return 0; }
static Vec mv(const Mat *m, const Vec *v) { Vec r; for (int i = 0; i < M; i++) { int s = 0; for (int j = 0; j < M; j++) s += m->a[i][j] * v->v[j]; r.v[i] = md(s); } return r; }
static int matinv(const Mat *m, Mat *out) {
    int W[M][2 * M];
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { W[i][j] = md(m->a[i][j]); W[i][M + j] = (i == j); }
    for (int c = 0; c < M; c++) {
        int piv = -1; for (int i = c; i < M; i++) if (W[i][c]) { piv = i; break; }
        if (piv < 0) return 0;
        if (piv != c) for (int j = 0; j < 2 * M; j++) { int t = W[c][j]; W[c][j] = W[piv][j]; W[piv][j] = t; }
        int iv = modinv(W[c][c]); for (int j = 0; j < 2 * M; j++) W[c][j] = md(W[c][j] * iv);
        for (int i = 0; i < M; i++) if (i != c && W[i][c]) { int f = W[i][c]; for (int j = 0; j < 2 * M; j++) W[i][j] = md(W[i][j] - f * W[c][j]); }
    }
    for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) out->a[i][j] = W[i][M + j];
    return 1;
}
static int rank4(Vec *c) {   // rank of 4 column vectors
    int W[M][M]; for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) W[i][j] = c[j].v[i];
    int rk = 0;
    for (int col = 0; col < M && rk < M; col++) {
        int piv = -1; for (int i = rk; i < M; i++) if (W[i][col]) { piv = i; break; }
        if (piv < 0) continue;
        for (int j = 0; j < M; j++) { int t = W[rk][j]; W[rk][j] = W[piv][j]; W[piv][j] = t; }
        int iv = modinv(W[rk][col]); for (int j = 0; j < M; j++) W[rk][j] = md(W[rk][j] * iv);
        for (int i = 0; i < M; i++) if (i != rk && W[i][col]) { int f = W[i][col]; for (int j = 0; j < M; j++) W[i][j] = md(W[i][j] - f * W[rk][j]); }
        rk++;
    }
    return rk;
}
static int is_zero(const Vec *v) { for (int i = 0; i < M; i++) if (v->v[i]) return 0; return 1; }
static void derive(int t) { Mat iv; matinv(&A[t], &iv); A[t + 1] = iv; d[t + 1] = mv(&iv, &d[t]); }
static Mat matmul(const Mat *x, const Mat *y) { Mat r; memset(&r, 0, sizeof r); for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) { int s = 0; for (int k = 0; k < M; k++) s += x->a[i][k] * y->a[k][j]; r.a[i][j] = md(s); } return r; }
// self-inverse line symbol s (exlift's convention): M_s^2 = I, M_s d_s = d_s: M_s = Pm diag(1, e1, e2, e3) Pm^-1 with Pm e_0 = d_s,
// e_i in {+1, -1} not all +1 (identity excluded)
static void random_invol(int t) {
    Mat Pm, Pi, Dg; memset(&Dg, 0, sizeof Dg);
    do { for (int j = 0; j < M; j++) for (int i = 0; i < M; i++) Pm.a[i][j] = rnd(P); for (int i = 0; i < M; i++) Pm.a[i][0] = d[t].v[i]; } while (!matinv(&Pm, &Pi));
    int all1 = 1; for (int j = 0; j < M; j++) { int e = (j == 0) ? 1 : (rnd(2) ? 1 : P - 1); Dg.a[j][j] = e; if (e != 1) all1 = 0; }
    if (all1) { int j = 1 + rnd(M - 1); Dg.a[j][j] = P - 1; }
    Mat t1 = matmul(&Pm, &Dg); A[t] = matmul(&t1, &Pi);
}
static void random_pair(int t) { Mat iv; do { for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) A[t].a[i][j] = rnd(P); } while (!matinv(&A[t], &iv)); do { for (int i = 0; i < M; i++) d[t].v[i] = rnd(P); } while (is_zero(&d[t])); derive(t); }
static int bad_words(void) {
    int bad = 0;
    for (int t1 = 0; t1 < K; t1++) for (int t2 = 0; t2 < K; t2++) { if (t2 == inv_[t1]) continue;
        for (int t3 = 0; t3 < K; t3++) { if (t3 == inv_[t2]) continue;
            for (int t4 = 0; t4 < K; t4++) { if (t4 == inv_[t3]) continue;
                Vec c[4]; c[3] = d[t4]; c[2] = mv(&A[t4], &d[t3]); Vec x = mv(&A[t3], &d[t2]); c[1] = mv(&A[t4], &x); x = mv(&A[t2], &d[t1]); x = mv(&A[t3], &x); c[0] = mv(&A[t4], &x);
                if (rank4(c) < M) bad++;
            } } }
    return bad;
}
static void save(const char *fn, int bad) {
    FILE *f = fopen(fn, "w"); if (!f) return;
    fprintf(f, "# bad=%d p=%d m=4 (univ.c; last line = pause symbol with M = I)\n", bad, P);
    for (int t = 0; t < K; t++) { fprintf(f, "A"); for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) fprintf(f, " %d", A[t].a[i][j]); fprintf(f, " d"); for (int i = 0; i < M; i++) fprintf(f, " %d", d[t].v[i]); fprintf(f, "\n"); }
    fprintf(f, "A"); for (int i = 0; i < M; i++) for (int j = 0; j < M; j++) fprintf(f, " %d", i == j); fprintf(f, " d"); Vec z; do { for (int i = 0; i < M; i++) z.v[i] = rnd(P); } while (is_zero(&z)); for (int i = 0; i < M; i++) fprintf(f, " %d", z.v[i]); fprintf(f, "\n");
    fclose(f);
}
int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: univ npairs seed steps restarts [--p P] [--out f] [--T0 t] [--T1 t]\n"); return 2; }
    int npairs = atoi(argv[1]); uint64_t seed = strtoull(argv[2], NULL, 10); long steps = atol(argv[3]); int restarts = atoi(argv[4]);
    const char *out = "univ_best.chart"; double T0 = 2.0, T1 = 0.05; int ninv = 0;
    for (int i = 5; i < argc; i++) { if (!strcmp(argv[i], "--p")) P = atoi(argv[++i]); else if (!strcmp(argv[i], "--out")) out = argv[++i]; else if (!strcmp(argv[i], "--T0")) T0 = atof(argv[++i]); else if (!strcmp(argv[i], "--T1")) T1 = atof(argv[++i]); else if (!strcmp(argv[i], "--inv")) ninv = atoi(argv[++i]); }
    K = 2 * npairs; for (int t = 0; t < K; t += 2) { inv_[t] = t + 1; inv_[t + 1] = t; }
    for (int i = 0; i < ninv; i++) { inv_[K] = K; K++; }   // self-inverse line symbols s (M^2 = I, M d = d), the (15,4) alphabet with --inv 1
    int nfree = npairs + ninv;
    rs ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    int bestall = 1 << 30; long hist[512] = {0};
    for (int r = 0; r < restarts; r++) {
        for (int t = 0; t < 2 * npairs; t += 2) random_pair(t);
        for (int t = 2 * npairs; t < K; t++) { do { for (int i = 0; i < M; i++) d[t].v[i] = rnd(P); } while (is_zero(&d[t])); random_invol(t); }
        int cur = bad_words(), best = cur; Mat sA[MAXK]; Vec sd[MAXK];
        for (long st = 0; st < steps; st++) {
            double Temp = T0 * pow(T1 / T0, (double)st / (double)(steps - 1));
            memcpy(sA, A, sizeof A); memcpy(sd, d, sizeof d);
            int f = rnd(nfree), mvk = rnd(10);
            if (f >= npairs) {   // involution symbol: redraw M fixing d, or move d (then redraw M)
                int t = 2 * npairs + (f - npairs);
                if (mvk < 6) random_invol(t);
                else { int j = rnd(M); d[t].v[j] = md(d[t].v[j] + 1 + rnd(P - 1)); if (is_zero(&d[t])) d[t] = sd[t]; random_invol(t); }
            } else {
            int t = 2 * f;
            if (mvk < 3) { Mat iv; int j = rnd(M); for (int tries = 0; tries < 50; tries++) { for (int i = 0; i < M; i++) A[t].a[i][j] = rnd(P); if (matinv(&A[t], &iv)) break; if (tries == 49) A[t] = sA[t]; } }
            else if (mvk < 6) { int i0 = rnd(M), j0 = rnd(M), l = 1 + rnd(P - 1); if (i0 == j0) j0 = (j0 + 1) % M; for (int i = 0; i < M; i++) A[t].a[i][j0] = md(A[t].a[i][j0] + l * A[t].a[i][i0]); }
            else if (mvk < 9) { int j = rnd(M); d[t].v[j] = md(d[t].v[j] + 1 + rnd(P - 1)); if (is_zero(&d[t])) d[t] = sd[t]; }
            else { do { for (int i = 0; i < M; i++) d[t].v[i] = rnd(P); } while (is_zero(&d[t])); }
            derive(t);
            }
            int v = bad_words();
            if (v <= cur || (xr() >> 11) * (1.0 / 9007199254740992.0) < exp((cur - v) / Temp)) { cur = v; if (v < best) { best = v; if (best < bestall) { bestall = best; save(out, best); } if (best == 0) break; } }
            else { memcpy(A, sA, sizeof A); memcpy(d, sd, sizeof d); }
        }
        if (best < 512) hist[best]++;
        printf("restart %d best bad words %d\n", r, best); fflush(stdout);
    }
    printf("SUMMARY univ npairs=%d p=%d words=%d bestall=%d hist", npairs, P, K * (K - 1) * (K - 1) * (K - 1), bestall);
    for (int b = 0; b < 512; b++) if (hist[b]) printf(" %d:%ld", b, hist[b]); printf("\n");
    return 0;
}
