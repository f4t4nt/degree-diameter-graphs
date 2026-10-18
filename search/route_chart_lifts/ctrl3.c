// ctrl3.c: controller enumerator for DIAMETER-3 route-chart lifts over F_3^3 (27-point cells; RESULTS.md section 2).
// Controllers = regular lifts of a k-vertex Schreier base (alpha0, beta0 permutations, sigma0 involution: links and
// semi-edges) over a group table G (the ctrl5.c skeleton), symbols a, a^-1, b, b^-1 [, s], optionally plus plain
// symbols p_i = right multiplication by central involutions z_i of G (--plain z1,z2,..; they commute with the deck
// action and the symbols, so the pencil {p t u, t p u, t u p} lands in one fibre).
// Lift condition at x = 0 (chart-independent necessary condition, exact for the walk types that matter):
//   every ordered pair (root, v) at controller distance 0, 2 or 3 (and the plain neighbours (root, root z_i)) needs
//   an NB 3-walk of line symbols (full-rank word) or, with plain symbols, a 2-walk (t, u) of line symbols from
//   (root z_i) (the pencil of 3 parallel planes), i.e. A^3 + A > 0 for the line sub-controller, relaxed by Z A^2.
//   Line-symbol neighbours (distance 1) are covered by the edge pencils (automatic for the coplanar chart).
// With 5 line symbols the symbol chart has >= 4 bad words t s pi(t)^-1 (7 patterns); --pi tests them (the per-edge
// chart of Mizuno's construction is not bound by them, so the default is anyword = A^3 + A > 0).
// Usage: ctrl3 k alpha0 beta0 sigma0|- group.tbl [--plain z1,z2,..] [--sample S seed] [--maxout M] [--near T]
//        [--ctrl PREFIX] [--edges PREFIX] [--pi] [--first v1,v2,...]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 1024
#define MAXK 16
#define MAXW 80

static int k, nG, N, K = 5, KL = 5;   // K = number of line symbols
static int a0[MAXK], b0[MAXK], s0[MAXK], a0i[MAXK], b0i[MAXK];
static int *T, *inv, *invols, ninvols = 0;
static int va[MAXK], vb[MAXK], vs[MAXK];
static int nbr[MAXN][8];
static int symInv[5] = {1, 0, 3, 2, 4};
static int words[MAXW][3], NW = 0;
static int nplain = 0, plainz[4];
static int dist[MAXK][MAXN];
static uint64_t wmask[MAXK][MAXN][2];
static uint8_t pencil[MAXK][MAXN];
static int usepi = 0;
static int PI[7][4] = {{0,1,2,3},{1,0,2,3},{0,1,3,2},{2,1,0,3},{3,1,2,0},{0,2,1,3},{0,3,2,1}};
static uint64_t badmask[7][2];

static void build_words(void) {
    NW = 0;
    for (int w0 = 0; w0 < K; w0++) for (int w1 = 0; w1 < K; w1++) for (int w2 = 0; w2 < K; w2++) {
        if (w1 == symInv[w0] || w2 == symInv[w1]) continue;
        words[NW][0] = w0; words[NW][1] = w1; words[NW][2] = w2;
        if (K == 5) for (int p = 0; p < 7; p++) if (w1 == 4 && w0 < 4 && w2 == symInv[PI[p][w0]]) badmask[p][NW >> 6] |= 1ULL << (NW & 63);
        NW++;
    }
}

static int build_and_check(void) {
    static int alpha[MAXN], beta[MAXN], alphai[MAXN], betai[MAXN], sig[MAXN];
    for (int c = 0; c < k; c++) for (int g = 0; g < nG; g++) {
        int v = c * nG + g;
        alpha[v] = a0[c] * nG + T[g * nG + va[c]];
        beta[v] = b0[c] * nG + T[g * nG + vb[c]];
        if (K == 5) {
            int d = s0[c];
            if (d == c) sig[v] = c * nG + T[g * nG + vs[c]];
            else if (c < d) sig[v] = d * nG + T[g * nG + vs[c]];
            else sig[v] = d * nG + T[g * nG + inv[vs[d]]];
        }
    }
    for (int v = 0; v < N; v++) { alphai[alpha[v]] = v; betai[beta[v]] = v; }
    for (int v = 0; v < N; v++) {
        nbr[v][0] = alpha[v]; nbr[v][1] = alphai[v]; nbr[v][2] = beta[v]; nbr[v][3] = betai[v]; if (K == 5) nbr[v][4] = sig[v];
        for (int j = 0; j < nplain; j++) nbr[v][K + j] = (v / nG) * nG + T[(v % nG) * nG + plainz[j]];
    }
    int Kt = K + nplain;
    for (int c = 0; c < k; c++) {   // simple at the roots (G-invariant)
        int v = c * nG;
        for (int i = 0; i < Kt; i++) { if (nbr[v][i] == v) return 0; for (int j = i + 1; j < Kt; j++) if (nbr[v][i] == nbr[v][j]) return 0; }
    }
    static int queue[MAXN];
    for (int c = 0; c < k; c++) {
        int *d = dist[c];
        for (int v = 0; v < N; v++) d[v] = -1;
        int s = c * nG, qh = 0, qt = 0; d[s] = 0; queue[qt++] = s;
        while (qh < qt) {
            int x = queue[qh++];
            if (d[x] >= 3) continue;
            for (int i = 0; i < Kt; i++) { int y = nbr[x][i]; if (d[y] < 0) { d[y] = d[x] + 1; queue[qt++] = y; } }
        }
        if (qt < N) return 0;
    }
    for (int c = 0; c < k; c++) {
        memset(wmask[c], 0, N * sizeof(wmask[c][0])); memset(pencil[c], 0, N);
        int s = c * nG;
        for (int w = 0; w < NW; w++) {
            int x = s;
            for (int i = 0; i < 3; i++) x = nbr[x][words[w][i]];
            wmask[c][x][w >> 6] |= 1ULL << (w & 63);
        }
        for (int j = 0; j < nplain; j++) {
            int z = nbr[s][K + j];
            for (int t = 0; t < K; t++) for (int u = 0; u < K; u++) { if (u == symInv[t]) continue; pencil[c][nbr[nbr[z][t]][u]] = 1; }
        }
    }
    return 1;
}

// pairs (root c, v) needing coverage without a good 3-word / pencil (pattern p = -1: any word); returns count
static int uncovered(int p) {
    int cnt = 0;
    for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) {
        int d = dist[c][v];
        if (d == 1) {   // line-symbol neighbours are covered by edge pencils; plain neighbours are not
            int isline = 0; for (int i = 0; i < K; i++) if (nbr[c * nG][i] == v) isline = 1;
            if (isline) continue;
        }
        uint64_t any = 0;
        for (int i = 0; i < 2; i++) any |= wmask[c][v][i] & ~(p < 0 ? 0 : badmask[p][i]);
        if (!any && !pencil[c][v]) cnt++;
    }
    return cnt;
}

static int parse_perm(const char *s, int *out) {
    int n = 0; const char *p = s;
    if (!strcmp(s, "-")) return -1;
    while (*p) { out[n++] = atoi(p); while (*p && *p != ',') p++; if (*p == ',') p++; }
    return n;
}
static uint64_t rng_state = 88172645463325252ULL;
static uint64_t xr(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }
static void print_volt(FILE *f) {
    fprintf(f, "va="); for (int c = 0; c < k; c++) fprintf(f, "%d%s", va[c], c + 1 < k ? "," : "");
    fprintf(f, " vb="); for (int c = 0; c < k; c++) fprintf(f, "%d%s", vb[c], c + 1 < k ? "," : "");
    if (K == 5) { fprintf(f, " vs="); for (int c = 0; c < k; c++) fprintf(f, "%d%s", vs[c], c + 1 < k ? "," : ""); }
    if (nplain) { fprintf(f, " plain="); for (int j = 0; j < nplain; j++) fprintf(f, "%d%s", plainz[j], j + 1 < nplain ? "," : ""); }
}
static void write_ctrl(const char *fn, const char *tag) {
    FILE *f = fopen(fn, "w"); int Kt = K + nplain;
    fprintf(f, "%d %d\n", N, Kt);
    for (int t = 0; t < K; t++) fprintf(f, "%d ", symInv[t]);
    for (int j = 0; j < nplain; j++) fprintf(f, "%d ", K + j);
    fprintf(f, "\n");
    for (int t = 0; t < K; t++) fprintf(f, "1 ");
    for (int j = 0; j < nplain; j++) fprintf(f, "0 ");
    fprintf(f, "\n");
    for (int v = 0; v < N; v++) { for (int t = 0; t < Kt; t++) fprintf(f, "%d ", nbr[v][t]); fprintf(f, "\n"); }
    fprintf(f, "%d", k); for (int c = 0; c < k; c++) fprintf(f, " %d", c * nG); fprintf(f, "\n");
    fprintf(f, "# %s ", tag); print_volt(f); fprintf(f, "\n");
    fclose(f);
}
static void write_edges(const char *fn) {
    FILE *f = fopen(fn, "w"); int Kt = K + nplain;
    for (int v = 0; v < N; v++) for (int t = 0; t < Kt; t++) if (nbr[v][t] > v) fprintf(f, "%d %d\n", v, nbr[v][t]);
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: ctrl3 k alpha0 beta0 sigma0|- group.tbl [options]\n"); return 2; }
    k = atoi(argv[1]);
    if (parse_perm(argv[2], a0) != k || parse_perm(argv[3], b0) != k) { fprintf(stderr, "perm length\n"); return 2; }
    int ns = parse_perm(argv[4], s0);
    if (ns == -1) K = 4; else if (ns != k) { fprintf(stderr, "sigma0 length\n"); return 2; }
    KL = K;
    if (K == 5) for (int c = 0; c < k; c++) if (s0[s0[c]] != c) { fprintf(stderr, "sigma0 not an involution\n"); return 2; }
    for (int c = 0; c < k; c++) { a0i[a0[c]] = c; b0i[b0[c]] = c; }
    FILE *f = fopen(argv[5], "r"); if (!f) { perror("group"); return 2; }
    if (fscanf(f, "%d", &nG) != 1) return 2;
    T = malloc(nG * nG * sizeof(int)); inv = malloc(nG * sizeof(int)); invols = malloc(nG * sizeof(int));
    for (int i = 0; i < nG * nG; i++) if (fscanf(f, "%d", &T[i]) != 1) return 2;
    fclose(f);
    for (int g = 0; g < nG; g++) for (int h = 0; h < nG; h++) if (T[g * nG + h] == 0) inv[g] = h;
    for (int g = 1; g < nG; g++) if (inv[g] == g) invols[ninvols++] = g;
    N = k * nG; if (N > MAXN) { fprintf(stderr, "N too big\n"); return 2; }
    long sample = 0; int maxout = 20, nearT = -1; uint64_t seed = 1; const char *ctrlprefix = NULL, *edgeprefix = NULL;
    int nfirst = 0, firstvals[256];
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "--sample")) { sample = atol(argv[++i]); seed = atoll(argv[++i]); }
        else if (!strcmp(argv[i], "--maxout")) maxout = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--near")) nearT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--ctrl")) ctrlprefix = argv[++i];
        else if (!strcmp(argv[i], "--edges")) edgeprefix = argv[++i];
        else if (!strcmp(argv[i], "--pi")) usepi = 1;
        else if (!strcmp(argv[i], "--first")) nfirst = parse_perm(argv[++i], firstvals);
        else if (!strcmp(argv[i], "--plain")) { nplain = parse_perm(argv[++i], plainz); for (int j = 0; j < nplain; j++) { int z = plainz[j]; if (inv[z] != z) { fprintf(stderr, "plain %d not an involution\n", z); return 2; } for (int g = 0; g < nG; g++) if (T[g * nG + z] != T[z * nG + g]) { fprintf(stderr, "plain %d not central\n", z); return 2; } } }
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    rng_state ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    build_words();
    int intree_a[MAXK] = {0}, intree_b[MAXK] = {0}, intree_s[MAXK] = {0}, seen[MAXK] = {0}, q[MAXK], qh = 0, qt = 0;
    seen[0] = 1; q[qt++] = 0;
    while (qh < qt) {
        int c = q[qh++];
        int tg[5] = {a0[c], a0i[c], b0[c], b0i[c], K == 5 ? s0[c] : c};
        for (int s = 0; s < K; s++) {
            int d = tg[s]; if (seen[d]) continue;
            seen[d] = 1; q[qt++] = d;
            if (s == 0) intree_a[c] = 1; else if (s == 1) intree_a[d] = 1; else if (s == 2) intree_b[c] = 1; else if (s == 3) intree_b[d] = 1;
            else { intree_s[c] = 1; intree_s[d] = 1; }
        }
    }
    if (qt < k) { fprintf(stderr, "base not connected\n"); return 2; }
    int nfree = 0, freeslot[3 * MAXK], domain[3 * MAXK];
    for (int c = 0; c < k; c++) if (!intree_a[c]) { freeslot[nfree] = c; domain[nfree++] = nG; }
    for (int c = 0; c < k; c++) if (!intree_b[c]) { freeslot[nfree] = k + c; domain[nfree++] = nG; }
    if (K == 5) for (int c = 0; c < k; c++) {
        if (s0[c] == c) { if (ninvols == 0) { fprintf(stderr, "semi-edge but G has no involution\n"); return 2; } freeslot[nfree] = 2 * k + c; domain[nfree++] = ninvols; }
        else if (c < s0[c] && !intree_s[c]) { freeslot[nfree] = 2 * k + c; domain[nfree++] = nG; }
    }
    memset(va, 0, sizeof va); memset(vb, 0, sizeof vb); memset(vs, 0, sizeof vs);
    double total = 1; for (int i = 0; i < nfree; i++) total *= (i == 0 && nfirst) ? nfirst : domain[i];
    fprintf(stderr, "k=%d K=%d+%d nG=%d N=%d involutions=%d free slots=%d (%.3g assignments)%s\n", k, K, nplain, nG, N, ninvols, nfree, total, sample ? " sampling" : "");
    long tried = 0, valid = 0, pass = 0, near = 0, written = 0; int best = 1 << 30;
    int idx[3 * MAXK]; memset(idx, 0, sizeof idx);
    long iter = 0;
    while (1) {
        if (sample) {
            if (iter++ >= sample) break;
            for (int i = 0; i < nfree; i++) idx[i] = (i == 0 && nfirst) ? (int)(xr() % nfirst) : (int)(xr() % domain[i]);
        } else {
            if (iter++ > 0) {
                int i = 0;
                while (i < nfree) { int dom = (i == 0 && nfirst) ? nfirst : domain[i]; if (++idx[i] < dom) break; idx[i] = 0; i++; }
                if (i == nfree) break;
            }
        }
        for (int i = 0; i < nfree; i++) {
            int s = freeslot[i]; int val = (i == 0 && nfirst) ? firstvals[idx[i]] : idx[i];
            if (s < k) va[s] = val; else if (s < 2 * k) vb[s - k] = val;
            else { int c = s - 2 * k; vs[c] = (s0[c] == c) ? invols[val] : val; }
        }
        tried++;
        if (!build_and_check()) continue;
        valid++;
        int u = uncovered(-1); int bp = -1;
        if (u == 0 && usepi && K == 5) { int bu = 1 << 30; for (int p = 0; p < 7; p++) { int x = uncovered(p); if (x < bu) { bu = x; bp = p; } } u = bu; }
        if (u < best) best = u;
        if (u == 0) {
            pass++;
            if (pass <= maxout) { printf("PASS N=%d base=%s/%s/%s %s ", N, argv[2], argv[3], argv[4], argv[5]); print_volt(stdout); if (usepi) printf(" pi=%d", bp); printf("\n"); fflush(stdout); }
            if (written < maxout && (ctrlprefix || edgeprefix)) {
                char fn[512];
                if (ctrlprefix) { snprintf(fn, sizeof fn, "%s_%ld.ctrl", ctrlprefix, written); char tag[256]; snprintf(tag, sizeof tag, "base=%s/%s/%s %s", argv[2], argv[3], argv[4], argv[5]); write_ctrl(fn, tag); }
                if (edgeprefix) { snprintf(fn, sizeof fn, "%s_%ld.edges", edgeprefix, written); write_edges(fn); }
                written++;
            }
        } else if (nearT >= 0 && u <= nearT) {
            near++;
            if (near <= maxout) { printf("NEAR uncovered=%d N=%d base=%s/%s/%s %s ", u, N, argv[2], argv[3], argv[4], argv[5]); print_volt(stdout); printf("\n"); fflush(stdout); }
        }
    }
    printf("SUMMARY k=%d K=%d+%d nG=%d N=%d base=%s/%s/%s tried=%ld valid(simple,diam<=3)=%ld pass=%ld near=%ld best_uncovered=%d\n", k, K, nplain, nG, N, argv[2], argv[3], argv[4], tried, valid, pass, near, best == (1 << 30) ? -1 : best);
    return 0;
}
