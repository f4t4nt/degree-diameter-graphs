// Controller search: regular lifts of a k-vertex Schreier base (two permutations alpha0, beta0 on k points)
// over a group G (multiplication table), filtered for the route-chart-lift diameter-4 conditions:
//   * simple 4-regular lift, diameter <= 4;
//   * every ordered pair at controller distance 1, 3 or 4 has a non-backtracking walk of length exactly 4 whose
//     word (over a, A=a^-1, b, B=b^-1) is controllable for the chart (word not in the chart's BAD set).
//   (pairs at distance 0 and 2 are covered by the chart-only conditions C0/C2, see chart_cond.py)
// Usage: ctrl_lift k alpha0 beta0 group.tbl badsets.txt [--sample S seed] [--maxout M] [--nearmiss T]
//   alpha0, beta0: comma-separated images, e.g. 1,2,3,0
//   group.tbl: text: n, then n*n integers T[g*n+h] = g*h (identity = 0)
//   badsets.txt: one chart per line: count followed by word ids (lexicographic index among the 108 reduced words)
// Enumeration: voltages on a spanning tree of the base are the identity; all other voltages run over G
// (exhaustive) or are sampled uniformly (--sample).  Survivors are printed as "HIT chart=<i> va=.. vb=..";
// near misses (<= T uncovered pairs for the best chart) as "NEAR ...".
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 512
#define MAXK 16
#define NW 108

static int k, nG, N;
static int a0[MAXK], b0[MAXK], a0i[MAXK], b0i[MAXK];
static int *T, *inv;
static int va[MAXK], vb[MAXK];
static int alpha[MAXN], beta[MAXN], alphai[MAXN], betai[MAXN];
static int nbr[MAXN][4];   // symbol s: 0=a (alpha), 1=A (alpha^-1), 2=b, 3=B
static int words[NW][4];
static int nbad = 0; static uint64_t badlo[64], badhi[64];
static int need0 = 0, need2 = 0;  // also require good NB-4 walks for pairs at distance 0 / 2
// exact coverage mode (--chart FILE, rclift chart format): F <= 256 fibre points, K = 4 charts (a, A, b, B)
static int exact = 0, F = 0, pF = 0, mF = 0, KF = 0;
static int *clin[4], ctsz[4], *ctset[4];
static int *addF;
static uint64_t cover[MAXN][4];
static int vaddF(int a, int b) { int r = 0, pw = 1; for (int i = 0; i < mF; i++) { int da = a % pF, db = b % pF; a /= pF; b /= pF; r += ((da + db) % pF) * pw; pw *= pF; } return r; }
static void load_chart(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror("chart"); exit(2); }
    if (fscanf(f, "%d %d %d %d", &F, &pF, &mF, &KF) != 4 || KF != 4 || F > 256) { fprintf(stderr, "bad chart header\n"); exit(2); }
    for (int k = 0; k < 4; k++) {
        clin[k] = malloc(F * sizeof(int));
        for (int x = 0; x < F; x++) if (fscanf(f, "%d", &clin[k][x]) != 1) exit(2);
        if (fscanf(f, "%d", &ctsz[k]) != 1) exit(2);
        ctset[k] = malloc(ctsz[k] * sizeof(int));
        for (int j = 0; j < ctsz[k]; j++) if (fscanf(f, "%d", &ctset[k][j]) != 1) exit(2);
    }
    fclose(f);
    addF = malloc(F * F * sizeof(int));
    for (int a = 0; a < F; a++) for (int b = 0; b < F; b++) addF[a * F + b] = vaddF(a, b);
    exact = 1;
}
// image of a point set (bitmask over F) under chart k: {lin[x] + tau}
static void chart_image(const uint64_t *in, int k, uint64_t *out) {
    out[0] = out[1] = out[2] = out[3] = 0;
    for (int x = 0; x < F; x++) if (in[x >> 6] >> (x & 63) & 1) {
        int lx = clin[k][x];
        for (int j = 0; j < ctsz[k]; j++) { int y = addF[lx * F + ctset[k][j]]; out[y >> 6] |= 1ULL << (y & 63); }
    }
}
// exact: number of (root fibre vertex, target vertex) pairs not fully covered by walks of length <= 4
static int exact_uncovered(void) {
    static uint64_t cur[340][4]; static int curv[340];
    static uint64_t nxt[340][4]; static int nxtv[340];
    int total = 0;
    uint64_t full[4] = {0, 0, 0, 0};
    for (int x = 0; x < F; x++) full[x >> 6] |= 1ULL << (x & 63);
    for (int c = 0; c < k; c++) {
        for (int v = 0; v < N; v++) cover[v][0] = cover[v][1] = cover[v][2] = cover[v][3] = 0;
        int s = c * nG;
        cover[s][0] = 1;  // the root (c, e) with x = 0 reaches itself
        int ncur = 1; cur[0][0] = 1; cur[0][1] = cur[0][2] = cur[0][3] = 0; curv[0] = s;
        for (int L = 1; L <= 4; L++) {
            int nn = 0;
            for (int i = 0; i < ncur; i++) for (int sym = 0; sym < 4; sym++) {
                int w = nbr[curv[i]][sym];
                chart_image(cur[i], sym, nxt[nn]); nxtv[nn] = w;
                for (int j = 0; j < 4; j++) cover[w][j] |= nxt[nn][j];
                nn++;
            }
            ncur = nn; memcpy(cur, nxt, nn * sizeof(cur[0])); memcpy(curv, nxtv, nn * sizeof(int));
        }
        for (int v = 0; v < N; v++) for (int j = 0; j < 4; j++) if ((cover[v][j] & full[j]) != full[j]) { total++; break; }
    }
    return total;
}
static int dist[MAXK][MAXN];
static uint64_t wlo[MAXK][MAXN], whi[MAXK][MAXN];

static void build_words(void) {
    int inv4[4] = {1, 0, 3, 2}, n = 0;
    for (int w0 = 0; w0 < 4; w0++) for (int w1 = 0; w1 < 4; w1++) for (int w2 = 0; w2 < 4; w2++) for (int w3 = 0; w3 < 4; w3++) {
        if (w1 == inv4[w0] || w2 == inv4[w1] || w3 == inv4[w2]) continue;
        words[n][0] = w0; words[n][1] = w1; words[n][2] = w2; words[n][3] = w3; n++;
    }
    if (n != NW) { fprintf(stderr, "word count %d\n", n); exit(1); }
}

// returns 0 if invalid (not simple / disconnected / diameter > 4), else fills dist and word masks
static int build_and_check(void) {
    for (int c = 0; c < k; c++) for (int g = 0; g < nG; g++) {
        int v = c * nG + g;
        alpha[v] = a0[c] * nG + T[g * nG + va[c]];
        beta[v] = b0[c] * nG + T[g * nG + vb[c]];
    }
    for (int v = 0; v < N; v++) { alphai[alpha[v]] = v; betai[beta[v]] = v; }
    // simplicity: 4 distinct neighbours, none equal to v; by fibre-transitivity check only roots (c, e)... but
    // the condition is G-invariant, so checking the k roots suffices
    for (int c = 0; c < k; c++) {
        int v = c * nG;
        int nb[4] = {alpha[v], alphai[v], beta[v], betai[v]};
        for (int i = 0; i < 4; i++) { if (nb[i] == v) return 0; for (int j = i + 1; j < 4; j++) if (nb[i] == nb[j]) return 0; }
    }
    for (int v = 0; v < N; v++) { nbr[v][0] = alpha[v]; nbr[v][1] = alphai[v]; nbr[v][2] = beta[v]; nbr[v][3] = betai[v]; }
    // BFS from the k roots
    static int queue[MAXN];
    for (int c = 0; c < k; c++) {
        int *d = dist[c];
        for (int v = 0; v < N; v++) d[v] = -1;
        int s = c * nG, qh = 0, qt = 0; d[s] = 0; queue[qt++] = s;
        while (qh < qt) {
            int x = queue[qh++];
            if (d[x] >= 4) continue;
            for (int i = 0; i < 4; i++) { int y = nbr[x][i]; if (d[y] < 0) { d[y] = d[x] + 1; queue[qt++] = y; } }
        }
        if (qt < N) return 0;  // some vertex at distance > 4 (or disconnected)
    }
    // NB-4 word masks from the roots
    for (int c = 0; c < k; c++) {
        memset(wlo[c], 0, N * sizeof(uint64_t)); memset(whi[c], 0, N * sizeof(uint64_t));
        int s = c * nG;
        for (int w = 0; w < NW; w++) {
            int x = s;
            for (int i = 0; i < 4; i++) x = nbr[x][words[w][i]];
            if (w < 64) wlo[c][x] |= 1ULL << w; else whi[c][x] |= 1ULL << (w - 64);
        }
    }
    return 1;
}

// number of ordered pairs (root c, v) at distance 1,3,4 without a good NB-4 word for chart j (-1: any word)
static int uncovered(int j) {
    int cnt = 0;
    uint64_t blo = j < 0 ? 0 : badlo[j], bhi = j < 0 ? 0 : badhi[j];
    for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) {
        int d = dist[c][v];
        if (d == 1 || d == 3 || d == 4 || (need0 && d == 0) || (need2 && d == 2)) {
            if (((wlo[c][v] & ~blo) | (whi[c][v] & ~bhi)) == 0) cnt++;
        }
    }
    return cnt * nG;  // each root represents nG vertices of its fibre
}

static int parse_perm(const char *s, int *out) {
    int n = 0; const char *p = s;
    while (*p) { out[n++] = atoi(p); while (*p && *p != ',') p++; if (*p == ',') p++; }
    return n;
}

static uint64_t rng_state = 88172645463325252ULL;
static uint64_t xr(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage\n"); return 2; }
    k = atoi(argv[1]);
    if (parse_perm(argv[2], a0) != k || parse_perm(argv[3], b0) != k) { fprintf(stderr, "perm length\n"); return 2; }
    for (int c = 0; c < k; c++) { a0i[a0[c]] = c; b0i[b0[c]] = c; }
    FILE *f = fopen(argv[4], "r"); if (!f) { perror("group"); return 2; }
    if (fscanf(f, "%d", &nG) != 1) return 2;
    T = malloc(nG * nG * sizeof(int)); inv = malloc(nG * sizeof(int));
    for (int i = 0; i < nG * nG; i++) if (fscanf(f, "%d", &T[i]) != 1) return 2;
    fclose(f);
    for (int g = 0; g < nG; g++) for (int h = 0; h < nG; h++) if (T[g * nG + h] == 0) inv[g] = h;
    N = k * nG; if (N > MAXN) { fprintf(stderr, "N too big\n"); return 2; }
    f = fopen(argv[5], "r"); if (!f) { perror("badsets"); return 2; }
    int cnt;
    while (fscanf(f, "%d", &cnt) == 1) {
        badlo[nbad] = badhi[nbad] = 0;
        for (int i = 0; i < cnt; i++) { int w; if (fscanf(f, "%d", &w) != 1) return 2; if (w < 64) badlo[nbad] |= 1ULL << w; else badhi[nbad] |= 1ULL << (w - 64); }
        nbad++;
    }
    fclose(f);
    long sample = 0; int maxout = 50, nearT = -1; uint64_t seed = 1;
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "--sample")) { sample = atol(argv[++i]); seed = atoll(argv[++i]); }
        else if (!strcmp(argv[i], "--maxout")) maxout = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--nearmiss")) nearT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--need0")) need0 = 1;
        else if (!strcmp(argv[i], "--chart")) load_chart(argv[++i]);
        else if (!strcmp(argv[i], "--need2")) need2 = 1;
    }
    rng_state ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    build_words();
    // spanning tree of the base: BFS over darts (a: c->a0[c], A: c->a0i[c], b, B); tree darts get identity voltage.
    // free voltage slots: va[c] for a-darts not in tree, vb[c] for b-darts not in tree.
    int intree_a[MAXK] = {0}, intree_b[MAXK] = {0}, seen[MAXK] = {0}, q[MAXK], qh = 0, qt = 0;
    seen[0] = 1; q[qt++] = 0;
    while (qh < qt) {
        int c = q[qh++];
        int tg[4] = {a0[c], a0i[c], b0[c], b0i[c]};
        for (int s = 0; s < 4; s++) {
            int d = tg[s]; if (seen[d]) continue;
            seen[d] = 1; q[qt++] = d;
            if (s == 0) intree_a[c] = 1; else if (s == 1) intree_a[d] = 1; else if (s == 2) intree_b[c] = 1; else intree_b[d] = 1;
        }
    }
    if (qt < k) { fprintf(stderr, "base not connected\n"); return 2; }
    int nfree = 0, *freeslot = malloc(2 * k * sizeof(int));
    for (int c = 0; c < k; c++) if (!intree_a[c]) freeslot[nfree++] = c;          // a-slot c
    for (int c = 0; c < k; c++) if (!intree_b[c]) freeslot[nfree++] = k + c;      // b-slot c
    memset(va, 0, sizeof va); memset(vb, 0, sizeof vb);
    double total = 1; for (int i = 0; i < nfree; i++) total *= nG;
    fprintf(stderr, "k=%d nG=%d N=%d free voltages=%d (%.3g assignments) charts=%d%s\n", k, nG, N, nfree, total, nbad, sample ? " sampling" : "");
    long tried = 0, valid = 0, hits = 0, near = 0, anyword_ok = 0;
    int best = 1 << 30;
    int idx[2 * MAXK]; memset(idx, 0, sizeof idx);
    long iter = 0;
    while (1) {
        if (sample) {
            if (iter++ >= sample) break;
            for (int i = 0; i < nfree; i++) idx[i] = xr() % nG;
        } else {
            if (iter++ > 0) {
                int i = 0;
                while (i < nfree) { if (++idx[i] < nG) break; idx[i] = 0; i++; }
                if (i == nfree) break;
            }
        }
        for (int i = 0; i < nfree; i++) { int s = freeslot[i]; if (s < k) va[s] = idx[i]; else vb[s - k] = idx[i]; }
        tried++;
        if (!build_and_check()) continue;
        valid++;
        if (exact) {
            int ue = exact_uncovered();
            if (ue < best) best = ue;
            if (ue == 0) {
                hits++;
                if (hits <= maxout) {
                    printf("HIT exact va=");
                    for (int c = 0; c < k; c++) printf("%d%s", va[c], c + 1 < k ? "," : "");
                    printf(" vb=");
                    for (int c = 0; c < k; c++) printf("%d%s", vb[c], c + 1 < k ? "," : "");
                    printf("\n"); fflush(stdout);
                }
            } else if (nearT >= 0 && ue <= nearT) {
                near++;
                if (near <= maxout) {
                    printf("NEAR exact uncovered=%d va=", ue);
                    for (int c = 0; c < k; c++) printf("%d%s", va[c], c + 1 < k ? "," : "");
                    printf(" vb=");
                    for (int c = 0; c < k; c++) printf("%d%s", vb[c], c + 1 < k ? "," : "");
                    printf("\n"); fflush(stdout);
                }
            }
            continue;
        }
        int u_any = uncovered(-1);
        if (u_any == 0) anyword_ok++;
        int bestj = -1, bestu = 1 << 30;
        for (int j = 0; j < nbad; j++) { int u = uncovered(j); if (u < bestu) { bestu = u; bestj = j; } }
        if (nbad == 0) { bestu = u_any; bestj = -1; }
        if (bestu < best) best = bestu;
        if (bestu == 0) {
            hits++;
            if (hits <= maxout) {
                printf("HIT chart=%d va=", bestj);
                for (int c = 0; c < k; c++) printf("%d%s", va[c], c + 1 < k ? "," : "");
                printf(" vb=");
                for (int c = 0; c < k; c++) printf("%d%s", vb[c], c + 1 < k ? "," : "");
                printf("\n"); fflush(stdout);
            }
        } else if (nearT >= 0 && bestu <= nearT) {
            near++;
            if (near <= maxout) {
                printf("NEAR uncovered=%d anyword=%d chart=%d va=", bestu, u_any, bestj);
                for (int c = 0; c < k; c++) printf("%d%s", va[c], c + 1 < k ? "," : "");
                printf(" vb=");
                for (int c = 0; c < k; c++) printf("%d%s", vb[c], c + 1 < k ? "," : "");
                printf("\n"); fflush(stdout);
            }
        }
    }
    printf("SUMMARY k=%d nG=%d N=%d tried=%ld valid(simple,diam<=4)=%ld anyword_ok=%ld hits=%ld near=%ld best_uncovered=%d\n", k, nG, N, tried, valid, anyword_ok, hits, near, best == (1 << 30) ? -1 : best);
    return 0;
}
