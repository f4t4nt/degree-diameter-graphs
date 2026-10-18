// Controller search for route-chart lifts with an odd-degree (5 symbols: a, a^-1, b, b^-1, s) or even-degree
// (4 symbols) alphabet.  Controllers = regular lifts of a k-vertex Schreier base (alpha0, beta0, sigma0) over a
// group G (multiplication table); sigma0 is an involution of the base: c <-> sigma0[c] is a link carrying one free
// voltage, a fixed point c is a semi-edge whose voltage runs over the involutions of G.
// Filters (in order, each cheap):
//   simple 5-regular lift, diameter <= 4 from the k roots (c, e);
//   anyword: every ordered pair (root, v) at distance 1, 3, 4 (plus 0 / 2 with --need0 / --need2) has a
//            non-backtracking walk of length exactly 4;
//   word filter (--bad FILE): for each chart's bad-word set, every such pair has an NB-4 walk with a good word;
//   exact (--chart FILE): full route-chart lift (fibre F = p^m points, K charts in rclift format) BFS from the k*F
//            roots (c, e, x) with early exit: 0 bad pairs = diameter-4 lift.
// Usage: ctrl5 k alpha0 beta0 sigma0 group.tbl [--bad badsets.txt] [--chart chart.cht] [--sample S seed]
//        [--maxout M] [--nearmiss T] [--need0] [--need2] [--first v1,v2,...] [--only v1,v2,...] [--dump FILE]
//   alpha0, beta0, sigma0: comma-separated images (sigma0 = "-" for the 4-symbol alphabet)
//   --first: restrict the first free voltage slot to the listed values (orbit representatives under Aut(G))
//   --only: evaluate the single assignment with these values in the free slots (group indices, in slot order:
//           a-slots, b-slots, s-slots of the base vertices not fixed by the spanning tree)
//   --dump: append one line per anyword-OK controller: voltages + word masks of every (root, v) pair (for chart
//           co-design: wordhits format "va.. vb.. vs.. | per pair: c v dist mask0 .. mask4")
// Words: reduced words of length 4 over K symbols, lexicographic index (K (K-1)^3 = 320 for K = 5, 108 for K = 4).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 1024
#define MAXK 16
#define MAXW 320
#define MAXF 256
#define MAXLIFT (MAXN * 81)

static int k, nG, N, K = 5;
static int a0[MAXK], b0[MAXK], s0[MAXK], a0i[MAXK], b0i[MAXK];
static int *T, *inv, *invols, ninvols = 0;
static int va[MAXK], vb[MAXK], vs[MAXK];
static int nbr[MAXN][5];       // symbol s: 0=a, 1=a^-1, 2=b, 3=b^-1, 4=s
static int symInv[5] = {1, 0, 3, 2, 4};
static int words[MAXW][4], NW = 0;
static int nbad = 0; static uint64_t badmask[64][5];
static int need0 = 0, need2 = 0;
static int dist[MAXK][MAXN];
static uint64_t wmask[MAXK][MAXN][5];
// exact lift
static int exact = 0, F = 0, pF = 0, mF = 0, KF = 0;
static int *clin[5], ctsz[5], *ctset[5];
static int *addF;
static int *ladj, *ldist; static long *lqueue; static int ldeg = 0;

static int vaddF(int a, int b) { int r = 0, pw = 1; for (int i = 0; i < mF; i++) { int da = a % pF, db = b % pF; a /= pF; b /= pF; r += ((da + db) % pF) * pw; pw *= pF; } return r; }
static void load_chart(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror("chart"); exit(2); }
    if (fscanf(f, "%d %d %d %d", &F, &pF, &mF, &KF) != 4 || KF != K || F > MAXF) { fprintf(stderr, "bad chart header (K=%d, need %d)\n", KF, K); exit(2); }
    for (int s = 0; s < K; s++) {
        clin[s] = malloc(F * sizeof(int));
        for (int x = 0; x < F; x++) if (fscanf(f, "%d", &clin[s][x]) != 1) exit(2);
        if (fscanf(f, "%d", &ctsz[s]) != 1) exit(2);
        ctset[s] = malloc(ctsz[s] * sizeof(int));
        for (int j = 0; j < ctsz[s]; j++) if (fscanf(f, "%d", &ctset[s][j]) != 1) exit(2);
        ldeg += ctsz[s];
    }
    fclose(f);
    addF = malloc(F * F * sizeof(int));
    for (int a = 0; a < F; a++) for (int b = 0; b < F; b++) addF[a * F + b] = vaddF(a, b);
    exact = 1;
}
// exact lift BFS from the k*F roots; returns number of uncovered (root, target) pairs (x nG = lift bad pairs),
// stopping early after `limit` uncovered pairs
static long exact_bad(long limit) {
    long n = (long)N * F;
    if (!ladj) { ladj = malloc((size_t)n * ldeg * sizeof(int)); ldist = malloc(n * sizeof(int)); lqueue = malloc(n * sizeof(long)); }
    // adjacency: vertex (v, x) index v*F + x; neighbours via each symbol
    for (int v = 0; v < N; v++) for (int x = 0; x < F; x++) {
        long u = (long)v * F + x; int *row = ladj + u * ldeg; int c = 0;
        for (int s = 0; s < K; s++) { int w = nbr[v][s]; int lx = clin[s][x]; for (int j = 0; j < ctsz[s]; j++) row[c++] = w * F + addF[lx * F + ctset[s][j]]; }
    }
    long bad = 0;
    for (int c = 0; c < k; c++) for (int x0 = 0; x0 < F; x0++) {
        long s = (long)c * nG * F + x0;
        for (long i = 0; i < n; i++) ldist[i] = -1;
        long qh = 0, qt = 0; ldist[s] = 0; lqueue[qt++] = s;
        while (qh < qt) {
            long x = lqueue[qh++];
            if (ldist[x] >= 4) continue;
            int *row = ladj + x * ldeg;
            for (int i = 0; i < ldeg; i++) { long y = row[i]; if (ldist[y] < 0) { ldist[y] = ldist[x] + 1; lqueue[qt++] = y; } }
        }
        bad += n - qt;
        if (bad > limit) return bad;
    }
    return bad;
}

static void build_words(void) {
    NW = 0;
    for (int w0 = 0; w0 < K; w0++) for (int w1 = 0; w1 < K; w1++) for (int w2 = 0; w2 < K; w2++) for (int w3 = 0; w3 < K; w3++) {
        if (w1 == symInv[w0] || w2 == symInv[w1] || w3 == symInv[w2]) continue;
        words[NW][0] = w0; words[NW][1] = w1; words[NW][2] = w2; words[NW][3] = w3; NW++;
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
    for (int v = 0; v < N; v++) { nbr[v][0] = alpha[v]; nbr[v][1] = alphai[v]; nbr[v][2] = beta[v]; nbr[v][3] = betai[v]; if (K == 5) nbr[v][4] = sig[v]; }
    // simplicity at the k roots (G-invariant condition)
    for (int c = 0; c < k; c++) {
        int v = c * nG;
        for (int i = 0; i < K; i++) { if (nbr[v][i] == v) return 0; for (int j = i + 1; j < K; j++) if (nbr[v][i] == nbr[v][j]) return 0; }
    }
    static int queue[MAXN];
    for (int c = 0; c < k; c++) {
        int *d = dist[c];
        for (int v = 0; v < N; v++) d[v] = -1;
        int s = c * nG, qh = 0, qt = 0; d[s] = 0; queue[qt++] = s;
        while (qh < qt) {
            int x = queue[qh++];
            if (d[x] >= 4) continue;
            for (int i = 0; i < K; i++) { int y = nbr[x][i]; if (d[y] < 0) { d[y] = d[x] + 1; queue[qt++] = y; } }
        }
        if (qt < N) return 0;
    }
    for (int c = 0; c < k; c++) {
        memset(wmask[c], 0, N * sizeof(wmask[c][0]));
        int s = c * nG;
        for (int w = 0; w < NW; w++) {
            int x = s;
            for (int i = 0; i < 4; i++) x = nbr[x][words[w][i]];
            wmask[c][x][w >> 6] |= 1ULL << (w & 63);
        }
    }
    return 1;
}

static int needed(int d) { return d == 1 || d == 3 || d == 4 || (need0 && d == 0) || (need2 && d == 2); }

// pairs (root c, v) that need an NB-4 walk and have no good word for chart j (-1: any word); returns count x nG
static int uncovered(int j) {
    int cnt = 0;
    for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) {
        if (!needed(dist[c][v])) continue;
        uint64_t any = 0;
        for (int i = 0; i < 5; i++) any |= wmask[c][v][i] & ~(j < 0 ? 0 : badmask[j][i]);
        if (!any) cnt++;
    }
    return cnt * nG;
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
    fprintf(f, "va=");
    for (int c = 0; c < k; c++) fprintf(f, "%d%s", va[c], c + 1 < k ? "," : "");
    fprintf(f, " vb=");
    for (int c = 0; c < k; c++) fprintf(f, "%d%s", vb[c], c + 1 < k ? "," : "");
    if (K == 5) { fprintf(f, " vs="); for (int c = 0; c < k; c++) fprintf(f, "%d%s", vs[c], c + 1 < k ? "," : ""); }
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: ctrl5 k alpha0 beta0 sigma0 group.tbl [options]\n"); return 2; }
    k = atoi(argv[1]);
    if (parse_perm(argv[2], a0) != k || parse_perm(argv[3], b0) != k) { fprintf(stderr, "perm length\n"); return 2; }
    int ns = parse_perm(argv[4], s0);
    if (ns == -1) K = 4; else if (ns != k) { fprintf(stderr, "sigma0 length\n"); return 2; }
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
    long sample = 0; int maxout = 50, nearT = -1; uint64_t seed = 1; const char *dumpfn = NULL;
    int nfirst = 0, firstvals[256], nonly = 0, onlyvals[3 * MAXK];
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "--sample")) { sample = atol(argv[++i]); seed = atoll(argv[++i]); }
        else if (!strcmp(argv[i], "--maxout")) maxout = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--nearmiss")) nearT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--need0")) need0 = 1;
        else if (!strcmp(argv[i], "--need2")) need2 = 1;
        else if (!strcmp(argv[i], "--dump")) dumpfn = argv[++i];
        else if (!strcmp(argv[i], "--first")) nfirst = parse_perm(argv[++i], firstvals);
        else if (!strcmp(argv[i], "--only")) nonly = parse_perm(argv[++i], onlyvals);
        else if (!strcmp(argv[i], "--chart")) load_chart(argv[++i]);
        else if (!strcmp(argv[i], "--bad")) {
            FILE *b = fopen(argv[++i], "r"); if (!b) { perror("bad"); return 2; }
            int cnt;
            while (fscanf(b, "%d", &cnt) == 1) {
                memset(badmask[nbad], 0, sizeof badmask[nbad]);
                for (int j = 0; j < cnt; j++) { int w; if (fscanf(b, "%d", &w) != 1) return 2; badmask[nbad][w >> 6] |= 1ULL << (w & 63); }
                nbad++;
            }
            fclose(b);
        }
    }
    rng_state ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    build_words();
    // spanning tree over darts a, A, b, B, s (links only); tree darts get the identity voltage
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
    // free slots: 0..k-1 a-slots, k..2k-1 b-slots, 2k.. s-slots (link: c < s0[c]; semi-edge: c == s0[c], involution domain)
    int nfree = 0, freeslot[3 * MAXK], domain[3 * MAXK];
    for (int c = 0; c < k; c++) if (!intree_a[c]) { freeslot[nfree] = c; domain[nfree++] = nG; }
    for (int c = 0; c < k; c++) if (!intree_b[c]) { freeslot[nfree] = k + c; domain[nfree++] = nG; }
    if (K == 5) for (int c = 0; c < k; c++) {
        if (s0[c] == c) { if (ninvols == 0) { fprintf(stderr, "semi-edge but G has no involution\n"); return 2; } freeslot[nfree] = 2 * k + c; domain[nfree++] = ninvols; }
        else if (c < s0[c] && !intree_s[c]) { freeslot[nfree] = 2 * k + c; domain[nfree++] = nG; }
    }
    memset(va, 0, sizeof va); memset(vb, 0, sizeof vb); memset(vs, 0, sizeof vs);
    double total = 1; for (int i = 0; i < nfree; i++) total *= (i == 0 && nfirst) ? nfirst : domain[i];
    fprintf(stderr, "k=%d K=%d nG=%d N=%d involutions=%d free slots=%d (%.3g assignments) badsets=%d%s%s\n", k, K, nG, N, ninvols, nfree, total, nbad, sample ? " sampling" : "", exact ? " exact" : "");
    FILE *dump = dumpfn ? fopen(dumpfn, "a") : NULL;
    long tried = 0, valid = 0, anyok = 0, hits = 0, near = 0, exact_hits = 0, exact_tried = 0;
    int best = 1 << 30; long best_exact = -1;
    int idx[3 * MAXK]; memset(idx, 0, sizeof idx);
    long iter = 0;
    if (nonly && nonly != nfree) { fprintf(stderr, "--only needs %d values\n", nfree); return 2; }
    while (1) {
        if (nonly) {
            if (iter++ > 0) break;
            for (int i = 0; i < nfree; i++) idx[i] = onlyvals[i];
        } else if (sample) {
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
            int s = freeslot[i]; int val = (i == 0 && nfirst && !nonly) ? firstvals[idx[i]] : idx[i];
            if (s < k) va[s] = val; else if (s < 2 * k) vb[s - k] = val;
            else { int c = s - 2 * k; vs[c] = (s0[c] == c && !nonly) ? invols[val] : val; }
        }
        tried++;
        if (!build_and_check()) continue;
        valid++;
        int u_any = uncovered(-1);
        if (u_any) { if (u_any < best && nbad == 0) best = u_any; continue; }
        anyok++;
        if (dump) {
            fprintf(dump, "base=%s/%s/%s ", argv[2], argv[3], argv[4]);
            print_volt(dump); fprintf(dump, " |");
            for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) if (needed(dist[c][v]) || dist[c][v] == 0 || dist[c][v] == 2) {
                int w2 = -1;
                if (dist[c][v] == 2) {   // index of a geodesic 2-word (t, u), u != t^-1, lexicographic among the K(K-1) reduced 2-words
                    int idx = 0;
                    for (int t = 0; t < K && w2 < 0; t++) for (int u = 0; u < K; u++) { if (u == symInv[t]) continue; if (nbr[nbr[c * nG][t]][u] == v) { w2 = idx; break; } idx++; }
                }
                fprintf(dump, " %d %d %d %d", c, v, dist[c][v], w2);
                for (int i = 0; i < 5; i++) fprintf(dump, " %llx", (unsigned long long)wmask[c][v][i]);
            }
            fprintf(dump, "\n"); fflush(dump);
        }
        int bestj = -1, bestu = 0;
        if (nbad) {
            bestu = 1 << 30;
            for (int j = 0; j < nbad; j++) { int u = uncovered(j); if (u < bestu) { bestu = u; bestj = j; } }
            if (bestu < best) best = bestu;
        }
        if (bestu == 0) {
            hits++;
            if (exact) {
                exact_tried++;
                long eb = exact_bad(0);
                if (best_exact < 0 || eb < best_exact) best_exact = eb;
                if (eb == 0) {
                    exact_hits++;
                    printf("EXACTHIT chart=%d ", bestj); print_volt(stdout); printf("\n"); fflush(stdout);
                } else if (hits <= maxout) {
                    printf("WORDHIT exactbad=%ld chart=%d ", eb, bestj); print_volt(stdout); printf("\n"); fflush(stdout);
                }
            } else if (hits <= maxout) {
                printf("HIT chart=%d ", bestj); print_volt(stdout); printf("\n"); fflush(stdout);
            }
        } else if (nearT >= 0 && bestu <= nearT) {
            near++;
            if (near <= maxout) { printf("NEAR uncovered=%d chart=%d ", bestu, bestj); print_volt(stdout); printf("\n"); fflush(stdout); }
        }
    }
    if (dump) fclose(dump);
    printf("SUMMARY k=%d K=%d nG=%d N=%d tried=%ld valid(simple,diam<=4)=%ld anyword_ok=%ld wordhits=%ld exact_tried=%ld exact_hits=%ld near=%ld best_uncovered=%d best_exact=%ld\n",
           k, K, nG, N, tried, valid, anyok, hits, exact_tried, exact_hits, near, best == (1 << 30) ? -1 : best, best_exact);
    return 0;
}
