// mixctrl: controller screen for route-chart lifts with MIXED alphabets over F_3^4
//   four line symbols a, a^-1, b, b^-1 (3 neighbours each) + P plain symbols (1 neighbour each), P = 1..4: degree 12 + P;
//   with --s also a fifth, self-inverse LINE symbol s (3 neighbours; the (15,4) alphabet): degree 15 + P.  --s needs k = 1:
//   the s voltage runs over the involutions of G (a third free slot), s non-backtracking means no "s s".
// Controllers: regular lifts over a group G (table) of a k-vertex Schreier base (alpha0, beta0) with voltages va[c], vb[c]
// (spanning-tree normalisation, exhaustive or sampled); the plain symbols are CENTRAL elements of G acting by right
// multiplication at every vertex:
//   inv z        a self-inverse plain symbol s: (c, g) -> (c, g z), z a central involution      (degree +1)
//   pair z       a plain pair p, p^-1: (c, g) -> (c, g z^{+-1}), z central of order >= 3        (degree +2)
// --plain takes a list of items: "inv 13" (degree 13), "pair 21" or the older spelling "z3 21" (14), "inv 13 13" (two
// involutions, 14), "pair 5 inv 13" (15), "pair 3 pair 9" (16).  --specs inv|pair|pairinv|pairpair enumerates ALL
// such specs of the group (central elements up to inversion, unordered pairs of classes) inside the (a, b) loop.
// Chart-independent coverage options of an ordered pair (root c, vertex v) with line distance dL (line symbols only):
//   dL = 0: own fibre (C0, chart only); dL = 2: C2 pencil of the geodesic 2-word (chart only);
//   otherwise NEEDED: a non-backtracking line 4-word reaching v (controllable -> 81 points) or a TILING: a line 3-word W
//   with W(c) z = v for a plain shift z (central, so the 4 insertions of the plain letter into W all reach W(c) z and,
//   for scalar plain charts M = +-I, are cosets of one hyperplane V_W: a tiling when the residues cover F_3).
// A needed pair with neither option is HARD (not coverable by 27-point pieces for any chart) unless v
// is a TWIN of a line neighbour, v = t(c) z, which the 32 backtracking one-plain words reach with 2-flats through d_t
// (a chart-only condition: SOFT, counted separately).  Plain-plain targets c z1 z2 are hard (only 3-point pieces).
// hard deficit 0 = admissible.  --dump writes one line per admissible controller:
//   "base=A0/B0 va=.. vb=.. plain=<spec> group=<file> | c v dL w2 m4hi m4lo m3[shift 0] m3[shift 1] ..."
//   (m4: 108-bit NB-4 word mask, m3: 36-bit NB-3 word masks of the tilings per distinct plain shift, w2 = geodesic
//   2-word index or -1).  k = 1 (Cayley) runs enumerate (a, b) up to a <-> a^-1, b <-> b^-1, a <-> b.
// Usage: mixctrl k alpha0 beta0 group.tbl (--plain ITEMS | --specs KIND) [--sample S seed] [--dump FILE] [--maxout M]
//        [--first v1,v2,..] [--stats] [--need2] [--show] [--nodiam] [--all-ab]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 1024
#define MAXK 16
#define MAXL 5
#define MAXP 4
#define MAXSYM (MAXL + MAXP)
#define NW4MAX 320
#define NW3MAX 80
#define M4W 5
#define M3W 2
#define MAXSPEC 512
static int KL = 4, NW4 = 108, NW3 = 36, use_s = 0, vs = 0;   // line symbols: 4, or 5 with --s

static int k, nG, N;
static int a0[MAXK], b0[MAXK], a0i[MAXK], b0i[MAXK];
static int *T, *inv;
static int va[MAXK], vb[MAXK];
static int nbr[MAXN][MAXSYM];          // 0 a, 1 a^-1, 2 b, 3 b^-1, 4.. plain symbols
static int symInv[MAXL] = {1, 0, 3, 2, 4};
static int w4[NW4MAX][4], w3[NW3MAX][3];
static int dL[MAXK][MAXN], dC[MAXK][MAXN];
static uint64_t m4[MAXK][MAXN][M4W];   // NB-4 word masks per (root, v) (108 bits, 320 with --s)
static uint64_t m3[MAXK][MAXN][M3W];   // NB-3 word masks per (root, v) (endpoint of the line 3-walk; 36 / 80 bits)
static inline int m4any(int c, int v) { int r = 0; for (int w = 0; w < M4W; w++) r |= (m4[c][v][w] != 0); return r; }
static inline int m4pop(int c, int v) { int r = 0; for (int w = 0; w < M4W; w++) r += __builtin_popcountll(m4[c][v][w]); return r; }
static unsigned char twin[MAXK][MAXN]; // v = z(t(c)) for a line symbol t and a plain shift z
static int need2 = 0, nodiam = 0;

// a plain spec: P symbols, each a central element; pairs occupy two consecutive symbols (z, z^-1)
typedef struct { int P; int sym[MAXP]; int nshift; int shift[MAXP]; char name[64]; } Spec;
static Spec specs[MAXSPEC]; static int nspec = 0;
static int *cen, ncen;                 // central elements != e
static int gorder(int z) { int o = 1, x = z; while (x != 0) { x = T[x * nG + z]; o++; } return o; }
static int is_central(int z) { for (int g = 0; g < nG; g++) if (T[g * nG + z] != T[z * nG + g]) return 0; return 1; }
static void spec_finish(Spec *s) {
    s->nshift = 0;
    for (int j = 0; j < s->P; j++) { int z = s->sym[j], seen = 0; for (int i = 0; i < s->nshift; i++) seen |= (s->shift[i] == z); if (!seen) s->shift[s->nshift++] = z; }
}
static int spec_add_item(Spec *s, const char *kind, int z, char *name, size_t nsz) {
    if (z <= 0 || z >= nG) { fprintf(stderr, "plain element %d out of range\n", z); return 0; }
    if (!is_central(z)) { fprintf(stderr, "plain element %d not central\n", z); return 0; }
    int o = gorder(z);
    if (!strcmp(kind, "inv")) {
        if (o != 2) { fprintf(stderr, "plain element %d not an involution\n", z); return 0; }
        if (s->P + 1 > MAXP) { fprintf(stderr, "too many plain symbols\n"); return 0; }
        s->sym[s->P++] = z;
        snprintf(name + strlen(name), nsz - strlen(name), "%sinv:%d", strlen(name) ? "," : "", z);
    } else {
        if (o < 3) { fprintf(stderr, "plain element %d is an involution: use inv\n", z); return 0; }
        if (s->P + 2 > MAXP) { fprintf(stderr, "too many plain symbols\n"); return 0; }
        s->sym[s->P++] = z; s->sym[s->P++] = inv[z];
        snprintf(name + strlen(name), nsz - strlen(name), "%spair:%d", strlen(name) ? "," : "", z);
    }
    return 1;
}
// a key that is invariant under automorphisms for one central element: order + numbers of 2nd/3rd/4th/6th roots
static int *rootcnt[7];
static void root_counts(void) {
    int ks[4] = {2, 3, 4, 6};
    for (int i = 0; i < 4; i++) { int kk = ks[i]; rootcnt[kk] = calloc(nG, sizeof(int));
        for (int g = 0; g < nG; g++) { int x = 0; for (int e = 0; e < kk; e++) x = T[x * nG + g]; rootcnt[kk][x]++; } }
}
static long ckey(int z) { return ((((long)gorder(z) * 1000 + rootcnt[2][z]) * 1000 + rootcnt[3][z]) * 1000 + rootcnt[4][z]) * 1000 + rootcnt[6][z]; }
// --specs KIND: every spec of the kind up to the invariant key of its elements (single element: ckey; two elements:
// the two ckeys with the orders of z w and z w^-1), one representative per key; pairs use z < z^-1.
static long speckeys[MAXSPEC];
static void add_spec(const char *k1, int z, const char *k2, int w, long key) {
    for (int i = 0; i < nspec; i++) if (speckeys[i] == key) return;
    if (nspec >= MAXSPEC) return;
    Spec *s = &specs[nspec]; memset(s, 0, sizeof *s);
    if (!spec_add_item(s, k1, z, s->name, sizeof s->name)) exit(2);
    if (k2 && !spec_add_item(s, k2, w, s->name, sizeof s->name)) exit(2);
    spec_finish(s); speckeys[nspec++] = key;
}
static long pairkey(int z, int w) {
    long kz = ckey(z), kw = ckey(w); long lo = kz < kw ? kz : kw, hi = kz < kw ? kw : kz;
    int o1 = gorder(T[z * nG + w]), o2 = gorder(T[z * nG + inv[w]]); int olo = o1 < o2 ? o1 : o2, ohi = o1 < o2 ? o2 : o1;
    return ((lo * 1000003L + hi) * 1000 + olo) * 1000 + ohi;
}
static void enumerate_specs(const char *kind) {
    root_counts();
    int one_inv = !strcmp(kind, "inv"), one_pair = !strcmp(kind, "pair"), pi = !strcmp(kind, "pairinv"), pp = !strcmp(kind, "pairpair");
    if (!one_inv && !one_pair && !pi && !pp) { fprintf(stderr, "--specs kind must be inv|pair|pairinv|pairpair\n"); exit(2); }
    for (int ci = 0; ci < ncen; ci++) {
        int z = cen[ci], oz = gorder(z);
        if (one_inv && oz == 2) add_spec("inv", z, NULL, 0, ckey(z));
        if (oz < 3 || z > inv[z]) continue;
        if (one_pair) add_spec("pair", z, NULL, 0, ckey(z));
        for (int cj = 0; cj < ncen && (pi || pp); cj++) {
            int w = cen[cj], ow = gorder(w);
            if (pi && ow == 2) add_spec("pair", z, "inv", w, pairkey(z, w));
            if (pp && ow >= 3 && w < inv[w] && w > z) add_spec("pair", z, "pair", w, pairkey(z, w));
        }
    }
}

static void build_words(void) {
    int n = 0;
    for (int x0 = 0; x0 < KL; x0++) for (int x1 = 0; x1 < KL; x1++) for (int x2 = 0; x2 < KL; x2++) for (int x3 = 0; x3 < KL; x3++) {
        if (x1 == symInv[x0] || x2 == symInv[x1] || x3 == symInv[x2]) continue;
        w4[n][0] = x0; w4[n][1] = x1; w4[n][2] = x2; w4[n][3] = x3; n++;
    }
    NW4 = n;
    n = 0;
    for (int x0 = 0; x0 < KL; x0++) for (int x1 = 0; x1 < KL; x1++) for (int x2 = 0; x2 < KL; x2++) {
        if (x1 == symInv[x0] || x2 == symInv[x1]) continue;
        w3[n][0] = x0; w3[n][1] = x1; w3[n][2] = x2; n++;
    }
    NW3 = n;
    if (NW4 != (use_s ? 320 : 108) || NW3 != (use_s ? 80 : 36)) { fprintf(stderr, "word counts %d %d\n", NW4, NW3); exit(1); }
}

// line part of a controller: simplicity of the line edges, line distances, word masks (independent of the plain spec)
static int build_line(void) {
    static int alpha[MAXN], beta[MAXN], alphai[MAXN], betai[MAXN];
    for (int c = 0; c < k; c++) for (int g = 0; g < nG; g++) {
        int v = c * nG + g;
        alpha[v] = a0[c] * nG + T[g * nG + va[c]];
        beta[v] = b0[c] * nG + T[g * nG + vb[c]];
    }
    for (int v = 0; v < N; v++) { alphai[alpha[v]] = v; betai[beta[v]] = v; }
    for (int v = 0; v < N; v++) {
        nbr[v][0] = alpha[v]; nbr[v][1] = alphai[v]; nbr[v][2] = beta[v]; nbr[v][3] = betai[v];
        // s: k = 1 semi-edge with an involution voltage; k = 2 link 0 <-> 1 with voltage vs on the dart 0 -> 1 (vs^-1 back)
        if (use_s) nbr[v][4] = (k == 1) ? T[v * nG + vs] : (v < nG ? nG + T[v * nG + vs] : T[(v - nG) * nG + inv[vs]]);
    }
    for (int c = 0; c < k; c++) {
        int v = c * nG;
        for (int i = 0; i < KL; i++) { if (nbr[v][i] == v) return 0; for (int j = i + 1; j < KL; j++) if (nbr[v][i] == nbr[v][j]) return 0; }
    }
    static int queue[MAXN];
    for (int c = 0; c < k; c++) {
        int *e = dL[c]; for (int v = 0; v < N; v++) e[v] = -1;
        int s = c * nG, qh = 0, qt = 0; e[s] = 0; queue[qt++] = s;
        while (qh < qt) { int x = queue[qh++]; for (int i = 0; i < KL; i++) { int y = nbr[x][i]; if (e[y] < 0) { e[y] = e[x] + 1; queue[qt++] = y; } } }
        memset(m4[c], 0, N * sizeof(m4[c][0])); memset(m3[c], 0, N * sizeof(m3[c][0]));
        for (int w = 0; w < NW4; w++) { int x = s; for (int i = 0; i < 4; i++) x = nbr[x][w4[w][i]]; m4[c][x][w >> 6] |= 1ULL << (w & 63); }
        for (int w = 0; w < NW3; w++) { int x = s; for (int i = 0; i < 3; i++) x = nbr[x][w3[w][i]]; m3[c][x][w >> 6] |= 1ULL << (w & 63); }
    }
    return 1;
}
// plain part: neighbours 4.., simplicity against the line edges, controller distance (diameter <= 4), twins
static const Spec *cur;
static int build_plain(const Spec *s) {
    cur = s;
    for (int c = 0; c < k; c++) for (int g = 0; g < nG; g++) { int v = c * nG + g; for (int j = 0; j < s->P; j++) nbr[v][KL + j] = c * nG + T[g * nG + s->sym[j]]; }
    // simplicity at the roots (G-invariant): plain neighbours distinct from the root and the line neighbours; two plain
    // symbols on the same element only for two involutions (double edge with two charts)
    for (int c = 0; c < k; c++) {
        int v = c * nG;
        for (int j = 0; j < s->P; j++) { int y = nbr[v][KL + j]; if (y == v) return 0; for (int i = 0; i < KL; i++) if (nbr[v][i] == y) return 0; }
        for (int j = 0; j < s->P; j++) for (int i = j + 1; i < s->P; i++) if (nbr[v][KL + j] == nbr[v][KL + i] && gorder(s->sym[j]) != 2) return 0;
    }
    static int queue[MAXN];
    int K = KL + s->P;
    for (int c = 0; c < k; c++) {
        int *d = dC[c]; for (int v = 0; v < N; v++) d[v] = -1;
        int src = c * nG, qh = 0, qt = 0; d[src] = 0; queue[qt++] = src;
        while (qh < qt) { int x = queue[qh++]; if (d[x] >= 4) continue; for (int i = 0; i < K; i++) { int y = nbr[x][i]; if (d[y] < 0) { d[y] = d[x] + 1; queue[qt++] = y; } } }
        if (qt < N && !nodiam) return 0;
        memset(twin[c], 0, N);
        for (int t = 0; t < KL; t++) for (int j = 0; j < s->P; j++) twin[c][nbr[nbr[src][t]][KL + j]] = 1;
    }
    return 1;
}
// tiling mask of (c, v) for shift index i: 3-words W with W(c) z = v  <=>  W(c) = v z^-1
static inline uint64_t tmask(int c, int v, int i, int w) {
    int cc = v / nG, g = v % nG;
    return m3[c][cc * nG + T[g * nG + inv[cur->shift[i]]]][w];
}
static inline int tmask_any(int c, int v) { for (int i = 0; i < cur->nshift; i++) for (int w = 0; w < M3W; w++) if (tmask(c, v, i, w)) return 1; return 0; }
static inline int tmask_pop(int c, int v) { int r = 0; for (int w = 0; w < M3W; w++) { uint64_t t = 0; for (int i = 0; i < cur->nshift; i++) t |= tmask(c, v, i, w); r += __builtin_popcountll(t); } return r; }
static int needed(int d) { return d != 0 && (d != 2 || need2); }

static int nstat_needed, nstat_tileonly, nstat_wordonly, nstat_single, nstat_soft;
static int deficit(void) {
    int cnt = 0; nstat_needed = nstat_tileonly = nstat_wordonly = nstat_single = nstat_soft = 0;
    for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) {
        if (!needed(dL[c][v])) continue;
        nstat_needed++;
        int a = m4any(c, v), t = tmask_any(c, v);
        if (!a && !t) { if (twin[c][v]) nstat_soft++; else cnt++; }
        else {
            if (!a) nstat_tileonly++;
            if (!t) nstat_wordonly++;
            if (m4pop(c, v) + tmask_pop(c, v) == 1) nstat_single++;
        }
    }
    return cnt;
}

static int parse_perm(const char *s, int *out) {
    int n = 0; const char *p = s;
    while (*p) { out[n++] = atoi(p); while (*p && *p != ',') p++; if (*p == ',') p++; }
    return n;
}
static int is_int(const char *s) { if (*s == '-' ) s++; if (!*s) return 0; while (*s) { if (*s < '0' || *s > '9') return 0; s++; } return 1; }
static uint64_t rng_state = 88172645463325252ULL;
static uint64_t xr(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static void print_volts(FILE *f) {
    fprintf(f, "va="); for (int c = 0; c < k; c++) fprintf(f, "%d%s", va[c], c + 1 < k ? "," : "");
    fprintf(f, " vb="); for (int c = 0; c < k; c++) fprintf(f, "%d%s", vb[c], c + 1 < k ? "," : "");
    if (use_s) fprintf(f, " vs=%d", vs);
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: mixctrl k alpha0 beta0 group.tbl (--plain ITEMS | --specs KIND) [options]\n"); return 2; }
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
    cen = malloc(nG * sizeof(int)); ncen = 0; for (int z = 1; z < nG; z++) if (is_central(z)) cen[ncen++] = z;
    long sample = 0; int maxout = 20; uint64_t seed = 1; const char *dumpfn = NULL; int stats = 0, show = 0, allab = 0;
    int nfirst = 0, firstvals[256];
    const char *speckind = NULL;
    for (int i = 5; i < argc; i++) {
        if (!strcmp(argv[i], "--plain")) {
            Spec *s = &specs[0]; memset(s, 0, sizeof *s);
            while (i + 1 < argc && strncmp(argv[i + 1], "--", 2)) {
                const char *kind = argv[++i];
                if (strcmp(kind, "inv") && strcmp(kind, "pair") && strcmp(kind, "z3")) { fprintf(stderr, "plain item kind %s\n", kind); return 2; }
                if (!strcmp(kind, "z3")) kind = "pair";
                if (i + 1 >= argc || !is_int(argv[i + 1])) { fprintf(stderr, "plain item needs an element\n"); return 2; }
                if (!spec_add_item(s, kind, atoi(argv[++i]), s->name, sizeof s->name)) return 2;
                if (!strcmp(kind, "inv") && i + 1 < argc && is_int(argv[i + 1])) { if (!spec_add_item(s, "inv", atoi(argv[++i]), s->name, sizeof s->name)) return 2; }
            }
            if (s->P == 0) { fprintf(stderr, "--plain needs items\n"); return 2; }
            spec_finish(s); nspec = 1;
        }
        else if (!strcmp(argv[i], "--specs")) speckind = argv[++i];
        else if (!strcmp(argv[i], "--sample")) { sample = atol(argv[++i]); seed = atoll(argv[++i]); }
        else if (!strcmp(argv[i], "--maxout")) maxout = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--dump")) dumpfn = argv[++i];
        else if (!strcmp(argv[i], "--first")) nfirst = parse_perm(argv[++i], firstvals);
        else if (!strcmp(argv[i], "--stats")) stats = 1;
        else if (!strcmp(argv[i], "--need2")) need2 = 1;
        else if (!strcmp(argv[i], "--show")) show = 1;
        else if (!strcmp(argv[i], "--nodiam")) nodiam = 1;
        else if (!strcmp(argv[i], "--all-ab")) allab = 1;
        else if (!strcmp(argv[i], "--s")) { use_s = 1; KL = 5; }
        else if (!strcmp(argv[i], "--soft2")) { /* older spelling: hard / soft is the only mode now */ }
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    if (use_s && k > 2) { fprintf(stderr, "--s needs k = 1 (semi-edge) or k = 2 (link)\n"); return 2; }
    if (speckind) { nspec = 0; enumerate_specs(speckind); }
    if (!nspec) { fprintf(stderr, "no plain spec (--plain or --specs); central elements: %d\n", ncen); printf("SUMMARY k=%d nG=%d N=%d nspecs=0\n", k, nG, N); return 0; }
    rng_state ^= seed * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    build_words();
    // spanning tree over the darts a, A, b, B; tree darts get the identity voltage
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
    int nfree = 0, freeslot[2 * MAXK];
    for (int c = 0; c < k; c++) if (!intree_a[c]) freeslot[nfree++] = c;
    for (int c = 0; c < k; c++) if (!intree_b[c]) freeslot[nfree++] = k + c;
    memset(va, 0, sizeof va); memset(vb, 0, sizeof vb);
    double total = 1; for (int i = 0; i < nfree; i++) total *= (i == 0 && nfirst) ? nfirst : nG;
    int reduce = (k == 1 && !allab && !sample && !nfirst);
    int invs[MAXN], ninv = 0; for (int g = 1; g < nG; g++) if (inv[g] == g) invs[ninv++] = g;
    int nsv = (k == 1) ? ninv : nG;   // s voltages: involutions (semi-edge) or every element (link)
    if (use_s) { total *= nsv; if (!nsv) { fprintf(stderr, "no involutions for s\n"); return 2; } }
    fprintf(stderr, "k=%d nG=%d N=%d central=%d specs=%d (P=%d..%d) free slots=%d (%.3g assignments)%s%s\n", k, nG, N, ncen, nspec, specs[0].P, specs[nspec - 1].P, nfree, total, sample ? " sampling" : "", reduce ? " a/b reduced" : "");
    FILE *dump = dumpfn ? fopen(dumpfn, "a") : NULL;
    long tried = 0, linevalid = 0;
    static long valid[MAXSPEC], ok[MAXSPEC], hist[MAXSPEC][16]; static int best[MAXSPEC], bestva[MAXSPEC][MAXK], bestvb[MAXSPEC][MAXK], bestsoft[MAXSPEC], bestvs[MAXSPEC];
    for (int s = 0; s < nspec; s++) { best[s] = 1 << 30; valid[s] = ok[s] = 0; }
    int idx[2 * MAXK]; memset(idx, 0, sizeof idx);
    long iter = 0;
    while (1) {
        if (sample) {
            if (iter++ >= sample) break;
            for (int i = 0; i < nfree; i++) idx[i] = (i == 0 && nfirst) ? (int)(xr() % nfirst) : (int)(xr() % nG);
        } else {
            if (iter++ > 0) {
                int i = 0;
                while (i < nfree) { int dom = (i == 0 && nfirst) ? nfirst : nG; if (++idx[i] < dom) break; idx[i] = 0; i++; }
                if (i == nfree) break;
            }
        }
        for (int i = 0; i < nfree; i++) {
            int s = freeslot[i]; int val = (i == 0 && nfirst) ? firstvals[idx[i]] : idx[i];
            if (s < k) va[s] = val; else vb[s - k] = val;
        }
        if (reduce && (inv[va[0]] < va[0] || inv[vb[0]] < vb[0] || vb[0] < va[0])) continue;
        for (int ii = 0; ii < (use_s ? (sample ? 1 : nsv) : 1); ii++) {
        vs = !use_s ? 0 : (k == 1 ? invs[sample ? (int)(xr() % ninv) : ii] : (sample ? (int)(xr() % nG) : ii));
        tried++;
        if (!build_line()) continue;
        linevalid++;
        for (int si = 0; si < nspec; si++) {
            Spec *s = &specs[si];
            if (!build_plain(s)) continue;
            valid[si]++;
            int def = deficit();
            if (def < best[si] || (def == best[si] && nstat_soft < bestsoft[si])) { best[si] = def; bestsoft[si] = nstat_soft; memcpy(bestva[si], va, sizeof va); memcpy(bestvb[si], vb, sizeof vb); bestvs[si] = vs; }
            hist[si][def < 15 ? def : 15]++;
            if (def) continue;
            ok[si]++;
            if (ok[si] <= maxout) {
                printf("OK "); print_volts(stdout);
                printf(" plain=%s needed=%d tileonly=%d wordonly=%d singletons=%d soft=%d\n", s->name, nstat_needed, nstat_tileonly, nstat_wordonly, nstat_single, nstat_soft);
                fflush(stdout);
            }
            if (dump) {
                fprintf(dump, "base=%s/%s ", argv[2], argv[3]); print_volts(dump);
                fprintf(dump, " plain=%s group=%s soft=%d |", s->name, argv[4], nstat_soft);
                for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) {
                    int d = dL[c][v];
                    if (d == 0) continue;
                    int w2 = -1;
                    if (d == 2) {
                        int ix = 0;
                        for (int t = 0; t < KL && w2 < 0; t++) for (int u = 0; u < KL; u++) { if (u == symInv[t]) continue; if (nbr[nbr[c * nG][t]][u] == v) { w2 = ix; break; } ix++; }
                    }
                    fprintf(dump, " %d %d %d %d", c, v, d, w2);
                    for (int w = (use_s ? M4W : 2) - 1; w >= 0; w--) fprintf(dump, " %llx", (unsigned long long)m4[c][v][w]);
                    for (int i = 0; i < s->nshift; i++) for (int w = (use_s ? M3W : 1) - 1; w >= 0; w--) fprintf(dump, " %llx", (unsigned long long)tmask(c, v, i, w));
                }
                fprintf(dump, "\n"); fflush(dump);
            }
        }
        }
    }
    if (dump) fclose(dump);
    for (int si = 0; si < nspec; si++) {
        Spec *s = &specs[si];
        if (show && best[si] < (1 << 30)) {
            memcpy(va, bestva[si], sizeof va); memcpy(vb, bestvb[si], sizeof vb); vs = bestvs[si]; build_line(); build_plain(s); deficit();
            printf("BEST plain=%s ", s->name); print_volts(stdout);
            printf(" deficit=%d soft=%d\n", best[si], nstat_soft);
            int tot[8] = {0}, unc[8] = {0}, wo[8] = {0}, to[8] = {0}, bo[8] = {0};
            for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) {
                int d = dL[c][v]; int di = d < 0 ? 7 : (d > 6 ? 6 : d);
                tot[di]++;
                if (!needed(d)) continue;
                int a = m4any(c, v), t = tmask_any(c, v);
                if (!a && !t) unc[di]++; else if (a && !t) wo[di]++; else if (!a && t) to[di]++; else bo[di]++;
            }
            for (int c = 0; c < k; c++) for (int v = 0; v < N; v++) {
                int d = dL[c][v]; if (!needed(d)) continue;
                int a = m4any(c, v), t = tmask_any(c, v);
                if (!a && !t) { int w3 = 0; for (int w = 0; w < M3W; w++) w3 += __builtin_popcountll(m3[c][v][w]);
                    printf("  %s c=%d v=%d dL=%d dC=%d walks3=%d\n", twin[c][v] ? "SOFT" : "DEFICIT", c, v, d, dC[c][v], w3); }
            }
            for (int d = 0; d < 8; d++) if (tot[d]) printf("  dL=%s%d: targets %d uncovered %d word-only %d tiling-only %d both %d\n", d == 7 ? "inf/" : "", d == 7 ? -1 : d, tot[d], unc[d], wo[d], to[d], bo[d]);
        }
        printf("SUMMARY k=%d nG=%d N=%d n=%d deg=%d plain=%s tried=%ld linevalid=%ld valid(simple,diam<=4)=%ld mixedOK=%ld best_deficit=%d best_soft=%d", k, nG, N, 81 * N, 3 * KL + s->P, s->name, tried, linevalid, valid[si], ok[si], best[si] == (1 << 30) ? -1 : best[si], best[si] == (1 << 30) ? -1 : bestsoft[si]);
        if (stats) { printf(" hist="); for (int i = 0; i < 16; i++) printf("%ld%s", hist[si][i], i < 15 ? "," : "+"); }
        printf("\n");
    }
    return 0;
}
