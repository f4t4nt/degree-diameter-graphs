/*
 * sweep_odd.c -- COMPLETE enumeration of the odd-degree Cayley graphs Cay(G, {t} u {a_1^{+-1}, ..., a_K^{+-1}}) on the
 * metacyclic group G = Z_N x|_s Z_M (N odd, M even, s^M = 1 mod N, any kernel), t = y^(M/2), degree 2K + 1.
 *
 * Conventions: element x^i y^j <-> index i + N j, product
 * (i1, j1)(i2, j2) = (i1 + s^j1 i2 mod N, j1 + j2 mod M), edges g ~ g a (right multiplication), BFS from the identity
 * (vertex-transitive: its eccentricity is the diameter).  defect = number of vertices at distance > D.
 *
 * Symmetry: A_t = { phi_{i,j,l} : x -> x^i, y -> x^j y^l | i unit mod N, l unit mod M with s^l = s mod N,
 * j s_{M/2} = 0 mod N (s_b = 1 + s + ... + s^(b-1)) } is a group of automorphisms of G fixing t (checked explicitly on
 * small groups with --check).  Connection sets are K-subsets of the "classes" {g, g^-1} (g != e, t, not an involution),
 * indexed in element order; a set is emitted only if it can be the lexicographically smallest member of its A_t-orbit:
 * with S = {c_0 < c_1 < ... < c_{K-1}}, stab_m = pointwise stabiliser of c_0..c_{m-1} in A_t (stab_0 = A_t) and
 * omin_m(c) = min of the stab_m-orbit of c, every c_k must satisfy omin_m(c_k) >= c_m for all m <= k (equality forced
 * at m = k).  These are necessary conditions for lexicographic minimality, so every orbit is represented (--check
 * verifies this by brute force: every A_t-orbit of K-subsets contains an emitted set).
 *
 * usage: sweep_odd N M s K D [--bound b] [--sample R seed] [--brute] [--check] [--list L] [--quiet]
 *                        [--c0 lo hi] [--c1 lo hi] [--count] [--c0screen R seed]
 *   default: the reduced sweep, tracking the exact minimum defect (bound = current best) and listing every set
 *            with defect <= L (default 0) as FOUND/NEAR lines "x1 y1 x2 y2 ..." (slot order; t implicit).
 *   --bound b   stop the BFS as soon as more than b vertices are unreached (b >= 0; default: adaptive, exact minimum)
 *   --sample R seed   R uniformly random K-subsets of classes, full BFS each: mean / quartiles / min of the defect
 *   --brute     evaluate ALL K-subsets without reduction (small groups only)
 *   --check     validate A_t (automorphism property on all element pairs) and the orbit coverage of the reduced
 *               enumeration against brute force (small groups only), then compare the minima of both enumerations
 *
 * Chunking: the reduced enumeration is a tree (level k chooses c_k); the chunk options restrict
 *   --c0 lo hi  the first slot to class indices c_0 in [lo, hi)      (default [0, classes))
 *   --c1 lo hi  the second slot to class indices c_1 in [lo, hi)     (default [0, classes))
 *   so the union of the chunks over a partition of the ranges is exactly the full reduced enumeration (each chunk
 *   applies the same minimality tests; nothing depends on the order or on other chunks).  After every completed
 *   c_0 subtree (or every completed c_1 subtree when the c_0 range is a single class) a line
 *   "PROGRESS level=<k> c=<c> candidates=<count> best=<b> time=<s>" is printed, so a run killed by timeout documents
 *   the exact completed prefix (resume with --c1 c+1 hi, resp. --c0 c+1 hi).
 *   --count     walk the reduced enumeration without any BFS (exact chunk sizes for planning)
 *   --c0screen R seed   for every A_t-orbit representative c_0: R random completions (K-1 further random classes),
 *               full BFS each; prints the defect distribution per c_0 (the "first-generator class" softness screen)
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef uint64_t u64;
typedef long long ll;

static int N, M, S, D, K, W, deg;
static ll n;
static u64 lastmask;
static int *rpow, *sb;
static int narcs, arc_x[16], arc_y[16], *arcT[16];
static u64 *vis, *cur, *nxt;
static uint8_t *rowflag_cur, *rowflag_nxt;
static unsigned *rowstamp, stampc = 0;
static long scalar_thresh;
static ll capacity[64];     /* capacity[r] = (deg-1) + (deg-1)^2 + ... + (deg-1)^r */
static long nbfs = 0;

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static int igcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

/* ---------- group, classes ---------- */
static int Cn, *cls_x, *cls_y, *elem2cls, ninvol;

static inline int inv_x(int x, int y) { int jj = (M - y) % M; return (int)((N - ((ll)rpow[jj] * x) % N) % N); }
static inline int inv_y(int y) { return (M - y) % M; }

static void build_classes(void) {
    elem2cls = malloc(n * sizeof(int));
    cls_x = malloc(n * sizeof(int)); cls_y = malloc(n * sizeof(int));
    Cn = 0; ninvol = 0;
    for (ll idx = 0; idx < n; idx++) elem2cls[idx] = -1;
    for (ll idx = 0; idx < n; idx++) {
        int x = (int)(idx % N), y = (int)(idx / N);
        if (idx == 0) continue;
        int ix = inv_x(x, y), iy = inv_y(y); ll idx2 = ix + (ll)N * iy;
        if (idx2 == idx) { ninvol++; continue; }           /* involution (t among them) */
        if (idx2 < idx) continue;                            /* already seen as the inverse */
        cls_x[Cn] = x; cls_y[Cn] = y; elem2cls[idx] = Cn; elem2cls[idx2] = Cn; Cn++;
    }
}

/* ---------- A_t ---------- */
static int nA, *A_i, *A_j, *A_l;

static void build_aut(void) {
    int h = M / 2;
    int cap = 1024; A_i = malloc(cap * sizeof(int)); A_j = malloc(cap * sizeof(int)); A_l = malloc(cap * sizeof(int)); nA = 0;
    for (int l = 1; l < M; l++) {
        if (igcd(l, M) != 1 || rpow[l] != rpow[1]) continue;
        for (int i = 1; i < N; i++) {
            if (igcd(i, N) != 1) continue;
            for (int j = 0; j < N; j++) {
                if (((ll)j * sb[h]) % N != 0) continue;
                if (nA == cap) { cap *= 2; A_i = realloc(A_i, cap * sizeof(int)); A_j = realloc(A_j, cap * sizeof(int)); A_l = realloc(A_l, cap * sizeof(int)); }
                A_i[nA] = i; A_j[nA] = j; A_l[nA] = l; nA++;
            }
        }
    }
}
/* image of the element (x, y) under phi_a: x^{i x + j s_y} y^{l y} */
static inline ll img_elem_g(int i, int j, int l, int x, int y) {
    ll xx = ((ll)i * x + (ll)j * sb[y]) % N;
    ll yy = ((ll)l * y) % M;
    return xx + (ll)N * yy;
}
static inline ll img_elem(int a, int x, int y) { return img_elem_g(A_i[a], A_j[a], A_l[a], x, y); }
static inline int img_cls(int a, int c) { return elem2cls[img_elem(a, cls_x[c], cls_y[c])]; }

/* generators of A_t: a greedy generating set of the i-maps ((Z/N)^*), of the l-maps (the subgroup {l : s^l = s} of
 * (Z/M)^*) and the translation j -> j0 generating J = {j : j s_{M/2} = 0}; verified by closure (composition
 * phi_2 o phi_1 = (i1 i2, i2 j1 + j2 s_{l1}, l1 l2)) to generate exactly |A_t| elements. */
static int ngen, G_i[64], G_j[64], G_l[64];
static void add_gen(int i, int j, int l) { G_i[ngen] = i; G_j[ngen] = j; G_l[ngen] = l; ngen++; }
static void greedy_gens(int mod, int (*member)(int), int which) {
    /* which = 0: i-maps (i, 0, 1); which = 1: l-maps (1, 0, l) */
    uint8_t *H = calloc(mod, 1); int *list = malloc(mod * sizeof(int)); int hs = 0, target = 0;
    for (int v = 1; v < mod; v++) if (member(v)) target++;
    H[1 % mod] = 1; list[hs++] = 1 % mod;
    for (int g = 2; g < mod && hs < target; g++) {
        if (!member(g) || H[g]) continue;
        if (which == 0) add_gen(g, 0, 1); else add_gen(1, 0, g);
        int old = hs;
        for (int t = 0; t < old; t++) { ll v = list[t]; for (;;) { v = (v * g) % mod; if (H[v]) break; H[v] = 1; list[hs++] = (int)v; } }
    }
    if (hs != target) { fprintf(stderr, "generator search failed (%d of %d)\n", hs, target); exit(3); }
    free(H); free(list);
}
static int member_i(int v) { return igcd(v, N) == 1; }
static int member_l(int v) { return igcd(v, M) == 1 && rpow[v] == rpow[1]; }
static void build_gens(void) {
    ngen = 0;
    greedy_gens(N, member_i, 0);
    greedy_gens(M, member_l, 1);
    int h = M / 2, g = sb[h] == 0 ? N : igcd(N, sb[h]);
    if (g > 1) add_gen(1, N / g, 1);
    /* closure check */
    size_t bits = (size_t)N * N * M;
    u64 *seen = calloc((bits + 63) / 64, 8);
    int *qi = malloc(nA * sizeof(int)), *qj = malloc(nA * sizeof(int)), *ql = malloc(nA * sizeof(int));
    long qh = 0, qt = 0;
    #define IDX(i, j, l) ((size_t)(i) + (size_t)N * ((size_t)(j) + (size_t)N * (size_t)(l)))
    #define SEEN(k) (seen[(k) >> 6] & (1ULL << ((k) & 63)))
    #define SETSEEN(k) (seen[(k) >> 6] |= (1ULL << ((k) & 63)))
    SETSEEN(IDX(1, 0, 1)); qi[qt] = 1; qj[qt] = 0; ql[qt] = 1; qt++;
    while (qh < qt) {
        int i1 = qi[qh], j1 = qj[qh], l1 = ql[qh]; qh++;
        for (int t = 0; t < ngen; t++) {
            int i2 = G_i[t], j2 = G_j[t], l2 = G_l[t];
            int i = (int)(((ll)i1 * i2) % N), j = (int)((((ll)i2 * j1) + (ll)j2 * sb[l1]) % N), l = (int)(((ll)l1 * l2) % M);
            size_t k = IDX(i, j, l);
            if (SEEN(k)) continue;
            SETSEEN(k);
            if (qt >= nA) { fprintf(stderr, "closure exceeds |A_t|: generators leave A_t?\n"); exit(3); }
            qi[qt] = i; qj[qt] = j; ql[qt] = l; qt++;
        }
    }
    if (qt != nA) { fprintf(stderr, "closure of the generators has %ld elements, |A_t| = %d\n", qt, nA); exit(3); }
    free(seen); free(qi); free(qj); free(ql);
    #undef IDX
    #undef SEEN
    #undef SETSEEN
}

/* orbit minima of the classes under A_t by union-find over the generators */
static int *uf;
static int uf_find(int a) { while (uf[a] != a) { uf[a] = uf[uf[a]]; a = uf[a]; } return a; }
static void compute_omin0(int *out) {
    uf = malloc(Cn * sizeof(int));
    for (int c = 0; c < Cn; c++) uf[c] = c;
    for (int t = 0; t < ngen; t++)
        for (int c = 0; c < Cn; c++) {
            int d = elem2cls[img_elem_g(G_i[t], G_j[t], G_l[t], cls_x[c], cls_y[c])];
            int ra = uf_find(c), rb = uf_find(d);
            if (ra != rb) { if (ra < rb) uf[rb] = ra; else uf[ra] = rb; }
        }
    int *first = malloc(Cn * sizeof(int));
    for (int c = 0; c < Cn; c++) first[c] = -1;
    for (int c = 0; c < Cn; c++) { int r = uf_find(c); if (first[r] < 0) first[r] = c; out[c] = first[r]; }
    free(first); free(uf);
}

static inline ll gmul(ll g, ll h) {
    ll x1 = g % N, y1 = g / N, x2 = h % N, y2 = h / N;
    return (x1 + ((ll)rpow[y1] * x2) % N) % N + (ll)N * ((y1 + y2) % M);
}

/* ---------- bit-parallel BFS (single fibre) ---------- */
static inline void rot_or(u64 *restrict dst, const u64 *restrict src, int t) {
    if (t == 0) { for (int i = 0; i < W; i++) dst[i] |= src[i]; return; }
    int q = t >> 6, s = t & 63;
    if (s == 0) {
        dst[W - 1] |= src[W - 1 - q] & lastmask;
        for (int i = W - 2; i >= q; i--) dst[i] |= src[i - q];
    } else {
        int sh = 64 - s;
        { int j = W - 1 - q; u64 v = src[j] << s; if (j >= 1) v |= src[j - 1] >> sh; dst[W - 1] |= v & lastmask; }
        for (int i = W - 2; i > q; i--) dst[i] |= (src[i - q] << s) | (src[i - q - 1] >> sh);
        if (q < W - 1) dst[q] |= src[0] << s;
    }
    int u = N - t, q2 = u >> 6, s2 = u & 63;
    if (s2 == 0) {
        for (int i = 0; i + q2 < W; i++) dst[i] |= src[i + q2];
    } else {
        int sh2 = 64 - s2, i;
        for (i = 0; i + q2 + 1 < W; i++) dst[i] |= (src[i + q2] >> s2) | (src[i + q2 + 1] << sh2);
        if (i + q2 < W) dst[i] |= src[i + q2] >> s2;
    }
}

static long level_scalar(void) {
    unsigned st = ++stampc;
    for (int y = 0; y < M; y++) {
        if (!rowflag_cur[y]) continue;
        u64 *dp[16]; int tt[16];
        for (int q = 0; q < narcs; q++) {
            int y2 = y + arc_y[q]; if (y2 >= M) y2 -= M;
            if (rowstamp[y2] != st) { memset(nxt + (size_t)y2 * W, 0, (size_t)W * 8); rowstamp[y2] = st; }
            dp[q] = nxt + (size_t)y2 * W; tt[q] = arcT[q][y];
        }
        const u64 *src = cur + (size_t)y * W;
        for (int w = 0; w < W; w++) {
            u64 bits = src[w];
            while (bits) {
                int b = __builtin_ctzll(bits); bits &= bits - 1;
                int x = (w << 6) + b;
                for (int q = 0; q < narcs; q++) { int x2 = x + tt[q]; if (x2 >= N) x2 -= N; dp[q][x2 >> 6] |= 1ULL << (x2 & 63); }
            }
        }
    }
    long cnt = 0;
    for (int d = 0; d < M; d++) {
        if (rowstamp[d] != st) { rowflag_nxt[d] = 0; continue; }
        u64 *nx = nxt + (size_t)d * W, *vs = vis + (size_t)d * W, any = 0;
        for (int w = 0; w < W; w++) { u64 v = nx[w] & ~vs[w]; nx[w] = v; vs[w] |= v; cnt += __builtin_popcountll(v); any |= v; }
        rowflag_nxt[d] = any != 0;
    }
    return cnt;
}

static long level_bp(void) {
    long cnt = 0;
    for (int d = 0; d < M; d++) {
        int any_src = 0;
        u64 *nx = nxt + (size_t)d * W;
        for (int q = 0; q < narcs; q++) {
            int ys = d - arc_y[q]; if (ys < 0) ys += M;
            if (!rowflag_cur[ys]) continue;
            if (!any_src) { memset(nx, 0, (size_t)W * 8); any_src = 1; }
            rot_or(nx, cur + (size_t)ys * W, arcT[q][ys]);
        }
        if (!any_src) { rowflag_nxt[d] = 0; continue; }
        u64 *vs = vis + (size_t)d * W, any = 0;
        for (int w = 0; w < W; w++) { u64 v = nx[w] & ~vs[w]; nx[w] = v; vs[w] |= v; cnt += __builtin_popcountll(v); any |= v; }
        rowflag_nxt[d] = any != 0;
    }
    return cnt;
}

/* last level: unreached count, row by row, stopping once it exceeds bound */
static long level_last(ll bound) {
    long unreached = 0;
    u64 *buf = nxt;
    for (int d = 0; d < M; d++) {
        int any_src = 0;
        for (int q = 0; q < narcs; q++) {
            int ys = d - arc_y[q]; if (ys < 0) ys += M;
            if (!rowflag_cur[ys]) continue;
            if (!any_src) { memset(buf, 0, (size_t)W * 8); any_src = 1; }
            rot_or(buf, cur + (size_t)ys * W, arcT[q][ys]);
        }
        const u64 *vs = vis + (size_t)d * W;
        if (!any_src) {
            for (int w = 0; w < W; w++) unreached += __builtin_popcountll(~vs[w] & (w == W - 1 ? lastmask : ~0ULL));
        } else {
            for (int w = 0; w < W; w++) unreached += __builtin_popcountll(~(vs[w] | buf[w]) & (w == W - 1 ? lastmask : ~0ULL));
        }
        if (unreached > bound) return unreached;
    }
    return unreached;
}

/* number of vertices at distance > D from the identity: exact when <= bound, otherwise some value > bound */
static long bfs(ll bound) {
    nbfs++;
    memset(vis, 0, (size_t)M * W * 8);
    memset(rowflag_cur, 0, M);
    memset(cur, 0, (size_t)W * 8);
    cur[0] = 1; vis[0] = 1; rowflag_cur[0] = 1;
    ll reached = 1, cntprev = 1;
    for (int lev = 1; lev <= D; lev++) {
        if (lev == D) return level_last(bound);
        long cnt = cntprev < scalar_thresh ? level_scalar() : level_bp();
        reached += cnt; cntprev = cnt;
        if (cnt == 0) return n - reached;
        if (reached + cnt * capacity[D - lev] < n - bound) return bound + 1;   /* (deg-1)-ary tree bound */
        u64 *t = cur; cur = nxt; nxt = t;
        uint8_t *tf = rowflag_cur; rowflag_cur = rowflag_nxt; rowflag_nxt = tf;
    }
    return n - reached;
}

static void set_arc(int a, int ax, int ay) { arc_x[a] = ax; arc_y[a] = ay; for (int y = 0; y < M; y++) arcT[a][y] = (int)(((ll)rpow[y] * ax) % N); }
static void set_slots(const int *sel) {
    set_arc(0, 0, M / 2);
    for (int k = 0; k < K; k++) {
        int x = cls_x[sel[k]], y = cls_y[sel[k]];
        set_arc(1 + 2 * k, x, y);
        set_arc(2 + 2 * k, inv_x(x, y), inv_y(y));
    }
}

/* ---------- plain BFS over the multiplication (cross-check of the bit-parallel engine) ---------- */
static long plain_defect(const int *sel) {
    uint8_t *dist = malloc(n); ll *queue = malloc(n * sizeof(ll));
    memset(dist, 255, n);
    ll gens[16]; int ng = 0;
    gens[ng++] = (ll)N * (M / 2);
    for (int k = 0; k < K; k++) { int x = cls_x[sel[k]], y = cls_y[sel[k]]; gens[ng++] = x + (ll)N * y; gens[ng++] = inv_x(x, y) + (ll)N * inv_y(y); }
    ll qh = 0, qt = 0; dist[0] = 0; queue[qt++] = 0;
    while (qh < qt) {
        ll g = queue[qh++];
        if (dist[g] >= D) continue;
        for (int a = 0; a < ng; a++) { ll h = gmul(g, gens[a]); if (dist[h] == 255) { dist[h] = dist[g] + 1; queue[qt++] = h; } }
    }
    long un = n - qt;
    free(dist); free(queue);
    return un;
}

/* ---------- reduced enumeration ---------- */
static int *stab[8], nstab[8], *omin[8], trivial[8];
static int sel[8];
static long best = -1, nbest = 0, ncand = 0, nfound = 0;
static int list_T = 0, quiet = 0, fixed_bound = -1;
static int brute_mode = 0;
static int c0lo = 0, c0hi = -1, c1lo = 0, c1hi = -1, count_mode = 0, progress_level = 0;
static double t_start;
static uint8_t *emitted = NULL;      /* --check: bitmap over K-subset ranks */
static ll *binom_tab = NULL;         /* binom_tab[c * (K+1) + k] = C(c, k) */

static ll subset_rank(const int *s) {   /* combinatorial number system, s strictly increasing */
    ll r = 0; for (int k = 0; k < K; k++) r += binom_tab[(ll)s[k] * (K + 1) + (k + 1)]; return r;
}

static void print_set(const char *tag, long d) {
    printf("%s %ld :", tag, d);
    for (int k = 0; k < K; k++) printf(" %d %d", cls_x[sel[k]], cls_y[sel[k]]);
    printf("\n"); fflush(stdout);
}

static void evaluate(void) {
    ncand++;
    if (emitted) emitted[subset_rank(sel)] = 1;
    if (count_mode) return;
    set_slots(sel);
    ll bound = fixed_bound >= 0 ? fixed_bound : (best < 0 ? n : best);
    if (list_T > bound) bound = list_T;
    long d = bfs(bound);
    if (d > bound) return;
    if (d <= list_T) { print_set(d == 0 ? "FOUND" : "NEAR", d); if (d == 0) nfound++; }
    if (best < 0 || d < best) { best = d; nbest = 1; if (!quiet) { printf("  best %ld after %ld candidates (%.1f s)\n", best, ncand, now() - t_start); fflush(stdout); } }
    else if (d == best) nbest++;
}

/* generating set of the subgroup given as an element list (indices into A): random elements, each time re-closing
 * under composition until the closure has exactly cnt elements (so the set provably generates the subgroup).
 * Returns the number of generators written to g[] (indices into A). */
static u64 rs_gen = 0x2545F4914F6CDD1DULL;
static inline u64 rnd_gen(void) { rs_gen ^= rs_gen << 13; rs_gen ^= rs_gen >> 7; rs_gen ^= rs_gen << 17; return rs_gen; }
static int *hk_i, *hk_j, *hk_l; static size_t hcap;   /* open-addressing hash set of triples */
static int gens_of_list(const int *list, int cnt, int *g, int maxg) {
    size_t cap = 16; while (cap < (size_t)cnt * 2) cap <<= 1;
    if (cap > hcap) { hcap = cap; hk_i = realloc(hk_i, cap * sizeof(int)); hk_j = realloc(hk_j, cap * sizeof(int)); hk_l = realloc(hk_l, cap * sizeof(int)); }
    int *qi = malloc(cnt * sizeof(int)), *qj = malloc(cnt * sizeof(int)), *ql = malloc(cnt * sizeof(int));
    int ng = 0; long csize = 1;
    for (;;) {
        /* closure of the current generators */
        for (size_t h = 0; h < cap; h++) hk_i[h] = -1;
        long qh = 0, qt = 0;
        #define HPUT(i, j, l, inserted) do { size_t h = ((size_t)(i) * 0x9E3779B1ULL ^ (size_t)(j) * 0x85EBCA77ULL ^ (size_t)(l) * 0xC2B2AE3DULL) & (cap - 1); \
            inserted = 0; for (;;) { if (hk_i[h] < 0) { hk_i[h] = (i); hk_j[h] = (j); hk_l[h] = (l); inserted = 1; break; } \
            if (hk_i[h] == (i) && hk_j[h] == (j) && hk_l[h] == (l)) break; h = (h + 1) & (cap - 1); } } while (0)
        int ins; HPUT(1, 0, 1, ins); qi[qt] = 1; qj[qt] = 0; ql[qt] = 1; qt++;
        while (qh < qt) {
            int i1 = qi[qh], j1 = qj[qh], l1 = ql[qh]; qh++;
            for (int t = 0; t < ng; t++) {
                int a = g[t], i2 = A_i[a], j2 = A_j[a], l2 = A_l[a];
                int i = (int)(((ll)i1 * i2) % N), j = (int)((((ll)i2 * j1) + (ll)j2 * sb[l1]) % N), l = (int)(((ll)l1 * l2) % M);
                HPUT(i, j, l, ins);
                if (ins) { if (qt >= cnt) { fprintf(stderr, "stabiliser closure exceeds its size\n"); exit(3); } qi[qt] = i; qj[qt] = j; ql[qt] = l; qt++; }
            }
        }
        csize = qt;
        if (csize == cnt) break;
        if (ng == maxg) { fprintf(stderr, "too many generators\n"); exit(3); }
        /* a random element of the list outside the closure */
        int pick = -1;
        for (int tries = 0; tries < 64 && pick < 0; tries++) {
            int t = (int)(rnd_gen() % (u64)cnt); int a = list[t];
            size_t h = ((size_t)A_i[a] * 0x9E3779B1ULL ^ (size_t)A_j[a] * 0x85EBCA77ULL ^ (size_t)A_l[a] * 0xC2B2AE3DULL) & (cap - 1);
            int found = 0;
            for (;;) { if (hk_i[h] < 0) break; if (hk_i[h] == A_i[a] && hk_j[h] == A_j[a] && hk_l[h] == A_l[a]) { found = 1; break; } h = (h + 1) & (cap - 1); }
            if (!found) pick = a;
        }
        if (pick < 0) {   /* scan */
            for (int t = 0; t < cnt && pick < 0; t++) {
                int a = list[t];
                size_t h = ((size_t)A_i[a] * 0x9E3779B1ULL ^ (size_t)A_j[a] * 0x85EBCA77ULL ^ (size_t)A_l[a] * 0xC2B2AE3DULL) & (cap - 1);
                int found = 0;
                for (;;) { if (hk_i[h] < 0) break; if (hk_i[h] == A_i[a] && hk_j[h] == A_j[a] && hk_l[h] == A_l[a]) { found = 1; break; } h = (h + 1) & (cap - 1); }
                if (!found) pick = a;
            }
        }
        if (pick < 0) { fprintf(stderr, "no element outside the closure although |closure| < cnt\n"); exit(3); }
        g[ng++] = pick;
        #undef HPUT
    }
    free(qi); free(qj); free(ql);
    return ng;
}

static void compute_omin(int k) {
    /* omin[k][c] = min over stab[k] of img_cls(c); large stabilisers through a verified generating set + union-find */
    if (nstab[k] <= 48) {
        for (int c = 0; c < Cn; c++) omin[k][c] = c;
        for (int t = 0; t < nstab[k]; t++) {
            int a = stab[k][t];
            if (A_i[a] == 1 && A_j[a] == 0 && A_l[a] == 1) continue;
            for (int c = 0; c < Cn; c++) { int v = img_cls(a, c); if (v < omin[k][c]) omin[k][c] = v; }
        }
        return;
    }
    int g[64]; int ng = gens_of_list(stab[k], nstab[k], g, 64);
    uf = malloc(Cn * sizeof(int));
    for (int c = 0; c < Cn; c++) uf[c] = c;
    for (int t = 0; t < ng; t++)
        for (int c = 0; c < Cn; c++) {
            int d = img_cls(g[t], c);
            int ra = uf_find(c), rb = uf_find(d);
            if (ra != rb) { if (ra < rb) uf[rb] = ra; else uf[ra] = rb; }
        }
    int *first = malloc(Cn * sizeof(int));
    for (int c = 0; c < Cn; c++) first[c] = -1;
    for (int c = 0; c < Cn; c++) { int r = uf_find(c); if (first[r] < 0) first[r] = c; omin[k][c] = first[r]; }
    free(first); free(uf);
}

static void progress(int k, int c) {
    printf("PROGRESS level=%d c=%d candidates=%ld best=%ld time=%.1f\n", k, c, ncand, best, now() - t_start);
    fflush(stdout);
}

static void recurse(int k) {
    int lo = k == 0 ? 0 : sel[k - 1] + 1, hi = Cn;
    if (k == 0) { if (c0lo > lo) lo = c0lo; if (c0hi < hi) hi = c0hi; }
    if (k == 1) { if (c1lo > lo) lo = c1lo; if (c1hi < hi) hi = c1hi; }
    for (int c = lo; c < hi; c++) {
        int ok = 1;
        for (int m = 0; m < k; m++) { if (trivial[m]) break; if (omin[m][c] < sel[m]) { ok = 0; break; } }
        if (!ok) continue;
        if (!trivial[k] && omin[k][c] != c) continue;
        sel[k] = c;
        if (k == K - 1) { evaluate(); continue; }
        /* stabiliser of c inside stab[k] */
        if (trivial[k]) { trivial[k + 1] = 1; nstab[k + 1] = 1; stab[k + 1][0] = stab[k][0]; }
        else {
            int m = 0;
            for (int t = 0; t < nstab[k]; t++) { int a = stab[k][t]; if (img_cls(a, c) == c) stab[k + 1][m++] = a; }
            nstab[k + 1] = m;
            trivial[k + 1] = (m == 1);
            if (!trivial[k + 1]) compute_omin(k + 1);
        }
        recurse(k + 1);
        if (k == progress_level && !emitted) progress(k, c);
    }
}

static void brute(void) {
    /* all K-subsets in lexicographic order */
    for (int k = 0; k < K; k++) sel[k] = k;
    if (Cn < K) return;
    for (;;) {
        evaluate();
        int k = K - 1;
        while (k >= 0 && sel[k] == Cn - K + k) k--;
        if (k < 0) break;
        sel[k]++;
        for (int m = k + 1; m < K; m++) sel[m] = sel[m - 1] + 1;
    }
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: sweep_odd N M s K D [--bound b] [--sample R seed] [--brute] [--check] [--list L] [--quiet]\n"); return 2; }
    N = atoi(argv[1]); M = atoi(argv[2]); S = atoi(argv[3]); K = atoi(argv[4]); D = atoi(argv[5]);
    long sampleR = 0, screenR = 0; unsigned seed = 1; int check = 0;
    for (int a = 6; a < argc; a++) {
        if (!strcmp(argv[a], "--bound")) fixed_bound = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--sample")) { sampleR = atol(argv[++a]); seed = (unsigned)atol(argv[++a]); }
        else if (!strcmp(argv[a], "--c0screen")) { screenR = atol(argv[++a]); seed = (unsigned)atol(argv[++a]); }
        else if (!strcmp(argv[a], "--brute")) brute_mode = 1;
        else if (!strcmp(argv[a], "--check")) check = 1;
        else if (!strcmp(argv[a], "--list")) list_T = atoi(argv[++a]);
        else if (!strcmp(argv[a], "--quiet")) quiet = 1;
        else if (!strcmp(argv[a], "--c0")) { c0lo = atoi(argv[++a]); c0hi = atoi(argv[++a]); }
        else if (!strcmp(argv[a], "--c1")) { c1lo = atoi(argv[++a]); c1hi = atoi(argv[++a]); }
        else if (!strcmp(argv[a], "--count")) count_mode = 1;
        else { fprintf(stderr, "unknown option %s\n", argv[a]); return 2; }
    }
    if (N < 3 || (N & 1) == 0 || M < 2 || (M & 1) || K < 1 || K > 6 || D < 1) { fprintf(stderr, "bad parameters\n"); return 2; }
    n = (ll)N * M; deg = 2 * K + 1; W = (N + 63) / 64; lastmask = (N % 64) ? ((1ULL << (N % 64)) - 1) : ~0ULL;
    rpow = malloc((M + 1) * sizeof(int)); sb = malloc((M + 1) * sizeof(int));
    rpow[0] = 1 % N; sb[0] = 0;
    for (int j = 1; j <= M; j++) { rpow[j] = (int)(((ll)rpow[j - 1] * S) % N); sb[j] = (int)((sb[j - 1] + rpow[j - 1]) % N); }
    if (rpow[M] != 1 % N) { fprintf(stderr, "s^M != 1 mod N\n"); return 2; }
    int ord = 1; while (rpow[ord] != 1 % N) ord++;
    int p = igcd(N, (int)((1 + (ll)rpow[M / 2]) % N)); /* number of involutions = gcd(N, 1 + s^(M/2)) */
    capacity[0] = 0; { ll pw = 1; for (int r = 1; r < 64; r++) { pw *= (deg - 1); capacity[r] = capacity[r - 1] + pw; if (capacity[r] > (ll)4e18 / deg) { for (int r2 = r + 1; r2 < 64; r2++) capacity[r2] = capacity[r]; break; } } }
    build_classes();
    if (c0hi < 0 || c0hi > Cn) c0hi = Cn;
    if (c1hi < 0 || c1hi > Cn) c1hi = Cn;
    if (sampleR == 0) { build_aut(); build_gens(); }
    vis = calloc((size_t)M * W, 8); cur = calloc((size_t)M * W, 8); nxt = calloc((size_t)M * W, 8);
    rowflag_cur = calloc(M, 1); rowflag_nxt = calloc(M, 1); rowstamp = calloc(M, sizeof(unsigned));
    for (int a = 0; a < 16; a++) arcT[a] = malloc(M * sizeof(int));
    narcs = 2 * K + 1;
    scalar_thresh = (long)M * W / 2; if (scalar_thresh < 8) scalar_thresh = 8;
    printf("G = Z_%d x|_%d Z_%d, n = %lld, ord(s) = %d (kernel %d), involutions p = %d (q = %d), degree %d, D = %d, classes %d, |A_t| = %d\n",
           N, S, M, n, ord, M / ord, p, N / p, deg, D, Cn, nA);
    if (sampleR == 0) printf("  A_t generated by %d maps (closure verified)\n", ngen);
    fflush(stdout);
    double t0 = now();
    if (sampleR > 0) {
        u64 rs = 0x9E3779B97F4A7C15ULL ^ ((u64)seed * 0x100000001B3ULL);
        long *vals = malloc(sampleR * sizeof(long)); double sum = 0;
        for (long r = 0; r < sampleR; r++) {
            for (int k = 0; k < K; k++) {
                int again;
                do { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; sel[k] = (int)(rs % (u64)Cn); again = 0; for (int m = 0; m < k; m++) if (sel[m] == sel[k]) again = 1; } while (again);
            }
            set_slots(sel);
            long d = bfs(n);
            vals[r] = d; sum += d;
        }
        /* sort */
        for (long i = 1; i < sampleR; i++) { long v = vals[i], j = i - 1; while (j >= 0 && vals[j] > v) { vals[j + 1] = vals[j]; j--; } vals[j + 1] = v; }
        printf("SAMPLE N=%d M=%d s=%d p=%d ord=%d n=%lld R=%ld mean %.2f q1 %ld median %ld q3 %ld min %ld max %ld  (%.2f s, %.1f us/BFS)\n", N, M, S, p, ord, n, sampleR, sum / sampleR,
               vals[sampleR / 4], vals[sampleR / 2], vals[3 * sampleR / 4], vals[0], vals[sampleR - 1], now() - t0, 1e6 * (now() - t0) / sampleR);
        return 0;
    }
    if (check) {
        /* automorphism property */
        int fails = 0;
        for (int a = 0; a < nA; a++) {
            if (img_elem(a, 0, M / 2) != (ll)N * (M / 2)) { printf("CHECK FAIL: phi_%d does not fix t\n", a); fails++; }
            uint8_t *seen = calloc(n, 1); int inj = 1;
            for (ll g = 0; g < n; g++) { ll im = img_elem(a, (int)(g % N), (int)(g / N)); if (seen[im]) { inj = 0; break; } seen[im] = 1; }
            free(seen);
            if (!inj) { printf("CHECK FAIL: phi_%d not injective\n", a); fails++; }
            int hom = 1;
            for (ll g = 0; g < n && hom; g++)
                for (ll h = 0; h < n; h++) {
                    ll gh = gmul(g, h);
                    if (img_elem(a, (int)(gh % N), (int)(gh / N)) != gmul(img_elem(a, (int)(g % N), (int)(g / N)), img_elem(a, (int)(h % N), (int)(h / N)))) { hom = 0; break; }
                }
            if (!hom) { printf("CHECK FAIL: phi_%d (i=%d j=%d l=%d) is not a homomorphism\n", a, A_i[a], A_j[a], A_l[a]); fails++; }
        }
        printf("check: %d maps verified as automorphisms fixing t, %d failures\n", nA, fails);
        /* orbit coverage */
        ll total = 1; for (int k = 0; k < K; k++) total = total * (Cn - k) / (k + 1);
        binom_tab = calloc((size_t)(Cn + 1) * (K + 1), sizeof(ll));
        for (int c = 0; c <= Cn; c++) for (int k = 0; k <= K; k++) {
            if (k == 0) binom_tab[(ll)c * (K + 1) + k] = 1;
            else if (c == 0) binom_tab[(ll)c * (K + 1) + k] = 0;
            else binom_tab[(ll)c * (K + 1) + k] = binom_tab[(ll)(c - 1) * (K + 1) + k - 1] + binom_tab[(ll)(c - 1) * (K + 1) + k];
        }
        emitted = calloc(total, 1);
        for (int k = 0; k < 8; k++) { stab[k] = malloc((nA + 1) * sizeof(int)); omin[k] = malloc((Cn + 1) * sizeof(int)); }
        for (int a = 0; a < nA; a++) stab[0][a] = a;
        nstab[0] = nA; trivial[0] = (nA == 1); if (!trivial[0]) compute_omin0(omin[0]);
        quiet = 1; best = -1; ncand = 0; nbest = 0;
        recurse(0);
        long red_best = best, red_cand = ncand, red_nbest = nbest;
        uint8_t *mark = emitted; emitted = NULL;
        uint8_t *done = calloc(total, 1);
        long orbits = 0, missing = 0; int s[8], im[8];
        for (int k = 0; k < K; k++) s[k] = k;
        if (Cn >= K) for (;;) {
            ll r = subset_rank(s);
            if (!done[r]) {
                orbits++; int hit = 0;
                for (int a = 0; a < nA; a++) {
                    for (int k = 0; k < K; k++) im[k] = img_cls(a, s[k]);
                    for (int p2 = 1; p2 < K; p2++) { int v = im[p2], q = p2 - 1; while (q >= 0 && im[q] > v) { im[q + 1] = im[q]; q--; } im[q + 1] = v; }
                    ll r2 = subset_rank(im); done[r2] = 1; if (mark[r2]) hit = 1;
                }
                if (!hit) { missing++; if (missing <= 5) { printf("CHECK FAIL: orbit of {"); for (int k = 0; k < K; k++) printf(" (%d,%d)", cls_x[s[k]], cls_y[s[k]]); printf(" } has no emitted member\n"); } }
            }
            int k = K - 1;
            while (k >= 0 && s[k] == Cn - K + k) k--;
            if (k < 0) break;
            s[k]++;
            for (int m = k + 1; m < K; m++) s[m] = s[m - 1] + 1;
        }
        printf("check: %lld subsets, %ld orbits under A_t, %ld emitted (redundancy %.3f), %ld orbits without an emitted member\n",
               total, orbits, red_cand, orbits ? (double)red_cand / orbits : 0.0, missing);
        /* brute-force minimum (bit-parallel) and plain-BFS cross-check of the minimisers */
        best = -1; ncand = 0; nbest = 0; list_T = 0;
        brute();
        printf("check: reduced minimum %ld (%ld sets), brute-force minimum %ld (%ld sets) over %ld sets: %s\n",
               red_best, red_nbest, best, nbest, ncand, red_best == best ? "AGREE" : "DISAGREE");
        if (red_best != best) fails++;
        /* plain BFS cross-check on 200 random subsets */
        u64 rs = 0x9E3779B97F4A7C15ULL; int mism = 0;
        for (int r = 0; r < 200 && Cn >= K; r++) {
            for (int k = 0; k < K; k++) { int again; do { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; sel[k] = (int)(rs % (u64)Cn); again = 0; for (int m = 0; m < k; m++) if (sel[m] == sel[k]) again = 1; } while (again); }
            set_slots(sel);
            long d1 = bfs(n), d2 = plain_defect(sel);
            if (d1 != d2) { mism++; if (mism <= 3) printf("CHECK FAIL: bit-parallel %ld vs plain %ld\n", d1, d2); }
        }
        printf("check: plain-BFS cross-check on 200 random sets: %d mismatches\n", mism);
        fails += mism + (missing > 0);
        printf("CHECK %s\n", fails ? "FAILED" : "PASSED");
        return fails ? 1 : 0;
    }
    if (brute_mode) {
        brute();
        printf("BRUTE: %ld candidates, minimum defect %ld (%ld sets), %ld with defect 0, %.2f s, %.2f us/candidate\n",
               ncand, best, nbest, nfound, now() - t0, ncand ? 1e6 * (now() - t0) / ncand : 0.0);
        return 0;
    }
    for (int k = 0; k < 8; k++) { stab[k] = malloc((nA + 1) * sizeof(int)); omin[k] = malloc((Cn + 1) * sizeof(int)); }
    for (int a = 0; a < nA; a++) stab[0][a] = a;
    nstab[0] = nA; trivial[0] = (nA == 1); if (!trivial[0]) compute_omin0(omin[0]);
    { long reps = 0; for (int c = 0; c < Cn; c++) if (trivial[0] || omin[0][c] == c) reps++; printf("  %ld A_t-orbits of classes (%.2f s)\n", reps, now() - t0); fflush(stdout); }
    if (screenR > 0) {
        /* per first-generator class: R random completions, full BFS each */
        u64 rs = 0x9E3779B97F4A7C15ULL ^ ((u64)seed * 0x100000001B3ULL);
        long *vals = malloc(screenR * sizeof(long));
        int *osz = calloc(Cn, sizeof(int));
        for (int c = 0; c < Cn; c++) osz[trivial[0] ? c : omin[0][c]]++;
        /* the chunk of c0 consists of the sets whose smallest class (after normalisation) is c0, i.e. all other
         * classes are > c0: the completions are drawn from (c0, Cn) so that the screen samples the chunk's own
         * population (the omin conditions only remove duplicates of sets already counted) */
        for (int c0 = 0; c0 < Cn; c0++) {
            if (!trivial[0] && omin[0][c0] != c0) continue;
            if (Cn - c0 - 1 < K - 1) continue;
            double sum = 0;
            for (long r = 0; r < screenR; r++) {
                sel[0] = c0;
                for (int k = 1; k < K; k++) {
                    int again;
                    do { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; sel[k] = c0 + 1 + (int)(rs % (u64)(Cn - c0 - 1)); again = 0; for (int m = 1; m < k; m++) if (sel[m] == sel[k]) again = 1; } while (again);
                }
                set_slots(sel);
                long d = bfs(n);
                vals[r] = d; sum += d;
            }
            for (long i = 1; i < screenR; i++) { long v = vals[i], j = i - 1; while (j >= 0 && vals[j] > v) { vals[j + 1] = vals[j]; j--; } vals[j + 1] = v; }
            printf("C0SCREEN c0=%d elem=(%d,%d) orbit=%d R=%ld mean %.2f q1 %ld median %ld q3 %ld min %ld\n", c0, cls_x[c0], cls_y[c0], osz[c0], screenR,
                   sum / screenR, vals[screenR / 4], vals[screenR / 2], vals[3 * screenR / 4], vals[0]);
        }
        printf("C0SCREEN done (%.2f s)\n", now() - t0);
        return 0;
    }
    progress_level = (K >= 3 && c0hi == c0lo + 1) ? 1 : 0;
    printf("  chunk: c0 in [%d, %d), c1 in [%d, %d), progress per level-%d subtree%s\n", c0lo, c0hi, c1lo, c1hi, progress_level, count_mode ? " (COUNT ONLY, no BFS)" : "");
    fflush(stdout);
    t_start = now();
    recurse(0);
    double dt = now() - t0;
    printf("RESULT N=%d M=%d s=%d p=%d ord=%d n=%lld K=%d D=%d c0=%d:%d c1=%d:%d candidates=%ld min_defect=%ld n_min=%ld found0=%ld time=%.2f us_per=%.2f%s\n",
           N, M, S, p, ord, n, K, D, c0lo, c0hi, c1lo, c1hi, ncand, best, nbest, nfound, dt, ncand ? 1e6 * dt / ncand : 0.0, count_mode ? " COUNT" : "");
    return 0;
}
