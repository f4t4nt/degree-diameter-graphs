/*
 * bitlift.c -- bit-parallel voltage-lift evaluator and exhaustive prefix sweeper for the metacyclic groups
 *   G(N, M, r, t) = <x, y | x^N = 1, y^M = x^t, y x y^-1 = x^r>   (r^M = 1 mod N, t (r - 1) = 0 mod N; --t T, default 0,
 *   i.e. the split group Z_N x|_r Z_M).
 * Element x^a y^b has index a + N b; product (x1,y1)(x2,y2) = (x1 + r^y1 x2 + t [y1 + y2 >= M], (y1 + y2) mod M); a dart
 * u -> v with voltage g joins the lift vertex (u, h) to (v, h g); lift vertex (v, h) has index v |G| + index(h).
 *
 * The BFS keeps one bitset per base vertex (M rows of N bits, row y = {(x, y)}): right multiplication by (x2, y2) maps
 * row y to row (y + y2) mod M rotated by r^y x2 + t [y + y2 >= M], so one level of the lift BFS is a handful of word
 * rotations per row instead of a random table lookup per dart.  Rows are stored doubled (v | v << N) so a rotation is
 * a single funnel shift per word.  Left multiplication by G is an automorphism transitive on each fibre, so one BFS per
 * base vertex decides the diameter; the defect is the number of lift vertices beyond distance D from the roots.
 *
 * Base spec (--base): "k=<vertices>;deg=<degree>;" followed by slots in order: loop:u (a loop at u, voltage g and g^-1),
 * link:u,v (a link u -> v with voltage g; link*:u,v is a fixed trivial link), semi:u (a semi-edge at u, voltage an
 * involution).  Example (shape 103, degree 5): "k=2;deg=5;loop:0;link*:0,1;link:0,1;link:0,1;loop:1".
 *
 * Modes: --mode prefix --prefix FILE (exhaustive over the file's lines: "v1 .. v_{f-1} lo hi" sweeps the last free slot
 * over the index range [lo, hi), "v1 .. v_f" is one explicit assignment, "v1 .. v_{f-1} -1 h1 h2 ..." lists explicit
 * values of the last free slot), --mode sa (simulated annealing with targeted moves; --T0 --T1 --steps --time --ptarget
 * --start), --eval "[v0,...]" (one assignment, full profiles).  Every candidate with defect <=
 * --report R is printed as "NEAR defect=d voltages=[..]"; a defect-0 candidate as "FOUND ..."; an exhaustive run ends
 * with "EXHAUSTIVE ... candidates=.. complete=1".  --out FILE appends one JSON line per event.
 *
 * Antipode instrumentation (groups with a unique central involution z): --aprec FILE records, per candidate, the BFS
 * level at which each root reaches (root, z); --require-ap rejects a candidate as soon as a root misses its antipode.
 *
 * Build: the row width is a compile-time constant NW = ceil(N/64) words:
 *   gcc -O3 -march=native -DNW=4 -o bitlift_nw4 bitlift.c -lm
 * The default kernels use AVX-512 (KERNEL 1 for NW <= 3 and M >= 8, KERNEL 2 for NW >= 4).  On a machine without
 * AVX-512 add -DKERNEL=0 for the portable scalar kernel (same results, slower).  -DMAXM=<max M> (default 128) and
 * -DMAXK=<max base vertices> (default 4) bound the static tables.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef NW
#define NW 2
#endif
#ifndef MAXM
#define MAXM 128
#endif
#ifndef MAXK
#define MAXK 4
#endif
#define MAXS 24
#define MAXDEG 8
#define MAXDART (MAXK * MAXDEG)

typedef uint32_t gel;

static int N, M, m, n;          /* |Z_N|, |Z_M|, |G| = N*M, lift order k*m */
static long r;
static long tt = 0;             /* y^M = x^tt (non-split parameter; 0 = split group) */
static int rp[MAXM];            /* r^y mod N */
static uint64_t lastmask;
static int qN, bN;              /* N = 64*qN + bN, bN > 0 */
static char groupdesc[256];

static double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* ---------------- group arithmetic ---------------- */
static inline int gx(gel g) { return (int)(g % (gel)N); }
static inline int gy(gel g) { return (int)(g / (gel)N); }
static inline gel gmk(long x, long y) { return (gel)(((x % N + N) % N) + (long)N * ((y % M + M) % M)); }
/* (x1,y1)(x2,y2) = (x1 + r^y1 x2 + tt [y1 + y2 >= M], (y1 + y2) mod M); the inverse of (x, y) is (x', y') with
 * y' = (M - y) mod M and x' = -r^y' (x + tt [y != 0]) (for y != 0 the product's y-sum is exactly M) */
static inline gel gmul(gel g, gel a) {
    int ys = gy(g) + gy(a);
    return gmk(gx(g) + (long)rp[gy(g)] * gx(a) + (ys >= M ? tt : 0), ys);
}
static inline gel ginv(gel g) { int y = (M - gy(g)) % M; return gmk(-(long)rp[y] * (gx(g) + (gy(g) ? tt : 0)), y); }

static void build_meta(void) {
    /* N a multiple of 64 is allowed: bN = 0 makes the doubled row a word-aligned copy (lastmask = all ones) */
    if ((N + 63) / 64 != NW) { fprintf(stderr, "binary compiled for NW=%d words, need %d for N=%d\n", NW, (N + 63) / 64, N); exit(1); }
    if (M > MAXM) { fprintf(stderr, "M too large\n"); exit(1); }
    m = N * M;
    rp[0] = 1 % N;
    for (int y = 1; y < M; y++) rp[y] = (int)(((long)rp[y - 1] * r) % N);
    if (((long)rp[M - 1] * r) % N != 1 % N) { fprintf(stderr, "r^M != 1 mod N\n"); exit(1); }
    tt = ((tt % N) + N) % N;
    if ((tt * ((r - 1) % N)) % N != 0) { fprintf(stderr, "t (r - 1) != 0 mod N: not a group presentation\n"); exit(1); }
    qN = N / 64; bN = N % 64;
    lastmask = bN ? ((1ULL << bN) - 1) : ~0ULL;
    for (int g = 0; g < m; g++) { gel i = ginv((gel)g); if (gmul((gel)g, i) != 0 || gmul(i, (gel)g) != 0) { fprintf(stderr, "inverse check failed\n"); exit(1); } }
    /* associativity on a pseudo-random sample (the t-term is the only non-trivial part of the group law) */
    { uint64_t s = 12345; for (int it = 0; it < 20000; it++) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL; gel a = (gel)((s >> 33) % m);
        s = s * 6364136223846793005ULL + 1442695040888963407ULL; gel b = (gel)((s >> 33) % m);
        s = s * 6364136223846793005ULL + 1442695040888963407ULL; gel c = (gel)((s >> 33) % m);
        if (gmul(gmul(a, b), c) != gmul(a, gmul(b, c))) { fprintf(stderr, "associativity failed\n"); exit(1); } } }
    if (tt) snprintf(groupdesc, sizeof groupdesc, "G(%d,%d,%ld,t=%ld)", N, M, r, tt);
    else snprintf(groupdesc, sizeof groupdesc, "Z%d x|_%ld Z%d", N, r, M);
}

/* antipode instrumentation: the unique involution z of G (must be central), the BFS level at which
 * (root, z) is first reached (apdist[root], -1 if beyond D), per-root miss counters, optional per-candidate records */
static gel zinv; static int have_z = 0;
static int apdist[MAXK];
static long ap_miss[MAXK], ap_allok = 0, ap_cand = 0;
static FILE *aprec = NULL;      /* --aprec FILE: one line per simple candidate: voltages defect apdist... */
static int require_ap = 0;      /* --require-ap: reject a candidate as soon as a root misses its antipode */
static void find_involution(void) {
    int cnt = 0;
    for (int g = 1; g < m; g++) if (ginv((gel)g) == (gel)g) { zinv = (gel)g; cnt++; }
    if (cnt != 1) { fprintf(stderr, "group has %d involutions; antipode instrumentation off\n", cnt); return; }
    /* central? */
    for (int g = 0; g < m; g++) if (gmul((gel)g, zinv) != gmul(zinv, (gel)g)) { fprintf(stderr, "involution not central?!\n"); exit(1); }
    have_z = 1;
    fprintf(stderr, "unique central involution z = x^%d y^%d (index %u)\n", gx(zinv), gy(zinv), (unsigned)zinv);
}

/* ---------------- base graph ---------------- */
enum { LINK = 0, LOOP = 1, SEMI = 2 };
static int deg = 5;
static int k, nslot;
static int st[MAXS], su[MAXS], sv[MAXS], sfixed[MAXS];
static int nfree, freeslot[MAXS];
static gel val[MAXS];
static char basedesc[512];

static void parse_base(const char *spec) {
    char buf[1024]; strncpy(buf, spec, sizeof buf - 1); buf[sizeof buf - 1] = 0;
    snprintf(basedesc, sizeof basedesc, "%s", spec);
    nslot = 0; k = 0; nfree = 0;
    char *save = NULL;
    for (char *tok = strtok_r(buf, ";", &save); tok; tok = strtok_r(NULL, ";", &save)) {
        if (!strncmp(tok, "k=", 2)) { k = atoi(tok + 2); continue; }
        if (!strncmp(tok, "deg=", 4)) { deg = atoi(tok + 4); continue; }
        int fixed = 0, type;
        char *colon = strchr(tok, ':');
        if (!colon) { fprintf(stderr, "bad token %s\n", tok); exit(1); }
        *colon = 0;
        size_t L = strlen(tok);
        if (L && tok[L - 1] == '*') { fixed = 1; tok[L - 1] = 0; }
        if (!strcmp(tok, "link")) type = LINK; else if (!strcmp(tok, "loop")) type = LOOP;
        else if (!strcmp(tok, "semi")) type = SEMI; else { fprintf(stderr, "bad type %s\n", tok); exit(1); }
        int u = 0, v = 0;
        if (type == LINK) { if (sscanf(colon + 1, "%d,%d", &u, &v) != 2) { fprintf(stderr, "bad link\n"); exit(1); } }
        else { u = v = atoi(colon + 1); }
        if (nslot >= MAXS) { fprintf(stderr, "too many slots\n"); exit(1); }
        st[nslot] = type; su[nslot] = u; sv[nslot] = v; sfixed[nslot] = fixed;
        val[nslot] = 0;
        if (!fixed) freeslot[nfree++] = nslot;
        nslot++;
    }
    if (k < 1 || k > MAXK) { fprintf(stderr, "bad k\n"); exit(1); }
    if (deg < 1 || deg > MAXDEG) { fprintf(stderr, "bad deg\n"); exit(1); }
    int dg[MAXK] = {0};
    for (int s = 0; s < nslot; s++) {
        if (su[s] < 0 || su[s] >= k || sv[s] < 0 || sv[s] >= k) { fprintf(stderr, "bad vertex in slot %d\n", s); exit(1); }
        if (st[s] == LINK) { dg[su[s]]++; dg[sv[s]]++; }
        else if (st[s] == LOOP) dg[su[s]] += 2;
        else dg[su[s]] += 1;
        if (sfixed[s] && st[s] != LINK) { fprintf(stderr, "only links may be fixed\n"); exit(1); }
    }
    for (int u = 0; u < k; u++) if (dg[u] != deg) { fprintf(stderr, "base vertex %d has degree %d, want %d\n", u, dg[u], deg); exit(1); }
}

/* darts: tail, head, slot, direction, voltage, and per-row rotation tables (target row, word offset, bit shift) */
static int ndart;
static int dtail[MAXDART], dhead[MAXDART], dslot[MAXDART], ddir[MAXDART];
static gel dvolt[MAXDART];
static uint16_t dry[MAXDART][MAXM]; static uint8_t dq[MAXDART][MAXM], db[MAXDART][MAXM];   /* dry: row index < MAXM (may exceed 255) */
static gel slot_cached[MAXS];
static int darts_built = 0;
static void dart_tables_kernel(int d);   /* kernel-specific tables derived from dry/dq/db */

static void dart_tables(int d) {
    int x2 = gx(dvolt[d]), y2 = gy(dvolt[d]);
    for (int y = 0; y < M; y++) {
        int s = (int)(((long)rp[y] * x2 + (y + y2 >= M ? tt : 0)) % N);
        int t = N - s;                     /* take bits [t, t+N) of v | v<<N: rotation left by s */
        dry[d][y] = (uint16_t)((y + y2) % M);
        dq[d][y] = (uint8_t)(t >> 6);
        db[d][y] = (uint8_t)(t & 63);
    }
}

static void setup_darts(void) {
    if (!darts_built) {
        ndart = 0;
        for (int s = 0; s < nslot; s++) {
            int u = su[s], v = sv[s];
            if (st[s] == LINK) {
                dtail[ndart] = u; dhead[ndart] = v; dslot[ndart] = s; ddir[ndart] = 0; ndart++;
                dtail[ndart] = v; dhead[ndart] = u; dslot[ndart] = s; ddir[ndart] = 1; ndart++;
            } else if (st[s] == LOOP) {
                dtail[ndart] = u; dhead[ndart] = u; dslot[ndart] = s; ddir[ndart] = 0; ndart++;
                dtail[ndart] = u; dhead[ndart] = u; dslot[ndart] = s; ddir[ndart] = 1; ndart++;
            } else { dtail[ndart] = u; dhead[ndart] = u; dslot[ndart] = s; ddir[ndart] = 0; ndart++; }
        }
        for (int s = 0; s < nslot; s++) slot_cached[s] = (gel)-1;
        darts_built = 1;
    }
    for (int d = 0; d < ndart; d++) {
        int s = dslot[d];
        if (slot_cached[s] == val[s]) continue;
        dvolt[d] = ddir[d] ? ginv(val[s]) : val[s];
        dart_tables(d);
        dart_tables_kernel(d);
    }
    for (int s = 0; s < nslot; s++) slot_cached[s] = val[s];
}

/* simple-graph check at (u,e) for every u (regular cover => suffices) */
static int is_simple(void) {
    for (int u = 0; u < k; u++) {
        uint64_t c[MAXDEG]; int t = 0;
        for (int d = 0; d < ndart; d++) {
            if (dtail[d] != u) continue;
            c[t] = (uint64_t)dhead[d] * m + dvolt[d];
            if (c[t] == (uint64_t)u * m) return 0;
            for (int s = 0; s < t; s++) if (c[s] == c[t]) return 0;
            t++;
        }
    }
    return 1;
}

/* ---------------- bit-parallel BFS ---------------- */
/* KERNEL 0: scalar, row layout.  KERNEL 1: AVX-512 column layout (8 rows per vector, per-lane funnel shifts;
 * NW <= 3, M >= 8).  KERNEL 2: AVX-512 row layout (vectors along the words of a row; NW >= 4). */
#ifndef KERNEL
#if NW <= 3
#define KERNEL 1
#else
#define KERNEL 2
#endif
#endif
#if KERNEL > 0
#include <immintrin.h>
#endif
static int D;
static long bfs_calls = 0;
static long levcnt[64];

static inline uint64_t shrd64(uint64_t lo, uint64_t hi, int b) {
    return (uint64_t)(((((unsigned __int128)hi) << 64) | lo) >> b);
}

#if KERNEL == 1
/* column layout: X[u][word][row], rows padded to MP = 8*MB; NFc has 8 overflow rows (targets y+y2 >= M) folded
 * back after the darts */
#define MPMAX MAXM
#define MBMAX (MAXM / 8)
static uint64_t Fc[MAXK][NW][MPMAX] __attribute__((aligned(64)));
static uint64_t Vc[MAXK][NW][MPMAX] __attribute__((aligned(64)));
static uint64_t NFc[MAXK][NW][MPMAX + 8] __attribute__((aligned(64)));
static uint64_t Wc[MAXK][2 * NW + 1][MPMAX] __attribute__((aligned(64)));
static uint64_t dbv[MAXDART][MBMAX][8] __attribute__((aligned(64)));   /* per-lane shift counts */
static __mmask8 dqm[MAXDART][MBMAX][NW];                                 /* lanes with word offset q' (q' >= 1) */
static int dry0[MAXDART][MBMAX];                                         /* target row of lane 0 */
static int MB, MP;

static void dart_tables_kernel(int d) {
    int y2 = gy(dvolt[d]);
    for (int kb = 0; kb < MB; kb++) {
        dry0[d][kb] = (8 * kb + y2) % M;
        for (int qq = 0; qq < NW; qq++) dqm[d][kb][qq] = 0;
        for (int l = 0; l < 8; l++) {
            int y = 8 * kb + l;
            if (y >= M) { dbv[d][kb][l] = 0; continue; }
            dbv[d][kb][l] = (uint64_t)db[d][y];
            dqm[d][kb][dq[d][y]] |= (__mmask8)(1 << l);
        }
    }
}

static int bfs_root(int root, int *uDm1) {
    bfs_calls++;
    memset(Fc, 0, sizeof(uint64_t) * NW * MPMAX * k);
    memset(Vc, 0, sizeof(uint64_t) * NW * MPMAX * k);
    Fc[root][0][0] = 1; Vc[root][0][0] = 1;
    long reached = 1; *uDm1 = -1;
    levcnt[0] = 1; apdist[root] = -1;
    const __m512i lmask = _mm512_set1_epi64((long long)lastmask);
    const __m512i vbN = _mm512_set1_epi64(bN), vbN2 = _mm512_set1_epi64(64 - bN);
    for (int lev = 1; lev <= D; lev++) {
        for (int u = 0; u < k; u++)
            for (int kb = 0; kb < MB; kb++) {
                __m512i f[NW], w[2 * NW + 1];
                for (int i = 0; i < NW; i++) f[i] = _mm512_load_si512((const void *)&Fc[u][i][8 * kb]);
                for (int i = 0; i < NW; i++) w[i] = f[i];
                for (int i = NW; i < 2 * NW + 1; i++) w[i] = _mm512_setzero_si512();
                for (int i = 0; i < NW; i++) {
                    w[i + qN] = _mm512_or_si512(w[i + qN], _mm512_sllv_epi64(f[i], vbN));
                    w[i + qN + 1] = _mm512_or_si512(w[i + qN + 1], _mm512_srlv_epi64(f[i], vbN2));
                }
                for (int i = 0; i < 2 * NW + 1; i++) _mm512_store_si512((void *)&Wc[u][i][8 * kb], w[i]);
            }
        memset(NFc, 0, sizeof(uint64_t) * NW * (MPMAX + 8) * k);
        for (int d = 0; d < ndart; d++) {
            int u = dtail[d], v = dhead[d];
            for (int kb = 0; kb < MB; kb++) {
                __m512i x[NW + 1];
                for (int c = 0; c <= NW; c++) {
                    x[c] = _mm512_load_si512((const void *)&Wc[u][c][8 * kb]);
                    for (int qq = 1; qq < NW; qq++)
                        x[c] = _mm512_mask_blend_epi64(dqm[d][kb][qq], x[c], _mm512_load_si512((const void *)&Wc[u][c + qq][8 * kb]));
                }
                __m512i bv = _mm512_load_si512((const void *)dbv[d][kb]);
                int ry0 = dry0[d][kb];
                for (int i = 0; i < NW; i++) {
                    __m512i o = _mm512_shrdv_epi64(x[i], x[i + 1], bv);
                    if (i == NW - 1) o = _mm512_and_si512(o, lmask);
                    uint64_t *t = &NFc[v][i][ry0];
                    _mm512_storeu_si512((void *)t, _mm512_or_si512(_mm512_loadu_si512((const void *)t), o));
                }
            }
        }
        for (int u = 0; u < k; u++)
            for (int i = 0; i < NW; i++)
                for (int y = 0; y < 8; y++) { NFc[u][i][y] |= NFc[u][i][y + M]; NFc[u][i][y + M] = 0; }
        __m512i vcnt = _mm512_setzero_si512();
        for (int u = 0; u < k; u++)
            for (int i = 0; i < NW; i++)
                for (int kb = 0; kb < MB; kb++) {
                    __m512i nf = _mm512_load_si512((const void *)&NFc[u][i][8 * kb]);
                    __m512i vv = _mm512_load_si512((const void *)&Vc[u][i][8 * kb]);
                    __m512i xx = _mm512_andnot_si512(vv, nf);
                    _mm512_store_si512((void *)&Vc[u][i][8 * kb], _mm512_or_si512(vv, xx));
                    _mm512_store_si512((void *)&Fc[u][i][8 * kb], xx);
                    vcnt = _mm512_add_epi64(vcnt, _mm512_popcnt_epi64(xx));
                }
        long cnt = _mm512_reduce_add_epi64(vcnt);
        levcnt[lev] = cnt;
        reached += cnt;
        if (have_z && apdist[root] < 0 && ((Vc[root][gx(zinv) >> 6][gy(zinv)] >> (gx(zinv) & 63)) & 1)) apdist[root] = lev;
        if (lev == D - 1) *uDm1 = (int)(n - reached);
        if (cnt == 0) { for (int l = lev + 1; l <= D; l++) levcnt[l] = 0; break; }
    }
    if (*uDm1 < 0) *uDm1 = (int)(n - reached);
    return (int)(n - reached);
}

static int unreached_list(uint32_t *buf) {
    int c = 0;
    for (int u = 0; u < k; u++)
        for (int y = 0; y < M; y++)
            for (int i = 0; i < NW; i++) {
                uint64_t x = ~Vc[u][i][y];
                if (i == NW - 1) x &= lastmask;
                while (x) { int bit = __builtin_ctzll(x); x &= x - 1; buf[c++] = (uint32_t)u * m + (gel)(64 * i + bit) + (gel)N * y; }
            }
    return c;
}

#else /* KERNEL 0 or 2: row layout */
#if KERNEL == 2
#define NWP (((NW + 7) / 8) * 8)
#else
#define NWP NW
#endif
static void dart_tables_kernel(int d) { (void)d; }
typedef uint64_t row_t[NWP];
typedef uint64_t drow_t[3 * NWP];
static row_t V[MAXK][MAXM] __attribute__((aligned(64)));
static row_t F[MAXK][MAXM] __attribute__((aligned(64)));
static row_t NF[MAXK][MAXM] __attribute__((aligned(64)));
static drow_t W[MAXK][MAXM] __attribute__((aligned(64)));
static uint64_t wmask[NWP] __attribute__((aligned(64)));
static uint8_t nz[MAXK][MAXM];

static int bfs_root(int root, int *uDm1) {
    bfs_calls++;
    memset(V, 0, sizeof(row_t) * MAXM * k);
    memset(F, 0, sizeof(row_t) * MAXM * k);
    memset(nz, 0, sizeof nz);
    V[root][0][0] = 1; F[root][0][0] = 1; nz[root][0] = 1;
    long reached = 1; *uDm1 = -1;
    levcnt[0] = 1; apdist[root] = -1;
    for (int lev = 1; lev <= D; lev++) {
        for (int u = 0; u < k; u++)
            for (int y = 0; y < M; y++) {
                if (!nz[u][y]) continue;
                uint64_t *w = W[u][y]; const uint64_t *v = F[u][y];
                for (int i = 0; i < NW; i++) w[i] = v[i];
                for (int i = NW; i < 3 * NWP; i++) w[i] = 0;
                for (int i = 0; i < NW; i++) { w[i + qN] |= v[i] << bN; if (bN) w[i + qN + 1] |= v[i] >> (64 - bN); }
            }
        memset(NF, 0, sizeof(row_t) * MAXM * k);
        for (int d = 0; d < ndart; d++) {
            int u = dtail[d], v = dhead[d];
            const uint16_t *ry = dry[d]; const uint8_t *rq = dq[d], *rb = db[d];
            for (int y = 0; y < M; y++) {
                if (!nz[u][y]) continue;
                const uint64_t *w = W[u][y] + rq[y];
                uint64_t *o = NF[v][ry[y]];
                int b = rb[y];
#if KERNEL == 2
                __m512i bv = _mm512_set1_epi64(b);
                for (int i = 0; i < NWP; i += 8) {
                    __m512i lo = _mm512_loadu_si512((const void *)(w + i)), hi = _mm512_loadu_si512((const void *)(w + i + 1));
                    __m512i r8 = _mm512_and_si512(_mm512_shrdv_epi64(lo, hi, bv), _mm512_load_si512((const void *)(wmask + i)));
                    _mm512_store_si512((void *)(o + i), _mm512_or_si512(_mm512_load_si512((const void *)(o + i)), r8));
                }
#else
                for (int i = 0; i < NW - 1; i++) o[i] |= shrd64(w[i], w[i + 1], b);
                o[NW - 1] |= shrd64(w[NW - 1], w[NW], b) & lastmask;
#endif
            }
        }
        long cnt = 0;
#if KERNEL == 2
        __m512i vcnt = _mm512_setzero_si512();
        for (int u = 0; u < k; u++)
            for (int y = 0; y < M; y++) {
                __m512i anyv = _mm512_setzero_si512();
                for (int i = 0; i < NWP; i += 8) {
                    __m512i nf = _mm512_load_si512((const void *)&NF[u][y][i]);
                    __m512i vv = _mm512_load_si512((const void *)&V[u][y][i]);
                    __m512i xx = _mm512_andnot_si512(vv, nf);
                    _mm512_store_si512((void *)&V[u][y][i], _mm512_or_si512(vv, xx));
                    _mm512_store_si512((void *)&F[u][y][i], xx);
                    anyv = _mm512_or_si512(anyv, xx);
                    vcnt = _mm512_add_epi64(vcnt, _mm512_popcnt_epi64(xx));
                }
                nz[u][y] = _mm512_test_epi64_mask(anyv, anyv) != 0;
            }
        cnt = _mm512_reduce_add_epi64(vcnt);
#else
        for (int u = 0; u < k; u++)
            for (int y = 0; y < M; y++) {
                uint64_t any = 0;
                for (int i = 0; i < NW; i++) {
                    uint64_t x = NF[u][y][i] & ~V[u][y][i];
                    V[u][y][i] |= x; F[u][y][i] = x; any |= x;
                    cnt += __builtin_popcountll(x);
                }
                nz[u][y] = any != 0;
            }
#endif
        levcnt[lev] = cnt;
        reached += cnt;
        if (have_z && apdist[root] < 0 && ((V[root][gy(zinv)][gx(zinv) >> 6] >> (gx(zinv) & 63)) & 1)) apdist[root] = lev;
        if (lev == D - 1) *uDm1 = (int)(n - reached);
        if (cnt == 0) { for (int l = lev + 1; l <= D; l++) levcnt[l] = 0; break; }
    }
    if (*uDm1 < 0) *uDm1 = (int)(n - reached);
    return (int)(n - reached);
}

/* vertices at distance > D from the last BFS root, as lift vertex codes u*m + g */
static int unreached_list(uint32_t *buf) {
    int c = 0;
    for (int u = 0; u < k; u++)
        for (int y = 0; y < M; y++)
            for (int i = 0; i < NW; i++) {
                uint64_t x = ~V[u][y][i];
                if (i == NW - 1) x &= lastmask;
                while (x) { int bit = __builtin_ctzll(x); x &= x - 1; buf[c++] = (uint32_t)u * m + (gel)(64 * i + bit) + (gel)N * y; }
            }
    return c;
}
#endif /* KERNEL */

static void kernel_init(void) {
#if KERNEL == 1
    if (M < 8) { fprintf(stderr, "column kernel needs M >= 8\n"); exit(1); }
    MB = (M + 7) / 8; MP = 8 * MB;
    if (MP + 8 > MPMAX + 8) { fprintf(stderr, "M too large for column kernel\n"); exit(1); }
#elif KERNEL == 2
    for (int i = 0; i < NWP; i++) wmask[i] = i < NW - 1 ? ~0ULL : (i == NW - 1 ? lastmask : 0);
#endif
}

/* total defect over all roots; roots after the first are skipped as soon as the sum exceeds thr (thr < 0: exact).
 * Returns the exact defect if <= thr (or thr < 0), else a value > thr.  Not-simple: 1<<40, *cost 1<<50. */
static long W_cost = 64;
static uint32_t *unrbuf[MAXK]; static int nunr_tmp[MAXK];
static long eval_full(long *cost, long thr, int want_unreached) {
    setup_darts();
    if (!is_simple()) { *cost = (long)1 << 50; return (long)1 << 40; }
    long tot = 0, c = 0;
    for (int rt = 0; rt < k; rt++) {
        int u1; int u = bfs_root(rt, &u1);
        if (want_unreached) nunr_tmp[rt] = unreached_list(unrbuf[rt]);
        tot += u; c += W_cost * u + u1;
        if (thr >= 0 && tot > thr) { *cost = (long)1 << 50; return tot; }
    }
    *cost = c;
    return tot;
}

/* ---------------- allowed sets ---------------- */
/* a non-tree link may carry the identity unless a tree link joins the same pair (double edge) */
static int par_tree(int slot) {
    for (int t = 0; t < nslot; t++)
        if (sfixed[t] && st[t] == LINK && su[t] == su[slot] && sv[t] == sv[slot]) return 1;
    return 0;
}
static int allowed(int slot, gel g) {
    switch (st[slot]) {
        case LINK: return g != 0 || !par_tree(slot);
        case LOOP: return g != 0 && ginv(g) != g;
        default:   return g != 0 && ginv(g) == g;
    }
}
static gel *allowlist[MAXS]; static int nallow[MAXS];
static void build_allow(void) {
    for (int s = 0; s < nslot; s++) {
        allowlist[s] = malloc(m * sizeof(gel)); nallow[s] = 0;
        for (int g = 0; g < m; g++) if (allowed(s, (gel)g)) allowlist[s][nallow[s]++] = (gel)g;
    }
}

/* ---------------- reporting ---------------- */
static FILE *out;
static const char *tag = "";
static void print_assignment(FILE *f) {
    fprintf(f, "[");
    for (int s = 0; s < nslot; s++) fprintf(f, "%s%d", s ? "," : "", (int)val[s]);
    fprintf(f, "]");
}
static int nfound = 0;
static void report_found(void) {
    nfound++;
    long cost; long tot = eval_full(&cost, -1, 0);
    printf("FOUND diam<=%d (defect %ld) group=%s base=%s voltages=", D, tot, groupdesc, basedesc);
    print_assignment(stdout); printf("\n"); fflush(stdout);
    if (out) {
        fprintf(out, "{\"event\":\"found\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"voltages\":",
                tag, D, n, groupdesc, basedesc);
        print_assignment(out); fprintf(out, "}\n"); fflush(out);
    }
}

/* ---------------- exhaustive (prefix file) ---------------- */
static long best_cap = 64;
static long best_defect = (long)1 << 40; static gel best_val[MAXS];
static long report_thr = -1;
static long ncand = 0, nsimplecand = 0;
static double t_start, t_limit = 1e18;
static int stop_on_found = 0, stopped = 0;
static long hist0[64];          /* log2 histogram of u_D(root 0) over simple candidates */

static void consider_candidate(void) {
    ncand++;
    if ((ncand & 4095) == 0 && now() - t_start > t_limit) stopped = 1;
    setup_darts();
    if (!is_simple()) return;
    nsimplecand++;
    long want = best_defect > report_thr ? best_defect : report_thr;
    int thr = want < best_cap ? (int)want : (int)best_cap;
    int u1;
    if (aprec) {
        /* statistics mode: every root is evaluated in full; one record per simple candidate */
        int tot = 0, allok = 1;
        for (int rt = 0; rt < k; rt++) {
            tot += bfs_root(rt, &u1);
            if (apdist[rt] < 0) { ap_miss[rt]++; allok = 0; }
        }
        ap_cand++; if (allok) ap_allok++;
        { int b = tot == 0 ? 0 : (int)(log2((double)tot) + 1); if (b > 63) b = 63; hist0[b]++; }
        for (int s = 0; s < nslot; s++) fprintf(aprec, "%d ", (int)val[s]);
        fprintf(aprec, "%d", tot);
        for (int rt = 0; rt < k; rt++) fprintf(aprec, " %d", apdist[rt]);
        fputc('\n', aprec);
        if (tot == 0) { report_found(); if (stop_on_found) exit(0); }
        if (tot < best_defect) { best_defect = tot; memcpy(best_val, val, sizeof val); }
        if (tot <= report_thr) { printf("NEAR defect=%d voltages=", tot); print_assignment(stdout); printf("\n"); fflush(stdout); }
        return;
    }
    int tot = bfs_root(0, &u1);
    { int b = tot == 0 ? 0 : (int)(log2((double)tot) + 1); if (b > 63) b = 63; hist0[b]++; }
    if (have_z) { ap_cand++; if (apdist[0] < 0) { ap_miss[0]++; if (require_ap) return; } }
    if (tot > thr) return;
    for (int rt = 1; rt < k; rt++) {
        tot += bfs_root(rt, &u1);
        if (have_z && apdist[rt] < 0) { ap_miss[rt]++; if (require_ap) return; }
        if (tot > thr) return;
    }
    if (tot == 0) { report_found(); if (stop_on_found) exit(0); }
    if (tot < best_defect) { best_defect = tot; memcpy(best_val, val, sizeof val); }
    if (tot <= report_thr) { printf("NEAR defect=%d voltages=", tot); print_assignment(stdout); printf("\n"); fflush(stdout); }
}

static void loop_rest(int from) {
    if (stopped) return;
    if (from == nfree) { consider_candidate(); return; }
    int s = freeslot[from];
    for (int i = 0; i < nallow[s] && !stopped; i++) { val[s] = allowlist[s][i]; loop_rest(from + 1); }
}

static void run_prefix(const char *prefixfile) {
    build_allow();
    t_start = now();
    FILE *f = fopen(prefixfile, "r"); if (!f) { perror(prefixfile); exit(1); }
    static char line[1 << 16]; long nprefix = 0;
    long vals[MAXS + 4];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        int cnt = 0; char *p = line;
        while (*p && cnt < MAXS + 3) {
            char *end; long v = strtol(p, &end, 10);
            if (end == p) break;
            vals[cnt++] = v; p = end;
            /* compact line "v1 .. v_{nfree-1} -1 h1 h2 ...": explicit values of the last free slot follow the -1 */
            if (cnt == nfree && v == -1) break;
        }
        if (cnt == 0) continue;
        nprefix++;
        if (cnt == nfree && vals[nfree - 1] == -1) {
            for (int i = 0; i < nfree - 1; i++) val[freeslot[i]] = (gel)vals[i];
            int s = freeslot[nfree - 1];
            while (*p && !stopped) {
                char *end; long g = strtol(p, &end, 10);
                if (end == p) break;
                p = end;
                if (g < 0 || g >= m || !allowed(s, (gel)g)) continue;
                val[s] = (gel)g;
                consider_candidate();
            }
        } else if (cnt == nfree + 1) {
            for (int i = 0; i < nfree - 1; i++) val[freeslot[i]] = (gel)vals[i];
            int s = freeslot[nfree - 1];
            long lo = vals[nfree - 1], hi = vals[nfree];
            if (lo < 0) lo = 0;
            if (hi > m) hi = m;
            for (long g = lo; g < hi && !stopped; g++) {
                if (!allowed(s, (gel)g)) continue;
                val[s] = (gel)g;
                consider_candidate();
            }
        } else if (cnt <= nfree) {
            for (int i = 0; i < cnt; i++) val[freeslot[i]] = (gel)vals[i];
            loop_rest(cnt);
        } else { fprintf(stderr, "bad prefix line (%d numbers, nfree=%d)\n", cnt, nfree); exit(1); }
        if (stopped) { printf("TIME LIMIT after %ld prefix lines\n", nprefix); break; }
    }
    fclose(f);
    double el = now() - t_start;
    printf("EXHAUSTIVE group=%s base=%s D=%d n=%d prefixes=%ld candidates=%ld simple=%ld bfs=%ld time=%.1fs rate=%.0f cand/s best_defect=%ld best=",
           groupdesc, basedesc, D, n, nprefix, ncand, nsimplecand, bfs_calls, el, ncand / (el > 0 ? el : 1), best_defect);
    memcpy(val, best_val, sizeof val); print_assignment(stdout); printf(" found=%d complete=%d\n", nfound, !stopped);
    printf("  u_D(root 0) histogram (log2 bins): ");
    for (int b = 0; b < 64; b++) if (hist0[b]) printf("[%s%ld]=%ld ", b ? "<2^" : "", b ? (long)b : 0L, hist0[b]);
    printf("\n");
    if (have_z) {
        printf("  antipode: candidates with a root-0 BFS %ld, root-0 antipode missed %ld (%.2f%%)", ap_cand, ap_miss[0], 100.0 * ap_miss[0] / (ap_cand ? ap_cand : 1));
        for (int rt = 1; rt < k; rt++) printf(", root-%d missed %ld (of the candidates whose BFS from root %d ran)", rt, ap_miss[rt], rt);
        if (aprec) printf("; all antipodes reached: %ld of %ld (%.2f%%)", ap_allok, ap_cand, 100.0 * ap_allok / (ap_cand ? ap_cand : 1));
        printf("\n");
    }
    if (out) {
        fprintf(out, "{\"event\":\"exhaustive\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"prefixes\":%ld,\"candidates\":%ld,\"simple\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, nprefix, ncand, nsimplecand, el, best_defect);
        print_assignment(out); fprintf(out, ",\"found\":%d,\"complete\":%d,\"ap_cand\":%ld,\"ap_miss0\":%ld,\"ap_allok\":%ld,\"require_ap\":%d,\"t\":%ld,\"evaluator\":\"bitlift\"}\n", nfound, !stopped, ap_cand, ap_miss[0], ap_allok, require_ap, tt); fflush(out);
    }
}

/* ---------------- SA ---------------- */
static uint64_t rng_s[2] = {0x9E3779B97F4A7C15ULL, 0xD1B54A32D192ED03ULL};
static inline uint64_t rotl(uint64_t x, int rr) { return (x << rr) | (x >> (64 - rr)); }
static uint64_t rnd(void) {
    uint64_t s0 = rng_s[0], s1 = rng_s[1], rr = s0 + s1;
    s1 ^= s0; rng_s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); rng_s[1] = rotl(s1, 36); return rr;
}
static double rndu(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static void seed_rng(uint64_t s) { rng_s[0] = s * 0x9E3779B97F4A7C15ULL + 1; rng_s[1] = (s ^ 0xDEADBEEFCAFEULL) * 0xBF58476D1CE4E5B9ULL + 7; for (int i = 0; i < 20; i++) rnd(); }
static void random_assignment(void) {
    for (int i = 0; i < nfree; i++) { int s = freeslot[i]; val[s] = allowlist[s][rnd() % nallow[s]]; }
}

static uint32_t *cur_unr[MAXK]; static int ncur_unr[MAXK];
static double ptarget = 0.5;
static int have_start = 0; static gel start_val[MAXS];
static void parse_start(const char *vs) {
    char buf[2048]; strncpy(buf, vs, sizeof buf - 1); buf[sizeof buf - 1] = 0;
    int s = 0; char *save = NULL;
    for (char *tok = strtok_r(buf, ",[] ", &save); tok && s < nslot; tok = strtok_r(NULL, ",[] ", &save)) start_val[s++] = (gel)atoi(tok);
    if (s != nslot) { fprintf(stderr, "--start needs %d voltages, got %d\n", nslot, s); exit(1); }
    have_start = 1;
}

/* darts out of base vertex cv, excluding slot s */
static int targeted_move(gel *newval) {
    int totu = 0; for (int rt = 0; rt < k; rt++) totu += ncur_unr[rt];
    if (!totu) return -1;
    int pick = (int)(rnd() % totu), root = 0;
    while (pick >= ncur_unr[root]) { pick -= ncur_unr[root]; root++; }
    uint32_t c = cur_unr[root][pick]; int tv = c / m; gel tg = (gel)(c - tv * m);
    setup_darts();
    for (int tries = 0; tries < 20; tries++) {
        int s = freeslot[rnd() % nfree];
        int l = 1 + (int)(rnd() % D), j = (int)(rnd() % l);
        int cv = root, dir = 0, ok = 1; gel p = 0, q = 0;
        for (int t = 0; t < l && ok; t++) {
            if (t == j) {
                if (st[s] == LINK) { if (cv == su[s]) { dir = 0; cv = sv[s]; } else if (cv == sv[s]) { dir = 1; cv = su[s]; } else ok = 0; }
                else if (cv != su[s]) ok = 0;
            } else {
                int cand[MAXDART], nc = 0;
                for (int d = 0; d < ndart; d++) if (dtail[d] == cv && dslot[d] != s) cand[nc++] = d;
                if (!nc) { ok = 0; break; }
                int d = cand[rnd() % nc];
                gel vlt = dvolt[d];
                if (t < j) p = gmul(p, vlt); else q = gmul(q, vlt);
                cv = dhead[d];
            }
        }
        if (!ok || cv != tv) continue;
        gel x = gmul(gmul(ginv(p), tg), ginv(q));
        gel a = dir ? ginv(x) : x;
        if (!allowed(s, a) || a == val[s]) continue;
        *newval = a; return s;
    }
    return -1;
}

static void run_sa(double total_time, double T0, double T1, long steps_per_restart) {
    build_allow();
    for (int rt = 0; rt < k; rt++) { cur_unr[rt] = malloc(n * sizeof(uint32_t)); }
    t_start = now();
    long gbest = (long)1 << 40; gel gbest_val[MAXS]; long restarts = 0, evals = 0, rejected_early = 0, ntarget = 0, ntarget_acc = 0;
    long hist_local[64] = {0};
    long curcost, curtot;
    while (now() - t_start < total_time) {
        if (have_start && restarts == 0) { memcpy(val, start_val, sizeof val); curtot = eval_full(&curcost, -1, 1); if (curtot >= ((long)1 << 40)) { fprintf(stderr, "bad --start\n"); exit(1); } }
        else do { random_assignment(); curtot = eval_full(&curcost, -1, 1); } while (curtot >= ((long)1 << 40));
        for (int rt = 0; rt < k; rt++) { ncur_unr[rt] = nunr_tmp[rt]; memcpy(cur_unr[rt], unrbuf[rt], nunr_tmp[rt] * sizeof(uint32_t)); }
        long lbest = curtot; gel lbest_val[MAXS]; memcpy(lbest_val, val, sizeof val);
        restarts++;
        double lt = log(T1 / T0);
        for (long step = 0; step < steps_per_restart; step++) {
            double T = T0 * exp(lt * (double)step / steps_per_restart);
            int s = -1; gel nv = 0; int targeted = 0;
            if (rndu() < ptarget) { s = targeted_move(&nv); if (s >= 0) { targeted = 1; ntarget++; } }
            if (s < 0) {
                s = freeslot[rnd() % nfree];
                do { nv = allowlist[s][rnd() % nallow[s]]; } while (nv == val[s] && nallow[s] > 1);
            }
            gel old = val[s];
            val[s] = nv;
            long thr = (long)((curcost + 8.0 * T) / W_cost);
            long cost; long tot = eval_full(&cost, thr, 1); evals++;
            if (tot > thr) { val[s] = old; rejected_early++; }
            else {
                if (tot == 0) { report_found(); if (stop_on_found) exit(0); }
                long dd = cost - curcost;
                if (dd <= 0 || rndu() < exp(-(double)dd / T)) {
                    curcost = cost; curtot = tot;
                    if (tot < lbest) { lbest = tot; memcpy(lbest_val, val, sizeof val); }
                    for (int rt = 0; rt < k; rt++) { ncur_unr[rt] = nunr_tmp[rt]; memcpy(cur_unr[rt], unrbuf[rt], nunr_tmp[rt] * sizeof(uint32_t)); }
                    if (targeted) ntarget_acc++;
                } else val[s] = old;
                if (curtot < gbest) { gbest = curtot; memcpy(gbest_val, val, sizeof val); }
            }
            if ((step & 255) == 0 && now() - t_start > total_time) break;
        }
        { int b = lbest == 0 ? 0 : (int)(log2((double)lbest) + 1); if (b > 63) b = 63; hist_local[b]++; }
        {
            gel sv2[MAXS]; memcpy(sv2, val, sizeof val); memcpy(val, lbest_val, sizeof val);
            printf("  restart %ld: local best defect %ld (global %ld), evals %ld (early-rejected %ld, targeted %ld acc %ld), %.0fs local=",
                   restarts, lbest, gbest, evals, rejected_early, ntarget, ntarget_acc, now() - t_start);
            print_assignment(stdout); printf("\n"); fflush(stdout);
            memcpy(val, sv2, sizeof val);
        }
    }
    double el = now() - t_start;
    printf("SA group=%s base=%s D=%d n=%d restarts=%ld evals=%ld time=%.1fs rate=%.0f/s best_defect=%ld best=",
           groupdesc, basedesc, D, n, restarts, evals, el, evals / (el > 0 ? el : 1), gbest);
    memcpy(val, gbest_val, sizeof val); print_assignment(stdout);
    printf("\n  local-best histogram (log2 bins): ");
    for (int b = 0; b < 64; b++) if (hist_local[b]) printf("[%s%ld]=%ld ", b ? "<2^" : "", b ? (long)b : 0L, hist_local[b]);
    printf("\n");
    if (out) {
        fprintf(out, "{\"event\":\"sa\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"restarts\":%ld,\"evals\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, restarts, evals, el, gbest);
        print_assignment(out); fprintf(out, ",\"T0\":%g,\"T1\":%g,\"steps\":%ld,\"W\":%ld,\"ptarget\":%g,\"evaluator\":\"bitlift\"}\n", T0, T1, steps_per_restart, W_cost, ptarget); fflush(out);
    }
}

/* evaluate one assignment: per-root level counts and the diameter */
static void run_eval(const char *vs) {
    char buf[2048]; strncpy(buf, vs, sizeof buf - 1); buf[sizeof buf - 1] = 0;
    int s = 0; char *save = NULL;
    for (char *tok = strtok_r(buf, ",[] ", &save); tok && s < nslot; tok = strtok_r(NULL, ",[] ", &save)) val[s++] = (gel)atoi(tok);
    if (s != nslot) { fprintf(stderr, "need %d voltages, got %d\n", nslot, s); exit(1); }
    setup_darts();
    printf("simple=%d\n", is_simple());
    int saveD = D; D = 60;
    int maxd = 0; long tot = 0;
    for (int rt = 0; rt < k; rt++) {
        int u1; int u = bfs_root(rt, &u1);
        int ecc = 0; for (int l = 0; l <= D; l++) if (levcnt[l]) ecc = l;
        printf("root %d: unreached=%d ecc=%d levels:", rt, u, ecc);
        for (int l = 0; l <= ecc; l++) printf(" %ld", levcnt[l]);
        if (have_z) printf("  antipode distance %d", apdist[rt]);
        printf("\n");
        if (ecc > maxd) maxd = ecc;
        tot += u;
    }
    D = saveD;
    printf("EVAL group=%s base=%s n=%d diameter=%d (unreached=%ld)\n", groupdesc, basedesc, n, tot ? -1 : maxd, tot);
}

int main(int argc, char **argv) {
    const char *base = NULL, *mode = "prefix", *prefix = NULL, *outfn = NULL, *evalv = NULL, *startv = NULL;
    double satime = 60, T0 = 2000, T1 = 20; long steps = 20000; uint64_t seed = 1;
    N = 0; M = 1; r = 1; D = 7;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--meta")) { N = atoi(argv[++i]); M = atoi(argv[++i]); r = atol(argv[++i]); }
        else if (!strcmp(argv[i], "--t")) tt = atol(argv[++i]);
        else if (!strcmp(argv[i], "--base")) base = argv[++i];
        else if (!strcmp(argv[i], "--D")) D = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--mode")) mode = argv[++i];
        else if (!strcmp(argv[i], "--prefix")) prefix = argv[++i];
        else if (!strcmp(argv[i], "--out")) outfn = argv[++i];
        else if (!strcmp(argv[i], "--tag")) tag = argv[++i];
        else if (!strcmp(argv[i], "--time")) { satime = atof(argv[++i]); t_limit = satime; }
        else if (!strcmp(argv[i], "--T0")) T0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T1")) T1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--steps")) steps = atol(argv[++i]);
        else if (!strcmp(argv[i], "--W")) W_cost = atol(argv[++i]);
        else if (!strcmp(argv[i], "--seed")) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--eval")) evalv = argv[++i];
        else if (!strcmp(argv[i], "--stop")) stop_on_found = 1;
        else if (!strcmp(argv[i], "--ptarget")) ptarget = atof(argv[++i]);
        else if (!strcmp(argv[i], "--report")) report_thr = atol(argv[++i]);
        else if (!strcmp(argv[i], "--cap")) best_cap = atol(argv[++i]);
        else if (!strcmp(argv[i], "--start")) startv = argv[++i];
        else if (!strcmp(argv[i], "--aprec")) { aprec = fopen(argv[++i], "w"); if (!aprec) { perror("aprec"); exit(1); } }
        else if (!strcmp(argv[i], "--require-ap")) require_ap = 1;
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); exit(1); }
    }
    if (!base) { fprintf(stderr, "need --base\n"); exit(1); }
    if (N <= 0) { fprintf(stderr, "need --meta N M r\n"); exit(1); }
    build_meta();
    find_involution();
    if ((aprec || require_ap) && !have_z) { fprintf(stderr, "antipode options need a unique involution\n"); exit(1); }
    kernel_init();
    parse_base(base);
    if (startv) parse_start(startv);
    n = k * m;
    for (int rt = 0; rt < k; rt++) unrbuf[rt] = malloc((n + 1) * sizeof(uint32_t));
    if (D < 2 || D > 60) { fprintf(stderr, "bad D\n"); exit(1); }
    if (outfn) out = fopen(outfn, "a");
    seed_rng(seed);
    fprintf(stderr, "group %s order %d, base %s, n=%d, free slots=%d, deg=%d, D=%d, mode=%s (bitlift NW=%d MAXM=%d)\n", groupdesc, m, basedesc, n, nfree, deg, D, mode, NW, MAXM);
    if (evalv) run_eval(evalv);
    else if (!strcmp(mode, "prefix") || !strcmp(mode, "exhaustive")) { if (!prefix) { fprintf(stderr, "need --prefix FILE\n"); exit(1); } run_prefix(prefix); }
    else if (!strcmp(mode, "sa")) run_sa(satime, T0, T1, steps);
    else { fprintf(stderr, "bad mode\n"); exit(1); }
    if (out) fclose(out);
    if (aprec) fclose(aprec);
    return 0;
}
