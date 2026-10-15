/*
 * lift.c -- bit-parallel voltage-lift / Cayley-graph search over metacyclic groups Z_N x|_r Z_M
 *           for graphs of degree d (default 10) and diameter <= D: exhaustive prefix / level-pattern sweeps,
 *           random sampling (--ylist: y-pattern restriction of random assignments), rsearch, sa / ils with
 *           targeted moves, and the completion modes sweep1 (complete single-slot sweeps around --init) and
 *           pair (targeted two-slot sweep around --init).
 *
 * Group: element (x, y) <-> index x + N*y, (x1,y1)(x2,y2) = (x1 + r^y1 x2, y1 + y2), r^M = 1 mod N.
 * Base graph: k vertices, slots link:u,v (optionally link*:u,v fixed to e, spanning tree) and loop:u.
 *        Lift vertex (v, g); a dart u->v with voltage a joins (u, g) to (v, g*a).
 *        A Cayley graph with five inverse pairs is "k=1;loop:0;loop:0;loop:0;loop:0;loop:0".
 *
 * BFS: a fibre is stored as M rows of N bits (row y = elements (., y)).  Right multiplication by
 *        a = (xa, ya) maps row y to row y + ya rotated by r^y * xa, so one BFS level is a handful of
 *        word-rotations per dart and nonzero row of the frontier.  Eccentricity is constant on fibres, so
 *        k BFS runs decide the diameter.
 *
 * Modes:
 *   exhaustive --prefix FILE : lines "E v0 .. v_{f-2} L" (the last free slot runs over the canonical elements
 *                              of level >= L, i.e. y in [0, (M-1)/2], x in [1, (N-1)/2] if y = 0; a complete
 *                              single-slot sweep for L = 0), "T v0 v1 ..." (explicit values) or "P v0 .." (odometer
 *                              over the remaining slots).
 *   random, sa, ils, --eval, --selftest.  sa and ils accept --init "[..]" and use targeted moves (k = 1 only):
 *   a new value s for a slot is aimed at an uncovered element u (distance > D from e) as s = w1^-1 u w2^-1
 *   for short words w1, w2 over the other generators, so that the walk w1 s w2 reaches u.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef uint64_t u64;

static int N, M, m;
static long r;
static int W;             /* 64-bit words per row */
static u64 lastmask;
static int *rp;           /* r^y mod N */
static char groupdesc[256];

static double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline int gmul(int g, int a) {
    int x1 = g % N, y1 = g / N, x2 = a % N, y2 = a / N;
    long x = (x1 + (long)rp[y1] * x2) % N;
    int y = y1 + y2; if (y >= M) y -= M;
    return (int)x + N * y;
}
static inline int ginv(int g) {
    int x = g % N, y = g / N;
    int yi = (M - y) % M;
    long xi = (N - ((long)rp[yi] * x) % N) % N;
    return (int)xi + N * yi;
}

static void build_group(void) {
    m = N * M;
    rp = malloc(M * sizeof(int));
    rp[0] = 1 % N;
    for (int y = 1; y < M; y++) rp[y] = (int)(((long)rp[y - 1] * r) % N);
    if (((long)rp[M - 1] * r) % N != 1 % N) { fprintf(stderr, "r^M != 1 mod N\n"); exit(1); }
    W = (N + 63) / 64;
    lastmask = (N % 64) ? ((1ULL << (N % 64)) - 1) : ~0ULL;
    snprintf(groupdesc, sizeof groupdesc, "Z%d x|_%ld Z%d", N, r, M);
    for (int g = 0; g < m; g += 1 + m / 997) {
        if (gmul(g, ginv(g)) != 0 || gmul(ginv(g), g) != 0) { fprintf(stderr, "inverse check failed\n"); exit(1); }
        int a = (g * 7 + 3) % m, b = (g * 13 + 5) % m;
        if (gmul(gmul(g, a), b) != gmul(g, gmul(a, b))) { fprintf(stderr, "associativity failed\n"); exit(1); }
    }
}

/* ---------------- base graph ---------------- */
#define MAXK 8
#define MAXS 32
#define MAXD 24
enum { LINK = 0, LOOP = 1 };
static int k, deg = 10;
static int nslot, st[MAXS], su[MAXS], sv[MAXS], sfixed[MAXS];
static int nfree, freeslot[MAXS];
static int val[MAXS];
static char basedesc[512];

static void parse_base(const char *spec) {
    char buf[1024]; strncpy(buf, spec, sizeof buf - 1); buf[sizeof buf - 1] = 0;
    snprintf(basedesc, sizeof basedesc, "%s", spec);
    nslot = 0; k = 0; nfree = 0;
    char *save = NULL;
    for (char *tok = strtok_r(buf, ";", &save); tok; tok = strtok_r(NULL, ";", &save)) {
        if (!strncmp(tok, "k=", 2)) { k = atoi(tok + 2); continue; }
        int fixed = 0, type;
        char *colon = strchr(tok, ':');
        if (!colon) { fprintf(stderr, "bad token %s\n", tok); exit(1); }
        *colon = 0;
        size_t L = strlen(tok);
        if (L && tok[L - 1] == '*') { fixed = 1; tok[L - 1] = 0; }
        if (!strcmp(tok, "link")) type = LINK; else if (!strcmp(tok, "loop")) type = LOOP;
        else { fprintf(stderr, "bad type %s (no semi-edges: the generating sets have no involutions)\n", tok); exit(1); }
        int u = 0, v = 0;
        if (type == LINK) { if (sscanf(colon + 1, "%d,%d", &u, &v) != 2) { fprintf(stderr, "bad link\n"); exit(1); } }
        else { u = v = atoi(colon + 1); }
        if (type == LINK && u == v) { fprintf(stderr, "link with u == v\n"); exit(1); }
        st[nslot] = type; su[nslot] = u; sv[nslot] = v; sfixed[nslot] = fixed;
        val[nslot] = 0;
        if (!fixed) freeslot[nfree++] = nslot;
        nslot++;
    }
    if (k < 1 || k > MAXK) { fprintf(stderr, "bad k\n"); exit(1); }
    int dg[MAXK] = {0};
    for (int s = 0; s < nslot; s++) {
        if (st[s] == LINK) { dg[su[s]]++; dg[sv[s]]++; } else dg[su[s]] += 2;
        if (sfixed[s] && st[s] != LINK) { fprintf(stderr, "only links may be fixed\n"); exit(1); }
    }
    for (int u = 0; u < k; u++) if (dg[u] != deg) { fprintf(stderr, "base vertex %d has degree %d, want %d\n", u, dg[u], deg); exit(1); }
}

/* darts: for each slot, one or two darts (from, to, voltage) */
static int ndart;
static int dfrom[2 * MAXS], dto[2 * MAXS], dvolt[2 * MAXS], dslot[2 * MAXS];
static int *dtab[2 * MAXS];  /* dtab[d][y] = r^y * x(a) mod N */
static int dya[2 * MAXS];

static void setup_darts(void) {
    ndart = 0;
    for (int s = 0; s < nslot; s++) {
        int a = val[s];
        dfrom[ndart] = su[s]; dto[ndart] = sv[s]; dslot[ndart] = s; dvolt[ndart++] = a;
        dfrom[ndart] = sv[s]; dto[ndart] = su[s]; dslot[ndart] = s; dvolt[ndart++] = ginv(a);
    }
    for (int d = 0; d < ndart; d++) {
        int xa = dvolt[d] % N; dya[d] = dvolt[d] / N;
        for (int y = 0; y < M; y++) dtab[d][y] = (int)(((long)rp[y] * xa) % N);
    }
}

/* simple-graph check at (u, e) for every u (regular cover => suffices) */
static int is_simple(void) {
    for (int u = 0; u < k; u++) {
        long c[MAXD]; int t = 0;
        for (int d = 0; d < ndart; d++) if (dfrom[d] == u) {
            long code = (long)dto[d] * m + dvolt[d];
            if (code == (long)u * m) return 0;             /* self-loop: voltage e on a loop or on a link to itself */
            for (int s = 0; s < t; s++) if (c[s] == code) return 0;
            c[t++] = code;
        }
    }
    return 1;
}

/* ---------------- bitset BFS ---------------- */
static int n, D;
static u64 *vis, *cur, *nxt;      /* k * M * W words each */
static uint8_t *rowflag, *rowflag2; /* k * M: row of cur is nonzero */
static int FW;                    /* words per fibre = M*W */
static long pw[40];               /* (d-1)^i */
static long bfs_calls = 0;

/* dst |= rotate(src, t): bit x of the N-bit row src goes to bit (x + t) mod N.  0 <= t < N. */
static inline void rot_or(u64 *restrict dst, const u64 *restrict src, int t) {
    if (t == 0) { for (int i = 0; i < W; i++) dst[i] |= src[i]; return; }
    int q = t >> 6, s = t & 63;
    /* part 1: out bit j = src bit j - t for t <= j < N (left shift by t, top word masked to N bits) */
    if (s == 0) {
        dst[W - 1] |= src[W - 1 - q] & lastmask;
        for (int i = W - 2; i >= q; i--) dst[i] |= src[i - q];
    } else {
        int sh = 64 - s;
        { int j = W - 1 - q; u64 v = src[j] << s; if (j >= 1) v |= src[j - 1] >> sh; dst[W - 1] |= v & lastmask; }
        for (int i = W - 2; i > q; i--) dst[i] |= (src[i - q] << s) | (src[i - q - 1] >> sh);
        if (q < W - 1) dst[q] |= src[0] << s;
    }
    /* part 2: out bit j = src bit j + (N - t) for 0 <= j < t (right shift by N - t; src has no bits >= N) */
    int u = N - t, q2 = u >> 6, s2 = u & 63;
    if (s2 == 0) {
        for (int i = 0; i + q2 < W; i++) dst[i] |= src[i + q2];
    } else {
        int sh2 = 64 - s2, i;
        for (i = 0; i + q2 + 1 < W; i++) dst[i] |= (src[i + q2] >> s2) | (src[i + q2 + 1] << sh2);
        if (i + q2 < W) dst[i] |= src[i + q2] >> s2;
    }
}

static void bfs_init(int root) {
    memset(vis, 0, (size_t)k * FW * sizeof(u64));
    memset(cur, 0, (size_t)k * FW * sizeof(u64));
    memset(rowflag, 0, (size_t)k * M);
    vis[(size_t)root * FW] = 1; cur[(size_t)root * FW] = 1;   /* element e = (0,0): row 0, bit 0 */
    rowflag[root * M] = 1;
}

/* one BFS level: nxt = (cur * darts) \ vis; vis |= nxt; swap; returns the number of new vertices */
static long bfs_level(void) {
    memset(nxt, 0, (size_t)k * FW * sizeof(u64));
    for (int d = 0; d < ndart; d++) {
        const u64 *src = cur + (size_t)dfrom[d] * FW;
        u64 *dst = nxt + (size_t)dto[d] * FW;
        const uint8_t *rf = rowflag + dfrom[d] * M;
        const int *tab = dtab[d]; int ya = dya[d];
        for (int y = 0; y < M; y++) {
            if (!rf[y]) continue;
            int y2 = y + ya; if (y2 >= M) y2 -= M;
            rot_or(dst + (size_t)y2 * W, src + (size_t)y * W, tab[y]);
        }
    }
    long cnt = 0;
    for (int f = 0; f < k; f++) for (int y = 0; y < M; y++) {
        u64 *nx = nxt + (size_t)f * FW + (size_t)y * W, *vs = vis + (size_t)f * FW + (size_t)y * W;
        u64 any = 0;
        for (int i = 0; i < W; i++) { u64 v = nx[i] & ~vs[i]; nx[i] = v; vs[i] |= v; cnt += __builtin_popcountll(v); any |= v; }
        rowflag2[f * M + y] = any != 0;
    }
    u64 *tmp = cur; cur = nxt; nxt = tmp;
    uint8_t *t2 = rowflag; rowflag = rowflag2; rowflag2 = t2;
    return cnt;
}

/* BFS from (root, e).  Returns the number of vertices at distance > D (exact unless prune made it stop early,
 * in which case the return value is a positive lower bound); *uDm1 = number at distance > D-1 (when exact). */
static int bfs(int root, int *uDm1, int prune) {
    bfs_calls++;
    bfs_init(root);
    long reached = 1;
    *uDm1 = -1;
    for (int lev = 1; lev <= D; lev++) {
        long cnt = bfs_level();
        reached += cnt;
        if (lev == D - 1) *uDm1 = (int)(n - reached);
        if (cnt == 0) { if (*uDm1 < 0) *uDm1 = (int)(n - reached); return (int)(n - reached); }
        if (reached == n) { if (*uDm1 < 0) *uDm1 = 0; return 0; }
        if (prune && lev < D) {
            long rem = n - reached;
            long bound = cnt * ((pw[D - lev + 1] - (deg - 1)) / (deg - 2));
            if (rem > bound) return (int)(rem - bound);
        }
    }
    if (*uDm1 < 0) *uDm1 = 0;
    return (int)(n - reached);
}

/* eccentricity of (root, e) (or -1 if disconnected); fills dist histogram if given */
static int eccentricity(int root, long *unreached, long *hist) {
    bfs_init(root);
    long reached = 1; int lev = 0;
    if (hist) hist[0] = 1;
    while (reached < n) {
        lev++;
        long cnt = bfs_level();
        reached += cnt;
        if (hist && lev < 40) hist[lev] = cnt;
        if (cnt == 0) { lev--; break; }
    }
    *unreached = n - reached;
    return *unreached ? -1 : lev;
}

static long Wc = 64;
static long eval_full(long *cost) {
    setup_darts();
    if (!is_simple()) { *cost = (long)1 << 50; return (long)1 << 40; }
    long tot = 0, c = 0;
    for (int rt = 0; rt < k; rt++) { int u1; int u = bfs(rt, &u1, 0); tot += u; c += Wc * u + u1; }
    *cost = c;
    return tot;
}

/* ---------------- plain reference BFS (self-test) ---------------- */
static int plain_unreached(int root, int Dmax) {
    uint8_t *dist = malloc(n); memset(dist, 0xff, n);
    int *qa = malloc(n * sizeof(int)), *qb = malloc(n * sizeof(int));
    int rootc = root * m; dist[rootc] = 0; qa[0] = rootc; int ca = 1; long reached = 1;
    for (int lev = 1; lev <= Dmax && ca; lev++) {
        int cb = 0;
        for (int i = 0; i < ca; i++) {
            int v = qa[i] / m, g = qa[i] % m;
            for (int d = 0; d < ndart; d++) if (dfrom[d] == v) {
                int nc = dto[d] * m + gmul(g, dvolt[d]);
                if (dist[nc] == 0xff) { dist[nc] = (uint8_t)lev; qb[cb++] = nc; }
            }
        }
        reached += cb; int *t = qa; qa = qb; qb = t; ca = cb;
    }
    free(dist); free(qa); free(qb);
    return (int)(n - reached);
}

/* ---------------- allowed sets, reporting ---------------- */
static int allowed(int slot, int g) {
    if (st[slot] == LOOP) return g != 0 && ginv(g) != g;
    return 1;
}
static int *allowlist[MAXS]; static int nallow[MAXS];
static void build_allow(void) {
    for (int s = 0; s < nslot; s++) {
        allowlist[s] = malloc(m * sizeof(int)); nallow[s] = 0;
        for (int g = 0; g < m; g++) if (allowed(s, g)) allowlist[s][nallow[s]++] = g;
    }
}

static FILE *out;
static const char *tag = "";
static void print_assignment(FILE *f) {
    fprintf(f, "[");
    for (int s = 0; s < nslot; s++) fprintf(f, "%s%d", s ? "," : "", val[s]);
    fprintf(f, "]");
}
static int nfound = 0;
static void report_found(void) {
    nfound++;
    long cost; long tot = eval_full(&cost);
    printf("FOUND diam<=%d (defect %ld) group=%s base=%s voltages=", D, tot, groupdesc, basedesc);
    print_assignment(stdout); printf("\n"); fflush(stdout);
    if (out) {
        fprintf(out, "{\"event\":\"found\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"voltages\":",
                tag, D, n, groupdesc, basedesc);
        print_assignment(out); fprintf(out, "}\n"); fflush(out);
    }
}

/* ---------------- exhaustive ---------------- */
static long best_defect = (long)1 << 40; static int best_val[MAXS];
static long ncand = 0;
static double t_start, t_limit = 1e18;
static int stop_on_found = 0, stopped = 0;
/* --near T prints every candidate whose exact defect is in 1..T ("NEAR defect=d voltages=[..]"),
 * so that a sweep records all its defect-2 states and not only the running best.  u0 is exact or a lower bound
 * (prune), so a candidate with u0 > max(best, T) cannot be a near state and costs nothing extra. */
static long near_T = 4;

static void consider_candidate(void) {
    ncand++;
    if ((ncand & 1023) == 0 && now() - t_start > t_limit) stopped = 1;
    setup_darts();
    if (!is_simple()) return;
    int u1;
    int u0 = bfs(0, &u1, 1);
    if (u0 == 0) {
        int ok = 1;
        for (int rt = 1; rt < k; rt++) if (bfs(rt, &u1, 1) > 0) { ok = 0; break; }
        if (ok) { report_found(); if (stop_on_found) exit(0); }
    }
    long thr = best_defect > near_T ? best_defect : near_T;
    if (u0 <= thr) {
        long tot = bfs(0, &u1, 0);
        for (int rt = 1; rt < k && tot <= thr; rt++) tot += bfs(rt, &u1, 0);
        if (tot < best_defect) { best_defect = tot; memcpy(best_val, val, sizeof val); }
        if (tot > 0 && tot <= near_T) { printf("NEAR defect=%ld voltages=", tot); print_assignment(stdout); printf("\n"); fflush(stdout); }
    }
}

static void loop_rest(int from) {
    if (stopped) return;
    if (from == nfree) { consider_candidate(); return; }
    int s = freeslot[from];
    for (int i = 0; i < nallow[s] && !stopped; i++) { val[s] = allowlist[s][i]; loop_rest(from + 1); }
}

static int igcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
static long nskipped = 0;   /* candidates skipped because the voltages cannot generate G (k = 1 only) */

/* canonical elements of level >= L: y in [L, (M-1)/2]; for y = 0, x in [1, (N-1)/2] (odd N, M; for even M the
 * row y = M/2 is also run in full, for even N the involution x = N/2, y = 0 is excluded by allowed()).
 * For Cayley graphs (k = 1) the y-coordinates of the generators must generate Z_M, otherwise the graph is
 * disconnected: levels y with gcd(y, y(other generators), M) > 1 are skipped. */
static void enumerate_last(int L) {
    int s = freeslot[nfree - 1];
    int H = M / 2;
    int gy = M;
    if (k == 1) for (int i = 0; i < nfree - 1; i++) gy = igcd(gy, val[freeslot[i]] / N);
    if (L == 0) {
        if (gy == 1) for (int x = 1; x <= (N - 1) / 2 && !stopped; x++) { val[s] = x; consider_candidate(); }
        else nskipped += (N - 1) / 2;
        L = 1;
    }
    for (int y = L; y <= H && !stopped; y++) {
        if (igcd(gy, y) != 1) { nskipped += N; continue; }
        for (int x = 0; x < N && !stopped; x++) { val[s] = x + N * y; if (allowed(s, val[s])) consider_candidate(); }
    }
}

/* Level pattern sweep (Cayley graphs, k = 1, three free slots): "Y ya yb yc" with levels
 * 0 <= y <= (M-1)/2.  a = (0, ya); b = (1, yb) with c = (x, yc) for every x in Z_N; and b = (0, yb) with
 * c in {(0, yc), (1, yc)}.  For prime N this covers every connection set whose generators have levels
 * {ya, yb, yc} up to the automorphisms x -> i x + j S_y (which fix the levels) and inversion: the pair (a, b) has
 * at most these representatives, and the stabiliser of ((0, ya), (1, yb)) is trivial.  (ya = 0 means a = (1, 0).)
 * Each pattern costs N + 2 candidates and is a complete, well-defined family. */
static long npattern_found = 0;
static void sweep_pattern(int ya, int yb, int yc) {
    if (nfree != 3 || k != 1) { fprintf(stderr, "Y lines need k = 1 and 3 free slots\n"); exit(1); }
    int sa = freeslot[0], sb = freeslot[1], sc = freeslot[2];
    int gy = igcd(igcd(ya, yb), igcd(yc, M));
    if (gy != 1) { nskipped += N + 2; return; }
    int f0 = nfound;
    val[sa] = ya == 0 ? 1 : N * ya;
    val[sb] = 1 + N * yb;
    for (int x = 0; x < N && !stopped; x++) { val[sc] = x + N * yc; consider_candidate(); }
    if (ya != 0) {
        val[sb] = N * yb;
        for (int x = 0; x < 2 && !stopped; x++) { val[sc] = x + N * yc; consider_candidate(); }
    }
    if (nfound > f0) npattern_found++;
    printf("  pattern %d %d %d: found %d so far, best defect %ld, %.0fs\n", ya, yb, yc, nfound, best_defect, now() - t_start); fflush(stdout);
}

static void run_exhaustive(const char *prefixfile) {
    build_allow();
    t_start = now();
    long nprefix = 0;
    if (prefixfile) {
        FILE *f = fopen(prefixfile, "r"); if (!f) { perror(prefixfile); exit(1); }
        char line[4096];
        while (fgets(line, sizeof line, f) && !stopped) {
            char kind = line[0];
            if (kind != 'E' && kind != 'T' && kind != 'P' && kind != 'Y') continue;
            char *p = line + 1; long v[64]; int nv = 0;
            while (nv < 64) { char *end; long x = strtol(p, &end, 10); if (end == p) break; v[nv++] = x; p = end; }
            nprefix++;
            if (kind == 'Y') {
                if (nv != 3) { fprintf(stderr, "Y line needs 3 levels\n"); exit(1); }
                sweep_pattern((int)v[0], (int)v[1], (int)v[2]);
            } else if (kind == 'E') {
                if (nv != nfree) { fprintf(stderr, "E line needs %d numbers\n", nfree); exit(1); }
                for (int i = 0; i < nfree - 1; i++) val[freeslot[i]] = (int)v[i];
                enumerate_last((int)v[nfree - 1]);
            } else if (kind == 'T') {
                if (nv != nfree) { fprintf(stderr, "T line needs %d numbers\n", nfree); exit(1); }
                for (int i = 0; i < nfree; i++) val[freeslot[i]] = (int)v[i];
                consider_candidate();
            } else {
                for (int i = 0; i < nv; i++) val[freeslot[i]] = (int)v[i];
                loop_rest(nv);
            }
        }
        fclose(f);
        if (stopped) printf("TIME LIMIT after %ld prefixes\n", nprefix);
        printf("prefixes processed: %ld\n", nprefix);
    } else {
        loop_rest(0);
        if (stopped) printf("TIME LIMIT\n");
    }
    double el = now() - t_start;
    printf("EXHAUSTIVE group=%s base=%s D=%d n=%d candidates=%ld skipped=%ld bfs=%ld time=%.1fs rate=%.0f cand/s best_defect=%ld best=",
           groupdesc, basedesc, D, n, ncand, nskipped, bfs_calls, el, ncand / (el > 0 ? el : 1), best_defect);
    memcpy(val, best_val, sizeof val); print_assignment(stdout); printf(" found=%d complete=%d\n", nfound, !stopped);
    if (out) {
        fprintf(out, "{\"event\":\"exhaustive\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"candidates\":%ld,\"skipped\":%ld,\"prefixes\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, ncand, nskipped, nprefix, el, best_defect);
        print_assignment(out); fprintf(out, ",\"found\":%d,\"complete\":%d}\n", nfound, !stopped); fflush(out);
    }
}

/* ---------------- random / SA / ILS ---------------- */
static uint64_t rng_s[2] = {0x9E3779B97F4A7C15ULL, 0xD1B54A32D192ED03ULL};
static inline uint64_t rotl(uint64_t x, int rr) { return (x << rr) | (x >> (64 - rr)); }
static uint64_t rnd(void) {
    uint64_t s0 = rng_s[0], s1 = rng_s[1], res = s0 + s1;
    s1 ^= s0; rng_s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); rng_s[1] = rotl(s1, 36); return res;
}
static double rndu(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static void seed_rng(uint64_t s) { rng_s[0] = s * 0x9E3779B97F4A7C15ULL + 1; rng_s[1] = (s ^ 0xDEADBEEFCAFEULL) * 0xBF58476D1CE4E5B9ULL + 7; for (int i = 0; i < 20; i++) rnd(); }

/* --ylist FILE: lines of nfree Z_M coordinates (y-patterns); random assignments then take a random line and
 * uniform Z_N coordinates (Cayley graphs, k = 1).  Without it every slot is uniform over its allowed set. */
static int *ylist = NULL, nylist = 0;
static void load_ylist(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
    int cap = 1024; ylist = malloc((size_t)cap * nfree * sizeof(int));
    char line[4096];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        if (nylist == cap) { cap *= 2; ylist = realloc(ylist, (size_t)cap * nfree * sizeof(int)); }
        char *p = line; int nv = 0;
        while (nv < nfree) { char *end; long y = strtol(p, &end, 10); if (end == p) break; ylist[nylist * nfree + nv++] = (int)(((y % M) + M) % M); p = end; }
        if (nv != nfree) { fprintf(stderr, "ylist line needs %d numbers: %s", nfree, line); exit(1); }
        nylist++;
    }
    fclose(f);
    fprintf(stderr, "ylist: %d y-patterns from %s\n", nylist, fn);
}

static void random_assignment(void) {
    if (ylist) {
        int L = (int)(rnd() % nylist);
        for (int i = 0; i < nfree; i++) {
            int s = freeslot[i];
            do { val[s] = (int)(rnd() % N) + N * ylist[L * nfree + i]; } while (!allowed(s, val[s]));
        }
        return;
    }
    for (int i = 0; i < nfree; i++) { int s = freeslot[i]; val[s] = allowlist[s][rnd() % nallow[s]]; }
}

/* Uncovered elements of the current state (k = 1): zero bits of vis after the last BFS of the accepted state.
 * A sample of at most UMAX, taken from a random starting word. */
#define UMAX 64
static int unc[UMAX]; static int nunc = 0;
static void extract_uncovered(void) {
    nunc = 0;
    if (k != 1) return;
    int start = (int)(rnd() % FW);
    for (int i = 0; i < FW && nunc < UMAX; i++) {
        int w = start + i; if (w >= FW) w -= FW;
        u64 z = ~vis[w];
        int y = w / W, wi = w % W;
        if (wi == W - 1) z &= lastmask;
        while (z && nunc < UMAX) {
            int b = __builtin_ctzll(z); z &= z - 1;
            unc[nunc++] = wi * 64 + b + N * y;
        }
    }
}

/* a value for `slot` aimed at a random uncovered element: s = w1^-1 u w2^-1, |w1| + |w2| <= J, words over the
 * darts of the other slots (both directions).  Falls back to a uniform value when there is nothing to aim at. */
static int targetJ = 3;
static int targeted_value(int slot) {
    if (nunc == 0 || k != 1) return allowlist[slot][rnd() % nallow[slot]];
    int u = unc[rnd() % nunc];
    int j = (int)(rnd() % (targetJ + 1));
    int j1 = (int)(rnd() % (j + 1)), j2 = j - j1;
    int others[2 * MAXS], no = 0;
    for (int d = 0; d < ndart; d++) if (dslot[d] != slot) others[no++] = dvolt[d];
    if (no == 0) return allowlist[slot][rnd() % nallow[slot]];
    int w1 = 0, w2 = 0;
    for (int t = 0; t < j1; t++) w1 = gmul(w1, others[rnd() % no]);
    for (int t = 0; t < j2; t++) w2 = gmul(w2, others[rnd() % no]);
    int s = gmul(gmul(ginv(w1), u), ginv(w2));
    if (!allowed(slot, s)) return allowlist[slot][rnd() % nallow[slot]];
    return s;
}

static void parse_assignment(const char *vs) {
    char buf[2048]; strncpy(buf, vs, sizeof buf - 1); buf[sizeof buf - 1] = 0;
    int s = 0; char *save = NULL;
    for (char *tok = strtok_r(buf, ",[] ", &save); tok && s < nslot; tok = strtok_r(NULL, ",[] ", &save)) val[s++] = atoi(tok);
    if (s != nslot) { fprintf(stderr, "need %d voltages, got %d\n", nslot, s); exit(1); }
}

static const char *initv = NULL;
static void initial_state(long *curcost, long *curtot) {
    if (initv) { parse_assignment(initv); *curtot = eval_full(curcost); initv = NULL; if (*curtot < ((long)1 << 40)) return; }
    do { random_assignment(); *curtot = eval_full(curcost); } while (*curtot >= ((long)1 << 40));
}

static void run_random(long samples) {
    build_allow();
    t_start = now();
    long hist[64] = {0}; long best = (long)1 << 40; long nsimple = 0;
    for (long i = 0; i < samples; i++) {
        random_assignment();
        long cost; long tot = eval_full(&cost);
        if (tot >= ((long)1 << 40)) continue;
        nsimple++;
        int b = tot == 0 ? 0 : (int)(log2((double)tot) + 1); if (b > 63) b = 63;
        hist[b]++;
        if (tot < best) { best = tot; memcpy(best_val, val, sizeof val); }
        if (tot == 0) { report_found(); if (stop_on_found) exit(0); }
        if (now() - t_start > t_limit) break;
    }
    double el = now() - t_start;
    printf("RANDOM group=%s base=%s D=%d n=%d samples=%ld simple=%ld time=%.1fs rate=%.0f/s best_defect=%ld best=",
           groupdesc, basedesc, D, n, samples, nsimple, el, nsimple / (el > 0 ? el : 1), best);
    memcpy(val, best_val, sizeof val); print_assignment(stdout); printf("\n  defect histogram (log2 bins): ");
    for (int b = 0; b < 64; b++) if (hist[b]) printf("[%s%ld]=%ld ", b ? "<2^" : "", b ? (long)b : 0L, hist[b]);
    printf("\n");
    if (out) {
        fprintf(out, "{\"event\":\"random\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"samples\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, nsimple, el, best);
        print_assignment(out); fprintf(out, "}\n"); fflush(out);
    }
}

/* Random candidates through the pruned decision BFS (consider_candidate: tree-bound prune, full BFS only
 * when the pruned bound does not exceed the best defect so far).  Sampling mode; with --ylist
 * it samples only the given y-patterns. */
static void run_rsearch(double total_time) {
    build_allow();
    t_start = now();
    t_limit = total_time;
    double t_report = t_start + 60;
    while (!stopped) {
        random_assignment();
        consider_candidate();
        if (now() > t_report) {
            double el = now() - t_start;
            printf("  %.0fs: candidates=%ld rate=%.0f/s best_defect=%ld bfs=%ld\n", el, ncand, ncand / el, best_defect, bfs_calls);
            fflush(stdout);
            t_report += 60;
        }
    }
    double el = now() - t_start;
    printf("RSEARCH group=%s base=%s D=%d n=%d candidates=%ld bfs=%ld time=%.1fs rate=%.0f cand/s best_defect=%ld best=",
           groupdesc, basedesc, D, n, ncand, bfs_calls, el, ncand / (el > 0 ? el : 1), best_defect);
    memcpy(val, best_val, sizeof val); print_assignment(stdout); printf(" found=%d\n", nfound);
    if (out) {
        fprintf(out, "{\"event\":\"rsearch\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"candidates\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, ncand, el, best_defect);
        print_assignment(out); fprintf(out, ",\"found\":%d}\n", nfound); fflush(out);
    }
}

/* Simulated annealing on cost = Wc * (unreached at distance > D) + (unreached at distance > D-1); one slot
 * changes per move, to a targeted value with probability ptarget, else uniform. */
static double ptarget = 0.5;
static long completeT = -1;      /* --complete T: after a restart whose local best defect is <= T, run the complete
                                    single-slot sweeps (sweep_slot) from that local best state */
static void sweep_slot(int s);
static void run_sa(double total_time, double T0, double T1, long steps_per_restart) {
    build_allow();
    t_start = now();
    long gbest = (long)1 << 40; int gbest_val[MAXS]; long restarts = 0, evals = 0, ntarget = 0, acc_t = 0, acc_r = 0;
    long curcost, curtot; int lbest_val[MAXS]; long ncompl = 0;
    while (now() - t_start < total_time) {
        initial_state(&curcost, &curtot);
        extract_uncovered();
        long lbest = curtot; memcpy(lbest_val, val, sizeof val);
        restarts++;
        double lt = log(T1 / T0);
        for (long step = 0; step < steps_per_restart; step++) {
            double T = T0 * exp(lt * (double)step / steps_per_restart);
            int s = freeslot[rnd() % nfree];
            int old = val[s], nv;
            int targeted = rndu() < ptarget;
            if (targeted) { nv = targeted_value(s); ntarget++; }
            else do { nv = allowlist[s][rnd() % nallow[s]]; } while (nv == old && nallow[s] > 1);
            if (nv == old) continue;
            val[s] = nv;
            long cost; long tot = eval_full(&cost); evals++;
            if (tot == 0) { report_found(); if (stop_on_found) exit(0); }
            long dd = cost - curcost;
            if (dd <= 0 || rndu() < exp(-(double)dd / T)) {
                curcost = cost; curtot = tot; if (tot <= lbest) { lbest = tot; memcpy(lbest_val, val, sizeof val); }
                extract_uncovered();
                if (targeted) acc_t++; else acc_r++;
            } else val[s] = old;
            if (curtot < gbest) { gbest = curtot; memcpy(gbest_val, val, sizeof val); }
            if ((step & 255) == 0 && now() - t_start > total_time) break;
        }
        printf("  restart %ld: local best defect %ld (global %ld), evals %ld (targeted %ld, accepted t/r %ld/%ld), %.0fs local=",
               restarts, lbest, gbest, evals, ntarget, acc_t, acc_r, now() - t_start);
        { int sv[MAXS]; memcpy(sv, val, sizeof val); memcpy(val, lbest_val, sizeof val); print_assignment(stdout); memcpy(val, sv, sizeof val); }
        printf("\n"); fflush(stdout);
        if (lbest <= completeT && now() - t_start < total_time) {
            /* exact completion: every single-generator change of the local best state */
            int sv[MAXS]; memcpy(sv, val, sizeof val);
            memcpy(val, lbest_val, sizeof val);
            long c0 = ncand; double t0 = now(); long bd = best_defect; best_defect = lbest;
            for (int i = 0; i < nfree && !stopped; i++) sweep_slot(freeslot[i]);
            ncompl++;
            printf("  completion %ld of the defect-%ld state: %ld candidates in %.0fs, best %ld%s\n", ncompl, lbest, ncand - c0, now() - t0,
                   best_defect, stopped ? " (time limit)" : "");
            fflush(stdout);
            if (best_defect < gbest) { gbest = best_defect; memcpy(gbest_val, best_val, sizeof val); }
            if (bd < best_defect) best_defect = bd;
            memcpy(val, sv, sizeof val);
            evals += ncand - c0;
        }
    }
    double el = now() - t_start;
    printf("SA group=%s base=%s D=%d n=%d restarts=%ld evals=%ld completions=%ld time=%.1fs rate=%.0f/s best_defect=%ld best=",
           groupdesc, basedesc, D, n, restarts, evals, ncompl, el, evals / (el > 0 ? el : 1), gbest);
    memcpy(val, gbest_val, sizeof val); print_assignment(stdout); printf("\n");
    if (out) {
        fprintf(out, "{\"event\":\"sa\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"restarts\":%ld,\"evals\":%ld,\"completions\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, restarts, evals, ncompl, el, gbest);
        print_assignment(out); fprintf(out, "}\n"); fflush(out);
    }
}

/* Iterated local search: coordinate descent where one slot sweep evaluates `sweep` uniform random values
 * (0 = every allowed value) plus, for each sampled uncovered element u, the targeted values w1^-1 u w2^-1 over all
 * words with |w1| + |w2| <= 2 in the other generators; the best value is taken.  Kicks of `kick` random slots at
 * a local optimum, restart after `patience` kicks without improving the local best.  Cost as in SA. */
static int cand_buf[1 << 20]; static int ncand_buf;
static void add_cand(int slot, int s) { if (ncand_buf < (1 << 20) && allowed(slot, s)) cand_buf[ncand_buf++] = s; }
static void targeted_candidates(int slot) {
    if (k != 1) return;
    int others[2 * MAXS], no = 0;
    for (int d = 0; d < ndart; d++) if (dslot[d] != slot) others[no++] = dvolt[d];
    for (int i = 0; i < nunc; i++) {
        int u = unc[i];
        add_cand(slot, u);
        for (int a = 0; a < no; a++) {
            int ia = ginv(others[a]);
            add_cand(slot, gmul(ia, u));            /* s = a^-1 u : walk a, s */
            add_cand(slot, gmul(u, ia));            /* s = u a^-1 : walk s, a */
            for (int b = 0; b < no; b++) {
                int ib = ginv(others[b]);
                add_cand(slot, gmul(gmul(ib, ia), u));    /* (ab)^-1 u */
                add_cand(slot, gmul(gmul(ia, u), ib));    /* a^-1 u b^-1 */
                add_cand(slot, gmul(u, gmul(ib, ia)));    /* u (ab)^-1 */
            }
        }
    }
}

static void run_ils(double total_time, int kick, int patience, long sweep) {
    build_allow();
    t_start = now();
    long gbest = (long)1 << 40; int gbest_val[MAXS]; long restarts = 0, evals = 0, sweeps = 0;
    long curcost, curtot;
    while (now() - t_start < total_time) {
        initial_state(&curcost, &curtot);
        extract_uncovered();
        restarts++;
        long lbest = curcost; int fails = 0;
        while (fails < patience && now() - t_start < total_time) {
            int improved = 1;
            while (improved && now() - t_start < total_time) {
                improved = 0;
                for (int i = 0; i < nfree; i++) {
                    int s = freeslot[(i + restarts) % nfree];
                    int old = val[s], bestv = old; long bestc = curcost, bestt = curtot;
                    setup_darts();
                    ncand_buf = 0;
                    if (sweep == 0) { for (int j = 0; j < nallow[s]; j++) cand_buf[ncand_buf++] = allowlist[s][j]; }
                    else {
                        targeted_candidates(s);
                        for (long j = 0; j < sweep; j++) add_cand(s, allowlist[s][rnd() % nallow[s]]);
                    }
                    for (int j = 0; j < ncand_buf; j++) {
                        if (cand_buf[j] == old) continue;
                        val[s] = cand_buf[j];
                        long cost; long tot = eval_full(&cost); evals++;
                        if (tot == 0) { report_found(); if (stop_on_found) exit(0); }
                        if (cost < bestc) { bestc = cost; bestt = tot; bestv = val[s]; }
                        if ((j & 63) == 0 && now() - t_start > total_time) break;
                    }
                    sweeps++;
                    val[s] = bestv;
                    if (bestc < curcost) {
                        curcost = bestc; curtot = bestt; improved = 1;
                        eval_full(&bestc); extract_uncovered();   /* refresh the uncovered sample of the new state */
                    }
                    if (curtot < gbest) { gbest = curtot; memcpy(gbest_val, val, sizeof val); }
                }
            }
            if (curcost < lbest) { lbest = curcost; fails = 0; } else fails++;
            printf("  restart %ld: local optimum defect %ld cost %ld (global best %ld), sweeps %ld evals %ld, %.0fs\n",
                   restarts, curtot, curcost, gbest, sweeps, evals, now() - t_start);
            fflush(stdout);
            /* kick */
            int saved[MAXS]; memcpy(saved, val, sizeof val);
            long savedcost = curcost, savedtot = curtot;
            for (int t = 0; t < 50; t++) {
                memcpy(val, saved, sizeof val);
                for (int q = 0; q < kick; q++) { int s = freeslot[rnd() % nfree]; val[s] = rndu() < ptarget ? targeted_value(s) : allowlist[s][rnd() % nallow[s]]; }
                curtot = eval_full(&curcost); evals++;
                if (curtot < ((long)1 << 40)) break;
            }
            if (curtot >= ((long)1 << 40)) { memcpy(val, saved, sizeof val); curcost = savedcost; curtot = savedtot; eval_full(&curcost); }
            extract_uncovered();
        }
    }
    double el = now() - t_start;
    printf("ILS group=%s base=%s D=%d n=%d restarts=%ld sweeps=%ld evals=%ld time=%.1fs rate=%.0f/s best_defect=%ld best=",
           groupdesc, basedesc, D, n, restarts, sweeps, evals, el, evals / (el > 0 ? el : 1), gbest);
    memcpy(val, gbest_val, sizeof val); print_assignment(stdout); printf("\n");
    if (out) {
        fprintf(out, "{\"event\":\"ils\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"restarts\":%ld,\"evals\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, restarts, evals, el, gbest);
        print_assignment(out); fprintf(out, "}\n"); fflush(out);
    }
}

/* ---------------- completion around a near miss (--init) ---------------- */

/* Complete single-slot sweeps: every canonical value (one per inverse pair) of one slot with the others fixed at
 * --init, for every free slot in turn.  Exact: if some single generator change of the start state gives
 * diameter <= D, it is found.  Candidates run through consider_candidate (pruned decision BFS). */
static void sweep_slot(int s) {
    int saved = val[s];
    int H = M / 2;
    int gy = M;
    if (k == 1) for (int i = 0; i < nfree; i++) if (freeslot[i] != s) gy = igcd(gy, val[freeslot[i]] / N);
    if (gy == 1) for (int x = 1; x <= (N - 1) / 2 && !stopped; x++) { val[s] = x; consider_candidate(); }
    else nskipped += (N - 1) / 2;
    for (int y = 1; y <= H && !stopped; y++) {
        if (igcd(gy, y) != 1) { nskipped += N; continue; }
        for (int x = 0; x < N && !stopped; x++) { val[s] = x + N * y; if (allowed(s, val[s])) consider_candidate(); }
    }
    val[s] = saved;
}

static void run_sweep1(void) {
    build_allow();
    if (!initv) { fprintf(stderr, "sweep1 needs --init\n"); exit(1); }
    parse_assignment(initv);
    t_start = now();
    long c0; long d0 = eval_full(&c0);
    printf("start state defect %ld cost %ld voltages=", d0, c0); print_assignment(stdout); printf("\n"); fflush(stdout);
    for (int i = 0; i < nfree && !stopped; i++) {
        int s = freeslot[i];
        long before = ncand; double t0 = now();
        sweep_slot(s);
        printf("  slot %d swept: %ld candidates, %.0fs, best defect so far %ld\n", s, ncand - before, now() - t0, best_defect);
        fflush(stdout);
    }
    double el = now() - t_start;
    printf("SWEEP1 group=%s base=%s D=%d n=%d candidates=%ld skipped=%ld time=%.1fs rate=%.0f cand/s best_defect=%ld best=",
           groupdesc, basedesc, D, n, ncand, nskipped, el, ncand / (el > 0 ? el : 1), best_defect);
    memcpy(val, best_val, sizeof val); print_assignment(stdout); printf(" found=%d complete=%d\n", nfound, !stopped);
    if (out) {
        fprintf(out, "{\"event\":\"sweep1\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"init\":\"%s\",\"candidates\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, initv, ncand, el, best_defect);
        print_assignment(out); fprintf(out, ",\"found\":%d,\"complete\":%d}\n", nfound, !stopped); fflush(out);
    }
}

/* All uncovered elements of the last full BFS (k = 1), up to UALL. */
#define UALL 4096
static int uall[UALL]; static int nuall = 0;
static void extract_uncovered_all(void) {
    nuall = 0;
    for (int w = 0; w < FW && nuall < UALL; w++) {
        u64 z = ~vis[w];
        int y = w / W, wi = w % W;
        if (wi == W - 1) z &= lastmask;
        while (z && nuall < UALL) { int b = __builtin_ctzll(z); z &= z - 1; uall[nuall++] = wi * 64 + b + N * y; }
    }
}

/* Targeted values for `slot` aimed at every element of list[0..nl): s = w1^-1 u w2^-1 over all words with
 * |w1| + |w2| <= J in the darts of the other slots (current values), deduplicated. */
static int cmp_int(const void *a, const void *b) { int x = *(const int *)a, y = *(const int *)b; return (x > y) - (x < y); }
static void targeted_candidates_J(int slot, const int *list, int nl, int J) {
    ncand_buf = 0;
    if (k != 1) return;
    int others[2 * MAXS], no = 0;
    for (int d = 0; d < ndart; d++) if (dslot[d] != slot) others[no++] = dvolt[d];
    /* words of length <= J over others: enumerate by odometer; left words w1 and right words w2 with |w1|+|w2| <= J */
    static int wl[1 << 16]; static int wlen[1 << 16]; int nw = 0;
    wl[nw] = 0; wlen[nw] = 0; nw++;
    for (int start = 0, L = 1; L <= J; L++) {
        int end = nw;
        for (int i = start; i < end; i++) for (int a = 0; a < no && nw < (1 << 16); a++) { wl[nw] = gmul(wl[i], others[a]); wlen[nw] = L; nw++; }
        start = end;
    }
    for (int t = 0; t < nl; t++) {
        int u = list[t];
        for (int i = 0; i < nw; i++) {
            int lu = gmul(ginv(wl[i]), u);
            for (int j = 0; j < nw; j++) { if (wlen[i] + wlen[j] > J) continue; add_cand(slot, gmul(lu, ginv(wl[j]))); }
        }
    }
    qsort(cand_buf, ncand_buf, sizeof(int), cmp_int);
    int m2 = 0;
    for (int i = 0; i < ncand_buf; i++) if (i == 0 || cand_buf[i] != cand_buf[i - 1]) cand_buf[m2++] = cand_buf[i];
    ncand_buf = m2;
}

/* Targeted two-slot sweep around --init: for every slot A and every targeted value a of A (aimed at all
 * uncovered elements of the start state, words up to J), evaluate; if the defect is <= T, take the new uncovered
 * set and try, for every other slot B, every targeted value of B aimed at it.  Reports FOUND / the best state. */
static int pairT = 16;
static void run_pair(int J) {
    build_allow();
    if (!initv) { fprintf(stderr, "pair needs --init\n"); exit(1); }
    parse_assignment(initv);
    t_start = now();
    long c0; long d0 = eval_full(&c0);
    extract_uncovered_all();
    printf("start state defect %ld (uncovered listed %d) voltages=", d0, nuall); print_assignment(stdout); printf("\n"); fflush(stdout);
    int u0[UALL]; int nu0 = nuall; memcpy(u0, uall, nu0 * sizeof(int));
    int start[MAXS]; memcpy(start, val, sizeof val);
    long evals = 0, good = 0;
    int *candA = malloc(sizeof(int) * (1 << 20));
    int *candB = malloc(sizeof(int) * (1 << 20));
    for (int ia = 0; ia < nfree && !stopped; ia++) {
        int A = freeslot[ia];
        memcpy(val, start, sizeof val); setup_darts();
        targeted_candidates_J(A, u0, nu0, J);
        int nA = ncand_buf; memcpy(candA, cand_buf, nA * sizeof(int));
        printf("  slot %d: %d targeted values\n", A, nA); fflush(stdout);
        for (int i = 0; i < nA && !stopped; i++) {
            if (candA[i] == start[A]) continue;
            memcpy(val, start, sizeof val); val[A] = candA[i];
            long c; long d = eval_full(&c); evals++;
            if (d >= ((long)1 << 40)) continue;
            if (d == 0) { report_found(); if (stop_on_found) exit(0); }
            if (d < best_defect) { best_defect = d; memcpy(best_val, val, sizeof val); }
            if (d > pairT) continue;
            good++;
            extract_uncovered_all();
            int u1[UALL]; int nu1 = nuall; memcpy(u1, uall, nu1 * sizeof(int));
            int stateA[MAXS]; memcpy(stateA, val, sizeof val);
            for (int ib = 0; ib < nfree && !stopped; ib++) {
                int B = freeslot[ib]; if (B == A) continue;
                memcpy(val, stateA, sizeof val); setup_darts();
                targeted_candidates_J(B, u1, nu1, J);
                int nB = ncand_buf; memcpy(candB, cand_buf, nB * sizeof(int));
                for (int j = 0; j < nB && !stopped; j++) {
                    if (candB[j] == stateA[B]) continue;
                    memcpy(val, stateA, sizeof val); val[B] = candB[j];
                    long c2; long d2 = eval_full(&c2); evals++;
                    if (d2 >= ((long)1 << 40)) continue;
                    if (d2 == 0) { report_found(); if (stop_on_found) exit(0); }
                    if (d2 < best_defect) { best_defect = d2; memcpy(best_val, val, sizeof val); }
                    if ((evals & 255) == 0 && now() - t_start > t_limit) stopped = 1;
                }
            }
            if ((i & 15) == 0) { printf("    %d/%d values of slot %d, evals %ld, good %ld, best %ld, %.0fs\n", i, nA, A, evals, good, best_defect, now() - t_start); fflush(stdout); }
            if (now() - t_start > t_limit) stopped = 1;
        }
    }
    double el = now() - t_start;
    printf("PAIR group=%s base=%s D=%d n=%d J=%d T=%d evals=%ld good=%ld time=%.1fs best_defect=%ld best=",
           groupdesc, basedesc, D, n, J, pairT, evals, good, el, best_defect);
    memcpy(val, best_val, sizeof val); print_assignment(stdout); printf(" found=%d complete=%d\n", nfound, !stopped);
    if (out) {
        fprintf(out, "{\"event\":\"pair\",\"tag\":\"%s\",\"D\":%d,\"n\":%d,\"group\":\"%s\",\"base\":\"%s\",\"init\":\"%s\",\"J\":%d,\"T\":%d,\"evals\":%ld,\"time\":%.1f,\"best_defect\":%ld,\"best\":",
                tag, D, n, groupdesc, basedesc, initv, J, pairT, evals, el, best_defect);
        print_assignment(out); fprintf(out, ",\"found\":%d,\"complete\":%d}\n", nfound, !stopped); fflush(out);
    }
}

static void run_eval(void) {
    setup_darts();
    printf("simple=%d\n", is_simple());
    int maxd = 0; long tot = 0;
    for (int rt = 0; rt < k; rt++) {
        long u; long hist[40] = {0}; int ecc = eccentricity(rt, &u, hist);
        printf("root %d: unreached=%ld ecc=%d dist distribution:", rt, u, ecc);
        for (int i = 0; i < 40 && (i == 0 || hist[i]); i++) printf(" %ld", hist[i]);
        printf("\n");
        if (ecc > maxd) maxd = ecc;
        tot += u;
    }
    printf("EVAL group=%s base=%s n=%d diameter=%d (unreached=%ld)\n", groupdesc, basedesc, n, tot ? -1 : maxd, tot);
}

static void run_selftest(long trials) {
    build_allow();
    int bad = 0;
    for (long t = 0; t < trials; t++) {
        random_assignment();
        setup_darts();
        if (!is_simple()) continue;
        for (int rt = 0; rt < k; rt++) {
            int u1; int u = bfs(rt, &u1, 0);
            int ref = plain_unreached(rt, D);
            int ref1 = plain_unreached(rt, D - 1);
            if (u != ref || u1 != ref1) { bad++; printf("MISMATCH root %d: bitset %d/%d plain %d/%d voltages=", rt, u, u1, ref, ref1); print_assignment(stdout); printf("\n"); }
            int up = bfs(rt, &u1, 1);
            if ((up == 0) != (ref == 0) || up > ref) { bad++; printf("PRUNE MISMATCH root %d: pruned %d exact %d\n", rt, up, ref); }
        }
        /* uncovered extraction: every sampled element must be unreached according to the plain BFS */
        if (k == 1) {
            int u1; bfs(0, &u1, 0); extract_uncovered();
            uint8_t *dist = malloc(n); memset(dist, 0xff, n);
            int *qa = malloc(n * sizeof(int)), *qb = malloc(n * sizeof(int));
            dist[0] = 0; qa[0] = 0; int ca = 1;
            for (int lev = 1; lev <= D && ca; lev++) {
                int cb = 0;
                for (int i = 0; i < ca; i++) for (int d = 0; d < ndart; d++) { int nc = gmul(qa[i], dvolt[d]); if (dist[nc] == 0xff) { dist[nc] = (uint8_t)lev; qb[cb++] = nc; } }
                int *tq = qa; qa = qb; qb = tq; ca = cb;
            }
            for (int i = 0; i < nunc; i++) if (dist[unc[i]] != 0xff) { bad++; printf("UNCOVERED MISMATCH: element %d has distance %d\n", unc[i], dist[unc[i]]); }
            free(dist); free(qa); free(qb);
        }
    }
    printf("SELFTEST group=%s base=%s D=%d trials=%ld mismatches=%d\n", groupdesc, basedesc, D, trials, bad);
}

int main(int argc, char **argv) {
    const char *base = NULL, *mode = "exhaustive", *prefix = NULL, *outfn = NULL, *evalv = NULL, *ylistfn = NULL;
    long samples = 10000; double satime = 60, T0 = 2000, T1 = 20; long steps = 20000; uint64_t seed = 1;
    long selftest = 0; int kick = 2, patience = 20; long sweep = 2000;
    N = 0; M = 1; r = 1; D = 7;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--meta")) { N = atoi(argv[++i]); M = atoi(argv[++i]); r = atol(argv[++i]); }
        else if (!strcmp(argv[i], "--base")) base = argv[++i];
        else if (!strcmp(argv[i], "--deg")) deg = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--D")) D = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--mode")) mode = argv[++i];
        else if (!strcmp(argv[i], "--prefix")) prefix = argv[++i];
        else if (!strcmp(argv[i], "--out")) outfn = argv[++i];
        else if (!strcmp(argv[i], "--tag")) tag = argv[++i];
        else if (!strcmp(argv[i], "--samples")) samples = atol(argv[++i]);
        else if (!strcmp(argv[i], "--time")) { satime = atof(argv[++i]); t_limit = satime; }
        else if (!strcmp(argv[i], "--T0")) T0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T1")) T1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--steps")) steps = atol(argv[++i]);
        else if (!strcmp(argv[i], "--W")) Wc = atol(argv[++i]);
        else if (!strcmp(argv[i], "--seed")) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--eval")) evalv = argv[++i];
        else if (!strcmp(argv[i], "--init")) initv = argv[++i];
        else if (!strcmp(argv[i], "--stop")) stop_on_found = 1;
        else if (!strcmp(argv[i], "--selftest")) selftest = atol(argv[++i]);
        else if (!strcmp(argv[i], "--kick")) kick = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--patience")) patience = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--sweep")) sweep = atol(argv[++i]);
        else if (!strcmp(argv[i], "--ptarget")) ptarget = atof(argv[++i]);
        else if (!strcmp(argv[i], "--J")) targetJ = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--pairT")) pairT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--complete")) completeT = atol(argv[++i]);
        else if (!strcmp(argv[i], "--near")) near_T = atol(argv[++i]);
        else if (!strcmp(argv[i], "--ylist")) ylistfn = argv[++i];
        else { fprintf(stderr, "unknown arg %s\n", argv[i]); exit(1); }
    }
    if (!base || N <= 0) { fprintf(stderr, "need --base and --meta N M r\n"); exit(1); }
    if (deg > MAXD || deg < 3) { fprintf(stderr, "bad degree\n"); exit(1); }
    build_group();
    parse_base(base);
    n = k * m;
    FW = M * W;
    vis = calloc((size_t)k * FW, sizeof(u64)); cur = calloc((size_t)k * FW, sizeof(u64)); nxt = calloc((size_t)k * FW, sizeof(u64));
    rowflag = calloc((size_t)k * M, 1); rowflag2 = calloc((size_t)k * M, 1);
    for (int d = 0; d < 2 * nslot; d++) dtab[d] = malloc(M * sizeof(int));
    pw[0] = 1; for (int i = 1; i < 40; i++) pw[i] = pw[i - 1] * (deg - 1);
    if (outfn) out = fopen(outfn, "a");
    seed_rng(seed);
    if (ylistfn) { if (k != 1) { fprintf(stderr, "--ylist needs k = 1\n"); exit(1); } load_ylist(ylistfn); }
    fprintf(stderr, "group %s order %d, base %s, n=%d, free slots=%d, D=%d, deg=%d, mode=%s\n", groupdesc, m, basedesc, n, nfree, D, deg, mode);
    if (evalv) { parse_assignment(evalv); run_eval(); }
    else if (selftest) run_selftest(selftest);
    else if (!strcmp(mode, "exhaustive")) run_exhaustive(prefix);
    else if (!strcmp(mode, "random")) run_random(samples);
    else if (!strcmp(mode, "rsearch")) run_rsearch(satime);
    else if (!strcmp(mode, "sa")) run_sa(satime, T0, T1, steps);
    else if (!strcmp(mode, "ils")) run_ils(satime, kick, patience, sweep);
    else if (!strcmp(mode, "sweep1")) run_sweep1();
    else if (!strcmp(mode, "pair")) run_pair(targetJ);
    else { fprintf(stderr, "bad mode\n"); exit(1); }
    if (out) fclose(out);
    return 0;
}
