/* tablefree_sweep.c: table-free evaluator and sweeper for G = K x| C_M, K = S x Z_q, where S is the direct product of the non-cyclic
 * blocks of a gf.py spec (F_p^e, F_p^2 x| C_k, p^{1+2}_+, Z_{p^a} x Z_{p^b}; tables of the blocks come from the TFG1 file
 * written by tfgroup.py) and Z_q is the CRT fusion of its cyclic blocks (alpha acts on it as x -> u x).  Elements
 * (s, x, m) = (s, x) y^m; (s1,x1,m1)(s2,x2,m2) = (s1 alpha_S^m1(s2), x1 + u^m1 x2 mod q, m1 + m2 mod M).  I/O uses gf.py's
 * indices (k M + m, k the mixed-radix block index in spec order), so sets agree with table_sweep.c where both apply.
 *
 * BFS: bit-parallel over rows (fibre, s, m) of q bits (bit x): right multiplication by a fixed element maps a row to a row
 * with a rotation by u^m x_g, so one level of BFS is (#active rows) x (#arcs) rotate-or operations on ceil(q/64) words.
 * Families:  lift = shape-103 two-vertex lift (loop a at 0, loop b at 1, links e, c2, c3 from 0 to 1; n = 2|G|, degree 5;
 *            two BFS, from (0, e) and (1, e), decide the diameter);  cay = Cay(G, {a, b, c}^+-) (degree 6; optional
 *            involution t = y^(M/2) with --tinv for degree 7).
 * Sweeps (complete, chunked by candidate index): lift in a level pattern (ma, mb, m2, m3) up to H = Inn(G) o B (B = the
 * block automorphisms of S commuting with alpha_S, times the units of Z_q), cay over all 3-pair sets up to H x inversions x
 * permutations exactly as table_sweep.c.  --mode cover brute-forces completeness on small groups.
 *
 * tf GROUP.tfg --family lift|cay --D D --mode info|eval|random|patscreen|exhaustive|cover|edges|selftest
 *    [--set a b c2 c3 | a b c] [--pat ma mb m2 m3] [--patfile F] [--samples K] [--seed s] [--chunk i N] [--near T]
 *    [--time S] [--tag TAG] [--out results.jsonl] [--dry] [--edges FILE] [--tinv] [--stabcap C] [--unitcap U]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static double now(void) { struct timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + 1e-6 * tv.tv_usec; }
static void die(const char *m) { fprintf(stderr, "FATAL: %s\n", m); exit(1); }

/* ---------------- group ---------------- */
typedef struct { int type, size, ngens, naut, q, u; uint16_t *T, *alpha; int *gens; uint16_t *aut; } Block;
static int nblk; static Block *blk; static char *spec;
static int M, nS = 1, q = 1, nK, n;             /* n = |G| */
static long long u = 1, *upow;                   /* upow[m] = u^m mod q */
static uint16_t *TS, *alS, *apS, *invS;          /* S table, alpha_S, apS[m*nS+s] = alpha_S^m(s), inverses */
static int *Sgens, nSgens;                       /* S indices of the block generators */
static int *autS, nautS = 1;                     /* automorphisms of S commuting with alpha_S (nautS x nS), [0] = id */
static int *radS; static int *sblk, nsblk;       /* non-cyc blocks: indices and S radices */
static int *cblk, ncblk; static long long *crt;  /* cyc blocks and CRT coefficients */
static long long *radAll;                        /* gf radices of all blocks */
static int *units, nunits;                       /* a generating set of Z_q^* (small) */
static int *allunits, nallunits;                 /* all units (for B lists) */
static int D = 9, family = 0;                    /* 1 = lift, 2 = cay */
static int use_tinv = 0, tinv = -1;
static const char *gfile;

static inline int es(int t) { return t / (q * M); }
static inline int ex(int t) { return (t / M) % q; }
static inline int em(int t) { return t % M; }
static inline int mk(int s, int x, int m) { return (s * q + x) * M + m; }
static inline int mul(int a, int b) {
    int s1 = a / (q * M), r1 = a % (q * M), x1 = r1 / M, m1 = r1 % M;
    int s2 = b / (q * M), r2 = b % (q * M), x2 = r2 / M, m2 = r2 % M;
    int s = TS[(size_t)s1 * nS + apS[(size_t)m1 * nS + s2]];
    int x = (int)((x1 + upow[m1] * x2) % q);
    int m = m1 + m2; if (m >= M) m -= M;
    return (s * q + x) * M + m;
}
static inline int inv(int a) {
    int s = a / (q * M), r = a % (q * M), x = r / M, m = r % M;
    int mi = (M - m) % M;
    int si = apS[(size_t)mi * nS + invS[s]];
    int xi = (int)((q - (upow[mi] * x) % q) % q);
    return (si * q + xi) * M + mi;
}
static inline int cnj(int c, int g) { return mul(mul(c, g), inv(c)); }
static inline int bapply(int j, int v, int g) {   /* beta_(sigma_j, v): (s, x, m) -> (sigma_j(s), v x, m) */
    int s = g / (q * M), r = g % (q * M), x = r / M, m = r % M;
    return (autS[(size_t)j * nS + s] * q + (int)(((long long)v * x) % q)) * M + m;
}
static inline int kpart(int g) { return g / M; }        /* internal K index s*q + x */

/* gf.py index <-> internal */
static int to_gf(int t) {
    int s = es(t), x = ex(t), m = em(t);
    long long k = 0;
    int si = 0;
    for (int b = 0; b < nblk; b++) {
        int d;
        if (blk[b].type == 0) d = x % blk[b].q;
        else { d = (s / radS[si]) % blk[b].size; si++; }
        k = k * blk[b].size + d;
    }
    return (int)(k * M + m);
}
static int from_gf(long long idx) {
    int m = (int)(idx % M); long long k = idx / M;
    int digs[16];
    for (int b = nblk - 1; b >= 0; b--) { digs[b] = (int)(k % blk[b].size); k /= blk[b].size; }
    int s = 0, si = 0; long long x = 0;
    for (int b = 0; b < nblk; b++) {
        if (blk[b].type == 0) { x = (x + (long long)digs[b] * crt[b]) % q; }
        else { s += digs[b] * radS[si]; si++; }
    }
    return mk(s, (int)x, m);
}

static long long egcd_inv(long long a, long long m) {   /* a^-1 mod m */
    long long g = m, x0 = 0, x1 = 1, aa = a % m;
    while (aa) { long long qq = g / aa, t = g - qq * aa; g = aa; aa = t; t = x0 - qq * x1; x0 = x1; x1 = t; }
    if (g != 1) die("CRT: moduli not coprime");
    return ((x0 % m) + m) % m;
}

static void load(const char *fn) {
    FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
    char magic[4]; int32_t hdr[2];
    if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "TFG1", 4)) die("bad magic");
    if (fread(hdr, 4, 2, f) != 2) die("bad header");
    M = hdr[0]; nblk = hdr[1]; if (nblk > 16) die("too many blocks");
    blk = calloc(nblk, sizeof(Block));
    for (int b = 0; b < nblk; b++) {
        int32_t h[4]; if (fread(h, 4, 4, f) != 4) die("short block header");
        Block *B = &blk[b]; B->type = h[0]; B->size = h[1]; B->ngens = h[2]; B->naut = h[3];
        if (B->type == 0) {
            int32_t qu[2]; if (fread(qu, 4, 2, f) != 2) die("short cyc"); B->q = qu[0]; B->u = qu[1];
            if (B->size != B->q) die("cyc size");
        } else {
            size_t sz = (size_t)B->size * B->size;
            B->T = malloc(sz * 2); B->alpha = malloc((size_t)B->size * 2); B->gens = malloc(B->ngens * sizeof(int));
            B->aut = malloc((size_t)B->naut * B->size * 2);
            if (fread(B->T, 2, sz, f) != sz || fread(B->alpha, 2, B->size, f) != (size_t)B->size) die("short table");
            if (fread(B->gens, 4, B->ngens, f) != (size_t)B->ngens) die("short gens");
            if (fread(B->aut, 2, (size_t)B->naut * B->size, f) != (size_t)B->naut * B->size) die("short aut");
        }
    }
    spec = calloc(1 << 16, 1); size_t got = fread(spec, 1, (1 << 16) - 1, f); spec[got] = 0;
    fclose(f);
    /* blocks: S = non-cyc (mixed radix, first most significant), Z_q = CRT of the cyc blocks */
    sblk = malloc(nblk * sizeof(int)); cblk = malloc(nblk * sizeof(int)); radS = malloc(nblk * sizeof(int));
    crt = calloc(nblk, sizeof(long long)); radAll = malloc(nblk * sizeof(long long));
    nsblk = ncblk = 0; long long qq = 1;
    for (int b = 0; b < nblk; b++) { if (blk[b].type == 0) { cblk[ncblk++] = b; qq *= blk[b].q; } else sblk[nsblk++] = b; }
    if (qq > 2000000) die("q too large");
    q = (int)qq;
    nS = 1; for (int i = nsblk - 1; i >= 0; i--) { radS[i] = nS; nS *= blk[sblk[i]].size; }
    if ((long long)nS * q * M >= (1LL << 31)) die("|G| too large");
    nK = nS * q; n = nK * M;
    long long r = 1; for (int b = nblk - 1; b >= 0; b--) { radAll[b] = r; r *= blk[b].size; }
    /* CRT coefficients and u */
    u = 0;
    for (int i = 0; i < ncblk; i++) {
        int b = cblk[i]; long long qi = blk[b].q, Mi = q / qi;
        crt[b] = qi == 1 ? 0 : (Mi * egcd_inv(Mi % qi, qi)) % q;
        u = (u + (long long)blk[b].u * crt[b]) % q;
    }
    if (q == 1) u = 0;
    upow = malloc(M * sizeof(long long)); upow[0] = 1 % q;
    for (int m = 1; m < M; m++) upow[m] = (upow[m - 1] * u) % q;
    if ((upow[M - 1] * u) % q != 1 % q) die("u^M != 1 mod q");
    /* S table and alpha_S */
    TS = malloc((size_t)nS * nS * 2); alS = malloc((size_t)nS * 2); invS = malloc((size_t)nS * 2);
    int *dig = malloc(nsblk * sizeof(int)), *dig2 = malloc(nsblk * sizeof(int));
    for (int s1 = 0; s1 < nS; s1++) {
        for (int i = 0; i < nsblk; i++) dig[i] = (s1 / radS[i]) % blk[sblk[i]].size;
        int al = 0; for (int i = 0; i < nsblk; i++) al += blk[sblk[i]].alpha[dig[i]] * radS[i];
        alS[s1] = (uint16_t)al;
        for (int s2 = 0; s2 < nS; s2++) {
            for (int i = 0; i < nsblk; i++) dig2[i] = (s2 / radS[i]) % blk[sblk[i]].size;
            int pr = 0; for (int i = 0; i < nsblk; i++) { Block *B = &blk[sblk[i]]; pr += B->T[(size_t)dig[i] * B->size + dig2[i]] * radS[i]; }
            TS[(size_t)s1 * nS + s2] = (uint16_t)pr;
            if (pr == 0) invS[s1] = (uint16_t)s2;
        }
    }
    apS = malloc((size_t)M * nS * 2);
    for (int s = 0; s < nS; s++) apS[s] = (uint16_t)s;
    for (int m = 1; m < M; m++) for (int s = 0; s < nS; s++) apS[(size_t)m * nS + s] = alS[apS[(size_t)(m - 1) * nS + s]];
    for (int s = 0; s < nS; s++) if (alS[apS[(size_t)(M - 1) * nS + s]] != s) die("alpha_S^M != id");
    /* S generators */
    Sgens = malloc(64 * sizeof(int)); nSgens = 0;
    for (int i = 0; i < nsblk; i++) for (int g = 0; g < blk[sblk[i]].ngens; g++) Sgens[nSgens++] = blk[sblk[i]].gens[g] * radS[i];
    /* automorphisms of S: product of the block lists while <= 4000 */
    long long tot = 1; int *use = calloc(nsblk, sizeof(int));
    for (int i = 0; i < nsblk; i++) { if (tot * blk[sblk[i]].naut <= 4000) { use[i] = 1; tot *= blk[sblk[i]].naut; } }
    nautS = (int)tot; autS = malloc((size_t)nautS * nS * sizeof(int));
    for (int j = 0; j < nautS; j++) {
        int rem = j; int ai[16];
        for (int i = nsblk - 1; i >= 0; i--) { int cnt = use[i] ? blk[sblk[i]].naut : 1; ai[i] = rem % cnt; rem /= cnt; }
        for (int s = 0; s < nS; s++) {
            int img = 0;
            for (int i = 0; i < nsblk; i++) { Block *B = &blk[sblk[i]]; int d = (s / radS[i]) % B->size; img += (use[i] ? B->aut[(size_t)ai[i] * B->size + d] : d) * radS[i]; }
            autS[(size_t)j * nS + s] = img;
        }
    }
    free(use);
    /* units of Z_q: all, and a small generating set */
    allunits = malloc(q * sizeof(int)); nallunits = 0;
    for (int v = 1; v < q; v++) { long long a = v, b = q; while (b) { long long t = a % b; a = b; b = t; } if (a == 1) allunits[nallunits++] = v; }
    if (q == 1) { allunits[nallunits++] = 0; }
    units = malloc(32 * sizeof(int)); nunits = 0;
    if (q > 1) {
        char *gen = calloc(q, 1); int *sub = malloc(q * sizeof(int)); int nsub = 1; sub[0] = 1 % q; gen[1 % q] = 1;
        uint64_t rs = 0x1234567ULL;
        for (int tries = 0; nsub < nallunits && tries < 100000; tries++) {
            rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
            int v = allunits[rs % nallunits];
            if (gen[v]) continue;
            units[nunits++] = v;
            /* close the subgroup */
            for (int i = 0; i < nsub; i++) { long long w = sub[i]; for (;;) { w = (w * v) % q; if (gen[w]) break; gen[w] = 1; sub[nsub++] = (int)w; } }
            for (int i = 0; i < nsub; i++) for (int j = 0; j < nunits; j++) { long long w = ((long long)sub[i] * units[j]) % q; if (!gen[w]) { gen[w] = 1; sub[nsub++] = (int)w; } }
            if (nunits >= 31) break;
        }
        if (nsub != nallunits) die("unit generators");
        free(gen); free(sub);
    }
    if (use_tinv) { if (M % 2) die("tinv needs M even"); tinv = mk(0, 0, M / 2); if (inv(tinv) != tinv) die("tinv not an involution"); }
    /* sanity: associativity sample, identity, inverses */
    uint64_t x = 88172645463325252ULL;
    for (int i = 0; i < 300000; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        int a = x % n, b = (x >> 20) % n, c = (x >> 41) % n;
        if (mul(mul(a, b), c) != mul(a, mul(b, c))) die("associativity");
        if (mul(a, inv(a)) != 0 || mul(0, a) != a || mul(a, 0) != a) die("inverse/identity");
        if (from_gf(to_gf(a)) != a) die("index map");
        for (int j = 1; j < nautS && j < 4; j++) if (bapply(j, 1 % q, mul(a, b)) != mul(bapply(j, 1 % q, a), bapply(j, 1 % q, b))) die("autS not an automorphism");
        if (nunits) { int v = units[i % nunits]; if (bapply(0, v, mul(a, b)) != mul(bapply(0, v, a), bapply(0, v, b))) die("unit aut"); }
    }
    free(dig); free(dig2);
}

/* ---------------- invariants: Z(G), involutions, k(G) by union-find over the generators' conjugations ---------------- */
static int *uf;
static int ufind(int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }
static void uunion(int a, int b) { a = ufind(a); b = ufind(b); if (a != b) { if (a < b) uf[b] = a; else uf[a] = b; } }
static int Ggens[80], nGgens;
static void group_gens(void) {
    nGgens = 0;
    for (int i = 0; i < nSgens; i++) Ggens[nGgens++] = mk(Sgens[i], 0, 0);
    if (q > 1) Ggens[nGgens++] = mk(0, 1, 0);
    if (M > 1) Ggens[nGgens++] = mk(0, 0, 1);
}
static int nclasses, Zsize, ninvol;
static void invariants(void) {
    group_gens();
    uf = malloc(n * sizeof(int)); for (int g = 0; g < n; g++) uf[g] = g;
    for (int i = 0; i < nGgens; i++) { int c = Ggens[i]; for (int g = 0; g < n; g++) uunion(g, cnj(c, g)); }
    nclasses = 0; for (int g = 0; g < n; g++) if (ufind(g) == g) nclasses++;
    Zsize = 0; for (int g = 0; g < n; g++) { int z = 1; for (int i = 0; i < nGgens && z; i++) if (cnj(Ggens[i], g) != g) z = 0; Zsize += z; }
    ninvol = 0; for (int g = 1; g < n; g++) if (mul(g, g) == 0) ninvol++;
    free(uf);
}

/* ---------------- bit-parallel BFS ---------------- */
static int W, R, F;                 /* words per row, rows, fibres */
static uint64_t topmask;
static uint64_t *vis, *fr, *nx;
static int *tgt, *rot, deg;         /* per row: deg slots */
static int *active, *touched, nactive; static char *tflag;
static long bfs_calls = 0;

static inline void rot_or(uint64_t *dst, const uint64_t *src, int t) {
    if (t == 0) { for (int j = 0; j < W; j++) dst[j] |= src[j]; return; }
    if (W == 1) { dst[0] |= ((src[0] << t) | (src[0] >> (q - t))) & topmask; return; }
    int ws = t >> 6, bs = t & 63;
    for (int j = W - 1; j >= ws; j--) {
        uint64_t v = src[j - ws] << bs;
        if (bs && j - ws - 1 >= 0) v |= src[j - ws - 1] >> (64 - bs);
        if (j == W - 1) v &= topmask;
        dst[j] |= v;
    }
    int t2 = q - t, ws2 = t2 >> 6, bs2 = t2 & 63;
    for (int j = 0; j + ws2 < W; j++) {
        uint64_t v = src[j + ws2] >> bs2;
        if (bs2 && j + ws2 + 1 < W) v |= src[j + ws2 + 1] << (64 - bs2);
        dst[j] |= v;
    }
}

static void bfs_alloc(void) {
    W = (q + 63) / 64; topmask = (q % 64) ? ((1ULL << (q % 64)) - 1) : ~0ULL;
    F = family == 1 ? 2 : 1; R = F * nS * M;
    deg = family == 1 ? 5 : (6 + use_tinv);
    vis = malloc((size_t)R * W * 8); fr = malloc((size_t)R * W * 8); nx = malloc((size_t)R * W * 8);
    tgt = malloc((size_t)R * deg * sizeof(int)); rot = malloc((size_t)R * deg * sizeof(int));
    active = malloc(R * sizeof(int)); touched = malloc(R * sizeof(int)); tflag = calloc(R, 1);
}
/* slot: row (f, s, m) times g into fibre f2 */
static inline void set_slot(int f, int s, int m, int k, int g, int f2) {
    int gs = es(g), gx = ex(g), gm = em(g);
    int s2 = TS[(size_t)s * nS + apS[(size_t)m * nS + gs]];
    int m2 = m + gm; if (m2 >= M) m2 -= M;
    int row = (f * nS + s) * M + m;
    tgt[(size_t)row * deg + k] = (f2 * nS + s2) * M + m2;
    rot[(size_t)row * deg + k] = (int)((upow[m] * gx) % q);
}
static int la, lb, lc2, lc3;        /* lift voltages (internal indices) */
static int ga, gb, gc;              /* cay generators */
static void setup_lift(void) {
    int ia = inv(la), ib = inv(lb), i2 = inv(lc2), i3 = inv(lc3);
    for (int s = 0; s < nS; s++) for (int m = 0; m < M; m++) {
        set_slot(0, s, m, 0, la, 0); set_slot(0, s, m, 1, ia, 0); set_slot(0, s, m, 2, 0, 1); set_slot(0, s, m, 3, lc2, 1); set_slot(0, s, m, 4, lc3, 1);
        set_slot(1, s, m, 0, lb, 1); set_slot(1, s, m, 1, ib, 1); set_slot(1, s, m, 2, 0, 0); set_slot(1, s, m, 3, i2, 0); set_slot(1, s, m, 4, i3, 0);
    }
}
static void setup_cay(void) {
    for (int s = 0; s < nS; s++) for (int m = 0; m < M; m++) {
        set_slot(0, s, m, 0, ga, 0); set_slot(0, s, m, 1, inv(ga), 0); set_slot(0, s, m, 2, gb, 0);
        set_slot(0, s, m, 3, inv(gb), 0); set_slot(0, s, m, 4, gc, 0); set_slot(0, s, m, 5, inv(gc), 0);
        if (use_tinv) set_slot(0, s, m, 6, tinv, 0);
    }
}
/* BFS from (root row, bit 0); returns #vertices beyond D (or a lower bound > bound_T when bound_T >= 0 and the defect
 * provably exceeds it); prof[0..D] optional */
static int bfs(int root_row, int bound_T, int *prof) {
    bfs_calls++;
    int N = R * q;                     /* vertices */
    memset(vis, 0, (size_t)R * W * 8); memset(fr, 0, (size_t)R * W * 8); memset(nx, 0, (size_t)R * W * 8);
    vis[(size_t)root_row * W] = 1; fr[(size_t)root_row * W] = 1;
    nactive = 1; active[0] = root_row;
    long reached = 1, fcount = 1;
    if (prof) prof[0] = 1;
    int out = deg - 1;
    for (int lev = 1; lev <= D; lev++) {
        if (bound_T >= 0) {
            /* upper bound on what the remaining levels can add: fcount*(d or d-1) * (1 + (d-1) + ... ) */
            long mx = fcount * (lev == 1 ? deg : out), term = mx;
            for (int l = lev + 1; l <= D; l++) { term *= out; mx += term; if (mx > N) break; }
            if (N - reached - mx > bound_T) return (int)(N - reached - mx);
        }
        int ntouched = 0;
        for (int i = 0; i < nactive; i++) {
            int r = active[i]; const uint64_t *src = fr + (size_t)r * W;
            const int *tg = tgt + (size_t)r * deg, *ro = rot + (size_t)r * deg;
            for (int k = 0; k < deg; k++) {
                int t = tg[k];
                rot_or(nx + (size_t)t * W, src, ro[k]);
                if (!tflag[t]) { tflag[t] = 1; touched[ntouched++] = t; }
            }
        }
        for (int i = 0; i < nactive; i++) memset(fr + (size_t)active[i] * W, 0, W * 8);
        nactive = 0; fcount = 0;
        for (int i = 0; i < ntouched; i++) {
            int r = touched[i]; tflag[r] = 0;
            uint64_t *vr = vis + (size_t)r * W, *nr = nx + (size_t)r * W, *fr_ = fr + (size_t)r * W;
            long cnt = 0;
            for (int j = 0; j < W; j++) { uint64_t nw = nr[j] & ~vr[j]; vr[j] |= nw; fr_[j] = nw; cnt += __builtin_popcountll(nw); nr[j] = 0; }
            if (cnt) { active[nactive++] = r; fcount += cnt; }
        }
        reached += fcount;
        if (prof) prof[lev] = (int)fcount;
        if (!nactive) break;
    }
    return (int)(N - reached);
}
static int lift_defect(int bound_T, int *p0, int *p1) {
    setup_lift();
    int d0 = bfs(0, bound_T, p0);
    if (bound_T >= 0 && d0 > bound_T) return d0;
    int d1 = bfs(nS * M, bound_T >= 0 ? bound_T - d0 : -1, p1);   /* root (1, e): row (1*nS + 0)*M + 0 */
    return d0 + d1;
}
static int cay_defect(int bound_T, int *prof) { setup_cay(); return bfs(0, bound_T, prof); }

/* plain BFS with mul() for the self-test */
static int plain_defect_lift(int root, int *prof) {
    int N2 = 2 * n; uint8_t *seen = calloc(N2, 1); int *qu = malloc(N2 * sizeof(int));
    int st = root * n; seen[st] = 1; qu[0] = st; int qh = 0, qt = 1, reached = 1;
    int ia = inv(la), ib = inv(lb), i2 = inv(lc2), i3 = inv(lc3);
    prof[0] = 1;
    for (int lev = 1; lev <= D; lev++) {
        int end = qt;
        while (qh < end) {
            int v = qu[qh++], g = v % n, w[5];
            if (v < n) { w[0] = mul(g, la); w[1] = mul(g, ia); w[2] = n + g; w[3] = n + mul(g, lc2); w[4] = n + mul(g, lc3); }
            else { w[0] = n + mul(g, lb); w[1] = n + mul(g, ib); w[2] = g; w[3] = mul(g, i2); w[4] = mul(g, i3); }
            for (int k = 0; k < 5; k++) if (!seen[w[k]]) { seen[w[k]] = 1; qu[qt++] = w[k]; }
        }
        reached = qt; prof[lev] = qt - end;
        if (qh == qt) break;
    }
    free(seen); free(qu);
    return N2 - reached;
}
static int plain_defect_cay(int *prof) {
    uint8_t *seen = calloc(n, 1); int *qu = malloc(n * sizeof(int));
    seen[0] = 1; qu[0] = 0; int qh = 0, qt = 1, reached = 1;
    int gens[7] = {ga, inv(ga), gb, inv(gb), gc, inv(gc), tinv}; int ng = 6 + use_tinv;
    prof[0] = 1;
    for (int lev = 1; lev <= D; lev++) {
        int end = qt;
        while (qh < end) { int v = qu[qh++]; for (int k = 0; k < ng; k++) { int w = mul(v, gens[k]); if (!seen[w]) { seen[w] = 1; qu[qt++] = w; } } }
        reached = qt; prof[lev] = qt - end;
        if (qh == qt) break;
    }
    free(seen); free(qu);
    return n - reached;
}

/* ---------------- RNG ---------------- */
static uint64_t rs[2] = {0x9E3779B97F4A7C15ULL, 0xD1B54A32D192ED03ULL};
static inline uint64_t rotl(uint64_t x, int r) { return (x << r) | (x >> (64 - r)); }
static uint64_t rnd(void) { uint64_t s0 = rs[0], s1 = rs[1], r = s0 + s1; s1 ^= s0; rs[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); rs[1] = rotl(s1, 36); return r; }
static void seed_rng(uint64_t seed) { rs[0] = seed * 0x9E3779B97F4A7C15ULL + 1; rs[1] = (seed ^ 0xDEADBEEFCAFEULL) * 0xBF58476D1CE4E5B9ULL + 7; for (int i = 0; i < 20; i++) rnd(); }
static inline int is_pair(int g) { return g != 0 && inv(g) != g; }

/* ---------------- reporting ---------------- */
static const char *tag = ""; static FILE *out = NULL;
static int near_T = 0; static double t_limit = 1e18, t_start;
static long ncand = 0, nmine = 0, nfound = 0; static int best_defect = 1 << 30, best_set[4] = {0, 0, 0, 0};
static int chunk_i = 0, chunk_N = 1, stopped = 0;
static int pat[4] = {-1, -1, -1, -1};

static void report(int *s, int d) {   /* s: internal indices; printed in gf indices */
    int g[4] = {to_gf(s[0]), to_gf(s[1]), family == 1 ? to_gf(s[2]) : (s[2] ? to_gf(s[2]) : 0), family == 1 ? to_gf(s[3]) : 0};
    if (d == 0) {
        nfound++;
        if (family == 1) printf("FOUND diam<=%d (defect 0) family=lift group=%s set=[%d,%d,%d,%d]\n", D, gfile, g[0], g[1], g[2], g[3]);
        else printf("FOUND diam<=%d (defect 0) family=cay group=%s set=[%d,%d,%d]\n", D, gfile, g[0], g[1], g[2]);
        fflush(stdout);
        if (out) { fprintf(out, "{\"event\":\"found\",\"tag\":\"%s\",\"family\":%d,\"D\":%d,\"n\":%d,\"group\":\"%s\",\"set\":[%d,%d,%d,%d]}\n", tag, family, D, family == 1 ? 2 * n : n, gfile, g[0], g[1], g[2], g[3]); fflush(out); }
    } else if (d <= near_T) {
        if (family == 1) printf("NEAR defect=%d set=[%d,%d,%d,%d]\n", d, g[0], g[1], g[2], g[3]);
        else printf("NEAR defect=%d set=[%d,%d,%d]\n", d, g[0], g[1], g[2]);
        fflush(stdout);
    }
    if (d < best_defect) { best_defect = d; memcpy(best_set, g, 4 * sizeof(int)); }
}

/* ---------------- structure for the sweeps ---------------- */
static int stabcap = 4000, unitcap = 2000;
/* B list for stabilisers: pairs (sigma j, unit v) over all sigma and a deterministic subset of the units */
static int *Blist_j, *Blist_v, nBlist;
static void build_Blist(void) {
    int nv = nallunits < unitcap ? nallunits : unitcap;
    nBlist = nautS * nv;
    Blist_j = malloc(nBlist * sizeof(int)); Blist_v = malloc(nBlist * sizeof(int));
    uint64_t x = 0xABCDEF123ULL; int k = 0;
    /* units: all when few, else a deterministic pseudo-random subset including 1 */
    int *vs = malloc(nv * sizeof(int));
    if (nv == nallunits) memcpy(vs, allunits, nv * sizeof(int));
    else { vs[0] = 1; for (int i = 1; i < nv; i++) { x ^= x << 13; x ^= x >> 7; x ^= x << 17; vs[i] = allunits[x % nallunits]; } }
    for (int j = 0; j < nautS; j++) for (int i = 0; i < nv; i++) { Blist_j[k] = j; Blist_v[k] = vs[i]; k++; }
    free(vs);
}

/* orbits of H (conj by all of G, all beta) on the K-coset of layer m (K indices 0..nK-1), merging inverses if mergeinv.
 * uf over nK.  Returns the number of orbits; omin[] gets the smallest element (K index) of each orbit; orb[k] = orbit idx */
static int coset_orbits(int m, int mergeinv, int *orb, int *omin) {
    for (int k = 0; k < nK; k++) uf[k] = k;
    for (int i = 0; i < nGgens; i++) { int c = Ggens[i]; for (int k = 0; k < nK; k++) uunion(k, kpart(cnj(c, k * M + m))); }
    for (int j = 1; j < nautS; j++) for (int k = 0; k < nK; k++) uunion(k, kpart(bapply(j, 1 % q, k * M + m)));
    for (int i = 0; i < nunits; i++) for (int k = 0; k < nK; k++) uunion(k, kpart(bapply(0, units[i], k * M + m)));
    if (mergeinv) for (int k = 0; k < nK; k++) uunion(k, kpart(inv(k * M + m)));
    int cnt = 0;
    for (int k = 0; k < nK; k++) if (ufind(k) == k) omin[cnt++] = k;
    int *idx = malloc(nK * sizeof(int));
    for (int i = 0; i < cnt; i++) idx[omin[i]] = i;
    for (int k = 0; k < nK; k++) orb[k] = idx[ufind(k)];
    free(idx);
    return cnt;
}
/* explicit stabiliser of a in H as pairs (c, beta=(j, v)): c beta(a) c^-1 = a; capped at stabcap elements (deterministic) */
static int *stC, *stJ, *stV, nst;
static int *cconj;    /* over K indices of a's coset */
static void stabiliser(int a) {
    int ma = em(a);
    for (int k = 0; k < nK; k++) cconj[k] = -1;
    for (int c = 0; c < n; c++) { int g = cnj(c, a); if (em(g) != ma) die("conj changes layer"); if (cconj[kpart(g)] < 0) cconj[kpart(g)] = c; }
    nst = 0;
    /* centraliser elements (beta = id) */
    for (int c = 0; c < n && nst < stabcap; c++) if (cnj(c, a) == a && c != 0) { stC[nst] = c; stJ[nst] = 0; stV[nst] = 1 % q; nst++; }
    int ncent = nst;
    /* beta's whose image is conjugate to a: (c0, beta) with c0 = cconj[beta(a)]^-1 ... c0 beta(a) c0^-1 = a */
    for (int i = 0; i < nBlist && nst < stabcap; i++) {
        int j = Blist_j[i], v = Blist_v[i];
        if (j == 0 && v == 1 % q) continue;
        int x = bapply(j, v, a); int c1 = cconj[kpart(x)];
        if (c1 < 0) continue;                 /* beta(a) not conjugate to a */
        /* c1 a c1^-1 = beta(a)  =>  c1^-1 beta(a) c1 = a: c0 = c1^-1 */
        int c0 = inv(c1);
        stC[nst] = c0; stJ[nst] = j; stV[nst] = v; nst++;
    }
    (void)ncent;
}
/* orbits on the K-coset of layer m under the stabiliser pairs (restricted to those fixing extra point p2 if p2 >= 0) */
static int stab_orbits(int m, int p2, int *orb, int *omin, long *nused) {
    long used = 0;
    for (int i = 0; i < nst; i++) if (p2 < 0 || cnj(stC[i], bapply(stJ[i], stV[i], p2)) == p2) used++;
    *nused = used;
    if (!used) { for (int k = 0; k < nK; k++) { orb[k] = k; omin[k] = k; } return nK; }
    for (int k = 0; k < nK; k++) uf[k] = k;
    for (int i = 0; i < nst; i++) {
        int c = stC[i], j = stJ[i], v = stV[i];
        if (p2 >= 0 && cnj(c, bapply(j, v, p2)) != p2) continue;
        for (int k = 0; k < nK; k++) uunion(k, kpart(cnj(c, bapply(j, v, k * M + m))));
    }
    int cnt = 0;
    for (int k = 0; k < nK; k++) if (ufind(k) == k) omin[cnt++] = k;
    int *idx = malloc(nK * sizeof(int));
    for (int i = 0; i < cnt; i++) idx[omin[i]] = i;
    for (int k = 0; k < nK; k++) orb[k] = idx[ufind(k)];
    free(idx);
    return cnt;
}

static void consider_lift(int a, int b, int c2, int c3) {
    long idx = ncand++;
    if (idx % chunk_N != chunk_i) return;
    nmine++;
    la = a; lb = b; lc2 = c2; lc3 = c3;
    int T = near_T > best_defect ? near_T : best_defect;
    int d = lift_defect(T, NULL, NULL);
    if (d <= T) { int s[4] = {a, b, c2, c3}; report(s, d); }
}
static void consider_cay(int a, int b, int c) {
    long idx = ncand++;
    if (idx % chunk_N != chunk_i) return;
    nmine++;
    ga = a; gb = b; gc = c;
    int T = near_T > best_defect ? near_T : best_defect;
    int d = cay_defect(T, NULL);
    if (d <= T) { int s[4] = {a, b, c, 0}; report(s, d); }
}

/* cover test bitsets (small groups) */
static uint8_t *covered; static int *pairidx, npair;
static inline void setbit(size_t bit) { covered[bit >> 3] |= 1 << (bit & 7); }
static inline int getbit(size_t bit) { return (covered[bit >> 3] >> (bit & 7)) & 1; }
static void mark_lift_full(int a, int b, int c2, int c3, int ma, int mb) {
    /* orbit under all (c, beta) in G x B (full B: all sigma x all units), a<->a^-1 if allowed, b<->b^-1 if allowed */
    for (int j = 0; j < nautS; j++) for (int vi = 0; vi < nallunits; vi++) for (int c = 0; c < n; c++) {
        int v = allunits[vi];
        int A = cnj(c, bapply(j, v, a)), Bb = cnj(c, bapply(j, v, b)), C2 = cnj(c, bapply(j, v, c2)), C3 = cnj(c, bapply(j, v, c3));
        int As[2] = {A, inv(A)}, Bs[2] = {Bb, inv(Bb)};
        for (int ia = 0; ia < (2 * ma % M == 0 ? 2 : 1); ia++) for (int ib = 0; ib < (2 * mb % M == 0 ? 2 : 1); ib++)
            setbit((((size_t)kpart(As[ia]) * nK + kpart(Bs[ib])) * nK + kpart(C2)) * nK + kpart(C3));
    }
}
static void mark_cay_full(int a, int b, int c) {
    for (int j = 0; j < nautS; j++) for (int vi = 0; vi < nallunits; vi++) for (int x = 0; x < n; x++) {
        int v = allunits[vi];
        int p[3] = {pairidx[cnj(x, bapply(j, v, a))], pairidx[cnj(x, bapply(j, v, b))], pairidx[cnj(x, bapply(j, v, c))]}; int t;
        if (p[0] > p[1]) { t = p[0]; p[0] = p[1]; p[1] = t; }
        if (p[1] > p[2]) { t = p[1]; p[1] = p[2]; p[2] = t; }
        if (p[0] > p[1]) { t = p[0]; p[0] = p[1]; p[1] = t; }
        setbit(((size_t)p[0] * npair + p[1]) * npair + p[2]);
    }
}

/* lift sweep over the pattern (ma, mb, m2, m3): a = H-orbit reps on coset ma; c2 = Stab(a)-orbit reps on coset m2;
 * c3 = Stab(a, c2)-orbit reps on coset m3; b = all of coset mb (b <= b^-1 when 2 mb = 0) */
static void run_lift_pattern(int dry) {
    int ma = pat[0], mb = pat[1], m2 = pat[2], m3 = pat[3];
    if (ma < 0) die("--pat ma mb m2 m3 required");
    group_gens(); build_Blist();
    uf = malloc(nK * sizeof(int)); cconj = malloc(nK * sizeof(int));
    stC = malloc(stabcap * sizeof(int)); stJ = malloc(stabcap * sizeof(int)); stV = malloc(stabcap * sizeof(int));
    int *orbA = malloc(nK * sizeof(int)), *minA = malloc(nK * sizeof(int));
    int *orb2 = malloc(nK * sizeof(int)), *min2 = malloc(nK * sizeof(int));
    int *orb3 = malloc(nK * sizeof(int)), *min3 = malloc(nK * sizeof(int));
    if (dry == 2) { if (nK > 40) die("cover needs nK <= 40"); size_t bits = (size_t)nK * nK * nK * nK; covered = calloc(bits / 8 + 1, 1); }
    int nA = coset_orbits(ma, (2 * ma) % M == 0, orbA, minA);
    /* b-coset: allowed b (non-identity, non-involution, b <= b^-1 when 2 mb = 0) */
    char *bok = malloc(nK);
    for (int kb = 0; kb < nK; kb++) { int b = kb * M + mb; bok[kb] = !(b == 0 || inv(b) == b || ((2 * mb) % M == 0 && inv(b) < b)); }
    long stab_total = 0, used_total = 0; int nA_used = 0;
    double t_last = now();
    for (int ia = 0; ia < nA && !stopped; ia++) {
        int a = minA[ia] * M + ma;
        if (a == 0 || inv(a) == a) continue;
        nA_used++;
        stabiliser(a); stab_total += nst;
        long used; int n2 = stab_orbits(m2, -1, orb2, min2, &used); used_total += used;
        for (int i2 = 0; i2 < n2 && !stopped; i2++) {
            int c2 = min2[i2] * M + m2;
            if (c2 == 0) continue;
            long used2; int n3 = stab_orbits(m3, c2, orb3, min3, &used2);
            for (int i3 = 0; i3 < n3; i3++) {
                int c3 = min3[i3] * M + m3;
                if (c3 == 0 || c3 == c2) continue;
                for (int kb = 0; kb < nK; kb++) {
                    if (!bok[kb]) continue;
                    int b = kb * M + mb;
                    if (dry) { ncand++; if (dry == 2) mark_lift_full(a, b, c2, c3, ma, mb); continue; }
                    consider_lift(a, b, c2, c3);
                }
            }
            if (now() - t_last > 60) {
                t_last = now();
                printf("PROGRESS a %d/%d c2 %d/%d candidates=%ld mine=%ld bfs=%ld best=%d %.0fs\n", ia, nA, i2, n2, ncand, nmine, bfs_calls, best_defect, now() - t_start); fflush(stdout);
            }
            if (now() - t_start > t_limit) stopped = 1;
        }
    }
    if (dry) {
        printf("COUNT family=lift pat=[%d,%d,%d,%d] n=%d |G|=%d nK=%d |autS|=%d |Blist|=%d aorbits=%d a_used=%d stab=%ld used=%ld candidates(total over all chunks)=%ld time=%.1fs\n",
               ma, mb, m2, m3, 2 * n, n, nK, nautS, nBlist, nA, nA_used, stab_total, used_total, ncand, now() - t_start);
        if (dry == 2) {
            long total = 0, miss = 0;
            for (int ka = 0; ka < nK; ka++) { int a = ka * M + ma; if (a == 0 || inv(a) == a) continue;
                for (int kb = 0; kb < nK; kb++) { int b = kb * M + mb; if (b == 0 || inv(b) == b) continue;
                    for (int k2 = 0; k2 < nK; k2++) { int c2 = k2 * M + m2; if (c2 == 0) continue;
                        for (int k3 = 0; k3 < nK; k3++) { int c3 = k3 * M + m3; if (c3 == 0 || c3 == c2) continue;
                            total++; if (!getbit((((size_t)ka * nK + kb) * nK + k2) * nK + k3)) miss++; } } } }
            printf("COVER pat=[%d,%d,%d,%d] total=%ld uncovered=%ld %s\n", ma, mb, m2, m3, total, miss, miss ? "FAIL" : "OK");
        }
        return;
    }
    double el = now() - t_start;
    printf("EXHAUSTIVE family=lift group=%s D=%d n=%d pat=[%d,%d,%d,%d] chunk=%d/%d candidates=%ld mine=%ld bfs=%ld time=%.1fs rate=%.0f cand/s best_defect=%d best=[%d,%d,%d,%d] found=%ld complete=%d\n",
           gfile, D, 2 * n, ma, mb, m2, m3, chunk_i, chunk_N, ncand, nmine, bfs_calls, el, nmine / (el > 0 ? el : 1), best_defect, best_set[0], best_set[1], best_set[2], best_set[3], nfound, !stopped);
    if (out) { fprintf(out, "{\"event\":\"exhaustive\",\"tag\":\"%s\",\"family\":1,\"D\":%d,\"n\":%d,\"group\":\"%s\",\"pat\":[%d,%d,%d,%d],\"chunk\":[%d,%d],\"candidates\":%ld,\"mine\":%ld,\"time\":%.1f,\"best_defect\":%d,\"best\":[%d,%d,%d,%d],\"found\":%ld,\"complete\":%d}\n",
                          tag, D, 2 * n, gfile, ma, mb, m2, m3, chunk_i, chunk_N, ncand, nmine, el, best_defect, best_set[0], best_set[1], best_set[2], best_set[3], nfound, !stopped); fflush(out); }
}

/* cay sweep (whole family) as table_sweep.c's cay7 with H = Inn(G) o B: types = H-orbits on pairs over G */
static int *type_, *typemin, ntype;
static int number_orbits_G(int *orb, int *omin) {
    int cnt = 0;
    for (int g = 0; g < n; g++) if (ufind(g) == g) omin[cnt++] = g;
    int *idx = malloc(n * sizeof(int));
    for (int i = 0; i < cnt; i++) idx[omin[i]] = i;
    for (int g = 0; g < n; g++) orb[g] = idx[ufind(g)];
    free(idx);
    return cnt;
}
static void run_cay(int dry) {
    group_gens(); build_Blist();
    uf = malloc(n * sizeof(int)); cconj = malloc(nK * sizeof(int));
    stC = malloc(stabcap * sizeof(int)); stJ = malloc(stabcap * sizeof(int)); stV = malloc(stabcap * sizeof(int));
    for (int g = 0; g < n; g++) uf[g] = g;
    for (int g = 0; g < n; g++) uunion(g, inv(g));
    for (int i = 0; i < nGgens; i++) { int c = Ggens[i]; for (int g = 0; g < n; g++) uunion(g, cnj(c, g)); }
    for (int j = 1; j < nautS; j++) for (int g = 0; g < n; g++) uunion(g, bapply(j, 1 % q, g));
    for (int i = 0; i < nunits; i++) for (int g = 0; g < n; g++) uunion(g, bapply(0, units[i], g));
    type_ = malloc(n * sizeof(int)); typemin = malloc(n * sizeof(int));
    ntype = number_orbits_G(type_, typemin);
    int *sorb = malloc(n * sizeof(int)), *somin = malloc(n * sizeof(int));
    if (dry == 2) {
        pairidx = malloc(n * sizeof(int)); npair = 0;
        for (int g = 0; g < n; g++) pairidx[g] = -1;
        for (int g = 0; g < n; g++) if (is_pair(g) && g < inv(g)) { pairidx[g] = pairidx[inv(g)] = npair++; }
        size_t bits = (size_t)npair * npair * npair; covered = calloc(bits / 8 + 1, 1);
    }
    double t_last = now(); long stab_total = 0;
    for (int ta = 0; ta < ntype && !stopped; ta++) {
        int a = typemin[ta];
        if (!is_pair(a)) continue;
        /* stabiliser of the pair {a, a^-1}: pairs (c, beta) with c beta(a) c^-1 in {a, a^-1} */
        int ma = em(a), ia = inv(a);
        for (int k = 0; k < nK; k++) cconj[k] = -1;
        for (int c = 0; c < n; c++) { int g = cnj(c, a); if (cconj[kpart(g)] < 0) cconj[kpart(g)] = c; }
        nst = 0;
        for (int c = 0; c < n && nst < stabcap; c++) { int g = cnj(c, a); if ((g == a || g == ia) && c != 0) { stC[nst] = c; stJ[nst] = 0; stV[nst] = 1 % q; nst++; } }
        for (int i = 0; i < nBlist && nst < stabcap; i++) {
            int j = Blist_j[i], v = Blist_v[i];
            if (j == 0 && v == 1 % q) continue;
            int x = bapply(j, v, a); int c1 = -1;
            if (em(x) == ma) c1 = cconj[kpart(x)];
            if (c1 < 0 && em(inv(x)) == ma) c1 = cconj[kpart(inv(x))];
            if (c1 < 0) continue;
            stC[nst] = inv(c1); stJ[nst] = j; stV[nst] = v; nst++;
        }
        stab_total += nst;
        for (int g = 0; g < n; g++) uf[g] = g;
        for (int g = 0; g < n; g++) uunion(g, inv(g));
        for (int i = 0; i < nst; i++) { int c = stC[i], j = stJ[i], v = stV[i]; for (int g = 0; g < n; g++) uunion(g, cnj(c, bapply(j, v, g))); }
        number_orbits_G(sorb, somin);
        for (int b = 0; b < n && !stopped; b++) {
            if (!is_pair(b) || b == a || b == ia || type_[b] < ta) continue;
            if (somin[sorb[b]] != b) continue;
            int ob = sorb[b], ib = inv(b);
            for (int c = b + 1; c < n; c++) {
                if (!is_pair(c) || c > inv(c) || c == a || c == ia || c == ib) continue;
                if (type_[c] < ta || sorb[c] < ob) continue;
                if (dry) { ncand++; if (dry == 2) mark_cay_full(a, b, c); continue; }
                consider_cay(a, b, c);
            }
            if (now() - t_last > 60) {
                t_last = now();
                printf("PROGRESS type %d/%d b=%d candidates=%ld mine=%ld bfs=%ld best=%d %.0fs\n", ta, ntype, b, ncand, nmine, bfs_calls, best_defect, now() - t_start); fflush(stdout);
            }
            if (now() - t_start > t_limit) stopped = 1;
        }
    }
    if (dry) {
        printf("COUNT family=cay n=%d |autS|=%d |Blist|=%d types=%d stab=%ld candidates(total over all chunks)=%ld time=%.1fs\n", n, nautS, nBlist, ntype, stab_total, ncand, now() - t_start);
        if (dry == 2) {
            long total = 0, miss = 0;
            for (int p0 = 0; p0 < npair; p0++) for (int p1 = p0 + 1; p1 < npair; p1++) for (int p2 = p1 + 1; p2 < npair; p2++) { total++; if (!getbit(((size_t)p0 * npair + p1) * npair + p2)) miss++; }
            printf("COVER total=%ld uncovered=%ld %s\n", total, miss, miss ? "FAIL" : "OK");
        }
        return;
    }
    double el = now() - t_start;
    printf("EXHAUSTIVE family=cay group=%s D=%d n=%d chunk=%d/%d candidates=%ld mine=%ld bfs=%ld time=%.1fs rate=%.0f cand/s best_defect=%d best=[%d,%d,%d,%d] found=%ld complete=%d\n",
           gfile, D, n, chunk_i, chunk_N, ncand, nmine, bfs_calls, el, nmine / (el > 0 ? el : 1), best_defect, best_set[0], best_set[1], best_set[2], best_set[3], nfound, !stopped);
    if (out) { fprintf(out, "{\"event\":\"exhaustive\",\"tag\":\"%s\",\"family\":2,\"D\":%d,\"n\":%d,\"group\":\"%s\",\"chunk\":[%d,%d],\"candidates\":%ld,\"mine\":%ld,\"time\":%.1f,\"best_defect\":%d,\"best\":[%d,%d,%d,%d],\"found\":%ld,\"complete\":%d}\n",
                          tag, D, n, gfile, chunk_i, chunk_N, ncand, nmine, el, best_defect, best_set[0], best_set[1], best_set[2], best_set[3], nfound, !stopped); fflush(out); }
}

/* random screens: whole family, or inside a level pattern */
static int rand_elem_layer(int m, int needpair) {
    for (;;) { int g = (int)(rnd() % nK) * M + m; if (g == 0) continue; if (needpair && inv(g) == g) continue; return g; }
}
static void random_screen(long K, const int *P, const char *label) {
    long hist[40] = {0}; int best = 1 << 30, bs[4] = {0, 0, 0, 0}; double suml = 0; int prof0[64], prof1[64];
    for (long i = 0; i < K; i++) {
        int d;
        if (family == 1) {
            int a, b, c2, c3;
            if (P) {   /* a negative layer means: random layer for that voltage (loop-layer screens) */
                a = rand_elem_layer(P[0] < 0 ? (int)(rnd() % M) : P[0], 1); b = rand_elem_layer(P[1] < 0 ? (int)(rnd() % M) : P[1], 1);
                c2 = rand_elem_layer(P[2] < 0 ? (int)(rnd() % M) : P[2], 0); do c3 = rand_elem_layer(P[3] < 0 ? (int)(rnd() % M) : P[3], 0); while (c3 == c2);
            }
            else { do a = rnd() % n; while (!is_pair(a)); do b = rnd() % n; while (!is_pair(b)); c2 = 1 + rnd() % (n - 1); do c3 = 1 + rnd() % (n - 1); while (c3 == c2); }
            la = a; lb = b; lc2 = c2; lc3 = c3;
            d = lift_defect(-1, prof0, prof1);
            if (d < best) { best = d; bs[0] = to_gf(a); bs[1] = to_gf(b); bs[2] = to_gf(c2); bs[3] = to_gf(c3); }
        } else {
            int a, b, c;
            do a = rnd() % n; while (!is_pair(a));
            do b = rnd() % n; while (!is_pair(b) || b == a || b == inv(a));
            do c = rnd() % n; while (!is_pair(c) || c == a || c == inv(a) || c == b || c == inv(b));
            ga = a; gb = b; gc = c;
            d = cay_defect(-1, prof0);
            if (d < best) { best = d; bs[0] = to_gf(a); bs[1] = to_gf(b); bs[2] = to_gf(c); }
        }
        int bin = 0; if (d > 0) { bin = 1; while ((1 << bin) <= d) bin++; }
        hist[bin]++;
        suml += d > 0 ? bin - 0.5 : 0;
    }
    long f64 = 0, f256 = 0, f1k = 0;
    for (int b = 0; b <= 6; b++) f64 += hist[b];
    for (int b = 0; b <= 8; b++) f256 += hist[b];
    for (int b = 0; b <= 10; b++) f1k += hist[b];
    printf("RANDOM family=%s group=%s n=%d D=%d %ssamples=%ld f64=%.4f f256=%.4f f1k=%.4f ml=%.3f best_defect=%d best=[%d,%d,%d,%d] hist=",
           family == 1 ? "lift" : "cay", gfile, family == 1 ? 2 * n : n, D, label, K, (double)f64 / K, (double)f256 / K, (double)f1k / K, suml / K, best, bs[0], bs[1], bs[2], bs[3]);
    for (int b = 0; b < 40; b++) if (hist[b]) printf("[%s%d]=%ld ", b ? "<2^" : "", b ? b : 0, hist[b]);
    printf("\n"); fflush(stdout);
}

static void write_edges(int *s, const char *fn) {
    FILE *f = fopen(fn, "w"); if (!f) { perror(fn); exit(1); }
    long m = 0;
    if (family == 1) {
        int a = s[0], b = s[1], c2 = s[2], c3 = s[3], ia = inv(a), ib = inv(b);
        for (int g = 0; g < n; g++) {
            int G_ = to_gf(g), w;
            w = to_gf(mul(g, a)); if (G_ < w) { fprintf(f, "%d %d\n", G_, w); m++; }
            w = to_gf(mul(g, ia)); if (G_ < w) { fprintf(f, "%d %d\n", G_, w); m++; }
            fprintf(f, "%d %d\n", G_, n + G_); m++;
            fprintf(f, "%d %d\n", G_, n + to_gf(mul(g, c2))); m++;
            fprintf(f, "%d %d\n", G_, n + to_gf(mul(g, c3))); m++;
            w = to_gf(mul(g, b)); if (G_ < w) { fprintf(f, "%d %d\n", n + G_, n + w); m++; }
            w = to_gf(mul(g, ib)); if (G_ < w) { fprintf(f, "%d %d\n", n + G_, n + w); m++; }
        }
    } else {
        int gens[7] = {s[0], inv(s[0]), s[1], inv(s[1]), s[2], inv(s[2]), tinv}; int ng = 6 + use_tinv;
        for (int g = 0; g < n; g++) for (int i = 0; i < ng; i++) { int G_ = to_gf(g), w = to_gf(mul(g, gens[i])); if (G_ < w) { fprintf(f, "%d %d\n", G_, w); m++; } }
    }
    fclose(f);
    printf("edges written: %ld to %s\n", m, fn);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: see header\n"); return 2; }
    gfile = argv[1];
    const char *mode = "info"; long long set[4] = {0, 0, 0, 0}; const char *efile = NULL, *patfile = NULL; long samples = 1000; uint64_t seed = 7; int dry = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--mode")) mode = argv[++i];
        else if (!strcmp(argv[i], "--family")) { const char *f = argv[++i]; family = !strcmp(f, "lift") ? 1 : !strcmp(f, "cay") ? 2 : 0; }
        else if (!strcmp(argv[i], "--D")) D = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--set")) { for (int k = 0; k < 4; k++) { if (i + 1 < argc && strncmp(argv[i + 1], "--", 2)) set[k] = atoll(argv[++i]); } }
        else if (!strcmp(argv[i], "--pat")) { for (int k = 0; k < 4; k++) pat[k] = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--patfile")) patfile = argv[++i];
        else if (!strcmp(argv[i], "--edges")) efile = argv[++i];
        else if (!strcmp(argv[i], "--samples")) samples = atol(argv[++i]);
        else if (!strcmp(argv[i], "--seed")) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--chunk")) { chunk_i = atoi(argv[++i]); chunk_N = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--near")) near_T = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--time")) t_limit = atof(argv[++i]);
        else if (!strcmp(argv[i], "--tag")) tag = argv[++i];
        else if (!strcmp(argv[i], "--out")) out = fopen(argv[++i], "a");
        else if (!strcmp(argv[i], "--dry")) dry = 1;
        else if (!strcmp(argv[i], "--tinv")) use_tinv = 1;
        else if (!strcmp(argv[i], "--stabcap")) stabcap = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--unitcap")) unitcap = atoi(argv[++i]);
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    if (!family) die("--family lift|cay required");
    t_start = now();
    load(gfile);
    bfs_alloc();
    seed_rng(seed);
    fprintf(stderr, "group %s: M=%d nS=%d q=%d u=%lld |G|=%d |autS|=%d units=%d (gens %d) W=%d rows=%d loaded in %.2fs\n", gfile, M, nS, q, u, n, nautS, nallunits, nunits, W, R, now() - t_start);
    if (!strcmp(mode, "info")) {
        invariants();
        printf("INFO group=%s M=%d nS=%d q=%d u=%lld n=%d involutions=%d centre=%d classes=%d cp=%.6f |autS|=%d units=%d\nspec=%s\n", gfile, M, nS, q, u, n, ninvol, Zsize, nclasses, (double)nclasses / n, nautS, nallunits, spec);
        if (pat[0] >= 0 && family == 1) run_lift_pattern(1);
        if (family == 2 && dry) run_cay(1);
    } else if (!strcmp(mode, "eval")) {
        int p0[64] = {0}, p1[64] = {0};
        if (family == 1) {
            la = from_gf(set[0]); lb = from_gf(set[1]); lc2 = from_gf(set[2]); lc3 = from_gf(set[3]);
            int d = lift_defect(-1, p0, p1);
            printf("EVAL family=lift set=[%lld,%lld,%lld,%lld] D=%d root0 profile=", set[0], set[1], set[2], set[3], D);
            for (int l = 0; l <= D; l++) printf("%d ", p0[l]);
            printf("| root1 profile=");
            for (int l = 0; l <= D; l++) printf("%d ", p1[l]);
            printf("| total defect=%d\n", d);
        } else {
            ga = from_gf(set[0]); gb = from_gf(set[1]); gc = from_gf(set[2]);
            int d = cay_defect(-1, p0);
            printf("EVAL family=cay set=[%lld,%lld,%lld] D=%d defect=%d profile=", set[0], set[1], set[2], D, d);
            for (int l = 0; l <= D; l++) printf("%d ", p0[l]);
            printf("\n");
        }
    } else if (!strcmp(mode, "selftest")) {
        /* bit-parallel BFS against the plain BFS on random sets */
        int bad = 0;
        for (int i = 0; i < 20; i++) {
            int p0[64] = {0}, p1[64] = {0}, q0[64] = {0}, q1[64] = {0};
            if (family == 1) {
                do la = rnd() % n; while (!is_pair(la)); do lb = rnd() % n; while (!is_pair(lb)); lc2 = 1 + rnd() % (n - 1); do lc3 = 1 + rnd() % (n - 1); while (lc3 == lc2);
                int d = lift_defect(-1, p0, p1), e0 = plain_defect_lift(0, q0), e1 = plain_defect_lift(1, q1);
                if (d != e0 + e1 || memcmp(p0, q0, sizeof(p0)) || memcmp(p1, q1, sizeof(p1))) bad++;
            } else {
                do ga = rnd() % n; while (!is_pair(ga)); do gb = rnd() % n; while (!is_pair(gb)); do gc = rnd() % n; while (!is_pair(gc));
                int d = cay_defect(-1, p0), e = plain_defect_cay(q0);
                if (d != e || memcmp(p0, q0, sizeof(p0))) bad++;
            }
        }
        printf("SELFTEST group=%s family=%s D=%d mismatches=%d/20 %s\n", gfile, family == 1 ? "lift" : "cay", D, bad, bad ? "FAIL" : "OK");
    } else if (!strcmp(mode, "random")) {
        if (pat[0] >= 0) { char lab[64]; snprintf(lab, 64, "pat=[%d,%d,%d,%d] ", pat[0], pat[1], pat[2], pat[3]); random_screen(samples, pat, lab); }
        else random_screen(samples, NULL, "");
    } else if (!strcmp(mode, "patscreen")) {
        FILE *pf = fopen(patfile, "r"); if (!pf) die("patfile");
        int P[4];
        while (fscanf(pf, "%d %d %d %d", &P[0], &P[1], &P[2], &P[3]) == 4) { char lab[64]; snprintf(lab, 64, "pat=[%d,%d,%d,%d] ", P[0], P[1], P[2], P[3]); random_screen(samples, P, lab); }
        fclose(pf);
    } else if (!strcmp(mode, "exhaustive")) {
        if (family == 1) run_lift_pattern(dry); else run_cay(dry);
    } else if (!strcmp(mode, "cover")) {
        if (family == 1) run_lift_pattern(2); else run_cay(2);
    } else if (!strcmp(mode, "edges")) {
        int s[4] = {from_gf(set[0]), from_gf(set[1]), from_gf(set[2]), family == 1 ? from_gf(set[3]) : 0};
        write_edges(s, efile);
    } else die("mode");
    return 0;
}
