/* table_sweep.c: complete sweeps over table groups G = K x| C_M (GRPF group files of gf.py; element (k, m) = k y^m has index
 * k*M + m) for two families, with the symmetry group H = Inn(C) o B where C is the conjugating subgroup (C_G(t) for the
 * one-involution Cayley family, G for the lifts) and B is a list of automorphisms of K commuting with alpha (an optional
 * AUTK file written by aut_k.py: units of cyclic blocks, GL centralisers, ...; every beta extends to G by y -> y and is
 * re-verified at load).  Any subgroup of Aut(G) gives a complete reduction; a larger B only removes redundancy.
 *
 *   --family cay7 : Cay(G, {t} u {a, b, c}^+-), t = y^(M/2) (all involutions of G are K-conjugate when |K| is odd),
 *                   degree 7; candidates (a, b, c) up to H x inversions x permutations exactly as the Inn(G)-reduced sweep with H in
 *                   place of Inn(G): a = the smallest element of each H-orbit on inverse pairs, b = the smallest element of
 *                   each Stab_H(a)-orbit on the remaining pairs with type(b) >= type(a), c = every remaining pair
 *                   representative with type(c) >= type(a) and sorb(c) >= sorb(b)  (type = H-orbit index, sorb = Stab_H(a)
 *                   orbit index, both numbered by smallest element).
 *   --family lift : shape-103 two-vertex lift, vertices (v, g), loop a at 0, loop b at 1, links e, c2, c3 from 0 to 1
 *                   ((0, g) ~ (0, g a^+-1), (1, g), (1, g c2), (1, g c3); (1, g) ~ (1, g b^+-1)), n = 2|G|, degree 5.
 *                   Symmetries: H diagonal, a <-> a^-1, b <-> b^-1, the S_3 of re-gauging fibre 1 ((c2, c3) -> (c3, c2),
 *                   (c2^-1, c3 c2^-1), (c3 c2^-1, c2^-1), (c3^-1, c2 c3^-1), (c2 c3^-1, c3^-1), with b conjugated), the vertex
 *                   swap (a, b, c2, c3) -> (b, a, c2^-1, c3^-1).  Candidates: a = smallest element of each H-orbit on pairs,
 *                   c2 = smallest element of each Stab_H(a)-orbit on G \ {e}, c3 any other non-identity element such that the
 *                   key (sorb(c2), sorb(c3)) is lexicographically <= the keys of the five other S_3 images, b = every pair
 *                   representative (b < b^-1) with type(b) >= type(a).  Two BFS (from (0, e) and (1, e)) decide the diameter
 *                   since G acts transitively on each fibre by left multiplication; defect = lift vertices beyond D from
 *                   the two roots.
 *
 *   table_sweep GRP [--aut AUTFILE] --family F --D D --mode info|random|exhaustive|eval|edges|cover [--samples K] [--seed s]
 *       [--chunk i N] [--near T] [--time S] [--tag TAG] [--out results.jsonl] [--dry] [--set a b c | a b c2 c3] [--edges FILE]
 *   Lines: FOUND / NEAR / PROGRESS / EXHAUSTIVE ... complete=1, RANDOM (screen), INFO.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static int nK, M, n, D = 6;
static uint16_t *TK, *alpha, *apow, *RT;   /* RT[s*n + g] = g*s */
static int *inv;
static char *spec;
static const char *gfile, *autfile = NULL;
static int family = 0;                     /* 1 = cay7, 2 = lift */
static uint16_t *B; static int nB = 0;     /* automorphisms of K (nB x nK), B[0] = identity */
static int *C; static int nC;              /* conjugating subgroup */
static int tinv = -1;                      /* the involution t = y^(M/2) for cay7 */

static double now(void) { struct timeval tv; gettimeofday(&tv, NULL); return tv.tv_sec + 1e-6 * tv.tv_usec; }
static void die(const char *m) { fprintf(stderr, "FATAL: %s\n", m); exit(1); }
static inline int mul(int a, int b) { return RT[(size_t)b * n + a]; }
static inline int cnj(int c, int g) { return RT[(size_t)inv[c] * n + RT[(size_t)g * n + c]]; }   /* c g c^-1 */
static inline int bapply(int bi, int g) { return B[(size_t)bi * nK + g / M] * M + g % M; }

static uint64_t rs[2] = {0x9E3779B97F4A7C15ULL, 0xD1B54A32D192ED03ULL};
static inline uint64_t rotl(uint64_t x, int r) { return (x << r) | (x >> (64 - r)); }
static uint64_t rnd(void) { uint64_t s0 = rs[0], s1 = rs[1], r = s0 + s1; s1 ^= s0; rs[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); rs[1] = rotl(s1, 36); return r; }

static void load(const char *fn) {
    FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
    char magic[4]; int32_t hdr[2];
    if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "GRPF", 4)) die("bad magic");
    if (fread(hdr, 4, 2, f) != 2) die("bad header");
    nK = hdr[0]; M = hdr[1]; n = nK * M;
    if (n >= 65536 || nK < 1 || M < 1) die("n must be < 65536");
    TK = malloc((size_t)nK * nK * 2); alpha = malloc((size_t)nK * 2);
    if (fread(TK, 2, (size_t)nK * nK, f) != (size_t)nK * nK) die("short TK");
    if (fread(alpha, 2, nK, f) != (size_t)nK) die("short alpha");
    spec = calloc(1 << 16, 1); size_t got = fread(spec, 1, (1 << 16) - 1, f); spec[got] = 0;
    fclose(f);
    for (int a = 0; a < nK; a++) { if (TK[a] != a || TK[(size_t)a * nK] != a) die("K: 0 is not the identity"); }
    char *seen = calloc(nK, 1);
    for (int a = 0; a < nK; a++) {
        memset(seen, 0, nK);
        for (int b = 0; b < nK; b++) { int c = TK[(size_t)a * nK + b]; if (c >= nK || seen[c]) die("K: row not a permutation"); seen[c] = 1; }
    }
    memset(seen, 0, nK);
    for (int a = 0; a < nK; a++) { if (alpha[a] >= nK || seen[alpha[a]]) die("alpha not a bijection"); seen[alpha[a]] = 1; }
    for (int a = 0; a < nK; a++) for (int b = 0; b < nK; b++)
        if (alpha[TK[(size_t)a * nK + b]] != TK[(size_t)alpha[a] * nK + alpha[b]]) die("alpha not a homomorphism");
    apow = malloc((size_t)M * nK * 2);
    for (int k = 0; k < nK; k++) apow[k] = k;
    for (int m = 1; m < M; m++) for (int k = 0; k < nK; k++) apow[(size_t)m * nK + k] = alpha[apow[(size_t)(m - 1) * nK + k]];
    for (int k = 0; k < nK; k++) if (alpha[apow[(size_t)(M - 1) * nK + k]] != k) die("alpha^M != id");
    free(seen);
    RT = malloc((size_t)n * n * 2); if (!RT) die("RT malloc");
    for (int s = 0; s < n; s++) {
        int k2 = s / M, m2 = s % M;
        uint16_t *row = RT + (size_t)s * n;
        for (int g = 0; g < n; g++) {
            int k1 = g / M, m1 = g % M;
            row[g] = (uint16_t)(TK[(size_t)k1 * nK + apow[(size_t)m1 * nK + k2]] * M + (m1 + m2) % M);
        }
    }
    inv = malloc(n * sizeof(int));
    for (int g = 0; g < n; g++) inv[g] = -1;
    for (int s = 0; s < n; s++) { if (RT[(size_t)s * n + 0] != s) die("G: identity"); }
    for (int g = 0; g < n; g++) for (int s = 0; s < n; s++) if (RT[(size_t)s * n + g] == 0) { inv[g] = s; break; }
    for (int g = 0; g < n; g++) { if (inv[g] < 0 || RT[(size_t)g * n + inv[g]] != 0) die("G: inverses"); }
    uint64_t x = 88172645463325252ULL;
    for (int i = 0; i < 2000000; i++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        int g = x % n, s = (x >> 20) % n, t = (x >> 40) % n;
        if (RT[(size_t)t * n + RT[(size_t)s * n + g]] != RT[(size_t)RT[(size_t)t * n + s] * n + g]) die("G: associativity");
    }
}

static void load_aut(const char *fn) {
    /* AUTK: magic, int32 nB, int32 nK, uint16 B[nB*nK]; identity prepended if missing; every beta verified */
    int nb = 0; uint16_t *Bf = NULL;
    if (fn) {
        FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
        char magic[4]; int32_t hdr[2];
        if (fread(magic, 1, 4, f) != 4 || memcmp(magic, "AUTK", 4)) die("aut: bad magic");
        if (fread(hdr, 4, 2, f) != 2) die("aut: bad header");
        nb = hdr[0]; if (hdr[1] != nK) die("aut: nK mismatch");
        Bf = malloc((size_t)nb * nK * 2);
        if (fread(Bf, 2, (size_t)nb * nK, f) != (size_t)nb * nK) die("aut: short");
        fclose(f);
    }
    B = malloc((size_t)(nb + 1) * nK * 2);
    for (int k = 0; k < nK; k++) B[k] = k;
    nB = 1;
    char *seen = calloc(nK, 1);
    uint64_t x = 0x1234567887654321ULL;
    for (int i = 0; i < nb; i++) {
        uint16_t *b = Bf + (size_t)i * nK;
        int isid = 1; for (int k = 0; k < nK; k++) if (b[k] != k) { isid = 0; break; }
        if (isid) continue;
        memset(seen, 0, nK);
        for (int k = 0; k < nK; k++) { if (b[k] >= nK || seen[b[k]]) die("aut: beta not a bijection"); seen[b[k]] = 1; }
        for (int k = 0; k < nK; k++) if (b[alpha[k]] != alpha[b[k]]) die("aut: beta does not commute with alpha");
        if (nK <= 600 || i < 3) {
            for (int a = 0; a < nK; a++) for (int c = 0; c < nK; c++) if (b[TK[(size_t)a * nK + c]] != TK[(size_t)b[a] * nK + b[c]]) die("aut: beta not a homomorphism");
        } else {
            for (int s = 0; s < 50000; s++) {
                x ^= x << 13; x ^= x >> 7; x ^= x << 17;
                int a = x % nK, c = (x >> 24) % nK;
                if (b[TK[(size_t)a * nK + c]] != TK[(size_t)b[a] * nK + b[c]]) die("aut: beta not a homomorphism (sample)");
            }
        }
        memcpy(B + (size_t)nB * nK, b, (size_t)nK * 2); nB++;
    }
    free(seen); free(Bf);
    /* G-level check on a sample: phi(g h) = phi(g) phi(h) with phi = beta extended by y -> y */
    for (int i = 1; i < nB; i++) for (int s = 0; s < 20000; s++) {
        x ^= x << 13; x ^= x >> 7; x ^= x << 17;
        int g = x % n, h = (x >> 24) % n;
        if (bapply(i, mul(g, h)) != mul(bapply(i, g), bapply(i, h))) die("aut: beta extension is not an automorphism of G");
    }
}

/* ---------------- BFS ---------------- */
static uint8_t *seen; static int *queue; static long bfs_calls = 0;
static int gens[7];           /* cay7: t a a^-1 b b^-1 c c^-1 */
static int la, lb, lc2, lc3;  /* lift voltages */

/* cay7: BFS from the identity over the 7 arcs; returns the defect (vertices beyond D), or a value > bound_T that is a
 * lower bound when the defect provably exceeds bound_T (bound_T >= 0). */
static int bfs7(int bound_T, int *prof) {
    bfs_calls++;
    memset(seen, 0, n);
    seen[0] = 1; queue[0] = 0;
    int qh = 0, qt = 1, reached = 1;
    const uint16_t *R[7];
    for (int i = 0; i < 7; i++) R[i] = RT + (size_t)gens[i] * n;
    if (prof) prof[0] = 1;
    for (int lev = 1; lev <= D; lev++) {
        int end = qt;
        if (lev == D && bound_T >= 0) {
            int rem = end - qh;
            if (n - reached - 6 * rem > bound_T) return n - reached - 6 * rem;
        }
        while (qh < end) {
            int v = queue[qh++];
            for (int i = 0; i < 7; i++) { int w = R[i][v]; if (!seen[w]) { seen[w] = 1; queue[qt++] = w; } }
            if (lev == D && bound_T >= 0 && (qh & 255) == 0) {
                int rem = end - qh; reached = qt;
                if (n - reached - 6 * rem > bound_T) return n - reached - 6 * rem;
            }
        }
        reached = qt;
        if (prof) prof[lev] = qt - end;
        if (qh == qt) break;
    }
    return n - reached;
}

/* lift: BFS from (root, e) over the lift (2n vertices, index v*n + g); returns the number of lift vertices beyond D
 * from that root, or a lower bound > bound_T.  prof[] optional. */
static int bfs_lift(int root, int bound_T, int *prof) {
    bfs_calls++;
    int N2 = 2 * n;
    memset(seen, 0, N2);
    int start = root * n;
    seen[start] = 1; queue[0] = start;
    int qh = 0, qt = 1, reached = 1;
    const uint16_t *Ra = RT + (size_t)la * n, *Rai = RT + (size_t)inv[la] * n;
    const uint16_t *Rb = RT + (size_t)lb * n, *Rbi = RT + (size_t)inv[lb] * n;
    const uint16_t *R2 = RT + (size_t)lc2 * n, *R3 = RT + (size_t)lc3 * n;
    const uint16_t *R2i = RT + (size_t)inv[lc2] * n, *R3i = RT + (size_t)inv[lc3] * n;
    if (prof) prof[0] = 1;
    for (int lev = 1; lev <= D; lev++) {
        int end = qt;
        if (lev == D && bound_T >= 0) {
            int rem = end - qh;
            if (N2 - reached - 4 * rem > bound_T) return N2 - reached - 4 * rem;
        }
        while (qh < end) {
            int v = queue[qh++];
            int g = v % n; int w;
            if (v < n) {
                w = Ra[g];      if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
                w = Rai[g];     if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
                w = n + g;      if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
                w = n + R2[g];  if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
                w = n + R3[g];  if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
            } else {
                w = n + Rb[g];  if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
                w = n + Rbi[g]; if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
                w = g;          if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
                w = R2i[g];     if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
                w = R3i[g];     if (!seen[w]) { seen[w] = 1; queue[qt++] = w; }
            }
            if (lev == D && bound_T >= 0 && (qh & 255) == 0) {
                int rem = end - qh; reached = qt;
                if (N2 - reached - 4 * rem > bound_T) return N2 - reached - 4 * rem;
            }
        }
        reached = qt;
        if (prof) prof[lev] = qt - end;
        if (qh == qt) break;
    }
    return N2 - reached;
}

static int lift_defect(int bound_T) {
    int d0 = bfs_lift(0, bound_T, NULL);
    if (bound_T >= 0 && d0 > bound_T) return d0;
    int d1 = bfs_lift(1, bound_T >= 0 ? bound_T - d0 : -1, NULL);
    return d0 + d1;
}

/* ---------------- structure: classes, H-orbits, stabilisers ---------------- */
static int *cls, ncls, *clsrep;
static void classes(void) {
    cls = malloc(n * sizeof(int)); for (int g = 0; g < n; g++) cls[g] = -1;
    ncls = 0; int cap = 1024; clsrep = malloc(cap * sizeof(int));
    for (int g = 0; g < n; g++) if (cls[g] < 0) {
        if (ncls == cap) { cap *= 2; clsrep = realloc(clsrep, cap * sizeof(int)); }
        clsrep[ncls] = g;
        for (int x = 0; x < n; x++) cls[cnj(x, g)] = ncls;
        ncls++;
    }
}

static int *uf;
static int ufind(int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; }
static void uunion(int a, int b) { a = ufind(a); b = ufind(b); if (a != b) { if (a < b) uf[b] = a; else uf[a] = b; } }
/* number orbits by their smallest element: out[g] = orbit index, omin[idx] = smallest element; returns #orbits */
static int number_orbits(int *out, int *omin) {
    int cnt = 0;
    for (int g = 0; g < n; g++) if (ufind(g) == g) { omin[cnt] = g; cnt++; }
    /* ufind root is the minimum element of the component because unions keep the smaller index as root */
    int *idx = malloc(n * sizeof(int));
    for (int i = 0; i < cnt; i++) idx[omin[i]] = i;
    for (int g = 0; g < n; g++) out[g] = idx[ufind(g)];
    free(idx);
    return cnt;
}

static int *type_, ntype, *typemin;       /* H-orbits on inverse pairs */
static void h_orbits(void) {
    uf = malloc(n * sizeof(int)); for (int g = 0; g < n; g++) uf[g] = g;
    for (int g = 0; g < n; g++) uunion(g, inv[g]);
    for (int i = 0; i < nC; i++) { int c = C[i]; for (int g = 0; g < n; g++) uunion(g, cnj(c, g)); }
    for (int bi = 1; bi < nB; bi++) for (int g = 0; g < n; g++) uunion(g, bapply(bi, g));
    type_ = malloc(n * sizeof(int)); typemin = malloc(n * sizeof(int));
    ntype = number_orbits(type_, typemin);
}

/* Stab_H({a, a^-1}) orbits on G (merge inverses if pairs = 1); returns the number of orbits; sorb[] / somin[] filled;
 * also returns the stabiliser's generator count through *ngen */
static int *cconj;                          /* cconj[g] = one c in C with c a c^-1 = g, or -1 */
/* explicit Stab_H({a, a^-1}) = { (c' c0(beta), beta) : beta in B', c' in stabC }, stabC = {c in C : c a c^-1 in {a, a^-1}} */
static int *stabC, nstabC, *stabB, *stabB0, nstabB;
static int *sorb2, *somin2;
static int stab_orbits(int a, int pairs, int *sorb, int *somin, long *ngen) {
    for (int g = 0; g < n; g++) uf[g] = g;
    if (pairs) for (int g = 0; g < n; g++) uunion(g, inv[g]);
    for (int g = 0; g < n; g++) cconj[g] = -1;
    for (int i = 0; i < nC; i++) { int c = C[i], g = cnj(c, a); if (cconj[g] < 0) cconj[g] = c; }
    long gens_ = 0;
    int ia = inv[a];
    nstabC = 0; nstabB = 0;
    for (int i = 0; i < nC; i++) {
        int c = C[i], g = cnj(c, a);
        if (g == a || g == ia) { gens_++; stabC[nstabC++] = c; for (int h = 0; h < n; h++) uunion(h, cnj(c, h)); }
    }
    stabB[nstabB] = 0; stabB0[nstabB] = 0; nstabB++;
    for (int bi = 1; bi < nB; bi++) {
        int x = bapply(bi, a);
        int c1 = cconj[x];
        if (c1 < 0) c1 = cconj[inv[x]];
        if (c1 >= 0) {
            int c0 = inv[c1]; gens_++; stabB[nstabB] = bi; stabB0[nstabB] = c0; nstabB++;
            for (int h = 0; h < n; h++) uunion(h, cnj(c0, bapply(bi, h)));
        }
    }
    *ngen = gens_;
    return number_orbits(sorb, somin);
}

/* orbits on G of the pointwise stabiliser of c2 inside the explicit Stab_H(a) (lift family) */
static int point_stab_orbits(int c2, int *sorb2, int *somin2, long *nelem) {
    for (int g = 0; g < n; g++) uf[g] = g;
    long cnt = 0;
    for (int j = 0; j < nstabB; j++) {
        int bi = stabB[j], c0 = stabB0[j];
        int x = bapply(bi, c2);
        for (int i = 0; i < nstabC; i++) {
            int c = mul(stabC[i], c0);                 /* c' c0 */
            if (cnj(c, x) != c2) continue;
            cnt++;
            if (bi == 0 && c == 0) continue;               /* the identity */
            /* the orbits of the subgroup generated by the first 96 stabiliser elements found: a subgroup of the pointwise
             * stabiliser keeps the reduction complete; the cap bounds the union-find cost for central c2 */
            if (cnt > 96) continue;
            for (int h = 0; h < n; h++) uunion(h, cnj(c, bapply(bi, h)));
        }
    }
    *nelem = cnt;
    return number_orbits(sorb2, somin2);
}

static int allowed(int g) { return g != 0 && inv[g] != g; }

static inline void lift_key(int x, int y, const int *sorb, int *k) {
    k[0] = sorb[x]; k[1] = sorb[y]; k[2] = sorb[mul(inv[x], y)]; k[3] = sorb[mul(y, inv[x])];
}
static inline int lift_key_less(int x, int y, const int *sorb, const int *k0) {
    int k[4]; lift_key(x, y, sorb, k);
    for (int i = 0; i < 4; i++) { if (k[i] < k0[i]) return 1; if (k[i] > k0[i]) return 0; }
    return 0;
}

/* ---------------- enumeration ---------------- */
static const char *tag = ""; static FILE *out = NULL;
static int near_T = 0; static double t_limit = 1e18, t_start;
static long ncand = 0, nmine = 0, nfound = 0; static int best_defect = 1 << 30, best_set[4] = {0, 0, 0, 0};
static int chunk_i = 0, chunk_N = 1, stopped = 0;

static void report(int *s, int d) {
    if (d == 0) {
        nfound++;
        if (family == 1) printf("FOUND diam<=%d (defect 0) family=cay7 group=%s set=[%d,%d,%d]\n", D, gfile, s[0], s[1], s[2]);
        else printf("FOUND diam<=%d (defect 0) family=lift group=%s set=[%d,%d,%d,%d]\n", D, gfile, s[0], s[1], s[2], s[3]);
        fflush(stdout);
        if (out) { fprintf(out, "{\"event\":\"found\",\"tag\":\"%s\",\"family\":%d,\"D\":%d,\"n\":%d,\"group\":\"%s\",\"set\":[%d,%d,%d,%d]}\n", tag, family, D, n, gfile, s[0], s[1], s[2], s[3]); fflush(out); }
    } else if (d <= near_T) {
        if (family == 1) printf("NEAR defect=%d set=[%d,%d,%d]\n", d, s[0], s[1], s[2]);
        else printf("NEAR defect=%d set=[%d,%d,%d,%d]\n", d, s[0], s[1], s[2], s[3]);
        fflush(stdout);
    }
    if (d < best_defect) { best_defect = d; memcpy(best_set, s, 4 * sizeof(int)); }
}

static void consider7(int a, int b, int c) {
    long idx = ncand++;
    if (idx % chunk_N != chunk_i) return;
    nmine++;
    gens[0] = tinv; gens[1] = a; gens[2] = inv[a]; gens[3] = b; gens[4] = inv[b]; gens[5] = c; gens[6] = inv[c];
    int T = near_T > best_defect ? near_T : best_defect;
    int d = bfs7(T, NULL);
    if (d <= T) { int s[4] = {a, b, c, 0}; report(s, d); }
}

static void consider_lift(int a, int b, int c2, int c3) {
    long idx = ncand++;
    if (idx % chunk_N != chunk_i) return;
    nmine++;
    la = a; lb = b; lc2 = c2; lc3 = c3;
    int T = near_T > best_defect ? near_T : best_defect;
    int d = lift_defect(T);
    if (d <= T) { int s[4] = {a, b, c2, c3}; report(s, d); }
}

/* cover test (small groups): bitsets over triples (cay7) / quadruples (lift) of the enumerated candidates' full orbits */
static uint8_t *covered; static int *pairidx, npair;
static inline void setbit(size_t bit) { covered[bit >> 3] |= 1 << (bit & 7); }
static inline int getbit(size_t bit) { return (covered[bit >> 3] >> (bit & 7)) & 1; }
static void mark7(int a, int b, int c) {
    int p[3] = {pairidx[a], pairidx[b], pairidx[c]}; int t;
    if (p[0] > p[1]) { t = p[0]; p[0] = p[1]; p[1] = t; }
    if (p[1] > p[2]) { t = p[1]; p[1] = p[2]; p[2] = t; }
    if (p[0] > p[1]) { t = p[0]; p[0] = p[1]; p[1] = t; }
    setbit(((size_t)p[0] * npair + p[1]) * npair + p[2]);
}
static void mark_orbit7(int a, int b, int c) {
    for (int bi = 0; bi < nB; bi++) for (int i = 0; i < nC; i++) {
        int x = C[i];
        mark7(cnj(x, bapply(bi, a)), cnj(x, bapply(bi, b)), cnj(x, bapply(bi, c)));
    }
}
static void mark_lift_tuple(int a, int b, int c2, int c3) {
    /* all 48 images under inversions, S_3 re-gauging, swap */
    int A[2] = {a, inv[a]}, Bb[2] = {b, inv[b]};
    for (int ia = 0; ia < 2; ia++) for (int ib = 0; ib < 2; ib++) for (int sw = 0; sw < 2; sw++) {
        int a_ = A[ia], b_ = Bb[ib], x2 = c2, x3 = c3;
        if (sw) { int t = a_; a_ = b_; b_ = t; x2 = inv[c2]; x3 = inv[c3]; }
        /* S_3 regauge: h in {e, x2^-1, x3^-1}: links (c h), loop b -> h^-1 b h */
        int hs[3] = {0, inv[x2], inv[x3]};
        for (int hi = 0; hi < 3; hi++) {
            int h = hs[hi];
            int l1 = mul(x2, h), l2 = mul(x3, h), l0 = h;         /* new links: h, x2 h, x3 h ; one of them is e */
            int bb = mul(mul(inv[h], b_), h);
            int links[3] = {l0, l1, l2}; int u[2], k = 0;
            for (int j = 0; j < 3; j++) if (links[j] != 0) { if (k < 2) u[k] = links[j]; k++; }
            if (k != 2) continue;
            setbit((((size_t)a_ * n + bb) * n + u[0]) * n + u[1]);
            setbit((((size_t)a_ * n + bb) * n + u[1]) * n + u[0]);
        }
    }
}
static void mark_orbit_lift(int a, int b, int c2, int c3) {
    for (int bi = 0; bi < nB; bi++) for (int i = 0; i < nC; i++) {
        int x = C[i];
        mark_lift_tuple(cnj(x, bapply(bi, a)), cnj(x, bapply(bi, b)), cnj(x, bapply(bi, c2)), cnj(x, bapply(bi, c3)));
    }
}

static void run_exhaustive(int dry) {
    /* dry: 1 = count only, 2 = cover test */
    h_orbits();
    int *sorb = malloc(n * sizeof(int)), *somin = malloc(n * sizeof(int));
    sorb2 = malloc(n * sizeof(int)); somin2 = malloc(n * sizeof(int));
    cconj = malloc(n * sizeof(int));
    stabC = malloc(nC * sizeof(int)); stabB = malloc((nB + 1) * sizeof(int)); stabB0 = malloc((nB + 1) * sizeof(int));
    long pstab_total = 0;
    if (dry == 2) {
        if (family == 1) {
            pairidx = malloc(n * sizeof(int)); npair = 0;
            for (int g = 0; g < n; g++) pairidx[g] = -1;
            for (int g = 0; g < n; g++) if (allowed(g) && g < inv[g]) { pairidx[g] = pairidx[inv[g]] = npair++; }
            size_t bits = (size_t)npair * npair * npair;
            covered = calloc(bits / 8 + 1, 1);
        } else {
            if (n > 80) die("cover mode (lift) needs n <= 80");
            size_t bits = (size_t)n * n * n * n;
            covered = calloc(bits / 8 + 1, 1);
        }
    }
    double t_last = now();
    long stabgens_total = 0;
    for (int ta = 0; ta < ntype && !stopped; ta++) {
        int a = typemin[ta];
        if (!allowed(a)) continue;
        long ng;
        int nso = stab_orbits(a, family == 1, sorb, somin, &ng);
        stabgens_total += ng;
        if (family == 1) {
            for (int b = 0; b < n && !stopped; b++) {
                if (!allowed(b) || b == a || b == inv[a] || type_[b] < ta) continue;
                if (somin[sorb[b]] != b) continue;
                int ob = sorb[b];
                for (int c = b + 1; c < n; c++) {
                    if (!allowed(c) || c > inv[c] || c == a || c == inv[a] || c == inv[b]) continue;
                    if (type_[c] < ta || sorb[c] < ob) continue;
                    if (dry) { ncand++; if (dry == 2) mark_orbit7(a, b, c); continue; }
                    consider7(a, b, c);
                }
                if (now() - t_last > 60) {
                    t_last = now();
                    printf("PROGRESS type %d/%d b=%d candidates=%ld mine=%ld bfs=%ld best=%d %.0fs\n", ta, ntype, b, ncand, nmine, bfs_calls, best_defect, now() - t_start); fflush(stdout);
                }
                if (now() - t_start > t_limit) stopped = 1;
            }
        } else {
            (void)nso;
            for (int c2 = 1; c2 < n && !stopped; c2++) {
                if (somin[sorb[c2]] != c2) continue;
                long nel;
                point_stab_orbits(c2, sorb2, somin2, &nel);
                pstab_total += nel;
                for (int c3 = 1; c3 < n; c3++) {
                    if (c3 == c2 || somin2[sorb2[c3]] != c3) continue;
                    /* Stab_H(a)-invariant key of an ordered pair (x, y): (sorb x, sorb y, sorb x^-1 y, sorb y x^-1); the
                     * candidate's key must be lexicographically <= the keys of its five other S_3 images */
                    int i2 = inv[c2], i3 = inv[c3];
                    int k0[4]; lift_key(c2, c3, sorb, k0);
                    int u = mul(c3, i2), w = mul(c2, i3);
                    if (lift_key_less(c3, c2, sorb, k0)) continue;
                    if (lift_key_less(i2, u, sorb, k0)) continue;
                    if (lift_key_less(u, i2, sorb, k0)) continue;
                    if (lift_key_less(i3, w, sorb, k0)) continue;
                    if (lift_key_less(w, i3, sorb, k0)) continue;
                    for (int b = 1; b < n; b++) {
                        if (!allowed(b) || b > inv[b] || type_[b] < ta) continue;
                        if (dry) { ncand++; if (dry == 2) mark_orbit_lift(a, b, c2, c3); continue; }
                        consider_lift(a, b, c2, c3);
                    }
                }
                if (now() - t_last > 60) {
                    t_last = now();
                    printf("PROGRESS type %d/%d c2=%d candidates=%ld mine=%ld bfs=%ld best=%d %.0fs\n", ta, ntype, c2, ncand, nmine, bfs_calls, best_defect, now() - t_start); fflush(stdout);
                }
                if (now() - t_start > t_limit) stopped = 1;
            }
        }
    }
    if (dry) {
        printf("COUNT family=%s n=%d |C|=%d |B|=%d types=%d stabgens=%ld pstab=%ld candidates(total over all chunks)=%ld time=%.1fs\n",
               family == 1 ? "cay7" : "lift", n, nC, nB, ntype, stabgens_total, pstab_total, ncand, now() - t_start);
        if (dry == 2) {
            long total = 0, miss = 0;
            if (family == 1) {
                for (int p0 = 0; p0 < npair; p0++) for (int p1 = p0 + 1; p1 < npair; p1++) for (int p2 = p1 + 1; p2 < npair; p2++) {
                    total++;
                    if (!getbit(((size_t)p0 * npair + p1) * npair + p2)) miss++;
                }
            } else {
                for (int a = 1; a < n; a++) { if (!allowed(a)) continue;
                    for (int b = 1; b < n; b++) { if (!allowed(b)) continue;
                        for (int c2 = 1; c2 < n; c2++) for (int c3 = 1; c3 < n; c3++) { if (c3 == c2) continue;
                            total++;
                            if (!getbit((((size_t)a * n + b) * n + c2) * n + c3)) miss++;
                        } } }
            }
            printf("COVER total=%ld uncovered=%ld %s\n", total, miss, miss ? "FAIL" : "OK");
        }
        return;
    }
    double el = now() - t_start;
    printf("EXHAUSTIVE family=%s group=%s D=%d n=%d chunk=%d/%d candidates=%ld mine=%ld bfs=%ld time=%.1fs rate=%.0f cand/s best_defect=%d best=[%d,%d,%d,%d] found=%ld complete=%d\n",
           family == 1 ? "cay7" : "lift", gfile, D, family == 1 ? n : 2 * n, chunk_i, chunk_N, ncand, nmine, bfs_calls, el, nmine / (el > 0 ? el : 1), best_defect, best_set[0], best_set[1], best_set[2], best_set[3], nfound, !stopped);
    if (out) {
        fprintf(out, "{\"event\":\"exhaustive\",\"tag\":\"%s\",\"family\":%d,\"D\":%d,\"n\":%d,\"group\":\"%s\",\"chunk\":[%d,%d],\"candidates\":%ld,\"mine\":%ld,\"time\":%.1f,\"best_defect\":%d,\"best\":[%d,%d,%d,%d],\"found\":%ld,\"complete\":%d}\n",
                tag, family, D, family == 1 ? n : 2 * n, gfile, chunk_i, chunk_N, ncand, nmine, el, best_defect, best_set[0], best_set[1], best_set[2], best_set[3], nfound, !stopped);
        fflush(out);
    }
}

static void run_random(long K, uint64_t seed) {
    rs[0] = seed * 0x9E3779B97F4A7C15ULL + 1; rs[1] = (seed ^ 0xDEADBEEFCAFEULL) * 0xBF58476D1CE4E5B9ULL + 7;
    for (int i = 0; i < 20; i++) rnd();
    int *reps = malloc(n * sizeof(int)); int nr = 0;
    for (int g = 0; g < n; g++) if (allowed(g) && g < inv[g]) reps[nr++] = g;
    long hist[40] = {0}; int best = 1 << 30, bs[4] = {0, 0, 0, 0}; double suml = 0;
    for (long i = 0; i < K; i++) {
        int d;
        if (family == 1) {
            int a, b, c;
            a = reps[rnd() % nr]; do b = reps[rnd() % nr]; while (b == a); do c = reps[rnd() % nr]; while (c == a || c == b);
            gens[0] = tinv; gens[1] = a; gens[2] = inv[a]; gens[3] = b; gens[4] = inv[b]; gens[5] = c; gens[6] = inv[c];
            d = bfs7(-1, NULL);
            if (d < best) { best = d; bs[0] = a; bs[1] = b; bs[2] = c; }
        } else {
            int a = reps[rnd() % nr], b = reps[rnd() % nr]; if (rnd() & 1) a = inv[a]; if (rnd() & 1) b = inv[b];
            int c2, c3; c2 = 1 + rnd() % (n - 1); do c3 = 1 + rnd() % (n - 1); while (c3 == c2);
            la = a; lb = b; lc2 = c2; lc3 = c3;
            d = lift_defect(-1);
            if (d < best) { best = d; bs[0] = a; bs[1] = b; bs[2] = c2; bs[3] = c3; }
        }
        int bin = 0; if (d > 0) { bin = 1; while ((1 << bin) <= d) bin++; }
        hist[bin]++;
        suml += d > 0 ? bin - 0.5 : 0;
    }
    long f64 = 0, f256 = 0;
    for (int b = 0; b <= 6; b++) f64 += hist[b];
    for (int b = 0; b <= 8; b++) f256 += hist[b];
    printf("RANDOM family=%s group=%s n=%d D=%d samples=%ld f64=%.4f f256=%.4f ml=%.2f best_defect=%d best=[%d,%d,%d,%d] classes=%d hist=",
           family == 1 ? "cay7" : "lift", gfile, family == 1 ? n : 2 * n, D, K, (double)f64 / K, (double)f256 / K, suml / K, best, bs[0], bs[1], bs[2], bs[3], ncls);
    for (int b = 0; b < 40; b++) if (hist[b]) printf("[%s%d]=%ld ", b ? "<2^" : "", b ? b : 0, hist[b]);
    printf("\n");
}

static void write_edges(int *s, const char *fn) {
    FILE *f = fopen(fn, "w"); if (!f) { perror(fn); exit(1); }
    long m = 0;
    if (family == 1) {
        gens[0] = tinv; gens[1] = s[0]; gens[2] = inv[s[0]]; gens[3] = s[1]; gens[4] = inv[s[1]]; gens[5] = s[2]; gens[6] = inv[s[2]];
        for (int g = 0; g < n; g++) for (int i = 0; i < 7; i++) { int w = RT[(size_t)gens[i] * n + g]; if (g < w) { fprintf(f, "%d %d\n", g, w); m++; } }
    } else {
        int a = s[0], b = s[1], c2 = s[2], c3 = s[3];
        for (int g = 0; g < n; g++) {
            int w;
            w = mul(g, a); if (g < w) { fprintf(f, "%d %d\n", g, w); m++; }
            w = mul(g, inv[a]); if (g < w) { fprintf(f, "%d %d\n", g, w); m++; }
            fprintf(f, "%d %d\n", g, n + g); m++;
            fprintf(f, "%d %d\n", g, n + mul(g, c2)); m++;
            fprintf(f, "%d %d\n", g, n + mul(g, c3)); m++;
            w = mul(g, b); if (g < w) { fprintf(f, "%d %d\n", n + g, n + w); m++; }
            w = mul(g, inv[b]); if (g < w) { fprintf(f, "%d %d\n", n + g, n + w); m++; }
        }
    }
    fclose(f);
    printf("edges written: %ld to %s\n", m, fn);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: see header\n"); return 2; }
    gfile = argv[1];
    const char *mode = "info"; int set[4] = {0, 0, 0, 0}; const char *efile = NULL; long samples = 1000; uint64_t seed = 7; int dry = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--mode")) mode = argv[++i];
        else if (!strcmp(argv[i], "--family")) { const char *f = argv[++i]; family = !strcmp(f, "cay7") ? 1 : !strcmp(f, "lift") ? 2 : 0; }
        else if (!strcmp(argv[i], "--aut")) autfile = argv[++i];
        else if (!strcmp(argv[i], "--D")) D = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--set")) { set[0] = atoi(argv[++i]); set[1] = atoi(argv[++i]); set[2] = atoi(argv[++i]); if (i + 1 < argc && argv[i + 1][0] != '-') set[3] = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--edges")) efile = argv[++i];
        else if (!strcmp(argv[i], "--samples")) samples = atol(argv[++i]);
        else if (!strcmp(argv[i], "--seed")) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--chunk")) { chunk_i = atoi(argv[++i]); chunk_N = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--near")) near_T = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--time")) t_limit = atof(argv[++i]);
        else if (!strcmp(argv[i], "--tag")) tag = argv[++i];
        else if (!strcmp(argv[i], "--out")) out = fopen(argv[++i], "a");
        else if (!strcmp(argv[i], "--dry")) dry = 1;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 2; }
    }
    if (!family) die("--family cay7|lift required");
    t_start = now();
    load(gfile);
    load_aut(autfile);
    seen = malloc(2 * (size_t)n); queue = malloc(2 * (size_t)n * sizeof(int));
    int ninv = 0; for (int g = 1; g < n; g++) if (inv[g] == g) ninv++;
    if (family == 1) {
        if (M % 2) die("cay7 needs M even");
        tinv = M / 2;
        if (inv[tinv] != tinv || tinv == 0) die("t = y^(M/2) is not an involution");
        /* all involutions must be conjugate to t (true for |K| odd; checked) */
        C = malloc(n * sizeof(int)); nC = 0;
        for (int g = 0; g < n; g++) if (cnj(g, tinv) == tinv) C[nC++] = g;
        if ((long)nC * (ninv) != n) die("involutions are not a single class of size n/|C_G(t)|");
    } else {
        C = malloc(n * sizeof(int)); nC = n;
        for (int g = 0; g < n; g++) C[g] = g;
    }
    classes();
    int Z = 0; for (int g = 0; g < n; g++) { int central = 1; for (int s = 0; s < n && central; s++) if (mul(s, g) != mul(g, s)) central = 0; Z += central; }
    fprintf(stderr, "group %s: nK=%d M=%d n=%d involutions=%d centre=%d classes=%d |C|=%d |B|=%d loaded in %.1fs\n", gfile, nK, M, n, ninv, Z, ncls, nC, nB, now() - t_start);
    if (!strcmp(mode, "info")) {
        long hist[8192] = {0};
        for (int g = 0; g < n; g++) { int o = 1, h = g; while (h != 0) { h = mul(h, g); o++; if (o > 8191) die("order"); } hist[o]++; }
        printf("INFO family=%s group=%s nK=%d M=%d n=%d involutions=%d centre=%d classes=%d cp=%.6f |C|=%d |B|=%d orders=", family == 1 ? "cay7" : "lift", gfile, nK, M, n, ninv, Z, ncls, (double)ncls / n, nC, nB);
        for (int o = 1; o < 8192; o++) if (hist[o]) printf("%d:%ld ", o, hist[o]);
        printf("\nspec=%s\n", spec);
        run_exhaustive(1);
    } else if (!strcmp(mode, "eval")) {
        int prof[64] = {0};
        if (family == 1) {
            gens[0] = tinv; gens[1] = set[0]; gens[2] = inv[set[0]]; gens[3] = set[1]; gens[4] = inv[set[1]]; gens[5] = set[2]; gens[6] = inv[set[2]];
            int d = bfs7(-1, prof);
            printf("EVAL family=cay7 t=%d set=[%d,%d,%d] inv=[%d,%d,%d] D=%d defect=%d profile=", tinv, set[0], set[1], set[2], inv[set[0]], inv[set[1]], inv[set[2]], D, d);
            for (int l = 0; l <= D; l++) printf("%d ", prof[l]);
            printf("\n");
        } else {
            la = set[0]; lb = set[1]; lc2 = set[2]; lc3 = set[3];
            int d0 = bfs_lift(0, -1, prof);
            printf("EVAL family=lift set=[%d,%d,%d,%d] D=%d root0 defect=%d profile=", la, lb, lc2, lc3, D, d0);
            for (int l = 0; l <= D; l++) printf("%d ", prof[l]);
            memset(prof, 0, sizeof(prof));
            int d1 = bfs_lift(1, -1, prof);
            printf("| root1 defect=%d profile=", d1);
            for (int l = 0; l <= D; l++) printf("%d ", prof[l]);
            printf("| total defect=%d\n", d0 + d1);
        }
    } else if (!strcmp(mode, "edges")) {
        write_edges(set, efile);
    } else if (!strcmp(mode, "random")) {
        run_random(samples, seed);
    } else if (!strcmp(mode, "exhaustive")) {
        run_exhaustive(dry);
    } else if (!strcmp(mode, "cover")) {
        run_exhaustive(2);
    } else die("mode");
    return 0;
}
