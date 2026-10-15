/*
 * cayley_sa_twoplane.c -- the large-order variant of cayley_sa.c: simulated annealing with targeted moves and an exact
 * targeted completion for Cayley graphs (and 2-vertex lifts) on G = Z_N x|_r Z_M, with a TWO-PLANE in-place BFS so that
 * n up to about 3e10 fits in 8 GB of memory (a three-set BFS needs frontier, next frontier and visited: 3 n bits).
 *
 * Group convention (the same as the verifiers): element x^i y^j, product (i1, j1)(i2, j2) = (i1 + r^j1 i2 mod N,
 * j1 + j2 mod M), edges g ~ g s.  Slots 0..k-1 are the k inverse pairs; with --invol the involution t = x^tx y^(M/2)
 * is added (degree 2k + 1) and printed as the last pair of every state line.
 *
 * Two-plane BFS.  A vertex set is stored as (fibres x M) rows of N bits (row = Z_M-layer y of one fibre).
 * Plane V = visited, plane P = (BFS level is odd) for visited vertices and 0 for unvisited ones.  The
 * frontier of level L is V & (P == L mod 2) up to older vertices of the same parity (levels L-2, L-4, ...),
 * which are harmless: all their neighbours are visited already, and the per-row cost of the rotate-or
 * does not depend on how many bits a row carries.  Level L -> L+1 is written in place, destination row by
 * destination row: tmp = OR over the arcs into this row of rotate(frontier of the source row); new = tmp &
 * ~V[row]; V[row] |= new; if L+1 is odd, P[row] |= new.  The new bits have parity (L+1) mod 2, so a row
 * that is written early in the pass and read as a source later in the same pass does not leak them: the
 * frontier mask selects parity L mod 2.  The last level only counts the unreached vertices row by row, aborts as soon
 * as the count exceeds the bound, and collects up to UMAX = 64 unreached ("far") vertices of an exact evaluation.
 * Right multiplication by s = (sx, sy) maps the layer y onto the layer y + sy by x -> x + r^y sx mod N,
 * i.e. by a CYCLIC ROTATION of the whole row; one BFS level is (#arcs x M) rotate-ors of N-bit rows.
 *
 * usage: cayley_sa_twoplane N M r D k c seconds seed [stall [temp [xonly_pct]]] [--invol] [--xunit] [--tmove pct] [--tx X]
 *        [--ptarget pct] [--J J] [--near T] [--tcomplete J secs] [--nbr R] [--screen R] [--init x1 y1 x2 y2 ...]
 *   c = 0: Cayley graph, k inverse pairs (degree 2k, + 1 with --invol); slots 0..k-1, the involution pair printed last.
 *   seconds = 0 with --init: evaluate the given set once (full sphere profile, far vertices) and exit.
 *   stall: SA restarts after this many evaluations without improvement; temp: Metropolis temperature on the defect;
 *   xonly_pct: share of random moves that change only the x-coordinate of a slot.
 *   --invol / --xunit / --tmove / --tx   as in cayley_sa.c (one involution, unit x-coordinates, involution moves).
 *   --ptarget pct      SA: pct % of the moves are targeted: a far vertex u of the last accepted state, a slot j, random
 *   [--J J]            non-backtracking walks w1, w2 over the arcs of the OTHER slots (t included) with |w1| + |w2| <= J
 *                      (J uniform in 0..J), s' = w1^-1 u w2^-1 or its inverse, resampled until the slot's y is unchanged
 *                      (x-only); u is then at distance <= |w1| + |w2| + 1.  Self-checked (the walk w1 s'^(+-1) w2 ends
 *                      at u) for the first 2000 moves.  t-moves (--tmove) are aimed at a far involution of the last
 *                      accepted state when there is one, random otherwise.
 *   --near T           print every new distinct accepted state with defect <= T ("NEAR d (... s, eval e, from d0, how) :").
 *   --nbr R            baseline for --tcomplete: R random single-slot x-changes (t-changes 1 in K+1) of the --init state, exact.
 *   --tcomplete J secs exact targeted completion of the --init state (defect d <= UMAX, far vertices v, one aim per
 *                      inverse class): for every slot j (and for t) EVERY x-value whose slot element or its inverse is
 *                      L v R for non-backtracking walks L, R over the other arcs with |L| + |R| <= J, deduplicated, with
 *                      the minimal |L| + |R| recorded; the candidates are evaluated in order of increasing level (random
 *                      order within a level) with bound = d, capped at secs seconds; the first HIST are evaluated
 *                      exactly (defect distribution of the neighbourhood); TIE / IMPROVED / FOUND lines.  The exact
 *                      evaluations also check that the aim is no longer far.  Before the sweep the popcounts of the rows
 *                      of the ball B_(D-1) are printed (the density of the weakest targeted condition "v = w s'^(+-1)
 *                      with w in B_(D-1)").
 *   --screen R         R random states, defect distribution (family softness screen).
 * Output lines "BEST defect ..." carry the full generator list; "FOUND" when the defect is 0.
 *
 * Build: gcc -O3 -march=native -o cayley_sa_twoplane cayley_sa_twoplane.c -lm   (portable C, no SIMD intrinsics)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef long long ll;
typedef unsigned long long ull;

static ll N, M, R, *rpow;
static int D, K, C, deg, nfib, nslots, nrows;
static ll W;                 /* words per row */
static ll nn;                /* total vertices */
static uint64_t *V, *P;      /* the two planes, nrows * W words each */
static uint64_t *tmp;        /* one row (W + 2 words), cache resident: the gathered next-level row */
static uint8_t *rowflag_cur, *rowflag_nxt;
static uint64_t lastmask;    /* valid bits of the last word of a row */

/* slots: gx[i], gy[i].  For the lift, link 0 is implicit (trivial). */
static ll gx[32], gy[32];

/* arcs (generator or its inverse, link or its inverse): from fibre ff, to fibre tf, y-step ay, shift table T[y] */
typedef struct { int ff, tf; ll ay, ax; ll *T; } arc_t;
static arc_t arcs[64];
static int narcs;

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static int invol = 0, xunit = 0; static ll tx = 0, tstep = 1, ninv = 1; static unsigned tmove_pct = 0;
static ll gcdll(ll a, ll b);
static inline int x_ok(ll x) { return !xunit || gcdll(x, N) == 1; }
static void random_tx(void) { if (ninv > 1) { ll nt; do nt = tstep * (ll)(rnd() % ninv); while (nt == tx); tx = nt; } }
static inline int tx_ok(ll x) { return x >= 0 && x < N && x % tstep == 0; }

static inline ll inv_x(ll sx, ll sy) { ll jj = (M - sy) % M; return (N - (rpow[jj] * sx) % N) % N; }
static inline ll inv_y(ll sy) { return (M - sy) % M; }
/* scalar group law on (x, y) pairs */
static inline void gmul2(ll ax, ll ay, ll bx, ll by, ll *cx, ll *cy) { *cx = (ax + (rpow[ay] * bx) % N) % N; *cy = (ay + by) % M; }

static void add_arc(int ff, int tf, ll ax, ll ay) {
    arc_t *a = &arcs[narcs++];
    a->ff = ff; a->tf = tf; a->ay = ay; a->ax = ax;
    for (ll y = 0; y < M; y++) a->T[y] = (rpow[y] * ax) % N;
}

static void build_arcs(void) {
    narcs = 0;
    if (C == 0) {
        for (int i = 0; i < K; i++) { add_arc(0, 0, gx[i], gy[i]); add_arc(0, 0, inv_x(gx[i], gy[i]), inv_y(gy[i])); }
        if (invol) add_arc(0, 0, tx, M / 2);   /* t = x^tx y^(M/2), its own inverse */
    } else {
        for (int i = 0; i < K; i++) { add_arc(0, 0, gx[i], gy[i]); add_arc(0, 0, inv_x(gx[i], gy[i]), inv_y(gy[i])); }
        for (int i = K; i < 2 * K; i++) { add_arc(1, 1, gx[i], gy[i]); add_arc(1, 1, inv_x(gx[i], gy[i]), inv_y(gy[i])); }
        add_arc(0, 1, 0, 0); add_arc(1, 0, 0, 0);
        for (int j = 2 * K; j < nslots; j++) { add_arc(0, 1, gx[j], gy[j]); add_arc(1, 0, inv_x(gx[j], gy[j]), inv_y(gy[j])); }
    }
}
/* Cayley mode: arcs 2i, 2i+1 are slot i (s_i, s_i^-1), arc 2K is t (slot K, its own inverse) */
static inline int arc_slot(int a) { return a < 2 * K ? a / 2 : K; }
static inline int inv_arc(int a) { return a < 2 * K ? (a ^ 1) : a; }

/* simplicity: all arcs leaving a fibre must be distinct non-identity group elements per target fibre */
static int simple_ok(void) {
    for (int a = 0; a < narcs; a++) {
        if (arcs[a].ff == arcs[a].tf && arcs[a].ax == 0 && arcs[a].ay == 0) return 0;   /* loop generator = identity */
        for (int b = a + 1; b < narcs; b++)
            if (arcs[a].ff == arcs[b].ff && arcs[a].tf == arcs[b].tf && arcs[a].ax == arcs[b].ax && arcs[a].ay == arcs[b].ay) return 0;
    }
    return 1;
}

#define MW(w) (ps ? (vs[w] & (ps[w] ^ pm)) : vs[w])

static inline void shift_or_m(uint64_t *dst, ll dpos, const uint64_t *vs, const uint64_t *ps, uint64_t pm, ll spos, ll len) {
    if (len <= 0) return;
    ll dend = dpos + len;
    ll dw0 = dpos >> 6, dw1 = (dend - 1) >> 6;
    ll full0 = dw0 + ((dpos & 63) != 0), full1 = dw1 - ((dend & 63) != 0);   /* full destination words */
    for (ll dw = dw0; dw <= dw1; dw++) {
        if (dw == full0 && full1 >= full0) {
            ll q = (dw << 6) - dpos + spos; ll sw = q >> 6; int sb = (int)(q & 63);
            ll cntw = full1 - full0 + 1;
            if (sb == 0) { for (ll t = 0; t < cntw; t++) dst[dw + t] |= MW(sw + t); }
            else { int rb = 64 - sb; for (ll t = 0; t < cntw; t++) dst[dw + t] |= (MW(sw + t) >> sb) | (MW(sw + t + 1) << rb); }
            dw = full1; continue;
        }
        ll lo = dw << 6; if (lo < dpos) lo = dpos;
        ll hi = (dw + 1) << 6; if (hi > dend) hi = dend;
        ll s0 = spos + (lo - dpos);
        ll sw = s0 >> 6; int sb = (int)(s0 & 63);
        uint64_t bits = MW(sw) >> sb;
        if (sb) bits |= MW(sw + 1) << (64 - sb);
        int nb = (int)(hi - lo);
        uint64_t mask = nb == 64 ? ~0ULL : ((1ULL << nb) - 1);
        dst[dw] |= (bits & mask) << (lo & 63);
    }
}

/* dst[(x + c) mod N] |= masked src[x] for all x in [0, N) */
static inline void rot_or_m(uint64_t *dst, const uint64_t *vs, const uint64_t *ps, uint64_t pm, ll c) {
    if (c == 0) { for (ll w = 0; w < W; w++) dst[w] |= MW(w); return; }
    shift_or_m(dst, c, vs, ps, pm, 0, N - c);
    shift_or_m(dst, 0, vs, ps, pm, N - c, c);
}

static ll sphere[64];
static int last_lev;
static double bfs_seconds;   /* wall time of the last bfs() */

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

/* ---------- far vertices ----------
 * The last level collects, when collect_far is set, up to UMAX unreached vertices coded x + N y (fibre 0 only: the
 * targeted options are for Cayley graphs) into cand_far[]; cand_far_ok says the list is complete (exact evaluation
 * with at most UMAX unreached vertices).  far[] / nfar / far_ok hold the list of the last accepted SA state. */
#define UMAX 64
static ll cand_far[UMAX], farv[UMAX]; static int cand_nfar = 0, nfar = 0, collect_far = 0, cand_far_ok = 0, far_ok = 0;
static ll *rowpop = NULL; static int want_rowpop = 0;   /* popcount of every row of V after level D-1 (the ball B_(D-1)) */
static long abort_rows = 0, abort_calls = 0; static int last_aborted = 0;           /* last-level aborts (statistics) */

static inline void collect_bits(uint64_t u, int tf, ll y, ll w) {
    if (tf != 0) return;
    while (u && cand_nfar < UMAX) { int b = __builtin_ctzll(u); u &= u - 1; cand_far[cand_nfar++] = (w << 6) + b + N * y; }
}

/* BFS from (root fibre, identity); returns #vertices beyond distance D in the whole graph (both fibres),
 * or a value > bound as soon as the Moore bound on the remaining levels, or the row-wise count of the last
 * level, shows that the result exceeds `bound`. */
static ll bfs(int rootfib, ll bound) {
    double t0 = now();
    size_t words = (size_t)nrows * W;
    memset(V, 0, words * 8); memset(P, 0, words * 8);
    memset(rowflag_cur, 0, nrows);
    ll rrow = (ll)rootfib * M;
    V[rrow * W] = 1; rowflag_cur[rrow] = 1;          /* level 0 is even: P stays 0 */
    ll reached = 1; sphere[0] = 1;
    cand_nfar = 0; cand_far_ok = 0; last_aborted = 0;
    int lev;
    for (lev = 1; lev <= D; lev++) {
        int fpar = (lev - 1) & 1;                    /* parity of the frontier level */
        int setP = lev & 1;                          /* new vertices are at an odd level: mark them in P */
        int last = (lev == D);                       /* last level: count only */
        memset(rowflag_nxt, 0, nrows);
        ll cnt = 0, unreached = 0;
        if (last) abort_calls++;
        for (int tf = 0; tf < nfib; tf++) {
            for (ll y2 = 0; y2 < M; y2++) {
                ll drow = (ll)tf * M + y2;
                uint64_t *vd = V + drow * W, *pd = P + drow * W;
                int any = 0;
                for (int a = 0; a < narcs; a++) {
                    if (arcs[a].tf != tf) continue;
                    ll y = y2 - arcs[a].ay; if (y < 0) y += M;
                    if (rowflag_cur[(ll)arcs[a].ff * M + y]) { any = 1; break; }
                }
                if (!any) {
                    if (last) {   /* no new vertex can land in this row: every unvisited vertex of the row is unreached */
                        for (ll w = 0; w < W; w++) {
                            uint64_t u = ~vd[w] & (w == W - 1 ? lastmask : ~0ULL);
                            unreached += __builtin_popcountll(u);
                            if (u && collect_far) collect_bits(u, tf, y2, w);
                        }
                        if (unreached > bound) { abort_rows++; last_aborted = 1; last_lev = D; bfs_seconds = now() - t0; return unreached; }
                    }
                    continue;
                }
                memset(tmp, 0, (size_t)(W + 2) * 8);
                for (int a = 0; a < narcs; a++) {
                    if (arcs[a].tf != tf) continue;
                    ll y = y2 - arcs[a].ay; if (y < 0) y += M;
                    ll srow = (ll)arcs[a].ff * M + y;
                    if (!rowflag_cur[srow]) continue;
                    const uint64_t *vs = V + srow * W, *ps = (srow < drow && !last) ? P + srow * W : NULL;
                    rot_or_m(tmp, vs, ps, fpar ? 0 : ~0ULL, arcs[a].T[y]);
                }
                ll c = 0;
                if (last) {
                    for (ll w = 0; w < W; w++) {
                        uint64_t mask = (w == W - 1) ? lastmask : ~0ULL;
                        c += __builtin_popcountll(tmp[w] & ~vd[w]);
                        uint64_t u = ~(vd[w] | tmp[w]) & mask;
                        unreached += __builtin_popcountll(u);
                        if (u && collect_far) collect_bits(u, tf, y2, w);
                    }
                    if (unreached > bound) { abort_rows++; last_aborted = 1; last_lev = D; bfs_seconds = now() - t0; return unreached; }
                } else if (setP) {
                    for (ll w = 0; w < W; w++) { uint64_t b = tmp[w] & ~vd[w]; vd[w] |= b; pd[w] |= b; c += __builtin_popcountll(b); }
                } else {
                    for (ll w = 0; w < W; w++) { uint64_t b = tmp[w] & ~vd[w]; vd[w] |= b; c += __builtin_popcountll(b); }
                }
                if (c) rowflag_nxt[drow] = 1;
                cnt += c;
            }
        }
        sphere[lev] = cnt; reached += cnt;
        if (last) {
            if (unreached != nn - reached) { printf("FAIL: last-level row count %lld != %lld\n", unreached, nn - reached); exit(3); }
            cand_far_ok = (unreached <= UMAX) && (cand_nfar == unreached);
        }
        if (cnt == 0) { lev++; break; }
        if (lev < D) {
            /* Moore bound on what the remaining levels can add */
            double m = (double)cnt, tot = 0;
            for (int i = lev + 1; i <= D; i++) { m *= (deg - 1); tot += m; }
            if ((double)(nn - reached) - tot > (double)bound) { last_lev = lev; bfs_seconds = now() - t0; return nn - reached - (ll)tot; }
        }
        if (want_rowpop && lev == D - 1) {
            for (ll row = 0; row < nrows; row++) { ll c = 0; const uint64_t *vr = V + row * W; for (ll w = 0; w < W; w++) c += __builtin_popcountll(vr[w]); rowpop[row] = c; }
        }
        uint8_t *tfl = rowflag_cur; rowflag_cur = rowflag_nxt; rowflag_nxt = tfl;
    }
    last_lev = lev - 1;
    for (int l = lev; l <= D; l++) sphere[l] = 0;
    bfs_seconds = now() - t0;
    return nn - reached;
}

static ll eval(ll bound, int verbose) {
    build_arcs();
    if (!simple_ok()) return (ll)1 << 50;
    ll total = 0;
    for (int f = 0; f < nfib; f++) {
        ll d = bfs(f, bound - total);
        if (verbose) { printf("    root fibre %d spheres:", f); for (int l = 0; l <= D; l++) printf(" %lld", sphere[l]); printf("  unreached %lld  (%.1f s)\n", d, bfs_seconds); fflush(stdout); }
        total += d;
        if (total > bound) return total;
    }
    return total;
}

static void print_far(const char *what) {
    printf("    %s (%d, x y invol):", what, cand_nfar);
    for (int i = 0; i < cand_nfar; i++) { ll x = cand_far[i] % N, y = cand_far[i] / N; printf("  %lld %lld %d", x, y, invol && y == M / 2 && tx_ok(x)); }
    printf("%s\n", cand_far_ok ? "" : "  (incomplete)");
}

/* ---------- Poisson yardstick on the Z_M quotient (non-backtracking walks of length <= D) ---------- */
static double yardstick(void) {
    int Q = nfib * (int)M;
    static double *cntA = NULL, *cntB = NULL, *Cy = NULL;
    if (!cntA) { cntA = malloc(sizeof(double) * Q * 64); cntB = malloc(sizeof(double) * Q * 64); Cy = malloc(sizeof(double) * Q); }
    int inv[64];
    for (int a = 0; a < narcs; a++) for (int b = 0; b < narcs; b++)
        if (arcs[a].ff == arcs[b].tf && arcs[a].tf == arcs[b].ff && ((arcs[a].ay + arcs[b].ay) % M) == 0 && ((arcs[a].ax + (rpow[arcs[a].ay] * arcs[b].ax) % N) % N) == 0) inv[a] = b;
    double Pp = 0;
    for (int rf = 0; rf < nfib; rf++) {
        for (int q = 0; q < Q; q++) Cy[q] = 0;
        memset(cntA, 0, sizeof(double) * Q * 64);
        for (int a = 0; a < narcs; a++) if (arcs[a].ff == rf) { int q = arcs[a].tf * M + arcs[a].ay % M; cntA[q * 64 + a] += 1; }
        for (int len = 1; len <= D; len++) {
            for (int q = 0; q < Q; q++) for (int a = 0; a < narcs; a++) Cy[q] += cntA[q * 64 + a];
            if (len == D) break;
            memset(cntB, 0, sizeof(double) * Q * 64);
            for (int q = 0; q < Q; q++) {
                int f = q / M; ll y = q % M;
                for (int a = 0; a < narcs; a++) {
                    double c = cntA[q * 64 + a]; if (c == 0) continue;
                    for (int b = 0; b < narcs; b++) {
                        if (arcs[b].ff != f || b == inv[a]) continue;
                        ll y2 = (y + arcs[b].ay) % M; int q2 = arcs[b].tf * M + y2;
                        cntB[q2 * 64 + b] += c;
                    }
                }
            }
            double *t = cntA; cntA = cntB; cntB = t;
        }
        for (int q = 0; q < Q; q++) Pp += (double)N * exp(-Cy[q] / (double)N);
    }
    return Pp;
}

static ll gcdll(ll a, ll b) { while (b) { ll t = a % b; a = b; b = t; } return a; }
static int ypattern_ok(void) {
    ll g = invol ? M / 2 : M;
    for (int i = 0; i < nslots; i++) g = gcdll(g, gy[i]);
    return g == 1;
}

static void random_slot(int i, int full) {
    do gx[i] = rnd() % N; while (!x_ok(gx[i]));
    if (full) {
        if (C == 0 || i < 2 * K) gy[i] = 1 + rnd() % (M - 1);   /* loops: avoid Z_N-elements (y = 0) */
        else gy[i] = rnd() % M;
    }
}

static void print_gens(void) {
    for (int i = 0; i < nslots; i++) printf(" %lld %lld", gx[i], gy[i]);
    if (invol) printf(" %lld %lld", tx, M / 2);
    printf("\n");
}

static void pick_ypattern(int tries) {
    double bestP = 1e300; ll by[32];
    for (int t = 0; t < tries; t++) {
        for (int i = 0; i < nslots; i++) random_slot(i, 1);
        if (!ypattern_ok()) continue;
        build_arcs();
        if (!simple_ok()) continue;
        double Pp = yardstick();
        if (Pp < bestP) { bestP = Pp; for (int i = 0; i < nslots; i++) by[i] = gy[i]; }
    }
    for (int i = 0; i < nslots; i++) { gy[i] = by[i]; gx[i] = rnd() % N; }
    printf("y-pattern chosen from %d samples: yardstick P = %.2f, y =", tries, bestP);
    for (int i = 0; i < nslots; i++) printf(" %lld", gy[i]);
    printf("\n"); fflush(stdout);
}

/* ---------- targeted moves (x-only, Cayley graphs) ---------- */
static int targetJ = 9; static long tgt_checked = 0; static unsigned ptarget_pct = 0;
/* slot < 0: random slot 0..K-1; ucode < 0: random far vertex of the last accepted state.  Returns 1 when a move was made. */
static int targeted_move_at(int slot, ll ucode) {
    if (ucode < 0) { if (!far_ok || nfar == 0) return 0; ucode = farv[rnd() % nfar]; }
    ll ux = ucode % N, uy = ucode / N;
    if (slot < 0) slot = (int)(rnd() % K);
    int dir = (int)(rnd() % 2);
    int others[64], no = 0;
    for (int a = 0; a < narcs; a++) if (arc_slot(a) != slot) others[no++] = a;
    long maxtries = 60L * M;
    for (long tr = 0; tr < maxtries; tr++) {
        int L = (int)(rnd() % (targetJ + 1)); int L1 = (int)(rnd() % (L + 1)), L2 = L - L1;
        ll w1x = 0, w1y = 0, w2x = 0, w2y = 0; int ok = 1, last = -1;
        for (int t = 0; t < L1; t++) {
            int arc = -1;
            for (int q = 0; q < 20; q++) { int c = others[rnd() % no]; if (c != last) { arc = c; break; } }
            if (arc < 0) { ok = 0; break; }
            gmul2(w1x, w1y, arcs[arc].ax, arcs[arc].ay, &w1x, &w1y); last = inv_arc(arc);
        }
        if (!ok) continue;
        last = -1;
        for (int t = 0; t < L2; t++) {
            int arc = -1;
            for (int q = 0; q < 20; q++) { int c = others[rnd() % no]; if (c != last) { arc = c; break; } }
            if (arc < 0) { ok = 0; break; }
            gmul2(w2x, w2y, arcs[arc].ax, arcs[arc].ay, &w2x, &w2y); last = inv_arc(arc);
        }
        if (!ok) continue;
        ll sx, sy, mx, my;
        gmul2(inv_x(w1x, w1y), inv_y(w1y), ux, uy, &mx, &my);
        gmul2(mx, my, inv_x(w2x, w2y), inv_y(w2y), &sx, &sy);
        if (dir) { ll ix = inv_x(sx, sy), iy = inv_y(sy); sx = ix; sy = iy; }
        if (sy != gy[slot]) continue;                           /* x-only: the slot keeps its y */
        if (!x_ok(sx) || sx == gx[slot]) continue;
        if (tgt_checked < 2000) {                              /* self-check: the walk w1 s'^(+-1) w2 ends at u */
            tgt_checked++;
            ll px, py, qx, qy;
            if (dir) { qx = inv_x(sx, sy); qy = inv_y(sy); } else { qx = sx; qy = sy; }
            gmul2(w1x, w1y, qx, qy, &px, &py);
            gmul2(px, py, w2x, w2y, &px, &py);
            if (px != ux || py != uy) { printf("FAIL: targeted move self-check (walk ends at (%lld,%lld), u = (%lld,%lld))\n", px, py, ux, uy); exit(3); }
        }
        gx[slot] = sx;
        return 1;
    }
    return 0;
}

/* involution move: aimed at a far involution of the last accepted state when there is one (returns 1), else random (2) */
static int tmove(void) {
    if (ninv <= 1) return 0;
    ll cands[UMAX]; int nc = 0;
    if (far_ok) for (int i = 0; i < nfar; i++) { ll x = farv[i] % N, y = farv[i] / N; if (y == M / 2 && tx_ok(x) && x != tx) cands[nc++] = x; }
    if (nc > 0) { tx = cands[rnd() % nc]; return 1; }
    random_tx(); return 2;
}

/* ---------- NEAR list (distinct accepted states with defect <= near_T) ---------- */
static long near_T = -1; static ll (*near_list)[33] = NULL; static int near_n = 0;
static int near_seen(void) {
    if (!near_list) near_list = malloc(sizeof(ll) * 33 * 4000);
    for (int t = 0; t < near_n; t++) { int same = near_list[t][32] == tx; for (int i = 0; i < nslots && same; i++) if (near_list[t][2 * i] != gx[i] || near_list[t][2 * i + 1] != gy[i]) same = 0; if (same) return 1; }
    if (near_n < 4000) { for (int i = 0; i < nslots; i++) { near_list[near_n][2 * i] = gx[i]; near_list[near_n][2 * i + 1] = gy[i]; } near_list[near_n][32] = tx; near_n++; }
    return 0;
}

/* ---------- exact targeted completion ---------- */
/* J-bounded ball over the arcs outside skip_slot (non-backtracking walks of length 0..J from e, with repeats) */
static ll *ballx[8], *bally[8]; static int *balll[8]; static long balln[8], ballcap[8];
static void ball_level_alloc(int l, long need) {
    if (ballcap[l] >= need) return;
    ballcap[l] = need + 16; ballx[l] = realloc(ballx[l], ballcap[l] * sizeof(ll)); bally[l] = realloc(bally[l], ballcap[l] * sizeof(ll)); balll[l] = realloc(balll[l], ballcap[l] * sizeof(int));
}
static void build_ball(int skip_slot, int J) {
    ball_level_alloc(0, 1); ballx[0][0] = 0; bally[0][0] = 0; balll[0][0] = -1; balln[0] = 1;
    for (int l = 1; l <= J; l++) {
        ball_level_alloc(l, balln[l - 1] * (long)narcs);
        long n = 0;
        for (long k = 0; k < balln[l - 1]; k++) {
            ll px = ballx[l - 1][k], py = bally[l - 1][k]; int last = balll[l - 1][k];
            for (int a = 0; a < narcs; a++) {
                if (arc_slot(a) == skip_slot) continue;
                if (last >= 0 && a == inv_arc(last)) continue;
                ll cx, cy; gmul2(px, py, arcs[a].ax, arcs[a].ay, &cx, &cy);
                ballx[l][n] = cx; bally[l][n] = cy; balll[l][n] = a; n++;
            }
        }
        balln[l] = n;
    }
}
typedef struct { int slot, level; ll x; } tcand_t;
static tcand_t *tc = NULL; static long ntc = 0, tccap = 0;
static unsigned *tstamp = NULL, tstamp_cur = 0;
static void tc_push(int slot, int level, ll x) {
    if (ntc == tccap) { tccap = tccap ? 2 * tccap : 4096; tc = realloc(tc, tccap * sizeof(tcand_t)); }
    tc[ntc].slot = slot; tc[ntc].level = level; tc[ntc].x = x; ntc++;
}
/* all x-values x' such that slot j := x^x' y^(y_j) (or its inverse) equals L v R with |L| + |R| <= J; slot K = t
 * (x' must be an involution x-coordinate).  Deduplicated per slot across aims by tstamp (one stamp per slot). */
static long tc_enumerate(int j, int J, ll vx, ll vy, long *per_level) {
    build_ball(j, J);
    ll yj = j < K ? gy[j] : M / 2, yjinv = (M - yj) % M, curx = j < K ? gx[j] : tx; long added = 0;
    for (int k = 0; k <= J; k++) for (int l = 0; l <= k; l++) {
        int m = k - l;
        for (long a = 0; a < balln[l]; a++) {
            ll px, py; gmul2(ballx[l][a], bally[l][a], vx, vy, &px, &py);            /* P = L v */
            for (long b = 0; b < balln[m]; b++) {
                ll sx, sy; gmul2(px, py, ballx[m][b], bally[m][b], &sx, &sy);        /* s = L v R */
                ll cx;
                if (sy == yj) cx = sx; else if (sy == yjinv) cx = inv_x(sx, sy); else continue;
                if (j == K) { if (!tx_ok(cx)) continue; } else if (!x_ok(cx)) continue;
                if (j < K && yj == yjinv) { ll ix = inv_x(cx, yj); if (ix < cx) cx = ix; }   /* y_j = M/2: x and its inverse are the same slot */
                if (cx == curx || tstamp[cx] == tstamp_cur) continue;
                if (tgt_checked < 2000) {   /* self-check: L^-1 s'^(+-1) R^-1 = v for the slot element s' = x^cx y^yj */
                    tgt_checked++;
                    ll qx = cx, qy = yj; if (sy != yj) { qx = inv_x(cx, yj); qy = inv_y(yj); }
                    ll wx, wy; gmul2(inv_x(ballx[l][a], bally[l][a]), inv_y(bally[l][a]), qx, qy, &wx, &wy);
                    gmul2(wx, wy, inv_x(ballx[m][b], bally[m][b]), inv_y(bally[m][b]), &wx, &wy);
                    if (wx != vx || wy != vy) { printf("FAIL: tcomplete self-check (slot %d, walk ends at (%lld,%lld), v = (%lld,%lld))\n", j, wx, wy, vx, vy); exit(3); }
                }
                tstamp[cx] = tstamp_cur; tc_push(j, k, cx); per_level[k]++; added++;
            }
        }
    }
    return added;
}
#define HIST 24
static int tcomplete(int J, double secs) {
    collect_far = 1; want_rowpop = 1; rowpop = malloc(sizeof(ll) * nrows);
    double t0 = now();
    ll start = eval((ll)1 << 50, 1);
    want_rowpop = 0;
    printf("start: defect %lld (%.1f s)\nBEST %lld :", start, now() - t0, start); print_gens();
    if (start == 0) { printf("FOUND N=%lld M=%lld r=%lld D=%d k=%d c=%d degree %d n=%lld + involution (last pair) :", N, M, R, D, K, C, deg, nn); print_gens(); return 0; }
    if (!cand_far_ok) { printf("tcomplete: start defect %lld > %d, no complete far-vertex list\n", start, UMAX); return 1; }
    print_far("far vertices");
    /* aims: one far vertex per inverse class */
    ll aimx[UMAX], aimy[UMAX]; int naim = 0;
    for (int i = 0; i < cand_nfar; i++) {
        ll x = cand_far[i] % N, y = cand_far[i] / N, ix = inv_x(x, y), iy = inv_y(y); int seen = 0;
        for (int q = 0; q < naim; q++) if ((aimx[q] == x && aimy[q] == y) || (aimx[q] == ix && aimy[q] == iy)) seen = 1;
        if (!seen) { aimx[naim] = x; aimy[naim] = y; naim++; }
    }
    printf("  aims (one per inverse class):"); for (int q = 0; q < naim; q++) printf("  x^%lld y^%lld%s", aimx[q], aimy[q], aimy[q] == M / 2 && tx_ok(aimx[q]) ? " (involution)" : ""); printf("\n");
    /* density of the weakest targeted condition: s' with v s'^-1 in B_(D-1) (row y_v - y_j of the ball) or v s' in B_(D-1) (row y_v + y_j) */
    printf("  ball B_%d of the start state: %lld vertices (%.1f %% of n); the 'sphere-%d' filter (v = w s'^(+-1), w in B_%d of the FULL set) admits per slot:\n", D - 1, nn - sphere[D] - start, 100.0 * (nn - sphere[D] - start) / nn, D - 1, D - 1);
    for (int q = 0; q < naim; q++) for (int j = 0; j <= (invol ? K : K - 1); j++) {
        ll yj = j < K ? gy[j] : M / 2, r1 = ((aimy[q] - yj) % M + M) % M, r2 = (aimy[q] + yj) % M;
        printf("    aim %d slot %s%d (y %lld): rows %lld / %lld of B_%d have %lld / %lld vertices -> %.1f %% / %.1f %% of the %lld x-values%s\n", q, j < K ? "" : "t=", j, yj, r1, r2, D - 1, rowpop[r1], rowpop[r2], 100.0 * rowpop[r1] / N, 100.0 * rowpop[r2] / N, N, j < K ? "" : " (t: only the involution x-coordinates among them)");
    }
    /* enumeration */
    tstamp = calloc(N, sizeof(unsigned));
    long per_level[8] = {0}; int nsl = invol ? K + 1 : K; long per_slot[33] = {0};
    double te = now();
    for (int j = 0; j < nsl; j++) {
        tstamp_cur++;
        for (int q = 0; q < naim; q++) per_slot[j] += tc_enumerate(j, J, aimx[q], aimy[q], per_level);
    }
    printf("  targeted set J = %d: %ld candidates (%.1f s); per slot:", J, ntc, now() - te);
    for (int j = 0; j < nsl; j++) printf(" %s%d:%ld", j < K ? "" : "t=", j, per_slot[j]);
    printf("; per level |L|+|R|:"); for (int k = 0; k <= J; k++) printf(" %d:%ld", k, per_level[k]); printf("\n");
    for (int l = 1; l <= J; l++) printf("    ball level %d over the arcs of the other slots: %ld walks\n", l, balln[l]);
    fflush(stdout);
    /* order: by level, random within a level */
    for (long i = 1; i < ntc; i++) { tcand_t v = tc[i]; long k = i - 1; while (k >= 0 && tc[k].level > v.level) { tc[k + 1] = tc[k]; k--; } tc[k + 1] = v; }
    for (long i = 0; i < ntc;) { long e = i; while (e < ntc && tc[e].level == tc[i].level) e++; for (long k = e - 1; k > i; k--) { long r = i + (long)(rnd() % (k - i + 1)); tcand_t t = tc[k]; tc[k] = tc[r]; tc[r] = t; } i = e; }
    ll sx[32], stx = tx; for (int i = 0; i < nslots; i++) sx[i] = gx[i];
    long done = 0, ties = 0, improved = 0, nexact = 0, exact_sum = 0; ll exact_vals[HIST]; long hist[8] = {0}, done_level[8] = {0};
    int found = 0; double ts = now();
    for (long i = 0; i < ntc && !found; i++) {
        if (now() - ts > secs) { printf("  time cap %.0f s reached after %ld of %ld candidates (level %d partly done)\n", secs, done, ntc, tc[i].level); break; }
        int j = tc[i].slot;
        if (j < K) gx[j] = tc[i].x; else tx = tc[i].x;
        int exact = nexact < HIST;
        ll val = eval(exact ? (ll)1 << 50 : start, 0);
        done++; done_level[tc[i].level]++;
        if (val == (ll)1 << 50) { printf("  cand %ld: slot %s%d level %d x %lld -> not simple\n", i + 1, j < K ? "" : "t=", j, tc[i].level, tc[i].x); }
        else {
            if (exact) {
                exact_vals[nexact++] = val; exact_sum += val;
                int bkt = val == 0 ? 0 : val <= 2 ? 1 : val <= 4 ? 2 : val <= 8 ? 3 : val <= 16 ? 4 : val <= 32 ? 5 : val <= 64 ? 6 : 7; hist[bkt]++;
                /* end-to-end check: the aim must be reached in the candidate graph */
                if (cand_far_ok) for (int f = 0; f < cand_nfar; f++) for (int q = 0; q < naim; q++) {
                    ll x = cand_far[f] % N, y = cand_far[f] / N;
                    if (x == aimx[q] && y == aimy[q]) { printf("FAIL: aim x^%lld y^%lld still far after the targeted change (slot %d x %lld)\n", x, y, j, tc[i].x); exit(3); }
                }
            }
            printf("  cand %ld: slot %s%d level %d x %lld -> %s%lld (%.1f s, %.0f s)\n", i + 1, j < K ? "" : "t=", j, tc[i].level, tc[i].x, val > start && !exact ? "> " : "", val > start && !exact ? start : val, bfs_seconds, now() - t0);
            if (val <= start) {
                if (val < start) { improved++; printf("IMPROVED %lld (slot %s%d level %d, %.0f s) :", val, j < K ? "" : "t=", j, tc[i].level, now() - t0); }
                else { ties++; printf("TIE %lld (slot %s%d level %d, %.0f s) :", val, j < K ? "" : "t=", j, tc[i].level, now() - t0); }
                print_gens();
                if (val == 0) {
                    printf("FOUND: re-evaluation with full profiles\n");
                    eval((ll)1 << 50, 1);
                    printf("FOUND N=%lld M=%lld r=%lld D=%d k=%d c=%d degree %d n=%lld + involution (last pair) :", N, M, R, D, K, C, deg, nn); print_gens();
                    found = 1;
                }
            }
        }
        fflush(stdout);
        if (!found) { gx[j] = sx[j]; tx = stx; }
    }
    double mean = nexact ? (double)exact_sum / nexact : 0;
    for (long i = 1; i < nexact; i++) { ll v = exact_vals[i]; long k = i - 1; while (k >= 0 && exact_vals[k] > v) { exact_vals[k + 1] = exact_vals[k]; k--; } exact_vals[k + 1] = v; }
    printf("TCOMPLETE done: %ld of %ld candidates evaluated (per level:", done, ntc);
    for (int k = 0; k <= J; k++) printf(" %d:%ld/%ld", k, done_level[k], per_level[k]);
    printf("), %ld ties, %ld improved, %s, %.0f s; last-level aborts %ld of %ld; exact sample of %ld: mean %.2f, min %lld, q25 %lld, median %lld, q75 %lld, max %lld; hist [0]=%ld [<=2]=%ld [<=4]=%ld [<=8]=%ld [<=16]=%ld [<=32]=%ld [<=64]=%ld [>64]=%ld\n",
           ties, improved, found ? "FOUND" : "no hit", now() - t0, abort_rows, abort_calls, nexact, mean, nexact ? exact_vals[0] : -1, nexact ? exact_vals[nexact / 4] : -1, nexact ? exact_vals[nexact / 2] : -1, nexact ? exact_vals[3 * nexact / 4] : -1, nexact ? exact_vals[nexact - 1] : -1,
           hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7]);
    if (!found) { for (int i = 0; i < nslots; i++) gx[i] = sx[i]; tx = stx; }
    printf("BEST %lld :", found ? 0 : start); print_gens(); fflush(stdout);
    return found ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc < 9) { fprintf(stderr, "usage: cayley_sa_twoplane N M r D k c seconds seed [stall [temp [xonly_pct]]] [--invol] [--xunit] [--tmove pct] [--tx X] [--ptarget pct] [--J J] [--near T] [--tcomplete J secs] [--boundtest R] [--nbr R] [--screen R] [--init x1 y1 ...]\n"); return 2; }
    N = atoll(argv[1]); M = atoll(argv[2]); R = atoll(argv[3]); D = atoi(argv[4]); K = atoi(argv[5]); C = atoi(argv[6]);
    int seconds = atoi(argv[7]); rs = (uint64_t)atoll(argv[8]) * 2654435761ULL + 88172645463325252ULL;
    int stall = 150; int ai = 9; double temp = 0.0;
    unsigned xonly_pct = 50;                                                              /* share of x-only moves */
    if (argc > ai && strncmp(argv[ai], "--", 2)) { stall = atoi(argv[ai]); ai++; }
    if (argc > ai && strncmp(argv[ai], "--", 2)) { temp = atof(argv[ai]); ai++; }       /* SA temperature (0 = pure descent) */
    if (argc > ai && strncmp(argv[ai], "--", 2)) { xonly_pct = (unsigned)atoi(argv[ai]); ai++; }
    long screen = 0, boundtest = 0, nbr = 0; int tcJ = -1; double tcsecs = 0; targetJ = D - 1;
    for (;;) {   /* options */
        if (argc > ai && !strcmp(argv[ai], "--invol")) { invol = 1; ai++; continue; }
        if (argc > ai && !strcmp(argv[ai], "--xunit")) { xunit = 1; ai++; continue; }
        if (argc > ai + 1 && !strcmp(argv[ai], "--tmove")) { tmove_pct = (unsigned)atoi(argv[ai + 1]); ai += 2; continue; }
        if (argc > ai + 1 && !strcmp(argv[ai], "--tx")) { tx = atoll(argv[ai + 1]); ai += 2; continue; }
        if (argc > ai + 1 && !strcmp(argv[ai], "--ptarget")) { ptarget_pct = (unsigned)atoi(argv[ai + 1]); ai += 2; continue; }
        if (argc > ai + 1 && !strcmp(argv[ai], "--J")) { targetJ = atoi(argv[ai + 1]); ai += 2; continue; }
        if (argc > ai + 1 && !strcmp(argv[ai], "--near")) { near_T = atol(argv[ai + 1]); ai += 2; continue; }
        if (argc > ai + 2 && !strcmp(argv[ai], "--tcomplete")) { tcJ = atoi(argv[ai + 1]); tcsecs = atof(argv[ai + 2]); ai += 3; continue; }
        if (argc > ai + 1 && !strcmp(argv[ai], "--boundtest")) { boundtest = atol(argv[ai + 1]); ai += 2; continue; }
        if (argc > ai + 1 && !strcmp(argv[ai], "--nbr")) { nbr = atol(argv[ai + 1]); ai += 2; continue; }
        break;
    }
    if (invol && (C != 0 || M % 2)) { printf("FAIL: --invol needs a Cayley graph (c = 0) and M even\n"); return 2; }
    if ((ptarget_pct > 0 || tcJ >= 0) && C != 0) { printf("FAIL: --ptarget / --tcomplete need a Cayley graph (c = 0)\n"); return 2; }
    if (tcJ > 6 || targetJ < 0 || targetJ > 12) { printf("FAIL: --tcomplete J in 0..6, --J in 0..12\n"); return 2; }
    if (argc > ai + 1 && !strcmp(argv[ai], "--screen")) { screen = atol(argv[ai + 1]); ai += 2; }
    if (argc > ai && strcmp(argv[ai], "--init")) { fprintf(stderr, "unexpected argument %s\n", argv[ai]); return 2; }
    nfib = C ? 2 : 1; nslots = C ? 2 * K + C - 1 : K; deg = C ? 2 * K + C : 2 * K; nrows = nfib * (int)M;
    if (invol) deg = 2 * K + 1;
    rpow = malloc(M * sizeof(ll)); rpow[0] = 1 % N;
    for (ll j = 1; j < M; j++) { rpow[j] = (rpow[j - 1] * R) % N; if (rpow[j] == 1 % N) { printf("FAIL: r has order %lld < M\n", j); return 1; } }
    if ((rpow[M - 1] * R) % N != 1 % N) { printf("FAIL: r^M != 1 mod N\n"); return 1; }
    for (int a = 0; a < 64; a++) arcs[a].T = malloc(M * sizeof(ll));
    W = (N + 63) / 64; nn = (ll)nrows * N; lastmask = (N % 64) ? ((1ULL << (N % 64)) - 1) : ~0ULL;
    if (invol) {
        ll h = rpow[M / 2]; ninv = gcdll(N, (1 + h) % N); tstep = N / ninv;
        if (tx % tstep) { printf("FAIL: --tx %lld is not an involution x-coordinate (0 mod %lld)\n", tx, tstep); return 2; }
        printf("involution t = x^%lld y^%lld (r^(M/2) = %lld mod N): %lld involutions x^i y^(M/2) (i = 0 mod %lld)%s; x units only: %d; t-moves %u %%\n",
               tx, M / 2, h, ninv, tstep, h == 1 % N ? " (CENTRAL: bad)" : "", xunit, tmove_pct);
    }
    size_t words = (size_t)nrows * W + 2;
    V = calloc(words, 8); P = calloc(words, 8);
    tmp = calloc((size_t)W + 2, 8);
    rowflag_cur = calloc(nrows, 1); rowflag_nxt = calloc(nrows, 1);
    if (!V || !P || !tmp || !rowflag_cur || !rowflag_nxt) { printf("alloc failed\n"); return 2; }
    printf("G = Z_%lld x|_%lld Z_%lld (|G| = %lld), %s, degree %d, D = %d, n = %lld, %d slots, %d s, seed %s, stall %d; two planes of %.2f GB each; targeted moves %u %% (J %d), near threshold %ld%s\n",
           N, R, M, N * M, C ? "2-vertex lift" : "Cayley graph", deg, D, nn, nslots, seconds, argv[8], stall, words * 8 / 1e9, ptarget_pct, targetJ, near_T, tcJ >= 0 ? ", targeted completion mode" : "");
    double moore = 1, m = 1; for (int i = 1; i <= D; i++) { m *= (i == 1 ? deg : deg - 1); moore += m; }
    printf("Moore ball %.0f, Moore/n = %.3f, n e^(-Moore/n) = %.2f\n", moore, moore / nn, nn * exp(-moore / nn));
    int have_init = 0;
    if (argc > ai && !strcmp(argv[ai], "--init")) {
        if (argc < ai + 1 + 2 * nslots) { printf("FAIL: --init needs %d pairs\n", nslots); return 2; }
        for (int i = 0; i < nslots; i++) { gx[i] = atoll(argv[ai + 1 + 2 * i]) % N; gy[i] = atoll(argv[ai + 2 + 2 * i]) % M; }
        have_init = 1; printf("initial set from the command line\n");
    }
    if (boundtest > 0) {   /* validation: random states in the --init pattern; bounded verdicts must agree with the exact value */
        if (!have_init) { printf("FAIL: --boundtest needs --init\n"); return 2; }
        long bad = 0, tested = 0, aborted0 = abort_rows; collect_far = 1;
        for (long t = 0; t < boundtest; t++) {
            if (t > 0) { for (int i = 0; i < nslots; i++) random_slot(i, 0); if (invol && ninv > 1 && (t & 1)) random_tx(); }
            build_arcs(); if (!ypattern_ok() || !simple_ok()) continue;
            ll full = eval((ll)1 << 50, 0); int nf = cand_nfar, fok = cand_far_ok; ll fl[UMAX]; memcpy(fl, cand_far, sizeof fl);
            if ((full <= UMAX) != (fok && nf == full)) { bad++; printf("FAR MISMATCH state %ld: defect %lld, %d collected, ok %d\n", t, full, nf, fok); }
            for (ll b = 0; b <= 48; b += (b < 8 ? 1 : 8)) {
                ll r = eval(b, 0); tested++;
                if ((r <= b) != (full <= b) || (r <= b && r != full)) { bad++; printf("BOUND MISMATCH state %ld bound %lld: exact %lld, bounded %lld\n", t, b, full, r); }
                if (r <= b && full <= UMAX && (cand_nfar != nf || memcmp(fl, cand_far, sizeof(ll) * nf))) { bad++; printf("FAR LIST MISMATCH state %ld bound %lld\n", t, b); }
            }
            printf("  state %ld: defect %lld, %d far vertices\n", t, full, nf);
        }
        printf("boundtest: %ld verdicts on %ld states, %ld mismatches, %ld last-level aborts\n", tested, boundtest, bad, abort_rows - aborted0);
        return bad ? 1 : 0;
    }
    if (nbr > 0) {   /* baseline for --tcomplete: R random single-slot x-changes (and t-changes, 1 in K+1) of the --init state, exact defects */
        if (!have_init) { printf("FAIL: --nbr needs --init\n"); return 2; }
        build_arcs(); double t0 = now(); ll start = eval((ll)1 << 50, 1);
        printf("start: defect %lld\nBEST %lld :", start, start); print_gens(); fflush(stdout);
        ll sx[32], stx = tx; for (int i = 0; i < nslots; i++) sx[i] = gx[i];
        ll vals[4096]; long got = 0, hist[8] = {0}; double sum = 0;
        while (got < nbr && got < 4096) {
            int j = (int)(rnd() % (invol && ninv > 1 ? K + 1 : K));
            if (j < K) random_slot(j, 0); else random_tx();
            build_arcs(); if (!simple_ok()) { gx[j < K ? j : 0] = sx[j < K ? j : 0]; tx = stx; continue; }
            ll v = eval((ll)1 << 50, 0); vals[got++] = v; sum += v;
            int bkt = v == 0 ? 0 : v <= 2 ? 1 : v <= 4 ? 2 : v <= 8 ? 3 : v <= 16 ? 4 : v <= 32 ? 5 : v <= 64 ? 6 : 7; hist[bkt]++;
            printf("  nbr %ld: slot %s%d x %lld -> %lld (%.1f s, %.0f s)\n", got, j < K ? "" : "t=", j, j < K ? gx[j] : tx, v, bfs_seconds, now() - t0);
            if (v <= start) { printf("%s %lld (random neighbour, slot %s%d) :", v < start ? "IMPROVED" : "TIE", v, j < K ? "" : "t=", j); print_gens(); }
            if (v == 0) { printf("FOUND N=%lld M=%lld r=%lld D=%d k=%d c=%d degree %d n=%lld + involution (last pair) :", N, M, R, D, K, C, deg, nn); print_gens(); fflush(stdout); return 0; }
            fflush(stdout);
            if (j < K) gx[j] = sx[j]; else tx = stx;
        }
        for (long i = 1; i < got; i++) { ll v = vals[i]; long k = i - 1; while (k >= 0 && vals[k] > v) { vals[k + 1] = vals[k]; k--; } vals[k + 1] = v; }
        printf("NBR done: %ld random single-slot neighbours of the defect-%lld state, %.0f s: mean %.2f, min %lld, q25 %lld, median %lld, q75 %lld, max %lld; hist [0]=%ld [<=2]=%ld [<=4]=%ld [<=8]=%ld [<=16]=%ld [<=32]=%ld [<=64]=%ld [>64]=%ld\n",
               got, start, now() - t0, got ? sum / got : 0, got ? vals[0] : -1, got ? vals[got / 4] : -1, got ? vals[got / 2] : -1, got ? vals[3 * got / 4] : -1, got ? vals[got - 1] : -1, hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7]);
        return 1;
    }
    if (tcJ >= 0) {
        if (!have_init) { printf("FAIL: --tcomplete needs --init\n"); return 2; }
        build_arcs();
        if (!ypattern_ok() || !simple_ok()) { printf("FAIL: start state is not a simple connected candidate\n"); return 1; }
        printf("yardstick P of the start pattern: %.2f\n", yardstick());
        return tcomplete(tcJ, tcsecs);
    }
    time_t t0 = time(NULL);
    if (screen > 0) {
        ll *def = malloc(sizeof(ll) * screen); long got = 0; double sum = 0; ll mn = (ll)1 << 50; ll mnx[32], mny[32];
        while (got < screen) {
            for (int i = 0; i < nslots; i++) random_slot(i, 1);
            if (invol) random_tx();
            build_arcs();
            if (!ypattern_ok() || !simple_ok()) continue;
            ll d = eval((ll)1 << 50, 0);
            def[got++] = d; sum += (double)d;
            printf("  screen %ld: defect %lld (%ld s):", got, d, (long)(time(NULL) - t0)); print_gens(); fflush(stdout);
            if (d < mn) { mn = d; for (int i = 0; i < nslots; i++) { mnx[i] = gx[i]; mny[i] = gy[i]; } }
            if (d == 0) { printf("FOUND N=%lld M=%lld r=%lld D=%d k=%d c=%d degree %d n=%lld%s :", N, M, R, D, K, C, deg, nn, invol ? " + involution (last pair)" : ""); print_gens(); fflush(stdout); }
        }
        long le8 = 0, le16 = 0, le32 = 0, le64 = 0;
        for (long i = 0; i < got; i++) { le8 += def[i] <= 8; le16 += def[i] <= 16; le32 += def[i] <= 32; le64 += def[i] <= 64; }
        for (long i = 1; i < got; i++) { ll v = def[i]; long j = i - 1; while (j >= 0 && def[j] > v) { def[j + 1] = def[j]; j--; } def[j + 1] = v; }
        printf("SCREEN N=%lld M=%lld r=%lld n=%lld samples=%ld time=%ld s: min %lld, median %lld, mean %.1f, q25 %lld, q75 %lld, <=8 %.2f%%, <=16 %.2f%%, <=32 %.2f%%, <=64 %.2f%%\n",
               N, M, R, nn, got, (long)(time(NULL) - t0), mn, def[got / 2], sum / got, def[got / 4], def[(3 * got) / 4],
               100.0 * le8 / got, 100.0 * le16 / got, 100.0 * le32 / got, 100.0 * le64 / got);
        printf("  best set:"); for (int i = 0; i < nslots; i++) printf(" %lld %lld", mnx[i], mny[i]); printf("\n");
        return 0;
    }
    if (!have_init) pick_ypattern(400);
    build_arcs();
    printf("yardstick P of the start pattern: %.2f\n", yardstick());
    collect_far = (ptarget_pct > 0 || tmove_pct > 0 || near_T >= 0 || seconds == 0);
    ll best = eval((ll)1 << 50, 1);
    if (collect_far && cand_far_ok) print_far("far vertices");
    memcpy(farv, cand_far, sizeof farv); nfar = cand_nfar; far_ok = cand_far_ok;
    ll bx[32], by[32], btx = tx; for (int i = 0; i < nslots; i++) { bx[i] = gx[i]; by[i] = gy[i]; }
    ll gbest = best; ll gbx[32], gby[32], gbtx = tx; for (int i = 0; i < nslots; i++) { gbx[i] = gx[i]; gby[i] = gy[i]; }
    printf("start: defect %lld (%ld s)\n", best, (long)(time(NULL) - t0)); printf("BEST %lld :", best); print_gens(); fflush(stdout);
    if (near_T >= 0 && best <= near_T) near_seen();
    long evals = 1, accepted = 0, since_improve = 0, kicks = 0, ntgt = 0, acc_tgt = 0, ntm = 0, acc_tm = 0, ntm_aimed = 0;
    double Pcur = yardstick();
    while (gbest > 0 && (ll)(time(NULL) - t0) < seconds) {
        ll ox[32], oy[32], otx = tx; for (int i = 0; i < nslots; i++) { ox[i] = gx[i]; oy[i] = gy[i]; }
        /* spend more evaluations before kicking when the defect is already small (the last steps are the hardest) */
        long stall_eff = stall * (best <= 2 ? 6 : best <= 4 ? 3 : best <= 8 ? 2 : 1);
        if (since_improve >= stall_eff) {
            kicks++; since_improve = 0;
            if (rnd() % 2) { for (int i = 0; i < nslots; i++) { gx[i] = bx[i]; gy[i] = by[i]; } tx = btx; }
            else { for (int i = 0; i < nslots; i++) { gx[i] = gbx[i]; gy[i] = gby[i]; } tx = gbtx; }
            if (invol && tmove_pct > 0) random_tx();
            if (kicks % 3 == 0) { pick_ypattern(400); printf("kick %ld: new y-pattern\n", kicks); }
            else { int nk = 2 + (int)(rnd() % 2); for (int t = 0; t < nk; t++) random_slot((int)(rnd() % nslots), (int)(rnd() % 2)); printf("kick %ld: %d slots re-randomised\n", kicks, nk); }
            build_arcs();
            if (!ypattern_ok() || !simple_ok()) { for (int i = 0; i < nslots; i++) { gx[i] = ox[i]; gy[i] = oy[i]; } tx = otx; since_improve = stall; continue; }
            best = eval((ll)1 << 50, 0); evals++;
            memcpy(farv, cand_far, sizeof farv); nfar = cand_nfar; far_ok = cand_far_ok;
            for (int i = 0; i < nslots; i++) { bx[i] = gx[i]; by[i] = gy[i]; } btx = tx;
            Pcur = yardstick();
            printf("  after kick: defect %lld (P %.2f), %ld s\n", best, Pcur, (long)(time(NULL) - t0)); fflush(stdout);
            continue;
        }
        int targeted = 0, tmoved = 0;
        if (invol && tmove_pct > 0 && (unsigned)(rnd() % 100) < tmove_pct) { tmoved = tmove(); if (tmoved) { ntm++; ntm_aimed += tmoved == 1; } }
        if (!tmoved && ptarget_pct > 0 && (unsigned)(rnd() % 100) < ptarget_pct) targeted = targeted_move_at(-1, -1);
        if (!tmoved && !targeted) {
            unsigned u = (unsigned)(rnd() % 100);
            if (u < xonly_pct) random_slot((int)(rnd() % nslots), 0);
            else if (u < xonly_pct + (100 - xonly_pct) * 3 / 5) random_slot((int)(rnd() % nslots), 1);
            else { int a = (int)(rnd() % nslots), b = (int)(rnd() % nslots); random_slot(a, (int)(rnd() % 2)); random_slot(b, (int)(rnd() % 2)); }
        }
        build_arcs();
        int ok = ypattern_ok() && simple_ok();
        if (ok) {
            int ychanged = 0; for (int i = 0; i < nslots; i++) if (gy[i] != oy[i]) ychanged = 1;
            if (ychanged) { double Pp = yardstick(); if (Pp > 2.0 * Pcur + 5) ok = 0; }
        }
        if (!ok) { for (int i = 0; i < nslots; i++) { gx[i] = ox[i]; gy[i] = oy[i]; } tx = otx; continue; }
        /* SA: a worse set is accepted with probability exp(-(val - best) / temp); bound the BFS accordingly */
        ll slack = 0;
        if (temp > 0) { double uu = (rnd() % 1000000 + 0.5) / 1000000.0; slack = (ll)(-temp * log(uu)); }
        ll val = eval(best + slack, 0); evals++; ntgt += targeted;
        if (val <= best + slack) {
            ll prev = best;
            accepted++; acc_tgt += targeted; acc_tm += tmoved != 0;
            memcpy(farv, cand_far, sizeof farv); nfar = cand_nfar; far_ok = cand_far_ok;
            const char *how = targeted ? ", targeted" : tmoved == 1 ? ", t aimed" : tmoved == 2 ? ", t random" : "";
            if (val < best) {
                since_improve = 0;
                printf("eval %ld (%ld s): defect %lld -> %lld%s\n", evals, (long)(time(NULL) - t0), best, val, how);
                printf("BEST %lld :", val); print_gens(); fflush(stdout);
            } else since_improve++;
            best = val; for (int i = 0; i < nslots; i++) { bx[i] = gx[i]; by[i] = gy[i]; } btx = tx;
            if (near_T >= 0 && val <= near_T && !near_seen()) { printf("NEAR %lld (%ld s, eval %ld, from %lld%s) :", val, (long)(time(NULL) - t0), evals, prev, how); print_gens(); fflush(stdout); }
            Pcur = yardstick();
            if (val < gbest) { gbest = val; for (int i = 0; i < nslots; i++) { gbx[i] = gx[i]; gby[i] = gy[i]; } gbtx = tx; }
        } else {
            since_improve++;
            for (int i = 0; i < nslots; i++) { gx[i] = ox[i]; gy[i] = oy[i]; } tx = otx;
            if (evals % 10 == 0) { printf("  eval %ld (%ld s): current defect %lld, last candidate %s%lld%s\n", evals, (long)(time(NULL) - t0), best, last_aborted ? "> " : "", last_aborted ? best + slack : val, targeted ? " (targeted)" : tmoved ? " (t-move)" : ""); fflush(stdout); }
        }
    }
    for (int i = 0; i < nslots; i++) { gx[i] = gbx[i]; gy[i] = gby[i]; } tx = gbtx;
    printf("done: %ld evaluations, %ld accepted (%ld of %ld targeted, %ld of %ld t-moves (%ld aimed)), %ld kicks, %ld s; global best defect %lld; last-level aborts %ld of %ld; %d near states\n",
           evals, accepted, acc_tgt, ntgt, acc_tm, ntm, ntm_aimed, kicks, (long)(time(NULL) - t0), gbest, abort_rows, abort_calls, near_n);
    printf("BEST %lld :", gbest); print_gens();
    if (gbest == 0) {
        printf("FOUND: re-evaluation with full profiles\n");
        eval((ll)1 << 50, 1);
        printf("FOUND N=%lld M=%lld r=%lld D=%d k=%d c=%d degree %d n=%lld%s :", N, M, R, D, K, C, deg, nn, invol ? " + involution (last pair)" : ""); print_gens();
    }
    return gbest == 0 ? 0 : 1;
}
