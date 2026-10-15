/*
 * cayley_sa.c -- simulated annealing with level patterns, targeted moves and exact completion sweeps for Cayley graphs
 * and 2-vertex voltage lifts on the metacyclic group G = Z_N x|_r Z_M = <x, y | x^N = y^M = 1, y x y^-1 = x^r>.
 *
 * Conventions: element x^i y^j has index i + N j; product (i1, j1)(i2, j2) = (i1 + r^j1 i2 mod N, j1 + j2 mod M);
 * edges g ~ g s for s in the connection set S u S^-1.  Slots 0..k-1 are the loops at base vertex 0 (Cayley graph: the k
 * inverse pairs), k..2k-1 the loops at base vertex 1, 2k..2k+c-2 the links 1..c-1 (link 0 carries the trivial voltage);
 * --init and --ypat list the slots in that order as "x y" pairs resp. y-values.  The level of x^i y^j is j (its
 * y-coordinate); a level pattern is the list of the slots' levels.  The defect of a state is the number of vertices at
 * distance > D from the root(s); a defect-0 state is a graph of diameter <= D ("FOUND").
 *
 * Evaluator: bit-parallel BFS over rows of N bits (one row per Z_M layer of a fibre); right multiplication by x^a y^b
 * maps row y to row y + b rotated by r^y a.  Small levels are expanded by scattering set bits, large levels destination
 * row by destination row, and the last level is only counted, row by row, with early exit once the unreached count
 * exceeds the bound.  Left multiplication by G is an automorphism transitive on each fibre, so one BFS per fibre decides
 * the diameter.
 *
 * usage: cayley_sa N M r D k c [--init x1 y1 ...] [options]
 *   c = 0: Cayley graph with k inverse pairs (degree 2k; 2k + 1 with --invol); c > 0: 2-vertex lift with k loops per base
 *   vertex and c links (degree 2k + c).  r = 0: the smallest r of multiplicative order M mod N is used.
 *
 * Odd degree (Cayley graphs only, M even):
 *   --invol        add the involution t = x^tx y^(M/2) to the connection set (degree 2K + 1).  G has
 *                  I = gcd(N, 1 + r^(M/2)) involutions x^i y^(M/2), i = 0 mod N/I, all conjugate; the centraliser of t
 *                  has index I.  Every printed state ends with the pair "tx M/2" after the K slot pairs, so a FOUND line
 *                  is the whole connection set.
 *   --tx X         start value of tx (must be an involution x-coordinate).
 *   --tmove pct    SA: pct % of the moves change tx, aimed at an unreached involution of the last accepted state when
 *                  there is one (that vertex becomes a neighbour of e), random otherwise.
 *   --xunit        every slot's x-coordinate must be a unit mod N (no generator commutes with t when N = p q).
 *   --xunitI       weaker: x must only be a unit mod I (no generator commutes with t; x = 0 mod a centralised prime is
 *                  allowed, which some known graphs need).
 *   --qyard Q      diagnostic: Poisson yardstick on the quotient Z_Q x| Z_M (Q | N), printed with --eval.
 *
 * Level patterns:
 *   --patgen K     sample K level patterns, print "PAT coinc P ratio y..." where coinc is the number of radius-<= 2
 *                  coincidences of the Z_M quotient (reduced words of length <= 2 from each root fibre landing on the
 *                  same quotient vertex, or on the root or a radius-1 vertex, beyond the forced commutations of loops in
 *                  one fibre), P the Poisson yardstick (expected unreached count from the pattern alone, summed over the
 *                  root fibres), ratio = P / its value for a uniform pattern.
 *   --patopt R     R restarts of a local search for a coincidence-free pattern of minimal yardstick.
 *   --ypat y1 ..   fixed level pattern (x-coordinates random, or from --init); implies --fixy.
 *   --fixy / --freey   SA never / may change a y-coordinate.
 *   --maxcoinc c   reject start patterns with more than c coincidences.
 *   --sample S --ypat ...   defect quantiles of S random x-assignments in the pattern (family softness screen).
 *
 * Simulated annealing:
 *   --sa secs T xonly_pct   Metropolis on the defect for secs seconds at temperature T (slack ~ T * Exp(1)); xonly_pct %
 *                  of the random moves change only the x-coordinate of a slot.
 *   --ptarget pct  share of targeted moves: the new slot value is s' = w1^-1 u w2^-1 for an unreached vertex u of the last
 *                  accepted state and short non-backtracking walks w1, w2 over the other arcs (|w1| + |w2| <= J,
 *                  --J J, default D - 1), so that u gets a walk of length |w1| + |w2| + 1 through the new generator.
 *   --walkT T2     plateau walk: temperature T2 once a state with defect <= the --near threshold has been accepted.
 *   --near T       print every distinct state with defect <= T ("NEAR").
 *   --complete T --compl_cap s   on every new distinct state with defect <= T run the exact x-only single-slot sweep of
 *                  all slots (bound = that defect: FOUND / IMPROVED / TIE lines), capped at s seconds each.  The first
 *                  improved state stops the sweep and is handed back to the SA (nested completion).
 *   --seed s       random seed.
 *
 * Exact completion sweeps around a given state (--init required):
 *   single-slot sweep (default when neither --sa nor --sweep2 is given): every value of every slot (--slots i,j,..;
 *                  --xonly: x-coordinates only; --yfrom y; --order yardstick [--from a] [--to b]: (slot, y) units in
 *                  increasing yardstick order; --seconds s cap; --thresh t: scalar/vector BFS threshold).
 *   --sweep2 secs modes   (Cayley graphs with --invol) bounded sweeps of 2-slot / t neighbourhoods aimed at the far
 *                  vertices of the start state (defect d), each mode capped at secs seconds; the first HIST candidates of
 *                  each mode are evaluated exactly (defect distribution of the neighbourhood), the rest with bound = d
 *                  (TIE / IMPROVED / FOUND lines).  Modes:
 *                    A  all I values of tx (S fixed)
 *                    B  tx := v.x for the first far vertex v, then every unit x of every slot
 *                    C  every tx, then for each slot j one targeted x_j aimed at v
 *                    D  every ordered slot pair (i, j): every unit x_i, then one targeted x_j aimed at v
 *                    E  every ordered slot pair (i, j): every unit x_i, then EVERY x_j whose slot element (or its
 *                       inverse) is L v R for walks L, R over the other arcs with |L| + |R| <= J (--J2 J, default 2);
 *                       the set is enumerated exactly and deduplicated; the aim cycles through the far vertices.
 *                    F  the complete 2-slot x-neighbourhood: every unordered slot pair (or the pairs of --pairs i-j,..),
 *                       all admissible x_i times all admissible x_j.  With --s2all the start values are not skipped and
 *                       with --s2bound b the bound is b instead of the start defect, so mode F is a COMPLETE sweep of
 *                       the whole family {level pattern, t, the fixed slots} from an arbitrary start state.
 *                  Every IMPROVED state met in a sweep gets its complete x-only neighbourhood at once.
 *   --s2steps lo hi   mode F runs only the outer-slot steps [lo, hi) of its pseudo-random cyclic order
 *                  x_i = x0 + step * stride (0 <= step < N); the chunks of a partition of [0, N) partition the family
 *                  exactly when every chunk uses the same --seed.
 *   --x1mod        mode F restricts the outer slot's x to 3 residue classes mod q = N/I plus the classes with a non-unit
 *                  orbit step (a cheap quotient by the residual automorphisms x -> x^u, y -> x^(I m) y; every
 *                  isomorphism class is still met, about 4/q of the states are evaluated).
 *   --quot         mode F with exactly one --pairs pair evaluates exactly one state per orbit of the stabiliser H of t and
 *                  the fixed slots inside {x -> x^u, y -> x^c y}, computed by brute force over (u, c); the candidate
 *                  count equals the number of H-orbits of the family (Burnside).  Preferred over --x1mod.
 *   --noprerow / --prerow_max b   the involution-row pre-filter (after level D-2 the row y = M/2 is finished with the
 *                  radius-2 ball and the candidate is rejected if that row alone has more than the bound unreached
 *                  vertices; exact) off / with bound b.
 *   --eval         evaluate the --init state (full sphere profile, far vertices) and exit.
 *   --selftest     compare the bit-parallel BFS with a plain queue BFS on random states.
 *   --prerowtest R check the pre-filter's verdicts against exact evaluations on R random states.
 *
 * Output lines: "BEST d : x1 y1 ..." (new best of the SA), "NEAR d ...", "TIE d (...) :", "IMPROVED d (...) :" and
 * "FOUND N=.. M=.. r=.. D=.. k=.. c=.. degree .. n=.. [+ involution (last pair)] : x1 y1 x2 y2 ..." for a defect-0 state.
 *
 * Build: gcc -O3 -march=native -o cayley_sa cayley_sa.c -lm   (portable C; no SIMD intrinsics; -march=native is optional)
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

static int N, M, D, K, C, deg, nfib, nslots, nrows, W;
static ll R, *rpow, nn;
static u64 lastmask;
static int gx[32], gy[32];

typedef struct { int ff, tf, ay, ax; int *T; } arc_t;
static arc_t arcs[64]; static int narcs;
static int inv_arc[64];          /* inverse arc of each arc (the involution arc is its own inverse) */
static int out_arcs[2][64], nout[2], in_arcs[2][64], nin[2];
static int invol = 0, xunit = 0; /* --invol: fixed involution generator y^(M/2); --xunit: slot x must be a unit mod N */
static int tx = 0, tstep = 1, ninv = 1;   /* involution t = x^tx y^(M/2), tx = 0 mod tstep = N/ninv (ninv involutions) */
static unsigned tmove_pct = 0;            /* --tmove pct: share of SA moves that change tx */
static inline int tx_ok(int x) { return x >= 0 && x < N && x % tstep == 0; }
static int qyard_Q = 0;          /* --qyard Q: quotient yardstick on Z_Q x| Z_M */

static u64 *vis, *cur, *nxt;
static uint8_t *rowflag_cur, *rowflag_nxt;
static unsigned *rowstamp; static unsigned stamp_counter = 0;
static long sphere[64];
static long scalar_thresh;
static long nbfs = 0;

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

static inline int inv_x(int sx, int sy) { int jj = (M - sy) % M; return (int)((N - (rpow[jj] * sx) % N) % N); }
static inline int inv_y(int sy) { return (M - sy) % M; }

static void set_arc(int a, int ff, int tf, int ax, int ay) {
    arcs[a].ff = ff; arcs[a].tf = tf; arcs[a].ax = ax; arcs[a].ay = ay;
    for (int y = 0; y < M; y++) arcs[a].T[y] = (int)((rpow[y] * ax) % N);
}
static int slot_is_loop(int i) { return C == 0 || i < 2 * K; }
static int slot_fibre(int i) { return (C == 0 || i < K) ? 0 : 1; }
/* arcs of slot i: indices 2i, 2i+1 (loop: s, s^-1 in its fibre; link: 0->1 with voltage s, 1->0 with s^-1); the
 * trivial link's arcs are narcs-2 (0->1) and narcs-1 (1->0).  The inverse of arc a is a ^ 1 in every case. */
static void build_slot_arcs(int i) {
    int f = slot_fibre(i);
    if (slot_is_loop(i)) { set_arc(2 * i, f, f, gx[i], gy[i]); set_arc(2 * i + 1, f, f, inv_x(gx[i], gy[i]), inv_y(gy[i])); }
    else { set_arc(2 * i, 0, 1, gx[i], gy[i]); set_arc(2 * i + 1, 1, 0, inv_x(gx[i], gy[i]), inv_y(gy[i])); }
}
static int arc_slot(int a) { return a < 2 * nslots ? a / 2 : -1; }
static void build_arcs(void) {
    narcs = 2 * nslots;
    for (int i = 0; i < nslots; i++) build_slot_arcs(i);
    for (int a = 0; a < narcs; a++) inv_arc[a] = a ^ 1;
    if (C > 0) { set_arc(narcs, 0, 1, 0, 0); set_arc(narcs + 1, 1, 0, 0, 0); inv_arc[narcs] = narcs + 1; inv_arc[narcs + 1] = narcs; narcs += 2; }
    if (invol) { set_arc(narcs, 0, 0, tx, M / 2); inv_arc[narcs] = narcs; narcs += 1; }
    nout[0] = nout[1] = nin[0] = nin[1] = 0;
    for (int a = 0; a < narcs; a++) { out_arcs[arcs[a].ff][nout[arcs[a].ff]++] = a; in_arcs[arcs[a].tf][nin[arcs[a].tf]++] = a; }
}
static int igcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }
/* rebuild only the involution arc (index 2 nslots in Cayley mode) after a change of tx */
static void build_t_arc(void) { set_arc(2 * nslots, 0, 0, tx, M / 2); }
/* --xunit: a slot's x-coordinate must be a unit mod N (for N = p q with --invol: no generator commutes with t) */
/* --xunitI (xunit = 2): only gcd(x, ninv) = 1 is required, i.e. no generator commutes with t (I | x <=> t a = a t);
 * x may be 0 mod a centralised prime (the known (7,7) graph on 56952 vertices has x = 0 mod 3 in two slots, 3 centralised) */
static inline int x_ok(int x) { return !xunit || (xunit == 2 ? igcd(x, ninv) == 1 : igcd(x, N) == 1); }

/* simplicity: no identity loop arc, no two arcs with the same (from, to, voltage) */
static int simple_ok(void) {
    for (int a = 0; a < narcs; a++) {
        if (arcs[a].ff == arcs[a].tf && arcs[a].ax == 0 && arcs[a].ay == 0) return 0;
        for (int b = a + 1; b < narcs; b++)
            if (arcs[a].ff == arcs[b].ff && arcs[a].tf == arcs[b].tf && arcs[a].ax == arcs[b].ax && arcs[a].ay == arcs[b].ay) return 0;
    }
    return 1;
}
static int simple_ok_slot(int i) {
    for (int a = 2 * i; a <= 2 * i + 1; a++) {
        if (arcs[a].ff == arcs[a].tf && arcs[a].ax == 0 && arcs[a].ay == 0) return 0;
        for (int b = 0; b < narcs; b++) {
            if (b == a) continue;
            if (arcs[a].ff == arcs[b].ff && arcs[a].tf == arcs[b].tf && arcs[a].ax == arcs[b].ax && arcs[a].ay == arcs[b].ay) return 0;
        }
    }
    return 1;
}
static int ypattern_ok(void) { int g = invol ? M / 2 : M; for (int i = 0; i < nslots; i++) g = igcd(g, gy[i]); return g == 1; }
static int xpattern_ok(void) { for (int i = 0; i < nslots; i++) if (!x_ok(gx[i])) return 0; return 1; }

/* dst |= rotate(src, t): bit x of the N-bit row src goes to bit (x + t) mod N; 0 <= t < N */
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
    unsigned st = ++stamp_counter;
    for (int s = 0; s < nrows; s++) {
        if (!rowflag_cur[s]) continue;
        int f = s / M, y = s - f * M, na = nout[f];
        u64 *dp[64]; int tt[64];
        for (int q = 0; q < na; q++) {
            arc_t *a = &arcs[out_arcs[f][q]];
            int y2 = y + a->ay; if (y2 >= M) y2 -= M;
            int d = a->tf * M + y2;
            if (rowstamp[d] != st) { memset(nxt + (size_t)d * W, 0, (size_t)W * 8); rowstamp[d] = st; }
            dp[q] = nxt + (size_t)d * W; tt[q] = a->T[y];
        }
        const u64 *src = cur + (size_t)s * W;
        for (int w = 0; w < W; w++) {
            u64 bits = src[w];
            while (bits) {
                int b = __builtin_ctzll(bits); bits &= bits - 1;
                int x = (w << 6) + b;
                for (int q = 0; q < na; q++) { int x2 = x + tt[q]; if (x2 >= N) x2 -= N; dp[q][x2 >> 6] |= 1ULL << (x2 & 63); }
            }
        }
    }
    long cnt = 0;
    for (int d = 0; d < nrows; d++) {
        if (rowstamp[d] != st) { rowflag_nxt[d] = 0; continue; }
        u64 *nx = nxt + (size_t)d * W, *vs = vis + (size_t)d * W, any = 0;
        for (int w = 0; w < W; w++) { u64 v = nx[w] & ~vs[w]; nx[w] = v; vs[w] |= v; cnt += __builtin_popcountll(v); any |= v; }
        rowflag_nxt[d] = any != 0;
    }
    return cnt;
}

static long level_bp(void) {
    long cnt = 0;
    for (int d = 0; d < nrows; d++) {
        int f = d / M, y = d - f * M, any_src = 0;
        u64 *nx = nxt + (size_t)d * W;
        for (int q = 0; q < nin[f]; q++) {
            arc_t *a = &arcs[in_arcs[f][q]];
            int ys = y - a->ay; if (ys < 0) ys += M;
            int s = a->ff * M + ys;
            if (!rowflag_cur[s]) continue;
            if (!any_src) { memset(nx, 0, (size_t)W * 8); any_src = 1; }
            rot_or(nx, cur + (size_t)s * W, a->T[ys]);
        }
        if (!any_src) { rowflag_nxt[d] = 0; continue; }
        u64 *vs = vis + (size_t)d * W, any = 0;
        for (int w = 0; w < W; w++) { u64 v = nx[w] & ~vs[w]; nx[w] = v; vs[w] |= v; cnt += __builtin_popcountll(v); any |= v; }
        rowflag_nxt[d] = any != 0;
    }
    return cnt;
}

/* ---------- uncovered vertices (for targeted moves) ----------
 * The last level collects, when collect_unc is set, up to UMAX vertices that stay unreached (fibre, x, y),
 * coded fibre * (N M) + x + N y, into cand_unc[root fibre]; they belong to the candidate being evaluated and are
 * copied to unc[] when the candidate is accepted.  Only exact evaluations (unreached <= bound) complete the list. */
#define UMAX 64
static ll cand_unc[2][UMAX], unc[2][UMAX];
static int cand_nunc[2], nunc[2];
static int collect_unc = 0, cur_root = 0;

static inline void collect_row(u64 bits, int f, int y, int w) {
    while (bits && cand_nunc[cur_root] < UMAX) {
        int b = __builtin_ctzll(bits); bits &= bits - 1;
        cand_unc[cur_root][cand_nunc[cur_root]++] = (ll)f * ((ll)N * M) + (ll)((w << 6) + b) + (ll)N * y;
    }
}

/* last level: count the vertices that stay unreached, row by row; stop as soon as the count exceeds bound.
 * *newcnt = number of vertices at distance exactly D (only when not stopped). */
static long level_last(ll bound, long *newcnt) {
    long unreached = 0, cnt = 0;
    u64 *buf = nxt;
    if (collect_unc) cand_nunc[cur_root] = 0;
    for (int d = 0; d < nrows; d++) {
        int f = d / M, y = d - f * M, any_src = 0;
        for (int q = 0; q < nin[f]; q++) {
            arc_t *a = &arcs[in_arcs[f][q]];
            int ys = y - a->ay; if (ys < 0) ys += M;
            int s = a->ff * M + ys;
            if (!rowflag_cur[s]) continue;
            if (!any_src) { memset(buf, 0, (size_t)W * 8); any_src = 1; }
            rot_or(buf, cur + (size_t)s * W, a->T[ys]);
        }
        const u64 *vs = vis + (size_t)d * W;
        if (!any_src) {
            for (int w = 0; w < W; w++) {
                u64 u = ~vs[w] & (w == W - 1 ? lastmask : ~0ULL);
                unreached += __builtin_popcountll(u);
                if (u && collect_unc) collect_row(u, f, y, w);
            }
        } else {
            for (int w = 0; w < W; w++) {
                u64 mask = (w == W - 1) ? lastmask : ~0ULL;
                cnt += __builtin_popcountll(buf[w] & ~vs[w]);
                u64 u = ~(vs[w] | buf[w]) & mask;
                unreached += __builtin_popcountll(u);
                if (u && collect_unc) collect_row(u, f, y, w);
            }
        }
        if (unreached > bound) { *newcnt = -1; return unreached; }
    }
    *newcnt = cnt;
    return unreached;
}

/* involution-row pre-filter.  In a one-involution Cayley graph (Cayley mode, --invol) the row y = M/2 holds the
 * ninv involutions, which are reached at half rate (pairing argument: an involution is reached by a walk w only together with t w), so for a random candidate it
 * carries a Poisson(ninv e^(-Moore/2n)) number of the unreached vertices (a few at the cells this was used for), while an
 * ordinary row carries ~0.3.  After level D-2 the vertices of that row within distance D are exactly
 * (B_(D-2) . B_2) n row = U_{b in B_2} rot(vis[row M/2 - b.y], r^(M/2 - b.y) b.x): |B_2| = 1 + deg + deg^2 rotations
 * instead of the two remaining full levels (2 deg M rotations each).  When the row alone has more than `bound`
 * unreached vertices the candidate is rejected there (exact: those vertices are unreached in the whole graph).
 * Used only when bound <= prerow_max (completion sweeps); the SA's bounds are too loose for it to pay. */
static int prerow_max = 4, prerow_on = 1; static long prerow_rejects = 0, prerow_calls = 0;
static int b2x[1 + 64 + 64 * 64], b2y[1 + 64 + 64 * 64], nb2 = 0;
static void build_b2(void) {   /* B_2 of the current connection set (with repeats; inverse pairs a a^-1 = e skipped) */
    nb2 = 0; b2x[nb2] = 0; b2y[nb2] = 0; nb2++;
    for (int a = 0; a < narcs; a++) { b2x[nb2] = arcs[a].ax; b2y[nb2] = arcs[a].ay; nb2++; }
    for (int a = 0; a < narcs; a++) for (int b = 0; b < narcs; b++) {
        if (inv_arc[a] == b) continue;
        int ay = arcs[a].ay, cy = ay + arcs[b].ay; if (cy >= M) cy -= M;
        b2x[nb2] = (int)((arcs[a].ax + (rpow[ay] * (ll)arcs[b].ax) % N) % N); b2y[nb2] = cy; nb2++;
    }
}
/* number of vertices of row `row` (fibre 0) beyond distance D, given vis = B_(D-2); stops counting above bound */
static long prerow_count(int row, ll bound) {
    u64 *buf = nxt;   /* scratch: nxt is free between levels (every level clears the rows it writes) */
    memset(buf, 0, (size_t)W * 8);
    for (int b = 0; b < nb2; b++) {
        int ys = row - b2y[b]; if (ys < 0) ys += M;
        rot_or(buf, vis + (size_t)ys * W, (int)((rpow[ys] * (ll)b2x[b]) % N));
    }
    long unreached = 0;
    for (int w = 0; w < W; w++) {
        u64 u = ~buf[w] & (w == W - 1 ? lastmask : ~0ULL);
        unreached += __builtin_popcountll(u);
        if (unreached > bound) break;
    }
    return unreached;
}

/* BFS from (fibre rootfib, identity); returns the number of vertices at distance > D in the whole graph, exactly
 * when that number is <= bound, otherwise some value > bound. */
static long bfs(int rootfib, ll bound) {
    nbfs++;
    cur_root = rootfib;
    if (collect_unc) cand_nunc[rootfib] = 0;
    size_t words = (size_t)nrows * W;
    memset(vis, 0, words * 8);
    memset(rowflag_cur, 0, nrows);
    int rrow = rootfib * M;
    memset(cur + (size_t)rrow * W, 0, (size_t)W * 8);
    cur[(size_t)rrow * W] = 1; vis[(size_t)rrow * W] = 1; rowflag_cur[rrow] = 1;
    long reached = 1, cntprev = 1; sphere[0] = 1;
    int prerow = prerow_on && invol && C == 0 && bound <= prerow_max && D >= 4;
    for (int lev = 1; lev <= D; lev++) {
        if (lev == D) { long nc; long un = level_last(bound, &nc); sphere[D] = nc; return un; }
        long cnt = cntprev < scalar_thresh ? level_scalar() : level_bp();
        sphere[lev] = cnt; reached += cnt; cntprev = cnt;
        if (cnt == 0) { for (int l = lev + 1; l <= D; l++) sphere[l] = 0; return nn - reached; }   /* disconnected: no uncovered sample */
        u64 *t = cur; cur = nxt; nxt = t;
        uint8_t *tf = rowflag_cur; rowflag_cur = rowflag_nxt; rowflag_nxt = tf;
        if (prerow && lev == D - 2) {
            prerow_calls++; build_b2();
            long u = prerow_count(M / 2, bound);
            if (u > bound) { prerow_rejects++; for (int l = lev + 1; l <= D; l++) sphere[l] = -1; return u; }
        }
    }
    return nn - reached;
}

static void list_unreached(int rootfib) {
    size_t words = (size_t)nrows * W;
    memset(vis, 0, words * 8); memset(rowflag_cur, 0, nrows);
    int rrow = rootfib * M;
    memset(cur + (size_t)rrow * W, 0, (size_t)W * 8);
    cur[(size_t)rrow * W] = 1; vis[(size_t)rrow * W] = 1; rowflag_cur[rrow] = 1;
    long cntprev = 1;
    for (int lev = 1; lev <= D; lev++) {
        long cnt = cntprev < scalar_thresh ? level_scalar() : level_bp();
        cntprev = cnt; if (cnt == 0) break;
        u64 *t = cur; cur = nxt; nxt = t;
        uint8_t *tf = rowflag_cur; rowflag_cur = rowflag_nxt; rowflag_nxt = tf;
    }
    printf("    unreached from root fibre %d (fibre x y):", rootfib);
    for (int d = 0; d < nrows; d++) {
        int f = d / M, y = d - f * M;
        for (int w = 0; w < W; w++) {
            u64 bits = ~vis[(size_t)d * W + w] & (w == W - 1 ? lastmask : ~0ULL);
            while (bits) { int b = __builtin_ctzll(bits); bits &= bits - 1; printf("  %d %d %d", f, (w << 6) + b, y); }
        }
    }
    printf("\n");
}

static long eval(ll bound, int verbose) {
    long total = 0;
    for (int f = 0; f < nfib; f++) {
        long d = bfs(f, bound - total);
        if (verbose) { printf("    root fibre %d spheres:", f); for (int l = 0; l <= D; l++) printf(" %ld", sphere[l]); printf("  unreached %ld\n", d); fflush(stdout); }
        total += d;
        if (total > bound) return total;
    }
    return total;
}

/* the K slot pairs followed by the involution pair "tx M/2" (the whole connection set, verifier order) */
static void print_gens(void) { for (int i = 0; i < nslots; i++) printf(" %d %d", gx[i], gy[i]); if (invol) printf(" %d %d", tx, M / 2); printf("\n"); }

/* ---------- Poisson yardstick on the Z_M quotient: expected unreached count from the y-pattern alone ---------- */
static double yardstick(void) {
    int Q = nfib * M;
    static double *cntA = NULL, *cntB = NULL, *Cy = NULL;
    if (!cntA) { cntA = malloc(sizeof(double) * Q * 64); cntB = malloc(sizeof(double) * Q * 64); Cy = malloc(sizeof(double) * Q); }
    double P = 0;
    for (int rf = 0; rf < nfib; rf++) {
        for (int q = 0; q < Q; q++) Cy[q] = 0;
        memset(cntA, 0, sizeof(double) * Q * 64);
        for (int a = 0; a < narcs; a++) if (arcs[a].ff == rf) { int q = arcs[a].tf * M + arcs[a].ay % M; cntA[q * 64 + a] += 1; }
        for (int len = 1; len <= D; len++) {
            for (int q = 0; q < Q; q++) for (int a = 0; a < narcs; a++) Cy[q] += cntA[q * 64 + a];
            if (len == D) break;
            memset(cntB, 0, sizeof(double) * Q * 64);
            for (int q = 0; q < Q; q++) {
                int f = q / M, y = q % M;
                for (int a = 0; a < narcs; a++) {
                    double c = cntA[q * 64 + a]; if (c == 0) continue;
                    for (int b = 0; b < narcs; b++) {
                        if (arcs[b].ff != f || b == inv_arc[a]) continue;
                        int y2 = (y + arcs[b].ay) % M, q2 = arcs[b].tf * M + y2;
                        cntB[q2 * 64 + b] += c;
                    }
                }
            }
            double *t = cntA; cntA = cntB; cntB = t;
        }
        for (int q = 0; q < Q; q++) P += (double)N * exp(-Cy[q] / (double)N);
    }
    return P;
}
static double moore_ball(void) { double m = 1, tot = 1; for (int i = 1; i <= D; i++) { m *= (i == 1 ? deg : deg - 1); tot += m; } return tot; }
static double uniform_yardstick(void) { return nfib * (double)nn * exp(-moore_ball() / (double)nn); }

/* Poisson yardstick on the quotient Z_Q x| Z_M (Q | N; x -> x mod Q), fibres of size N/Q: non-backtracking walk
 * counts of length <= D per quotient element (tracked per last arc), P_Q = sum (N/Q) exp(-count / (N/Q)).  Q = 1 is
 * yardstick(); Q = p (N = p q, --invol) sees the involution class and the centraliser of t.  Cayley graphs only. */
static double qyardstick(int Q) {
    if (nfib != 1 || N % Q) return -1;
    size_t V = (size_t)Q * M;
    float *cA = calloc(V * narcs, sizeof(float)), *cB = calloc(V * narcs, sizeof(float)); double *Cy = calloc(V, sizeof(double));
    if (!cA || !cB || !Cy) { printf("qyardstick: alloc failed\n"); exit(2); }
    int *rq = malloc(M * sizeof(int)); for (int j = 0; j < M; j++) rq[j] = (int)(rpow[j] % Q);
    for (int a = 0; a < narcs; a++) { int i = arcs[a].ax % Q, j = arcs[a].ay; cA[((size_t)i * M + j) * narcs + a] += 1; }
    for (int len = 1; len <= D; len++) {
        for (size_t v = 0; v < V; v++) for (int a = 0; a < narcs; a++) Cy[v] += cA[v * narcs + a];
        if (len == D) break;
        memset(cB, 0, V * narcs * sizeof(float));
        for (size_t v = 0; v < V; v++) {
            int i = (int)(v / M), j = (int)(v % M);
            for (int a = 0; a < narcs; a++) {
                float c = cA[v * narcs + a]; if (c == 0) continue;
                for (int b = 0; b < narcs; b++) {
                    if (b == inv_arc[a]) continue;
                    int i2 = (int)((i + (ll)rq[j] * arcs[b].ax) % Q), j2 = j + arcs[b].ay; if (j2 >= M) j2 -= M;
                    cB[((size_t)i2 * M + j2) * narcs + b] += c;
                }
            }
        }
        float *t = cA; cA = cB; cB = t;
    }
    double fib = (double)N / Q, P = 0;
    for (size_t v = 0; v < V; v++) P += fib * exp(-Cy[v] / fib);
    free(cA); free(cB); free(Cy); free(rq);
    return P;
}

/* Radius-<= 2 coincidences of the quotient (both root fibres): reduced words of length 1 or 2 whose quotient endpoints
 * coincide, or hit the root / a radius-1 vertex, beyond the forced commutations ab = ba of two loops in one fibre.
 * 0 = "full sphere" (quotient profile 1, deg, deg(deg-1) - (commutations)). */
static int sphere2_coinc(int *s2_out) {
    static uint8_t *mark = NULL; if (!mark) mark = malloc((size_t)nfib * M);
    int coinc = 0, s2 = 0;
    for (int rf = 0; rf < nfib; rf++) {
        memset(mark, 0, (size_t)nfib * M);
        mark[rf * M] = 1;                       /* root */
        int F1 = 0, E1 = 0;
        for (int a = 0; a < narcs; a++) if (arcs[a].ff == rf) { int q = arcs[a].tf * M + arcs[a].ay; F1++; if (!mark[q]) { mark[q] = 1; E1++; } }
        coinc += F1 - E1;
        /* mark value 1 = distance <= 1; now the 2-words */
        int F2 = 0, E2 = 0;
        for (int a = 0; a < narcs; a++) {
            if (arcs[a].ff != rf) continue;
            for (int b = 0; b < narcs; b++) {
                if (arcs[b].ff != arcs[a].tf || b == inv_arc[a]) continue;
                int loops_same = (arcs[a].ff == arcs[a].tf) && (arcs[b].ff == arcs[b].tf);
                if (loops_same && a > b) continue;      /* ab = ba in the quotient: count the unordered pair once */
                F2++;
                int q = arcs[b].tf * M + (arcs[a].ay + arcs[b].ay) % M;
                if (mark[q] == 0) { mark[q] = 2; E2++; }
            }
        }
        coinc += F2 - E2; s2 += E2;
    }
    if (s2_out) *s2_out = s2;
    return coinc;
}

/* ---------- plain scalar BFS for the self-test (independent of the bit-parallel code) ---------- */
static ll gmul(ll g, ll a) { ll i1 = g % N, j1 = g / N, i2 = a % N, j2 = a / N; return (i1 + (rpow[j1] * i2) % N) % N + N * ((j1 + j2) % M); }
static uint8_t *plain_dist = NULL;
static long plain_bfs(int rootfib, long *sph) {
    ll G = (ll)N * M; ll n = nfib * G;
    if (!plain_dist) plain_dist = malloc(n);
    uint8_t *dist = plain_dist; memset(dist, 0xff, n);
    ll *qa = malloc(n * sizeof(ll)), *qb = malloc(n * sizeof(ll));
    qa[0] = rootfib * G; dist[qa[0]] = 0; ll ca = 1, reached = 1; sph[0] = 1;
    for (int lev = 1; lev <= D; lev++) {
        ll cb = 0;
        for (ll i = 0; i < ca; i++) {
            ll v = qa[i]; int f = (int)(v / G); ll g = v % G;
            for (int q = 0; q < nout[f]; q++) {
                arc_t *a = &arcs[out_arcs[f][q]];
                ll w = a->tf * G + gmul(g, a->ax + (ll)N * a->ay);
                if (dist[w] == 0xff) { dist[w] = (uint8_t)lev; qb[cb++] = w; }
            }
        }
        sph[lev] = cb; reached += cb; ll *t = qa; qa = qb; qb = t; ca = cb;
    }
    free(qa); free(qb);
    return n - reached;
}

static u64 rs = 88172645463325252ULL;
static inline u64 rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

/* scalar group law on (x, y) pairs */
static inline void gmul2(ll ax, ll ay, ll bx, ll by, ll *cx, ll *cy) { *cx = (ax + (rpow[ay] * bx) % N) % N; *cy = (ay + by) % M; }

/* Targeted move (fibre-aware): pick an uncovered vertex u of root fibre fr, a slot, a direction of its arc
 * (a -> b), random non-backtracking walks w1 (fibre fr -> a) and w2 (b -> fibre of u) over the arcs of the OTHER
 * slots with |w1| + |w2| <= targetJ, and set the slot to s' = w1^-1 u w2^-1 (inverted when the arc is used
 * backwards).  keep_y: resample until the slot's y-coordinate is unchanged.  Returns 0 if nothing could be aimed at. */
static int targetJ = 6;
static long tgt_checked = 0;
/* slot < 0 = random slot; ucode < 0 = random uncovered vertex of a random root fibre (else the given vertex code,
 * root fibre 0: Cayley graphs only). */
static int targeted_move_at(int keep_y, int slot, ll ucode) {
    int fr = 0;
    if (ucode < 0) {
        int roots[2], nr = 0;
        for (int f = 0; f < nfib; f++) if (nunc[f] > 0) roots[nr++] = f;
        if (nr == 0) return 0;
        fr = roots[rnd() % nr];
        ucode = unc[fr][rnd() % nunc[fr]];
    }
    ll NM = (ll)N * M;
    int fu = (int)(ucode / NM); ll uidx = ucode % NM; ll ux = uidx % N, uy = uidx / N;
    if (slot < 0) slot = (int)(rnd() % nslots);
    int a, b;
    if (slot_is_loop(slot)) a = b = slot_fibre(slot); else { a = 0; b = 1; }
    int dir = (int)(rnd() % 2);
    if (dir) { int t = a; a = b; b = t; }
    int others[64], no = 0;
    for (int i = 0; i < narcs; i++) if (arc_slot(i) != slot) others[no++] = i;
    long maxtries = keep_y ? 60L * M : 400;
    for (long tr = 0; tr < maxtries; tr++) {
        int L = (int)(rnd() % (targetJ + 1)); int L1 = (int)(rnd() % (L + 1)), L2 = L - L1;
        ll w1x = 0, w1y = 0, w2x = 0, w2y = 0; int f = fr, ok = 1, last = -1;
        for (int t = 0; t < L1; t++) {
            int arc = -1;
            for (int q = 0; q < 20; q++) { int c = others[rnd() % no]; if (arcs[c].ff == f && c != last) { arc = c; break; } }
            if (arc < 0) { ok = 0; break; }
            gmul2(w1x, w1y, arcs[arc].ax, arcs[arc].ay, &w1x, &w1y); f = arcs[arc].tf; last = inv_arc[arc];
        }
        if (!ok || f != a) continue;
        f = b; last = -1;
        for (int t = 0; t < L2; t++) {
            int arc = -1;
            for (int q = 0; q < 20; q++) { int c = others[rnd() % no]; if (arcs[c].ff == f && c != last) { arc = c; break; } }
            if (arc < 0) { ok = 0; break; }
            gmul2(w2x, w2y, arcs[arc].ax, arcs[arc].ay, &w2x, &w2y); f = arcs[arc].tf; last = inv_arc[arc];
        }
        if (!ok || f != fu) continue;
        ll sx, sy, mx, my;
        gmul2(inv_x((int)w1x, (int)w1y), inv_y((int)w1y), ux, uy, &mx, &my);
        gmul2(mx, my, inv_x((int)w2x, (int)w2y), inv_y((int)w2y), &sx, &sy);
        if (dir) { ll ix = inv_x((int)sx, (int)sy), iy = inv_y((int)sy); sx = ix; sy = iy; }
        if (keep_y && sy != gy[slot]) continue;
        if (slot_is_loop(slot) && sy == 0) continue;            /* loops stay out of Z_N */
        if (!x_ok((int)sx)) continue;
        if (sx == gx[slot] && sy == gy[slot]) continue;
        if (tgt_checked < 2000) {                              /* self-check: the walk w1 s'^(+-1) w2 ends at u */
            tgt_checked++;
            ll px, py, qx, qy;
            if (dir) { qx = inv_x((int)sx, (int)sy); qy = inv_y((int)sy); } else { qx = sx; qy = sy; }
            gmul2(w1x, w1y, qx, qy, &px, &py);
            gmul2(px, py, w2x, w2y, &px, &py);
            if (px != ux || py != uy) { printf("FAIL: targeted move self-check (walk ends at (%lld,%lld), u = (%lld,%lld))\n", px, py, ux, uy); exit(3); }
        }
        gx[slot] = (int)sx; gy[slot] = (int)sy;
        return 1;
    }
    return 0;
}
static int targeted_move(int keep_y) { return targeted_move_at(keep_y, -1, -1); }

/* involution move.  If the last accepted state's uncovered sample holds an involution (y = M/2, x = 0 mod tstep)
 * other than the current t, set tx to its x-coordinate (it becomes a neighbour of e); otherwise a random tx.  Returns
 * 1 if targeted, 2 if random, 0 if impossible (ninv = 1). */
static int invol_move(void) {
    if (ninv <= 1) return 0;
    int cands[UMAX], nc = 0;
    for (int i = 0; i < nunc[0]; i++) {
        ll code = unc[0][i]; ll idx = code % ((ll)N * M); int ux = (int)(idx % N), uy = (int)(idx / N);
        if (uy == M / 2 && tx_ok(ux) && ux != tx) cands[nc++] = ux;
    }
    if (nc > 0) { tx = cands[rnd() % nc]; build_t_arc(); return 1; }
    int nt; do nt = tstep * (int)(rnd() % ninv); while (nt == tx);
    tx = nt; build_t_arc(); return 2;
}

static int selftest(int trials) {
    int bad = 0;
    for (int t = 0; t < trials; t++) {
        for (int i = 0; i < nslots; i++) { do gx[i] = (int)(rnd() % N); while (!x_ok(gx[i])); gy[i] = slot_is_loop(i) ? 1 + (int)(rnd() % (M - 1)) : (int)(rnd() % M); }
        if (invol) tx = tstep * (int)(rnd() % ninv);     /* a random involution each trial */
        build_arcs();
        if (!ypattern_ok() || !simple_ok()) { t--; continue; }
        for (int f = 0; f < nfib; f++) {
            long sp[64], sq[64];
            long ref = plain_bfs(f, sp);
            collect_unc = 1;
            long got = bfs(f, (ll)1 << 50);
            collect_unc = 0;
            memcpy(sq, sphere, sizeof sq);
            int ok = got == ref; for (int l = 0; l <= D; l++) if (sq[l] != sp[l]) ok = 0;
            if (!ok) { bad++; printf("MISMATCH trial %d fibre %d: plain %ld bp %ld; spheres", t, f, ref, got); for (int l = 0; l <= D; l++) printf(" %ld/%ld", sp[l], sq[l]); printf("\n"); }
            /* uncovered sample: every collected vertex is unreached in the plain BFS, and the count is right */
            ll G = (ll)N * M;
            for (int i = 0; i < cand_nunc[f]; i++) {
                ll code = cand_unc[f][i]; int ff = (int)(code / G); ll idx = code % G;
                if (plain_dist[ff * G + idx] != 0xff) { bad++; printf("UNCOVERED MISMATCH trial %d fibre %d: vertex (%d, %lld, %lld) has distance %d\n", t, f, ff, idx % N, idx / N, plain_dist[ff * G + idx]); }
            }
            if ((ref < UMAX && cand_nunc[f] != ref) || (ref >= UMAX && cand_nunc[f] != UMAX)) { bad++; printf("UNCOVERED COUNT trial %d fibre %d: %d collected, %ld unreached\n", t, f, cand_nunc[f], ref); }
            for (int k = 0; k < 4; k++) {
                ll b = (ll)(rnd() % (2 * (ref + 2)));
                long gb = bfs(f, b);
                if ((gb <= b && gb != ref) || (gb > b && ref <= b)) { bad++; printf("BOUND MISMATCH trial %d fibre %d bound %lld: exact %ld got %ld\n", t, f, b, ref, gb); }
            }
        }
        long save = scalar_thresh; long sp0[64], sp1[64];
        scalar_thresh = 0; long a0 = bfs(0, (ll)1 << 50); memcpy(sp0, sphere, sizeof sp0);
        scalar_thresh = (long)1 << 40; long a1 = bfs(0, (ll)1 << 50); memcpy(sp1, sphere, sizeof sp1);
        scalar_thresh = save;
        if (a0 != a1 || memcmp(sp0, sp1, sizeof(long) * (D + 1))) { bad++; printf("PATH MISMATCH trial %d: bp-only %ld scalar-only %ld\n", t, a0, a1); }
        /* targeted moves: 50 moves from this state, each self-checked inside targeted_move */
        collect_unc = 1; eval((ll)1 << 50, 0); collect_unc = 0;
        memcpy(unc, cand_unc, sizeof unc); memcpy(nunc, cand_nunc, sizeof nunc);
        int sx[32], sy[32]; memcpy(sx, gx, sizeof gx); memcpy(sy, gy, sizeof gy);
        int nt = 0, nk = 0;
        for (int q = 0; q < 50; q++) {
            int keep = q & 1;
            if (targeted_move(keep)) {
                nt++;
                int ych = 0; for (int i = 0; i < nslots; i++) if (gy[i] != sy[i]) ych = 1;
                if (keep && ych) { bad++; printf("KEEP-Y MISMATCH trial %d\n", t); }
                nk += keep;
            }
            memcpy(gx, sx, sizeof gx); memcpy(gy, sy, sizeof gy);
        }
        build_arcs();
        printf("  trial %d: defects per fibre ok, %d targeted moves (%d keep-y), coinc %d\n", t, nt, nk, sphere2_coinc(NULL));
    }
    printf("selftest: %d trials, %d mismatches\n", trials, bad);
    return bad;
}

/* ---------- random patterns ---------- */
static void random_slot(int i, int full) {
    do gx[i] = (int)(rnd() % N); while (!x_ok(gx[i]));
    if (full) gy[i] = slot_is_loop(i) ? 1 + (int)(rnd() % (M - 1)) : (int)(rnd() % M);
}
static double pick_ypattern(int tries) {
    double bestP = 1e300; int by[32];
    for (int t = 0; t < tries; t++) {
        for (int i = 0; i < nslots; i++) random_slot(i, 1);
        if (!ypattern_ok()) continue;
        build_arcs();
        if (!simple_ok()) continue;
        double P = yardstick();
        if (P < bestP) { bestP = P; for (int i = 0; i < nslots; i++) by[i] = gy[i]; }
    }
    for (int i = 0; i < nslots; i++) { gy[i] = by[i]; do gx[i] = (int)(rnd() % N); while (!x_ok(gx[i])); }
    build_arcs();
    printf("y-pattern chosen from %d samples: yardstick P = %.2f, y =", tries, bestP);
    for (int i = 0; i < nslots; i++) printf(" %d", gy[i]);
    printf("\n"); fflush(stdout);
    return bestP;
}
static void patgen(long K) {
    double U = uniform_yardstick();
    printf("uniform yardstick %.3f (Moore ball %.0f, n %lld, %d root%s)\n", U, moore_ball(), nn, nfib, nfib > 1 ? "s" : "");
    long kept = 0;
    for (long t = 0; t < K; t++) {
        for (int i = 0; i < nslots; i++) random_slot(i, 1);
        if (!ypattern_ok()) continue;
        build_arcs();
        if (!simple_ok()) continue;
        int s2; int co = sphere2_coinc(&s2);
        double P = yardstick();
        printf("PAT %d %.3f %.4f :", co, P, P / U); for (int i = 0; i < nslots; i++) printf(" %d", gy[i]); printf("\n");
        kept++;
    }
    printf("patgen: %ld patterns printed of %ld samples\n", kept, K);
}

/* --patopt R: R hill-climbing restarts on the coincidence count (one level changed at a time, first improvement),
 * printing each local minimum as "PAT coinc P ratio : y..." (yardstick computed once per local minimum). */
static void patopt(long R) {
    double U = uniform_yardstick();
    printf("uniform yardstick %.3f (Moore ball %.0f, n %lld, %d root%s)\n", U, moore_ball(), nn, nfib, nfib > 1 ? "s" : "");
    for (long t = 0; t < R; t++) {
        do { for (int i = 0; i < nslots; i++) random_slot(i, 1); build_arcs(); } while (!ypattern_ok() || !simple_ok());
        int cur = sphere2_coinc(NULL), bestc_all = 1 << 30, besty_all[32];
        /* iterated local search: hill-climb (best single-level change per slot), then kick one slot; 40 kicks */
        for (int kick = 0; kick <= 40 && bestc_all > 0; kick++) {
            if (kick > 0) {
                int saved[32]; memcpy(saved, gy, sizeof gy);
                int i = (int)(rnd() % nslots); gy[i] = slot_is_loop(i) ? 1 + (int)(rnd() % (M - 1)) : (int)(rnd() % M); build_slot_arcs(i);
                if (!ypattern_ok() || !simple_ok_slot(i)) { memcpy(gy, saved, sizeof gy); build_arcs(); continue; }
                cur = sphere2_coinc(NULL);
            }
            int improved = 1;
            while (improved && cur > 0) {
                improved = 0;
                for (int i = 0; i < nslots && cur > 0; i++) {
                    int oy = gy[i], besty = oy, bestc = cur;
                    for (int y = slot_is_loop(i) ? 1 : 0; y < M; y++) {
                        if (y == oy) continue;
                        gy[i] = y; build_slot_arcs(i);
                        if (!ypattern_ok() || !simple_ok_slot(i)) continue;
                        int c = sphere2_coinc(NULL);
                        if (c < bestc) { bestc = c; besty = y; }
                    }
                    gy[i] = besty; build_slot_arcs(i);
                    if (bestc < cur) { cur = bestc; improved = 1; }
                }
            }
            if (cur < bestc_all) { bestc_all = cur; memcpy(besty_all, gy, sizeof gy); }
        }
        memcpy(gy, besty_all, sizeof gy); build_arcs();
        if (!ypattern_ok() || !simple_ok()) continue;
        double P = yardstick();
        printf("PAT %d %.3f %.4f :", bestc_all, P, P / U); for (int i = 0; i < nslots; i++) printf(" %d", gy[i]); printf("\n");
        fflush(stdout);
    }
}

/* ---------- exact x-only single-slot sweep around the current state (bound = start defect) ---------- */
/* Order: when the sweep is capped (seconds < 1e9) the slots are visited round-robin and x in a pseudo-random cyclic
 * order (x = x0 + step * stride mod N, stride coprime to N) so that a capped sweep samples every slot evenly instead of
 * the first x-values of slot 0 only; an uncapped sweep is exhaustive in both orders. */
/* Nested completion.  The first IMPROVED state (defect < start) stops the sweep and is handed back
 * in imp_x / imp_y / imp_val; sa() adopts it, continues the chain from it and sweeps it in turn (descent through the
 * x-only neighbourhoods; with an involution in S the defect can be odd, so defect-1 states occur and must be swept). */
static int imp_x[32], imp_y[32], imp_tx = 0; static long imp_val = -1; static int nest_stop = 1;
static long xsweep(ll start_defect, double seconds, long *ncands, int *found) {
    double t0 = now();
    int sx[32], sy[32], stx = tx; memcpy(sx, gx, sizeof gx); memcpy(sy, gy, sizeof gy);
    ll best = start_defect; long cands = 0, ties = 0, improved = 0; *found = 0; imp_val = -1;
    int save_collect = collect_unc; collect_unc = 0;
    int capped = seconds < 1e9;
    int stride = 1, x0[32]; for (int i = 0; i < nslots; i++) x0[i] = capped ? (int)(rnd() % N) : 0;
    if (capped) { do stride = 1 + (int)(rnd() % (N - 1)); while (igcd(stride, N) != 1); }
    /* with --tmove the involution is swept first (all ninv values of tx; cheap), as "slot t" */
    if (invol && tmove_pct > 0 && ninv > 1) {
        for (int k = 0; k < ninv && !*found; k++) {
            int x = k * tstep;
            if (x == stx) continue;
            if ((cands & 255) == 0 && now() - t0 > seconds) { printf("  xsweep: time limit in the t-sweep after %ld candidates\n", cands); goto done; }
            tx = x; build_t_arc();
            if (!simple_ok()) { tx = stx; build_t_arc(); continue; }
            long v = eval(best, 0); cands++;
            if (v <= best) {
                if (v < best) { improved++; printf("IMPROVED %ld (slot t -> %d %d, %.0f s) :", v, x, M / 2, now() - t0); }
                else { ties++; printf("TIE %ld (slot t -> %d %d, %.0f s) :", v, x, M / 2, now() - t0); }
                print_gens(); fflush(stdout);
                if (v == 0) {
                    printf("FOUND: re-evaluation with full profiles\n");
                    eval((ll)1 << 50, 1);
                    printf("FOUND N=%d M=%d r=%lld D=%d k=%d c=%d degree %d n=%lld%s :", N, M, R, D, K, C, deg, nn, invol ? " + involution (last pair)" : ""); print_gens();
                    fflush(stdout); *found = 1; break;
                }
                if (v < best && nest_stop) {
                    memcpy(imp_x, gx, sizeof gx); memcpy(imp_y, gy, sizeof gy); imp_tx = tx; imp_val = v;
                    tx = stx; build_t_arc();
                    printf("  xsweep: improved state handed to the chain after %ld candidates\n", cands);
                    goto done;
                }
            }
            tx = stx; build_t_arc();
        }
    }
    for (long step = 0; step < N && !*found; step++) {
        for (int i = 0; i < nslots && !*found; i++) {
            int x = (int)((x0[i] + step * (ll)stride) % N);
            if (x == sx[i] || !x_ok(x)) continue;
            if ((cands & 255) == 0 && now() - t0 > seconds) { printf("  xsweep: time limit at step %ld slot %d after %ld candidates\n", step, i, cands); goto done; }
            gx[i] = x; build_slot_arcs(i);
            if (!simple_ok_slot(i)) { gx[i] = sx[i]; build_slot_arcs(i); continue; }
            long v = eval(best, 0); cands++;
            if (v <= best) {
                if (v < best) { improved++; printf("IMPROVED %ld (slot %d -> %d %d, %.0f s) :", v, i, x, sy[i], now() - t0); }
                else { ties++; printf("TIE %ld (slot %d -> %d %d, %.0f s) :", v, i, x, sy[i], now() - t0); }
                print_gens(); fflush(stdout);
                if (v == 0) {
                    printf("FOUND: re-evaluation with full profiles\n");
                    eval((ll)1 << 50, 1);
                    printf("FOUND N=%d M=%d r=%lld D=%d k=%d c=%d degree %d n=%lld%s :", N, M, R, D, K, C, deg, nn, invol ? " + involution (last pair)" : ""); print_gens();
                    fflush(stdout); *found = 1; break;
                }
                if (v < best && nest_stop) {
                    memcpy(imp_x, gx, sizeof gx); memcpy(imp_y, gy, sizeof gy); imp_tx = tx; imp_val = v;
                    gx[i] = sx[i]; build_slot_arcs(i);
                    printf("  xsweep: improved state handed to the chain after %ld candidates\n", cands);
                    goto done;
                }
            }
            gx[i] = sx[i]; build_slot_arcs(i);   /* restore the slot (the completed state is kept when found) */
        }
    }
done:
    if (!*found) { memcpy(gx, sx, sizeof gx); memcpy(gy, sy, sizeof gy); tx = stx; build_arcs(); }
    collect_unc = save_collect;
    *ncands = cands;
    printf("  xsweep done: %ld candidates, %ld ties, %ld improved, %.0f s\n", cands, ties, improved, now() - t0);
    return best;
}

/* ---------- simulated annealing ---------- */
static int near_n = 0; static int near_list[4000][33];   /* slot pairs + tx */
static int near_seen(void) {
    for (int t = 0; t < near_n; t++) { int same = near_list[t][32] == tx; for (int i = 0; i < nslots && same; i++) if (near_list[t][2 * i] != gx[i] || near_list[t][2 * i + 1] != gy[i]) same = 0; if (same) return 1; }
    if (near_n < 4000) { for (int i = 0; i < nslots; i++) { near_list[near_n][2 * i] = gx[i]; near_list[near_n][2 * i + 1] = gy[i]; } near_list[near_n][32] = tx; near_n++; }
    return 0;
}
static int fixy = 0; static unsigned ptarget_pct = 0; static long completeT = -1; static double compl_cap = 600;
static int maxcoinc = 1 << 30;   /* --maxcoinc c: y-changing moves must keep the radius-2 coincidence count <= c */
static double walkT = -1;        /* --walkT T2: temperature after the first accepted state with defect <= near_T (plateau walk) */

static int sa(double seconds, double temp, unsigned xonly_pct, long near_T) {
    double t0 = now();
    collect_unc = ptarget_pct > 0 || tmove_pct > 0;   /* t-moves aim at unreached involutions too */
    long best = eval((ll)1 << 50, 1);
    memcpy(unc, cand_unc, sizeof unc); memcpy(nunc, cand_nunc, sizeof nunc);
    printf("start: defect %ld\nBEST %ld :", best, best); print_gens(); fflush(stdout);
    long gbest = best; int gbx[32], gby[32], gbtx = tx; for (int i = 0; i < nslots; i++) { gbx[i] = gx[i]; gby[i] = gy[i]; }
    long evals = 1, accepted = 0, ntgt = 0, acc_tgt = 0, ncompl = 0, compl_cands = 0; double next_report = 60, Pcur = yardstick();
    long hist[8] = {0}, ntm = 0, acc_tm = 0, ntm_tgt = 0;
    if (best <= near_T) near_seen();
    while (gbest > 0 && now() - t0 < seconds) {
        int ox[32], oy[32], otx = tx; for (int i = 0; i < nslots; i++) { ox[i] = gx[i]; oy[i] = gy[i]; }
        int targeted = 0, tmoved = 0;
        /* involution move (targeted at an unreached involution of the current state when there is one) */
        if (invol && tmove_pct > 0 && (unsigned)(rnd() % 100) < tmove_pct) { tmoved = invol_move(); if (tmoved) { ntm++; ntm_tgt += tmoved == 1; } }
        if (!tmoved && ptarget_pct > 0 && (unsigned)(rnd() % 100) < ptarget_pct) targeted = targeted_move(fixy);
        if (!targeted && !tmoved) {
            unsigned u = fixy ? 0 : (unsigned)(rnd() % 100);
            if (u < xonly_pct) random_slot((int)(rnd() % nslots), 0);
            else if (u < xonly_pct + (100 - xonly_pct) * 3 / 5) random_slot((int)(rnd() % nslots), 1);
            else { int a = (int)(rnd() % nslots), b = (int)(rnd() % nslots); random_slot(a, (int)(rnd() % 2)); random_slot(b, (int)(rnd() % 2)); }
        }
        build_arcs();
        int ok = ypattern_ok() && simple_ok();
        int ychanged = 0; for (int i = 0; i < nslots; i++) if (gy[i] != oy[i]) ychanged = 1;
        double Pnew = Pcur;
        if (ok && ychanged) { if (sphere2_coinc(NULL) > maxcoinc) ok = 0; else { Pnew = yardstick(); if (Pnew > 2.0 * Pcur + 5) ok = 0; } }
        if (!ok) { for (int i = 0; i < nslots; i++) { gx[i] = ox[i]; gy[i] = oy[i]; } tx = otx; build_arcs(); continue; }
        ll slack = 0;
        if (temp > 0) { double uu = (rnd() % 1000000 + 0.5) / 1000000.0; slack = (ll)(-temp * log(uu)); }
        long val = eval(best + slack, 0); evals++; ntgt += targeted;
        if (val <= best + slack) { int bkt = val == 0 ? 0 : val <= 2 ? 1 : val <= 4 ? 2 : val <= 8 ? 3 : val <= 16 ? 4 : val <= 32 ? 5 : val <= 64 ? 6 : 7; hist[bkt]++; }
        if ((evals & 127) == 0) {
            double el = now() - t0;
            if (el >= next_report) { printf("  status %.0f s: %ld evaluations (%.0f/s), %ld accepted (%ld of %ld targeted, %ld of %ld t-moves (%ld aimed)), current defect %ld, global best %ld, %d near states\n", el, evals, evals / el, accepted, acc_tgt, ntgt, acc_tm, ntm, ntm_tgt, best, gbest, near_n); fflush(stdout); next_report = el + 60; }
        }
        if (val <= best + slack) {
            long prev = best;
            accepted++; acc_tgt += targeted; acc_tm += tmoved != 0; best = val; Pcur = Pnew;
            memcpy(unc, cand_unc, sizeof unc); memcpy(nunc, cand_nunc, sizeof nunc);
            int is_new_near = 0;
            const char *how = targeted ? ", targeted" : tmoved == 1 ? ", t aimed" : tmoved == 2 ? ", t random" : "";
            if (val <= near_T && !near_seen()) { is_new_near = 1; printf("NEAR %ld (%.0f s, eval %ld, from %ld%s) :", val, now() - t0, evals, prev, how); print_gens(); fflush(stdout); }
            if (val <= near_T && walkT >= 0 && temp != walkT) { printf("  plateau walk: T %.2f -> %.2f at eval %ld (%.0f s)\n", temp, walkT, evals, now() - t0); temp = walkT; fflush(stdout); }
            if (val < gbest) {
                gbest = val; gbtx = tx; for (int i = 0; i < nslots; i++) { gbx[i] = gx[i]; gby[i] = gy[i]; }
                printf("eval %ld (%.0f s): global best defect %ld%s\nBEST %ld :", evals, now() - t0, val, how[0] ? how + 1 : "", val); print_gens(); fflush(stdout);
            }
            if (val > 0 && val <= completeT && is_new_near && now() - t0 < seconds) {
                long cur = val; int hit = 0;
                for (int depth = 0; depth < 12 && now() - t0 < seconds; depth++) {
                    ncompl++;
                    printf("  completion %ld: x-only sweep of the defect-%ld state (cap %.0f s%s)\n", ncompl, cur, compl_cap, depth ? ", nested" : ""); fflush(stdout);
                    long nc; int found;
                    int cx[32], cy[32], ctx = tx; memcpy(cx, gx, sizeof gx); memcpy(cy, gy, sizeof gy);
                    xsweep(cur, compl_cap, &nc, &found); compl_cands += nc;
                    if (found) { hit = 1; break; }
                    if (imp_val < 0 || imp_val >= cur) { memcpy(gx, cx, sizeof gx); memcpy(gy, cy, sizeof gy); tx = ctx; build_arcs(); break; }
                    /* adopt the improved state: the chain moves there (fresh uncovered sample) and sweeps it in turn */
                    memcpy(gx, imp_x, sizeof gx); memcpy(gy, imp_y, sizeof gy); tx = imp_tx; build_arcs();
                    long chk = eval((ll)1 << 50, 0); evals++;
                    if (chk != imp_val) { printf("FAIL: nested re-evaluation %ld != %ld\n", chk, imp_val); exit(3); }
                    memcpy(unc, cand_unc, sizeof unc); memcpy(nunc, cand_nunc, sizeof nunc);
                    cur = chk; best = cur; near_seen();
                    printf("NEAR %ld (%.0f s, eval %ld, nested from sweep) :", cur, now() - t0, evals); print_gens(); fflush(stdout);
                    if (cur < gbest) { gbest = cur; gbtx = tx; for (int i = 0; i < nslots; i++) { gbx[i] = gx[i]; gby[i] = gy[i]; } printf("eval %ld (%.0f s): global best defect %ld (nested)\nBEST %ld :", evals, now() - t0, cur, cur); print_gens(); fflush(stdout); }
                }
                if (hit) { gbest = 0; gbtx = tx; for (int i = 0; i < nslots; i++) { gbx[i] = gx[i]; gby[i] = gy[i]; } break; }
                /* the chain continues from the (last adopted) near state; its uncovered sample is in unc[] */
            }
        } else {
            for (int i = 0; i < nslots; i++) { gx[i] = ox[i]; gy[i] = oy[i]; }
            tx = otx;
            build_arcs();
        }
    }
    for (int i = 0; i < nslots; i++) { gx[i] = gbx[i]; gy[i] = gby[i]; }
    tx = gbtx;
    build_arcs();
    printf("done: %ld evaluations, %ld accepted (%ld of %ld targeted, %ld of %ld t-moves (%ld aimed)), %ld completions (%ld candidates), %.0f s; global best defect %ld, %d near states\n",
           evals, accepted, acc_tgt, ntgt, acc_tm, ntm, ntm_tgt, ncompl, compl_cands, now() - t0, gbest, near_n);
    printf("accepted-state defect histogram [0]=%ld [<=2]=%ld [<=4]=%ld [<=8]=%ld [<=16]=%ld [<=32]=%ld [<=64]=%ld [>64]=%ld\n", hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7]);
    printf("BEST %ld :", gbest); print_gens();
    if (gbest == 0) {
        printf("FOUND: re-evaluation with full profiles\n");
        collect_unc = 0;
        eval((ll)1 << 50, 1);
        printf("FOUND N=%d M=%d r=%lld D=%d k=%d c=%d degree %d n=%lld%s :", N, M, R, D, K, C, deg, nn, invol ? " + involution (last pair)" : ""); print_gens();
    }
    fflush(stdout);
    return gbest == 0 ? 0 : 1;
}

/* ---------- involution / 2-slot neighbourhood sweeps of a near state, aimed at its first far vertex ---------- */
#define HIST 1500
static long s2_vals[HIST], s2_nh, s2_cands, s2_ties, s2_improved; static ll s2_start; static double s2_t0;
static int s2_cmp(const void *a, const void *b) { long x = *(const long *)a, y = *(const long *)b; return x < y ? -1 : x > y; }
/* evaluate the current state: exactly for the first HIST candidates of the mode (defect distribution), with bound =
 * start defect afterwards; prints TIE / IMPROVED / FOUND; returns 1 on a hit */
static long s2_nested = 0, s2_nested_cands = 0;   /* exact x-only sweeps of the improved states met in a 2-slot sweep */
static int s2_J = 2;                                /* --J2: walk budget |L| + |R| of mode E */
static unsigned s2_pairmask = 0;                    /* --pairs i-j,...: mode F restricted to these unordered slot pairs (bit i*8+j, i < j) */
static int s2_all = 0;                              /* --s2all: mode F does not skip the start values of the swept slots (whole-family sweeps) */
static long s2_step_lo = 0, s2_step_hi = -1;        /* --s2steps lo hi: mode F runs the outer-slot steps [lo, hi) only (chunks of one family) */
static int s2_x1mod = 0;                            /* --x1mod: mode F quotients the family by its residual automorphisms (outer x in 3 classes mod q = N/I) */
static int s2_quot = 0;                             /* --quot: mode F quotients the family by the EXACT stabiliser (computed by brute force) */
static long s2_bound = -1;                          /* --s2bound b: evaluate with bound b instead of the start defect (whole-family sweeps from an arbitrary start state) */

/* --quot.  The maps phi_{u,c}: x -> x^u, y -> x^c y (u a unit mod N, c mod N) are automorphisms of G iff
 * c S_M = 0 mod N (S_k = 1 + r + ... + r^(k-1)); phi(x^a y^b) = x^(u a + c S_b) y^b.  Mode F keeps the levels, t and the
 * non-swept slots fixed, so the group H of those phi that fix t = x^tx y^(M/2) and every non-swept slot acts on the
 * family and preserves the defect (isomorphic Cayley graphs).  H is computed by brute force over (u, c) (no assumption
 * on r - 1 being a unit, on I, or on x_0 being a unit: groups with a central prime have a smaller H than q = N/I).
 * The sweep then evaluates exactly one state per H-orbit: the outer slot's x only when it is the minimum of its
 * H-orbit, and the inner slot's x only when it is the minimum of its orbit under the stabiliser of the outer x in H
 * (the pairs of an H-orbit with the outer x fixed form one orbit of that stabiliser).  So the candidate count equals
 * the number of H-orbits of the family (Burnside check). */
static ll *Sk = NULL;                               /* Sk[b] = S_b mod N */
static int *qh_u = NULL, *qh_c = NULL; static long qh_n = 0;       /* H as (u, c) pairs */
static int *qs_u = NULL, *qs_c = NULL; static long qs_n = 0;       /* stabiliser of the current outer x in H */
static void quot_build(int outer, int inner) {
    if (!Sk) { Sk = malloc((M + 1) * sizeof(ll)); Sk[0] = 0; for (int b = 0; b < M; b++) Sk[b + 1] = (Sk[b] + rpow[b]) % N; }
    long cap = 1024; qh_u = realloc(qh_u, cap * sizeof(int)); qh_c = realloc(qh_c, cap * sizeof(int)); qh_n = 0;
    ll SM = Sk[M], St = invol ? Sk[M / 2] : 0;
    for (int c = 0; c < N; c++) {
        if ((c * SM) % N) continue;                              /* (x^c y)^M = x^(c S_M) must be 1 */
        for (int u = 1; u < N; u++) {
            if (igcd(u, N) != 1) continue;
            if (invol && (((ll)u * tx + c * St) % N) != tx % N) continue;     /* fixes t */
            int ok = 1;
            for (int k = 0; k < nslots && ok; k++) {
                if (k == outer || k == inner) continue;
                if ((((ll)u * gx[k] + c * Sk[gy[k]]) % N) != gx[k]) ok = 0;    /* fixes the non-swept slot k */
            }
            if (!ok) continue;
            if (qh_n == cap) { cap *= 2; qh_u = realloc(qh_u, cap * sizeof(int)); qh_c = realloc(qh_c, cap * sizeof(int)); }
            qh_u[qh_n] = u; qh_c[qh_n] = c; qh_n++;
        }
    }
    qs_u = realloc(qs_u, (qh_n + 1) * sizeof(int)); qs_c = realloc(qs_c, (qh_n + 1) * sizeof(int));
}
/* 1 iff x is the minimum of its orbit under the (u, c) list; fills the stabiliser list when stab != 0 */
static int quot_rep(int x, int level, const int *lu, const int *lc, long ln, int stab) {
    int ismin = 1; if (stab) qs_n = 0;
    for (long h = 0; h < ln; h++) {
        int y = (int)((((ll)lu[h] * x) + (ll)lc[h] * Sk[level]) % N);
        if (y < x) { ismin = 0; if (!stab) return 0; }
        if (stab && y == x) { qs_u[qs_n] = lu[h]; qs_c[qs_n] = lc[h]; qs_n++; }
    }
    return ismin;
}
static int s2_eval(const char *what) {
    long v = s2_nh < HIST ? eval((ll)1 << 50, 0) : eval(s2_start, 0);
    if (s2_nh < HIST) s2_vals[s2_nh++] = v;
    s2_cands++;
    if (v <= s2_start) {
        if (v < s2_start) { s2_improved++; printf("IMPROVED %ld (%s, %.0f s) :", v, what, now() - s2_t0); }
        else { s2_ties++; printf("TIE %ld (%s, %.0f s) :", v, what, now() - s2_t0); }
        print_gens(); fflush(stdout);
        if (v == 0) {
            printf("FOUND: re-evaluation with full profiles\n");
            eval((ll)1 << 50, 1);
            printf("FOUND N=%d M=%d r=%lld D=%d k=%d c=%d degree %d n=%lld + involution (last pair) :", N, M, R, D, K, C, deg, nn); print_gens();
            fflush(stdout); return 1;
        }
        if (v < s2_start) {   /* an improved state gets its complete x-only neighbourhood at once (cheap: K phi(N) candidates) */
            long nc; int found; int save_nest = nest_stop; nest_stop = 0;
            printf("  nested: exact x-only sweep of the improved defect-%ld state\n", v); fflush(stdout);
            xsweep(v, 1e18, &nc, &found);
            nest_stop = save_nest; s2_nested++; s2_nested_cands += nc;
            if (found) return 1;   /* xsweep keeps the completed state in gx/gy/tx */
        }
    }
    return 0;
}

/* J-bounded ball over the arcs outside slot j (non-backtracking walks of length 0..J from e, with repeats). */
static int *ballx[8], *bally[8], *balll[8]; static long balln[8], ballcap[8];
static void ball_level_alloc(int l, long need) {
    if (ballcap[l] >= need) return;
    ballcap[l] = need + 16; ballx[l] = realloc(ballx[l], ballcap[l] * sizeof(int)); bally[l] = realloc(bally[l], ballcap[l] * sizeof(int)); balll[l] = realloc(balll[l], ballcap[l] * sizeof(int));
}
static void build_ball(int skip_slot, int J) {
    ball_level_alloc(0, 1); ballx[0][0] = 0; bally[0][0] = 0; balll[0][0] = -1; balln[0] = 1;
    for (int l = 1; l <= J; l++) {
        ball_level_alloc(l, balln[l - 1] * (long)narcs);
        long n = 0;
        for (long k = 0; k < balln[l - 1]; k++) {
            int px = ballx[l - 1][k], py = bally[l - 1][k], last = balll[l - 1][k];
            for (int a = 0; a < narcs; a++) {
                if (arc_slot(a) == skip_slot) continue;
                if (last >= 0 && a == inv_arc[last]) continue;
                int cy = py + arcs[a].ay; if (cy >= M) cy -= M;
                ballx[l][n] = (int)((px + (rpow[py] * (ll)arcs[a].ax) % N) % N); bally[l][n] = cy; balll[l][n] = a; n++;
            }
        }
        balln[l] = n;
    }
}
/* all x-coordinates x_j such that slot j := x^x_j y^(y_j) (or its inverse) equals L v R for walks L, R over the
 * other arcs with |L| + |R| <= J, i.e. the walk L^-1 s_j^(+-1) R^-1 of length <= J + 1 ends at v.  Deduplicated
 * (stamp array), units only, the current value excluded.  Returns the count; the values are in tlist[]. */
static unsigned *tstamp = NULL, tstamp_cur = 0; static int *tlist = NULL;
static int targeted_set(int j, int J, int vx, int vy) {
    if (!tstamp) { tstamp = calloc(N, sizeof(unsigned)); tlist = malloc(N * sizeof(int)); }
    if (++tstamp_cur == 0) { memset(tstamp, 0, N * sizeof(unsigned)); tstamp_cur = 1; }
    build_ball(j, J);
    int yj = gy[j], yjinv = (M - yj) % M, cnt = 0, curx = gx[j];
    for (int l = 0; l <= J; l++) for (long a = 0; a < balln[l]; a++) {
        ll px, py; gmul2(ballx[l][a], bally[l][a], vx, vy, &px, &py);            /* P = L v */
        for (int m = 0; m + l <= J; m++) for (long b = 0; b < balln[m]; b++) {
            ll sx, sy; gmul2(px, py, ballx[m][b], bally[m][b], &sx, &sy);        /* s = L v R */
            int cx;
            if (sy == yj) cx = (int)sx; else if (sy == yjinv) cx = inv_x((int)sx, (int)sy); else continue;
            if (cx == curx || !x_ok(cx) || tstamp[cx] == tstamp_cur) continue;
            tstamp[cx] = tstamp_cur; tlist[cnt++] = cx;
        }
    }
    return cnt;
}
static void s2_summary(char mode) {
    qsort(s2_vals, s2_nh, sizeof(long), s2_cmp);
    double mean = 0, lq = 0; long nq = s2_nh / 4 > 0 ? s2_nh / 4 : 1, hist[8] = {0};
    for (long i = 0; i < s2_nh; i++) { mean += s2_vals[i]; long v = s2_vals[i]; hist[v == 0 ? 0 : v <= 2 ? 1 : v <= 4 ? 2 : v <= 8 ? 3 : v <= 16 ? 4 : v <= 32 ? 5 : v <= 64 ? 6 : 7]++; }
    for (long i = 0; i < nq; i++) lq += s2_vals[i];
    if (s2_nh) { mean /= s2_nh; lq /= nq; }
    printf("SWEEP2 %c done: %ld candidates, %ld ties, %ld improved, %.0f s; exact sample of %ld: mean %.2f, lower-quartile mean %.2f, min %ld, q25 %ld, median %ld, q75 %ld, max %ld; hist [0]=%ld [<=2]=%ld [<=4]=%ld [<=8]=%ld [<=16]=%ld [<=32]=%ld [<=64]=%ld [>64]=%ld\n",
           mode, s2_cands, s2_ties, s2_improved, now() - s2_t0, s2_nh, mean, lq, s2_nh ? s2_vals[0] : -1, s2_nh ? s2_vals[s2_nh / 4] : -1, s2_nh ? s2_vals[s2_nh / 2] : -1, s2_nh ? s2_vals[3 * s2_nh / 4] : -1, s2_nh ? s2_vals[s2_nh - 1] : -1,
           hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7]);
    fflush(stdout);
}
static int sweep2(double secs, const char *modes) {
    collect_unc = 1;
    long start = eval((ll)1 << 50, 1);
    collect_unc = 0;
    memcpy(unc, cand_unc, sizeof unc); memcpy(nunc, cand_nunc, sizeof nunc);
    printf("start: defect %ld\nBEST %ld :", start, start); print_gens(); fflush(stdout);
    if (start == 0) { printf("FOUND N=%d M=%d r=%lld D=%d k=%d c=%d degree %d n=%lld + involution (last pair) :", N, M, R, D, K, C, deg, nn); print_gens(); return 0; }
    int have_far = start <= UMAX;   /* modes A and F do not aim at a far vertex, so any start defect is fine for them */
    if (!have_far) { printf("sweep2: start defect %ld > %d, no complete far-vertex list: only modes A and F can run\n", start, UMAX); nunc[0] = 0; }
    ll NM = (ll)N * M, vcode = -1; int vx = 0, vy = 0;
    printf("  far vertices (x y%s):", invol ? " invol" : "");
    for (int i = 0; i < nunc[0]; i++) {
        ll idx = unc[0][i] % NM; int ux = (int)(idx % N), uy = (int)(idx / N), isinv = uy == M / 2 && tx_ok(ux);
        printf("  %d %d %d", ux, uy, isinv);
        if (vcode < 0 || (isinv && !(vy == M / 2 && tx_ok(vx)))) { vcode = unc[0][i]; vx = ux; vy = uy; }
    }
    printf("\n  aim: v = x^%d y^%d (%s); modes %s, %.0f s each, exact sample %d per mode\n", vx, vy, vy == M / 2 && tx_ok(vx) ? "involution" : "not an involution", modes, secs, HIST);
    fflush(stdout);
    int sx[32], sy[32], stx = tx; memcpy(sx, gx, sizeof gx); memcpy(sy, gy, sizeof gy);
    int found = 0; char what[128];
    if (s2_bound >= 0 && s2_bound < start) { printf("  bound %ld (--s2bound) instead of the start defect %ld\n", s2_bound, start); start = s2_bound; }
    for (const char *mp = modes; *mp && !found; mp++) {
        char mode = *mp;
        if (!have_far && mode != 'A' && mode != 'F') { printf("SWEEP2 %c skipped: no far-vertex list\n", mode); continue; }
        s2_nh = s2_cands = s2_ties = s2_improved = 0; s2_start = start; s2_t0 = now();
        int capped = 0;
        if (mode == 'A') {                      /* every involution, S fixed */
            for (int k = 0; k < ninv && !found; k++) {
                int x = k * tstep; if (x == stx) continue;
                if ((s2_cands & 31) == 0 && now() - s2_t0 > secs) { capped = 1; break; }
                tx = x; build_t_arc();
                if (simple_ok()) { snprintf(what, sizeof what, "slot t -> %d %d", x, M / 2); found = s2_eval(what); }
                if (!found) { tx = stx; build_t_arc(); }
            }
        } else if (mode == 'B') {               /* t := v, then every unit x of every slot */
            if (!(vy == M / 2 && tx_ok(vx)) || vx == stx) { printf("SWEEP2 B skipped: v is not an involution other than t\n"); continue; }
            tx = vx; build_t_arc();
            if (!simple_ok()) { printf("SWEEP2 B skipped: t := v is not simple\n"); tx = stx; build_t_arc(); continue; }
            int stride, x0[32]; for (int i = 0; i < nslots; i++) x0[i] = (int)(rnd() % N);
            do stride = 1 + (int)(rnd() % (N - 1)); while (igcd(stride, N) != 1);
            for (long step = 0; step < N && !found && !capped; step++) {
                for (int i = 0; i < nslots && !found; i++) {
                    int x = (int)((x0[i] + step * (ll)stride) % N);
                    if (x == sx[i] || !x_ok(x)) continue;
                    if ((s2_cands & 31) == 0 && now() - s2_t0 > secs) { capped = 1; break; }
                    gx[i] = x; build_slot_arcs(i);
                    if (simple_ok_slot(i)) { snprintf(what, sizeof what, "t := v %d %d, slot %d -> %d %d", vx, M / 2, i, x, sy[i]); found = s2_eval(what); }
                    if (!found) { gx[i] = sx[i]; build_slot_arcs(i); }
                }
            }
            if (!found) { tx = stx; build_t_arc(); }
        } else if (mode == 'C') {               /* every involution, then one targeted x per slot aimed at v */
            for (int k = 0; k < ninv && !found && !capped; k++) {
                int x = k * tstep;
                tx = x; build_t_arc();
                if (!simple_ok()) { tx = stx; build_t_arc(); continue; }
                for (int j = 0; j < nslots && !found; j++) {
                    if ((s2_cands & 31) == 0 && now() - s2_t0 > secs) { capped = 1; break; }
                    if (!targeted_move_at(1, j, vcode)) continue;
                    build_slot_arcs(j);
                    if (simple_ok_slot(j)) { snprintf(what, sizeof what, "slot t -> %d %d, slot %d -> %d %d aimed", x, M / 2, j, gx[j], gy[j]); found = s2_eval(what); }
                    if (!found) { gx[j] = sx[j]; gy[j] = sy[j]; build_slot_arcs(j); }
                }
                if (!found) { tx = stx; build_t_arc(); }
            }
        } else if (mode == 'D') {               /* every ordered slot pair (i, j): every unit x_i, then one targeted x_j aimed at v */
            int stride, x0[32]; for (int i = 0; i < nslots; i++) x0[i] = (int)(rnd() % N);
            do stride = 1 + (int)(rnd() % (N - 1)); while (igcd(stride, N) != 1);
            for (long step = 0; step < N && !found && !capped; step++) {
                for (int i = 0; i < nslots && !found && !capped; i++) {
                    int x = (int)((x0[i] + step * (ll)stride) % N);
                    if (x == sx[i] || !x_ok(x)) continue;
                    gx[i] = x; build_slot_arcs(i);
                    if (!simple_ok_slot(i)) { gx[i] = sx[i]; build_slot_arcs(i); continue; }
                    for (int j = 0; j < nslots && !found; j++) {
                        if (j == i) continue;
                        if ((s2_cands & 31) == 0 && now() - s2_t0 > secs) { capped = 1; break; }
                        if (!targeted_move_at(1, j, vcode)) continue;
                        build_slot_arcs(j);
                        if (simple_ok_slot(j)) { snprintf(what, sizeof what, "slot %d -> %d %d, slot %d -> %d %d aimed", i, x, sy[i], j, gx[j], gy[j]); found = s2_eval(what); }
                        if (!found) { gx[j] = sx[j]; gy[j] = sy[j]; build_slot_arcs(j); }
                    }
                    if (!found) { gx[i] = sx[i]; build_slot_arcs(i); }
                }
            }
        } else if (mode == 'E') {               /* every ordered slot pair (i, j): every unit x_i, then EVERY x_j whose
                                                 * slot is L v R with |L| + |R| <= J2 (exact targeted set, deduplicated); the
                                                 * aim cycles through the far vertices with the x_i index */
            int stride, x0[32]; for (int i = 0; i < nslots; i++) x0[i] = (int)(rnd() % N);
            do stride = 1 + (int)(rnd() % (N - 1)); while (igcd(stride, N) != 1);
            long nxi = 0, tsum = 0, tmax = 0, tmin = 1L << 40; double tlast = s2_t0;
            for (long step = 0; step < N && !found && !capped; step++) {
                ll acode = unc[0][step % nunc[0]]; ll aidx = acode % NM; int ax = (int)(aidx % N), ay = (int)(aidx / N);
                for (int i = 0; i < nslots && !found && !capped; i++) {
                    int x = (int)((x0[i] + step * (ll)stride) % N);
                    if (x == sx[i] || !x_ok(x)) continue;
                    if (now() - s2_t0 > secs) { capped = 1; break; }
                    gx[i] = x; build_slot_arcs(i);
                    if (!simple_ok_slot(i)) { gx[i] = sx[i]; build_slot_arcs(i); continue; }
                    nxi++;
                    for (int j = 0; j < nslots && !found; j++) {
                        if (j == i) continue;
                        int nt = targeted_set(j, s2_J, ax, ay);
                        tsum += nt; if (nt > tmax) tmax = nt; if (nt < tmin) tmin = nt;
                        for (int q = 0; q < nt && !found; q++) {
                            gx[j] = tlist[q]; build_slot_arcs(j);
                            if (simple_ok_slot(j)) { snprintf(what, sizeof what, "slot %d -> %d %d, slot %d -> %d %d aimed at %d %d", i, x, sy[i], j, gx[j], sy[j], ax, ay); found = s2_eval(what); }
                            if (!found) { gx[j] = sx[j]; build_slot_arcs(j); }
                        }
                    }
                    if (!found) { gx[i] = sx[i]; build_slot_arcs(i); }
                    if (now() - tlast > 60) { tlast = now(); printf("  E status %.0f s: %ld x_i values (%ld (i, j) pairs), %ld candidates (%.0f/s), targeted set per (i, x_i, j): mean %.1f min %ld max %ld; prerow %ld of %ld rejected; %ld ties %ld improved\n", tlast - s2_t0, nxi, nxi * (nslots - 1), s2_cands, s2_cands / (tlast - s2_t0), tsum / (double)(nxi * (nslots - 1)), tmin, tmax, prerow_rejects, prerow_calls, s2_ties, s2_improved); fflush(stdout); }
                }
            }
            printf("  E: %ld x_i values, %ld (i, x_i, j) triples, targeted set mean %.1f min %ld max %ld (J %d); prerow %ld of %ld rejected; nested x-sweeps %ld (%ld candidates)\n", nxi, nxi * (nslots - 1), nxi ? tsum / (double)(nxi * (nslots - 1)) : 0.0, tmin, tmax, s2_J, prerow_rejects, prerow_calls, s2_nested, s2_nested_cands);
        } else if (mode == 'F') {               /* the complete 2-slot x-neighbourhood: every unordered pair (i < j), every
                                                 * unit x_i (random order), every unit x_j */
            int stride, x0[32]; for (int i = 0; i < nslots; i++) x0[i] = (int)(rnd() % N);
            do stride = 1 + (int)(rnd() % (N - 1)); while (igcd(stride, N) != 1);
            long nxi = 0; double tlast = s2_t0;
            /* --s2steps lo hi restricts the outer loop to the steps [lo, hi) (x_i = x0 + step * stride visits every
             * residue once as step runs over 0..N-1, so the chunks of a partition of [0, N) partition the family exactly;
             * x0 and stride depend only on --seed, hence are identical in every chunk) */
            long step_lo = s2_step_lo, step_hi = s2_step_hi < 0 || s2_step_hi > N ? N : s2_step_hi;
            if (step_lo > 0 || step_hi < N) printf("  F: steps [%ld, %ld) of %d (--s2steps)\n", step_lo, step_hi, N);
            /* --x1mod: quotient by the residual automorphisms phi_{u,c}: x -> x^u, y -> x^c y with c = 0 mod ninv (they fix
             * t = y^(M/2) and the levels; u = 1 - c S_{y_0} restores x_0 = 1; S_b = 1 + r + ... + r^(b-1)).  They move the
             * outer slot's x along x + c (S_{y_i} - S_{y_0} x) mod N, i.e. through every residue mod q = N/ninv when
             * d = ninv (S_{y_i} - S_{y_0} x) is a unit mod q; so keeping only x = 1, 2 or 3 mod q (three classes, because
             * the orbit member with x = 1 may be invalid: at most two members are, the one with u = 0 mod q and the one
             * whose inner-slot image is 0 mod q) plus every x whose d is not a unit mod q meets every orbit at least once:
             * ~4/q of the states.  Off for q <= 4. */
            long qres = N / ninv, Sy0 = 0, Syi = 0;
            if (s2_x1mod && qres > 4) { for (int k = 0; k < gy[0]; k++) Sy0 = (Sy0 + rpow[k]) % N; printf("  F: --x1mod: outer slot restricted to x = 1, 2, 3 mod q = %ld plus the non-unit-step classes\n", qres); }
            /* --quot needs exactly one unordered pair (the stabiliser depends on which slots are fixed) */
            int q_outer = -1, q_inner = -1; long q_outer_reps = 0;
            if (s2_quot) {
                int npairs = 0;
                for (int i = 0; i < nslots; i++) for (int j = i + 1; j < nslots; j++) if (s2_pairmask & (1u << (i * 8 + j))) { npairs++; q_outer = i; q_inner = j; }
                if (npairs != 1) { printf("FAIL: --quot needs --pairs with exactly one pair\n"); exit(2); }
                quot_build(q_outer, q_inner);
                printf("  F: --quot: stabiliser of t and the fixed slots in {x -> x^u, y -> x^c y} has order %ld (q = N/I = %ld); one candidate per orbit\n", qh_n, qres);
            }
            for (long step = step_lo; step < step_hi && !found && !capped; step++) {
                for (int i = 0; i < nslots && !found && !capped; i++) {
                    int x = (int)((x0[i] + step * (ll)stride) % N);
                    if ((x == sx[i] && !s2_all) || !x_ok(x)) continue;
                    if (s2_x1mod && qres > 4) {
                        Syi = 0; for (int k = 0; k < gy[i]; k++) Syi = (Syi + rpow[k]) % N;
                        long d = (long)(((Syi - (Sy0 * (ll)x) % N + N) % N) % qres);
                        if (x % qres != 1 && x % qres != 2 && x % qres != 3 && igcd((int)d, (int)qres) == 1) continue;
                    }
                    if (s2_quot && i == q_outer) { if (!quot_rep(x, gy[i], qh_u, qh_c, qh_n, 1)) continue; q_outer_reps++; }
                    if (now() - s2_t0 > secs) { capped = 1; break; }
                    gx[i] = x; build_slot_arcs(i);
                    if (!simple_ok_slot(i)) { gx[i] = sx[i]; build_slot_arcs(i); continue; }
                    { int any = 0; for (int j = i + 1; j < nslots; j++) if (!s2_pairmask || (s2_pairmask & (1u << (i * 8 + j)))) any = 1; if (!any) { gx[i] = sx[i]; build_slot_arcs(i); continue; } }
                    nxi++;
                    for (int j = i + 1; j < nslots && !found; j++) {
                        if (s2_pairmask && !(s2_pairmask & (1u << (i * 8 + j)))) continue;   /* --pairs i-j,...: only these unordered pairs */
                        for (int xj = 0; xj < N && !found; xj++) {
                            if ((xj == sx[j] && !s2_all) || !x_ok(xj)) continue;
                            if (s2_quot && qs_n > 1 && !quot_rep(xj, gy[j], qs_u, qs_c, qs_n, 0)) continue;
                            gx[j] = xj; build_slot_arcs(j);
                            if (simple_ok_slot(j)) { snprintf(what, sizeof what, "slot %d -> %d %d, slot %d -> %d %d", i, x, sy[i], j, xj, sy[j]); found = s2_eval(what); }
                            if (!found) { gx[j] = sx[j]; build_slot_arcs(j); }
                        }
                    }
                    if (!found) { gx[i] = sx[i]; build_slot_arcs(i); }
                    if (now() - tlast > 60) { tlast = now(); printf("  F status %.0f s: %ld x_i values, %ld candidates (%.0f/s); prerow %ld of %ld rejected; %ld ties %ld improved\n", tlast - s2_t0, nxi, s2_cands, s2_cands / (tlast - s2_t0), prerow_rejects, prerow_calls, s2_ties, s2_improved); fflush(stdout); }
                }
            }
            printf("  F: %ld x_i values; prerow %ld of %ld rejected; nested x-sweeps %ld (%ld candidates)\n", nxi, prerow_rejects, prerow_calls, s2_nested, s2_nested_cands);
            if (s2_quot) printf("  F: --quot: %ld outer-slot orbit representatives, |H| = %ld, %ld candidates = H-orbits of the family\n", q_outer_reps, qh_n, s2_cands);
        } else { printf("SWEEP2: unknown mode %c\n", mode); continue; }
        if (capped) printf("  sweep2 %c: time limit after %ld candidates\n", mode, s2_cands);
        s2_summary(mode);
    }
    if (!found) { memcpy(gx, sx, sizeof gx); memcpy(gy, sy, sizeof gy); tx = stx; build_arcs(); }
    return found ? 0 : 1;
}

/* ---------- the full single-slot sweep ---------- */
typedef struct { int slot, y; double P; } unit_t;
static int cmp_unit(const void *a, const void *b) { double x = ((const unit_t *)a)->P, y = ((const unit_t *)b)->P; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: cayley_sa N M r D k c [--init x1 y1 ...] [options]\n"); return 2; }
    N = atoi(argv[1]); M = atoi(argv[2]); R = atoll(argv[3]); D = atoi(argv[4]); K = atoi(argv[5]); C = atoi(argv[6]);
    nfib = C ? 2 : 1; nslots = C ? 2 * K + C - 1 : K; deg = C ? 2 * K + C : 2 * K; nrows = nfib * M;
    for (int i = 7; i < argc; i++) { if (!strcmp(argv[i], "--invol")) invol = 1; if (!strcmp(argv[i], "--xunit")) xunit = 1; if (!strcmp(argv[i], "--xunitI")) xunit = 2; if (!strcmp(argv[i], "--qyard")) qyard_Q = atoi(argv[i + 1]); }
    if (invol && C) { printf("FAIL: --invol is for Cayley graphs (c = 0)\n"); return 2; }
    if (M % 2 == 0 && !invol) { printf("FAIL: M must be odd (or use --invol)\n"); return 2; }
    if (M % 2 && invol) { printf("FAIL: --invol needs M even\n"); return 2; }
    if (invol) deg = 2 * K + 1;
    rpow = malloc(M * sizeof(ll));
    if (R == 0) {
        for (ll r = 2; r < N; r++) {
            ll p = 1; int ord = 0;
            for (int j = 1; j <= M; j++) { p = (p * r) % N; if (p == 1) { ord = j; break; } }
            if (ord == M) { R = r; break; }
        }
        if (R == 0) { printf("FAIL: no element of order %d mod %d\n", M, N); return 2; }
    }
    rpow[0] = 1 % N;
    for (int j = 1; j < M; j++) { rpow[j] = (rpow[j - 1] * R) % N; if (rpow[j] == 1 % N) { printf("FAIL: r has order %d < M\n", j); return 1; } }
    if ((rpow[M - 1] * R) % N != 1 % N) { printf("FAIL: r^M != 1 mod N\n"); return 1; }
    for (int a = 0; a < 64; a++) arcs[a].T = malloc(M * sizeof(int));
    W = (N + 63) / 64; nn = (ll)nrows * N; lastmask = (N % 64) ? ((1ULL << (N % 64)) - 1) : ~0ULL;
    size_t words = (size_t)nrows * W + 2;
    vis = calloc(words, 8); cur = calloc(words, 8); nxt = calloc(words, 8);
    rowflag_cur = calloc(nrows, 1); rowflag_nxt = calloc(nrows, 1); rowstamp = calloc(nrows, sizeof(unsigned));
    scalar_thresh = (long)nrows * W / 2;
    targetJ = D - 1;
    if (invol) { ll h = rpow[M / 2]; ninv = igcd(N, (int)((1 + h) % N)); tstep = N / ninv; }   /* involution x-coordinates = multiples of tstep */

    int have_init = 0, have_ypat = 0, do_eval = 0, do_selftest = 0, order_yard = 0, from_u = 0, to_u = 1 << 30, yfrom = 0, xonly = 0;
    double seconds = 1e18; int slot_list[32], nsl = 0; long patgen_K = 0, patopt_R = 0;
    double sa_seconds = 0, sa_T = 0; unsigned sa_xonly = 100; long near_T = 4; int freey = 0; long sample_S = 0;
    double sweep2_secs = 0; const char *sweep2_modes = "ABCD"; long prerow_test = 0;
    for (int i = 7; i < argc; i++) if (!strcmp(argv[i], "--freey")) freey = 1;
    for (int i = 7; i < argc; i++) {
        if (!strcmp(argv[i], "--sa")) { sa_seconds = atof(argv[i + 1]); sa_T = atof(argv[i + 2]); sa_xonly = (unsigned)atoi(argv[i + 3]); i += 3; continue; }
        if (!strcmp(argv[i], "--tx")) { tx = atoi(argv[++i]); if (!invol || !tx_ok(tx)) { printf("FAIL: --tx %d is not an involution x-coordinate (needs --invol, 0 mod %d)\n", tx, tstep); return 2; } continue; }
        if (!strcmp(argv[i], "--tmove")) { tmove_pct = (unsigned)atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], "--sweep2")) { sweep2_secs = atof(argv[i + 1]); sweep2_modes = argv[i + 2]; i += 2; continue; }
        if (!strcmp(argv[i], "--prerowtest")) { prerow_test = atol(argv[++i]); continue; }
        if (!strcmp(argv[i], "--J2")) { s2_J = atoi(argv[++i]); if (s2_J < 0 || s2_J > 6) { printf("FAIL: --J2 in 0..6\n"); return 2; } continue; }
        if (!strcmp(argv[i], "--pairs")) { char *p = argv[++i]; while (*p) { int a = (int)strtol(p, &p, 10); if (*p == '-') p++; int b = (int)strtol(p, &p, 10); if (a > b) { int t = a; a = b; b = t; } if (a < 0 || b > 7 || a == b) { printf("FAIL: --pairs needs i-j with 0 <= i < j <= 7\n"); return 2; } s2_pairmask |= 1u << (a * 8 + b); if (*p == ',') p++; } continue; }
        if (!strcmp(argv[i], "--s2bound")) { s2_bound = atol(argv[++i]); continue; }
        if (!strcmp(argv[i], "--s2all")) { s2_all = 1; continue; }
        if (!strcmp(argv[i], "--s2steps")) { s2_step_lo = atol(argv[i + 1]); s2_step_hi = atol(argv[i + 2]); i += 2; if (s2_step_lo < 0 || s2_step_hi <= s2_step_lo) { printf("FAIL: --s2steps lo hi needs 0 <= lo < hi\n"); return 2; } continue; }
        if (!strcmp(argv[i], "--x1mod")) { s2_x1mod = 1; continue; }
        if (!strcmp(argv[i], "--quot")) { s2_quot = 1; continue; }
        if (!strcmp(argv[i], "--noprerow")) { prerow_on = 0; continue; }
        if (!strcmp(argv[i], "--prerow_max")) { prerow_max = atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], "--seed")) { rs = (u64)atoll(argv[++i]) * 2654435761ULL + 88172645463325252ULL; continue; }
        if (!strcmp(argv[i], "--near")) { near_T = atol(argv[++i]); continue; }
        if (!strcmp(argv[i], "--ptarget")) { ptarget_pct = (unsigned)atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], "--J")) { targetJ = atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], "--fixy")) { fixy = 1; continue; }
        if (!strcmp(argv[i], "--maxcoinc")) { maxcoinc = atoi(argv[++i]); continue; }
        if (!strcmp(argv[i], "--walkT")) { walkT = atof(argv[++i]); continue; }
        if (!strcmp(argv[i], "--invol") || !strcmp(argv[i], "--xunit") || !strcmp(argv[i], "--xunitI")) continue;
        if (!strcmp(argv[i], "--qyard")) { ++i; continue; }
        if (!strcmp(argv[i], "--sample")) { sample_S = atol(argv[++i]); continue; }
        if (!strcmp(argv[i], "--freey")) { fixy = 0; continue; }
        if (!strcmp(argv[i], "--complete")) { completeT = atol(argv[++i]); continue; }
        if (!strcmp(argv[i], "--compl_cap")) { compl_cap = atof(argv[++i]); continue; }
        if (!strcmp(argv[i], "--patgen")) { patgen_K = atol(argv[++i]); continue; }
        if (!strcmp(argv[i], "--patopt")) { patopt_R = atol(argv[++i]); continue; }
        if (!strcmp(argv[i], "--ypat")) {
            if (argc < i + 1 + nslots) { printf("FAIL: --ypat needs %d values\n", nslots); return 2; }
            for (int s = 0; s < nslots; s++) gy[s] = ((atoi(argv[i + 1 + s]) % M) + M) % M;
            i += nslots; have_ypat = 1; if (!freey) fixy = 1; continue;
        }
        if (!strcmp(argv[i], "--init")) {
            if (argc < i + 1 + 2 * nslots) { printf("FAIL: --init needs %d pairs\n", nslots); return 2; }
            for (int s = 0; s < nslots; s++) { gx[s] = atoi(argv[i + 1 + 2 * s]) % N; gy[s] = atoi(argv[i + 2 + 2 * s]) % M; }
            i += 2 * nslots; have_init = 1;
        } else if (!strcmp(argv[i], "--slots")) { char *p = argv[++i]; while (*p) { slot_list[nsl++] = (int)strtol(p, &p, 10); if (*p == ',') p++; } }
        else if (!strcmp(argv[i], "--yfrom")) yfrom = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--seconds")) seconds = atof(argv[++i]);
        else if (!strcmp(argv[i], "--thresh")) scalar_thresh = atol(argv[++i]);
        else if (!strcmp(argv[i], "--order")) { if (!strcmp(argv[++i], "yardstick")) order_yard = 1; }
        else if (!strcmp(argv[i], "--from")) from_u = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--to")) to_u = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--eval")) do_eval = 1;
        else if (!strcmp(argv[i], "--xonly")) xonly = 1;
        else if (!strcmp(argv[i], "--selftest")) do_selftest = 1;
        else { printf("FAIL: unknown option %s\n", argv[i]); return 2; }
    }
    printf("G = Z_%d x|_%lld Z_%d (|G| = %lld), %s, degree %d, D = %d, n = %lld, %d slots, W = %d, scalar threshold %ld\n",
           N, R, M, (ll)N * M, C ? "2-vertex lift" : "Cayley graph", deg, D, nn, nslots, W, scalar_thresh);
    printf("Moore ball %.0f, Moore/n = %.3f, uniform yardstick %.2f (%d root%s)\n", moore_ball(), moore_ball() / nn, uniform_yardstick(), nfib, nfib > 1 ? "s" : "");
    if (invol) {
        ll h = rpow[M / 2];
        int cidx = igcd(N, (int)((h - 1 + N) % N));
        printf("involution t = x^%d y^%d (r^(M/2) = %lld mod N): %d involutions x^i y^(M/2) (i = 0 mod %d), centraliser of t has index %d%s; x units only: %d; t-moves %u %%\n",
               tx, M / 2, h, ninv, tstep, N / cidx, h == 1 % N ? " (CENTRAL: bad)" : "", xunit, tmove_pct);
    }
    if (do_selftest) { build_arcs(); return selftest(12) ? 1 : 0; }
    if (prerow_test > 0) {   /* random states in the pattern (or --init), bounds 0..4 with the pre-filter on: same verdicts as the exact evaluation */
        if (!have_ypat && !have_init) { printf("FAIL: --prerowtest needs --ypat or --init\n"); return 2; }
        long bad = 0, tested = 0, rej = 0; double t_full = 0, t_pre = 0, t_nopre = 0;
        for (long t = 0; t < prerow_test; t++) {
            if (!have_init || t > 0) { for (int i = 0; i < nslots; i++) random_slot(i, 0); }
            build_arcs(); if (!ypattern_ok() || !simple_ok()) continue;
            double a = now(); long full = eval((ll)1 << 50, 0); t_full += now() - a;
            for (ll b = 0; b <= 4; b++) {
                prerow_on = 1; double c = now(); long r = eval(b, 0); t_pre += now() - c;
                prerow_on = 0; c = now(); long r2 = eval(b, 0); t_nopre += now() - c;
                prerow_on = 1; tested++;
                if ((r <= b) != (full <= b) || (r <= b && r != full) || (r2 <= b) != (full <= b)) { bad++; printf("MISMATCH state %ld bound %lld: full %ld, prerow %ld, noprerow %ld\n", t, b, full, r, r2); }
            }
        }
        printf("prerowtest: %ld verdicts, %ld mismatches; prerow rejected %ld of %ld calls; mean time per eval: full %.1f us, bounded with prerow %.1f us, bounded without %.1f us\n", tested, bad, prerow_rejects, prerow_calls, 1e6 * t_full / (tested / 5.0), 1e6 * t_pre / tested, 1e6 * t_nopre / tested);
        (void)rej;
        return bad ? 1 : 0;
    }
    if (sweep2_secs > 0) {
        if (!have_init || !invol || C) { printf("FAIL: --sweep2 needs --init and --invol (Cayley graphs)\n"); return 2; }
        build_arcs();
        if (!ypattern_ok() || !simple_ok()) { printf("FAIL: start state is not a simple connected candidate\n"); return 1; }
        return sweep2(sweep2_secs, sweep2_modes);
    }
    if (patgen_K > 0) { patgen(patgen_K); return 0; }
    if (patopt_R > 0) { patopt(patopt_R); return 0; }
    if (sample_S > 0) {   /* --sample S --ypat ...: S random x-assignments in the fixed pattern, defect quantiles */
        if (!have_ypat) { printf("FAIL: --sample needs --ypat\n"); return 2; }
        long *vals = malloc(sample_S * sizeof(long)); long got = 0; double t0 = now(); long hist[8] = {0};
        while (got < sample_S) {
            for (int i = 0; i < nslots; i++) random_slot(i, 0);
            build_arcs(); if (!ypattern_ok() || !simple_ok()) continue;
            long v = eval((ll)1 << 50, 0); vals[got++] = v;
            int bkt = v == 0 ? 0 : v <= 2 ? 1 : v <= 4 ? 2 : v <= 8 ? 3 : v <= 16 ? 4 : v <= 32 ? 5 : v <= 64 ? 6 : 7; hist[bkt]++;
            if (v <= near_T) { printf("NEAR %ld :", v); print_gens(); fflush(stdout); }
        }
        /* sort */
        for (long i = 1; i < got; i++) { long v = vals[i], j = i - 1; while (j >= 0 && vals[j] > v) { vals[j + 1] = vals[j]; j--; } vals[j + 1] = v; }
        double mean = 0, lq = 0; for (long i = 0; i < got; i++) mean += vals[i]; mean /= got;
        long nq = got / 4 > 0 ? got / 4 : 1; for (long i = 0; i < nq; i++) lq += vals[i]; lq /= nq;
        printf("SAMPLE %ld random states in %.1f s: mean %.2f, lower-quartile mean %.2f, min %ld, q25 %ld, median %ld, q75 %ld, max %ld; hist [0]=%ld [<=2]=%ld [<=4]=%ld [<=8]=%ld [<=16]=%ld [<=32]=%ld [<=64]=%ld [>64]=%ld; uniform P %.2f, pattern P %.2f\n",
               got, now() - t0, mean, lq, vals[0], vals[got / 4], vals[got / 2], vals[3 * got / 4], vals[got - 1], hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7], uniform_yardstick(), yardstick());
        return 0;
    }
    if (sa_seconds > 0) {
        if (have_ypat && !have_init) { for (int i = 0; i < nslots; i++) random_slot(i, 0); build_arcs(); int tries = 0; while ((!ypattern_ok() || !simple_ok()) && tries++ < 1000) { for (int i = 0; i < nslots; i++) random_slot(i, 0); build_arcs(); } }
        else if (!have_init) pick_ypattern(400);
        build_arcs();
        if (!ypattern_ok() || !simple_ok()) { printf("FAIL: start state is not a simple connected candidate\n"); return 1; }
        int s2; int co = sphere2_coinc(&s2);
        printf("SA %.0f s, T = %.2f, x-only moves %u %%, fixed y %d, targeted %u %% (J %d), near threshold %ld, complete <= %ld (cap %.0f s); start pattern: yardstick P = %.2f (ratio %.3f), radius-2 coincidences %d (sphere %d), y =",
               sa_seconds, sa_T, sa_xonly, fixy, ptarget_pct, targetJ, near_T, completeT, compl_cap, yardstick(), yardstick() / uniform_yardstick(), co, s2);
        for (int i = 0; i < nslots; i++) printf(" %d", gy[i]);
        printf("\n"); fflush(stdout);
        if (qyard_Q > 1) { double PQ = qyardstick(qyard_Q); printf("quotient yardstick of the start state on Z_%d x| Z_%d: P_Q = %.2f (ratio %.3f)\n", qyard_Q, M, PQ, PQ / uniform_yardstick()); fflush(stdout); }
        return sa(sa_seconds, sa_T, sa_xonly, near_T);
    }
    if (!have_init) { printf("FAIL: --init required\n"); return 2; }
    build_arcs();
    if (!ypattern_ok() || !simple_ok()) { printf("FAIL: start state is not a simple connected candidate\n"); return 1; }
    { int s2; int co = sphere2_coinc(&s2); printf("yardstick P of the start pattern: %.2f (ratio %.3f), radius-2 coincidences %d (sphere %d)\n", yardstick(), yardstick() / uniform_yardstick(), co, s2); }
    if (qyard_Q > 1) { double t0 = now(); double PQ = qyardstick(qyard_Q); printf("quotient yardstick on Z_%d x| Z_%d (fibre %d): P_Q = %.2f (ratio %.3f to uniform %.2f; %.1f s)\n", qyard_Q, M, N / qyard_Q, PQ, PQ / uniform_yardstick(), uniform_yardstick(), now() - t0); }
    if (xunit && !xpattern_ok()) printf("note: start state has a non-unit x-coordinate\n");
    double t0 = now();
    long start = eval((ll)1 << 50, 1);
    printf("start: defect %ld (%.3f s)\nBEST %ld :", start, now() - t0, start); print_gens(); fflush(stdout);
    if (do_eval) { if (start > 0 && start <= 64) for (int f = 0; f < nfib; f++) list_unreached(f); return 0; }

    if (nsl == 0) for (int i = 0; i < nslots; i++) slot_list[nsl++] = i;
    int sx[32], sy[32]; for (int i = 0; i < nslots; i++) { sx[i] = gx[i]; sy[i] = gy[i]; }
    unit_t *units = malloc(sizeof(unit_t) * nslots * M); int nunits = 0;
    for (int u = 0; u < nsl; u++) {
        int i = slot_list[u], ymax = slot_is_loop(i) ? (M - 1) / 2 : M - 1, ymin = (u == 0 ? yfrom : 0);
        if (xonly) { ymin = ymax = sy[i]; }
        for (int y = ymin; y <= ymax; y++) {
            gy[i] = y;
            double P = 0;
            if (order_yard) { if (ypattern_ok()) { build_slot_arcs(i); P = yardstick(); } else P = 1e300; }
            gy[i] = sy[i]; build_slot_arcs(i);
            units[nunits].slot = i; units[nunits].y = y; units[nunits].P = P; nunits++;
        }
    }
    if (order_yard) {
        qsort(units, nunits, sizeof(unit_t), cmp_unit);
        printf("units ordered by yardstick (%d): index slot y P\n", nunits);
        for (int u = 0; u < nunits; u++) printf("  U %d %d %d %.2f\n", u, units[u].slot, units[u].y, units[u].P);
    }
    if (to_u > nunits) to_u = nunits;
    ll best = start; long cands = 0, ties = 0, improved = 0, skipped = 0; double tlast = t0;
    int units_done = 0;
    for (int u = from_u; u < to_u; u++) {
        int i = units[u].slot, y = units[u].y;
        if (now() - t0 > seconds) {
            printf("RESUME slot %d y %d (unit index %d): time limit after %d units, %ld candidates, %.0f s\n", i, y, u, units_done, cands, now() - t0);
            break;
        }
        gy[i] = y;
        if (!ypattern_ok()) { skipped += N; gy[i] = sy[i]; units_done++; continue; }
        int xlo = 0, xhi = N - 1;
        if (slot_is_loop(i) && y == 0 && !xonly) { xlo = 1; xhi = (N - 1) / 2; }
        for (int x = xlo; x <= xhi; x++) {
            if (x == sx[i] && y == sy[i]) continue;
            gx[i] = x; build_slot_arcs(i);
            if (!simple_ok_slot(i)) { skipped++; continue; }
            long v = eval(best, 0); cands++;
            if (v <= best) {
                if (v < best) { improved++; printf("IMPROVED %ld (slot %d -> %d %d, %.0f s) :", v, i, x, y, now() - t0); print_gens(); }
                else { ties++; printf("TIE %ld (slot %d -> %d %d, %.0f s) :", v, i, x, y, now() - t0); print_gens(); }
                fflush(stdout);
                if (v == 0) {
                    printf("FOUND: re-evaluation with full profiles\n");
                    eval((ll)1 << 50, 1);
                    printf("FOUND N=%d M=%d r=%lld D=%d k=%d c=%d degree %d n=%lld%s :", N, M, R, D, K, C, deg, nn, invol ? " + involution (last pair)" : ""); print_gens();
                    fflush(stdout);
                    return 0;
                }
            }
        }
        gx[i] = sx[i]; gy[i] = sy[i]; build_slot_arcs(i);
        units_done++;
        double tn = now();
        if (tn - tlast > 30 || u + 1 == to_u || units[u + 1 < nunits ? u + 1 : u].slot != i) {
            printf("  unit %d (slot %d, y %d) done: %ld candidates, %ld skipped, %ld ties, %ld improved, %.0f s, %.0f cand/s\n",
                   u, i, y, cands, skipped, ties, improved, tn - t0, cands / (tn - t0 + 1e-9));
            fflush(stdout); tlast = tn;
        }
    }
    printf("done: %d of %d units [%d, %d), %ld candidates, %ld skipped, %ld BFS, %ld ties, %ld improved, %.0f s; best %lld\n",
           units_done, to_u - from_u, from_u, to_u, cands, skipped, nbfs, ties, improved, now() - t0, best);
    return 1;
}
