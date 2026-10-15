/* Simulated annealing for "pattern A" extensions of the generalized hexagon H(3):
 *   subdivide a set S of m incidence edges (new vertex x_e adjacent to both ends of e),
 *   optionally add h "hub" vertices with no old neighbours,
 *   add a graph F on the new vertices with F-degree 2 at subdivision vertices and 4 at hubs,
 *   so the result is 4-regular on 728 + m + h vertices.
 * Cost = number of unordered vertex pairs at distance > 6 (exact, see evaluate()).
 *
 * build:  gcc -O2 -o orbit_sa orbit_sa.c -lm
 * usage:  ./orbit_sa H3.edges m h seed seconds        (H3.edges: the incidence graph of H(3), one "u v" per line, or an adjacency list) [--T0 x] [--T1 y] [--Tre z] [--cycles k] [--pS p] [--pR p]
 *            [--focus q] [--mind d] [--fixS file] [--fixF file] [--initS file] [--initF file]
 *            [--eval] [--tag name] [--outdir dir] [--keepgoing] [--restart-best]
 * Features: small State (fast memcpy), --Tre (reheat temperature for cycles >= 1),
 * --focus (defect-weighted choice of the new vertex moved), --mind (pairwise edge-distance
 * constraint on S), best state written whenever cost <= 12, SIGTERM/SIGINT end the run cleanly.
 * S file: lines "p l" (an edge of H, any order); F file: lines "i j" (indices of new vertices).
 *
 * Exactness of the cost: let G0 = H with S subdivided (no F edges). Old-old pairs with
 * d_G0 <= 6 are fine. Those with d_G0 > 6 are exactly the pairs all of whose H-geodesics
 * are damaged: distance-4 pairs with >= 3 S-edges on the geodesic, distance-5 pairs with
 * >= 2, distance-6 pairs with an S-edge on each of the 4 geodesics (any other H-path is
 * longer by >= 2 because the girth is 12, so it has length >= 8).  For such a pair (a,b),
 * d_G(a,b) = min_i d_G0(a,x_i) + d_G(x_i,b) over subdivision vertices x_i (x_i = first new
 * vertex on a shortest path).  New-old and new-new distances come from BFS in G from each
 * new vertex.
 */
#include <math.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NOLD 728
#define MAXNEW 48
#define MAXN (NOLD + MAXNEW)
#define NE 1456
#define NPAIR (NOLD * NOLD)
#define UNREACH 255
/* MAXBOO bounds the old-old candidate list per state.  Random states have nBoo ~ 200-400; the
 * State is memcpy'd on every move, so the bound is kept small (v1 used 70000 ints = 560 KB per
 * copy, which dominated the running time).  States exceeding it get cost BIGCOST and are rejected. */
#define MAXBOO 8192
#define MAXBAD 8192
#define MAXDEF 64
#define BIGCOST 1000000

static uint16_t adjH[NOLD][4];
static uint8_t degH[NOLD];
static uint8_t distH[NOLD][NOLD];
static uint16_t ep[NE], el[NE];
static int edge_of[NOLD][4];
static int nE = 0;
static int *P4, *P5, *P6;
static int off4[NE + 1], off5[NE + 1], off6[NE + 1];
static int *near2;           /* edges at edge-distance <= 2 from e (incl. e), CSR */
static int offnear[NE + 1];
static uint8_t edist[NE][NE];

/* Persistent per-pair damage counters for the CURRENT state's S (old-old "boo" candidates):
   c4g/c5g = number of S-edges on the unique geodesic of a distance-4/5 pair, c6g = four 4-bit
   counters (one per geodesic) for distance-6 pairs, inBooG = membership in cur's boo list.
   A trial S-move applies +-1 edge deltas (about 4k entries) and is committed or undone. */
static uint32_t stamp[NPAIR];
static uint8_t c4g[NPAIR], c5g[NPAIR], inBooG[NPAIR];
static uint16_t c6g[NPAIR];
static int touched[MAXNEW * 2000];
static int ntouched;
static uint32_t gen = 1;

typedef struct {
    int S[MAXNEW];
    uint8_t inS[NE];
    uint16_t fadj[MAXNEW][4];
    uint8_t fdeg[MAXNEW];
    int fea[2 * MAXNEW], feb[2 * MAXNEW], nF;
    uint16_t adjG0[MAXN][4];
    uint8_t degG0[MAXN];
    uint8_t d0new[MAXNEW][MAXN];     /* G0 distances from each new vertex (BFS limit 6; UNREACH beyond) */
    uint8_t dGnew[MAXNEW][MAXNEW];   /* exact G distances between new vertices (capped at UNREACH) */
    uint16_t bad[MAXBAD];            /* old vertices with d0new > 6, grouped by new vertex */
    int badoff[MAXNEW + 1];
    uint16_t booA[MAXBOO], booB[MAXBOO];
    int nBoo, booOverflow;
    int cost, c_no, c_nn, c_oo;
    int def[MAXNEW];             /* new-old + new-new defects involving each new vertex */
    int ndef;                    /* sum of def[] */
    /* explicit defect lists (first MAXDEF of each kind) for targeted repair moves */
    uint16_t dno_i[MAXDEF], dno_v[MAXDEF], dnn_i[MAXDEF], dnn_j[MAXDEF], doo_a[MAXDEF], doo_b[MAXDEF];
    int n_dno, n_dnn, n_doo;
} State;

static int ms, hb, nnew, N;
static int tdeg[MAXNEW];

/* ---------- RNG ---------- */
static uint64_t rng_s;
static inline uint64_t rnd64(void) {
    uint64_t z = (rng_s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static inline int rndint(int n) { return (int)(rnd64() % (uint64_t)n); }
static inline double rnd01(void) { return (rnd64() >> 11) * (1.0 / 9007199254740992.0); }

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ---------- loading H ---------- */
static void load_H(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }
    char line[4096];
    memset(degH, 0, sizeof degH);
    static uint8_t seen[NOLD][NOLD];
    memset(seen, 0, sizeof seen);
    while (fgets(line, sizeof line, f)) {
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == 0) continue;
        char *end;
        long v = strtol(p, &end, 10);
        if (end == p) continue;
        p = end;
        for (;;) {
            long w = strtol(p, &end, 10);
            if (end == p) break;
            p = end;
            if (v == w) continue;
            if (!seen[v][w]) {
                seen[v][w] = seen[w][v] = 1;
                if (degH[v] >= 4 || degH[w] >= 4) { fprintf(stderr, "degree > 4 in H\n"); exit(1); }
                adjH[v][degH[v]++] = (uint16_t)w;
                adjH[w][degH[w]++] = (uint16_t)v;
            }
        }
    }
    fclose(f);
    for (int v = 0; v < NOLD; v++)
        if (degH[v] != 4) { fprintf(stderr, "vertex %d has degree %d\n", v, degH[v]); exit(1); }
    /* edges */
    for (int u = 0; u < NOLD; u++)
        for (int k = 0; k < 4; k++) {
            int w = adjH[u][k];
            if (u < w) {
                ep[nE] = u; el[nE] = w;
                edge_of[u][k] = nE;
                for (int kk = 0; kk < 4; kk++) if (adjH[w][kk] == u) edge_of[w][kk] = nE;
                nE++;
            }
        }
    if (nE != NE) { fprintf(stderr, "expected %d edges, got %d\n", NE, nE); exit(1); }
}

static void bfs_H(int s, uint8_t *dist) {
    static int q[NOLD];
    memset(dist, 255, NOLD);
    int qh = 0, qt = 0;
    dist[s] = 0; q[qt++] = s;
    while (qh < qt) {
        int u = q[qh++];
        for (int k = 0; k < 4; k++) {
            int w = adjH[u][k];
            if (dist[w] == 255) { dist[w] = dist[u] + 1; q[qt++] = w; }
        }
    }
}

static int walk_edges(int from, int to, const uint8_t *dto, int *out) {
    /* unique geodesic from `from` to `to` (d <= 5); dto = distH[to]; returns edges used */
    int cnt = 0, v = from;
    while (v != to) {
        int found = -1;
        for (int k = 0; k < 4; k++) {
            int w = adjH[v][k];
            if (dto[w] == dto[v] - 1) { if (found >= 0) { fprintf(stderr, "non-unique geodesic\n"); exit(1); } found = k; }
        }
        if (found < 0) { fprintf(stderr, "walk failed\n"); exit(1); }
        out[cnt++] = edge_of[v][found];
        v = adjH[v][found];
    }
    return cnt;
}

static void precompute(void) {
    for (int a = 0; a < NOLD; a++) bfs_H(a, distH[a]);
    int *cnt4 = calloc(NE, sizeof(int)), *cnt5 = calloc(NE, sizeof(int)), *cnt6 = calloc(NE, sizeof(int));
    int eds[8];
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            off4[0] = off5[0] = off6[0] = 0;
            for (int e = 0; e < NE; e++) {
                off4[e + 1] = off4[e] + cnt4[e];
                off5[e + 1] = off5[e] + cnt5[e];
                off6[e + 1] = off6[e] + cnt6[e];
            }
            P4 = malloc(sizeof(int) * off4[NE]);
            P5 = malloc(sizeof(int) * off5[NE]);
            P6 = malloc(sizeof(int) * off6[NE]);
            memset(cnt4, 0, NE * sizeof(int)); memset(cnt5, 0, NE * sizeof(int)); memset(cnt6, 0, NE * sizeof(int));
        }
        for (int a = 0; a < NOLD; a++) {
            for (int b = a + 1; b < NOLD; b++) {
                int d = distH[a][b];
                int pid = a * NOLD + b;
                if (d <= 3) continue;
                if (d <= 5) {
                    int c = walk_edges(b, a, distH[a], eds);
                    if (c != d) { fprintf(stderr, "bad walk\n"); exit(1); }
                    for (int t = 0; t < c; t++) {
                        int e = eds[t];
                        if (d == 4) { if (pass) P4[off4[e] + cnt4[e]] = pid; cnt4[e]++; }
                        else { if (pass) P5[off5[e] + cnt5[e]] = pid; cnt5[e]++; }
                    }
                } else if (d == 6) {
                    for (int g = 0; g < 4; g++) {
                        int w0 = adjH[a][g];
                        if (distH[w0][b] != 5) { fprintf(stderr, "distance-6 structure violated\n"); exit(1); }
                        int c = walk_edges(b, w0, distH[w0], eds);
                        eds[c++] = edge_of[a][g];
                        for (int t = 0; t < c; t++) {
                            int e = eds[t];
                            if (pass) P6[off6[e] + cnt6[e]] = (pid << 2) | g;
                            cnt6[e]++;
                        }
                    }
                } else { fprintf(stderr, "diameter > 6 in H\n"); exit(1); }
            }
        }
    }
    fprintf(stderr, "precomp: P4 %d P5 %d P6 %d entries per edge (edge 0): %d %d %d\n",
            off4[NE], off5[NE], off6[NE], cnt4[0], cnt5[0], cnt6[0]);
    free(cnt4); free(cnt5); free(cnt6);
    /* edge distances and near lists */
    int tot = 0;
    for (int e = 0; e < NE; e++) {
        for (int f = 0; f < NE; f++) {
            int d = distH[ep[e]][ep[f]];
            int t;
            t = distH[ep[e]][el[f]]; if (t < d) d = t;
            t = distH[el[e]][ep[f]]; if (t < d) d = t;
            t = distH[el[e]][el[f]]; if (t < d) d = t;
            edist[e][f] = (uint8_t)d;
            if (d <= 2) tot++;
        }
    }
    near2 = malloc(sizeof(int) * tot);
    tot = 0;
    for (int e = 0; e < NE; e++) {
        offnear[e] = tot;
        for (int f = 0; f < NE; f++) if (edist[e][f] <= 2) near2[tot++] = f;
    }
    offnear[NE] = tot;
}

/* ---------- state construction ---------- */
static void build_G0(State *st) {
    for (int u = 0; u < NOLD; u++) {
        memcpy(st->adjG0[u], adjH[u], sizeof adjH[u]);
        st->degG0[u] = 4;
    }
    for (int i = 0; i < nnew; i++) st->degG0[NOLD + i] = 0;
    for (int i = 0; i < ms; i++) {
        int e = st->S[i];
        int p = ep[e], l = el[e], x = NOLD + i;
        for (int k = 0; k < 4; k++) if (st->adjG0[p][k] == l) st->adjG0[p][k] = (uint16_t)x;
        for (int k = 0; k < 4; k++) if (st->adjG0[l][k] == p) st->adjG0[l][k] = (uint16_t)x;
        st->adjG0[x][0] = (uint16_t)p; st->adjG0[x][1] = (uint16_t)l; st->degG0[x] = 2;
    }
}

static void bfs_lim(const uint16_t (*adj)[4], const uint8_t *deg, int s, uint8_t *dist, int limit) {
    static int q[MAXN];
    memset(dist, UNREACH, N);
    int qh = 0, qt = 0;
    dist[s] = 0; q[qt++] = s;
    while (qh < qt) {
        int u = q[qh++];
        int du = dist[u];
        if (du >= limit) continue;
        int d = deg[u];
        for (int k = 0; k < d; k++) {
            int w = adj[u][k];
            if (dist[w] == UNREACH) { dist[w] = (uint8_t)(du + 1); q[qt++] = w; }
        }
    }
}

static void edge_delta(int e, int delta) {
    for (int t = off4[e]; t < off4[e + 1]; t++) c4g[P4[t]] = (uint8_t)(c4g[P4[t]] + delta);
    for (int t = off5[e]; t < off5[e + 1]; t++) c5g[P5[t]] = (uint8_t)(c5g[P5[t]] + delta);
    for (int t = off6[e]; t < off6[e + 1]; t++) {
        int pid = P6[t] >> 2, g = P6[t] & 3;
        c6g[pid] = (uint16_t)(c6g[pid] + delta * (1 << (4 * g)));
    }
}

static inline int qual(int pid) {
    uint16_t m = c6g[pid];
    return c5g[pid] >= 2 || c4g[pid] >= 3 || ((m & 0xF) && (m & 0xF0) && (m & 0xF00) && (m & 0xF000));
}

/* record the pairs whose counters edge e touches (deduplicated by stamp == gen) */
static void touch_edge(int e) {
    for (int t = off4[e]; t < off4[e + 1]; t++) { int pid = P4[t]; if (stamp[pid] != gen) { stamp[pid] = gen; touched[ntouched++] = pid; } }
    for (int t = off5[e]; t < off5[e + 1]; t++) { int pid = P5[t]; if (stamp[pid] != gen) { stamp[pid] = gen; touched[ntouched++] = pid; } }
    for (int t = off6[e]; t < off6[e + 1]; t++) { int pid = P6[t] >> 2; if (stamp[pid] != gen) { stamp[pid] = gen; touched[ntouched++] = pid; } }
}

static void boo_append(State *st, int pid) {
    if (st->nBoo >= MAXBOO) { st->booOverflow = 1; return; }
    st->booA[st->nBoo] = (uint16_t)(pid / NOLD); st->booB[st->nBoo] = (uint16_t)(pid % NOLD); st->nBoo++;
}

/* from scratch: reset the global counters to st->S and build st's boo list (initial state, restarts) */
static void boo_init(State *st) {
    memset(c4g, 0, sizeof c4g); memset(c5g, 0, sizeof c5g); memset(c6g, 0, sizeof c6g); memset(inBooG, 0, sizeof inBooG);
    gen++; ntouched = 0;
    for (int i = 0; i < ms; i++) { edge_delta(st->S[i], +1); touch_edge(st->S[i]); }
    st->nBoo = 0; st->booOverflow = 0;
    for (int t = 0; t < ntouched; t++) {
        int pid = touched[t];
        if (qual(pid)) { boo_append(st, pid); inBooG[pid] = (uint8_t)!st->booOverflow; }
    }
}

/* trial S-move e_old -> e_new from cur: update the global counters and build tr's boo list.
   Must be followed by boo_commit() (move accepted) or boo_undo() (rejected). */
static void boo_move(const State *cur, State *tr, int e_old, int e_new) {
    edge_delta(e_old, -1); edge_delta(e_new, +1);
    gen++; ntouched = 0;
    touch_edge(e_old); touch_edge(e_new);
    tr->nBoo = 0; tr->booOverflow = 0;
    for (int t = 0; t < cur->nBoo; t++) {
        int pid = cur->booA[t] * NOLD + cur->booB[t];
        if (stamp[pid] == gen && !qual(pid)) continue;   /* dropped out */
        boo_append(tr, pid);
    }
    for (int t = 0; t < ntouched; t++) {
        int pid = touched[t];
        if (!inBooG[pid] && qual(pid)) boo_append(tr, pid);  /* newly qualified */
    }
}
static void boo_commit(void) { for (int t = 0; t < ntouched; t++) inBooG[touched[t]] = (uint8_t)qual(touched[t]); }
static void boo_undo(int e_old, int e_new) { edge_delta(e_old, +1); edge_delta(e_new, -1); }

/* G0 distances from every new vertex (depth 6) and the lists of old vertices beyond 6 ("bad"):
   only those can be new-old defects.  Hubs have no G0 neighbours: every old vertex is bad. */
static void compute_d0new(State *st) {
    int nb = 0;
    for (int i = 0; i < nnew; i++) {
        if (i < ms) bfs_lim(st->adjG0, st->degG0, NOLD + i, st->d0new[i], 6);
        else { memset(st->d0new[i], UNREACH, N); st->d0new[i][NOLD + i] = 0; }
        st->badoff[i] = nb;
        const uint8_t *d = st->d0new[i];
        for (int v = 0; v < NOLD; v++)
            if (d[v] > 6) { if (nb >= MAXBAD) { st->booOverflow = 1; break; } st->bad[nb++] = (uint16_t)v; }
    }
    st->badoff[nnew] = nb;
}

/* Exact cost from d0new alone (no BFS in G):
   - dGnew = all-pairs shortest paths on the new vertices with w(i,j) = 1 if (i,j) in F, else
     d_G0(x_i,x_j): every G path between new vertices splits at its new-vertex visits into F edges
     and G0 segments, so Floyd-Warshall on this weighted graph gives d_G exactly;
   - d_G(x_i, v) = min_l dGnew(i,l) + d_G0(x_l, v)  (l = last new vertex before v);
   - d_G(a, b) = min_{i,l} d_G0(a, x_i) + dGnew(i,l) + d_G0(x_l, b) for old-old pairs with d_G0 > 6
     (every path of length <= 6 uses an F edge, so it visits new vertices; i = first, l = last).
   d0new values beyond 6 are UNREACH; all sums below stay > 6 then, which is all that matters. */
static int evaluate(State *st) {
    static int d[MAXNEW][MAXNEW];
    int n = nnew;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) d[i][j] = st->d0new[i][NOLD + j];
        d[i][i] = 0;
        for (int k = 0; k < st->fdeg[i]; k++) d[i][st->fadj[i][k]] = 1;
    }
    for (int k = 0; k < n; k++)
        for (int i = 0; i < n; i++) {
            int dik = d[i][k];
            if (dik >= UNREACH) continue;
            for (int j = 0; j < n; j++) {
                int v = dik + d[k][j];
                if (v < d[i][j]) d[i][j] = v;
            }
        }
    int cno = 0, cnn = 0, coo = 0;
    st->n_dno = st->n_dnn = st->n_doo = 0;
    for (int i = 0; i < n; i++) st->def[i] = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            st->dGnew[i][j] = (uint8_t)(d[i][j] > UNREACH ? UNREACH : d[i][j]);
            if (j > i && d[i][j] > 6) {
                cnn++; st->def[i]++; st->def[j]++;
                if (st->n_dnn < MAXDEF) { st->dnn_i[st->n_dnn] = (uint16_t)i; st->dnn_j[st->n_dnn] = (uint16_t)j; st->n_dnn++; }
            }
        }
    for (int i = 0; i < n; i++) {
        int c = 0;
        for (int t = st->badoff[i]; t < st->badoff[i + 1]; t++) {
            int v = st->bad[t];
            int best = UNREACH;
            for (int l = 0; l < n; l++) {
                int w = d[i][l] + st->d0new[l][v];
                if (w < best) best = w;
            }
            if (best > 6) {
                c++;
                if (st->n_dno < MAXDEF) { st->dno_i[st->n_dno] = (uint16_t)i; st->dno_v[st->n_dno] = (uint16_t)v; st->n_dno++; }
            }
        }
        cno += c; st->def[i] += c;
    }
    for (int t = 0; t < st->nBoo; t++) {
        int a = st->booA[t], b = st->booB[t];
        int best = UNREACH;
        for (int i = 0; i < n; i++) {
            int da = st->d0new[i][a];
            if (da > 5) continue;
            for (int l = 0; l < n; l++) {
                int w = da + d[i][l] + st->d0new[l][b];
                if (w < best) best = w;
            }
        }
        if (best > 6) {
            coo++;
            if (st->n_doo < MAXDEF) { st->doo_a[st->n_doo] = (uint16_t)a; st->doo_b[st->n_doo] = (uint16_t)b; st->n_doo++; }
        }
    }
    st->ndef = 0;
    for (int i = 0; i < n; i++) st->ndef += st->def[i];
    st->c_no = cno; st->c_nn = cnn; st->c_oo = coo;
    st->cost = st->booOverflow ? BIGCOST : cno + cnn + coo;
    return st->cost;
}

/* from scratch (initial state, restart from best): also resets the global damage counters to st */
static void full_rebuild(State *st) {
    build_G0(st);
    boo_init(st);
    compute_d0new(st);
    evaluate(st);
}

/* tr = cur with S[i] replaced (e_old -> e_new, already written into tr->S/inS) */
static void rebuild_smove(const State *cur, State *tr, int e_old, int e_new) {
    build_G0(tr);
    boo_move(cur, tr, e_old, e_new);
    compute_d0new(tr);
    evaluate(tr);
}

/* exact check by all-pairs BFS on the explicit graph; returns #pairs > 6 */
static int exact_cost(const State *st, int *maxdeg, int *mindeg) {
    static uint16_t adjG[MAXN][4];
    static uint8_t degG[MAXN];
    memcpy(adjG, st->adjG0, sizeof(uint16_t) * 4 * N);
    memcpy(degG, st->degG0, N);
    for (int i = 0; i < nnew; i++) {
        int x = NOLD + i;
        for (int k = 0; k < st->fdeg[i]; k++) adjG[x][degG[x]++] = (uint16_t)(NOLD + st->fadj[i][k]);
    }
    *maxdeg = 0; *mindeg = 99;
    for (int v = 0; v < N; v++) { if (degG[v] > *maxdeg) *maxdeg = degG[v]; if (degG[v] < *mindeg) *mindeg = degG[v]; }
    static uint8_t dist[MAXN];
    int bad = 0;
    for (int s = 0; s < N; s++) {
        bfs_lim(adjG, degG, s, dist, 20);
        for (int t = s + 1; t < N; t++) bad += (dist[t] > 6);
    }
    return bad;
}

static void write_edges(const State *st, const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }
    for (int u = 0; u < N; u++)
        for (int k = 0; k < st->degG0[u]; k++) {
            int w = st->adjG0[u][k];
            if (u < w) fprintf(f, "%d %d\n", u, w);
        }
    for (int i = 0; i < nnew; i++)
        for (int k = 0; k < st->fdeg[i]; k++) {
            int j = st->fadj[i][k];
            if (i < j) fprintf(f, "%d %d\n", NOLD + i, NOLD + j);
        }
    fclose(f);
}

static void print_state(const State *st, FILE *f) {
    fprintf(f, "S:");
    for (int i = 0; i < ms; i++) fprintf(f, " (%d,%d)", ep[st->S[i]], el[st->S[i]]);
    fprintf(f, "\nF:");
    for (int i = 0; i < nnew; i++)
        for (int k = 0; k < st->fdeg[i]; k++) if (i < st->fadj[i][k]) fprintf(f, " (%d,%d)", i, st->fadj[i][k]);
    fprintf(f, "\n");
}

/* ---------- F helpers ---------- */
static int f_has(const State *st, int a, int b) {
    for (int k = 0; k < st->fdeg[a]; k++) if (st->fadj[a][k] == b) return 1;
    return 0;
}
static void f_add(State *st, int a, int b) {
    st->fadj[a][st->fdeg[a]++] = (uint16_t)b;
    st->fadj[b][st->fdeg[b]++] = (uint16_t)a;
    st->fea[st->nF] = a < b ? a : b; st->feb[st->nF] = a < b ? b : a; st->nF++;
}
static void f_del(State *st, int a, int b) {
    for (int k = 0; k < st->fdeg[a]; k++) if (st->fadj[a][k] == b) { st->fadj[a][k] = st->fadj[a][--st->fdeg[a]]; break; }
    for (int k = 0; k < st->fdeg[b]; k++) if (st->fadj[b][k] == a) { st->fadj[b][k] = st->fadj[b][--st->fdeg[b]]; break; }
    int lo = a < b ? a : b, hi = a < b ? b : a;
    for (int t = 0; t < st->nF; t++) if (st->fea[t] == lo && st->feb[t] == hi) { st->fea[t] = st->fea[st->nF - 1]; st->feb[t] = st->feb[st->nF - 1]; st->nF--; break; }
}

static void random_F(State *st) {
    /* configuration model with rejection */
    for (;;) {
        for (int i = 0; i < nnew; i++) st->fdeg[i] = 0;
        st->nF = 0;
        int stubs[4 * MAXNEW], ns = 0;
        for (int i = 0; i < nnew; i++) for (int k = 0; k < tdeg[i]; k++) stubs[ns++] = i;
        int ok = 1;
        while (ns > 0 && ok) {
            int tries = 0;
            for (;;) {
                int ia = rndint(ns), ib = rndint(ns);
                int a = stubs[ia], b = stubs[ib];
                if (a != b && !f_has(st, a, b)) {
                    f_add(st, a, b);
                    if (ia > ib) { int t = ia; ia = ib; ib = t; }
                    stubs[ib] = stubs[ns - 1]; ns--;
                    stubs[ia] = stubs[ns - 1]; ns--;
                    break;
                }
                if (++tries > 1000) { ok = 0; break; }
            }
        }
        if (ok) return;
    }
}

static void random_S(State *st) {
    memset(st->inS, 0, sizeof st->inS);
    for (int i = 0; i < ms; i++) {
        int e;
        do { e = rndint(NE); } while (st->inS[e]);
        st->inS[e] = 1; st->S[i] = e;
    }
}

static int read_pairs(const char *path, int *A, int *B, int max) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }
    int n = 0; char line[256];
    while (fgets(line, sizeof line, f)) {
        int a, b;
        if (line[0] == '#') continue;
        if (sscanf(line, "%d %d", &a, &b) == 2) { if (n >= max) break; A[n] = a; B[n] = b; n++; }
    }
    fclose(f);
    return n;
}

static int find_edge(int u, int v) {
    for (int k = 0; k < 4; k++) if (adjH[u][k] == v) return edge_of[u][k];
    return -1;
}

/* Pick a new vertex: with probability pfocus (and if any new-old/new-new defect exists) proportional
   to its defect count, otherwise uniformly among the first n new vertices. */
static double pfocus = 0.0;
static int pick_new(const State *t, int n) {
    if (pfocus > 0 && t->ndef > 0 && rnd01() < pfocus) {
        int r = rndint(t->ndef);
        for (int i = 0; i < nnew; i++) {
            r -= t->def[i];
            if (r < 0) return i < n ? i : rndint(n);
        }
    }
    return rndint(n);
}

/* Pick an F edge: with probability pfocus one incident to a defect-weighted new vertex. */
static int pick_fedge(const State *t) {
    if (pfocus > 0 && t->ndef > 0 && rnd01() < pfocus) {
        int i = pick_new(t, nnew);
        if (t->fdeg[i] > 0) {
            int j = t->fadj[i][rndint(t->fdeg[i])];
            int lo = i < j ? i : j, hi = i < j ? j : i;
            for (int k = 0; k < t->nF; k++) if (t->fea[k] == lo && t->feb[k] == hi) return k;
        }
    }
    return rndint(t->nF);
}

/* random F move on t (relocate with probability prel, else 2-swap); returns 1 if a move was
   applied and evaluated, 0 if the drawn move was invalid (t may then be inconsistent: discard it). */
static int f_move(State *t, double prel) {
    if (rnd01() < prel) {
        /* relocate a degree-2 vertex i: (a,i),(i,b) -> (a,b); (c,d) -> (c,i),(i,d) */
        int i = pick_new(t, ms);
        if (t->fdeg[i] != 2) return 0;
        int a = t->fadj[i][0], b = t->fadj[i][1];
        if (f_has(t, a, b)) return 0;
        int k = rndint(t->nF);
        int c = t->fea[k], d = t->feb[k];
        if (c == i || d == i || (c == a && d == b) || (c == b && d == a)) return 0;
        f_del(t, a, i); f_del(t, i, b); f_add(t, a, b);
        f_del(t, c, d);
        if (f_has(t, c, i) || f_has(t, i, d)) return 0;
        f_add(t, c, i); f_add(t, i, d);
        evaluate(t);
        return 1;
    }
    /* 2-swap: (a,b),(c,d) -> (a,c),(b,d) or (a,d),(b,c) */
    int t1 = pick_fedge(t), t2 = rndint(t->nF);
    if (t1 == t2) return 0;
    int a = t->fea[t1], b = t->feb[t1], c = t->fea[t2], d = t->feb[t2];
    if (rnd01() < 0.5) { int x = c; c = d; d = x; }
    if (a == c || a == d || b == c || b == d || f_has(t, a, c) || f_has(t, b, d)) return 0;
    f_del(t, a, b); f_del(t, c, d); f_add(t, a, c); f_add(t, b, d);
    evaluate(t);
    return 1;
}

/* ---------- targeted repair moves (WalkSAT style) ----------
   Pick one current defect and propose a move that repairs it (other defects may appear; the
   Metropolis test decides).  New-old (x_i,v): either make x_i F-adjacent to a new vertex l with
   d_G0(x_l,v) <= 5, or relocate an S-edge lying on the damaged H-geodesic from the nearer end of
   e_i to v.  New-new (x_i,x_j): make x_i F-adjacent to l with d_G(x_l,x_j) <= 5, or relocate e_i/e_j.
   Old-old (a,b): relocate an S-edge on one of the damaged geodesics, or make x_i (near a) F-adjacent
   to x_l (near b) with d_G0(a,x_i) + 1 + d_G0(x_l,b) <= 6. */
static int mind;
static int s_ok(const State *st, int skip, int e);

/* 2-swap that creates the F edge (i,l): (i,a),(l,b) -> (i,l),(a,b) */
static int f_make_edge(State *t, int i, int l) {
    if (i == l || f_has(t, i, l) || t->fdeg[i] == 0 || t->fdeg[l] == 0) return 0;
    int a = t->fadj[i][rndint(t->fdeg[i])], b = t->fadj[l][rndint(t->fdeg[l])];
    if (a == b || a == l || b == i || f_has(t, a, b)) return 0;
    f_del(t, i, a); f_del(t, l, b); f_add(t, i, l); f_add(t, a, b);
    return 1;
}

/* S indices of the S-edges on the unique H-geodesic from u to v (d(u,v) <= 5); returns count */
static int s_on_geodesic(const State *t, int u, int v, int *out) {
    int cnt = 0, x = u;
    const uint8_t *dto = distH[v];
    while (x != v) {
        int found = -1;
        for (int k = 0; k < 4; k++) if (dto[adjH[x][k]] == dto[x] - 1) { found = k; break; }
        if (found < 0) return cnt;
        int e = edge_of[x][found];
        if (t->inS[e]) for (int j = 0; j < ms; j++) if (t->S[j] == e) { out[cnt++] = j; break; }
        x = adjH[x][found];
    }
    return cnt;
}

/* relocate S-edge j to a random admissible edge; performs the trial rebuild */
static int s_relocate(const State *cur, State *t, int j, int *e_old, int *e_new) {
    int e;
    if (rnd01() < 0.5) { int old = t->S[j]; e = near2[offnear[old] + rndint(offnear[old + 1] - offnear[old])]; }
    else e = rndint(NE);
    if (t->inS[e] || !s_ok(t, j, e)) return 0;
    *e_old = t->S[j]; *e_new = e;
    t->inS[*e_old] = 0; t->S[j] = e; t->inS[e] = 1;
    rebuild_smove(cur, t, *e_old, *e_new);
    return 1;
}

static int repair_move(const State *cur, State *t, int *smove, int *e_old, int *e_new) {
    int tot = t->n_dno + t->n_dnn + t->n_doo;
    if (tot == 0) return 0;
    int r = rndint(tot);
    int cand[MAXNEW * MAXNEW], nc = 0, sidx[8];
    *smove = 0;
    if (r < t->n_dno) {
        int i = t->dno_i[r], v = t->dno_v[r];
        if (rnd01() < 0.5 || i >= ms) {
            for (int l = 0; l < nnew; l++) if (l != i && t->d0new[l][v] <= 5 && !f_has(t, i, l)) cand[nc++] = l;
            if (nc == 0 || !f_make_edge(t, i, cand[rndint(nc)])) return 0;
            evaluate(t);
            return 1;
        }
        int p = ep[t->S[i]], l = el[t->S[i]];
        int u = distH[p][v] < distH[l][v] ? p : l;
        int ns = s_on_geodesic(t, u, v, sidx);
        if (ns == 0) return 0;
        *smove = 1;
        return s_relocate(cur, t, sidx[rndint(ns)], e_old, e_new);
    }
    r -= t->n_dno;
    if (r < t->n_dnn) {
        int i = t->dnn_i[r], j = t->dnn_j[r];
        if (rnd01() < 0.5) { int x = i; i = j; j = x; }
        if (rnd01() < 0.5 || i >= ms) {
            for (int l = 0; l < nnew; l++) if (l != i && t->dGnew[l][j] <= 5 && !f_has(t, i, l)) cand[nc++] = l;
            if (nc == 0 || !f_make_edge(t, i, cand[rndint(nc)])) return 0;
            evaluate(t);
            return 1;
        }
        *smove = 1;
        return s_relocate(cur, t, i, e_old, e_new);
    }
    r -= t->n_dnn;
    int a = t->doo_a[r], b = t->doo_b[r];
    if (rnd01() < 0.5) {
        for (int i = 0; i < nnew; i++) {
            if (t->d0new[i][a] > 4) continue;
            for (int l = 0; l < nnew; l++)
                if (l != i && t->d0new[i][a] + 1 + t->d0new[l][b] <= 6 && !f_has(t, i, l)) cand[nc++] = i * MAXNEW + l;
        }
        if (nc == 0) return 0;
        int c = cand[rndint(nc)];
        if (!f_make_edge(t, c / MAXNEW, c % MAXNEW)) return 0;
        evaluate(t);
        return 1;
    }
    int u = a, ns;
    if (distH[a][b] == 6) u = adjH[a][rndint(4)];   /* pick one of the 4 geodesics */
    ns = s_on_geodesic(t, u, b, sidx);
    if (u != a) { int e = find_edge(a, u); if (t->inS[e]) for (int j = 0; j < ms; j++) if (t->S[j] == e) { sidx[ns++] = j; break; } }
    if (ns == 0) return 0;
    *smove = 1;
    return s_relocate(cur, t, sidx[rndint(ns)], e_old, e_new);
}

/* ---------- exhaustive 2-move search from the initial state (--exh2) ----------
   Enumerates every F move (2-swaps and relocations) on top of (a) the state itself and (b) every
   admissible single S-move; reports the minimum cost found and stops at the first cost-0 state.
   Complete for the neighbourhood {F} + {S x F}; S x S pairs are not enumerated. */
static volatile int stop_flag = 0;
static long exh_evals;
static int exh_best, exh_prune = 8;   /* S-moves raising the cost by more than exh_prune are not combined with F moves */
static State *exh_tmp;
static int exh_apply_f_moves(const State *base, State *out) {
    /* returns 1 if a cost-0 state was found (copied to out) */
    for (int t1 = 0; t1 < base->nF; t1++)
        for (int t2 = t1 + 1; t2 < base->nF; t2++)
            for (int var = 0; var < 2; var++) {
                int a = base->fea[t1], b = base->feb[t1], c = base->fea[t2], d = base->feb[t2];
                if (var) { int x = c; c = d; d = x; }
                if (a == c || a == d || b == c || b == d || f_has(base, a, c) || f_has(base, b, d)) continue;
                memcpy(exh_tmp, base, sizeof(State));
                f_del(exh_tmp, a, b); f_del(exh_tmp, c, d); f_add(exh_tmp, a, c); f_add(exh_tmp, b, d);
                evaluate(exh_tmp); exh_evals++;
                if (exh_tmp->cost < exh_best) exh_best = exh_tmp->cost;
                if (exh_tmp->cost == 0) { memcpy(out, exh_tmp, sizeof(State)); return 1; }
            }
    for (int i = 0; i < ms; i++) {
        if (base->fdeg[i] != 2) continue;
        int a = base->fadj[i][0], b = base->fadj[i][1];
        if (f_has(base, a, b)) continue;
        for (int k = 0; k < base->nF; k++) {
            int c = base->fea[k], d = base->feb[k];
            if (c == i || d == i || (c == a && d == b) || (c == b && d == a)) continue;
            memcpy(exh_tmp, base, sizeof(State));
            f_del(exh_tmp, a, i); f_del(exh_tmp, i, b); f_add(exh_tmp, a, b); f_del(exh_tmp, c, d);
            if (f_has(exh_tmp, c, i) || f_has(exh_tmp, i, d)) continue;
            f_add(exh_tmp, c, i); f_add(exh_tmp, i, d);
            evaluate(exh_tmp); exh_evals++;
            if (exh_tmp->cost < exh_best) exh_best = exh_tmp->cost;
            if (exh_tmp->cost == 0) { memcpy(out, exh_tmp, sizeof(State)); return 1; }
        }
    }
    return 0;
}

static int exhaustive2(State *cur, State *trial, State *out, double t_start) {
    exh_tmp = calloc(1, sizeof(State));
    exh_evals = 0; exh_best = cur->cost;
    printf("exh2: F moves on the initial state\n"); fflush(stdout);
    if (exh_apply_f_moves(cur, out)) return 1;
    printf("exh2: F-only best %d after %ld evaluations; now S x F\n", exh_best, exh_evals); fflush(stdout);
    for (int i = 0; i < ms; i++) {
        for (int e = 0; e < NE; e++) {
            if (cur->inS[e] || !s_ok(cur, i, e)) continue;
            memcpy(trial, cur, sizeof(State));
            int e_old = trial->S[i];
            trial->inS[e_old] = 0; trial->S[i] = e; trial->inS[e] = 1;
            rebuild_smove(cur, trial, e_old, e);
            exh_evals++;
            if (trial->cost < exh_best) exh_best = trial->cost;
            int found = trial->cost == 0 ? (memcpy(out, trial, sizeof(State)), 1) : 0;
            if (!found && trial->cost <= cur->cost + exh_prune) found = exh_apply_f_moves(trial, out);
            boo_undo(e_old, e);
            if (found) return 1;
            if (stop_flag) return 0;
        }
        printf("exh2: S index %d done, t=%.0f, evaluations %ld, best %d\n", i, now_sec() - t_start, exh_evals, exh_best);
        fflush(stdout);
    }
    return 0;
}

/* ---------- batch evaluation (--Flist FILE) ----------
   FILE has lines "S: p l p l ..." (m flags of H) and "F: i j i j ..." (a 2-regular graph on the m
   new vertices, F-degree 2 everywhere).  F lines before the first S line form the COMMON list,
   evaluated for every S; F lines after an S line are evaluated for that S only.
   For every S prints the number of F evaluated, the minimum cost, how many F reach it, and the
   best F; a cost-0 F is re-checked by exact_cost() and written as an edge file. */
#define FL_MAX_COMMON 20000
static int fl_common[FL_MAX_COMMON][2 * MAXNEW], fl_common_n[FL_MAX_COMMON], fl_ncommon = 0;

static int fl_parse_F(const char *p, int *pairs) {
    int n = 0;
    for (;;) {
        char *end; long a = strtol(p, &end, 10); if (end == p) break; p = end;
        long b = strtol(p, &end, 10); if (end == p) break; p = end;
        if (n >= 2 * MAXNEW) return -1;
        pairs[n++] = (int)a; pairs[n++] = (int)b;
    }
    return n / 2;
}

static int fl_set_F(State *st, const int *pairs, int np) {
    for (int i = 0; i < nnew; i++) st->fdeg[i] = 0;
    st->nF = 0;
    for (int t = 0; t < np; t++) {
        int a = pairs[2 * t], b = pairs[2 * t + 1];
        if (a < 0 || b < 0 || a >= nnew || b >= nnew || a == b || f_has(st, a, b) || st->fdeg[a] >= 4 || st->fdeg[b] >= 4) return 0;
        f_add(st, a, b);
    }
    for (int i = 0; i < nnew; i++) if (st->fdeg[i] != tdeg[i]) return 0;
    return 1;
}

static int run_flist(const char *path, const char *outdir, const char *tag) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); return 1; }
    State *st = calloc(1, sizeof(State));
    static char line[1 << 16];
    static char curS[1 << 14], bestF[1 << 14];
    int nS = 0, nF = 0, best = BIGCOST, nbest = 0, have_S = 0, nzero_total = 0, nbad = 0;
    long total_F = 0;
    int pairs[2 * MAXNEW];
    double t0 = now_sec();
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'S') {
            if (have_S) {
                printf("S#%d [%s] F %d min %d (x%d) bestF [%s]\n", nS, curS, nF, best, nbest, bestF);
                if (nS % 500 == 0) { fprintf(stderr, "Flist: %d S done, %ld F evaluations, %.0f s\n", nS, total_F, now_sec() - t0); }
                fflush(stdout);
            }
            nS++; nF = 0; best = BIGCOST; nbest = 0; bestF[0] = 0;
            memset(st->inS, 0, sizeof st->inS);
            int k = 0;
            char *p = line + 2;
            for (;;) {
                char *end; long a = strtol(p, &end, 10); if (end == p) break; p = end;
                long b = strtol(p, &end, 10); if (end == p) break; p = end;
                int e = find_edge((int)a, (int)b);
                if (e < 0) { fprintf(stderr, "(%ld,%ld) is not an edge of H\n", a, b); return 1; }
                if (k >= ms) { fprintf(stderr, "S line has more than m = %d flags\n", ms); return 1; }
                st->S[k++] = e; st->inS[e] = 1;
            }
            if (k != ms) { fprintf(stderr, "S line has %d flags, m = %d\n", k, ms); return 1; }
            snprintf(curS, sizeof curS, "%.16000s", line + 2);
            curS[strcspn(curS, "\n")] = 0;
            for (int i = 0; i < nnew; i++) st->fdeg[i] = 0;
            st->nF = 0;
            full_rebuild(st);
            have_S = 1;
            /* common F list */
            for (int c = 0; c < fl_ncommon; c++) {
                if (!fl_set_F(st, fl_common[c], fl_common_n[c])) { nbad++; continue; }
                evaluate(st);
                nF++; total_F++;
                if (st->cost < best) {
                    best = st->cost; nbest = 0;
                    int pos = 0; bestF[0] = 0;
                    for (int t = 0; t < fl_common_n[c] && pos < 16000; t++) pos += snprintf(bestF + pos, sizeof bestF - pos, "%d %d ", fl_common[c][2 * t], fl_common[c][2 * t + 1]);
                }
                if (st->cost == best) nbest++;
                if (st->cost == 0) {
                    int mx, mn, ex = exact_cost(st, &mx, &mn);
                    nzero_total++;
                    char path2[1024];
                    snprintf(path2, sizeof path2, "%s/%s_n%d_S%d_zero%d.edges", outdir, tag, N, nS, nzero_total);
                    write_edges(st, path2);
                    printf("ZERO COST in Flist: S#%d exact=%d maxdeg %d mindeg %d n=%d wrote %s\n", nS, ex, mx, mn, N, path2);
                    print_state(st, stdout);
                    fflush(stdout);
                }
            }
        } else if (line[0] == 'F') {
            int np = fl_parse_F(line + 2, pairs);
            if (np < 0) { fprintf(stderr, "F line too long (skipped)\n"); nbad++; continue; }
            if (!have_S) {
                if (fl_ncommon >= FL_MAX_COMMON) { fprintf(stderr, "too many common F lines\n"); return 1; }
                memcpy(fl_common[fl_ncommon], pairs, sizeof(int) * 2 * np);
                fl_common_n[fl_ncommon++] = np;
                continue;
            }
            if (!fl_set_F(st, pairs, np)) { nbad++; continue; }
            evaluate(st);
            nF++; total_F++;
            if (st->cost < best) { best = st->cost; nbest = 0; snprintf(bestF, sizeof bestF, "%.16000s", line + 2); bestF[strcspn(bestF, "\n")] = 0; }
            if (st->cost == best) nbest++;
            if (st->cost == 0) {
                int mx, mn, ex = exact_cost(st, &mx, &mn);
                nzero_total++;
                char path2[1024];
                snprintf(path2, sizeof path2, "%s/%s_n%d_S%d_zero%d.edges", outdir, tag, N, nS, nzero_total);
                write_edges(st, path2);
                printf("ZERO COST in Flist: S#%d exact=%d maxdeg %d mindeg %d n=%d wrote %s\n", nS, ex, mx, mn, N, path2);
                print_state(st, stdout);
                fflush(stdout);
            }
        }
    }
    if (have_S) printf("S#%d [%s] F %d min %d (x%d) bestF [%s]\n", nS, curS, nF, best, nbest, bestF);
    printf("Flist done: %d S blocks, %d common F, %ld F evaluations, %d invalid F skipped, %d cost-0 states, %.1f s\n",
           nS, fl_ncommon, total_F, nbad, nzero_total, now_sec() - t0);
    fclose(f);
    free(st);
    return 0;
}

/* ---------- main ---------- */
static void on_signal(int sig) { (void)sig; stop_flag = 1; }

static void write_best_file(const State *best, const char *outdir, const char *tag, int ex) {
    char path[1024];
    snprintf(path, sizeof path, "%s/%s_best.txt", outdir, tag);
    FILE *f = fopen(path, "w");
    if (f) { fprintf(f, "# cost %d exact %d n %d\n", best->cost, ex, N); print_state(best, f); fclose(f); }
}

/* S edges must be pairwise at edge-distance >= mind (all known cost-0 states have min distance 2) */
static int s_ok(const State *st, int skip, int e) {
    if (mind <= 0) return 1;
    for (int j = 0; j < ms; j++) if (j != skip && edist[st->S[j]][e] < mind) return 0;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s H3_adjlist m h seed seconds [opts]\n", argv[0]); return 1; }
    const char *hpath = argv[1];
    ms = atoi(argv[2]); hb = atoi(argv[3]);
    uint64_t seed = strtoull(argv[4], 0, 10);
    double seconds = atof(argv[5]);
    double T0 = 2.0, T1 = 0.15, pS = 0.3, pR = 0.15, Tre = -1, pRep = 0.0;
    int cycles = 1, do_eval = 0, keepgoing = 0, restart_best = 0, check_every = 0, exh2 = 0;
    long nchecks = 0;
    const char *fixS = 0, *fixF = 0, *initS = 0, *initF = 0, *tag = "run", *outdir = ".", *flist = 0;
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "--Flist")) flist = argv[++i];   /* batch evaluation of S/F lines, see run_flist() */
        else if (!strcmp(argv[i], "--T0")) T0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--T1")) T1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--Tre")) Tre = atof(argv[++i]);    /* start temperature of cycles >= 1 (default T0) */
        else if (!strcmp(argv[i], "--focus")) pfocus = atof(argv[++i]); /* prob. of defect-focused vertex choice */
        else if (!strcmp(argv[i], "--mind")) mind = atoi(argv[++i]);  /* min pairwise edge-distance in S */
        else if (!strcmp(argv[i], "--check")) check_every = atoi(argv[++i]); /* self-check every k iterations */
        else if (!strcmp(argv[i], "--prep")) pRep = atof(argv[++i]);   /* prob. of a targeted repair move */
        else if (!strcmp(argv[i], "--cycles")) cycles = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--pS")) pS = atof(argv[++i]);
        else if (!strcmp(argv[i], "--pR")) pR = atof(argv[++i]);
        else if (!strcmp(argv[i], "--fixS")) fixS = argv[++i];
        else if (!strcmp(argv[i], "--fixF")) fixF = argv[++i];
        else if (!strcmp(argv[i], "--initS")) initS = argv[++i];   /* seed S with <= m edges, rest random */
        else if (!strcmp(argv[i], "--initF")) initF = argv[++i];   /* seed F on the first k new vertices; others inserted */
        else if (!strcmp(argv[i], "--eval")) do_eval = 1;
        else if (!strcmp(argv[i], "--exh2")) exh2 = 1;   /* exhaustive {F} + {S x F} 2-move search from the initial state */
        else if (!strcmp(argv[i], "--exh2prune")) exh_prune = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--tag")) tag = argv[++i];
        else if (!strcmp(argv[i], "--outdir")) outdir = argv[++i];
        else if (!strcmp(argv[i], "--keepgoing")) keepgoing = 1;
        else if (!strcmp(argv[i], "--restart-best")) restart_best = 1; /* each cycle restarts from the best state */
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (fixS) initS = fixS;
    if (fixF) initF = fixF;
    nnew = ms + hb; N = NOLD + nnew;
    if (nnew > MAXNEW || ms < 1) { fprintf(stderr, "bad m/h\n"); return 1; }
    for (int i = 0; i < nnew; i++) tdeg[i] = i < ms ? 2 : 4;
    rng_s = seed * 0x9E3779B97F4A7C15ULL + 12345;
    double t_start = now_sec();
    load_H(hpath);
    precompute();
    fprintf(stderr, "precompute done in %.2fs\n", now_sec() - t_start);
    if (flist) return run_flist(flist, outdir, tag);

    State *cur = calloc(1, sizeof(State)), *trial = calloc(1, sizeof(State)), *best = calloc(1, sizeof(State));
    if (initS) {
        int A[MAXNEW], B[MAXNEW];
        int n = read_pairs(initS, A, B, MAXNEW);
        if (n > ms || (fixS && n != ms)) { fprintf(stderr, "S file has %d edges, m = %d\n", n, ms); return 1; }
        memset(cur->inS, 0, sizeof cur->inS);
        for (int i = 0; i < n; i++) {
            int e = find_edge(A[i], B[i]);
            if (e < 0) { fprintf(stderr, "(%d,%d) is not an edge of H\n", A[i], B[i]); return 1; }
            if (cur->inS[e]) { fprintf(stderr, "duplicate S edge\n"); return 1; }
            cur->S[i] = e; cur->inS[e] = 1;
        }
        for (int i = n; i < ms; i++) {
            int e;
            do { e = rndint(NE); } while (cur->inS[e]);
            cur->inS[e] = 1; cur->S[i] = e;
        }
        printf("S seeded with %d edges from %s, %d random\n", n, initS, ms - n);
    } else random_S(cur);
    if (initF) {
        int A[2 * MAXNEW], B[2 * MAXNEW];
        int n = read_pairs(initF, A, B, 2 * MAXNEW);
        for (int i = 0; i < nnew; i++) cur->fdeg[i] = 0;
        cur->nF = 0;
        for (int i = 0; i < n; i++) {
            if (A[i] < 0 || B[i] < 0 || A[i] >= nnew || B[i] >= nnew || A[i] == B[i] || f_has(cur, A[i], B[i])) { fprintf(stderr, "bad F edge (%d,%d)\n", A[i], B[i]); return 1; }
            f_add(cur, A[i], B[i]);
        }
        /* insert the remaining new vertices (F-degree 0) into random F edges */
        for (int i = 0; i < nnew; i++) {
            if (cur->fdeg[i] == tdeg[i]) continue;
            if (cur->fdeg[i] != 0) { fprintf(stderr, "initF: new vertex %d has F-degree %d, want 0 or %d\n", i, cur->fdeg[i], tdeg[i]); return 1; }
            int need = tdeg[i] / 2;  /* 1 edge split for a subdivision vertex, 2 for a hub */
            for (int tries = 0; need > 0; tries++) {
                if (tries > 10000 || cur->nF == 0) { fprintf(stderr, "initF: cannot insert vertex %d\n", i); return 1; }
                int t = rndint(cur->nF);
                int a = cur->fea[t], b = cur->feb[t];
                if (a == i || b == i || f_has(cur, a, i) || f_has(cur, i, b)) continue;
                f_del(cur, a, b); f_add(cur, a, i); f_add(cur, i, b);
                need--;
            }
            printf("F: inserted new vertex %d\n", i);
        }
        for (int i = 0; i < nnew; i++) if (cur->fdeg[i] != tdeg[i]) { fprintf(stderr, "F: new vertex %d has F-degree %d, want %d\n", i, cur->fdeg[i], tdeg[i]); return 1; }
    } else random_F(cur);
    full_rebuild(cur);
    printf("init cost %d (new-old %d, new-new %d, old-old %d) nBoo %d\n", cur->cost, cur->c_no, cur->c_nn, cur->c_oo, cur->nBoo);
    if (do_eval) {
        int mx, mn;
        int ex = exact_cost(cur, &mx, &mn);
        printf("exact all-pairs BFS: pairs > 6 = %d, max degree %d, min degree %d, n = %d\n", ex, mx, mn, N);
        print_state(cur, stdout);
        return 0;
    }
    if (mind > 0) {
        int bad = 0;
        for (int i = 0; i < ms; i++) for (int j = i + 1; j < ms; j++) bad += (edist[cur->S[i]][cur->S[j]] < mind);
        printf("initial S has %d pairs at edge-distance < %d\n", bad, mind);
    }
    memcpy(best, cur, sizeof(State));
    signal(SIGTERM, on_signal); signal(SIGINT, on_signal);
    if (exh2) {
        int found = exhaustive2(cur, trial, best, t_start);
        int mx, mn, ex = exact_cost(best, &mx, &mn);
        printf("exh2 done: found=%d evaluations=%ld best=%d exact=%d time=%.1f\n", found, exh_evals, found ? best->cost : exh_best, ex, now_sec() - t_start);
        if (found && ex == 0) {
            char path[1024];
            snprintf(path, sizeof path, "%s/%s_n%d_seed%llu_exh2.edges", outdir, tag, N, (unsigned long long)seed);
            write_edges(best, path);
            printf("ZERO COST: exact all-pairs check pairs>6 = %d, maxdeg %d, mindeg %d, n = %d\nwrote %s\n", ex, mx, mn, N, path);
        }
        print_state(best, stdout);
        write_best_file(best, outdir, tag, ex);
        return 0;
    }
    long it = 0, acc = 0;
    double t_last = now_sec();
    double cyc_len = seconds / cycles;
    int done = 0, last_cyc = 0;
    if (Tre < 0) Tre = T0;
    while (!done && !stop_flag) {
        double el_t = now_sec() - t_start;
        if (el_t >= seconds) break;
        int cyc = (int)(el_t / cyc_len);
        double frac = (el_t - cyc * cyc_len) / cyc_len;
        double Tstart = cyc == 0 ? T0 : Tre;
        double T = Tstart * pow(T1 / Tstart, frac);
        if (cyc != last_cyc) {
            last_cyc = cyc;
            if (restart_best) {
                memcpy(cur, best, sizeof(State));
                full_rebuild(cur);   /* resets the global damage counters to cur->S */
                printf("cycle %d: restart from best (%d)\n", cyc, best->cost);
            }
        }
        it++;
        double r = rnd01();
        memcpy(trial, cur, sizeof(State));
        int moved = 0, smove = 0, e_old = -1, e_new = -1;
        if (pRep > 0 && r < pRep && !fixS && !fixF) {
            moved = repair_move(cur, trial, &smove, &e_old, &e_new);
            if (!moved) smove = 0;
        } else if (r < pRep + pS * (1 - pRep) && !fixS) {
            int i = pick_new(trial, ms);
            int e;
            if (rnd01() < 0.5) {
                int old = trial->S[i];
                e = near2[offnear[old] + rndint(offnear[old + 1] - offnear[old])];
            } else e = rndint(NE);
            if (!trial->inS[e] && s_ok(trial, i, e)) {
                e_old = trial->S[i]; e_new = e;
                trial->inS[e_old] = 0; trial->S[i] = e; trial->inS[e] = 1;
                rebuild_smove(cur, trial, e_old, e_new);
                moved = 1; smove = 1;
            }
        } else if (!fixF) {
            moved = f_move(trial, pR / (1.0 - pS));
        }
        if (!moved) continue;
        int delta = trial->cost - cur->cost;
        int accept = (delta <= 0 || rnd01() < exp(-delta / T));
        if (smove) { if (accept) boo_commit(); else boo_undo(e_old, e_new); }
        if (accept) {
            State *tmp = cur; cur = trial; trial = tmp;
            acc++;
            if (cur->cost < best->cost) {
                memcpy(best, cur, sizeof(State));
                printf("t=%.1f it=%ld cyc=%d T=%.3f best=%d (no %d nn %d oo %d) nBoo=%d\n", now_sec() - t_start, it, cyc, T,
                       best->cost, best->c_no, best->c_nn, best->c_oo, best->nBoo);
                if (best->cost <= 6) print_state(best, stdout);
                fflush(stdout);
                if (best->cost <= 12) write_best_file(best, outdir, tag, -1);  /* survive a kill */
                if (best->cost == 0) {
                    int mx, mn;
                    int ex = exact_cost(best, &mx, &mn);
                    printf("ZERO COST: exact all-pairs check pairs>6 = %d, maxdeg %d, mindeg %d, n = %d\n", ex, mx, mn, N);
                    char path[1024];
                    snprintf(path, sizeof path, "%s/%s_n%d_seed%llu_it%ld.edges", outdir, tag, N, (unsigned long long)seed, it);
                    write_edges(best, path);
                    printf("wrote %s\n", path);
                    fflush(stdout);
                    if (!keepgoing) { done = 1; break; }
                }
            }
        }
        if (check_every > 0 && it % check_every == 0) {
            /* consistency check of the incremental boo counters: rebuild cur from scratch and compare */
            int cost_inc = cur->cost, nboo_inc = cur->nBoo;
            full_rebuild(cur);
            int mx, mn, ex = exact_cost(cur, &mx, &mn);
            if (cur->cost != cost_inc || cur->nBoo != nboo_inc || ex != cur->cost) {
                printf("CHECK FAILED at it=%ld: incremental cost %d nBoo %d, scratch cost %d nBoo %d, exact %d\n",
                       it, cost_inc, nboo_inc, cur->cost, cur->nBoo, ex);
                fflush(stdout);
                return 2;
            }
            nchecks++;
        }
        if (now_sec() - t_last > 30) {
            t_last = now_sec();
            printf("t=%.0f it=%ld acc=%ld cyc=%d T=%.3f cur=%d best=%d\n", t_last - t_start, it, acc, cyc, T, cur->cost, best->cost);
            fflush(stdout);
        }
    }
    int mx, mn;
    int ex = exact_cost(best, &mx, &mn);
    printf("FINAL best=%d (no %d nn %d oo %d) exact=%d iterations=%ld accepted=%ld time=%.1f n=%d maxdeg=%d mindeg=%d\n",
           best->cost, best->c_no, best->c_nn, best->c_oo, ex, it, acc, now_sec() - t_start, N, mx, mn);
    if (check_every > 0) printf("self-checks passed: %ld\n", nchecks);
    print_state(best, stdout);
    write_best_file(best, outdir, tag, ex);
    return 0;
}
