/* hexagon_sa.c -- simulated annealing for generalized-hexagon insertions at q = 4 ((5,6) cell); --freeze and
 * APARTMENT GADGETS (--apts) are optional structured flag sets:
 *   --apts K : S slots 4k..4k+3 (k < K) hold the flags at positions 0,3,6,9 of an apartment (12-cycle) of
 *              H(4) in the order a=pos 0, b=pos 3, c=pos 6, d=pos 9, and F contains the frozen chords
 *              (a,c), (b,d) (antipodal flags).  Such a group is self-repairing (0 far pairs on its own).
 *              Group moves replace the 4 flags of a group by another gadget (any apartment, any of the 3
 *              phases, any of the 8 orientations), chosen by the min-conflicts proxy over an enumeration of
 *              gadgets and evaluated with the exact cost; local moves try the 24 placements of the current
 *              apartment.  Slots >= 4K are ordinary free flags (single-flag relocation as in sa5).
 *   Gadget tables: opp[e][j] = the 4096 flags opposite to e; gtab/htab[e][j] = the two other gadget flags
 *   (middle edges of the two 5-geodesics closing the apartment through e and opp[e][j]).
 *
 * Simulated annealing for "pattern A" extensions of the generalized hexagon H(4) (q = 4):
 *   H = incidence graph of H(4): 2730 vertices, 5-regular, bipartite, girth 12, diameter 6.
 *   Subdivide a set S of m incidence edges (new vertex x_e adjacent to both ends of e),
 *   optionally add h "hub" vertices with no old neighbours,
 *   add a graph F on the new vertices with F-degree 3 at subdivision vertices and 5 at hubs,
 *   so the result G is 5-regular on 2730 + m + h vertices.
 * Cost = number of unordered vertex pairs at distance > 6 (exact).
 *
 * build:  gcc -O3 -march=native -o hexagon_sa hexagon_sa.c -lm
 * usage:  ./hexagon_sa H4_incidence.edges m h seed seconds [--T0 x] [--T1 y] [--cycles k] [--pS p] [--pR p]
 *            [--fixS file] [--fixF file] [--initS file] [--initF file] [--eval] [--tag name]
 *            [--outdir dir] [--keepgoing] [--restart-best] [--spread d] [--validate K] [--bench]
 * S file: lines "p l" (an edge of H, any order); F file: lines "i j" (indices of new vertices).
 *
 * Exact cost (G0 = H with S subdivided, no F edges; all H-geodesics are unique up to distance 5
 * and there are exactly 5 geodesics between vertices at distance 6; girth 12 makes every other
 * path at least 2 longer):
 *  - d0[i][v] = d_G0(x_i, v) capped at 7 = 1 + min over the ends u of e_i of (d_H(u,v) + #S-edges on
 *    the u-v geodesic), the correction being nonzero only on the "shadow" of each S edge seen
 *    from u (the subtree beyond that edge in the ball of radius 5 around u, which is a tree).
 *  - D[i][j] = d_G(x_i, x_j): Floyd-Warshall on the new vertices with weights min(1 if F-edge,
 *    d0[i][x_j]) (every x_i-x_j path decomposes at the new vertices it visits).
 *  - new-old: d_G(x_i, v) = min_j D[i][j] + d0[j][v], so v is fine iff it lies in the union over j of
 *    the G0-balls of radius 6 - D[i][j] around x_j (bitsets).
 *  - old-old: pairs with d_G0 > 6 ("boo" pairs) are exactly: distance-4 pairs with >= 3 S-edges on
 *    the geodesic, distance-5 pairs with >= 2, distance-6 pairs with an S-edge on each of the 5
 *    geodesics; maintained incrementally from per-edge geodesic lists.  For such a pair,
 *    d_G(a,b) = min_{i,j} d0[i][a] + D[i][j] + d0[j][b], and only i with d0[i][a] <= 4 matter.
 * Reference implementations (BFS from new vertices; all-pairs BFS) are used by --validate and
 * every cost-0 state is re-checked by all-pairs BFS before its edge list is written.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NOLD 2730
#define DEG 5
#define NE 6825
#define MAXNEW 128
#define MAXN (NOLD + MAXNEW)
#define NWORDS ((MAXN + 63) / 64)
#define OLDWORDS (NOLD / 64)          /* 42 full words of old vertices */
#define OLDTAILMASK ((1ULL << (NOLD % 64)) - 1)   /* 42 bits of word 42 */
#define FAR 7
#define UNREACH 255
#define NPAIR (NOLD * NOLD)
#define MAXBOO (1 << 20)
#define MAXF (MAXNEW * DEG / 2)
#define NGEO 5
#define C6MASK 0x1249   /* low bit of each of the five 3-bit fields */

static uint16_t adjH[NOLD][DEG];
static uint16_t edge_of[NOLD][DEG];
static uint16_t ep[NE], el[NE];
static uint8_t distH[NOLD][NOLD];
static int *P4, *P5, *P6;
static int off4[NE + 1], off5[NE + 1], off6[NE + 1];
static int *near2, offnear[NE + 1];
#define EDN 2729                      /* edges within edge-distance 4 of an edge (incl. itself) */
static uint16_t *edlist;              /* [NE][EDN], classes t = 0..4 consecutive */
static int ed_off[NE][6];
static uint8_t (*base)[NOLD];
static uint64_t (*basebits)[7][NWORDS];

/* boo counters: persistent, reflect the current S */
static uint8_t *c4, *c5;
static uint16_t *c6;
static int *pos, *boo, nBoo;

typedef struct {
    uint8_t d0[MAXNEW][MAXN];
    uint64_t bits[MAXNEW][7][NWORDS];
    uint8_t nearcnt[NOLD];
    uint8_t nearid[NOLD][MAXNEW];
    uint8_t neard[NOLD][MAXNEW];
} SDer;

typedef struct {
    uint8_t fadj[MAXNEW][DEG];
    uint8_t fdeg[MAXNEW];
    uint8_t fea[MAXF], feb[MAXF];
    int nF;
} Fgraph;

static int ms, hb, nnew, N, nn32;
static int S[MAXNEW];
static uint8_t inS[NE];
static int sidx[NE];   /* index in S of an S edge, -1 otherwise */
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
    static uint8_t degH[NOLD];
    memset(degH, 0, sizeof degH);
    char line[256];
    int nedges = 0;
    while (fgets(line, sizeof line, f)) {
        int u, w;
        if (line[0] == '#') continue;
        if (sscanf(line, "%d %d", &u, &w) != 2) continue;
        if (u < 0 || w < 0 || u >= NOLD || w >= NOLD || u == w) { fprintf(stderr, "bad edge %d %d\n", u, w); exit(1); }
        if (degH[u] >= DEG || degH[w] >= DEG) { fprintf(stderr, "degree > %d in H\n", DEG); exit(1); }
        adjH[u][degH[u]++] = (uint16_t)w;
        adjH[w][degH[w]++] = (uint16_t)u;
        nedges++;
    }
    fclose(f);
    if (nedges != NE) { fprintf(stderr, "expected %d edges, got %d\n", NE, nedges); exit(1); }
    for (int v = 0; v < NOLD; v++) {
        if (degH[v] != DEG) { fprintf(stderr, "vertex %d has degree %d\n", v, degH[v]); exit(1); }
        /* sort neighbours for determinism */
        for (int a = 0; a < DEG; a++)
            for (int b = a + 1; b < DEG; b++)
                if (adjH[v][b] < adjH[v][a]) { uint16_t t = adjH[v][a]; adjH[v][a] = adjH[v][b]; adjH[v][b] = t; }
    }
    int nE = 0;
    for (int u = 0; u < NOLD; u++)
        for (int k = 0; k < DEG; k++) {
            int w = adjH[u][k];
            if (u < w) {
                ep[nE] = (uint16_t)u; el[nE] = (uint16_t)w;
                edge_of[u][k] = (uint16_t)nE;
                for (int kk = 0; kk < DEG; kk++) if (adjH[w][kk] == u) edge_of[w][kk] = (uint16_t)nE;
                nE++;
            }
        }
    if (nE != NE) { fprintf(stderr, "edge count mismatch\n"); exit(1); }
}

static void bfs_H(int s, uint8_t *dist) {
    static int q[NOLD];
    memset(dist, 255, NOLD);
    int qh = 0, qt = 0;
    dist[s] = 0; q[qt++] = s;
    while (qh < qt) {
        int u = q[qh++];
        for (int k = 0; k < DEG; k++) {
            int w = adjH[u][k];
            if (dist[w] == 255) { dist[w] = (uint8_t)(dist[u] + 1); q[qt++] = w; }
        }
    }
}

static int walk_edges(int from, int to, const uint8_t *dto, int *out) {
    /* unique geodesic from `from` to `to` (d <= 5); dto = distH[to]; returns the edges used */
    int cnt = 0, v = from;
    while (v != to) {
        int found = -1;
        for (int k = 0; k < DEG; k++) {
            int w = adjH[v][k];
            if (dto[w] == dto[v] - 1) { if (found >= 0) { fprintf(stderr, "non-unique geodesic\n"); exit(1); } found = k; }
        }
        if (found < 0) { fprintf(stderr, "walk failed\n"); exit(1); }
        out[cnt++] = edge_of[v][found];
        v = adjH[v][found];
    }
    return cnt;
}

static int find_edge(int u, int v) {
    for (int k = 0; k < DEG; k++) if (adjH[u][k] == v) return edge_of[u][k];
    return -1;
}

static void precompute(void) {
    for (int a = 0; a < NOLD; a++) bfs_H(a, distH[a]);
    for (int a = 0; a < NOLD; a++) for (int b = 0; b < NOLD; b++)
        if (distH[a][b] > 6) { fprintf(stderr, "diameter > 6 in H\n"); exit(1); }
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
            P4 = malloc(sizeof(int) * (size_t)off4[NE]);
            P5 = malloc(sizeof(int) * (size_t)off5[NE]);
            P6 = malloc(sizeof(int) * (size_t)off6[NE]);
            if (!P4 || !P5 || !P6) { fprintf(stderr, "alloc failed\n"); exit(1); }
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
                } else {
                    for (int g = 0; g < DEG; g++) {
                        int w0 = adjH[a][g];
                        if (distH[w0][b] != 5) { fprintf(stderr, "distance-6 structure violated\n"); exit(1); }
                        int c = walk_edges(b, w0, distH[w0], eds);
                        eds[c++] = edge_of[a][g];
                        for (int t = 0; t < c; t++) {
                            int e = eds[t];
                            if (pass) P6[off6[e] + cnt6[e]] = (pid << 3) | g;
                            cnt6[e]++;
                        }
                    }
                }
            }
        }
    }
    fprintf(stderr, "precomp: P4 %d P5 %d P6 %d entries; per edge (edge 0): %d %d %d\n",
            off4[NE], off5[NE], off6[NE], cnt4[0], cnt5[0], cnt6[0]);
    free(cnt4); free(cnt5); free(cnt6);
    /* base distances and ball bitsets per edge */
    base = malloc((size_t)NE * NOLD);
    basebits = calloc((size_t)NE, sizeof *basebits);
    if (!base || !basebits) { fprintf(stderr, "alloc failed\n"); exit(1); }
    for (int e = 0; e < NE; e++) {
        int p = ep[e], l = el[e];
        for (int v = 0; v < NOLD; v++) {
            int d = distH[p][v] < distH[l][v] ? distH[p][v] : distH[l][v];
            d += 1;
            if (d > 6) { fprintf(stderr, "base > 6?\n"); exit(1); }
            base[e][v] = (uint8_t)d;
            for (int r = d; r <= 6; r++) basebits[e][r][v >> 6] |= 1ULL << (v & 63);
        }
    }
    /* edge-distance classes: for each edge e, the edges f at edge-distance t = 0..4 (min endpoint
       distance), as one CSR with per-class offsets ed_off[e][t]..ed_off[e][t+1]; the classes have the
       same size for every edge (9 incl. e, 32, 128, 512, 2048 for t = 0..4; 2729 within distance 4).
       near2 (edges within distance 2, incl. e) is the prefix of classes 0..2 (169 edges). */
    {
        static int cls_size[5] = { 9, 32, 128, 512, 2048 };   /* class 0 includes e itself */
        edlist = malloc(sizeof(uint16_t) * (size_t)NE * EDN);
        if (!edlist) { fprintf(stderr, "alloc failed\n"); exit(1); }
        int fill[5];
        for (int e = 0; e < NE; e++) {
            int p = ep[e], l = el[e];
            ed_off[e][0] = 0;
            for (int t = 0; t < 5; t++) { ed_off[e][t + 1] = ed_off[e][t] + cls_size[t]; fill[t] = 0; }
            if (ed_off[e][5] != EDN) { fprintf(stderr, "EDN mismatch\n"); exit(1); }
            for (int f = 0; f < NE; f++) {
                int d = distH[p][ep[f]], t;
                t = distH[p][el[f]]; if (t < d) d = t;
                t = distH[l][ep[f]]; if (t < d) d = t;
                t = distH[l][el[f]]; if (t < d) d = t;
                if (d <= 4) {
                    if (fill[d] >= cls_size[d]) { fprintf(stderr, "edge-distance class %d overflow at edge %d\n", d, e); exit(1); }
                    edlist[(size_t)e * EDN + ed_off[e][d] + fill[d]++] = (uint16_t)f;
                }
            }
            for (int t = 0; t < 5; t++) if (fill[t] != cls_size[t]) { fprintf(stderr, "edge-distance class %d has %d edges at edge %d\n", t, fill[t], e); exit(1); }
        }
    }
    near2 = malloc(sizeof(int) * (size_t)NE * 169);
    for (int e = 0; e < NE; e++) {
        offnear[e] = e * 169;
        for (int k = 0; k < 169; k++) near2[e * 169 + k] = edlist[(size_t)e * EDN + k];
    }
    offnear[NE] = NE * 169;
    fprintf(stderr, "edge-distance classes built (near2: %d entries per edge)\n", offnear[1] - offnear[0]);
    /* boo counters */
    c4 = calloc(NPAIR, 1); c5 = calloc(NPAIR, 1); c6 = calloc(NPAIR, sizeof(uint16_t));
    pos = malloc(sizeof(int) * (size_t)NPAIR); boo = malloc(sizeof(int) * MAXBOO);
    if (!c4 || !c5 || !c6 || !pos || !boo) { fprintf(stderr, "alloc failed\n"); exit(1); }
    for (int i = 0; i < NPAIR; i++) pos[i] = -1;
    nBoo = 0;
    for (int e = 0; e < NE; e++) sidx[e] = -1;
}

/* ---------- boo pairs (old-old pairs with d_G0 > 6), incremental ---------- */
static inline int is_boo(int pid) {
    uint16_t c = c6[pid];
    return c4[pid] >= 3 || c5[pid] >= 2 || (((c | (c >> 1) | (c >> 2)) & C6MASK) == C6MASK);
}
static inline void boo_update(int pid, int before) {
    int after = is_boo(pid);
    if (after && !before) {
        if (nBoo >= MAXBOO) { fprintf(stderr, "MAXBOO exceeded\n"); exit(1); }
        pos[pid] = nBoo; boo[nBoo++] = pid;
    } else if (!after && before) {
        int k = pos[pid];
        boo[k] = boo[--nBoo]; pos[boo[k]] = k; pos[pid] = -1;
    }
}
static void boo_delta(int e, int sign) {
    for (int t = off4[e]; t < off4[e + 1]; t++) {
        int pid = P4[t]; int b = is_boo(pid);
        c4[pid] = (uint8_t)(c4[pid] + sign); boo_update(pid, b);
    }
    for (int t = off5[e]; t < off5[e + 1]; t++) {
        int pid = P5[t]; int b = is_boo(pid);
        c5[pid] = (uint8_t)(c5[pid] + sign); boo_update(pid, b);
    }
    for (int t = off6[e]; t < off6[e + 1]; t++) {
        int pid = P6[t] >> 3, g = P6[t] & 7; int b = is_boo(pid);
        c6[pid] = (uint16_t)(c6[pid] + sign * (1 << (3 * g))); boo_update(pid, b);
    }
}
static void set_S_edge(int i, int e) {
    /* replace S[i] by e (e not in S) and update the boo structure */
    int old = S[i];
    boo_delta(old, -1); boo_delta(e, +1);
    inS[old] = 0; inS[e] = 1; S[i] = e;
    sidx[old] = -1; sidx[e] = i;
}
static void install_S(const int *newS) {
    /* replace the whole S (boo structure brought in sync) */
    for (int i = 0; i < ms; i++) if (S[i] >= 0 && inS[S[i]]) { boo_delta(S[i], -1); inS[S[i]] = 0; sidx[S[i]] = -1; }
    for (int i = 0; i < ms; i++) {
        S[i] = newS[i];
        if (inS[S[i]]) { fprintf(stderr, "duplicate S edge\n"); exit(1); }
        inS[S[i]] = 1; sidx[S[i]] = i; boo_delta(S[i], +1);
    }
}

/* ---------- apartment groups: slots 4k..4k+3 ---------- */
static int napt = 0;                       /* number of apartment groups (slots 0..4*napt-1) */
static inline int in_group(int i, int k) { return i >= 4 * k && i < 4 * k + 4; }
static void set_group(int k, const int *x) {
    /* replace the four flags of group k by x[0..3] (slots 4k..4k+3); x may reuse the group's own flags
       (re-orientation / phase change), so remove all four first, then add all four */
    for (int t = 0; t < 4; t++) {
        int old = S[4 * k + t];
        boo_delta(old, -1); inS[old] = 0; sidx[old] = -1; S[4 * k + t] = -1;
    }
    for (int t = 0; t < 4; t++) {
        int e = x[t];
        if (e < 0 || e >= NE || inS[e]) { fprintf(stderr, "set_group: bad/duplicate flag %d in group %d\n", e, k); exit(1); }
        boo_delta(e, +1); inS[e] = 1; sidx[e] = 4 * k + t; S[4 * k + t] = e;
    }
}
/* undo buffer of the last applied S move (single flag: su_n = 1; group: su_n = 4) */
static int su_n = 0, su_i[4], su_e[4];
static void s_undo(void) {
    if (su_n == 1) set_S_edge(su_i[0], su_e[0]);
    else if (su_n == 4) set_group(su_i[0] / 4, su_e);
    su_n = 0;
}

/* ---------- S-derived data: d0 rows, ball bitsets, near lists ---------- */
static int16_t corrbuf[2][NOLD];
static int touched[2][NOLD], ntouched[2];

static void shadow(int u, int a, int b, int16_t *corr, int *touch, int *nt) {
    /* add 1 to corr[v] for every v in the subtree beyond b (away from a) within distance 5 of u */
    static int stack[4096];
    int sp = 0;
    stack[sp++] = b;
    while (sp > 0) {
        int w = stack[--sp];
        if (corr[w] == 0) touch[(*nt)++] = w;
        corr[w]++;
        int dw = distH[u][w];
        if (dw < 5)
            for (int k = 0; k < DEG; k++) {
                int x = adjH[w][k];
                if (distH[u][x] == dw + 1) stack[sp++] = x;
            }
    }
    (void)a;
}

static void compute_row(SDer *sd, int i) {
    uint8_t *row = sd->d0[i];
    uint64_t (*bits)[NWORDS] = sd->bits[i];
    if (i >= ms) {   /* hub */
        memset(row, FAR, N);
        memset(bits, 0, sizeof(uint64_t) * 7 * NWORDS);
        row[NOLD + i] = 0;
        for (int r = 0; r <= 6; r++) bits[r][(NOLD + i) >> 6] |= 1ULL << ((NOLD + i) & 63);
        return;
    }
    int e = S[i];
    int ends[2] = { ep[e], el[e] };
    memcpy(row, base[e], NOLD);
    memcpy(bits, basebits[e], sizeof(uint64_t) * 7 * NWORDS);
    for (int s = 0; s < 2; s++) {
        int u = ends[s];
        ntouched[s] = 0;
        for (int j = 0; j < ms; j++) {
            if (j == i) continue;
            int f = S[j];
            int a = ep[f], b = el[f];
            int ta = distH[u][a], tb = distH[u][b];
            if (ta > tb) { int t = a; a = b; b = t; t = ta; ta = tb; tb = t; }
            if (tb > 5) continue;
            shadow(u, a, b, corrbuf[s], touched[s], &ntouched[s]);
        }
    }
    for (int s = 0; s < 2; s++) {
        for (int t = 0; t < ntouched[s]; t++) {
            int v = touched[s][t];
            int dp = distH[ends[0]][v] + corrbuf[0][v];
            int dl = distH[ends[1]][v] + corrbuf[1][v];
            int d = 1 + (dp < dl ? dp : dl);
            if (d > 6) d = FAR;
            int old = row[v];
            if (d > old) {
                for (int r = old; r < d && r <= 6; r++) bits[r][v >> 6] &= ~(1ULL << (v & 63));
                row[v] = (uint8_t)d;
            }
        }
    }
    for (int s = 0; s < 2; s++) {
        for (int t = 0; t < ntouched[s]; t++) corrbuf[s][touched[s][t]] = 0;
    }
    /* new-vertex columns */
    for (int j = 0; j < nnew; j++) {
        int x = NOLD + j;
        int d;
        if (j == i) d = 0;
        else if (j < ms) {
            int f = S[j];
            int a = row[ep[f]], b = row[el[f]];
            d = 1 + (a < b ? a : b);
            if (d > 6) d = FAR;
        } else d = FAR;
        row[x] = (uint8_t)d;
        if (d <= 6) for (int r = d; r <= 6; r++) bits[r][x >> 6] |= 1ULL << (x & 63);
    }
    for (int x = NOLD + nnew; x < MAXN; x++) row[x] = FAR;
}

static void build_near(SDer *sd) {
    memset(sd->nearcnt, 0, NOLD);
    for (int i = 0; i < nnew; i++) {
        const uint64_t *b4 = sd->bits[i][4];
        for (int w = 0; w <= OLDWORDS; w++) {
            uint64_t word = b4[w];
            if (w == OLDWORDS) word &= OLDTAILMASK;
            while (word) {
                int bit = __builtin_ctzll(word);
                word &= word - 1;
                int v = w * 64 + bit;
                int c = sd->nearcnt[v];
                sd->nearid[v][c] = (uint8_t)i;
                sd->neard[v][c] = sd->d0[i][v];
                sd->nearcnt[v] = (uint8_t)(c + 1);
            }
        }
    }
}

static void build_sder(SDer *sd) {
    for (int i = 0; i < nnew; i++) compute_row(sd, i);
    build_near(sd);
}

/* ---------- F helpers ---------- */
static int f_has(const Fgraph *F, int a, int b) {
    for (int k = 0; k < F->fdeg[a]; k++) if (F->fadj[a][k] == b) return 1;
    return 0;
}
static void f_add(Fgraph *F, int a, int b) {
    F->fadj[a][F->fdeg[a]++] = (uint8_t)b;
    F->fadj[b][F->fdeg[b]++] = (uint8_t)a;
    F->fea[F->nF] = (uint8_t)(a < b ? a : b); F->feb[F->nF] = (uint8_t)(a < b ? b : a); F->nF++;
}
static void f_del(Fgraph *F, int a, int b) {
    for (int k = 0; k < F->fdeg[a]; k++) if (F->fadj[a][k] == b) { F->fadj[a][k] = F->fadj[a][--F->fdeg[a]]; break; }
    for (int k = 0; k < F->fdeg[b]; k++) if (F->fadj[b][k] == a) { F->fadj[b][k] = F->fadj[b][--F->fdeg[b]]; break; }
    int lo = a < b ? a : b, hi = a < b ? b : a;
    for (int t = 0; t < F->nF; t++)
        if (F->fea[t] == lo && F->feb[t] == hi) { F->fea[t] = F->fea[F->nF - 1]; F->feb[t] = F->feb[F->nF - 1]; F->nF--; break; }
}

static uint8_t frozen[MAXNEW][MAXNEW];   /* F edges that must not be removed (--freeze file / chords) */
static void add_chords(Fgraph *F) {
    /* the frozen chords of the apartment groups: (4k, 4k+2), (4k+1, 4k+3) */
    for (int k = 0; k < napt; k++) {
        f_add(F, 4 * k, 4 * k + 2); f_add(F, 4 * k + 1, 4 * k + 3);
        frozen[4 * k][4 * k + 2] = frozen[4 * k + 2][4 * k] = 1;
        frozen[4 * k + 1][4 * k + 3] = frozen[4 * k + 3][4 * k + 1] = 1;
    }
}
static void random_F(Fgraph *F) {
    /* random F of the prescribed degrees; when apartment groups exist it contains their chords */
    for (;;) {
        for (int i = 0; i < nnew; i++) F->fdeg[i] = 0;
        F->nF = 0;
        add_chords(F);
        int stubs[DEG * MAXNEW], ns = 0;
        for (int i = 0; i < nnew; i++) for (int k = F->fdeg[i]; k < tdeg[i]; k++) stubs[ns++] = i;
        int ok = 1;
        while (ns > 0 && ok) {
            int tries = 0;
            for (;;) {
                int ia = rndint(ns), ib = rndint(ns);
                int a = stubs[ia], b = stubs[ib];
                if (a != b && !f_has(F, a, b)) {
                    f_add(F, a, b);
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

/* undo log for F moves */
static int ulog_a[8], ulog_b[8], ulog_op[8], ulog_n;   /* op 1 = added, 0 = deleted */
static void lf_add(Fgraph *F, int a, int b) { f_add(F, a, b); ulog_a[ulog_n] = a; ulog_b[ulog_n] = b; ulog_op[ulog_n++] = 1; }
static void lf_del(Fgraph *F, int a, int b) { f_del(F, a, b); ulog_a[ulog_n] = a; ulog_b[ulog_n] = b; ulog_op[ulog_n++] = 0; }
static void f_undo(Fgraph *F) {
    for (int t = ulog_n - 1; t >= 0; t--) {
        if (ulog_op[t]) f_del(F, ulog_a[t], ulog_b[t]); else f_add(F, ulog_a[t], ulog_b[t]);
    }
    ulog_n = 0;
}

/* random F move (relocate with probability prel, else 2-swap); returns 1 if applied */
static int f_move(Fgraph *F, double prel) {
    ulog_n = 0;
    if (rnd01() < prel) {
        /* relocate two of i's edges: (i,a),(i,b),(c,d) -> (a,b),(c,i),(i,d) */
        int i = rndint(nnew);
        if (F->fdeg[i] < 2) return 0;
        int ka = rndint(F->fdeg[i]), kb = rndint(F->fdeg[i] - 1);
        if (kb >= ka) kb++;
        int a = F->fadj[i][ka], b = F->fadj[i][kb];
        if (f_has(F, a, b)) return 0;
        int k = rndint(F->nF);
        int c = F->fea[k], d = F->feb[k];
        if (c == i || d == i || f_has(F, c, i) || f_has(F, d, i)) return 0;
        if (frozen[i][a] || frozen[i][b] || frozen[c][d]) return 0;
        lf_del(F, i, a); lf_del(F, i, b); lf_add(F, a, b);
        lf_del(F, c, d); lf_add(F, c, i); lf_add(F, i, d);
        return 1;
    }
    int t1 = rndint(F->nF), t2 = rndint(F->nF);
    if (t1 == t2) return 0;
    int a = F->fea[t1], b = F->feb[t1], c = F->fea[t2], d = F->feb[t2];
    if (rnd01() < 0.5) { int x = c; c = d; d = x; }
    if (a == c || a == d || b == c || b == d || f_has(F, a, c) || f_has(F, b, d)) return 0;
    if (frozen[a][b] || frozen[c][d]) return 0;
    lf_del(F, a, b); lf_del(F, c, d); lf_add(F, a, c); lf_add(F, b, d);
    return 1;
}

/* ---------- fast exact evaluation ---------- */
/* result of an evaluation, including the violations (for targeted moves) */
typedef struct {
    int no, nn, oo;
    uint8_t D[MAXNEW][MAXNEW];                 /* new-new distances (7 = far) */
    uint64_t uncov[MAXNEW][OLDWORDS + 1];      /* old vertices at distance > 6 from x_i */
    int nocnt[MAXNEW];
    int viol_nn_a[MAXNEW * MAXNEW / 2], viol_nn_b[MAXNEW * MAXNEW / 2];
    int *viol_oo;                              /* unrepaired boo pairs (pids), capacity MAXBOO */
} Eval;
static Eval evA, evB;
static Eval *ev_cur = &evA, *ev_trial = &evB;

static int evaluate(const SDer *sd, const Fgraph *F, Eval *ev) {
    uint8_t (*D)[MAXNEW] = ev->D;
    /* D: new-new distances */
    for (int i = 0; i < nnew; i++) {
        memcpy(D[i], sd->d0[i] + NOLD, nnew);
        for (int j = nnew; j < nn32; j++) D[i][j] = FAR;
        for (int k = 0; k < F->fdeg[i]; k++) D[i][F->fadj[i][k]] = 1;
    }
    for (int k = 0; k < nnew; k++) {
        const uint8_t *Dk = D[k];
        for (int i = 0; i < nnew; i++) {
            int dik = D[i][k];
            if (dik >= 6) continue;
            uint8_t *Di = D[i];
            uint8_t dd = (uint8_t)dik;
            for (int j = 0; j < nn32; j++) {
                uint8_t v = (uint8_t)(dd + Dk[j]);
                Di[j] = Di[j] < v ? Di[j] : v;
            }
        }
    }
    int cnn = 0;
    for (int i = 0; i < nnew; i++)
        for (int j = i + 1; j < nnew; j++)
            if (D[i][j] > 6) { ev->viol_nn_a[cnn] = i; ev->viol_nn_b[cnn] = j; cnn++; }
    /* new-old via union of balls */
    int cno = 0;
    for (int i = 0; i < nnew; i++) {
        uint64_t *acc = ev->uncov[i];
        memcpy(acc, sd->bits[i][6], sizeof(uint64_t) * (OLDWORDS + 1));
        for (int j = 0; j < nnew; j++) {
            if (j == i) continue;
            int dij = D[i][j];
            if (dij > 5) continue;
            const uint64_t *b = sd->bits[j][6 - dij];
            for (int w = 0; w <= OLDWORDS; w++) acc[w] |= b[w];
        }
        int covered = 0;
        for (int w = 0; w < OLDWORDS; w++) { acc[w] = ~acc[w]; covered += __builtin_popcountll(acc[w]); }
        acc[OLDWORDS] = ~acc[OLDWORDS] & OLDTAILMASK;
        covered += __builtin_popcountll(acc[OLDWORDS]);
        ev->nocnt[i] = covered;   /* acc now holds the uncovered old vertices */
        cno += covered;
    }
    /* old-old boo pairs */
    int coo = 0;
    for (int t = 0; t < nBoo; t++) {
        int pid = boo[t];
        int a = pid / NOLD, b = pid % NOLD;
        int na = sd->nearcnt[a], nb = sd->nearcnt[b];
        int ok = 0;
        for (int ia = 0; ia < na && !ok; ia++) {
            int i = sd->nearid[a][ia], da = sd->neard[a][ia];
            const uint8_t *Di = D[i];
            for (int jb = 0; jb < nb; jb++) {
                if (da + Di[sd->nearid[b][jb]] + sd->neard[b][jb] <= 6) { ok = 1; break; }
            }
        }
        if (!ok) ev->viol_oo[coo++] = pid;
    }
    ev->no = cno; ev->nn = cnn; ev->oo = coo;
    return cno + cnn + coo;
}

/* ---------- targeted moves: pick a violation of the current state and propose a repair ---------- */

/* 2-swap creating the F edge (i,k): (i,a),(k,b) -> (i,k),(a,b).  All valid (a,b) are tried and the
   one with the lowest cost is kept (evaluated into ev_trial, cost in tm_cost); returns 1 if applied */
static int tm_cost;
static int kF = 9, kS = 4;   /* candidates examined per targeted move */
static const SDer *tm_sd;
static int f_create_edge(Fgraph *F, int i, int k) {
    ulog_n = 0;
    if (i == k || f_has(F, i, k) || F->fdeg[i] == 0 || F->fdeg[k] == 0) return 0;
    int ca[DEG * DEG], cb[DEG * DEG], nc = 0;
    for (int x = 0; x < F->fdeg[i]; x++)
        for (int y = 0; y < F->fdeg[k]; y++) {
            int a = F->fadj[i][x], b = F->fadj[k][y];
            if (a == k || b == i || a == b || f_has(F, a, b)) continue;
            if (frozen[i][a] || frozen[k][b]) continue;
            ca[nc] = a; cb[nc] = b; nc++;
        }
    if (nc == 0) return 0;
    int besta = -1, bestb = -1, bestc = 1 << 30;
    int tries = nc < kF ? nc : kF;
    for (int t = 0; t < tries; t++) {
        int idx = (tries == nc) ? t : rndint(nc);
        int a = ca[idx], b = cb[idx];
        ulog_n = 0;
        lf_del(F, i, a); lf_del(F, k, b); lf_add(F, i, k); lf_add(F, a, b);
        int c = evaluate(tm_sd, F, ev_trial);
        f_undo(F);
        if (c < bestc) { bestc = c; besta = a; bestb = b; }
    }
    ulog_n = 0;
    lf_del(F, i, besta); lf_del(F, k, bestb); lf_add(F, i, k); lf_add(F, besta, bestb);
    tm_cost = evaluate(tm_sd, F, ev_trial);
    return 1;
}

/* move S edge e_old to the best of kS candidate positions (near2 / uniform); the move is applied,
   sd_trial built and evaluated into ev_trial (cost in tm_cost); returns 1 if applied */
static SDer *tm_sd_trial;
static int random_target_edge(int e_old);
static int mc_mode = 0;              /* S moves use the min-conflicts relocation */
static long tabu_until[MAXNEW], tabu_len = 0, n_smoves = 0;
static int s_relocate_mc(const Fgraph *F, int i, int *s_i, int *s_old);
static int group_relocate_mc(const Fgraph *F, int k, int rslot, int rv, int rR);
static int group_tabu(int k) {
    if (tabu_len <= 0) return 0;
    for (int t = 0; t < 4; t++) if (tabu_until[4 * k + t] > n_smoves) return 1;
    return 0;
}
static int s_move_best(const Fgraph *F, int e_old, int *s_i, int *s_old) {
    int i = sidx[e_old];
    if (i < 0) return 0;
    if (i < 4 * napt) {   /* apartment member: the whole group moves */
        if (group_tabu(i / 4)) return 0;
        return group_relocate_mc(F, i / 4, -1, 0, 0);
    }
    if (tabu_len > 0 && tabu_until[i] > n_smoves) return 0;
    if (mc_mode) return s_relocate_mc(F, i, s_i, s_old);
    int beste = -1, bestc = 1 << 30;
    for (int t = 0; t < kS; t++) {
        int e = random_target_edge(e_old);
        if (inS[e]) continue;
        set_S_edge(i, e);
        build_sder(tm_sd_trial);
        int c = evaluate(tm_sd_trial, F, ev_trial);
        set_S_edge(i, e_old);
        if (c < bestc) { bestc = c; beste = e; }
    }
    if (beste < 0) return 0;
    set_S_edge(i, beste);
    build_sder(tm_sd_trial);
    tm_cost = evaluate(tm_sd_trial, F, ev_trial);
    *s_i = i; *s_old = e_old;
    su_n = 1; su_i[0] = i; su_e[0] = e_old;
    return 2;
}

/* ---------- F-only distances, F pre-optimisation, min-conflicts relocation ---------- */
static uint8_t DF[MAXNEW][MAXNEW];
static void compute_DF(const Fgraph *F) {
    static int q[MAXNEW];
    for (int s = 0; s < nnew; s++) {
        uint8_t *d = DF[s];
        memset(d, 99, MAXNEW);
        int qh = 0, qt = 0;
        d[s] = 0; q[qt++] = s;
        while (qh < qt) {
            int u = q[qh++];
            for (int k = 0; k < F->fdeg[u]; k++) {
                int w = F->fadj[u][k];
                if (d[w] == 99) { d[w] = (uint8_t)(d[u] + 1); q[qt++] = w; }
            }
        }
    }
}
/* F quality: pairs at F-distance >= 7 count 100000, at distance 6 count 100, at distance 5 count w5
   (fewer far pairs = fewer constraints on the embedding) */
static int w5 = 0;   /* optional weight of distance-5 pairs in the F objective (smoother landscape) */
static int f_quality(const Fgraph *F) {
    compute_DF(F);
    int c = 0;
    for (int i = 0; i < nnew; i++)
        for (int j = i + 1; j < nnew; j++) {
            int d = DF[i][j];
            c += d >= 7 ? 100000 : d == 6 ? 100 : d == 5 ? w5 : 0;
        }
    return c;
}
static double optF_per = 5.0;   /* length of one restart */
static int optF_keep = 0;       /* restarts start from the initial F instead of a random cubic graph */
static void optimise_F(Fgraph *F, double seconds) {
    /* restarts of optF_per s (random cubic graph + annealing on 2-swaps/relocations); keeps the best F */
    double t_all = now_sec();
    int nrest = seconds >= 2 * optF_per ? (int)(seconds / optF_per) : 1;
    double per = seconds / nrest;
    Fgraph bestF = *F, initF = *F;
    int bestq = f_quality(F);
    long it = 0;
    for (int r = 0; r < nrest; r++) {
        if (r > 0) { if (optF_keep) *F = initF; else random_F(F); }
        double t0 = now_sec();
        int cur = f_quality(F);
        while (now_sec() - t0 < per) {
            it++;
            double frac = (now_sec() - t0) / per;
            double T = 300.0 * pow(0.05, frac);   /* in units where one distance-6 pair costs 100 */
            if (!f_move(F, 0.2)) continue;
            int c = f_quality(F);
            int delta = c - cur;
            if (delta <= 0 || rnd01() < exp(-delta / T)) cur = c;
            else f_undo(F);
        }
        if (cur < bestq) { bestq = cur; bestF = *F; }
    }
    *F = bestF;
    compute_DF(F);
    int n6 = 0, n7 = 0, diam = 0;
    for (int i = 0; i < nnew; i++)
        for (int j = i + 1; j < nnew; j++) { n6 += DF[i][j] == 6; n7 += DF[i][j] >= 7; if (DF[i][j] > diam) diam = DF[i][j]; }
    printf("optimise_F: %d restarts, %ld moves, %.1f s, quality %d, F-diameter %d, pairs at F-distance 6: %d, >= 7: %d\n",
           nrest, it, now_sec() - t_all, bestq, diam, n6, n7);
}

/* proxy cost of placing new vertex i on edge e given the other S edges and the F-only distances:
   sum over k != i of 4^(4 - t) when t = edge-distance(e, e_k) <= 4 and DF[i][k] > t + 1 (shadows that
   F cannot repair directly), plus 1 when t >= 5 and DF[i][k] >= 7 (new-new pair out of reach). */
static int proxy[NE];
static int proxy_excl_lo = -1, proxy_excl_hi = -1;   /* slots skipped by the proxy (the moving group) */
static void relocation_proxy(int i) {
    int base_pen = 0;
    for (int k = 0; k < nnew; k++) if (k != i && !(k >= proxy_excl_lo && k <= proxy_excl_hi) && DF[i][k] >= 7) base_pen++;
    for (int e = 0; e < NE; e++) proxy[e] = base_pen;
    static const int shadow_w[5] = { 256, 64, 16, 4, 1 };
    for (int k = 0; k < nnew; k++) {
        if (k == i || k >= ms) continue;
        if (k >= proxy_excl_lo && k <= proxy_excl_hi) continue;
        int f = S[k];
        int dfk = DF[i][k];
        const uint16_t *L = edlist + (size_t)f * EDN;
        for (int t = 0; t <= 4; t++) {
            int w = (dfk > t + 1) ? shadow_w[t] : 0;
            if (dfk >= 7) w -= 1;   /* within distance 4: the new-new pair is reachable in G0 */
            if (w == 0) continue;
            for (int x = ed_off[f][t]; x < ed_off[f][t + 1]; x++) proxy[L[x]] += w;
        }
    }
}

/* relocate new vertex i (< ms) to the best of the kS lowest-proxy edges (ties broken randomly); the
   move is applied, sd_trial built and evaluated into ev_trial (cost in tm_cost); returns 2, or 0 */
static int reloc_v = -1, reloc_R = 0;   /* when reloc_v >= 0: only edges with an end within reloc_R of reloc_v */
static int s_relocate_mc(const Fgraph *F, int i, int *s_i, int *s_old) {
    if (i < 0 || i >= ms) return 0;
    relocation_proxy(i);
    /* select kS candidates: the kS smallest proxy values, random among ties, excluding S edges */
    int cand[64], ncand = 0;
    int K = kS < 64 ? kS : 64;
    uint32_t keys[64];
    for (int e = 0; e < NE; e++) {
        if (inS[e]) continue;
        if (reloc_v >= 0) {
            int d = distH[ep[e]][reloc_v] < distH[el[e]][reloc_v] ? distH[ep[e]][reloc_v] : distH[el[e]][reloc_v];
            if (d > reloc_R) continue;
        }
        uint32_t key = ((uint32_t)proxy[e] << 12) | (uint32_t)(rnd64() & 4095);
        if (ncand < K) { cand[ncand] = e; keys[ncand] = key; ncand++; }
        else {
            int worst = 0;
            for (int c = 1; c < ncand; c++) if (keys[c] > keys[worst]) worst = c;
            if (key < keys[worst]) { cand[worst] = e; keys[worst] = key; }
        }
    }
    if (ncand == 0) return 0;
    int e_old = S[i];
    int beste = -1, bestc = 1 << 30;
    for (int c = 0; c < ncand; c++) {
        set_S_edge(i, cand[c]);
        build_sder(tm_sd_trial);
        int cost = evaluate(tm_sd_trial, F, ev_trial);
        set_S_edge(i, e_old);
        if (cost < bestc) { bestc = cost; beste = cand[c]; }
    }
    set_S_edge(i, beste);
    build_sder(tm_sd_trial);
    tm_cost = evaluate(tm_sd_trial, F, ev_trial);
    *s_i = i; *s_old = e_old;
    su_n = 1; su_i[0] = i; su_e[0] = e_old;
    return 2;
}

/* ---------- apartment gadgets ---------- */
/* An apartment is a 12-cycle v0..v11 of H (girth 12 = twice the diameter, so it is isometric).  The gadget
   with flags at positions 0,3,6,9 is determined by any one of its antipodal flag pairs: for opposite flags
   e = (p,L), f = (p',L') (d(p,p') = d(L,L') = 6, d(L,p') = d(p,L') = 5) the cycle is p, L, ..., p', L', ..., p
   with the unique 5-geodesics L -> p' and L' -> p, and the flags at positions 3 and 9 are the middle edges
   of those geodesics.  Slots: a = pos 0, b = pos 3, c = pos 6, d = pos 9; chords (a,c), (b,d). */
#define NOPP 4096
static int edge_dist(int e, int f);
static uint16_t *mid5;      /* [NOLD*NOLD] middle edge of the unique geodesic between vertices at distance 5 */
static uint16_t *opp;       /* [NE][NOPP] flags opposite to e */
static uint16_t *gtab, *htab;   /* [NE][NOPP]: gadget flags at positions 3 and 9 for e at 0, opp[e][j] at 6 */
static int topT = 16;       /* edges per slot examined by the gadget enumeration */
static double pLocal = 0.3; /* share of random group moves that are local (same apartment: phase / orientation) */
static long nG = 0;         /* accepted group moves */
static double unlock_t = -1; /* when >= 0: after this many seconds the groups are dissolved into free flags
                                (the chords stay in F); the gadget structure then only serves as the start */

static int step_towards(int cur, int target) {
    for (int k = 0; k < DEG; k++) if (distH[adjH[cur][k]][target] == distH[cur][target] - 1) return adjH[cur][k];
    fprintf(stderr, "step_towards failed\n"); exit(1);
}
static void build_apartment_tables(void) {
    double t0 = now_sec();
    mid5 = malloc(sizeof(uint16_t) * (size_t)NPAIR);
    opp = malloc(sizeof(uint16_t) * (size_t)NE * NOPP);
    gtab = malloc(sizeof(uint16_t) * (size_t)NE * NOPP);
    htab = malloc(sizeof(uint16_t) * (size_t)NE * NOPP);
    if (!mid5 || !opp || !gtab || !htab) { fprintf(stderr, "alloc failed\n"); exit(1); }
    memset(mid5, 0xff, sizeof(uint16_t) * (size_t)NPAIR);
    long n5 = 0;
    for (int u = 0; u < NOLD; u++)
        for (int v = 0; v < NOLD; v++) {
            if (distH[u][v] != 5) continue;
            int v2 = step_towards(u, v), v3 = step_towards(v2, v), v4 = step_towards(v3, v);
            mid5[(size_t)u * NOLD + v] = (uint16_t)find_edge(v3, v4);
            n5++;
        }
    for (int u = 0; u < NOLD; u++)
        for (int v = u + 1; v < NOLD; v++)
            if (distH[u][v] == 5 && mid5[(size_t)u * NOLD + v] != mid5[(size_t)v * NOLD + u]) { fprintf(stderr, "mid5 not symmetric\n"); exit(1); }
    for (int e = 0; e < NE; e++) {
        int c = 0;
        for (int f = 0; f < NE; f++)
            if (edge_dist(e, f) == 5) {
                if (c >= NOPP) { fprintf(stderr, "more than %d opposite flags\n", NOPP); exit(1); }
                opp[(size_t)e * NOPP + c] = (uint16_t)f;
                /* e = (p,L), f = (p',L'): points < lines, so ep = point, el = line */
                gtab[(size_t)e * NOPP + c] = mid5[(size_t)el[e] * NOLD + ep[f]];
                htab[(size_t)e * NOPP + c] = mid5[(size_t)el[f] * NOLD + ep[e]];
                c++;
            }
        if (c != NOPP) { fprintf(stderr, "flag %d has %d opposite flags\n", e, c); exit(1); }
    }
    fprintf(stderr, "apartment tables: %ld ordered pairs at distance 5, %d x %d opposite pairs, %.1fs\n", n5, NE, NOPP, now_sec() - t0);
}
/* the 12 flags of the apartment through the opposite flags e (position 0) and f (position 6), in cyclic order */
static void apartment_flags(int e, int f, int *fl) {
    int v[12];
    v[0] = ep[e]; v[1] = el[e]; v[6] = ep[f]; v[7] = el[f];
    for (int i = 2; i <= 5; i++) v[i] = step_towards(v[i - 1], v[6]);
    for (int i = 8; i <= 11; i++) v[i] = step_towards(v[i - 1], v[0]);
    for (int i = 0; i < 12; i++) {
        fl[i] = find_edge(v[i], v[(i + 1) % 12]);
        if (fl[i] < 0) { fprintf(stderr, "apartment_flags: not a cycle\n"); exit(1); }
    }
}
/* check that slots 4k..4k+3 hold a gadget: (a,c), (b,d) opposite and {b,d} = the middle flags of the two
   geodesics closing the apartment through a and c */
static int check_group(int k) {
    int a = S[4 * k], b = S[4 * k + 1], c = S[4 * k + 2], d = S[4 * k + 3];
    if (edge_dist(a, c) != 5 || edge_dist(b, d) != 5) return 0;
    int g = mid5[(size_t)el[a] * NOLD + ep[c]], h = mid5[(size_t)el[c] * NOLD + ep[a]];
    return (b == g && d == h) || (b == h && d == g);
}

typedef struct { int x[4]; } Placement;   /* flags for slots a, b, c, d of a group */
/* the 8 orientations of the gadget with antipodal pairs (e,f), (g,h) */
static void gadget_orientations(int e, int f, int g, int h, Placement *pl) {
    int n = 0;
    int pairs[2][2] = { { e, f }, { g, h } };
    for (int pa = 0; pa < 2; pa++)
        for (int oa = 0; oa < 2; oa++)
            for (int ob = 0; ob < 2; ob++) {
                int *A = pairs[pa], *B = pairs[1 - pa];
                pl[n].x[0] = A[oa]; pl[n].x[2] = A[1 - oa];
                pl[n].x[1] = B[ob]; pl[n].x[3] = B[1 - ob];
                n++;
            }
}
static inline int flag_free_for(int e, int k) { return !inS[e] || sidx[e] / 4 == k; }

/* candidate gadget placements for group k by the min-conflicts proxy: for each slot s, the topT edges by
   proxy (optionally restricted to edges within rR of vertex rv when s == rslot) and every opposite flag
   as its chord partner; the two other gadget flags fill the other chord in both orders.  Keeps the K best
   scores (random tie-break, duplicates removed).  Returns the number of candidates. */
static int prox[4][NE];
static int group_candidates(int k, int rslot, int rv, int rR, Placement *cand, int K) {
    proxy_excl_lo = 4 * k; proxy_excl_hi = 4 * k + 3;
    for (int t = 0; t < 4; t++) { relocation_proxy(4 * k + t); memcpy(prox[t], proxy, sizeof(int) * NE); }
    proxy_excl_lo = proxy_excl_hi = -1;
    uint64_t keys[64];
    int nc = 0;
    if (K > 64) K = 64;
    for (int s = 0; s < 4; s++) {
        if (rslot >= 0 && s != rslot) continue;
        int sp = s ^ 2, s2 = (s + 1) & 3, s3 = (s + 3) & 3;   /* chord partner; the other chord's slots */
        /* topT edges for slot s */
        int top[64]; uint32_t tk[64]; int nt = 0;
        int T = topT < 64 ? topT : 64;
        for (int e = 0; e < NE; e++) {
            if (!flag_free_for(e, k)) continue;
            if (rslot >= 0) {
                int d = distH[ep[e]][rv] < distH[el[e]][rv] ? distH[ep[e]][rv] : distH[el[e]][rv];
                if (d > rR) continue;
            }
            uint32_t key = ((uint32_t)prox[s][e] << 12) | (uint32_t)(rnd64() & 4095);
            if (nt < T) { top[nt] = e; tk[nt] = key; nt++; }
            else {
                int worst = 0;
                for (int c = 1; c < nt; c++) if (tk[c] > tk[worst]) worst = c;
                if (key < tk[worst]) { top[worst] = e; tk[worst] = key; }
            }
        }
        for (int t = 0; t < nt; t++) {
            int e = top[t];
            const uint16_t *O = opp + (size_t)e * NOPP, *G = gtab + (size_t)e * NOPP, *Hh = htab + (size_t)e * NOPP;
            for (int j = 0; j < NOPP; j++) {
                int f = O[j], g = G[j], h = Hh[j];
                if (!flag_free_for(f, k) || !flag_free_for(g, k) || !flag_free_for(h, k)) continue;
                int base_sc = prox[s][e] + prox[sp][f];
                for (int o = 0; o < 2; o++) {
                    int xg = o ? h : g, xh = o ? g : h;
                    int sc = base_sc + prox[s2][xg] + prox[s3][xh];
                    uint64_t key = ((uint64_t)sc << 20) | (rnd64() & 0xfffff);
                    int worst = 0;
                    if (nc == K) {
                        for (int c = 1; c < nc; c++) if (keys[c] > keys[worst]) worst = c;
                        if (key >= keys[worst]) continue;
                    }
                    Placement p;
                    p.x[s] = e; p.x[sp] = f; p.x[s2] = xg; p.x[s3] = xh;
                    if (p.x[0] == S[4 * k] && p.x[1] == S[4 * k + 1] && p.x[2] == S[4 * k + 2] && p.x[3] == S[4 * k + 3]) continue;   /* identity */
                    int dup = 0;
                    for (int c = 0; c < nc && !dup; c++) dup = !memcmp(&cand[c], &p, sizeof p);
                    if (dup) continue;
                    if (nc < K) { cand[nc] = p; keys[nc] = key; nc++; }
                    else { cand[worst] = p; keys[worst] = key; }
                }
            }
        }
    }
    return nc;
}

/* evaluate the candidate placements of group k exactly, apply the best, return 2 (tm_cost set, undo buffer
   filled) or 0 when there is none */
static int group_apply_best(const Fgraph *F, int k, const Placement *cand, int nc) {
    if (nc == 0) return 0;
    int old[4];
    for (int t = 0; t < 4; t++) old[t] = S[4 * k + t];
    int best = -1, bestc = 1 << 30;
    for (int c = 0; c < nc; c++) {
        set_group(k, cand[c].x);
        build_sder(tm_sd_trial);
        int cost = evaluate(tm_sd_trial, F, ev_trial);
        set_group(k, old);
        if (cost < bestc) { bestc = cost; best = c; }
    }
    set_group(k, cand[best].x);
    build_sder(tm_sd_trial);
    tm_cost = evaluate(tm_sd_trial, F, ev_trial);
    su_n = 4;
    for (int t = 0; t < 4; t++) { su_i[t] = 4 * k + t; su_e[t] = old[t]; }
    return 2;
}
static int group_relocate_mc(const Fgraph *F, int k, int rslot, int rv, int rR) {
    Placement cand[64];
    int nc = group_candidates(k, rslot, rv, rR, cand, kS);
    return group_apply_best(F, k, cand, nc);
}
/* local move: the 23 other placements of the current apartment (3 phases x 8 orientations) */
static int group_local_move(const Fgraph *F, int k) {
    int fl[12];
    apartment_flags(S[4 * k], S[4 * k + 2], fl);
    Placement cand[24];
    int nc = 0;
    for (int ph = 0; ph < 3; ph++) {
        Placement o[8];
        gadget_orientations(fl[ph], fl[ph + 6], fl[ph + 3], fl[ph + 9], o);
        for (int i = 0; i < 8; i++) {
            int ok = 1, same = 1;
            for (int t = 0; t < 4; t++) { if (!flag_free_for(o[i].x[t], k)) ok = 0; if (o[i].x[t] != S[4 * k + t]) same = 0; }
            if (ok && !same) cand[nc++] = o[i];
        }
    }
    return group_apply_best(F, k, cand, nc);
}
/* random gadgets for the apartment groups and random flags for the free slots, pairwise edge-distance
   >= spread across groups while possible */
static void random_S_apts(int *out, int spread, int nfixed) {
    /* out[0..nfixed-1] are given (nfixed a multiple of 4 or >= 4*napt); the rest is generated */
    static uint8_t used[NE];
    memset(used, 0, sizeof used);
    for (int i = 0; i < nfixed; i++) used[out[i]] = 1;
    for (int k = nfixed / 4; k < napt; k++) {
        for (int tries = 0;; tries++) {
            int e = rndint(NE), j = rndint(NOPP);
            int x[4] = { e, gtab[(size_t)e * NOPP + j], opp[(size_t)e * NOPP + j], htab[(size_t)e * NOPP + j] };
            int ok = 1;
            for (int t = 0; t < 4 && ok; t++) if (used[x[t]]) ok = 0;
            if (ok && spread > 0 && tries < 20000)
                for (int t = 0; t < 4 && ok; t++)
                    for (int i = 0; i < 4 * k && ok; i++) if (edge_dist(out[i], x[t]) < spread) ok = 0;
            if (!ok) continue;
            for (int t = 0; t < 4; t++) { out[4 * k + t] = x[t]; used[x[t]] = 1; }
            break;
        }
    }
    int i0 = nfixed > 4 * napt ? nfixed : 4 * napt;
    for (int i = i0; i < ms; i++) {
        for (int tries = 0;; tries++) {
            int e = rndint(NE);
            if (used[e]) continue;
            int ok = 1;
            if (spread > 0 && tries < 20000)
                for (int j = 0; j < i; j++) if (edge_dist(out[j], e) < spread) { ok = 0; break; }
            if (!ok) continue;
            used[e] = 1; out[i] = e;
            break;
        }
    }
}

/* S edges on the unique H-geodesic from v to u (d_H(u,v) <= 5); returns count */
static int s_edges_on_geodesic(int u, int v, int *out) {
    int cnt = 0, cur = v;
    while (cur != u) {
        int found = -1;
        for (int k = 0; k < DEG; k++) if (distH[u][adjH[cur][k]] == distH[u][cur] - 1) { found = k; break; }
        if (found < 0) return cnt;
        int e = edge_of[cur][found];
        if (inS[e]) out[cnt++] = e;
        cur = adjH[cur][found];
    }
    return cnt;
}

static int random_target_edge(int e_old) {
    int e;
    if (rnd01() < 0.5) e = near2[offnear[e_old] + rndint(offnear[e_old + 1] - offnear[e_old])];
    else e = rndint(NE);
    return e;
}

/* Propose a targeted move for the current state (violations recorded by the last evaluate of it).
   Returns 0 (nothing to do / invalid), 1 (F move applied, undo with f_undo), 2 (S move applied:
   *s_i and *s_old tell which S entry changed and its former edge). pSfix = probability of an S repair. */
/* try f_create_edge(i, k) for every candidate k and keep the best (exhaustive variant for the endgame) */
static int f_create_edge_best(Fgraph *F, int i, const int *cand, int nc) {
    int bestk = -1, bestc = 1 << 30, besta = -1, bestb = -1;
    for (int c = 0; c < nc; c++) {
        int k = cand[c];
        if (!f_create_edge(F, i, k)) continue;
        if (tm_cost < bestc) {
            bestc = tm_cost; bestk = k;
            /* remember the swap partners: the log holds del(i,a) del(k,b) add(i,k) add(a,b) */
            besta = ulog_b[0]; bestb = ulog_b[1];
        }
        f_undo(F);
    }
    if (bestk < 0) return 0;
    ulog_n = 0;
    lf_del(F, i, besta); lf_del(F, bestk, bestb); lf_add(F, i, bestk); lf_add(F, besta, bestb);
    tm_cost = evaluate(tm_sd, F, ev_trial);
    return 1;
}
static int exhaustive_k = 0;   /* set by the main loop in the endgame */

static int targeted_move(const SDer *sd, SDer *sd_trial, Fgraph *F, const Eval *ev, double pSfix, int *s_i, int *s_old) {
    int cost = ev->no + ev->nn + ev->oo;
    if (cost <= 0) return 0;
    tm_sd = sd; tm_sd_trial = sd_trial;
    int r = rndint(cost);
    int cand[MAXNEW], nc = 0;
    int sc[64], nsc = 0;
    if (r < ev->no) {
        /* new-old violation (x_i, v) */
        int i = 0;
        while (r >= ev->nocnt[i]) { r -= ev->nocnt[i]; i++; }
        int v = -1;
        for (int w = 0; w <= OLDWORDS; w++) {
            int pc = __builtin_popcountll(ev->uncov[i][w]);
            if (r < pc) {
                uint64_t word = ev->uncov[i][w];
                while (r--) word &= word - 1;
                v = w * 64 + __builtin_ctzll(word);
                break;
            }
            r -= pc;
        }
        if (v < 0) return 0;
        if (rnd01() >= pSfix) {
            /* F repair: new F edge (i,k) with d0[k][v] <= 5 */
            for (int k = 0; k < nnew; k++)
                if (k != i && ((sd->bits[k][5][v >> 6] >> (v & 63)) & 1) && !f_has(F, i, k)) cand[nc++] = k;
            if (nc == 0) return 0;
            return exhaustive_k ? f_create_edge_best(F, i, cand, nc) : f_create_edge(F, i, cand[rndint(nc)]);
        }
        /* S repair: move a blocking S edge (or e_i itself), or (mc mode, half of the time) bring a
           vertex at F-distance 1 or 2 from x_i to an edge within distance 5 - DF of v so that it covers v */
        if (i >= ms) return 0;
        if (mc_mode && rnd01() < 0.5) {
            int dmax = rnd01() < 0.5 ? 1 : 2;
            for (int k = 0; k < ms; k++) if (k != i && DF[i][k] >= 1 && DF[i][k] <= dmax) cand[nc++] = k;
            if (nc == 0) return 0;
            int k = cand[rndint(nc)];
            if (k < 4 * napt) {   /* bring the group of k so that slot k covers v */
                if (group_tabu(k / 4)) return 0;
                return group_relocate_mc(F, k / 4, k % 4, v, 5 - DF[i][k]);
            }
            if (tabu_len > 0 && tabu_until[k] > n_smoves) return 0;
            reloc_v = v; reloc_R = 5 - DF[i][k];
            int res = s_relocate_mc(F, k, s_i, s_old);
            reloc_v = -1;
            return res;
        }
        sc[nsc++] = S[i];
        int ends[2] = { ep[S[i]], el[S[i]] };
        for (int s = 0; s < 2; s++)
            if (distH[ends[s]][v] <= 5) nsc += s_edges_on_geodesic(ends[s], v, sc + nsc);
        return s_move_best(F, sc[rndint(nsc)], s_i, s_old);
    }
    r -= ev->no;
    if (r < ev->nn) {
        int i = ev->viol_nn_a[r], j = ev->viol_nn_b[r];
        if (rnd01() < 0.5) { int t = i; i = j; j = t; }
        if (rnd01() >= pSfix) {
            for (int k = 0; k < nnew; k++)
                if (k != i && ev->D[k][j] <= 5 && !f_has(F, i, k)) cand[nc++] = k;
            if (nc == 0) return 0;
            return exhaustive_k ? f_create_edge_best(F, i, cand, nc) : f_create_edge(F, i, cand[rndint(nc)]);
        }
        if (i >= ms) return 0;
        return s_move_best(F, S[i], s_i, s_old);
    }
    r -= ev->nn;
    if (r >= ev->oo) return 0;
    int pid = ev->viol_oo[r];
    int a = pid / NOLD, b = pid % NOLD;
    if (rnd01() >= pSfix) {
        /* F repair: edge (i,j) with i near a, j near b, da + 1 + db <= 6 */
        int ci[MAXNEW], cj[MAXNEW];
        int na = sd->nearcnt[a], nb = sd->nearcnt[b];
        for (int ia = 0; ia < na; ia++)
            for (int jb = 0; jb < nb; jb++) {
                int i = sd->nearid[a][ia], j = sd->nearid[b][jb];
                if (i != j && sd->neard[a][ia] + 1 + sd->neard[b][jb] <= 6 && !f_has(F, i, j) && nc < MAXNEW) { ci[nc] = i; cj[nc] = j; nc++; }
            }
        if (nc == 0) return 0;
        if (exhaustive_k) {
            int bestc = 1 << 30, bi = -1, bj = -1, ba = -1, bb = -1;
            for (int c = 0; c < nc; c++) {
                if (!f_create_edge(F, ci[c], cj[c])) continue;
                if (tm_cost < bestc) { bestc = tm_cost; bi = ci[c]; bj = cj[c]; ba = ulog_b[0]; bb = ulog_b[1]; }
                f_undo(F);
            }
            if (bi < 0) return 0;
            ulog_n = 0;
            lf_del(F, bi, ba); lf_del(F, bj, bb); lf_add(F, bi, bj); lf_add(F, ba, bb);
            tm_cost = evaluate(tm_sd, F, ev_trial);
            return 1;
        }
        int t = rndint(nc);
        return f_create_edge(F, ci[t], cj[t]);
    }
    /* S repair: move one of the S edges on the a-b geodesic(s) */
    int d = distH[a][b];
    if (d <= 5) nsc = s_edges_on_geodesic(a, b, sc);
    else {
        for (int g = 0; g < DEG; g++) {
            int w0 = adjH[a][g];
            if (inS[edge_of[a][g]]) sc[nsc++] = edge_of[a][g];
            nsc += s_edges_on_geodesic(w0, b, sc + nsc);
        }
    }
    if (nsc == 0) return 0;
    return s_move_best(F, sc[rndint(nsc)], s_i, s_old);
}

/* ---------- reference implementations ---------- */
static uint16_t adjG[MAXN][DEG];
static uint8_t degG[MAXN];

static void build_G(const Fgraph *F, int withF) {
    for (int u = 0; u < NOLD; u++) { memcpy(adjG[u], adjH[u], sizeof adjH[u]); degG[u] = DEG; }
    for (int i = 0; i < nnew; i++) degG[NOLD + i] = 0;
    for (int i = 0; i < ms; i++) {
        int e = S[i], p = ep[e], l = el[e], x = NOLD + i;
        for (int k = 0; k < DEG; k++) if (adjG[p][k] == l) adjG[p][k] = (uint16_t)x;
        for (int k = 0; k < DEG; k++) if (adjG[l][k] == p) adjG[l][k] = (uint16_t)x;
        adjG[x][0] = (uint16_t)p; adjG[x][1] = (uint16_t)l; degG[x] = 2;
    }
    if (withF)
        for (int i = 0; i < nnew; i++) {
            int x = NOLD + i;
            for (int k = 0; k < F->fdeg[i]; k++) {
                if (degG[x] >= DEG) { fprintf(stderr, "degree overflow at new %d\n", i); exit(1); }
                adjG[x][degG[x]++] = (uint16_t)(NOLD + F->fadj[i][k]);
            }
        }
}

static void bfs_G(int s, uint8_t *dist, int limit) {
    static int q[MAXN];
    memset(dist, UNREACH, N);
    int qh = 0, qt = 0;
    dist[s] = 0; q[qt++] = s;
    while (qh < qt) {
        int u = q[qh++];
        int du = dist[u];
        if (du >= limit) continue;
        for (int k = 0; k < degG[u]; k++) {
            int w = adjG[u][k];
            if (dist[w] == UNREACH) { dist[w] = (uint8_t)(du + 1); q[qt++] = w; }
        }
    }
}

/* the method of the q=3 code: BFS from every new vertex in G, boo pairs via d0 + dG */
static uint8_t dGref[MAXNEW][MAXN], d0ref[MAXNEW][MAXN];
static int evaluate_bfs(const Fgraph *F, int *pno, int *pnn, int *poo) {
    build_G(F, 0);
    for (int i = 0; i < nnew; i++) bfs_G(NOLD + i, d0ref[i], 6);
    build_G(F, 1);
    for (int i = 0; i < nnew; i++) bfs_G(NOLD + i, dGref[i], 6);
    int cno = 0, cnn = 0, coo = 0;
    for (int i = 0; i < nnew; i++) {
        for (int v = 0; v < NOLD; v++) cno += (dGref[i][v] > 6);
        for (int j = i + 1; j < nnew; j++) cnn += (dGref[i][NOLD + j] > 6);
    }
    for (int t = 0; t < nBoo; t++) {
        int a = boo[t] / NOLD, b = boo[t] % NOLD, best = 99;
        for (int i = 0; i < nnew; i++) {
            int v = d0ref[i][a] + dGref[i][b];
            if (v < best) best = v;
        }
        coo += (best > 6);
    }
    *pno = cno; *pnn = cnn; *poo = coo;
    return cno + cnn + coo;
}

/* all-pairs BFS on the explicit graph; returns #pairs > 6 */
static int exact_cost(const Fgraph *F, int *maxdeg, int *mindeg) {
    build_G(F, 1);
    *maxdeg = 0; *mindeg = 99;
    for (int v = 0; v < N; v++) { if (degG[v] > *maxdeg) *maxdeg = degG[v]; if (degG[v] < *mindeg) *mindeg = degG[v]; }
    static uint8_t dist[MAXN];
    int bad = 0;
    for (int s = 0; s < N; s++) {
        bfs_G(s, dist, 30);
        for (int t = s + 1; t < N; t++) bad += (dist[t] > 6);
    }
    return bad;
}

static void write_edges(const Fgraph *F, const char *path) {
    build_G(F, 1);
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }
    for (int u = 0; u < N; u++)
        for (int k = 0; k < degG[u]; k++) {
            int w = adjG[u][k];
            if (u < w) fprintf(f, "%d %d\n", u, w);
        }
    fclose(f);
}

static void print_state(const Fgraph *F, FILE *f) {
    fprintf(f, "S:");
    for (int i = 0; i < ms; i++) fprintf(f, " (%d,%d)", ep[S[i]], el[S[i]]);
    fprintf(f, "\nF:");
    for (int i = 0; i < nnew; i++)
        for (int k = 0; k < F->fdeg[i]; k++) if (i < F->fadj[i][k]) fprintf(f, " (%d,%d)", i, F->fadj[i][k]);
    fprintf(f, "\n");
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

static int edge_dist(int e, int f) {
    int d = distH[ep[e]][ep[f]], t;
    t = distH[ep[e]][el[f]]; if (t < d) d = t;
    t = distH[el[e]][ep[f]]; if (t < d) d = t;
    t = distH[el[e]][el[f]]; if (t < d) d = t;
    return d;
}

static void random_S(int *out, int spread) {
    /* spread > 0: greedy random with pairwise edge-distance >= spread while possible */
    static uint8_t used[NE];
    memset(used, 0, sizeof used);
    for (int i = 0; i < ms; i++) {
        int e, tries = 0;
        for (;;) {
            e = rndint(NE);
            if (used[e]) continue;
            int ok = 1;
            if (spread > 0 && tries < 20000)
                for (int j = 0; j < i; j++) if (edge_dist(out[j], e) < spread) { ok = 0; break; }
            if (ok) break;
            tries++;
        }
        used[e] = 1; out[i] = e;
    }
}

/* ---------- validation ---------- */
static void validate(int K, int spread_mode) {
    Fgraph F;
    static SDer sd;
    static uint8_t dist[MAXN];
    static uint8_t *mark;
    mark = calloc(NPAIR, 1);
    int newS[MAXNEW];
    for (int k = 0; k < K; k++) {
        /* S: uniform, or clustered (edges drawn from near2 lists of chosen edges) */
        int clustered = (k % 2 == 1);
        if (!clustered) random_S(newS, spread_mode);
        else {
            static uint8_t used[NE];
            memset(used, 0, sizeof used);
            newS[0] = rndint(NE); used[newS[0]] = 1;
            for (int i = 1; i < ms; i++) {
                int e;
                do {
                    int src = newS[rndint(i)];
                    e = (rnd01() < 0.8) ? near2[offnear[src] + rndint(offnear[src + 1] - offnear[src])] : rndint(NE);
                } while (used[e]);
                used[e] = 1; newS[i] = e;
            }
        }
        install_S(newS);
        random_F(&F);
        build_sder(&sd);
        /* (1) rows vs BFS in G0 */
        build_G(&F, 0);
        int rowerr = 0;
        for (int i = 0; i < nnew; i++) {
            bfs_G(NOLD + i, dist, 6);
            for (int v = 0; v < N; v++) {
                int a = dist[v] > 6 ? FAR : dist[v], b = sd.d0[i][v];
                if (a != b) { if (rowerr < 5) fprintf(stderr, "row mismatch i=%d v=%d bfs=%d fast=%d\n", i, v, a, b); rowerr++; }
                /* bits consistency */
                for (int r = 0; r <= 6; r++) {
                    int inb = (sd.bits[i][r][v >> 6] >> (v & 63)) & 1;
                    if (inb != (b <= r)) { if (rowerr < 5) fprintf(stderr, "bits mismatch i=%d v=%d r=%d\n", i, v, r); rowerr++; }
                }
            }
        }
        /* (2) boo list vs BFS in G0 from old vertices */
        for (int t = 0; t < nBoo; t++) mark[boo[t]] = 1;
        int nbad0 = 0, booerr = 0;
        for (int a = 0; a < NOLD; a++) {
            bfs_G(a, dist, 6);
            for (int b = a + 1; b < NOLD; b++) {
                int far = dist[b] > 6;
                nbad0 += far;
                if (far != mark[a * NOLD + b]) booerr++;
            }
        }
        for (int t = 0; t < nBoo; t++) mark[boo[t]] = 0;
        /* (3) costs */
        int cf = evaluate(&sd, &F, ev_cur);
        int fno = ev_cur->no, fnn = ev_cur->nn, foo = ev_cur->oo;
        /* violation records must be consistent with the counts */
        { int s = 0; for (int i = 0; i < nnew; i++) s += ev_cur->nocnt[i]; if (s != fno) { fprintf(stderr, "nocnt mismatch\n"); exit(2); } }
        int bno, bnn, boo_;
        int cb = evaluate_bfs(&F, &bno, &bnn, &boo_);
        int mx, mn;
        int ce = exact_cost(&F, &mx, &mn);
        printf("validate %d: %s m=%d h=%d nBoo=%d (bfs %d, errors %d) rowerr=%d | fast %d (no %d nn %d oo %d) | bfs %d (no %d nn %d oo %d) | exact %d deg %d..%d\n",
               k, clustered ? "clustered" : "uniform  ", ms, hb, nBoo, nbad0, booerr, rowerr, cf, fno, fnn, foo, cb, bno, bnn, boo_, ce, mn, mx);
        fflush(stdout);
        if (rowerr || booerr || cf != cb || cf != ce) { fprintf(stderr, "VALIDATION FAILED\n"); exit(2); }
    }
    printf("validation passed (%d instances)\n", K);
    free(mark);
}

/* ---------- --exhx: exact 1-, 2- and 3-move neighbourhoods of a low-cost state ----------
   Level 1: every single-flag relocation (as --exh1) and, with --swaps, every pairwise label swap
   S[i] <-> S[j] (same flag set, different F roles; not reachable by one or two relocations because
   the intermediate state would put two flags on one edge) and optionally every 3-cycle of labels.
   Level 2: from every level-1 state with cost <= c0 + slack1 (lowest cost first) every "fixing move"
   of its remaining violations: (a) every relocation of an involved flag (the far new vertex x_k, the
   S edges blocking the geodesics from the ends of S[k] to the far old vertex w, the S edges blocking
   the near-miss paths through another new vertex x_l with D[k][l] + 1 + d_H(end of S[l], w) <= 6,
   both ends of a new-new pair, the S edges on the geodesics of an old-old pair), and (b) every
   relocation of a flag j with D[k][j] <= 4 to an edge with an end within 5 - D[k][j] of w (a shortcut
   x_k .. x_j - u .. w; the analogous bridges for new-new and old-old pairs).  Level 3: the same fixing
   moves from every level-2 state with cost <= slack2.  Phases 2 and 3 stop at the wall-clock budget
   (the seconds argument).  Cost-0 states are checked by all-pairs BFS and written as edge lists. */
#define XC_MAX 800000
#define XW ((NE + 63) / 64)
static uint64_t *xc_bits;
static int xc_slot[XC_MAX], xc_edge[XC_MAX], xc_n, xc_ninv;
static void xc_add(int slot, int e) {
    if (slot < 0 || slot >= ms || e < 0 || e >= NE || inS[e] || xc_n >= XC_MAX) return;
    uint64_t *w = xc_bits + (size_t)slot * XW + (e >> 6), m = 1ULL << (e & 63);
    if (*w & m) return;
    *w |= m; xc_slot[xc_n] = slot; xc_edge[xc_n] = e; xc_n++;
}
static void xc_clear(void) {
    for (int t = 0; t < xc_n; t++) xc_bits[(size_t)xc_slot[t] * XW + (xc_edge[t] >> 6)] = 0;
    xc_n = 0;
}
static void xc_mark(uint8_t *mark, int *inv, int *ninv, int k) {
    if (k >= 0 && k < ms && !mark[k]) { mark[k] = 1; inv[(*ninv)++] = k; }
}
/* candidate fixing moves for the violations recorded in ev (evaluated on the current S) */
static void xc_fixing(const Eval *ev) {
    static uint8_t mark[MAXNEW];
    memset(mark, 0, sizeof mark);
    int inv[MAXNEW], ninv = 0, sc[64], nsc;
    xc_clear();
    for (int k = 0; k < ms; k++) {
        if (ev->nocnt[k] == 0) continue;
        xc_mark(mark, inv, &ninv, k);
        for (int w = 0; w <= OLDWORDS; w++) {
            uint64_t word = ev->uncov[k][w];
            while (word) {
                int v = w * 64 + __builtin_ctzll(word); word &= word - 1;
                /* blockers: S edges on the geodesics from the ends of S[k] (l = k) and on the near-miss
                   paths x_k .. x_l - end .. v through other new vertices */
                for (int l = 0; l < ms; l++) {
                    int dkl = l == k ? 0 : ev->D[k][l];
                    if (dkl > 4) continue;
                    int ends[2] = { ep[S[l]], el[S[l]] };
                    for (int s = 0; s < 2; s++) {
                        if (distH[ends[s]][v] > 5 || dkl + 1 + distH[ends[s]][v] > 6) continue;
                        nsc = s_edges_on_geodesic(ends[s], v, sc);
                        for (int q = 0; q < nsc; q++) xc_mark(mark, inv, &ninv, sidx[sc[q]]);
                    }
                }
                /* shortcut providers: x_j at an edge with an end u, D[k][j] + 1 + d(u, v) <= 6 */
                for (int j = 0; j < ms; j++) {
                    if (j == k) continue;
                    int R = 5 - ev->D[k][j];
                    if (R < 1) continue;
                    for (int e = 0; e < NE; e++) {
                        if (inS[e]) continue;
                        int d = distH[ep[e]][v] < distH[el[e]][v] ? distH[ep[e]][v] : distH[el[e]][v];
                        if (d <= R) xc_add(j, e);
                    }
                }
            }
        }
    }
    for (int q = 0; q < ev->nn; q++) {
        int ka = ev->viol_nn_a[q], kb = ev->viol_nn_b[q];
        xc_mark(mark, inv, &ninv, ka); xc_mark(mark, inv, &ninv, kb);
        for (int side = 0; side < 2; side++) {
            int a = side ? kb : ka, b = side ? ka : kb;
            if (b >= ms) continue;
            int ends[2] = { ep[S[b]], el[S[b]] };
            for (int j = 0; j < ms; j++) {
                if (j == a || j == b) continue;
                int R = 4 - ev->D[a][j];   /* D[a][j] + 1 + d(u, end of S[b]) + 1 <= 6 */
                if (R < 0) continue;
                for (int e = 0; e < NE; e++) {
                    if (inS[e]) continue;
                    int d = 99;
                    for (int s = 0; s < 2; s++) { int t = distH[ep[e]][ends[s]]; if (t < d) d = t; t = distH[el[e]][ends[s]]; if (t < d) d = t; }
                    if (d <= R) xc_add(j, e);
                }
            }
        }
    }
    for (int q = 0; q < ev->oo; q++) {
        int pid = ev->viol_oo[q], a = pid / NOLD, b = pid % NOLD;
        nsc = 0;
        if (distH[a][b] <= 5) nsc = s_edges_on_geodesic(a, b, sc);
        else for (int g = 0; g < DEG; g++) { if (inS[edge_of[a][g]]) sc[nsc++] = edge_of[a][g]; nsc += s_edges_on_geodesic(adjH[a][g], b, sc + nsc); }
        for (int x = 0; x < nsc; x++) xc_mark(mark, inv, &ninv, sidx[sc[x]]);
        /* a flag on an edge e = (u1, u2) with d(a, u) + d(u', b) <= 4 bridges a and b in 6 */
        for (int e = 0; e < NE; e++) {
            if (inS[e]) continue;
            int u[2] = { ep[e], el[e] }, ok = 0;
            for (int s = 0; s < 2 && !ok; s++) for (int t = 0; t < 2 && !ok; t++) if (distH[a][u[s]] + distH[u[t]][b] <= 4) ok = 1;
            if (ok) for (int j = 0; j < ms; j++) xc_add(j, e);
        }
    }
    for (int x = 0; x < ninv; x++) for (int e = 0; e < NE; e++) xc_add(inv[x], e);
    xc_ninv = ninv;
}
typedef struct { int kind, a, b; } Op;   /* kind 0: slot a -> edge b; kind 1: swap the labels of slots a and b */
typedef struct { Op op[3]; int nop, c; } XState;
static void swap_S(int i, int j) {
    int ei = S[i], ej = S[j];
    S[i] = ej; S[j] = ei; sidx[ej] = i; sidx[ei] = j;   /* same edge set: inS and the boo structure are unchanged */
}
static int op_apply(const Op *op, int *undo) {
    if (op->kind == 0) { if (inS[op->b]) return 0; *undo = S[op->a]; set_S_edge(op->a, op->b); }
    else swap_S(op->a, op->b);
    return 1;
}
static void op_undo(const Op *op, int undo) {
    if (op->kind == 0) set_S_edge(op->a, undo); else swap_S(op->a, op->b);
}
static void op_print(const Op *op) {
    if (op->kind == 0) printf("slot %d -> (%d,%d)", op->a, ep[op->b], el[op->b]);
    else printf("swap slots %d <-> %d", op->a, op->b);
}
static long xz_count = 0;
static void xz_report(const Fgraph *F, const XState *st, const char *outdir, const char *tag) {
    int mx, mn;
    int ex = exact_cost(F, &mx, &mn);
    printf("ZERO COST (%d ops:", st->nop);
    for (int t = 0; t < st->nop; t++) { printf(" "); op_print(&st->op[t]); if (t + 1 < st->nop) printf(","); }
    printf("): exact all-pairs check pairs>6 = %d, maxdeg %d, mindeg %d, n = %d\n", ex, mx, mn, N);
    char path[1024];
    snprintf(path, sizeof path, "%s/%s_exhx_n%d_L%d_%ld.edges", outdir, tag, N, st->nop, ++xz_count);
    write_edges(F, path);
    printf("wrote %s\n", path);
    print_state(F, stdout);
    fflush(stdout);
}
static int xs_cmp(const void *x, const void *y) { return ((const XState *)x)->c - ((const XState *)y)->c; }
#define XS_MAX 65536
static XState xl1[XS_MAX], xl2[XS_MAX];
/* evaluate every fixing move of the current state (whose violations are in ev_trial); records states with
   cost <= keep into out (if out != NULL); returns the best cost seen */
static int xc_expand(SDer *sd_trial, const Fgraph *F, const XState *base, int keep, XState *out, int *nout, long *neval, int *hist,
                     double budget, double t_start, const char *outdir, const char *tag) {
    xc_fixing(ev_trial);
    int best = 1 << 30, last = base->nop ? base->op[base->nop - 1].a : -1, lastkind = base->nop ? base->op[base->nop - 1].kind : -1;
    for (int t = 0; t < xc_n; t++) {
        int j = xc_slot[t], e = xc_edge[t];
        if (lastkind == 0 && j == last) continue;   /* moving the last-moved flag again is a different shorter sequence */
        if (inS[e]) continue;
        int old = S[j];
        set_S_edge(j, e); build_sder(sd_trial);
        int c = evaluate(sd_trial, F, ev_trial);
        (*neval)++; hist[c < 7 ? c : 7]++;
        if (c < best) best = c;
        XState st = *base; st.op[st.nop].kind = 0; st.op[st.nop].a = j; st.op[st.nop].b = e; st.nop++; st.c = c;
        if (c == 0) xz_report(F, &st, outdir, tag);
        if (out && c <= keep && *nout < XS_MAX) out[(*nout)++] = st;
        set_S_edge(j, old);
        if ((t & 1023) == 1023 && now_sec() - t_start > budget) { printf("  (budget reached inside the candidate scan after %d of %d)\n", t + 1, xc_n); break; }
    }
    return best;
}
static void exhx_run(SDer *sd_trial, Fgraph *F, int c0, int slack1, int slack2, int max1, int max2, int swaps, double budget, double t_start,
                     const char *outdir, const char *tag) {
    xc_bits = calloc((size_t)ms * XW, sizeof(uint64_t));
    if (!xc_bits) { fprintf(stderr, "alloc failed\n"); exit(1); }
    long neval = 0; int hist[8] = { 0 }, bestc = c0, n1 = 0, n2 = 0;
    /* phase 1: all single relocations */
    for (int i = 0; i < ms; i++) {
        int old = S[i];
        for (int e = 0; e < NE; e++) {
            if (inS[e]) continue;
            set_S_edge(i, e); build_sder(sd_trial);
            int c = evaluate(sd_trial, F, ev_trial);
            neval++; hist[c < 7 ? c : 7]++;
            if (c < bestc) { bestc = c; printf("  slot %d -> (%d,%d): cost %d\n", i, ep[e], el[e], c); fflush(stdout); }
            XState st = { { { 0, i, e }, { 0, 0, 0 }, { 0, 0, 0 } }, 1, c };
            if (c == 0) xz_report(F, &st, outdir, tag);
            if (c <= c0 + slack1 && n1 < XS_MAX) xl1[n1++] = st;
            set_S_edge(i, old);
        }
        if (i % 10 == 9 || i == ms - 1) {
            printf("exhx phase 1: slot %d done, %ld evaluations, best %d, hist 0..6,7+: %d %d %d %d %d %d %d %d, %d kept, %.0fs\n", i, neval, bestc,
                   hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7], n1, now_sec() - t_start);
            fflush(stdout);
        }
    }
    printf("exhx phase 1 done: %ld evaluations, min cost %d, hist 0..6,7+: %d %d %d %d %d %d %d %d, %d states with cost <= %d, %.0fs\n", neval, bestc,
           hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7], n1, c0 + slack1, now_sec() - t_start);
    fflush(stdout);
    /* phase S: label swaps (and 3-cycles) */
    if (swaps >= 2) {
        int sh[8] = { 0 }; long nsw = 0; int sbest = 1 << 30;
        for (int i = 0; i < ms; i++)
            for (int j = i + 1; j < ms; j++) {
                swap_S(i, j); build_sder(sd_trial);
                int c = evaluate(sd_trial, F, ev_trial);
                nsw++; sh[c < 7 ? c : 7]++;
                if (c < sbest) sbest = c;
                XState st = { { { 1, i, j }, { 0, 0, 0 }, { 0, 0, 0 } }, 1, c };
                if (c == 0) xz_report(F, &st, outdir, tag);
                if (c <= c0 + slack1 && n1 < XS_MAX) xl1[n1++] = st;
                swap_S(i, j);
            }
        printf("exhx swaps: %ld pairwise label swaps, min cost %d, hist 0..6,7+: %d %d %d %d %d %d %d %d, %.0fs\n", nsw, sbest,
               sh[0], sh[1], sh[2], sh[3], sh[4], sh[5], sh[6], sh[7], now_sec() - t_start);
        fflush(stdout);
        if (swaps >= 3) {
            memset(sh, 0, sizeof sh); nsw = 0; sbest = 1 << 30;
            for (int i = 0; i < ms; i++)
                for (int j = i + 1; j < ms; j++)
                    for (int k = i + 1; k < ms; k++) {
                        if (k == j) continue;
                        /* 3-cycle i -> j -> k -> i of the labels: S[i], S[j], S[k] = S[k], S[i], S[j] */
                        swap_S(i, j); swap_S(i, k); build_sder(sd_trial);
                        int c = evaluate(sd_trial, F, ev_trial);
                        nsw++; sh[c < 7 ? c : 7]++;
                        if (c < sbest) sbest = c;
                        XState st = { { { 1, i, j }, { 1, i, k }, { 0, 0, 0 } }, 2, c };
                        if (c == 0) xz_report(F, &st, outdir, tag);
                        swap_S(i, k); swap_S(i, j);
                    }
            printf("exhx 3-cycles: %ld label 3-cycles, min cost %d, hist 0..6,7+: %d %d %d %d %d %d %d %d, %.0fs\n", nsw, sbest,
                   sh[0], sh[1], sh[2], sh[3], sh[4], sh[5], sh[6], sh[7], now_sec() - t_start);
            fflush(stdout);
        }
    }
    /* phase 2: fixing moves from the kept level-1 states, lowest cost first */
    qsort(xl1, n1, sizeof(XState), xs_cmp);
    int nexp1 = n1 < max1 ? n1 : max1;
    long neval2 = 0; int h2[8] = { 0 }, done1 = 0;
    for (int t = 0; t < nexp1; t++) {
        if (now_sec() - t_start > budget) { printf("exhx phase 2: budget reached after %d states\n", t); break; }
        int undo = -1;
        if (!op_apply(&xl1[t].op[0], &undo)) continue;
        build_sder(sd_trial);
        int c1 = evaluate(sd_trial, F, ev_trial);
        long before = neval2; int hb4[8]; memcpy(hb4, h2, sizeof hb4);
        int n2b = n2;
        int best2 = xc_expand(sd_trial, F, &xl1[t], slack2, xl2, &n2, &neval2, h2, budget, t_start, outdir, tag);
        printf("exhx phase 2: state %d/%d (", t + 1, nexp1); op_print(&xl1[t].op[0]);
        printf("; cost %d, %d involved flags, %d candidates): best %d, hist 0..6,7+: %d %d %d %d %d %d %d %d, %d level-2 states kept, %ld evaluations, %.0fs\n",
               c1, xc_ninv, xc_n, best2, h2[0] - hb4[0], h2[1] - hb4[1], h2[2] - hb4[2], h2[3] - hb4[3], h2[4] - hb4[4], h2[5] - hb4[5], h2[6] - hb4[6], h2[7] - hb4[7],
               n2 - n2b, neval2 - before, now_sec() - t_start);
        fflush(stdout);
        op_undo(&xl1[t].op[0], undo);
        done1++;
    }
    printf("exhx phase 2 done: %d of %d kept level-1 states expanded, %ld evaluations, hist 0..6,7+: %d %d %d %d %d %d %d %d, %d level-2 states with cost <= %d, %.0fs\n",
           done1, n1, neval2, h2[0], h2[1], h2[2], h2[3], h2[4], h2[5], h2[6], h2[7], n2, slack2, now_sec() - t_start);
    fflush(stdout);
    /* phase 3: fixing moves from the kept level-2 states */
    qsort(xl2, n2, sizeof(XState), xs_cmp);
    int nexp2 = n2 < max2 ? n2 : max2;
    long neval3 = 0; int h3[8] = { 0 }, done2 = 0;
    for (int t = 0; t < nexp2; t++) {
        if (now_sec() - t_start > budget) { printf("exhx phase 3: budget reached after %d states\n", t); break; }
        int undo[2] = { -1, -1 };
        if (!op_apply(&xl2[t].op[0], &undo[0])) continue;
        if (!op_apply(&xl2[t].op[1], &undo[1])) { op_undo(&xl2[t].op[0], undo[0]); continue; }
        build_sder(sd_trial);
        int c2 = evaluate(sd_trial, F, ev_trial);
        long before = neval3; int hb4[8]; memcpy(hb4, h3, sizeof hb4);
        int best3 = xc_expand(sd_trial, F, &xl2[t], -1, NULL, NULL, &neval3, h3, budget, t_start, outdir, tag);
        printf("exhx phase 3: state %d/%d (", t + 1, nexp2); op_print(&xl2[t].op[0]); printf(", "); op_print(&xl2[t].op[1]);
        printf("; cost %d, %d involved flags, %d candidates): best %d, hist 0..6,7+: %d %d %d %d %d %d %d %d, %ld evaluations, %.0fs\n",
               c2, xc_ninv, xc_n, best3, h3[0] - hb4[0], h3[1] - hb4[1], h3[2] - hb4[2], h3[3] - hb4[3], h3[4] - hb4[4], h3[5] - hb4[5], h3[6] - hb4[6], h3[7] - hb4[7],
               neval3 - before, now_sec() - t_start);
        fflush(stdout);
        op_undo(&xl2[t].op[1], undo[1]); op_undo(&xl2[t].op[0], undo[0]);
        done2++;
    }
    printf("exhx FINAL: phase 1 %ld evaluations (min %d); phase 2 %d/%d states, %ld evaluations; phase 3 %d/%d states, %ld evaluations; cost-0 states %ld; %.0fs\n",
           neval, bestc, done1, n1, neval2, done2, n2, neval3, xz_count, now_sec() - t_start);
    build_sder(sd_trial);
    int cend = evaluate(sd_trial, F, ev_trial);
    printf("exhx: base state restored, cost %d (expected %d)%s\n", cend, c0, cend == c0 ? "" : "  <-- APPLY/UNDO BUG");
    fflush(stdout);
}

/* ---------- main ---------- */
int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s H4_incidence.edges m h seed seconds [opts]\n", argv[0]); return 1; }
    const char *hpath = argv[1];
    ms = atoi(argv[2]); hb = atoi(argv[3]);
    uint64_t seed = strtoull(argv[4], 0, 10);
    double seconds = atof(argv[5]);
    double T0 = 2.0, T1 = 0.15, pS = 0.3, pR = 0.15, pT = 0.5, pSfix = 0.3, optF = 0.0;
    int cycles = 1, do_eval = 0, keepgoing = 0, restart_best = 0, nvalidate = 0, bench = 0, spread = 0, exh1 = 0;
    int endgame_th = 12, kS_end = 16, exh2_slack = -1, exh2_max = 40;
    int exhx_slack = -1, exhx_slack2 = 1, exhx_max1 = 100, exhx_max2 = 300, exhx_swaps = 2;
    const char *fixS = 0, *fixF = 0, *initS = 0, *initF = 0, *tag = "run", *outdir = "."; const char *freezeF = 0;
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "--T0")) T0 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--apts")) napt = atoi(argv[++i]);     /* number of apartment groups (slots 0..4K-1) */
        else if (!strcmp(argv[i], "--topT")) topT = atoi(argv[++i]);     /* edges per slot in the gadget enumeration */
        else if (!strcmp(argv[i], "--pLocal")) pLocal = atof(argv[++i]); /* share of local (phase/orientation) group moves */
        else if (!strcmp(argv[i], "--unlock")) unlock_t = atof(argv[++i]); /* seconds after which apartment groups become free flags */
        else if (!strcmp(argv[i], "--exh1")) exh1 = 1;                   /* exhaustive single-flag move neighbourhood of the initial state */
        else if (!strcmp(argv[i], "--exh2")) exh2_slack = atoi(argv[++i]); /* + second move of an involved flag from every 1-move state with cost <= c0 + slack */
        else if (!strcmp(argv[i], "--exh2max")) exh2_max = atoi(argv[++i]); /* cap on the number of 1-move states expanded */
        else if (!strcmp(argv[i], "--exhx")) exhx_slack = atoi(argv[++i]);   /* exact 1/2/3-move neighbourhoods: level-1 states with cost <= c0 + slack are expanded */
        else if (!strcmp(argv[i], "--exhx2")) exhx_slack2 = atoi(argv[++i]); /* level-2 states with cost <= this are expanded (phase 3) */
        else if (!strcmp(argv[i], "--exhxmax")) exhx_max1 = atoi(argv[++i]); /* cap on expanded level-1 states */
        else if (!strcmp(argv[i], "--exhx3max")) exhx_max2 = atoi(argv[++i]); /* cap on expanded level-2 states */
        else if (!strcmp(argv[i], "--swaps")) exhx_swaps = atoi(argv[++i]);   /* 0: no label swaps, 2: pairwise, 3: also 3-cycles */
        else if (!strcmp(argv[i], "--T1")) T1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--cycles")) cycles = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--pS")) pS = atof(argv[++i]);
        else if (!strcmp(argv[i], "--pR")) pR = atof(argv[++i]);
        else if (!strcmp(argv[i], "--pT")) pT = atof(argv[++i]);       /* fraction of targeted moves */
        else if (!strcmp(argv[i], "--pSfix")) pSfix = atof(argv[++i]); /* targeted moves that move an S edge */
        else if (!strcmp(argv[i], "--kF")) kF = atoi(argv[++i]);       /* candidates per targeted F move */
        else if (!strcmp(argv[i], "--kS")) kS = atoi(argv[++i]);       /* candidates per targeted S move */
        else if (!strcmp(argv[i], "--mc")) mc_mode = 1;                 /* min-conflicts S relocation */
        else if (!strcmp(argv[i], "--optF")) optF = atof(argv[++i]);   /* seconds of F pre-optimisation */
        else if (!strcmp(argv[i], "--optFper")) optF_per = atof(argv[++i]); /* seconds per F restart */
        else if (!strcmp(argv[i], "--optFkeep")) optF_keep = 1;       /* F restarts from the initial (seeded) F */
        else if (!strcmp(argv[i], "--w5")) w5 = atoi(argv[++i]);       /* weight of distance-5 pairs in the F objective */
        else if (!strcmp(argv[i], "--tabu")) tabu_len = atol(argv[++i]);
        else if (!strcmp(argv[i], "--endgame")) endgame_th = atoi(argv[++i]); /* cost threshold for intensification */
        else if (!strcmp(argv[i], "--kSend")) kS_end = atoi(argv[++i]);       /* relocation candidates in the endgame */
        else if (!strcmp(argv[i], "--freeze")) freezeF = argv[++i];
        else if (!strcmp(argv[i], "--fixS")) fixS = argv[++i];
        else if (!strcmp(argv[i], "--fixF")) fixF = argv[++i];
        else if (!strcmp(argv[i], "--initS")) initS = argv[++i];
        else if (!strcmp(argv[i], "--initF")) initF = argv[++i];
        else if (!strcmp(argv[i], "--eval")) do_eval = 1;
        else if (!strcmp(argv[i], "--tag")) tag = argv[++i];
        else if (!strcmp(argv[i], "--outdir")) outdir = argv[++i];
        else if (!strcmp(argv[i], "--keepgoing")) keepgoing = 1;
        else if (!strcmp(argv[i], "--restart-best")) restart_best = 1;
        else if (!strcmp(argv[i], "--spread")) spread = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--validate")) nvalidate = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--bench")) bench = 1;
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (fixS) initS = fixS;
    if (fixF) initF = fixF;
    nnew = ms + hb; N = NOLD + nnew; nn32 = (nnew + 31) / 32 * 32;
    if (nnew > MAXNEW || ms < 1 || nn32 > MAXNEW) { fprintf(stderr, "bad m/h\n"); return 1; }
    for (int i = 0; i < nnew; i++) tdeg[i] = i < ms ? DEG - 2 : DEG;
    { int sum = 0; for (int i = 0; i < nnew; i++) sum += tdeg[i]; if (sum % 2) { fprintf(stderr, "odd F degree sum: m must be even (when h is even)\n"); return 1; } }
    rng_s = seed * 0x9E3779B97F4A7C15ULL + 12345;
    double t_start = now_sec();
    load_H(hpath);
    precompute();
    evA.viol_oo = malloc(sizeof(int) * MAXBOO);
    evB.viol_oo = malloc(sizeof(int) * MAXBOO);
    if (!evA.viol_oo || !evB.viol_oo) { fprintf(stderr, "alloc failed\n"); return 1; }
    fprintf(stderr, "precompute done in %.2fs\n", now_sec() - t_start);
    memset(inS, 0, sizeof inS);
    for (int i = 0; i < ms; i++) S[i] = -1;
    memset(frozen, 0, sizeof frozen);
    if (napt < 0 || 4 * napt > ms) { fprintf(stderr, "bad --apts %d for m = %d\n", napt, ms); return 1; }
    if (napt > 0) build_apartment_tables();

    if (nvalidate > 0) { validate(nvalidate, spread); return 0; }

    static SDer sdA, sdB;
    SDer *sd_cur = &sdA, *sd_trial = &sdB;
    Fgraph F, bestF;
    int newS[MAXNEW], bestS[MAXNEW];

    if (initS) {
        int A[MAXNEW], B[MAXNEW];
        int n = read_pairs(initS, A, B, MAXNEW);
        if (n > ms || (fixS && n != ms)) { fprintf(stderr, "S file has %d edges, m = %d\n", n, ms); return 1; }
        static uint8_t used[NE];
        memset(used, 0, sizeof used);
        for (int i = 0; i < n; i++) {
            int e = find_edge(A[i], B[i]);
            if (e < 0) { fprintf(stderr, "(%d,%d) is not an edge of H\n", A[i], B[i]); return 1; }
            if (used[e]) { fprintf(stderr, "duplicate S edge\n"); return 1; }
            newS[i] = e; used[e] = 1;
        }
        if (napt > 0) {
            /* the file holds the first n slots (whole groups); random gadgets / flags fill the rest */
            if (n < 4 * napt && n % 4) { fprintf(stderr, "initS: %d flags do not fill whole apartment groups\n", n); return 1; }
            random_S_apts(newS, spread, n);
        } else
            for (int i = n; i < ms; i++) {
                int e;
                do { e = rndint(NE); } while (used[e]);
                used[e] = 1; newS[i] = e;
            }
        printf("S seeded with %d edges from %s, %d random\n", n, initS, ms - n);
    } else if (napt > 0) random_S_apts(newS, spread, 0);
    else random_S(newS, spread);
    /* install: S currently empty */
    install_S(newS);
    for (int k = 0; k < napt; k++)
        if (!check_group(k)) { fprintf(stderr, "slots %d..%d do not form an apartment gadget (a,b,c,d = positions 0,3,6,9)\n", 4 * k, 4 * k + 3); return 1; }
    if (napt > 0) printf("%d apartment groups verified (gadget flags at positions 0,3,6,9; chords (4k,4k+2),(4k+1,4k+3) frozen)\n", napt);

    if (initF) {
        int A[MAXF], B[MAXF];
        int n = read_pairs(initF, A, B, MAXF);
        for (int i = 0; i < nnew; i++) F.fdeg[i] = 0;
        F.nF = 0;
        for (int i = 0; i < n; i++) {
            if (A[i] < 0 || B[i] < 0 || A[i] >= nnew || B[i] >= nnew || A[i] == B[i] || f_has(&F, A[i], B[i])) { fprintf(stderr, "bad F edge (%d,%d)\n", A[i], B[i]); return 1; }
            f_add(&F, A[i], B[i]);
        }
        for (int k = 0; k < napt; k++) {
            /* the chords of every apartment group must be present; they are frozen */
            if (!f_has(&F, 4 * k, 4 * k + 2) || !f_has(&F, 4 * k + 1, 4 * k + 3)) { fprintf(stderr, "initF: chords of group %d missing\n", k); return 1; }
            frozen[4 * k][4 * k + 2] = frozen[4 * k + 2][4 * k] = 1;
            frozen[4 * k + 1][4 * k + 3] = frozen[4 * k + 3][4 * k + 1] = 1;
        }
        for (int i = 0; i < nnew; i++) {
            if (F.fdeg[i] == tdeg[i]) continue;
            if (F.fdeg[i] != 0) { fprintf(stderr, "initF: new vertex %d has F-degree %d, want 0 or %d\n", i, F.fdeg[i], tdeg[i]); return 1; }
            /* insert into random F edges: (a,b) -> (a,i),(i,b) gives 2; a third edge needs a pair of such plus one more edge... */
            int need = tdeg[i];
            for (int tries = 0; need >= 2; tries++) {
                if (tries > 10000 || F.nF == 0) { fprintf(stderr, "initF: cannot insert vertex %d\n", i); return 1; }
                int t = rndint(F.nF);
                int a = F.fea[t], b = F.feb[t];
                if (a == i || b == i || f_has(&F, a, i) || f_has(&F, i, b) || frozen[a][b]) continue;
                f_del(&F, a, b); f_add(&F, a, i); f_add(&F, i, b);
                need -= 2;
            }
            if (need == 1) {
                /* odd remainder: pair up with another deficient vertex */
                int j = -1;
                for (int jj = i + 1; jj < nnew; jj++) if (F.fdeg[jj] == 0 || F.fdeg[jj] == tdeg[jj] - 1) { j = jj; break; }
                if (j < 0) { fprintf(stderr, "initF: cannot complete degrees\n"); return 1; }
                if (F.fdeg[j] == 0) {
                    int need2 = tdeg[j] - 1;
                    for (int tries = 0; need2 >= 2; tries++) {
                        if (tries > 10000) { fprintf(stderr, "initF: cannot insert vertex %d\n", j); return 1; }
                        int t = rndint(F.nF);
                        int a = F.fea[t], b = F.feb[t];
                        if (a == j || b == j || a == i || b == i || f_has(&F, a, j) || f_has(&F, j, b) || frozen[a][b]) continue;
                        f_del(&F, a, b); f_add(&F, a, j); f_add(&F, j, b);
                        need2 -= 2;
                    }
                }
                if (f_has(&F, i, j)) { fprintf(stderr, "initF: completion failed\n"); return 1; }
                f_add(&F, i, j);
            }
            printf("F: inserted new vertex %d\n", i);
        }
        for (int i = 0; i < nnew; i++) if (F.fdeg[i] != tdeg[i]) { fprintf(stderr, "F: new vertex %d has F-degree %d, want %d\n", i, F.fdeg[i], tdeg[i]); return 1; }
        if (freezeF) {
            int FA[MAXF], FB[MAXF];
            int nfz = read_pairs(freezeF, FA, FB, MAXF);
            for (int i = 0; i < nfz; i++) {
                if (!f_has(&F, FA[i], FB[i])) { fprintf(stderr, "frozen edge (%d,%d) not in F\n", FA[i], FB[i]); return 1; }
                frozen[FA[i]][FB[i]] = frozen[FB[i]][FA[i]] = 1;
            }
            printf("%d F edges frozen\n", nfz);
        }
    } else random_F(&F);
    if (optF > 0 && !fixF) optimise_F(&F, optF);
    compute_DF(&F);
    tm_sd_trial = sd_trial;

    build_sder(sd_cur);
    int cur_cost = evaluate(sd_cur, &F, ev_cur);
    printf("init cost %d (new-old %d, new-new %d, old-old %d) nBoo %d\n", cur_cost, ev_cur->no, ev_cur->nn, ev_cur->oo, nBoo);
    fflush(stdout);
    if (exhx_slack >= 0) {
        exhx_run(sd_trial, &F, cur_cost, exhx_slack, exhx_slack2, exhx_max1, exhx_max2, exhx_swaps, seconds, t_start, outdir, tag);
        return 0;
    }
    if (exh2_slack >= 0) {
        /* phase 1: every single-flag move of the initial state (as --exh1), collecting the states with cost
           <= c0 + slack; phase 2: from each collected state, every move of every flag involved in one of
           its remaining violations (the flag itself, the S edges blocking the geodesic, both ends of a
           new-new pair, the S edges on the geodesics of an old-old pair).  Cost-0 states are checked by
           all-pairs BFS and written. */
        int c0 = cur_cost, bestc = cur_cost, hist[8] = { 0 };
        long neval = 0, nzero = 0;
        typedef struct { int slot, e, c; } Mv;
        static Mv low[65536]; int nlow = 0;
        for (int i = 0; i < ms; i++) {
            int old = S[i];
            for (int e = 0; e < NE; e++) {
                if (inS[e]) continue;
                set_S_edge(i, e); build_sder(sd_trial);
                int c = evaluate(sd_trial, &F, ev_trial);
                neval++;
                hist[c < 7 ? c : 7]++;
                if (c < bestc) { bestc = c; printf("  slot %d -> edge (%d,%d): cost %d\n", i, ep[e], el[e], c); fflush(stdout); }
                if (c == 0) {
                    int mx, mn;
                    int ex = exact_cost(&F, &mx, &mn);
                    printf("ZERO COST: exact all-pairs check pairs>6 = %d, maxdeg %d, mindeg %d, n = %d\n", ex, mx, mn, N);
                    char path[1024];
                    snprintf(path, sizeof path, "%s/%s_exh1_n%d_slot%d_edge%d.edges", outdir, tag, N, i, e);
                    write_edges(&F, path);
                    printf("wrote %s\n", path);
                    print_state(&F, stdout);
                    fflush(stdout);
                    nzero++;
                }
                if (c <= c0 + exh2_slack && nlow < 65536) { low[nlow].slot = i; low[nlow].e = e; low[nlow].c = c; nlow++; }
                set_S_edge(i, old);
            }
            printf("exh2 phase 1: slot %d done, %ld evaluations, best %d, hist 0..6,7+: %d %d %d %d %d %d %d %d, %d low states, %.0fs\n", i, neval, bestc,
                   hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7], nlow, now_sec() - t_start);
            fflush(stdout);
        }
        printf("exh2 phase 1 done: %ld evaluations, min cost %d, cost-0 states %ld, %d states with cost <= %d\n", neval, bestc, nzero, nlow, c0 + exh2_slack);
        fflush(stdout);
        /* phase 2, lowest-cost states first */
        for (int a = 0; a < nlow; a++) for (int b = a + 1; b < nlow; b++) if (low[b].c < low[a].c) { Mv t = low[a]; low[a] = low[b]; low[b] = t; }
        int nexp = nlow < exh2_max ? nlow : exh2_max;
        long neval2 = 0;
        for (int t = 0; t < nexp; t++) {
            int i = low[t].slot, old = S[i];
            set_S_edge(i, low[t].e); build_sder(sd_trial);
            int c1 = evaluate(sd_trial, &F, ev_trial);
            /* involved slots */
            int inv[MAXNEW], ninv = 0; static uint8_t mark[MAXNEW]; memset(mark, 0, sizeof mark);
            for (int k = 0; k < nnew; k++) {
                if (ev_trial->nocnt[k] == 0) continue;
                if (!mark[k] && k < ms) { mark[k] = 1; inv[ninv++] = k; }
                for (int w = 0; w <= OLDWORDS; w++) {
                    uint64_t word = ev_trial->uncov[k][w];
                    while (word) {
                        int v = w * 64 + __builtin_ctzll(word); word &= word - 1;
                        if (k >= ms) continue;
                        int ends[2] = { ep[S[k]], el[S[k]] };
                        for (int s2 = 0; s2 < 2; s2++) {
                            int sc[8], nsc = 0;
                            if (distH[ends[s2]][v] <= 5) nsc = s_edges_on_geodesic(ends[s2], v, sc);
                            for (int q = 0; q < nsc; q++) { int kk = sidx[sc[q]]; if (kk >= 0 && !mark[kk]) { mark[kk] = 1; inv[ninv++] = kk; } }
                        }
                    }
                }
            }
            for (int q = 0; q < ev_trial->nn; q++) {
                int ka = ev_trial->viol_nn_a[q], kb = ev_trial->viol_nn_b[q];
                if (ka < ms && !mark[ka]) { mark[ka] = 1; inv[ninv++] = ka; }
                if (kb < ms && !mark[kb]) { mark[kb] = 1; inv[ninv++] = kb; }
            }
            for (int q = 0; q < ev_trial->oo; q++) {
                int pid = ev_trial->viol_oo[q], a = pid / NOLD, b = pid % NOLD, sc[64], nsc = 0;
                if (distH[a][b] <= 5) nsc = s_edges_on_geodesic(a, b, sc);
                else for (int g = 0; g < DEG; g++) { if (inS[edge_of[a][g]]) sc[nsc++] = edge_of[a][g]; nsc += s_edges_on_geodesic(adjH[a][g], b, sc + nsc); }
                for (int x = 0; x < nsc; x++) { int kk = sidx[sc[x]]; if (kk >= 0 && !mark[kk]) { mark[kk] = 1; inv[ninv++] = kk; } }
            }
            int best2 = c1;
            for (int x = 0; x < ninv; x++) {
                int s2 = inv[x], old2 = S[s2];
                for (int e2 = 0; e2 < NE; e2++) {
                    if (inS[e2]) continue;
                    set_S_edge(s2, e2); build_sder(sd_trial);
                    int c = evaluate(sd_trial, &F, ev_trial);
                    neval2++;
                    if (c < best2) best2 = c;
                    if (c == 0) {
                        int mx, mn;
                        int ex = exact_cost(&F, &mx, &mn);
                        printf("ZERO COST (2 moves: slot %d -> (%d,%d), slot %d -> (%d,%d)): exact all-pairs check pairs>6 = %d, maxdeg %d, mindeg %d, n = %d\n",
                               i, ep[low[t].e], el[low[t].e], s2, ep[e2], el[e2], ex, mx, mn, N);
                        char path[1024];
                        snprintf(path, sizeof path, "%s/%s_exh2_n%d_%d_%d.edges", outdir, tag, N, t, e2);
                        write_edges(&F, path);
                        printf("wrote %s\n", path);
                        print_state(&F, stdout);
                        fflush(stdout);
                        nzero++;
                    }
                    set_S_edge(s2, old2);
                }
            }
            printf("exh2 phase 2: state %d/%d (slot %d -> (%d,%d), cost %d, %d involved slots): best 2-move cost %d, %ld evaluations, %.0fs\n",
                   t + 1, nexp, i, ep[low[t].e], el[low[t].e], c1, ninv, best2, neval2, now_sec() - t_start);
            fflush(stdout);
            set_S_edge(i, old);
        }
        printf("exh2 FINAL: phase 1 %ld evaluations (min %d), phase 2 %d states, %ld evaluations, cost-0 states %ld\n", neval, bestc, nexp, neval2, nzero);
        return 0;
    }
    if (exh1) {
        /* every single-flag move (slot i -> edge e) of the initial state, exact cost; cost-0 states are
           checked by all-pairs BFS and written */
        int bestc = cur_cost, hist[8] = { 0 };
        long neval = 0, nzero = 0;
        for (int i = 0; i < ms; i++) {
            int old = S[i];
            for (int e = 0; e < NE; e++) {
                if (inS[e]) continue;
                set_S_edge(i, e); build_sder(sd_trial);
                int c = evaluate(sd_trial, &F, ev_trial);
                neval++;
                hist[c < 7 ? c : 7]++;
                if (c < bestc) { bestc = c; printf("  slot %d -> edge (%d,%d): cost %d\n", i, ep[e], el[e], c); fflush(stdout); }
                if (c == 0) {
                    int mx, mn;
                    int ex = exact_cost(&F, &mx, &mn);
                    printf("ZERO COST: exact all-pairs check pairs>6 = %d, maxdeg %d, mindeg %d, n = %d\n", ex, mx, mn, N);
                    char path[1024];
                    snprintf(path, sizeof path, "%s/%s_exh1_n%d_slot%d_edge%d.edges", outdir, tag, N, i, e);
                    write_edges(&F, path);
                    printf("wrote %s\n", path);
                    print_state(&F, stdout);
                    fflush(stdout);
                    nzero++;
                }
                set_S_edge(i, old);
            }
            printf("exh1: slot %d done, %ld evaluations, best %d, cost histogram 0..6,7+: %d %d %d %d %d %d %d %d, %.0fs\n", i, neval, bestc,
                   hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], hist[6], hist[7], now_sec() - t_start);
            fflush(stdout);
        }
        printf("exh1 FINAL: %ld evaluations, min cost over the 1-move neighbourhood %d, cost-0 states %ld\n", neval, bestc, nzero);
        return 0;
    }
    if (do_eval) {
        int mx, mn;
        int ex = exact_cost(&F, &mx, &mn);
        int bno, bnn, boo_;
        int cb = evaluate_bfs(&F, &bno, &bnn, &boo_);
        printf("bfs-method cost %d (no %d nn %d oo %d); exact all-pairs BFS: pairs > 6 = %d, max degree %d, min degree %d, n = %d\n", cb, bno, bnn, boo_, ex, mx, mn, N);
        print_state(&F, stdout);
        /* details of the violations */
        compute_DF(&F);
        int diamF = 0, n6 = 0, n7 = 0;
        for (int i = 0; i < nnew; i++) for (int j = i + 1; j < nnew; j++) { if (DF[i][j] > diamF) diamF = DF[i][j]; n6 += DF[i][j] == 6; n7 += DF[i][j] >= 7; }
        printf("F: diameter %d, pairs at F-distance 6: %d, >= 7: %d\n", diamF, n6, n7);
        for (int t = 0; t < ev_cur->nn; t++) {
            int i = ev_cur->viol_nn_a[t], j = ev_cur->viol_nn_b[t];
            printf("  new-new violation (%d,%d): DF=%d edge-distance=%d\n", i, j, DF[i][j], (i < ms && j < ms) ? edge_dist(S[i], S[j]) : -1);
        }
        for (int i = 0; i < nnew; i++) {
            for (int w = 0; w <= OLDWORDS; w++) {
                uint64_t word = ev_cur->uncov[i][w];
                while (word) {
                    int v = w * 64 + __builtin_ctzll(word);
                    word &= word - 1;
                    printf("  new-old violation (x%d, %d):", i, v);
                    if (i < ms) {
                        int ends[2] = { ep[S[i]], el[S[i]] };
                        for (int s = 0; s < 2; s++) {
                            int sc[8], nsc = 0;
                            if (distH[ends[s]][v] <= 5) nsc = s_edges_on_geodesic(ends[s], v, sc);
                            printf(" d_H(%s,v)=%d blocked by", s ? "l" : "p", distH[ends[s]][v]);
                            for (int q = 0; q < nsc; q++) printf(" x%d(DF %d,t %d)", sidx[sc[q]], DF[i][sidx[sc[q]]], edge_dist(S[i], sc[q]));
                            printf(";");
                        }
                    }
                    printf("\n");
                }
            }
        }
        for (int t = 0; t < ev_cur->oo; t++) {
            int pid = ev_cur->viol_oo[t];
            printf("  old-old violation (%d,%d): d_H=%d\n", pid / NOLD, pid % NOLD, distH[pid / NOLD][pid % NOLD]);
        }
        /* edge-distance / F-distance profile of S */
        int prof[6][10];
        memset(prof, 0, sizeof prof);
        for (int i = 0; i < ms; i++) for (int j = i + 1; j < ms; j++) { int td = edge_dist(S[i], S[j]); if (td > 5) td = 5; int df = DF[i][j] > 9 ? 9 : DF[i][j]; prof[td][df]++; }
        printf("pairs of S edges by (edge-distance t, F-distance): rows t=0..5(+), cols DF=0..9\n");
        for (int td = 0; td <= 5; td++) { printf("  t=%d:", td); for (int df = 0; df <= 9; df++) printf(" %4d", prof[td][df]); printf("\n"); }
        return 0;
    }
    if (bench) {
        double t0 = now_sec(); int cnt = 0;
        for (int k = 0; k < 2000; k++) { if (f_move(&F, pR)) { evaluate(sd_cur, &F, ev_trial); f_undo(&F); cnt++; } }
        double t1 = now_sec();
        printf("bench: %d F moves+eval in %.3fs = %.0f/s\n", cnt, t1 - t0, cnt / (t1 - t0));
        t0 = now_sec(); cnt = 0;
        for (int k = 0; k < 200; k++) {
            int i = rndint(ms), e = rndint(NE);
            if (inS[e]) continue;
            int old = S[i];
            set_S_edge(i, e); build_sder(sd_trial); evaluate(sd_trial, &F, ev_trial);
            set_S_edge(i, old); cnt++;
        }
        t1 = now_sec();
        printf("bench: %d S moves (incl. undo) in %.3fs = %.0f/s\n", cnt, t1 - t0, cnt / (t1 - t0));
        t0 = now_sec();
        for (int k = 0; k < 50; k++) build_sder(sd_trial);
        t1 = now_sec();
        printf("bench: build_sder %.3f ms\n", (t1 - t0) / 50 * 1000);
        t0 = now_sec();
        for (int k = 0; k < 200; k++) evaluate(sd_cur, &F, ev_trial);
        t1 = now_sec();
        printf("bench: evaluate %.3f ms\n", (t1 - t0) / 200 * 1000);
        t0 = now_sec(); cnt = 0; int nT = 0;
        for (int k = 0; k < 2000; k++) {
            int si, so;
            int kind = targeted_move(sd_cur, sd_trial, &F, ev_cur, pSfix, &si, &so);
            if (kind == 0) continue;
            nT++;
            if (kind == 1) f_undo(&F);
            else s_undo();
        }
        t1 = now_sec();
        printf("bench: %d targeted moves (of 2000 proposals, kF %d kS %d) in %.3fs = %.0f/s\n", nT, kF, kS, t1 - t0, nT / (t1 - t0));
        /* delta histograms of proposed moves from the current state */
        int hist[4][16];
        memset(hist, 0, sizeof hist);
        for (int kind = 0; kind < 4; kind++) {
            for (int k = 0; k < 2000; k++) {
                int c = -1, si, so;
                if (kind == 0) { if (!f_move(&F, 0.0)) continue; c = evaluate(sd_cur, &F, ev_trial); f_undo(&F); }
                else if (kind == 1) { if (!f_move(&F, 1.0)) continue; c = evaluate(sd_cur, &F, ev_trial); f_undo(&F); }
                else if (kind == 2) {
                    int i = rndint(ms), e = random_target_edge(S[i]);
                    if (inS[e] || i < 4 * napt) continue;
                    int old = S[i];
                    set_S_edge(i, e); build_sder(sd_trial); c = evaluate(sd_trial, &F, ev_trial); set_S_edge(i, old);
                } else {
                    int kk = targeted_move(sd_cur, sd_trial, &F, ev_cur, pSfix, &si, &so);
                    if (kk == 0) continue;
                    c = tm_cost;
                    if (kk == 1) f_undo(&F); else s_undo();
                }
                int d = c - cur_cost;
                int bin = d < 0 ? 0 : d == 0 ? 1 : d <= 2 ? 2 : d <= 5 ? 3 : d <= 10 ? 4 : d <= 20 ? 5 : d <= 50 ? 6 : 7;
                hist[kind][bin]++;
            }
            const char *names[4] = { "2-swap", "relocate", "S-move", "targeted" };
            printf("delta histogram %-9s: <0:%d =0:%d 1-2:%d 3-5:%d 6-10:%d 11-20:%d 21-50:%d >50:%d\n", names[kind],
                   hist[kind][0], hist[kind][1], hist[kind][2], hist[kind][3], hist[kind][4], hist[kind][5], hist[kind][6], hist[kind][7]);
        }
        return 0;
    }

    memcpy(bestS, S, sizeof(int) * ms); bestF = F; int best_cost = cur_cost;
    int best_no = ev_cur->no, best_nn = ev_cur->nn, best_oo = ev_cur->oo;
    long it = 0, acc = 0, nS = 0, nFm = 0, nT = 0, accT = 0, aNeg = 0, aZero = 0, aPos = 0;
    int kS_base = kS;
    double t_last = now_sec();
    double cyc_len = seconds / cycles;
    int done = 0, last_cyc = 0;
    while (!done) {
        double el_t = now_sec() - t_start;
        if (el_t >= seconds) break;
        int cyc = (int)(el_t / cyc_len);
        double frac = (el_t - cyc * cyc_len) / cyc_len;
        double T = T0 * pow(T1 / T0, frac);
        if (cyc != last_cyc) {
            last_cyc = cyc;
            if (restart_best) {
                install_S(bestS); F = bestF; compute_DF(&F); build_sder(sd_cur); cur_cost = evaluate(sd_cur, &F, ev_cur);
                printf("cycle %d: restart from best (%d)\n", cyc, cur_cost);
            }
        }
        it++;
        if (unlock_t >= 0 && napt > 0 && el_t >= unlock_t) {
            printf("t=%.1f unlocking %d apartment groups: all %d flags now move individually (chords stay in F)\n", el_t, napt, ms);
            fflush(stdout);
            napt = 0;
        }
        double r = rnd01();
        int trial_cost;
        tm_sd_trial = sd_trial;
        if (cur_cost <= endgame_th) { exhaustive_k = 1; kS = kS_end; } else { exhaustive_k = 0; kS = kS_base; }
        if (r < pT && cur_cost > 0) {
            /* targeted repair of a current violation */
            int si = -1, so = -1;
            int kind = targeted_move(sd_cur, sd_trial, &F, ev_cur, fixS ? 0.0 : pSfix, &si, &so);
            if (kind == 0) continue;
            trial_cost = tm_cost;   /* the applied candidate was evaluated into ev_trial (and sd_trial if kind 2) */
            if (kind == 2) nS++; else nFm++;
            nT++;
            int delta = trial_cost - cur_cost;
            if (delta <= 0 || rnd01() < exp(-delta / T)) {
                if (delta < 0) aNeg++; else if (delta == 0) aZero++; else aPos++;
                if (kind == 2) { SDer *tmp = sd_cur; sd_cur = sd_trial; sd_trial = tmp; n_smoves++; for (int t = 0; t < su_n; t++) tabu_until[su_i[t]] = n_smoves + tabu_len; if (su_n == 4) nG++; }
                else compute_DF(&F);
                Eval *te = ev_cur; ev_cur = ev_trial; ev_trial = te;
                cur_cost = trial_cost; acc++; accT++;
            } else {
                if (kind == 2) s_undo(); else f_undo(&F);
                continue;
            }
        } else if (r < pT + (1.0 - pT) * pS && !fixS) {
            int i = rndint(ms), old = S[i];
            if (i < 4 * napt) {
                /* apartment group move: local (phase / orientation of the same apartment) or global */
                if (group_tabu(i / 4)) continue;
                int res = rnd01() < pLocal ? group_local_move(&F, i / 4) : group_relocate_mc(&F, i / 4, -1, 0, 0);
                if (!res) continue;
                trial_cost = tm_cost;
            } else if (mc_mode) {
                int si, so;
                if (tabu_len > 0 && tabu_until[i] > n_smoves) continue;
                if (!s_relocate_mc(&F, i, &si, &so)) continue;
                trial_cost = tm_cost;
            } else {
                int e;
                if (rnd01() < 0.5) e = near2[offnear[old] + rndint(offnear[old + 1] - offnear[old])];
                else e = rndint(NE);
                if (inS[e]) continue;
                set_S_edge(i, e);
                build_sder(sd_trial);
                trial_cost = evaluate(sd_trial, &F, ev_trial);
                su_n = 1; su_i[0] = i; su_e[0] = old;
            }
            nS++;
            int delta = trial_cost - cur_cost;
            if (delta <= 0 || rnd01() < exp(-delta / T)) {
                if (delta < 0) aNeg++; else if (delta == 0) aZero++; else aPos++;
                SDer *tmp = sd_cur; sd_cur = sd_trial; sd_trial = tmp;
                Eval *te = ev_cur; ev_cur = ev_trial; ev_trial = te;
                cur_cost = trial_cost; acc++;
                n_smoves++; for (int t = 0; t < su_n; t++) tabu_until[su_i[t]] = n_smoves + tabu_len;
                if (su_n == 4) nG++;
            } else {
                s_undo();
                continue;
            }
        } else if (!fixF) {
            if (!f_move(&F, pR / (1.0 - pS))) continue;
            trial_cost = evaluate(sd_cur, &F, ev_trial);
            nFm++;
            int delta = trial_cost - cur_cost;
            if (delta <= 0 || rnd01() < exp(-delta / T)) {
                Eval *te = ev_cur; ev_cur = ev_trial; ev_trial = te;
                cur_cost = trial_cost; acc++;
                compute_DF(&F);
            } else {
                f_undo(&F);
                continue;
            }
        } else continue;
        if (cur_cost < best_cost) {
            best_cost = cur_cost; best_no = ev_cur->no; best_nn = ev_cur->nn; best_oo = ev_cur->oo;
            memcpy(bestS, S, sizeof(int) * ms); bestF = F;
            printf("t=%.1f it=%ld cyc=%d T=%.3f best=%d (no %d nn %d oo %d) nBoo=%d\n", now_sec() - t_start, it, cyc, T,
                   best_cost, best_no, best_nn, best_oo, nBoo);
            if (best_cost <= 6) print_state(&F, stdout);
            fflush(stdout);
            if (best_cost == 0) {
                int mx, mn;
                int ex = exact_cost(&F, &mx, &mn);
                printf("ZERO COST: exact all-pairs check pairs>6 = %d, maxdeg %d, mindeg %d, n = %d\n", ex, mx, mn, N);
                char path[1024];
                snprintf(path, sizeof path, "%s/%s_n%d_seed%llu_it%ld.edges", outdir, tag, N, (unsigned long long)seed, it);
                write_edges(&F, path);
                printf("wrote %s\n", path);
                fflush(stdout);
                if (!keepgoing) { done = 1; break; }
            }
        }
        if (now_sec() - t_last > 30) {
            t_last = now_sec();
            printf("t=%.0f it=%ld acc=%ld (-%ld =%ld +%ld) S=%ld G=%ld F=%ld T=%ld accT=%ld cyc=%d T=%.3f cur=%d best=%d nBoo=%d\n", t_last - t_start, it, acc, aNeg, aZero, aPos, nS, nG, nFm, nT, accT, cyc, T, cur_cost, best_cost, nBoo);
            fflush(stdout);
        }
    }
    install_S(bestS); F = bestF;
    build_sder(sd_cur);
    int bc = evaluate(sd_cur, &F, ev_cur);
    int mx, mn;
    int ex = exact_cost(&F, &mx, &mn);
    printf("FINAL best=%d (no %d nn %d oo %d) exact=%d iterations=%ld accepted=%ld time=%.1f n=%d maxdeg=%d mindeg=%d\n",
           bc, ev_cur->no, ev_cur->nn, ev_cur->oo, ex, it, acc, now_sec() - t_start, N, mx, mn);
    print_state(&F, stdout);
    char path[1024];
    snprintf(path, sizeof path, "%s/%s_best.txt", outdir, tag);
    FILE *f = fopen(path, "w");
    if (f) { fprintf(f, "# cost %d exact %d n %d\n", bc, ex, N); print_state(&F, f); fclose(f); }
    /* seed files in the "p l" / "i j" format accepted by --initS/--initF */
    snprintf(path, sizeof path, "%s/%s_best_S.txt", outdir, tag);
    f = fopen(path, "w");
    if (f) { fprintf(f, "# cost %d exact %d n %d\n", bc, ex, N); for (int i = 0; i < ms; i++) fprintf(f, "%d %d\n", ep[S[i]], el[S[i]]); fclose(f); }
    snprintf(path, sizeof path, "%s/%s_best_F.txt", outdir, tag);
    f = fopen(path, "w");
    if (f) {
        fprintf(f, "# cost %d exact %d n %d\n", bc, ex, N);
        for (int i = 0; i < nnew; i++)
            for (int k = 0; k < F.fdeg[i]; k++) if (i < F.fadj[i][k]) fprintf(f, "%d %d\n", i, F.fadj[i][k]);
        fclose(f);
    }
    return 0;
}
