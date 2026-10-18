// grow: greedy iterated vertex insertion with "absorbed" edge deletions, diameter D kept exactly.
//
// One step on the current graph G (max degree Delta, diameter D):
//   1. crit_D(e) for every edge (number of pairs at distance exactly D all of whose geodesics use e); an edge with
//      crit_D = 0 is "absorbable": if it is deleted and the new vertex x is joined to both endpoints, every pair that
//      depended on it (distance <= D-1) gets the path through x of length +1 <= D, and no pair at distance D did.
//   2. Greedy maximum coverage of V by radius-(D-1) balls, items = absorbable edges (both endpoints, 2 slots) and
//      vertices with a free slot (1 slot), at most Delta slots: x must reach every vertex within D.
//   3. Apply, all-pairs BFS; far pairs between x and earlier new vertices are repaired by new-new edges when the
//      degrees allow; accept at 0 far pairs, else try the next first item (up to TRIALS first items).
// usage: grow EDGES Delta D maxsteps OUTPREFIX [seed] [CMAX] [TRIALS]
//   writes OUTPREFIX_<n>.edges after every accepted step and OUTPREFIX.log (deleted edges, new neighbourhoods).
//   CMAX (default 0): edges with crit_D <= CMAX are items; with CMAX > 0 the deletions may create far pairs, which the
//   exact all-pairs check after the insertion still gates (a step is accepted only at 0 far pairs).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned long long u64;
static int n, D, Delta, n0;
static int **nb, *deg, *cap;      // adjacency as dynamic arrays (cap = capacity)

static void add_edge(int a, int b) {
    if (deg[a] == cap[a]) { cap[a] = cap[a] ? 2 * cap[a] : 8; nb[a] = realloc(nb[a], cap[a] * sizeof(int)); }
    if (deg[b] == cap[b]) { cap[b] = cap[b] ? 2 * cap[b] : 8; nb[b] = realloc(nb[b], cap[b] * sizeof(int)); }
    nb[a][deg[a]++] = b; nb[b][deg[b]++] = a;
}
static int has_edge(int a, int b) { for (int i = 0; i < deg[a]; i++) if (nb[a][i] == b) return 1; return 0; }
static void del_edge(int a, int b) {
    for (int i = 0; i < deg[a]; i++) if (nb[a][i] == b) { nb[a][i] = nb[a][--deg[a]]; break; }
    for (int i = 0; i < deg[b]; i++) if (nb[b][i] == a) { nb[b][i] = nb[b][--deg[b]]; break; }
}

static int *dist, *queue;
// BFS from s up to radius R (R < 0: unbounded); returns number reached
static int bfs(int s, int R) {
    for (int i = 0; i < n; i++) dist[i] = -1;
    int qh = 0, qt = 0;
    dist[s] = 0; queue[qt++] = s;
    while (qh < qt) {
        int x = queue[qh++];
        if (R >= 0 && dist[x] >= R) continue;
        for (int i = 0; i < deg[x]; i++) { int y = nb[x][i]; if (dist[y] < 0) { dist[y] = dist[x] + 1; queue[qt++] = y; } }
    }
    return qt;
}

// all-pairs far check; returns count and fills far pairs (up to maxfar)
static long far_pairs(int *fu, int *fv, long maxfar) {
    long far = 0;
    for (int s = 0; s < n; s++) {
        bfs(s, -1);
        for (int t = s + 1; t < n; t++)
            if (dist[t] < 0 || dist[t] > D) { if (far < maxfar) { fu[far] = s; fv[far] = t; } far++; }
    }
    return far;
}

// crit_D per edge: edges identified by (a,b) a<b through a hash of position; we store crit in a map keyed by a*n+b
// (n <= 2e4 so long keys are fine).  Uses the dependency-intersection method of gtool.c.
typedef struct { long key; long c; } CritEnt;
static CritEnt *ctab; static long csize;
static long cfind(long key) {
    long h = (key * 11400714819323198485ULL) >> 20; h &= csize - 1;
    while (ctab[h].key != -1 && ctab[h].key != key) h = (h + 1) & (csize - 1);
    return h;
}
static void compute_crit(void) {
    long m = 0;
    for (int v = 0; v < n; v++) m += deg[v];
    m /= 2;
    csize = 1; while (csize < 4 * m) csize <<= 1;
    free(ctab);
    ctab = malloc(csize * sizeof(CritEnt));
    for (long i = 0; i < csize; i++) { ctab[i].key = -1; ctab[i].c = 0; }
    for (int v = 0; v < n; v++) for (int i = 0; i < deg[v]; i++) if (v < nb[v][i]) { long h = cfind((long)v * n + nb[v][i]); ctab[h].key = (long)v * n + nb[v][i]; }
    long *deps = malloc((size_t)n * D * sizeof(long));
    int *ndeps = malloc(n * sizeof(int));
    long *tmp = malloc((D + 1) * sizeof(long));
    for (int s = 0; s < n; s++) {
        for (int i = 0; i < n; i++) dist[i] = -1;
        int qh = 0, qt = 0;
        dist[s] = 0; queue[qt++] = s; ndeps[s] = 0;
        while (qh < qt) {
            int x = queue[qh++];
            if (dist[x] >= D) continue;
            for (int i = 0; i < deg[x]; i++) {
                int y = nb[x][i];
                long e = x < y ? (long)x * n + y : (long)y * n + x;
                if (dist[y] < 0) {
                    dist[y] = dist[x] + 1; queue[qt++] = y;
                    long *dy = deps + (size_t)y * D, *dx = deps + (size_t)x * D; int k = 0, placed = 0;
                    for (int j = 0; j < ndeps[x]; j++) { if (!placed && e < dx[j]) { dy[k++] = e; placed = 1; } dy[k++] = dx[j]; }
                    if (!placed) dy[k++] = e;
                    ndeps[y] = k;
                } else if (dist[y] == dist[x] + 1) {
                    long *dy = deps + (size_t)y * D, *dx = deps + (size_t)x * D; int k = 0;
                    for (int j = 0; j < ndeps[y]; j++) {
                        long v = dy[j]; int in = (v == e);
                        for (int l = 0; l < ndeps[x] && !in; l++) if (dx[l] == v) in = 1;
                        if (in) tmp[k++] = v;
                    }
                    memcpy(dy, tmp, k * sizeof(long));
                    ndeps[y] = k;
                }
            }
        }
        for (int i = 0; i < qt; i++) {
            int v = queue[i];
            if (dist[v] != D) continue;
            long *dv = deps + (size_t)v * D;
            for (int j = 0; j < ndeps[v]; j++) ctab[cfind(dv[j])].c++;
        }
    }
    free(deps); free(ndeps); free(tmp);
}
static long crit_of(int a, int b) { long key = a < b ? (long)a * n + b : (long)b * n + a; return ctab[cfind(key)].c / 2; }

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: grow EDGES Delta D maxsteps OUTPREFIX [seed]\n"); return 2; }
    Delta = atoi(argv[2]); D = atoi(argv[3]);
    int maxsteps = atoi(argv[4]);
    const char *prefix = argv[5];
    unsigned seed = argc > 6 ? atoi(argv[6]) : 1;
    long CMAX = argc > 7 ? atol(argv[7]) : 0;
    int TRIALS = argc > 8 ? atoi(argv[8]) : 60;
    srand(seed);
    // read
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("open"); return 2; }
    int a, b, maxv = -1;
    size_t ecap = 1 << 16, m = 0;
    int *eu = malloc(ecap * sizeof(int)), *ev = malloc(ecap * sizeof(int));
    while (fscanf(f, "%d %d", &a, &b) == 2) {
        if (m == ecap) { ecap *= 2; eu = realloc(eu, ecap * sizeof(int)); ev = realloc(ev, ecap * sizeof(int)); }
        eu[m] = a; ev[m] = b; m++;
        if (a > maxv) maxv = a;
        if (b > maxv) maxv = b;
    }
    fclose(f);
    n0 = n = maxv + 1;
    int NMAX = n + maxsteps + 1;
    nb = calloc(NMAX, sizeof(int *)); deg = calloc(NMAX, sizeof(int)); cap = calloc(NMAX, sizeof(int));
    for (size_t i = 0; i < m; i++) add_edge(eu[i], ev[i]);
    dist = malloc(NMAX * sizeof(int)); queue = malloc(NMAX * sizeof(int));
    char path[1024];
    snprintf(path, sizeof path, "%s.log", prefix);
    FILE *log = fopen(path, "a");
    int W = (NMAX + 63) / 64;
    long maxfar = 1 << 16;
    int *fu = malloc(maxfar * sizeof(int)), *fv = malloc(maxfar * sizeof(int));
    for (int step = 0; step < maxsteps; step++) {
        clock_t t0 = clock();
        compute_crit();
        // items: absorbable edges (crit_D == 0) and free-slot vertices
        int nitems = 0, icap = 1 << 12;
        int *ia = malloc(icap * sizeof(int)), *ib = malloc(icap * sizeof(int));   // ib = -1 for a vertex item
        for (int v = 0; v < n; v++) {
            if (deg[v] < Delta) { if (nitems == icap) { icap *= 2; ia = realloc(ia, icap * sizeof(int)); ib = realloc(ib, icap * sizeof(int)); } ia[nitems] = v; ib[nitems++] = -1; }
            for (int i = 0; i < deg[v]; i++) {
                int w = nb[v][i];
                if (v < w && crit_of(v, w) <= CMAX) {
                    if (nitems == icap) { icap *= 2; ia = realloc(ia, icap * sizeof(int)); ib = realloc(ib, icap * sizeof(int)); }
                    ia[nitems] = v; ib[nitems++] = w;
                }
            }
        }
        int nabs = 0; for (int i = 0; i < nitems; i++) if (ib[i] >= 0) nabs++;
        fprintf(stderr, "step %d: n=%d absorbable edges %d, free-slot vertices %d (crit %.1f s)\n", step + 1, n, nabs, nitems - nabs, (double)(clock() - t0) / CLOCKS_PER_SEC);
        if (nitems == 0) break;
        // balls of the item vertices
        u64 *balls = calloc((size_t)nitems * W, sizeof(u64));
        for (int i = 0; i < nitems; i++) {
            u64 *bl = balls + (size_t)i * W;
            int reached = bfs(ia[i], D - 1);
            for (int j = 0; j < reached; j++) bl[queue[j] >> 6] |= 1ULL << (queue[j] & 63);
            if (ib[i] >= 0) { reached = bfs(ib[i], D - 1); for (int j = 0; j < reached; j++) bl[queue[j] >> 6] |= 1ULL << (queue[j] & 63); }
        }
        // greedy cover trials with different first items
        int accepted = 0;
        char *used = calloc(nitems, 1);
        for (int trial = 0; trial < TRIALS && !accepted; trial++) {
            u64 *cov = calloc(W, sizeof(u64));
            memset(used, 0, nitems);
            int slots = 0, chosen[64], nchosen = 0;
            long covered = 0;
            while (slots < Delta) {
                int best = -1; double bestscore = -1;
                for (int i = 0; i < nitems; i++) {
                    if (used[i]) continue;
                    int cost = ib[i] >= 0 ? 2 : 1;
                    if (slots + cost > Delta) continue;
                    // endpoint already chosen as a vertex item or in another chosen edge: skip (slots must be distinct)
                    int clash = 0;
                    for (int c = 0; c < nchosen && !clash; c++) {
                        int ca = ia[chosen[c]], cb = ib[chosen[c]];
                        if (ia[i] == ca || ia[i] == cb || (ib[i] >= 0 && (ib[i] == ca || ib[i] == cb))) clash = 1;
                    }
                    if (clash) continue;
                    u64 *bl = balls + (size_t)i * W;
                    long gain = 0;
                    for (int w = 0; w < W; w++) gain += __builtin_popcountll(bl[w] & ~cov[w]);
                    double score = (double)gain / cost + (trial > 0 && nchosen == 0 ? (rand() % 1000) * 0.001 * gain : 0);
                    if (nchosen == 0 && trial > 0) score = rand() % 100000;   // randomise the first item after trial 0
                    if (score > bestscore) { bestscore = score; best = i; }
                }
                if (best < 0) break;
                used[best] = 1; chosen[nchosen++] = best;
                slots += ib[best] >= 0 ? 2 : 1;
                u64 *bl = balls + (size_t)best * W;
                long gain = 0;
                for (int w = 0; w < W; w++) { gain += __builtin_popcountll(bl[w] & ~cov[w]); cov[w] |= bl[w]; }
                covered += gain;
                if (covered == n) break;
            }
            free(cov);
            if (covered < n) { if (trial == 0) fprintf(stderr, "  trial %d: greedy cover leaves %ld uncovered with %d slots\n", trial, n - covered, slots); continue; }
            // apply tentatively
            int x = n;
            n++;
            deg[x] = 0;
            int dels[64][2], ndel = 0;
            for (int c = 0; c < nchosen; c++) {
                int ca = ia[chosen[c]], cb = ib[chosen[c]];
                if (cb >= 0) { del_edge(ca, cb); dels[ndel][0] = ca; dels[ndel][1] = cb; ndel++; add_edge(x, ca); add_edge(x, cb); }
                else add_edge(x, ca);
            }
            long far = far_pairs(fu, fv, maxfar);
            int extra[64], nextra = 0;
            if (far > 0 && far <= 64) {
                // repair new-new far pairs by edges when degrees allow
                int ok = 1;
                for (long i = 0; i < far && ok; i++) {
                    int u = fu[i], v = fv[i];
                    if (!((u >= n0 && v == x) || (v >= n0 && u == x))) ok = 0;
                    else { int o = (u == x) ? v : u; if (deg[o] >= Delta || deg[x] + nextra >= Delta) ok = 0; else extra[nextra++] = o; }
                }
                if (ok) { for (int i = 0; i < nextra; i++) add_edge(x, extra[i]); far = far_pairs(fu, fv, maxfar); }
                else nextra = 0;
            }
            if (far == 0) {
                accepted = 1;
                fprintf(log, "step %d: new vertex %d (degree %d) deleted", step + 1, x, deg[x]);
                for (int i = 0; i < ndel; i++) fprintf(log, " %d-%d", dels[i][0], dels[i][1]);
                fprintf(log, " ; N(x) =");
                for (int i = 0; i < deg[x]; i++) fprintf(log, " %d", nb[x][i]);
                fprintf(log, " ; trial %d\n", trial);
                fflush(log);
                fprintf(stderr, "  accepted trial %d: n=%d, x degree %d, %d deletions, %d new-new edges (%.1f s)\n", trial, n, deg[x], ndel, nextra, (double)(clock() - t0) / CLOCKS_PER_SEC);
                snprintf(path, sizeof path, "%s_%d.edges", prefix, n);
                FILE *o = fopen(path, "w");
                for (int v = 0; v < n; v++) for (int i = 0; i < deg[v]; i++) if (v < nb[v][i]) fprintf(o, "%d %d\n", v, nb[v][i]);
                fclose(o);
                if (step > 0) { snprintf(path, sizeof path, "%s_%d.edges", prefix, n - 1); remove(path); }
            } else {
                if (trial < 3) fprintf(stderr, "  trial %d: %ld far pairs (first %d-%d)\n", trial, far, fu[0], fv[0]);
                // undo
                for (int i = 0; i < nextra; i++) del_edge(x, extra[i]);
                for (int c = nchosen - 1; c >= 0; c--) {
                    int ca = ia[chosen[c]], cb = ib[chosen[c]];
                    if (cb >= 0) { del_edge(x, ca); del_edge(x, cb); add_edge(ca, cb); } else del_edge(x, ca);
                }
                n--;
            }
        }
        free(balls); free(ia); free(ib); free(used);
        if (!accepted) { fprintf(stderr, "step %d: no accepted trial\n", step + 1); break; }
    }
    fprintf(stderr, "final n=%d\n", n);
    fclose(log);
    return 0;
}
