// gtool: BFS utilities on an edge list (vertices 0..n-1, one "u v" per line).
//   gtool far  EDGES D [DELFILE]     list unordered pairs at distance > D after deleting the edges in DELFILE
//   gtool crit EDGES D               per edge e: for d = 1..D the number of unordered pairs (u,v) with
//                                    dist(u,v) = d all of whose geodesics use e ("critical pairs"); plus the
//                                    distance distribution from every vertex (summed) and the max |B_{D-1}|
//   gtool ecc  EDGES D [DELFILE]     per vertex: eccentricity and |B_{D-1}(v)| after deletions
//   gtool cover EDGES R K [CANDFILE] greedy maximum coverage of V by K balls B_R(c), centres c from CANDFILE
//                                    (one vertex per line; default all vertices): the uncovered count after every
//                                    pick (an upper bound on what K neighbours of one new vertex can reach in R+1)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, *deg, *off, *adj, *eid;   // eid[arc position] = edge id
static size_t m;
static int *eu, *ev;

static void read_graph(const char *path, const char *delpath) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(2); }
    size_t cap = 1 << 16;
    eu = malloc(cap * sizeof(int)); ev = malloc(cap * sizeof(int));
    int a, b, maxv = -1;
    m = 0;
    while (fscanf(f, "%d %d", &a, &b) == 2) {
        if (m == cap) { cap *= 2; eu = realloc(eu, cap * sizeof(int)); ev = realloc(ev, cap * sizeof(int)); }
        eu[m] = a < b ? a : b; ev[m] = a < b ? b : a; m++;
        if (a > maxv) maxv = a;
        if (b > maxv) maxv = b;
    }
    fclose(f);
    n = maxv + 1;
    if (delpath) {
        FILE *g = fopen(delpath, "r");
        if (!g) { perror(delpath); exit(2); }
        size_t ndel = 0;
        while (fscanf(g, "%d %d", &a, &b) == 2) {
            int lo = a < b ? a : b, hi = a < b ? b : a, found = 0;
            for (size_t i = 0; i < m; i++)
                if (eu[i] == lo && ev[i] == hi) { eu[i] = eu[m - 1]; ev[i] = ev[m - 1]; m--; found = 1; break; }
            if (!found) { fprintf(stderr, "deletion %d %d not an edge\n", a, b); exit(3); }
            ndel++;
        }
        fclose(g);
        fprintf(stderr, "deleted %zu edges\n", ndel);
    }
    deg = calloc(n, sizeof(int));
    for (size_t i = 0; i < m; i++) { deg[eu[i]]++; deg[ev[i]]++; }
    off = malloc((n + 1) * sizeof(int));
    off[0] = 0;
    for (int i = 0; i < n; i++) off[i + 1] = off[i] + deg[i];
    adj = malloc(off[n] * sizeof(int));
    eid = malloc(off[n] * sizeof(int));
    int *fill = calloc(n, sizeof(int));
    for (size_t i = 0; i < m; i++) {
        adj[off[eu[i]] + fill[eu[i]]] = ev[i]; eid[off[eu[i]] + fill[eu[i]]++] = (int)i;
        adj[off[ev[i]] + fill[ev[i]]] = eu[i]; eid[off[ev[i]] + fill[ev[i]]++] = (int)i;
    }
    free(fill);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: gtool far|crit|ecc EDGES D [DELFILE]\n"); return 2; }
    const char *mode = argv[1];
    int D = atoi(argv[3]);
    if (!strcmp(mode, "cover")) {
        read_graph(argv[2], NULL);
        int R = D, K = atoi(argv[4]);
        int *cand = malloc(n * sizeof(int)), nc = 0;
        if (argc > 5) {
            FILE *g = fopen(argv[5], "r");
            if (!g) { perror(argv[5]); return 2; }
            int c;
            while (fscanf(g, "%d", &c) == 1) cand[nc++] = c;
            fclose(g);
        } else for (int i = 0; i < n; i++) cand[nc++] = i;
        int W = (n + 63) / 64;
        unsigned long long *balls = calloc((size_t)nc * W, sizeof(unsigned long long));
        int *dist = malloc(n * sizeof(int)), *queue = malloc(n * sizeof(int));
        for (int ci = 0; ci < nc; ci++) {
            int s = cand[ci];
            for (int i = 0; i < n; i++) dist[i] = -1;
            int qh = 0, qt = 0;
            dist[s] = 0; queue[qt++] = s;
            unsigned long long *b = balls + (size_t)ci * W;
            while (qh < qt) {
                int x = queue[qh++];
                b[x >> 6] |= 1ULL << (x & 63);
                if (dist[x] >= R) continue;
                for (int i = off[x]; i < off[x + 1]; i++) {
                    int y = adj[i];
                    if (dist[y] < 0) { dist[y] = dist[x] + 1; queue[qt++] = y; }
                }
            }
        }
        unsigned long long *cov = calloc(W, sizeof(unsigned long long));
        char *used = calloc(nc, 1);
        long covered = 0;
        printf("# greedy cover: n %d candidates %d R %d K %d\n", n, nc, R, K);
        for (int k = 0; k < K; k++) {
            int best = -1;
            long bestgain = -1;
            for (int ci = 0; ci < nc; ci++) {
                if (used[ci]) continue;
                unsigned long long *b = balls + (size_t)ci * W;
                long gain = 0;
                for (int w = 0; w < W; w++) gain += __builtin_popcountll(b[w] & ~cov[w]);
                if (gain > bestgain) { bestgain = gain; best = ci; }
            }
            if (best < 0) break;
            used[best] = 1;
            unsigned long long *b = balls + (size_t)best * W;
            for (int w = 0; w < W; w++) cov[w] |= b[w];
            covered += bestgain;
            printf("pick %d: centre %d gain %ld covered %ld uncovered %ld\n", k + 1, cand[best], bestgain, covered, (long)n - covered);
        }
        return 0;
    }
    if (!strcmp(mode, "dist")) {
        // gtool dist EDGES R SRCFILE [DELFILE]: for every source s in SRCFILE one line "s u1 d1 u2 d2 ..." (d <= R)
        read_graph(argv[2], argc > 5 ? argv[5] : NULL);
        int R = D;
        FILE *g = fopen(argv[4], "r");
        if (!g) { perror(argv[4]); return 2; }
        int *dist = malloc(n * sizeof(int)), *queue = malloc(n * sizeof(int));
        int s;
        while (fscanf(g, "%d", &s) == 1) {
            for (int i = 0; i < n; i++) dist[i] = -1;
            int qh = 0, qt = 0;
            dist[s] = 0; queue[qt++] = s;
            while (qh < qt) {
                int x = queue[qh++];
                if (dist[x] >= R) continue;
                for (int i = off[x]; i < off[x + 1]; i++) {
                    int y = adj[i];
                    if (dist[y] < 0) { dist[y] = dist[x] + 1; queue[qt++] = y; }
                }
            }
            printf("%d", s);
            for (int i = 1; i < qt; i++) printf(" %d %d", queue[i], dist[queue[i]]);
            printf("\n");
        }
        return 0;
    }
    read_graph(argv[2], argc > 4 ? argv[4] : NULL);
    int *dist = malloc(n * sizeof(int)), *queue = malloc(n * sizeof(int));
    if (!strcmp(mode, "far") || !strcmp(mode, "ecc")) {
        long far = 0;
        int diam = 0;
        for (int s = 0; s < n; s++) {
            for (int i = 0; i < n; i++) dist[i] = -1;
            int qh = 0, qt = 0, ecc = 0;
            long ball = 0;
            dist[s] = 0; queue[qt++] = s;
            while (qh < qt) {
                int x = queue[qh++];
                if (dist[x] > ecc) ecc = dist[x];
                if (dist[x] <= D - 1) ball++;
                for (int i = off[x]; i < off[x + 1]; i++) {
                    int y = adj[i];
                    if (dist[y] < 0) { dist[y] = dist[x] + 1; queue[qt++] = y; }
                }
            }
            if (!strcmp(mode, "ecc")) printf("%d %d %ld\n", s, qt < n ? -1 : ecc, ball);
            else
                for (int t = s + 1; t < n; t++)
                    if (dist[t] < 0 || dist[t] > D) { printf("%d %d %d\n", s, t, dist[t]); far++; }
            if (ecc > diam) diam = ecc;
        }
        if (!strcmp(mode, "far")) fprintf(stderr, "n %d m %zu diameter %d far_pairs %ld\n", n, m, diam, far);
        return 0;
    }
    if (strcmp(mode, "crit")) { fprintf(stderr, "bad mode\n"); return 2; }
    // criticality: deps[v] = intersection over predecessors p of (deps[p] u {edge(p,v)}), as sorted small arrays
    long *crit = calloc((size_t)m * (D + 1), sizeof(long));
    int *deps = malloc((size_t)n * D * sizeof(int));   // up to D edge ids per vertex
    int *ndeps = malloc(n * sizeof(int));
    long *ddist = calloc(n + 1, sizeof(long));
    long maxball = 0;
    int *tmp = malloc((D + 1) * sizeof(int));
    for (int s = 0; s < n; s++) {
        for (int i = 0; i < n; i++) dist[i] = -1;
        int qh = 0, qt = 0;
        dist[s] = 0; queue[qt++] = s; ndeps[s] = 0;
        long ball = 1;
        while (qh < qt) {
            int x = queue[qh++];
            if (dist[x] >= D) continue;
            for (int i = off[x]; i < off[x + 1]; i++) {
                int y = adj[i];
                if (dist[y] < 0) {
                    dist[y] = dist[x] + 1; queue[qt++] = y;
                    if (dist[y] <= D - 1) ball++;
                    // first predecessor: deps[y] = deps[x] u {e}
                    int *dy = deps + (size_t)y * D, *dx = deps + (size_t)x * D, k = 0, e = eid[i];
                    // keep sorted
                    int placed = 0;
                    for (int j = 0; j < ndeps[x]; j++) {
                        if (!placed && e < dx[j]) { dy[k++] = e; placed = 1; }
                        dy[k++] = dx[j];
                    }
                    if (!placed) dy[k++] = e;
                    ndeps[y] = k;
                } else if (dist[y] == dist[x] + 1) {
                    // another predecessor: intersect deps[y] with deps[x] u {e}
                    int *dy = deps + (size_t)y * D, *dx = deps + (size_t)x * D, e = eid[i], k = 0;
                    for (int j = 0; j < ndeps[y]; j++) {
                        int v = dy[j], in = (v == e);
                        for (int l = 0; l < ndeps[x] && !in; l++) if (dx[l] == v) in = 1;
                        if (in) tmp[k++] = v;
                    }
                    memcpy(dy, tmp, k * sizeof(int));
                    ndeps[y] = k;
                }
            }
        }
        for (int i = 0; i < qt; i++) {
            int v = queue[i];
            ddist[dist[v]]++;
            if (v == s) continue;
            int *dv = deps + (size_t)v * D;
            for (int j = 0; j < ndeps[v]; j++) crit[(size_t)dv[j] * (D + 1) + dist[v]]++;
        }
        ddist[D + 1 > n ? n : D + 1] += 0;
        long reached = qt;
        if (reached < n) ddist[0] += 0;  // (vertices beyond D are not in queue: counted below)
        if (ball > maxball) maxball = ball;
        (void)reached;
    }
    printf("# n %d m %zu D %d\n", n, m, D);
    printf("# distance distribution (summed over sources, distances 0..%d; beyond: %ld): ", D, (long)n * n - 0);
    long tot = 0;
    for (int d = 0; d <= D; d++) { printf("%ld ", ddist[d]); tot += ddist[d]; }
    printf("| beyond_D %ld\n", (long)n * n - tot);
    printf("# max |B_%d(v)| = %ld\n", D - 1, maxball);
    printf("# edge u v crit_1 .. crit_%d (unordered pairs whose every geodesic uses the edge)\n", D);
    for (size_t e = 0; e < m; e++) {
        printf("%d %d", eu[e], ev[e]);
        for (int d = 1; d <= D; d++) printf(" %ld", crit[e * (D + 1) + d] / 2);
        printf("\n");
    }
    return 0;
}
