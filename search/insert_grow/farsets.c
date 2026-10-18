// farsets: for every edge e of a pool, the unordered pairs (u, v) with d_G(u, v) <= D and d_{G - e}(u, v) > D
// ("the far pairs created by deleting e alone").  Candidates are the pairs all of whose geodesics use e (dependency
// sets along the BFS DAG as in gtool crit); each candidate is re-checked by a BFS in G - e from its source.
// usage: farsets EDGES D POOLFILE        POOLFILE: one "u v" per line (must be edges)
// output: one line per pool edge "u v : u1 v1 u2 v2 ..." (far pairs, u_i < v_i)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, *deg, *off, *adj, *eid;
static size_t m;
static int *eu, *ev;

typedef struct { int *a; int len, cap; } Vec;
static void push(Vec *v, int x) {
    if (v->len == v->cap) { v->cap = v->cap ? 2 * v->cap : 16; v->a = realloc(v->a, v->cap * sizeof(int)); }
    v->a[v->len++] = x;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: farsets EDGES D POOLFILE\n"); return 2; }
    int D = atoi(argv[2]);
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 2; }
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
    // pool
    int *pidx = malloc(m * sizeof(int));
    for (size_t i = 0; i < m; i++) pidx[i] = -1;
    FILE *g = fopen(argv[3], "r");
    if (!g) { perror(argv[3]); return 2; }
    int np = 0, pcap = 1024;
    int *pe = malloc(pcap * sizeof(int));
    while (fscanf(g, "%d %d", &a, &b) == 2) {
        int lo = a < b ? a : b, hi = a < b ? b : a, found = -1;
        for (int i = off[lo]; i < off[lo + 1]; i++) if (adj[i] == hi) { found = eid[i]; break; }
        if (found < 0) { fprintf(stderr, "pool edge %d %d not an edge\n", a, b); return 3; }
        if (pidx[found] >= 0) continue;
        if (np == pcap) { pcap *= 2; pe = realloc(pe, pcap * sizeof(int)); }
        pidx[found] = np; pe[np++] = found;
    }
    fclose(g);
    Vec *cands = calloc(np, sizeof(Vec));   // pairs (s, v) with s < v, flattened
    int *dist = malloc(n * sizeof(int)), *queue = malloc(n * sizeof(int));
    int *deps = malloc((size_t)n * D * sizeof(int)), *ndeps = malloc(n * sizeof(int)), *tmp = malloc((D + 1) * sizeof(int));
    for (int s = 0; s < n; s++) {
        for (int i = 0; i < n; i++) dist[i] = -1;
        int qh = 0, qt = 0;
        dist[s] = 0; queue[qt++] = s; ndeps[s] = 0;
        while (qh < qt) {
            int x = queue[qh++];
            if (dist[x] >= D) continue;
            for (int i = off[x]; i < off[x + 1]; i++) {
                int y = adj[i];
                if (dist[y] < 0) {
                    dist[y] = dist[x] + 1; queue[qt++] = y;
                    int *dy = deps + (size_t)y * D, *dx = deps + (size_t)x * D, k = 0, e = eid[i], placed = 0;
                    for (int j = 0; j < ndeps[x]; j++) { if (!placed && e < dx[j]) { dy[k++] = e; placed = 1; } dy[k++] = dx[j]; }
                    if (!placed) dy[k++] = e;
                    ndeps[y] = k;
                } else if (dist[y] == dist[x] + 1) {
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
        for (int i = 1; i < qt; i++) {
            int v = queue[i];
            if (v < s) continue;
            int *dv = deps + (size_t)v * D;
            for (int j = 0; j < ndeps[v]; j++) if (pidx[dv[j]] >= 0) { push(&cands[pidx[dv[j]]], s); push(&cands[pidx[dv[j]]], v); }
        }
    }
    // re-check each candidate in G - e
    for (int p = 0; p < np; p++) {
        int e = pe[p];
        printf("%d %d :", eu[e], ev[e]);
        Vec *c = &cands[p];
        // group by source: candidates were appended in increasing s already (outer loop over s)
        int i = 0;
        while (i < c->len) {
            int s = c->a[i];
            for (int k = 0; k < n; k++) dist[k] = -1;
            int qh = 0, qt = 0;
            dist[s] = 0; queue[qt++] = s;
            while (qh < qt) {
                int x = queue[qh++];
                if (dist[x] >= D) continue;
                for (int t = off[x]; t < off[x + 1]; t++) {
                    if (eid[t] == e) continue;
                    int y = adj[t];
                    if (dist[y] < 0) { dist[y] = dist[x] + 1; queue[qt++] = y; }
                }
            }
            while (i < c->len && c->a[i] == s) {
                int v = c->a[i + 1];
                if (dist[v] < 0 || dist[v] > D) printf(" %d %d", s, v);
                i += 2;
            }
        }
        printf("\n");
    }
    return 0;
}
