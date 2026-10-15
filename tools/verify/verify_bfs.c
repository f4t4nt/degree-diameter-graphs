// Independent plain-BFS verifier for an edge list: order, degrees, simplicity, diameter,
// number of unordered pairs at distance > D.  Lines starting with '#' are ignored.
// Usage: verify_bfs <edges file> <D>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s edges D\n", argv[0]); return 2; }
    int D = atoi(argv[2]);
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("open"); return 2; }
    size_t cap = 1 << 20, m = 0;
    int *eu = malloc(cap * sizeof(int)), *ev = malloc(cap * sizeof(int));
    int a, b, maxv = -1;
    char line[4096];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "%d %d", &a, &b) != 2) continue;
        if (m == cap) { cap *= 2; eu = realloc(eu, cap * sizeof(int)); ev = realloc(ev, cap * sizeof(int)); }
        if (a == b) { printf("LOOP at %d\n", a); return 1; }
        eu[m] = a; ev[m] = b; m++;
        if (a > maxv) maxv = a;
        if (b > maxv) maxv = b;
    }
    fclose(f);
    int n = maxv + 1;
    int *deg = calloc(n, sizeof(int));
    for (size_t i = 0; i < m; i++) { deg[eu[i]]++; deg[ev[i]]++; }
    int *off = malloc((n + 1) * sizeof(int));
    off[0] = 0;
    for (int i = 0; i < n; i++) off[i + 1] = off[i] + deg[i];
    int *adj = malloc(off[n] * sizeof(int));
    int *fill = calloc(n, sizeof(int));
    for (size_t i = 0; i < m; i++) {
        adj[off[eu[i]] + fill[eu[i]]++] = ev[i];
        adj[off[ev[i]] + fill[ev[i]]++] = eu[i];
    }
    // simplicity: duplicate edges
    long dup = 0;
    for (int v = 0; v < n; v++)
        for (int i = off[v]; i < off[v + 1]; i++)
            for (int j = i + 1; j < off[v + 1]; j++)
                if (adj[i] == adj[j]) dup++;
    int mind = 1 << 30, maxd = 0, isolated = 0;
    for (int v = 0; v < n; v++) { if (deg[v] < mind) mind = deg[v]; if (deg[v] > maxd) maxd = deg[v]; if (deg[v] == 0) isolated++; }
    int *dist = malloc(n * sizeof(int)), *queue = malloc(n * sizeof(int));
    long far = 0, unreachable = 0;
    int diam = 0;
    for (int s = 0; s < n; s++) {
        for (int i = 0; i < n; i++) dist[i] = -1;
        int qh = 0, qt = 0;
        dist[s] = 0; queue[qt++] = s;
        while (qh < qt) {
            int x = queue[qh++];
            for (int i = off[x]; i < off[x + 1]; i++) {
                int y = adj[i];
                if (dist[y] < 0) { dist[y] = dist[x] + 1; queue[qt++] = y; }
            }
        }
        for (int t = s + 1; t < n; t++) {
            if (dist[t] < 0) unreachable++;
            else { if (dist[t] > diam) diam = dist[t]; if (dist[t] > D) far++; }
        }
    }
    printf("file: %s\norder: %d\nedges: %zu\nmin degree: %d\nmax degree: %d\nisolated: %d\nduplicate edges: %ld\n"
           "diameter: %d\nunreachable pairs: %ld\nunordered pairs at distance > %d: %ld\n",
           argv[1], n, m, mind, maxd, isolated, dup, diam, unreachable, D, far);
    return 0;
}
