// Route-chart lift evaluator (Mizuno / URSA construction, generalised).
// Input text format (one file):
//   N F p m D
//   for each controller vertex c (N lines):  deg  d_1 k_1  d_2 k_2 ... (target vertex, chart id)
//   K
//   for each chart k (K blocks):  one line with F integers lin[x] (the linear part as a permutation of [0,F)),
//                                 one line: t  tau_1 ... tau_t  (translation set, t = p^r elements)
// Vertex (c, x) has neighbours (d, lin_k[x] (+) tau) for every out-edge (d, k) of c and every tau in T_k,
// where (+) is coordinate-wise addition mod p of base-p digit vectors of length m.
// Output: order, degree range, duplicate/loop count, diameter (or >D), number of ordered pairs at distance > D,
// distance distribution from vertex 0.  Options: --edges FILE writes the edge list (u < v, one per line),
// --quiet prints one summary line, --maxbad B stops after B bad pairs.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int N, F, p, m, D, K;
static int *odeg, **otgt, **ochart;
static int **lin, *tsz, **tset;
static int *add; // F x F addition table (digit-wise mod p)

static int vadd(int a, int b) {
    int r = 0, pw = 1;
    for (int i = 0; i < m; i++) { int da = a % p, db = b % p; a /= p; b /= p; r += ((da + db) % p) * pw; pw *= p; }
    return r;
}

int main(int argc, char **argv) {
    const char *efile = NULL; int quiet = 0; long maxbad = -1;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--edges")) efile = argv[++i];
        else if (!strcmp(argv[i], "--quiet")) quiet = 1;
        else if (!strcmp(argv[i], "--maxbad")) maxbad = atol(argv[++i]);
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("open"); return 2; }
    if (fscanf(f, "%d %d %d %d %d", &N, &F, &p, &m, &D) != 5) { fprintf(stderr, "bad header\n"); return 2; }
    odeg = malloc(N * sizeof(int)); otgt = malloc(N * sizeof(int *)); ochart = malloc(N * sizeof(int *));
    for (int c = 0; c < N; c++) {
        if (fscanf(f, "%d", &odeg[c]) != 1) return 2;
        otgt[c] = malloc(odeg[c] * sizeof(int)); ochart[c] = malloc(odeg[c] * sizeof(int));
        for (int j = 0; j < odeg[c]; j++) if (fscanf(f, "%d %d", &otgt[c][j], &ochart[c][j]) != 2) return 2;
    }
    if (fscanf(f, "%d", &K) != 1) return 2;
    lin = malloc(K * sizeof(int *)); tsz = malloc(K * sizeof(int)); tset = malloc(K * sizeof(int *));
    for (int k = 0; k < K; k++) {
        lin[k] = malloc(F * sizeof(int));
        for (int x = 0; x < F; x++) if (fscanf(f, "%d", &lin[k][x]) != 1) return 2;
        if (fscanf(f, "%d", &tsz[k]) != 1) return 2;
        tset[k] = malloc(tsz[k] * sizeof(int));
        for (int j = 0; j < tsz[k]; j++) if (fscanf(f, "%d", &tset[k][j]) != 1) return 2;
    }
    fclose(f);
    add = malloc((size_t)F * F * sizeof(int));
    for (int a = 0; a < F; a++) for (int b = 0; b < F; b++) add[(size_t)a * F + b] = vadd(a, b);
    long n = (long)N * F;
    // degree of each lift vertex
    int *deg = calloc(n, sizeof(int));
    for (int c = 0; c < N; c++) for (int j = 0; j < odeg[c]; j++) for (int x = 0; x < F; x++) deg[(long)c * F + x] += tsz[ochart[c][j]];
    long *off = malloc((n + 1) * sizeof(long)); off[0] = 0;
    for (long v = 0; v < n; v++) off[v + 1] = off[v] + deg[v];
    int *adj = malloc(off[n] * sizeof(int));
    long *fill = calloc(n, sizeof(long));
    for (int c = 0; c < N; c++)
        for (int j = 0; j < odeg[c]; j++) {
            int d = otgt[c][j], k = ochart[c][j];
            for (int x = 0; x < F; x++) {
                long u = (long)c * F + x;
                for (int t = 0; t < tsz[k]; t++) {
                    long w = (long)d * F + add[(size_t)lin[k][x] * F + tset[k][t]];
                    adj[off[u] + fill[u]++] = w;
                }
            }
        }
    // symmetry check: every arc must have its reverse; loops and duplicates
    long loops = 0, dups = 0, asym = 0;
    int mind = 1 << 30, maxd = 0;
    for (long v = 0; v < n; v++) {
        if (deg[v] < mind) mind = deg[v];
        if (deg[v] > maxd) maxd = deg[v];
        for (long i = off[v]; i < off[v + 1]; i++) {
            long w = adj[i];
            if (w == v) loops++;
            for (long j = i + 1; j < off[v + 1]; j++) if (adj[j] == w) dups++;
            int found = 0;
            for (long j = off[w]; j < off[w + 1]; j++) if (adj[j] == v) { found = 1; break; }
            if (!found) asym++;
        }
    }
    // BFS from every vertex
    int *dist = malloc(n * sizeof(int)); long *queue = malloc(n * sizeof(long));
    long bad = 0, unreach = 0; int diam = 0;
    long hist0[64] = {0};
    for (long s = 0; s < n; s++) {
        for (long i = 0; i < n; i++) dist[i] = -1;
        long qh = 0, qt = 0; dist[s] = 0; queue[qt++] = s;
        while (qh < qt) {
            long x = queue[qh++];
            if (dist[x] >= D) continue;  // only need distances up to D
            for (long i = off[x]; i < off[x + 1]; i++) { long y = adj[i]; if (dist[y] < 0) { dist[y] = dist[x] + 1; queue[qt++] = y; } }
        }
        long far = n - qt;  // vertices not reached within D
        bad += far;
        if (far && diam <= D) diam = D + 1;
        if (s == 0) { for (long i = 0; i < n; i++) if (dist[i] >= 0) hist0[dist[i]]++; hist0[D + 1] = far; }
        if (!far) { for (long i = 0; i < qt; i++) if (dist[queue[i]] > diam) diam = dist[queue[i]]; }
        if (maxbad >= 0 && bad > maxbad) break;
    }
    if (quiet) {
        printf("n=%ld deg=%d..%d loops=%ld dups=%ld asym=%ld diam=%s%d badpairs=%ld\n", n, mind, maxd, loops, dups, asym, diam > D ? ">" : "", diam > D ? D : diam, bad);
    } else {
        printf("order: %ld\nmin degree: %d\nmax degree: %d\nloops: %ld\nduplicate arcs: %ld\nasymmetric arcs: %ld\n", n, mind, maxd, loops, dups, asym);
        if (diam > D) printf("diameter: > %d\n", D); else printf("diameter: %d\n", diam);
        printf("ordered pairs at distance > %d: %ld\n", D, bad);
        printf("distance distribution from vertex 0:");
        for (int i = 0; i <= D + 1; i++) printf(" %ld", hist0[i]);
        printf("\n");
    }
    if (efile) {
        FILE *g = fopen(efile, "w");
        for (long v = 0; v < n; v++) for (long i = off[v]; i < off[v + 1]; i++) if (adj[i] > v) fprintf(g, "%ld %d\n", v, adj[i]);
        fclose(g);
    }
    return 0;
}
