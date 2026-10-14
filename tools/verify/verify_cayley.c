// Independent verifier for a Cayley graph of a metacyclic group Z_N x|_r Z_M given as an edge list.
// Usage: verify_cayley <edges file> <D> <N> <M> <r> <gen1> <gen2> ...  (generator index = i + N*j for x^i y^j)
// Checks: r^M = 1 mod N; the file's edge set equals {g ~ g*s, g*s^-1} exactly (so the graph IS the
// Cayley graph, hence vertex-transitive); then BFS from the identity gives the diameter, cross-checked
// by BFS from 50 pseudo-random vertices.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef long long ll;
static ll N, M, R;
static ll *rpow;  // r^j mod N for j < M

static inline ll mul(ll a, ll b) {  // elements as index i + N*j
    ll i1 = a % N, j1 = a / N, i2 = b % N, j2 = b / N;
    ll i = (i1 + (rpow[j1] * i2) % N) % N;
    ll j = (j1 + j2) % M;
    return i + N * j;
}
static ll inv(ll a) {  // brute: find b with mul(a,b)==0 using structure: (i,j)^-1 = (-r^{-j} i, -j)
    ll i = a % N, j = a / N;
    ll jj = (M - j) % M;
    ll ii = (N - (rpow[jj] * i) % N) % N;  // r^{-j} = r^{M-j}
    return ii + N * jj;
}
static int cmpll(const void *x, const void *y) { ll a = *(const ll *)x, b = *(const ll *)y; return (a > b) - (a < b); }

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage\n"); return 2; }
    const char *path = argv[1]; int D = atoi(argv[2]); N = atoll(argv[3]); M = atoll(argv[4]); R = atoll(argv[5]);
    int k = argc - 6; ll gens[16];
    for (int t = 0; t < k; t++) gens[t] = atoll(argv[6 + t]);
    ll n = N * M;
    rpow = malloc(M * sizeof(ll)); rpow[0] = 1 % N;
    for (ll j = 1; j < M; j++) rpow[j] = (rpow[j - 1] * R) % N;
    if ((rpow[M - 1] * R) % N != 1 % N) { printf("FAIL: r^M != 1 mod N\n"); return 1; }
    // expected adjacency
    int deg = 2 * k;
    ll *exp_adj = malloc((size_t)n * deg * sizeof(ll));
    for (ll g = 0; g < n; g++) {
        for (int t = 0; t < k; t++) {
            exp_adj[g * deg + 2 * t] = mul(g, gens[t]);
            exp_adj[g * deg + 2 * t + 1] = mul(g, inv(gens[t]));
        }
        qsort(exp_adj + g * deg, deg, sizeof(ll), cmpll);
    }
    for (int t = 0; t < k; t++) {  // generators must not be involutions or identity, and distinct pairs
        if (gens[t] == 0 || inv(gens[t]) == gens[t]) { printf("FAIL: generator %lld is identity or involution\n", gens[t]); return 1; }
    }
    // read file
    FILE *f = fopen(path, "r"); if (!f) { perror("open"); return 2; }
    ll *cnt = calloc(n, sizeof(ll)); ll *adj = malloc((size_t)n * deg * sizeof(ll));
    ll a, b, m = 0;
    while (fscanf(f, "%lld %lld", &a, &b) == 2) {
        if (a < 0 || b < 0 || a >= n || b >= n || a == b) { printf("FAIL: bad edge %lld %lld\n", a, b); return 1; }
        if (cnt[a] >= deg || cnt[b] >= deg) { printf("FAIL: degree exceeds %d at edge %lld %lld\n", deg, a, b); return 1; }
        adj[a * deg + cnt[a]++] = b; adj[b * deg + cnt[b]++] = a; m++;
    }
    fclose(f);
    for (ll g = 0; g < n; g++) {
        if (cnt[g] != deg) { printf("FAIL: vertex %lld has degree %lld\n", g, cnt[g]); return 1; }
        qsort(adj + g * deg, deg, sizeof(ll), cmpll);
        if (memcmp(adj + g * deg, exp_adj + g * deg, deg * sizeof(ll)) != 0) { printf("FAIL: neighbourhood of %lld differs from the Cayley graph\n", g); return 1; }
        for (int t = 1; t < deg; t++) if (adj[g * deg + t] == adj[g * deg + t - 1]) { printf("FAIL: duplicate edge at %lld\n", g); return 1; }
    }
    printf("edge set equals Cay(Z_%lld x|_%lld Z_%lld, S) with |S| = %d: n = %lld, m = %lld\n", N, R, M, deg, n, m);
    // BFS from identity and from samples
    int *dist = malloc(n * sizeof(int)); ll *queue = malloc(n * sizeof(ll));
    ll seeds[51]; seeds[0] = 0; unsigned long long st = 88172645463325252ULL;
    for (int s = 1; s <= 50; s++) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; seeds[s] = (ll)(st % (unsigned long long)n); }
    int diam_id = -1;
    for (int s = 0; s <= 50; s++) {
        for (ll i = 0; i < n; i++) dist[i] = -1;
        ll qh = 0, qt = 0; dist[seeds[s]] = 0; queue[qt++] = seeds[s]; int ecc = 0;
        while (qh < qt) { ll x = queue[qh++]; for (int t = 0; t < deg; t++) { ll y = adj[x * deg + t]; if (dist[y] < 0) { dist[y] = dist[x] + 1; if (dist[y] > ecc) ecc = dist[y]; queue[qt++] = y; } } }
        if (qt != n) { printf("FAIL: disconnected from %lld\n", seeds[s]); return 1; }
        if (s == 0) {
            diam_id = ecc;
            ll prof[64] = {0}; for (ll i = 0; i < n; i++) prof[dist[i]]++;
            printf("eccentricity of identity: %d; profile:", ecc); for (int d = 0; d <= ecc; d++) printf(" %lld", prof[d]); printf("\n");
        } else if (ecc != diam_id) { printf("FAIL: eccentricity %d at vertex %lld differs (transitivity violated?)\n", ecc, seeds[s]); return 1; }
    }
    printf("50 sampled vertices agree. Vertex-transitive Cayley graph => diameter = %d (%s D = %d)\n", diam_id, diam_id <= D ? "satisfies" : "EXCEEDS", D);
    return diam_id <= D ? 0 : 1;
}
