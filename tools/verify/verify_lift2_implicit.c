// Independent verifier for a 2-vertex voltage lift (bi-Cayley graph) over G = Z_N x|_r Z_M, given implicitly.
// Usage: verify_lift2_implicit <D> <N> <M> <r> <k0> <k1> <c> i j ... (k0 loop voltages at base vertex 0,
//        then k1 loop voltages at base vertex 1, then c link voltages 0->1, each as a pair i j for x^i y^j)
// Vertices (v, g), index v*|G| + (i + N j). Edges: (0,g)~(0,g a^{+-1}), (1,g)~(1,g b^{+-1}), (0,g)~(1,g l).
// Left multiplication by any h in G is an automorphism transitive on each fibre, so BFS from (0,e) and (1,e)
// determine all eccentricities; both must be <= D. Also checks simplicity (distinct neighbours) at every vertex
// of a sample, r^M = 1 with r of order M, and that the two fibre eccentricities are reproduced from random sources.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef long long ll;
static ll N, M, G;
static ll *rpow;
static inline ll mulg(ll a, ll gi, ll gj) {
    ll i1 = a % N, j1 = a / N;
    ll i = (i1 + (rpow[j1] * gi) % N) % N;
    ll j = j1 + gj; if (j >= M) j -= M;
    return i + N * j;
}
static void invg(ll gi, ll gj, ll *ii, ll *jj) { *jj = (M - gj) % M; *ii = (N - (rpow[*jj] * gi) % N) % N; }

static int deg0, deg1;
static ll s0i[64], s0j[64], s1i[64], s1j[64], li[64], lj[64], linvi[64], linvj[64]; static int c;

static int neighbours(ll v, ll *out) {  // returns count
    int cnt = 0;
    if (v < G) { for (int t = 0; t < deg0 - c; t++) out[cnt++] = mulg(v, s0i[t], s0j[t]); for (int t = 0; t < c; t++) out[cnt++] = G + mulg(v, li[t], lj[t]); }
    else { ll g = v - G; for (int t = 0; t < deg1 - c; t++) out[cnt++] = G + mulg(g, s1i[t], s1j[t]); for (int t = 0; t < c; t++) out[cnt++] = mulg(g, linvi[t], linvj[t]); }
    return cnt;
}

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage\n"); return 2; }
    int D = atoi(argv[1]); N = atoll(argv[2]); M = atoll(argv[3]); ll r = atoll(argv[4]);
    int k0 = atoi(argv[5]), k1 = atoi(argv[6]); c = atoi(argv[7]);
    if (argc != 8 + 2 * (k0 + k1 + c)) { fprintf(stderr, "wrong number of voltages\n"); return 2; }
    G = N * M;
    rpow = malloc(M * sizeof(ll)); rpow[0] = 1 % N;
    for (ll j = 1; j < M; j++) { rpow[j] = (rpow[j - 1] * r) % N; if (rpow[j] == 1 % N) { printf("FAIL: r has order %lld < M\n", j); return 1; } }
    if ((rpow[M - 1] * r) % N != 1 % N) { printf("FAIL: r^M != 1\n"); return 1; }
    int p = 8;
    for (int t = 0; t < k0; t++) { ll i = atoll(argv[p++]) % N, j = atoll(argv[p++]) % M; ll ii, jj; invg(i, j, &ii, &jj); s0i[2*t] = i; s0j[2*t] = j; s0i[2*t+1] = ii; s0j[2*t+1] = jj; }
    for (int t = 0; t < k1; t++) { ll i = atoll(argv[p++]) % N, j = atoll(argv[p++]) % M; ll ii, jj; invg(i, j, &ii, &jj); s1i[2*t] = i; s1j[2*t] = j; s1i[2*t+1] = ii; s1j[2*t+1] = jj; }
    for (int t = 0; t < c; t++) { ll i = atoll(argv[p++]) % N, j = atoll(argv[p++]) % M; li[t] = i; lj[t] = j; invg(i, j, &linvi[t], &linvj[t]); }
    deg0 = 2 * k0 + c; deg1 = 2 * k1 + c;
    if (deg0 != deg1) printf("note: fibre degrees %d and %d differ\n", deg0, deg1);
    ll n = 2 * G;
    // simplicity: distinct neighbours of (0,e) and (1,e) (by transitivity this holds everywhere)
    ll nb[128];
    for (int v = 0; v < 2; v++) { int cnt = neighbours(v == 0 ? 0 : G, nb); for (int a = 0; a < cnt; a++) { if (nb[a] == (v == 0 ? 0 : G)) { printf("FAIL: loop\n"); return 1; } for (int b = a + 1; b < cnt; b++) if (nb[a] == nb[b]) { printf("FAIL: multi-edge or involutory loop voltage at fibre %d\n", v); return 1; } } }
    printf("lift over Z_%lld x|_%lld Z_%lld (|G| = %lld): n = %lld, fibre degrees %d/%d, simple\n", N, r, M, G, n, deg0, deg1);
    uint8_t *dist = malloc(n); uint32_t *queue = malloc(n * sizeof(uint32_t));
    if (!dist || !queue) { printf("alloc failed\n"); return 2; }
    int ecc_f[2] = {-1, -1}; unsigned long long st = 0x9E3779B97F4A7C15ULL;
    ll srcs[6] = {0, G, 0, 0, 0, 0};
    for (int s = 2; s < 6; s++) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; srcs[s] = (ll)(st % (unsigned long long)n); }
    for (int s = 0; s < 6; s++) {
        ll src = srcs[s]; memset(dist, 255, n);
        ll qh = 0, qt = 0; dist[src] = 0; queue[qt++] = (uint32_t)src; int ecc = 0; ll prof[64] = {0}; prof[0] = 1;
        while (qh < qt) { ll x = queue[qh++]; int dx = dist[x]; int cnt = neighbours(x, nb);
            for (int t = 0; t < cnt; t++) { ll y = nb[t]; if (dist[y] == 255) { dist[y] = dx + 1; if (dx + 1 > ecc) ecc = dx + 1; prof[dx + 1]++; queue[qt++] = (uint32_t)y; } } }
        printf("source %lld (fibre %d): reached %lld of %lld, eccentricity %d, profile:", src, src < G ? 0 : 1, qt, n, ecc);
        for (int d = 0; d <= ecc; d++) printf(" %lld", prof[d]); printf("\n");
        if (qt != n) { printf("FAIL: disconnected\n"); return 1; }
        int f = src < G ? 0 : 1;
        if (s < 2) ecc_f[f] = ecc; else if (ecc != ecc_f[f]) { printf("FAIL: eccentricity differs within fibre %d\n", f); return 1; }
    }
    int diam = ecc_f[0] > ecc_f[1] ? ecc_f[0] : ecc_f[1];
    printf("fibre eccentricities %d and %d; G acts transitively on each fibre by left multiplication, so diameter = %d: %s D = %d\n", ecc_f[0], ecc_f[1], diam, diam <= D ? "satisfies" : "EXCEEDS", D);
    return diam <= D ? 0 : 1;
}
