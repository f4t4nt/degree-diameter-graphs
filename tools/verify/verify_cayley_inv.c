// Independent verifier for a Cayley graph of G = Z_N x|_r Z_M given only by (N, M, r, S): no edge list.
// Usage: verify_cayley_implicit <D> <N> <M> <r> <nsources> i1 j1 i2 j2 ... (generators x^i y^j)
// Checks r^M = 1 (mod N) and r has order exactly M, that S u S^-1 has 2k distinct non-identity elements
// (so the Cayley graph is simple and 2k-regular), then BFS from the identity and from nsources-1 further
// vertices over the implicitly defined graph. Vertex-transitivity => the identity's eccentricity is the diameter.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef long long ll;
static ll N, M;
static ll *rpow;

static inline ll mul_idx(ll a, ll gi, ll gj) {  // a = i + N*j times generator (gi, gj)
    ll i1 = a % N, j1 = a / N;
    ll i = (i1 + (rpow[j1] * gi) % N) % N;
    ll j = j1 + gj; if (j >= M) j -= M;
    return i + N * j;
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage\n"); return 2; }
    int D = atoi(argv[1]); N = atoll(argv[2]); M = atoll(argv[3]); ll r = atoll(argv[4]); int nsrc = atoi(argv[5]);
    int k = (argc - 6) / 2; ll gi[32], gj[32];
    for (int t = 0; t < k; t++) { gi[t] = atoll(argv[6 + 2 * t]) % N; gj[t] = atoll(argv[7 + 2 * t]) % M; }
    ll n = N * M;
    rpow = malloc(M * sizeof(ll)); rpow[0] = 1 % N;
    for (ll j = 1; j < M; j++) { rpow[j] = (rpow[j - 1] * r) % N; if (rpow[j] == 1 % N) { printf("FAIL: r has order %lld < M\n", j); return 1; } }
    if ((rpow[M - 1] * r) % N != 1 % N) { printf("FAIL: r^M != 1 mod N\n"); return 1; }
    // inverses: (i,j)^-1 = (-r^{M-j} i, M-j)
    ll si[64], sj[64]; int cnt = 0;
    for (int t = 0; t < k; t++) {
        ll jj = (M - gj[t]) % M; ll ii = (N - (rpow[jj] * gi[t]) % N) % N;
        si[cnt] = gi[t]; sj[cnt] = gj[t]; cnt++; si[cnt] = ii; sj[cnt] = jj; cnt++;
    }
    // dedup S u S^-1 (involutions appear once)
    { int w = 0; for (int a = 0; a < cnt; a++) { if (si[a] == 0 && sj[a] == 0) { printf("FAIL: identity in S\n"); return 1; } int dup = 0; for (int b = 0; b < w; b++) if (si[a] == si[b] && sj[a] == sj[b]) dup = 1; if (!dup) { si[w] = si[a]; sj[w] = sj[a]; w++; } } cnt = w; }
    printf("G = Z_%lld x|_%lld Z_%lld (order %lld), %d generators, S u S^-1 has %d distinct elements: simple %d-regular Cayley graph\n", N, r, M, n, k, cnt, cnt);
    uint8_t *dist = malloc(n); uint32_t *queue = malloc(n * sizeof(uint32_t));
    if (!dist || !queue) { printf("alloc failed\n"); return 2; }
    unsigned long long st = 0x9E3779B97F4A7C15ULL; int diam0 = -1;
    for (int s = 0; s < nsrc; s++) {
        ll src = 0;
        if (s > 0) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; src = (ll)(st % (unsigned long long)n); }
        memset(dist, 255, n);
        ll qh = 0, qt = 0; dist[src] = 0; queue[qt++] = (uint32_t)src; int ecc = 0; ll prof[64] = {0}; prof[0] = 1;
        while (qh < qt) {
            ll x = queue[qh++]; int dx = dist[x];
            for (int t = 0; t < cnt; t++) {
                ll y = mul_idx(x, si[t], sj[t]);
                if (dist[y] == 255) { dist[y] = dx + 1; if (dx + 1 > ecc) ecc = dx + 1; prof[dx + 1]++; queue[qt++] = (uint32_t)y; }
            }
        }
        printf("source %lld: reached %lld of %lld, eccentricity %d, profile:", src, qt, n, ecc);
        for (int d = 0; d <= ecc; d++) printf(" %lld", prof[d]); printf("\n");
        if (qt != n) { printf("FAIL: disconnected\n"); return 1; }
        if (s == 0) diam0 = ecc; else if (ecc != diam0) { printf("FAIL: eccentricities differ\n"); return 1; }
    }
    printf("Cayley graph (vertex-transitive): diameter %d on %lld vertices, degree %d: %s D = %d\n", diam0, n, cnt, diam0 <= D ? "satisfies" : "EXCEEDS", D);
    return diam0 <= D ? 0 : 1;
}
