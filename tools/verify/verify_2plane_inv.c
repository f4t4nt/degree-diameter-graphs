// Two-bit-plane independent verifier for implicitly defined Cayley graphs / 2-vertex lifts over
// G = Z_N x|_r Z_M, for n up to ~3e10 under 8 GB (2 bits per vertex).
// States per vertex: (P0,P1) = (0,0) unvisited, (1,0) current frontier, (0,1) next frontier, (1,1) done.
// Usage: verify_2plane_implicit cay  <D> <N> <M> <r> i1 j1 ...
//        verify_2plane_implicit lift <D> <N> <M> <r> <k0> <k1> <c> i j ...
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef unsigned long long ull;
static ull N, M, G; static ull *rpow;
static inline ull mulg(ull a, ull gi, ull gj) { ull i1 = a % N, j1 = a / N; ull i = (i1 + (rpow[j1] * gi) % N) % N; ull j = j1 + gj; if (j >= M) j -= M; return i + N * j; }
static void invg(ull gi, ull gj, ull *ii, ull *jj) { *jj = (M - gj) % M; *ii = (N - (rpow[*jj] * gi) % N) % N; }
static int lift, deg0, deg1, c;
static ull s0i[64], s0j[64], s1i[64], s1j[64], li[64], lj[64], linvi[64], linvj[64];
static int neighbours(ull v, ull *out) {
    int cnt = 0;
    if (!lift) { for (int t = 0; t < deg0; t++) out[cnt++] = mulg(v, s0i[t], s0j[t]); return cnt; }
    if (v < G) { for (int t = 0; t < deg0 - c; t++) out[cnt++] = mulg(v, s0i[t], s0j[t]); for (int t = 0; t < c; t++) out[cnt++] = G + mulg(v, li[t], lj[t]); }
    else { ull g = v - G; for (int t = 0; t < deg1 - c; t++) out[cnt++] = G + mulg(g, s1i[t], s1j[t]); for (int t = 0; t < c; t++) out[cnt++] = mulg(g, linvi[t], linvj[t]); }
    return cnt;
}
#define GET(b, i) ((b[(i) >> 6] >> ((i) & 63)) & 1ULL)
#define SET(b, i) (b[(i) >> 6] |= 1ULL << ((i) & 63))

static int bfs(ull src, ull n, ull *P0, ull *P1, ull *prof) {
    ull words = (n + 63) / 64;
    memset(P0, 0, words * 8); memset(P1, 0, words * 8);
    SET(P0, src); ull reached = 1; prof[0] = 1; int level = 0; ull nb[128];
    while (1) {
        ull added = 0;
        for (ull w = 0; w < words; w++) {
            ull bits = P0[w] & ~P1[w];  // current frontier
            while (bits) { int b = __builtin_ctzll(bits); bits &= bits - 1; ull v = w * 64 + b;
                int cnt = neighbours(v, nb);
                for (int t = 0; t < cnt; t++) { ull y = nb[t]; ull wy = y >> 6, by = 1ULL << (y & 63);
                    if (!((P0[wy] | P1[wy]) & by)) { P1[wy] |= by; added++; } } }
        }
        if (!added) break;
        level++; prof[level] = added; reached += added;
        for (ull w = 0; w < words; w++) { ull cur = P0[w] & ~P1[w], nxt = ~P0[w] & P1[w]; P0[w] |= nxt; P1[w] = (P1[w] | cur) & ~nxt; }
        if (level > 60) break;
    }
    if (reached != n) { printf("FAIL: reached %llu of %llu\n", reached, n); return -1; }
    return level;
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage\n"); return 2; }
    lift = strcmp(argv[1], "lift") == 0; int D = atoi(argv[2]); N = strtoull(argv[3], 0, 10); M = strtoull(argv[4], 0, 10); ull r = strtoull(argv[5], 0, 10);
    G = N * M; rpow = malloc(M * sizeof(ull)); rpow[0] = 1 % N;
    for (ull j = 1; j < M; j++) { rpow[j] = (rpow[j - 1] * r) % N; if (rpow[j] == 1 % N) { printf("FAIL: r has order %llu < M\n", j); return 1; } }
    if ((rpow[M - 1] * r) % N != 1 % N) { printf("FAIL: r^M != 1\n"); return 1; }
    // overflow guard for rpow[j1]*gi: needs N < 2^32 for the product to fit in 64 bits
    if (N >= (1ULL << 32)) { printf("FAIL: N too large for 64-bit products\n"); return 2; }
    int p = 6; ull n;
    if (!lift) {
        int k = (argc - 6) / 2; deg0 = 2 * k;
        for (int t = 0; t < k; t++) { ull i = strtoull(argv[p++], 0, 10) % N, j = strtoull(argv[p++], 0, 10) % M; ull ii, jj; invg(i, j, &ii, &jj); s0i[2*t] = i; s0j[2*t] = j; s0i[2*t+1] = ii; s0j[2*t+1] = jj; }
        { int w = 0; for (int a = 0; a < deg0; a++) { if (s0i[a] == 0 && s0j[a] == 0) { printf("FAIL: identity in S\n"); return 1; } int dup = 0; for (int b = 0; b < w; b++) if (s0i[a] == s0i[b] && s0j[a] == s0j[b]) dup = 1; if (!dup) { s0i[w] = s0i[a]; s0j[w] = s0j[a]; w++; } } deg0 = w; }  // involutions appear once
        n = G; printf("Cay(Z_%llu x|_%llu Z_%llu, S), |S u S^-1| = %d, n = %llu\n", N, r, M, deg0, n);
    } else {
        int k0 = atoi(argv[6]), k1 = atoi(argv[7]); c = atoi(argv[8]); p = 9;
        for (int t = 0; t < k0; t++) { ull i = strtoull(argv[p++], 0, 10) % N, j = strtoull(argv[p++], 0, 10) % M; ull ii, jj; invg(i, j, &ii, &jj); s0i[2*t] = i; s0j[2*t] = j; s0i[2*t+1] = ii; s0j[2*t+1] = jj; }
        for (int t = 0; t < k1; t++) { ull i = strtoull(argv[p++], 0, 10) % N, j = strtoull(argv[p++], 0, 10) % M; ull ii, jj; invg(i, j, &ii, &jj); s1i[2*t] = i; s1j[2*t] = j; s1i[2*t+1] = ii; s1j[2*t+1] = jj; }
        for (int t = 0; t < c; t++) { ull i = strtoull(argv[p++], 0, 10) % N, j = strtoull(argv[p++], 0, 10) % M; li[t] = i; lj[t] = j; invg(i, j, &linvi[t], &linvj[t]); }
        deg0 = 2 * k0 + c; deg1 = 2 * k1 + c; n = 2 * G;
        ull nb[128]; for (int v = 0; v < 2; v++) { int cnt = neighbours(v ? G : 0, nb); for (int a = 0; a < cnt; a++) { if (nb[a] == (v ? G : 0)) { printf("FAIL: loop\n"); return 1; } for (int b = a + 1; b < cnt; b++) if (nb[a] == nb[b]) { printf("FAIL: multi-edge\n"); return 1; } } }
        printf("lift over Z_%llu x|_%llu Z_%llu, fibre degrees %d/%d, n = %llu\n", N, r, M, deg0, deg1, n);
    }
    ull words = (n + 63) / 64; ull *P0 = malloc(words * 8), *P1 = malloc(words * 8);
    if (!P0 || !P1) { printf("alloc failed\n"); return 2; }
    ull prof[64]; int diam = 0;
    for (int s = 0; s < (lift ? 2 : 1); s++) {
        memset(prof, 0, sizeof prof); int e = bfs(s ? G : 0, n, P0, P1, prof); if (e < 0) return 1;
        printf("source %s: eccentricity %d, profile:", s ? "(1,e)" : "identity", e); for (int d = 0; d <= e; d++) printf(" %llu", prof[d]); printf("\n"); fflush(stdout);
        if (e > diam) diam = e;
    }
    printf("%s: diameter %d on %llu vertices: %s D = %d\n", lift ? "lift (two fibre BFS)" : "Cayley graph (vertex-transitive)", diam, n, diam <= D ? "satisfies" : "EXCEEDS", D);
    return diam <= D ? 0 : 1;
}
