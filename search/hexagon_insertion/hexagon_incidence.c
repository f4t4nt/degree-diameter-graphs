/* Incidence graph of the split Cayley hexagon H(q), q in {3,4,5,7,8,9}.
 *
 * Points: points of the parabolic quadric Q(6,q): x0 x4 + x1 x5 + x2 x6 = x3^2 in PG(6,q).
 * Lines: lines of Q(6,q) whose Grassmann coordinates p_ij = x_i y_j - x_j y_i satisfy Tits' six
 * conditions (Van Maldeghem, Generalized Polygons 2.4.13):
 *     p12 = p34, p54 = p32, p20 = p35, p65 = p30, p01 = p36, p46 = p31.
 * Explicit field tables so that q = 7, 8, 9 take seconds.
 * GF(8) = GF(2)[t]/(t^3+t+1), GF(9) = GF(3)[i]/(i^2+1), GF(4) = GF(2)[w]/(w^2+w+1).
 *
 * The output edge file lists "p L" with points 0..P-1 and lines P..2P-1.  The program verifies:
 * P = (q^6-1)/(q-1) points and lines, q+1 points per line and lines per point, bipartite, girth 12
 * (BFS from a point and a line), and the distance distribution 1, q+1, (q+1)q, ..., q^5 from a
 * sample of vertices (every vertex with --all).
 *
 * build: gcc -O2 -o hexagon_incidence hexagon_incidence.c
 * usage: ./hexagon_incidence q out.edges [--all] [--points file]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int q, chr;               /* field size, characteristic */
static int add_t[16][16], mul_t[16][16], neg_t[16], inv_t[16];

static void build_field(void) {
    if (q == 3 || q == 5 || q == 7) {
        chr = q;
        for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) { add_t[a][b] = (a + b) % q; mul_t[a][b] = (a * b) % q; }
    } else if (q == 4 || q == 8) {
        chr = 2;
        int poly = (q == 4) ? 0x7 : 0xB;   /* w^2+w+1, t^3+t+1 */
        int deg = (q == 4) ? 2 : 3;
        for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
            add_t[a][b] = a ^ b;
            int r = 0, x = a;
            for (int i = 0; i < deg; i++) { if (b >> i & 1) r ^= x; x <<= 1; if (x >> deg & 1) x ^= poly; }
            mul_t[a][b] = r;
        }
    } else if (q == 9) {
        chr = 3;
        /* element a + 3b <-> a + b i, i^2 = -1 */
        for (int x = 0; x < 9; x++) for (int y = 0; y < 9; y++) {
            int a = x % 3, b = x / 3, c = y % 3, d = y / 3;
            add_t[x][y] = (a + c) % 3 + 3 * ((b + d) % 3);
            int re = ((a * c - b * d) % 3 + 3) % 3, im = (a * d + b * c) % 3;
            mul_t[x][y] = re + 3 * im;
        }
    } else { fprintf(stderr, "unsupported q\n"); exit(1); }
    for (int a = 0; a < q; a++) {
        for (int b = 0; b < q; b++) { if (add_t[a][b] == 0) neg_t[a] = b; if (mul_t[a][b] == 1) inv_t[a] = b; }
    }
}
#define ADD(a, b) add_t[a][b]
#define MUL(a, b) mul_t[a][b]
#define SUB(a, b) add_t[a][neg_t[b]]

static int quad(const int *x) {
    int s = ADD(MUL(x[0], x[4]), MUL(x[1], x[5]));
    s = ADD(s, MUL(x[2], x[6]));
    return SUB(s, MUL(x[3], x[3]));
}
static int bil(const int *x, const int *y) {
    /* polar form Q(x+y) - Q(x) - Q(y) */
    int s = 0;
    s = ADD(s, MUL(x[0], y[4])); s = ADD(s, MUL(x[4], y[0]));
    s = ADD(s, MUL(x[1], y[5])); s = ADD(s, MUL(x[5], y[1]));
    s = ADD(s, MUL(x[2], y[6])); s = ADD(s, MUL(x[6], y[2]));
    int t = MUL(x[3], y[3]); t = ADD(t, t);
    return SUB(s, t);
}
static int gr(const int *x, const int *y, int i, int j) { return SUB(MUL(x[i], y[j]), MUL(x[j], y[i])); }
static int tits(const int *x, const int *y) {
    return gr(x, y, 1, 2) == gr(x, y, 3, 4) && gr(x, y, 5, 4) == gr(x, y, 3, 2) && gr(x, y, 2, 0) == gr(x, y, 3, 5)
        && gr(x, y, 6, 5) == gr(x, y, 3, 0) && gr(x, y, 0, 1) == gr(x, y, 3, 6) && gr(x, y, 4, 6) == gr(x, y, 3, 1);
}

static int P;
static int (*pts)[7];
static int *index_of;   /* base-q code -> point index */
static int code_of(const int *x) { int c = 0; for (int i = 0; i < 7; i++) c = c * q + x[i]; return c; }
static void normalise(int *x) {
    int f = -1;
    for (int i = 0; i < 7; i++) if (x[i]) { f = x[i]; break; }
    int inv = inv_t[f];
    for (int i = 0; i < 7; i++) x[i] = MUL(inv, x[i]);
}

/* line hash set: canonical sorted (q+1)-tuples */
static int *lines; static int nlines, lcap;
static int *htab; static long hsize;
static uint64_t hash_line(const int *l) { uint64_t h = 1469598103934665603ULL; for (int i = 0; i <= q; i++) { h ^= (uint64_t)l[i] + 0x9E37; h *= 1099511628211ULL; } return h; }
static int cmpint(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }
static int cmp_line_idx(const void *A, const void *B) {
    const int *la = lines + (long)*(const int *)A * (q + 1), *lb = lines + (long)*(const int *)B * (q + 1);
    for (int i = 0; i <= q; i++) if (la[i] != lb[i]) return la[i] - lb[i];
    return 0;
}
static int add_line(int *l) {
    qsort(l, q + 1, sizeof(int), cmpint);
    uint64_t h = hash_line(l) % (uint64_t)hsize;
    for (;;) {
        int id = htab[h];
        if (id < 0) break;
        if (!memcmp(lines + (long)id * (q + 1), l, sizeof(int) * (q + 1))) return 0;
        h = (h + 1) % (uint64_t)hsize;
    }
    if (nlines >= lcap) { lcap = lcap ? 2 * lcap : 1024; lines = realloc(lines, sizeof(int) * (long)lcap * (q + 1)); }
    memcpy(lines + (long)nlines * (q + 1), l, sizeof(int) * (q + 1));
    htab[h] = nlines++;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s q out.edges [--all] [--points file]\n", argv[0]); return 1; }
    q = atoi(argv[1]);
    const char *out = argv[2], *ptsfile = 0;
    int all = 0;
    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "--all")) all = 1;
        else if (!strcmp(argv[i], "--points")) ptsfile = argv[++i];
    }
    build_field();
    double t0 = (double)clock() / CLOCKS_PER_SEC;
    long Q7 = 1; for (int i = 0; i < 7; i++) Q7 *= q;
    index_of = malloc(sizeof(int) * Q7);
    for (long c = 0; c < Q7; c++) index_of[c] = -1;
    long expectP = (Q7 / q - 1) / (q - 1);   /* (q^6-1)/(q-1) */
    pts = malloc(sizeof(int) * 7 * (expectP + 16));
    P = 0;
    int x[7];
    for (long c = 0; c < Q7; c++) {
        long t = c;
        for (int i = 6; i >= 0; i--) { x[i] = (int)(t % q); t /= q; }
        int f = -1;
        for (int i = 0; i < 7; i++) if (x[i]) { f = x[i]; break; }
        if (f != 1) continue;            /* normalised: first nonzero coordinate is 1 (skips 0) */
        if (quad(x) != 0) continue;
        if (P >= expectP + 16) { fprintf(stderr, "too many points\n"); return 1; }
        memcpy(pts[P], x, sizeof x); index_of[c] = P; P++;
    }
    printf("q=%d: %d points on Q(6,q) (expected %ld)\n", q, P, expectP);
    if (P != expectP) { fprintf(stderr, "point count mismatch\n"); return 1; }
    /* lines of H(q): for each pair of collinear points satisfying Tits' conditions */
    hsize = 4L * expectP + 1024;
    htab = malloc(sizeof(int) * hsize);
    for (long i = 0; i < hsize; i++) htab[i] = -1;
    int l[16], y[7];
    long npairs = 0;
    for (int a = 0; a < P; a++) {
        const int *xa = pts[a];
        for (int b = a + 1; b < P; b++) {
            const int *xb = pts[b];
            if (bil(xa, xb) != 0) continue;
            if (!tits(xa, xb)) continue;
            npairs++;
            l[0] = a; l[1] = b;
            int n = 2;
            for (int lam = 1; lam < q; lam++) {
                for (int i = 0; i < 7; i++) y[i] = ADD(xa[i], MUL(lam, xb[i]));
                normalise(y);
                int id = index_of[code_of(y)];
                if (id < 0) { fprintf(stderr, "point on line not on quadric?\n"); return 1; }
                l[n++] = id;
            }
            add_line(l);
        }
        if (a % 5000 == 4999) { fprintf(stderr, "  %d/%d points processed, %d lines so far\n", a + 1, P, nlines); }
    }
    printf("%d hexagon lines (expected %d) from %ld Tits pairs, %.1fs\n", nlines, P, npairs, (double)clock() / CLOCKS_PER_SEC - t0);
    if (nlines != P) { fprintf(stderr, "line count mismatch\n"); return 1; }
    /* sort lines lexicographically (by sorted point tuple) for determinism */
    {
        int *idx = malloc(sizeof(int) * nlines);
        for (int i = 0; i < nlines; i++) idx[i] = i;
        qsort(idx, nlines, sizeof(int), cmp_line_idx);
        int *sorted = malloc(sizeof(int) * (long)nlines * (q + 1));
        for (int i = 0; i < nlines; i++) memcpy(sorted + (long)i * (q + 1), lines + (long)idx[i] * (q + 1), sizeof(int) * (q + 1));
        free(lines); lines = sorted; free(idx);
    }
    /* incidence graph */
    int NV = 2 * P, D = q + 1;
    int *adj = malloc(sizeof(int) * (long)NV * D), *deg = calloc(NV, sizeof(int));
    for (int j = 0; j < nlines; j++)
        for (int k = 0; k <= q; k++) {
            int p = lines[(long)j * (q + 1) + k], L = P + j;
            if (deg[p] >= D || deg[L] >= D) { fprintf(stderr, "degree overflow at point %d / line %d\n", p, j); return 1; }
            adj[(long)p * D + deg[p]++] = L; adj[(long)L * D + deg[L]++] = p;
        }
    for (int v = 0; v < NV; v++) if (deg[v] != D) { fprintf(stderr, "vertex %d has degree %d\n", v, deg[v]); return 1; }
    printf("incidence graph: %d vertices, %d-regular, %d edges\n", NV, D, NV * D / 2);
    FILE *f = fopen(out, "w");
    if (!f) { perror(out); return 1; }
    for (int p = 0; p < P; p++) for (int k = 0; k < D; k++) fprintf(f, "%d %d\n", p, adj[(long)p * D + k]);
    fclose(f);
    printf("wrote %s\n", out);
    if (ptsfile) {
        FILE *g = fopen(ptsfile, "w");
        for (int p = 0; p < P; p++) { fprintf(g, "%d", p); for (int i = 0; i < 7; i++) fprintf(g, " %d", pts[p][i]); fprintf(g, "\n"); }
        fclose(g);
    }
    /* BFS checks */
    int *dist = malloc(sizeof(int) * NV), *par = malloc(sizeof(int) * NV), *queue = malloc(sizeof(int) * NV);
    long expect[7]; expect[0] = 1; expect[1] = D; for (int i = 2; i <= 6; i++) expect[i] = expect[i - 1] * q;
    expect[6] = 1; for (int i = 0; i < 5; i++) expect[6] *= q;   /* q^5 */
    int nsrc = all ? NV : 20, bad = 0, girth = 1 << 30, maxecc = 0, uniq_ok = 1, geo6_ok = 1;
    for (int s_i = 0; s_i < nsrc; s_i++) {
        int s = all ? s_i : (s_i % 2 ? P + (s_i * 7919) % P : (s_i * 7919) % P);
        for (int v = 0; v < NV; v++) dist[v] = -1;
        int *npaths = queue;   /* reuse: number of geodesics, computed in a second pass */
        int qh = 0, qt = 0;
        static int *qq; if (!qq) qq = malloc(sizeof(int) * NV);
        dist[s] = 0; par[s] = -1; qq[qt++] = s;
        while (qh < qt) {
            int u = qq[qh++];
            for (int k = 0; k < D; k++) {
                int w = adj[(long)u * D + k];
                if (dist[w] < 0) { dist[w] = dist[u] + 1; par[w] = u; qq[qt++] = w; }
                else if (w != par[u] && dist[w] >= dist[u]) { int c = dist[u] + dist[w] + 1; if (c < girth) girth = c; }
            }
        }
        long cnt[8] = {0};
        for (int v = 0; v < NV; v++) { if (dist[v] < 0 || dist[v] > 7) { bad++; continue; } cnt[dist[v]]++; if (dist[v] > maxecc) maxecc = dist[v]; }
        for (int d = 0; d <= 6; d++) if (cnt[d] != expect[d]) bad++;
        if (cnt[7]) bad++;
        /* geodesic counts (unique for d <= 5, exactly q+1 at distance 6) */
        for (int v = 0; v < NV; v++) npaths[v] = 0;
        npaths[s] = 1;
        for (int i = 0; i < qt; i++) {
            int u = qq[i];
            for (int k = 0; k < D; k++) { int w = adj[(long)u * D + k]; if (dist[w] == dist[u] + 1) npaths[w] += npaths[u]; }
        }
        for (int v = 0; v < NV; v++) {
            if (dist[v] >= 1 && dist[v] <= 5 && npaths[v] != 1) uniq_ok = 0;
            if (dist[v] == 6 && npaths[v] != D) geo6_ok = 0;
        }
        if (s_i < 2) {
            printf("from vertex %d: distribution", s);
            for (int d = 0; d <= 6; d++) printf(" %ld", cnt[d]);
            printf(" (expected"); for (int d = 0; d <= 6; d++) printf(" %ld", expect[d]); printf(")\n");
        }
    }
    printf("BFS from %d vertices: distribution mismatches %d, max eccentricity %d, girth %d, unique geodesics d<=5: %s, %d geodesics at d=6: %s\n",
           nsrc, bad, maxecc, girth, uniq_ok ? "yes" : "NO", D, geo6_ok ? "yes" : "NO");
    printf("total time %.1fs\n", (double)clock() / CLOCKS_PER_SEC - t0);
    return bad || girth != 12 || !uniq_ok || !geo6_ok;
}
