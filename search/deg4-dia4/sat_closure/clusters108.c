/* clusters108.c: C port of the chaining-only cluster growth of clusters108.py
 * (Growth with --kmax=0), same windows, bounds and state semantics, for the complete
 * enumeration of clusters at n = 108, cost <= 12.
 *
 * A cluster is a chaining-closed multiset of gadget windows (every half-edge demand a(t -> m) >
 * a(m -> t), a(t -> t) odd, is met inside the multiset) grown from its minimum window (the seed)
 * by adding windows of id >= seed; slots whose complement vector is provided by no window of the
 * table are "declared" inexact (at most 2 X per cluster, X = excess budget); all other unmatched
 * slots are left to the combine stage. Bounds (necessary conditions for the whole member) as in
 * Growth.bounds_ok: cost, size, slots <= 32, deficit with affordable quads, excess used, chaining
 * cost bound, per-type window bound, demanded deficit, slot-complement bound with tolerance.
 *
 * Difference from the Python (option --seedfilter restores the Python behaviour for validation):
 * the Python decides "complement provided by no window" over the windows of id >= seed only, so a
 * slot whose complements all have smaller ids than the seed was declared (or killed the state)
 * although another cluster of the member can provide it; here the test is over the whole table.
 *
 * Usage: clusters108 WINDOWS.txt N CMAX8 X MAXSIZE OUT [--seedfilter] [--text] [--seeds=lo,hi]
 *        [--log=FILE] [--maxsec=S] [--progress=K]
 * Output: one record per cluster. Binary: u8 nwin, nwin x (u16 id, u8 cnt), u8 ndecl, u8 nres,
 * then ndecl x (u8 j, u8 f0, u8 f1). Text (--text): "W id:cnt ... D ndecl R nres".
 * Windows file: export_windows.py. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NT 20
#define MAXW 4096
#define MAXSZ 40
#define MAXD 16
#define BIG 1000000000

typedef struct {
    int id, t, nx, quad, x[3], cost8, deficit, nh, hm[4], hk[4], ns, sj[4], sf0[4], sf1[4], ex[6], hh[3];
} Win;

static Win win[MAXW];
static int NW;
static int DEG[NT], DEFI[NT], HC[NT][NT];
static int minw[NT][NT], maxmult[NT][NT];
static int cheapest_quad8 = BIG;
static double minsc[3][16][16];   /* min cost per slot of a window providing (j, b) */
static int maxprov[3][16][16];    /* max window id providing (j, b), -1 if none */
static int *byhalf[NT][NT];       /* windows of type m carrying t, ascending ids */
static int nbyhalf[NT][NT];

static int N, Dtot, cmax8, X, maxsize, seedfilter, textout;
static FILE *out, *logf;

typedef struct {
    uint8_t n;
    uint16_t id[MAXSZ];
    uint8_t cnt[MAXSZ];
    uint8_t nd;                    /* distinct declared vectors */
    uint8_t dj[MAXD], df0[MAXD], df1[MAXD], dc[MAXD];
} State;

#define MAXVEC 192
typedef struct {
    int a[NT][NT];
    int sv[3][16][16];
    int nvec[3];
    uint8_t vf0[3][MAXVEC], vf1[3][MAXVEC];
    int slots[3], ex[6], hh[3], deficit, cost8, G, quads;
} Ana;

/* ------------------------------------------------------------------ hash set of 128-bit keys */
typedef struct { uint64_t *k1, *k2; size_t cap, n; } HSet;

static inline uint64_t mix64(uint64_t z) {
    z ^= z >> 30; z *= 0xbf58476d1ce4e5b9ULL; z ^= z >> 27; z *= 0x94d049bb133111ebULL; z ^= z >> 31;
    return z;
}

static void hs_init(HSet *h, size_t cap) {
    h->cap = cap; h->n = 0;
    h->k1 = calloc(cap, sizeof(uint64_t)); h->k2 = calloc(cap, sizeof(uint64_t));
    if (!h->k1 || !h->k2) { fprintf(stderr, "out of memory (hash set)\n"); exit(2); }
}

static void hs_free(HSet *h) { free(h->k1); free(h->k2); h->k1 = h->k2 = NULL; h->cap = h->n = 0; }

static int hs_insert_raw(HSet *h, uint64_t a, uint64_t b) {
    size_t mask = h->cap - 1, i = (size_t)(a ^ (b * 0x9E3779B97F4A7C15ULL)) & mask;
    for (;;) {
        if (h->k1[i] == 0 && h->k2[i] == 0) { h->k1[i] = a; h->k2[i] = b; h->n++; return 1; }
        if (h->k1[i] == a && h->k2[i] == b) return 0;
        i = (i + 1) & mask;
    }
}

static void hs_grow(HSet *h) {
    HSet g; hs_init(&g, h->cap * 2);
    for (size_t i = 0; i < h->cap; i++)
        if (h->k1[i] || h->k2[i]) hs_insert_raw(&g, h->k1[i], h->k2[i]);
    hs_free(h); *h = g;
}

static int hs_insert(HSet *h, uint64_t a, uint64_t b) {
    if (a == 0 && b == 0) a = 1;
    if (h->n * 5 >= h->cap * 3) hs_grow(h);
    return hs_insert_raw(h, a, b);
}

static void state_hash(const State *s, uint64_t *a, uint64_t *b) {
    uint64_t h1 = 0x243F6A8885A308D3ULL, h2 = 0x13198A2E03707344ULL;
    for (int i = 0; i < s->n; i++) {
        uint64_t v = ((uint64_t)s->id[i] << 8) | s->cnt[i];
        h1 = mix64(h1 ^ (v * 0x9E3779B97F4A7C15ULL + (uint64_t)i));
        h2 = mix64(h2 + v + 0xA54FF53A5F1D36F1ULL) * 0xD6E8FEB86659FD93ULL ^ (h2 >> 29);
    }
    h1 = mix64(h1 ^ 0xFFFFFFFFULL); h2 = mix64(h2 + 0x510E527FADE682D1ULL);
    for (int i = 0; i < s->nd; i++) {
        uint64_t v = ((uint64_t)s->dj[i] << 24) | ((uint64_t)s->df0[i] << 16) | ((uint64_t)s->df1[i] << 8) | s->dc[i];
        h1 = mix64(h1 ^ (v * 0x9E3779B97F4A7C15ULL + 77 + (uint64_t)i));
        h2 = mix64(h2 + v + 0x9B05688C2B3E6C1FULL) * 0xD6E8FEB86659FD93ULL ^ (h2 >> 29);
    }
    *a = h1; *b = h2;
}

/* ------------------------------------------------------------------ windows */
static void load_windows(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }
    int nt;
    if (fscanf(f, "NT %d", &nt) != 1 || nt != NT) { fprintf(stderr, "bad NT\n"); exit(1); }
    for (int t = 0; t < NT; t++) {
        int tt, x0, x1, x2;
        if (fscanf(f, "%d %d %d %d %d %d", &tt, &x0, &x1, &x2, &DEG[t], &DEFI[t]) != 6 || tt != t) { fprintf(stderr, "bad type line\n"); exit(1); }
    }
    for (int t = 0; t < NT; t++)
        for (int m = 0; m < NT; m++)
            if (fscanf(f, "%d", &HC[t][m]) != 1) { fprintf(stderr, "bad HC\n"); exit(1); }
    if (fscanf(f, " NW %d", &NW) != 1 || NW > MAXW) { fprintf(stderr, "bad NW\n"); exit(1); }
    for (int k = 0; k < NW; k++) {
        Win *w = &win[k];
        if (fscanf(f, "%d %d %d %d %d %d %d %d %d %d", &w->id, &w->t, &w->nx, &w->quad, &w->x[0], &w->x[1], &w->x[2], &w->cost8, &w->deficit, &w->nh) != 10 || w->id != k) {
            fprintf(stderr, "bad window %d\n", k); exit(1);
        }
        for (int i = 0; i < w->nh; i++) if (fscanf(f, "%d %d", &w->hm[i], &w->hk[i]) != 2) exit(1);
        if (fscanf(f, "%d", &w->ns) != 1) exit(1);
        for (int i = 0; i < w->ns; i++) if (fscanf(f, "%d %d %d", &w->sj[i], &w->sf0[i], &w->sf1[i]) != 3) exit(1);
        for (int i = 0; i < 6; i++) if (fscanf(f, "%d", &w->ex[i]) != 1) exit(1);
        for (int i = 0; i < 3; i++) if (fscanf(f, "%d", &w->hh[i]) != 1) exit(1);
    }
    fclose(f);
    /* derived tables */
    for (int t = 0; t < NT; t++) for (int m = 0; m < NT; m++) { minw[t][m] = BIG; maxmult[t][m] = 0; nbyhalf[t][m] = 0; }
    for (int j = 0; j < 3; j++) for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) { minsc[j][a][b] = BIG; maxprov[j][a][b] = -1; }
    for (int k = 0; k < NW; k++) {
        Win *w = &win[k];
        if (w->quad) { if (w->cost8 < cheapest_quad8) cheapest_quad8 = w->cost8; }
        else {
            for (int i = 0; i < w->nh; i++) {
                int m = w->hm[i];
                if (w->cost8 < minw[w->t][m]) minw[w->t][m] = w->cost8;
                if (w->hk[i] > maxmult[w->t][m]) maxmult[w->t][m] = w->hk[i];
                nbyhalf[w->t][m]++;
            }
        }
        if (w->nx > 0) {
            double c = (double)w->cost8 / w->nx;
            for (int i = 0; i < w->ns; i++) {
                int j = w->sj[i], a = w->sf0[i], b = w->sf1[i];
                if (c < minsc[j][a][b]) minsc[j][a][b] = c;
                if (k > maxprov[j][a][b]) maxprov[j][a][b] = k;
            }
        }
    }
    for (int t = 0; t < NT; t++) for (int m = 0; m < NT; m++) { byhalf[t][m] = malloc(sizeof(int) * (nbyhalf[t][m] + 1)); nbyhalf[t][m] = 0; }
    for (int k = 0; k < NW; k++) {
        Win *w = &win[k];
        if (w->quad) continue;
        for (int i = 0; i < w->nh; i++) byhalf[w->t][w->hm[i]][nbyhalf[w->t][w->hm[i]]++] = k;   /* ascending ids */
    }
}

/* ------------------------------------------------------------------ analysis */
static inline void ana_touch(Ana *A, int j, int f0, int f1) {
    /* register a vector in the per-block list (duplicates allowed; consumers use stamps) */
    if (A->nvec[j] >= MAXVEC) { fprintf(stderr, "vector list overflow\n"); exit(2); }
    A->vf0[j][A->nvec[j]] = f0; A->vf1[j][A->nvec[j]] = f1; A->nvec[j]++;
}

static void ana_add(Ana *A, const Win *w, int cnt) {
    if (!w->quad)
        for (int i = 0; i < w->nh; i++) A->a[w->t][w->hm[i]] += w->hk[i] * cnt;
    for (int i = 0; i < w->ns; i++) {
        A->sv[w->sj[i]][w->sf0[i]][w->sf1[i]] += cnt;
        if (cnt > 0) ana_touch(A, w->sj[i], w->sf0[i], w->sf1[i]);
    }
    for (int j = 0; j < 3; j++) { A->slots[j] += cnt * w->x[j]; A->hh[j] += cnt * w->hh[j]; }
    for (int k = 0; k < 6; k++) A->ex[k] += cnt * w->ex[k];
    A->deficit += cnt * w->deficit; A->cost8 += cnt * w->cost8; A->G += cnt;
    if (w->quad) A->quads += cnt;
}

static void analyse(const State *s, Ana *A) {
    memset(A, 0, sizeof(*A));
    for (int i = 0; i < s->n; i++) ana_add(A, &win[s->id[i]], s->cnt[i]);
    for (int i = 0; i < s->nd; i++) A->sv[s->dj[i]][s->df0[i]][s->df1[i]] -= s->dc[i];
}

static inline int ndecl_of(const State *s) { int n = 0; for (int i = 0; i < s->nd; i++) n += s->dc[i]; return n; }

static inline double partner_cost(int j, int f0, int f1) {
    int b0 = 4 - f0, b1 = 4 - f1;
    if (b0 < 0 || b1 < 0) return BIG;
    return minsc[j][b0][b1];
}

static int cmp_desc(const void *p, const void *q) {
    double a = *(const double *)p, b = *(const double *)q;
    return (a < b) - (a > b);
}

static int bounds_ok(const Ana *A, int used) {
    int cost8 = A->cost8;
    if (cost8 > cmax8 || A->G > maxsize) return 0;
    if (A->slots[0] > 32 || A->slots[1] > 32 || A->slots[2] > 32) return 0;
    int qmax = cheapest_quad8 < BIG ? (cmax8 - cost8) / cheapest_quad8 : 0;
    if (A->deficit > Dtot + qmax) return 0;
    int budget = cmax8 - cost8;
    if (used > Dtot - (cost8 + 7) / 8) return 0;
    long long chain8 = 0, need_def = 0, type_bound = 0;
    for (int m = 0; m < NT; m++) {
        int D = A->a[m][m] % 2, nwin = D;
        long long minc = D ? minw[m][m] : BIG;
        if (D) chain8 += HC[m][m];
        for (int t = 0; t < NT; t++) {
            if (t == m) continue;
            int d = A->a[t][m] - A->a[m][t];
            if (d > 0) {
                chain8 += d * HC[t][m];
                int mm = maxmult[m][t];
                if (mm == 0) return 0;
                int need = (d + mm - 1) / mm;
                if (need > nwin) nwin = need;
                if (minw[m][t] < minc) minc = minw[m][t];
                D += d;
            }
        }
        if (D) {
            int need = (D + DEG[m] - 1) / DEG[m];
            if (need > nwin) nwin = need;
            type_bound += nwin * minc;
            need_def += DEFI[m] * nwin;
        }
    }
    if (chain8 > budget || type_bound > budget) return 0;
    if (need_def > Dtot - A->deficit + qmax) return 0;
    int t1 = X - used, t2 = Dtot - (cost8 + 7) / 8 - used;
    int tol = 2 * (t1 < t2 ? t1 : t2);
    if (tol < 0) return 0;
    double costs[256]; int nc = 0;
    static int stamp[3][16][16]; static int stampv = 0;
    stampv++;
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < A->nvec[j]; k++) {
            int f0 = A->vf0[j][k], f1 = A->vf1[j][k];
            if (stamp[j][f0][f1] == stampv) continue;
            stamp[j][f0][f1] = stampv;
            if (f0 == 2 && f1 == 2) continue;
            int cnt = A->sv[j][f0][f1];
            if (cnt <= 0) continue;
            int other = (4 - f0 >= 0 && 4 - f1 >= 0) ? A->sv[j][4 - f0][4 - f1] : 0;
            if (cnt > other) {
                double pc = partner_cost(j, f0, f1);
                for (int r = 0; r < cnt - other && nc < 256; r++) costs[nc++] = pc;
            }
        }
    }
    if (nc > tol) {
        qsort(costs, nc, sizeof(double), cmp_desc);
        double s = 0;
        for (int k = tol; k < nc; k++) s += costs[k];
        if (s > budget + 1e-9) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ growth */
static State *stack; static size_t stack_n, stack_cap;
static HSet visited;
static long long st_states, st_clusters, st_dead_end, st_pruned_child, st_children, st_revisit, st_dead_seed, st_seeds;
static long long st_clusters_decl[MAXD];
static long long st_by_size[MAXSZ + 1];
static time_t t_start; static int maxsec = 0; static int stop_flag = 0; static long long maxstates = 0; static int incomplete_seeds = 0;

static void push_state(const State *s) {
    if (stack_n == stack_cap) {
        stack_cap = stack_cap ? stack_cap * 2 : 1024;
        stack = realloc(stack, stack_cap * sizeof(State));
        if (!stack) { fprintf(stderr, "out of memory (stack)\n"); exit(2); }
    }
    stack[stack_n++] = *s;
}

static void state_add_window(State *s, int wid) {
    int i = 0;
    while (i < s->n && s->id[i] < wid) i++;
    if (i < s->n && s->id[i] == wid) { s->cnt[i]++; return; }
    if (s->n >= MAXSZ) { fprintf(stderr, "state overflow\n"); exit(2); }
    for (int k = s->n; k > i; k--) { s->id[k] = s->id[k - 1]; s->cnt[k] = s->cnt[k - 1]; }
    s->id[i] = wid; s->cnt[i] = 1; s->n++;
}

static void state_add_decl(State *s, int j, int f0, int f1) {
    int key = (j << 8) | (f0 << 4) | f1, i = 0;
    while (i < s->nd && ((s->dj[i] << 8) | (s->df0[i] << 4) | s->df1[i]) < key) i++;
    if (i < s->nd && s->dj[i] == j && s->df0[i] == f0 && s->df1[i] == f1) { s->dc[i]++; return; }
    if (s->nd >= MAXD) { fprintf(stderr, "decl overflow\n"); exit(2); }
    for (int k = s->nd; k > i; k--) { s->dj[k] = s->dj[k - 1]; s->df0[k] = s->df0[k - 1]; s->df1[k] = s->df1[k - 1]; s->dc[k] = s->dc[k - 1]; }
    s->dj[i] = j; s->df0[i] = f0; s->df1[i] = f1; s->dc[i] = 1; s->nd++;
}

static void record(const State *s, const Ana *A, int nres) {
    int ndecl = ndecl_of(s);
    st_clusters++;
    if (ndecl < MAXD) st_clusters_decl[ndecl]++;
    if (A->G <= MAXSZ) st_by_size[A->G]++;
    if (textout) {
        fprintf(out, "W");
        for (int i = 0; i < s->n; i++) fprintf(out, " %d:%d", s->id[i], s->cnt[i]);
        fprintf(out, " D %d R %d\n", ndecl, nres + ndecl);
    } else {
        uint8_t hdr = s->n;
        fwrite(&hdr, 1, 1, out);
        for (int i = 0; i < s->n; i++) { uint16_t id = s->id[i]; uint8_t c = s->cnt[i]; fwrite(&id, 2, 1, out); fwrite(&c, 1, 1, out); }
        uint8_t nd = ndecl, nr = (uint8_t)(nres + ndecl);
        fwrite(&nd, 1, 1, out); fwrite(&nr, 1, 1, out);
        for (int i = 0; i < s->nd; i++) for (int c = 0; c < s->dc[i]; c++) { uint8_t v[3] = { s->dj[i], s->df0[i], s->df1[i] }; fwrite(v, 1, 3, out); }
    }
}

static int provider_exists(int seed, int j, int b0, int b1) {
    if (b0 < 0 || b1 < 0) return 0;
    int mp = maxprov[j][b0][b1];
    if (mp < 0) return 0;
    return seedfilter ? (mp >= seed) : 1;
}

static void grow(int seed) {
    State s0; memset(&s0, 0, sizeof(s0));
    s0.n = 1; s0.id[0] = seed; s0.cnt[0] = 1;
    Ana A;
    analyse(&s0, &A);
    if (!bounds_ok(&A, 0)) { st_dead_seed++; return; }
    hs_init(&visited, 1 << 12);
    uint64_t h1, h2;
    state_hash(&s0, &h1, &h2); hs_insert(&visited, h1, h2);
    stack_n = 0; push_state(&s0);
    /* per-seed start offsets into byhalf lists: first index with id >= seed */
    static int start[NT][NT];
    for (int t = 0; t < NT; t++) for (int m = 0; m < NT; m++) {
        int lo = 0, hi = nbyhalf[t][m];
        while (lo < hi) { int mid = (lo + hi) / 2; if (byhalf[t][m][mid] < seed) lo = mid + 1; else hi = mid; }
        start[t][m] = lo;
    }
    static int stamp[3][16][16]; static int stampv = 0;
    long long seed_states0 = st_states;
    while (stack_n) {
        State s = stack[--stack_n];
        st_states++;
        if ((st_states & 1023) == 0 && maxsec && (int)(time(NULL) - t_start) > maxsec) { stop_flag = 1; break; }
        if (maxstates && st_states - seed_states0 > maxstates) { incomplete_seeds++; fprintf(logf, "seed %d: state cap %lld exceeded, seed INCOMPLETE\n", seed, maxstates); fflush(logf); break; }
        analyse(&s, &A);
        int ndecl = ndecl_of(&s);
        /* demands: (kind, key, noptions); kind 0 = half-edge (m, t): windows of type m carrying t;
           kind 1 = slot (j, f0, f1): declare only */
        int best_kind = -1, best_n = BIG, best_m = 0, best_t = 0, best_j = 0, best_f0 = 0, best_f1 = 0;
        for (int t = 0; t < NT; t++) {
            if (A.a[t][t] % 2) {
                int n = nbyhalf[t][t] - start[t][t];
                if (n < best_n) { best_n = n; best_kind = 0; best_m = t; best_t = t; }
            }
            for (int m = 0; m < NT; m++) {
                if (m != t && A.a[t][m] > A.a[m][t]) {
                    int n = nbyhalf[m][t] - start[m][t];
                    if (n < best_n) { best_n = n; best_kind = 0; best_m = m; best_t = t; }
                }
            }
        }
        int nres = 0;
        stampv++;
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < A.nvec[j]; k++) {
                int f0 = A.vf0[j][k], f1 = A.vf1[j][k];
                if (stamp[j][f0][f1] == stampv) continue;
                stamp[j][f0][f1] = stampv;
                if (f0 == 2 && f1 == 2) continue;
                int cnt = A.sv[j][f0][f1];
                if (cnt <= 0) continue;
                int c0 = 4 - f0, c1 = 4 - f1;
                int other = (c0 >= 0 && c1 >= 0) ? A.sv[j][c0][c1] : 0;
                if (cnt <= other) continue;
                if (provider_exists(seed, j, c0, c1)) { nres += cnt - other; continue; }
                int n = (ndecl < 2 * X) ? 1 : 0;
                if (n < best_n) { best_n = n; best_kind = 1; best_j = j; best_f0 = f0; best_f1 = f1; }
            }
        }
        if (best_kind < 0) { record(&s, &A, nres); continue; }
        if (best_n == 0) { st_dead_end++; continue; }
        if (best_kind == 1) {
            st_children++;
            State c = s; state_add_decl(&c, best_j, best_f0, best_f1);
            Ana B = A; B.sv[best_j][best_f0][best_f1] -= 1;
            if (!bounds_ok(&B, (ndecl + 2) / 2)) { st_pruned_child++; continue; }
            state_hash(&c, &h1, &h2);
            if (!hs_insert(&visited, h1, h2)) { st_revisit++; continue; }
            push_state(&c);
            continue;
        }
        int *lst = byhalf[best_m][best_t];
        for (int k = start[best_m][best_t]; k < nbyhalf[best_m][best_t]; k++) {
            int wid = lst[k];
            st_children++;
            Ana B = A; ana_add(&B, &win[wid], 1);
            if (!bounds_ok(&B, (ndecl + 1) / 2)) { st_pruned_child++; continue; }
            State c = s; state_add_window(&c, wid);
            state_hash(&c, &h1, &h2);
            if (!hs_insert(&visited, h1, h2)) { st_revisit++; continue; }
            push_state(&c);
        }
    }
    hs_free(&visited);
}

static void print_stats(FILE *f, const char *tag) {
    fprintf(f, "%s seeds %lld states %lld clusters %lld dead_end %lld pruned_child %lld children %lld revisit %lld dead_seed %lld elapsed %ld s\n",
            tag, st_seeds, st_states, st_clusters, st_dead_end, st_pruned_child, st_children, st_revisit, st_dead_seed, (long)(time(NULL) - t_start));
    fprintf(f, "%s clusters by decl:", tag);
    for (int i = 0; i < MAXD; i++) if (st_clusters_decl[i]) fprintf(f, " %d:%lld", i, st_clusters_decl[i]);
    fprintf(f, "\n%s clusters by size:", tag);
    for (int i = 0; i <= MAXSZ; i++) if (st_by_size[i]) fprintf(f, " %d:%lld", i, st_by_size[i]);
    fprintf(f, "\n");
    fflush(f);
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s WINDOWS.txt N CMAX8 X MAXSIZE OUT [--seedfilter] [--text] [--seeds=lo,hi] [--log=FILE] [--maxsec=S] [--progress=K] [--maxstates=M]\n", argv[0]); return 1; }
    N = atoi(argv[2]); cmax8 = atoi(argv[3]); X = atoi(argv[4]); maxsize = atoi(argv[5]);
    Dtot = 3 * N - 96;
    int seed_lo = 0, seed_hi = -1, progress = 100;
    logf = stderr;
    for (int i = 7; i < argc; i++) {
        if (!strcmp(argv[i], "--seedfilter")) seedfilter = 1;
        else if (!strcmp(argv[i], "--text")) textout = 1;
        else if (!strncmp(argv[i], "--seeds=", 8)) sscanf(argv[i] + 8, "%d,%d", &seed_lo, &seed_hi);
        else if (!strncmp(argv[i], "--log=", 6)) { logf = fopen(argv[i] + 6, "a"); if (!logf) { perror("log"); return 1; } }
        else if (!strncmp(argv[i], "--maxsec=", 9)) maxsec = atoi(argv[i] + 9);
        else if (!strncmp(argv[i], "--progress=", 11)) progress = atoi(argv[i] + 11);
        else if (!strncmp(argv[i], "--maxstates=", 12)) maxstates = atoll(argv[i] + 12);
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    load_windows(argv[1]);
    if (seed_hi < 0 || seed_hi > NW) seed_hi = NW;
    out = fopen(argv[6], textout ? "w" : "wb");
    if (!out) { perror(argv[6]); return 1; }
    t_start = time(NULL);
    fprintf(logf, "clusters108.c: N=%d Dtot=%d cmax8=%d X=%d maxsize=%d windows=%d seeds=[%d,%d) seedfilter=%d\n", N, Dtot, cmax8, X, maxsize, NW, seed_lo, seed_hi, seedfilter);
    int last_seed = seed_lo;
    for (int seed = seed_lo; seed < seed_hi; seed++) {
        grow(seed);
        st_seeds++;
        last_seed = seed;
        if (stop_flag) { fprintf(logf, "time limit reached during seed %d\n", seed); break; }
        if (progress && (seed + 1) % progress == 0) { char tag[64]; snprintf(tag, sizeof tag, "seed %d:", seed); print_stats(logf, tag); }
    }
    fclose(out);
    char tag[96]; snprintf(tag, sizeof tag, "done (last seed %d%s, %d seeds over the state cap):", last_seed, stop_flag ? ", INCOMPLETE" : "", incomplete_seeds);
    print_stats(logf, tag);
    return (stop_flag || incomplete_seeds) ? 3 : 0;
}
