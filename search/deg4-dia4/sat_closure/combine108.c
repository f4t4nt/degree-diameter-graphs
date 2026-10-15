/* combine108.c: shapes whose gadget window multiset is m1 c1 + m2 c2 for clusters c1, c2 of a
 * cluster file written by clusters108.c (binary format), i.e. the pair combine of
 * pairs108.py in C for large cluster lists, plus the single-cluster case m c.
 *
 * Shape conditions (all necessary, as in pairs108.try_shape): total deficit D = 3N - 96 (quads
 * -1), total cost in the accepted set, G = |W| with u = (N - G) / 2 units integral, slots
 * (32 - 2u) per block, slot imbalance ||sum sigma||_1 <= 2 e where e = D - cost is the excess,
 * declared slots <= 2 e, six rows >= 0, hub-to-hub sum + 2u >= 16 per block pair, (2,2)-slot
 * parity per block at e = 0. The partner of (c1, m1, m2, cost) has a determined (deficit, cost,
 * G, slots) and a slot imbalance sigma2 with ||m1 sigma1 + m2 sigma2||_1 <= 2 e: clusters are
 * indexed by (deficit, cost, G, slots, sigma) for e = 0 and, for e > 0, by the same key with sigma
 * restricted to each of 5 coordinate groups (a residual of L1 norm <= 4 is nonzero on at most 4
 * coordinates, so at least one group is matched exactly). Survivors are written as window
 * multisets ("u G cost8 | wid:cnt ...") for shapes_from_c.py (canonical form under block
 * relabelling, slot-vector matching, realisation; the parts of pairs108/windows108 kept in
 * Python).
 *
 * Usage: combine108 WINDOWS.txt N ACCEPT8(comma list) MAXM OUT.txt CLUSTERS.bin [CLUSTERS2.bin ...]
 *        [--log=FILE] */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "clusters_io.h"

#define NGROUP 5
static int kgroup[MAXK];

static int N, Dtot, maxm, nacc; static int accept8[8]; static int cmax8;
static FILE *out;
static long long st_single_cand, st_pair_cand, st_survive, st_fail_imb, st_fail_decl, st_fail_rows, st_fail_hh, st_fail_par, st_lookups;

/* ------------------------------------------------------------------ hash index */
typedef struct { uint64_t *key; int32_t *val; int32_t *next; size_t cap, n; } HMap;   /* key -> chain of cluster ids */
static HMap tab[NGROUP + 1];     /* NGROUP: full sigma; 0..NGROUP-1: restricted */

static void hm_init(HMap *h, size_t cap) {
    h->cap = cap; h->n = 0;
    h->key = calloc(cap, sizeof(uint64_t)); h->val = malloc(cap * sizeof(int32_t)); h->next = malloc(ncl * sizeof(int32_t));
    if (!h->key || !h->val || !h->next) { fprintf(stderr, "oom (index)\n"); exit(2); }
    for (size_t i = 0; i < cap; i++) h->val[i] = -1;
}

static void hm_insert(HMap *h, uint64_t key, int32_t id) {
    size_t mask = h->cap - 1, i = key & mask;
    for (;;) {
        if (h->val[i] < 0) { h->key[i] = key; h->val[i] = id; h->next[id] = -1; h->n++; return; }
        if (h->key[i] == key) { h->next[id] = h->val[i]; h->val[i] = id; return; }
        i = (i + 1) & mask;
    }
}

static int32_t hm_find(const HMap *h, uint64_t key) {
    size_t mask = h->cap - 1, i = key & mask;
    for (;;) {
        if (h->val[i] < 0) return -1;
        if (h->key[i] == key) return h->val[i];
        i = (i + 1) & mask;
    }
}

static uint64_t key_of(int deficit, int cost8, int G, const int *slots, const int *sig, int group) {
    uint64_t h = 0x9E3779B97F4A7C15ULL;
    h = mix64(h ^ (uint64_t)(deficit + 1000)); h = mix64(h ^ (uint64_t)(cost8 + 7));
    h = mix64(h ^ (uint64_t)G); h = mix64(h ^ ((uint64_t)slots[0] << 20 | (uint64_t)slots[1] << 10 | (uint64_t)slots[2]));
    for (int k = 0; k < NK; k++) {
        if (group < NGROUP && kgroup[k] != group) continue;
        if (sig[k]) h = mix64(h ^ ((uint64_t)(k + 1) << 32 | (uint64_t)(uint32_t)(sig[k] + 1000)));
    }
    return h | 1;
}

static void build_index(void) {
    size_t cap = 1; while (cap < ncl * 2) cap <<= 1;
    for (int g = 0; g <= NGROUP; g++) hm_init(&tab[g], cap);
    int sig[MAXK];
    for (size_t i = 0; i < ncl; i++) {
        const Cl *c = &cl[i];
        cluster_sig(c, sig);
        int slots[3] = { c->slots[0], c->slots[1], c->slots[2] };
        for (int g = 0; g <= NGROUP; g++) hm_insert(&tab[g], key_of(c->deficit, c->cost8, c->G, slots, sig, g), (int32_t)i);
    }
}

/* ------------------------------------------------------------------ shape checks */
static int try_shape(const int *ids, const int *mult, int nparts) {
    int cost8 = 0, deficit = 0, G = 0, slots[3] = {0, 0, 0}, decl = 0, ex[6] = {0}, hh[3] = {0, 0, 0}, sv22[3] = {0, 0, 0};
    static int sig[MAXK]; memset(sig, 0, sizeof(int) * NK);
    for (int p = 0; p < nparts; p++) {
        const Cl *c = &cl[ids[p]]; int m = mult[p];
        cost8 += m * c->cost8; deficit += m * c->deficit; G += m * c->G; decl += m * c->decl;
        for (int j = 0; j < 3; j++) { slots[j] += m * c->slots[j]; hh[j] += m * c->hh[j]; sv22[j] += m * c->sv22[j]; }
        for (int k = 0; k < 6; k++) ex[k] += m * c->ex[k];
        for (int s = 0; s < c->nsig; s++) sig[(uint8_t)sarena[c->soff + 2 * s]] += m * sarena[c->soff + 2 * s + 1];
    }
    int ok = 0; for (int a = 0; a < nacc; a++) if (accept8[a] == cost8) ok = 1;
    if (!ok || deficit != Dtot || (N - G) % 2) return 0;
    int u = (N - G) / 2, g = 32 - 2 * u;
    if (slots[0] != g || slots[1] != g || slots[2] != g) return 0;
    int excess = Dtot - cost8 / 8, l1 = 0;
    for (int k = 0; k < NK; k++) l1 += sig[k] < 0 ? -sig[k] : sig[k];
    if (l1 > 2 * excess) { st_fail_imb++; return 0; }
    if (decl > 2 * excess) { st_fail_decl++; return 0; }
    for (int k = 0; k < 6; k++) if (ex[k] < 0) { st_fail_rows++; return 0; }
    for (int j = 0; j < 3; j++) if (hh[j] + 2 * u < 16) { st_fail_hh++; return 0; }
    if (excess == 0) for (int j = 0; j < 3; j++) if ((sv22[j] + 2 * u) % 2) { st_fail_par++; return 0; }
    st_survive++;
    /* window multiset */
    static int wc[MAXW]; static int touched[MAXW]; int nt = 0;
    for (int p = 0; p < nparts; p++) {
        const Cl *c = &cl[ids[p]];
        for (int i = 0; i < c->nwin; i++) { int wid = warena[c->woff + 2 * i]; if (!wc[wid]) touched[nt++] = wid; wc[wid] += mult[p] * warena[c->woff + 2 * i + 1]; }
    }
    fprintf(out, "%d %d %d %d |", u, G, cost8, nparts);
    for (int i = 0; i < nt; i++) fprintf(out, " %d:%d", touched[i], wc[touched[i]]);
    fprintf(out, " | parts");
    for (int p = 0; p < nparts; p++) fprintf(out, " %d*%d", mult[p], ids[p]);
    fprintf(out, "\n");
    for (int i = 0; i < nt; i++) wc[touched[i]] = 0;
    return 1;
}

static int divisible_target(const int *s1, int m1, int m2, int *t2, int group) {
    /* t2 = -m1 s1 / m2 on the coordinates of the group; 0 if not integral there */
    for (int k = 0; k < NK; k++) {
        if (group < NGROUP && kgroup[k] != group) { t2[k] = 0; continue; }
        int v = -m1 * s1[k];
        if (v % m2) return 0;
        t2[k] = v / m2;
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: %s WINDOWS.txt N ACCEPT8 MAXM OUT.txt CLUSTERS.bin [...] [--log=FILE]\n", argv[0]); return 1; }
    logf = stderr;
    for (int i = 6; i < argc; i++) if (!strncmp(argv[i], "--log=", 6)) { logf = fopen(argv[i] + 6, "a"); if (!logf) { perror("log"); return 1; } }
    N = atoi(argv[2]); Dtot = 3 * N - 96; maxm = atoi(argv[4]);
    { char *s = strdup(argv[3]), *p = strtok(s, ","); nacc = 0; cmax8 = 0; while (p) { accept8[nacc++] = atoi(p); if (accept8[nacc - 1] > cmax8) cmax8 = accept8[nacc - 1]; p = strtok(NULL, ","); } }
    time_t t0 = time(NULL);
    load_windows(argv[1]);
    for (int c = 0; c < NK; c++) kgroup[c] = c % NGROUP;
    int emax0 = 0; for (int a = 0; a < nacc; a++) if (Dtot - accept8[a] / 8 > emax0) emax0 = Dtot - accept8[a] / 8;
    for (int i = 6; i < argc; i++) if (strncmp(argv[i], "--", 2)) load_clusters(argv[i], cmax8, maxm == 0 ? 255 : 2 * emax0);
    fprintf(logf, "%zu clusters loaded, %ld s\n", ncl, (long)(time(NULL) - t0));
    if (maxm == 0) {
        /* statistics only: distribution by deficit, cost, size, balance, declared slots */
        static long long by_def[64], by_cost[128], by_size[64], by_decl[16], bal_by_def[64], by_l1[64];
        for (size_t i = 0; i < ncl; i++) {
            const Cl *c = &cl[i];
            by_def[c->deficit + 8]++; by_cost[c->cost8]++; by_size[c->G]++; by_decl[c->decl]++;
            int l1 = 0; for (int s = 0; s < c->nsig; s++) { int v = sarena[c->soff + 2 * s + 1]; l1 += v < 0 ? -v : v; }
            if (l1 == 0) bal_by_def[c->deficit + 8]++;
            by_l1[l1 > 63 ? 63 : l1]++;
        }
        fprintf(logf, "by deficit (count, self-balanced):"); for (int d = 0; d < 64; d++) if (by_def[d]) fprintf(logf, " %d:%lld/%lld", d - 8, by_def[d], bal_by_def[d]); fprintf(logf, "\n");
        fprintf(logf, "by cost8:"); for (int d = 0; d < 128; d++) if (by_cost[d]) fprintf(logf, " %d:%lld", d, by_cost[d]); fprintf(logf, "\n");
        fprintf(logf, "by size:"); for (int d = 0; d < 64; d++) if (by_size[d]) fprintf(logf, " %d:%lld", d, by_size[d]); fprintf(logf, "\n");
        fprintf(logf, "by decl:"); for (int d = 0; d < 16; d++) if (by_decl[d]) fprintf(logf, " %d:%lld", d, by_decl[d]); fprintf(logf, "\n");
        fprintf(logf, "by sigma L1:"); for (int d = 0; d < 64; d++) if (by_l1[d]) fprintf(logf, " %d:%lld", d, by_l1[d]); fprintf(logf, "\n");
        return 0;
    }
    out = fopen(argv[5], "w"); if (!out) { perror(argv[5]); return 1; }
    build_index();
    fprintf(logf, "index built, %ld s\n", (long)(time(NULL) - t0));
    int emax = Dtot - cmax8 / 8; for (int a = 0; a < nacc; a++) if (Dtot - accept8[a] / 8 > emax) emax = Dtot - accept8[a] / 8;
    /* singles and multiples m c */
    int s1[MAXK], t2[MAXK];
    for (size_t i = 0; i < ncl; i++) {
        const Cl *c = &cl[i];
        for (int m = 1; m <= maxm; m++) {
            if (m * c->G > N || m * c->cost8 > cmax8) break;
            if (c->deficit > 0 && m * c->deficit > Dtot) break;
            if (m * c->deficit != Dtot) continue;
            st_single_cand++;
            int id = (int)i; try_shape(&id, &m, 1);
        }
    }
    fprintf(logf, "singles done: candidates %lld survive %lld, %ld s\n", st_single_cand, st_survive, (long)(time(NULL) - t0));
    /* pairs m1 c1 + m2 c2, c1 < c2 */
    for (size_t i = 0; i < ncl; i++) {
        const Cl *c1 = &cl[i];
        cluster_sig(c1, s1);
        for (int m1 = 1; m1 <= maxm; m1++) {
            if (m1 * c1->G > N || m1 * c1->cost8 > cmax8) break;
            if (c1->deficit > 0 && m1 * c1->deficit > Dtot + 2) break;
            for (int m2 = 1; m2 <= maxm; m2++) {
                int d2m = Dtot - m1 * c1->deficit;
                if (d2m % m2) continue;
                int d2 = d2m / m2;
                for (int a = 0; a < nacc; a++) {
                    int c8m = accept8[a] - m1 * c1->cost8;
                    if (c8m < 0 || c8m % m2) continue;
                    int cst2 = c8m / m2;
                    int excess = Dtot - accept8[a] / 8;
                    for (int G = m1 * c1->G + m2; G <= N; G++) {
                        if ((N - G) % 2 || (G - m1 * c1->G) % m2) continue;
                        int G2 = (G - m1 * c1->G) / m2;
                        int g = 32 - (N - G), slots2[3], bad = 0;
                        for (int j = 0; j < 3; j++) { int v = g - m1 * c1->slots[j]; if (v < 0 || v % m2) { bad = 1; break; } slots2[j] = v / m2; }
                        if (bad) continue;
                        /* candidate partners: exact sigma (excess 0) or one of the 5 group tables */
                        int glo = excess == 0 ? NGROUP : 0, ghi = excess == 0 ? NGROUP : NGROUP - 1;
                        for (int grp = glo; grp <= ghi; grp++) {
                            if (!divisible_target(s1, m1, m2, t2, grp)) continue;
                            uint64_t key = key_of(d2, cst2, G2, slots2, t2, grp);
                            st_lookups++;
                            for (int32_t j2 = hm_find(&tab[grp], key); j2 >= 0; j2 = tab[grp].next[j2]) {
                                if ((size_t)j2 <= i) continue;        /* c1 < c2 */
                                const Cl *c2 = &cl[j2];
                                if (c2->deficit != d2 || c2->cost8 != cst2 || c2->G != G2 || c2->slots[0] != slots2[0] || c2->slots[1] != slots2[1] || c2->slots[2] != slots2[2]) continue;
                                if (grp < NGROUP) {
                                    /* a partner can be found through several clean groups: count it once, at the first clean group */
                                    int sig2[MAXK]; cluster_sig(c2, sig2);
                                    int first = -1;
                                    for (int gg = 0; gg < NGROUP && first < 0; gg++) {
                                        int clean = 1;
                                        for (int k = 0; k < NK && clean; k++) if (kgroup[k] == gg && m1 * s1[k] + m2 * sig2[k] != 0) clean = 0;
                                        if (clean) first = gg;
                                    }
                                    if (first != grp) continue;
                                }
                                st_pair_cand++;
                                int ids[2] = { (int)i, j2 }, mult[2] = { m1, m2 };
                                try_shape(ids, mult, 2);
                            }
                        }
                    }
                }
            }
        }
        if ((i + 1) % 1000000 == 0) fprintf(logf, "pairs: %zu clusters, lookups %lld candidates %lld survive %lld, %ld s\n", i + 1, st_lookups, st_pair_cand, st_survive, (long)(time(NULL) - t0));
    }
    fclose(out);
    fprintf(logf, "done: clusters %zu single_cand %lld pair_cand %lld lookups %lld survive %lld imbalance_fail %lld decl_fail %lld rows_fail %lld hubhub_fail %lld parity_fail %lld, %ld s\n",
            ncl, st_single_cand, st_pair_cand, st_lookups, st_survive, st_fail_imb, st_fail_decl, st_fail_rows, st_fail_hh, st_fail_par, (long)(time(NULL) - t0));
    return 0;
}
