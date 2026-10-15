/* combine_groups.c: shapes made of any number of clusters (the group combine of
 * combine108.py in C, with the cost ordering that makes it finite for large cluster
 * lists).
 *
 * Ordering. Rank the clusters by (cost desc, deficit desc, id). In a shape of total cost T in
 * {10, 12} whose clusters are c_1 (max cost) and c_2, ..., c_k, every c_i (i >= 2) has cost
 * <= min(cost(c_1), T - cost(c_1)) <= 6: all non-root clusters are "cheap" (cost <= 6, option
 * --cheap8, default 48 eighths). Hence the shape is the root group (the root plus cheap clusters
 * of rank >= root) together with groups of cheap clusters of rank >= root.
 *
 * Phase A (groups): from every cluster as seed, forced growth on the slot imbalance: while the
 * effective imbalance sigma(U) + tokens is nonzero, pick a nonzero coordinate (shortest raw fixer
 * list) and branch over the cheap clusters of rank >= seed with the opposite sign there that
 * still fit (cost, deficit <= D + 2, G <= N), or tokenise one unit (at most 2 e_max tokens: a
 * token is a unit of imbalance that survives in the member, which has at most 2 e inexact slots;
 * completeness as in combine108.py: in a member either a remaining cluster reduces the
 * coordinate or the unit of imbalance survives). A state with zero effective imbalance is a group.
 * States are memoised per seed.
 *
 * Phase B (shapes): multisets of groups with deficit D, cost accepted, (N - G) even, slots
 * (32 - 2u) per block, tokens <= 2 e, by coarse signature (deficit, cost, G, slots, tokens);
 * every multiset of signatures passing the totals is expanded to group multisets and checked:
 * exact imbalance ||sum sigma||_1 <= 2 e, declared slots <= 2 e, rows >= 0, hub-to-hub, (2,2)
 * parity at e = 0. Survivors (any number of clusters; singles and pairs are a cross-check of
 * combine108.c) are written as "u G cost8 nparts | wid:cnt ... | parts m*id ..." for
 * shapes_from_c.py. An expansion larger than --maxexp combinations is skipped and counted
 * (reported as incompleteness).
 *
 * Usage: combine_groups WINDOWS.txt N ACCEPT8 OUT.txt CLUSTERS.bin [...] [--log=FILE]
 *        [--cheap8=48] [--maxexp=20000000] [--groups=FILE] */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "clusters_io.h"

#define MAXP 24          /* clusters per group */
#define MAXTOK 8

static int N, Dtot, nacc, accept8[8], cmax8, emax, cheap8 = 48;
static long long maxexp = 20000000; static size_t progress = 20000;
static FILE *out;
static int *rank_of, *order;

/* fixers: cheap clusters by (coordinate, sign), sorted by cost8 asc */
static int *fix[MAXK][2], nfix[MAXK][2];

typedef struct {
    uint8_t n; int32_t id[MAXP]; uint8_t mult[MAXP];
    uint8_t ntok; uint8_t tokk[MAXTOK]; int8_t toks[MAXTOK];
} GState;

typedef struct {
    int32_t goff;         /* offset into garena: n pairs (cluster id, mult) */
    int16_t cost8, deficit; uint8_t n, ntok, G, slots[3];
    int32_t seed;
} Group;
static Group *groups; static size_t ngroups, capgroups;
static int32_t *garena; static size_t ngar, capgar;

/* hash set (128-bit) */
typedef struct { uint64_t *k1, *k2; size_t cap, n; } HSet;
static void hs_init(HSet *h, size_t cap) { h->cap = cap; h->n = 0; h->k1 = calloc(cap, 8); h->k2 = calloc(cap, 8); if (!h->k1 || !h->k2) { fprintf(stderr, "oom\n"); exit(2); } }
static void hs_free(HSet *h) { free(h->k1); free(h->k2); h->k1 = h->k2 = NULL; h->cap = h->n = 0; }
static int hs_insert_raw(HSet *h, uint64_t a, uint64_t b) {
    size_t mask = h->cap - 1, i = (size_t)(a ^ (b * 0x9E3779B97F4A7C15ULL)) & mask;
    for (;;) {
        if (h->k1[i] == 0 && h->k2[i] == 0) { h->k1[i] = a; h->k2[i] = b; h->n++; return 1; }
        if (h->k1[i] == a && h->k2[i] == b) return 0;
        i = (i + 1) & mask;
    }
}
static int hs_insert(HSet *h, uint64_t a, uint64_t b) {
    if (a == 0 && b == 0) a = 1;
    if (h->n * 5 >= h->cap * 3) { HSet g; hs_init(&g, h->cap * 2); for (size_t i = 0; i < h->cap; i++) if (h->k1[i] || h->k2[i]) hs_insert_raw(&g, h->k1[i], h->k2[i]); hs_free(h); *h = g; }
    return hs_insert_raw(h, a, b);
}

static void gstate_hash(const GState *s, uint64_t *a, uint64_t *b) {
    uint64_t h1 = 0x243F6A8885A308D3ULL, h2 = 0x13198A2E03707344ULL;
    for (int i = 0; i < s->n; i++) {
        uint64_t v = ((uint64_t)(uint32_t)s->id[i] << 8) | s->mult[i];
        h1 = mix64(h1 ^ (v * 0x9E3779B97F4A7C15ULL + (uint64_t)i)); h2 = mix64(h2 + v + 0xA54FF53A5F1D36F1ULL) * 0xD6E8FEB86659FD93ULL ^ (h2 >> 29);
    }
    h1 = mix64(h1 ^ 0xFFFFFFFFULL); h2 = mix64(h2 + 0x510E527FADE682D1ULL);
    for (int i = 0; i < s->ntok; i++) {
        uint64_t v = ((uint64_t)s->tokk[i] << 8) | (uint8_t)s->toks[i];
        h1 = mix64(h1 ^ (v * 0x9E3779B97F4A7C15ULL + 77 + (uint64_t)i)); h2 = mix64(h2 + v + 0x9B05688C2B3E6C1FULL) * 0xD6E8FEB86659FD93ULL ^ (h2 >> 29);
    }
    *a = h1; *b = h2;
}

static int cmp_rank(const void *p, const void *q) {
    const Cl *a = &cl[*(const int *)p], *b = &cl[*(const int *)q];
    if (a->cost8 != b->cost8) return b->cost8 - a->cost8;
    if (a->deficit != b->deficit) return b->deficit - a->deficit;
    return *(const int *)p - *(const int *)q;
}

static int cmp_cost_idx(const void *p, const void *q) {
    const Cl *a = &cl[*(const int *)p], *b = &cl[*(const int *)q];
    if (a->cost8 != b->cost8) return a->cost8 - b->cost8;
    return *(const int *)p - *(const int *)q;
}

static void gstate_add(GState *s, int id) {
    int i = 0;
    while (i < s->n && s->id[i] < id) i++;
    if (i < s->n && s->id[i] == id) { s->mult[i]++; return; }
    if (s->n >= MAXP) { fprintf(stderr, "group overflow\n"); exit(2); }
    for (int k = s->n; k > i; k--) { s->id[k] = s->id[k - 1]; s->mult[k] = s->mult[k - 1]; }
    s->id[i] = id; s->mult[i] = 1; s->n++;
}

static void gstate_token(GState *s, int k, int sign) {
    int i = 0;
    while (i < s->ntok && (s->tokk[i] < k || (s->tokk[i] == k && s->toks[i] < sign))) i++;
    if (s->ntok >= MAXTOK) { fprintf(stderr, "token overflow\n"); exit(2); }
    for (int q = s->ntok; q > i; q--) { s->tokk[q] = s->tokk[q - 1]; s->toks[q] = s->toks[q - 1]; }
    s->tokk[i] = k; s->toks[i] = sign; s->ntok++;
}

static void record_group(const GState *s, int cost8, int deficit, int G, const int *slots, int seed) {
    if (ngroups == capgroups) { capgroups = capgroups ? capgroups * 2 : 1 << 16; groups = realloc(groups, capgroups * sizeof(Group)); if (!groups) { fprintf(stderr, "oom\n"); exit(2); } }
    if (ngar + 2 * s->n > capgar) { capgar = capgar ? capgar * 2 : 1 << 18; if (ngar + 2 * s->n > capgar) capgar = ngar + 2 * s->n; garena = realloc(garena, capgar * sizeof(int32_t)); if (!garena) { fprintf(stderr, "oom\n"); exit(2); } }
    Group *g = &groups[ngroups++];
    g->goff = (int32_t)ngar; g->n = s->n; g->ntok = s->ntok; g->cost8 = cost8; g->deficit = deficit; g->G = G; g->seed = seed;
    for (int j = 0; j < 3; j++) g->slots[j] = slots[j];
    for (int i = 0; i < s->n; i++) { garena[ngar++] = s->id[i]; garena[ngar++] = s->mult[i]; }
}

static long long st_states, st_dead, st_tok_states;
static GState *stack; static size_t stack_n, stack_cap;
static void push(const GState *s) { if (stack_n == stack_cap) { stack_cap = stack_cap ? stack_cap * 2 : 1024; stack = realloc(stack, stack_cap * sizeof(GState)); if (!stack) exit(2); } stack[stack_n++] = *s; }

static void grow_groups(int seed) {
    HSet visited; hs_init(&visited, 1 << 10);
    GState s0; memset(&s0, 0, sizeof s0); s0.n = 1; s0.id[0] = seed; s0.mult[0] = 1;
    uint64_t h1, h2; gstate_hash(&s0, &h1, &h2); hs_insert(&visited, h1, h2);
    stack_n = 0; push(&s0);
    int sig[MAXK];
    int rseed = rank_of[seed];
    while (stack_n) {
        GState s = stack[--stack_n];
        st_states++;
        memset(sig, 0, sizeof(int) * NK);
        int cost8 = 0, deficit = 0, G = 0, slots[3] = {0, 0, 0};
        for (int i = 0; i < s.n; i++) {
            const Cl *c = &cl[s.id[i]]; int m = s.mult[i];
            cost8 += m * c->cost8; deficit += m * c->deficit; G += m * c->G;
            for (int j = 0; j < 3; j++) slots[j] += m * c->slots[j];
            cluster_sig_add(c, m, sig);
        }
        for (int i = 0; i < s.ntok; i++) sig[s.tokk[i]] -= s.toks[i];
        /* the coordinate with the fewest addable fixers (a coordinate without any is a dead end
           unless a token remains) */
        int best = -1, bestn = 1 << 30, tokens_left = 2 * emax - s.ntok;
        for (int k = 0; k < NK && bestn > 0; k++) {
            if (!sig[k]) continue;
            int *lst = fix[k][sig[k] > 0 ? 1 : 0], nl = nfix[k][sig[k] > 0 ? 1 : 0], n = 0;
            for (int q = 0; q < nl && n < bestn; q++) {
                int id = lst[q]; const Cl *c = &cl[id];
                if (cost8 + c->cost8 > cmax8) break;
                if (rank_of[id] < rseed) continue;
                int nd = deficit + c->deficit;
                if ((c->deficit >= 0 && nd > Dtot + 2) || nd < -2) continue;
                if (G + c->G > N || slots[0] + c->slots[0] > 32 || slots[1] + c->slots[1] > 32 || slots[2] + c->slots[2] > 32) continue;
                n++;
            }
            if (n < bestn) { bestn = n; best = k; }
        }
        if (best < 0) { record_group(&s, cost8, deficit, G, slots, seed); continue; }
        if (bestn == 0 && tokens_left <= 0) { st_dead++; continue; }
        int sign = sig[best] > 0 ? 1 : -1;
        int *lst = fix[best][sign > 0 ? 1 : 0], nl = nfix[best][sign > 0 ? 1 : 0];
        int any = 0;
        for (int q = 0; q < nl; q++) {
            int id = lst[q]; const Cl *c = &cl[id];
            if (cost8 + c->cost8 > cmax8) break;
            if (rank_of[id] < rseed) continue;
            int nd = deficit + c->deficit;
            if ((c->deficit >= 0 && nd > Dtot + 2) || nd < -2) continue;
            if (G + c->G > N || slots[0] + c->slots[0] > 32 || slots[1] + c->slots[1] > 32 || slots[2] + c->slots[2] > 32) continue;
            GState t = s; gstate_add(&t, id);
            gstate_hash(&t, &h1, &h2);
            if (!hs_insert(&visited, h1, h2)) continue;
            push(&t); any = 1;
        }
        if (s.ntok < 2 * emax) {
            GState t = s; gstate_token(&t, best, sign);
            gstate_hash(&t, &h1, &h2);
            if (hs_insert(&visited, h1, h2)) { push(&t); any = 1; st_tok_states++; }
        }
        if (!any) st_dead++;
    }
    hs_free(&visited);
}

/* ------------------------------------------------------------------ phase B */
typedef struct { int16_t deficit, cost8; uint8_t G, ntok, slots[3]; int32_t first, count; } Sig;
static Sig *sigs; static int nsigs;
static int32_t *gsorted;   /* group indices sorted by signature */
static long long st_coarse, st_fine, st_survive, st_skipped_exp, st_fail_imb, st_fail_decl, st_fail_rows, st_fail_hh, st_fail_par, st_dup;
static HSet outset;

static int cmp_group_sig(const void *p, const void *q) {
    const Group *a = &groups[*(const int32_t *)p], *b = &groups[*(const int32_t *)q];
    if (a->deficit != b->deficit) return b->deficit - a->deficit;
    if (a->cost8 != b->cost8) return a->cost8 - b->cost8;
    if (a->G != b->G) return a->G - b->G;
    for (int j = 0; j < 3; j++) if (a->slots[j] != b->slots[j]) return a->slots[j] - b->slots[j];
    return a->ntok - b->ntok;
}

static void emit(const int32_t *gl, int ng) {
    /* fine checks on a multiset of groups, then output of the window multiset */
    static int clid[MAXP * 8], clm[MAXP * 8]; int ncl_ = 0;
    for (int i = 0; i < ng; i++) {
        const Group *g = &groups[gl[i]];
        for (int q = 0; q < g->n; q++) {
            int id = garena[g->goff + 2 * q], m = garena[g->goff + 2 * q + 1];
            int f = -1; for (int z = 0; z < ncl_; z++) if (clid[z] == id) { f = z; break; }
            if (f < 0) { clid[ncl_] = id; clm[ncl_] = m; ncl_++; } else clm[f] += m;
        }
    }
    st_fine++;
    int cost8 = 0, deficit = 0, G = 0, decl = 0, ex[6] = {0}, hh[3] = {0, 0, 0}, sv22[3] = {0, 0, 0};
    static int sig[MAXK]; memset(sig, 0, sizeof(int) * NK);
    for (int z = 0; z < ncl_; z++) {
        const Cl *c = &cl[clid[z]]; int m = clm[z];
        cost8 += m * c->cost8; deficit += m * c->deficit; G += m * c->G; decl += m * c->decl;
        for (int j = 0; j < 3; j++) { hh[j] += m * c->hh[j]; sv22[j] += m * c->sv22[j]; }
        for (int k = 0; k < 6; k++) ex[k] += m * c->ex[k];
        cluster_sig_add(c, m, sig);
    }
    int u = (N - G) / 2, excess = Dtot - cost8 / 8, l1 = 0;
    for (int k = 0; k < NK; k++) l1 += sig[k] < 0 ? -sig[k] : sig[k];
    if (l1 > 2 * excess) { st_fail_imb++; return; }
    if (decl > 2 * excess) { st_fail_decl++; return; }
    for (int k = 0; k < 6; k++) if (ex[k] < 0) { st_fail_rows++; return; }
    for (int j = 0; j < 3; j++) if (hh[j] + 2 * u < 16) { st_fail_hh++; return; }
    if (excess == 0) for (int j = 0; j < 3; j++) if ((sv22[j] + 2 * u) % 2) { st_fail_par++; return; }
    /* window multiset, deduplicated */
    static int wc[MAXW], touched[MAXW]; int nt = 0;
    for (int z = 0; z < ncl_; z++) {
        const Cl *c = &cl[clid[z]];
        for (int i = 0; i < c->nwin; i++) { int wid = warena[c->woff + 2 * i]; if (!wc[wid]) touched[nt++] = wid; wc[wid] += clm[z] * warena[c->woff + 2 * i + 1]; }
    }
    /* sort touched for a canonical digest */
    for (int i = 1; i < nt; i++) { int v = touched[i], j = i; while (j > 0 && touched[j - 1] > v) { touched[j] = touched[j - 1]; j--; } touched[j] = v; }
    uint64_t h1 = 0x6A09E667F3BCC908ULL, h2 = 0xBB67AE8584CAA73BULL;
    for (int i = 0; i < nt; i++) { uint64_t v = ((uint64_t)touched[i] << 8) | wc[touched[i]]; h1 = mix64(h1 ^ (v + i)); h2 = mix64(h2 + v * 0x9E3779B97F4A7C15ULL); }
    if (!hs_insert(&outset, h1, h2)) { st_dup++; for (int i = 0; i < nt; i++) wc[touched[i]] = 0; return; }
    st_survive++;
    int nparts = 0; for (int z = 0; z < ncl_; z++) nparts += clm[z];
    fprintf(out, "%d %d %d %d |", u, G, cost8, nparts);
    for (int i = 0; i < nt; i++) fprintf(out, " %d:%d", touched[i], wc[touched[i]]);
    fprintf(out, " | parts");
    for (int z = 0; z < ncl_; z++) fprintf(out, " %d*%d", clm[z], clid[z]);
    fprintf(out, "\n");
    for (int i = 0; i < nt; i++) wc[touched[i]] = 0;
}

/* expand a multiset of signatures (chosen[] with multiplicities) into group multisets */
static int32_t exp_gl[64]; static int exp_ng;
static void expand(const int *chosen, const int *mult, int nch, int pos) {
    if (pos == nch) { emit(exp_gl, exp_ng); return; }
    const Sig *S = &sigs[chosen[pos]]; int m = mult[pos];
    /* combinations with replacement of m groups from S->count */
    int idx[16];
    for (int i = 0; i < m; i++) idx[i] = 0;
    for (;;) {
        for (int i = 0; i < m; i++) exp_gl[exp_ng + i] = gsorted[S->first + idx[i]];
        exp_ng += m; expand(chosen, mult, nch, pos + 1); exp_ng -= m;
        int i = m - 1;
        while (i >= 0 && idx[i] == S->count - 1) i--;
        if (i < 0) break;
        idx[i]++;
        for (int q = i + 1; q < m; q++) idx[q] = idx[i];
    }
}

static long long comb_count(long long n, int k) { long long r = 1; for (int i = 1; i <= k; i++) { r = r * (n + i - 1) / i; if (r > (1LL << 50)) return r; } return r; }

static void knapsack(int start, int *chosen, int *mult, int nch, int deficit, int cost8, int G, const int *slots, int tok) {
    st_coarse++;
    if (nch && deficit == Dtot && (N - G) % 2 == 0) {
        int ok = 0; for (int a = 0; a < nacc; a++) if (accept8[a] == cost8) ok = 1;
        if (ok) {
            int u = (N - G) / 2, g = 32 - 2 * u, excess = Dtot - cost8 / 8;
            if (slots[0] == g && slots[1] == g && slots[2] == g && tok <= 2 * excess) {
                long long total = 1;
                for (int i = 0; i < nch; i++) { total *= comb_count(sigs[chosen[i]].count, mult[i]); if (total > maxexp) break; }
                if (total > maxexp) { st_skipped_exp++; fprintf(logf, "skipped expansion of %lld combinations (signatures %d)\n", total, nch); }
                else { exp_ng = 0; expand(chosen, mult, nch, 0); }
            }
        }
    }
    for (int k = start; k < nsigs; k++) {
        const Sig *S = &sigs[k];
        int nc = cost8 + S->cost8; if (nc > cmax8) continue;
        int nd = deficit + S->deficit; if ((S->deficit >= 0 && nd > Dtot + 2) || nd < -2) continue;
        int ns[3] = { slots[0] + S->slots[0], slots[1] + S->slots[1], slots[2] + S->slots[2] };
        if (ns[0] > 32 || ns[1] > 32 || ns[2] > 32 || G + S->G > N) continue;
        if (tok + S->ntok > 2 * (Dtot - (nc + 7) / 8)) continue;
        if (nch && chosen[nch - 1] == k) { mult[nch - 1]++; knapsack(k, chosen, mult, nch, nd, nc, G + S->G, ns, tok + S->ntok); mult[nch - 1]--; }
        else { chosen[nch] = k; mult[nch] = 1; knapsack(k, chosen, mult, nch + 1, nd, nc, G + S->G, ns, tok + S->ntok); }
    }
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s WINDOWS.txt N ACCEPT8 OUT.txt CLUSTERS.bin [...] [--log=FILE] [--cheap8=48] [--maxexp=M] [--groups=FILE]\n", argv[0]); return 1; }
    logf = stderr; const char *gpath = NULL;
    for (int i = 5; i < argc; i++) {
        if (!strncmp(argv[i], "--log=", 6)) { logf = fopen(argv[i] + 6, "a"); if (!logf) { perror("log"); return 1; } }
        else if (!strncmp(argv[i], "--cheap8=", 9)) cheap8 = atoi(argv[i] + 9);
        else if (!strncmp(argv[i], "--maxexp=", 9)) maxexp = atoll(argv[i] + 9);
        else if (!strncmp(argv[i], "--groups=", 9)) gpath = argv[i] + 9;
        else if (!strncmp(argv[i], "--progress=", 11)) progress = (size_t)atoll(argv[i] + 11);
    }
    N = atoi(argv[2]); Dtot = 3 * N - 96;
    { char *s = strdup(argv[3]), *p = strtok(s, ","); nacc = 0; cmax8 = 0; while (p) { accept8[nacc++] = atoi(p); if (accept8[nacc - 1] > cmax8) cmax8 = accept8[nacc - 1]; p = strtok(NULL, ","); } }
    emax = 0; for (int a = 0; a < nacc; a++) if (Dtot - accept8[a] / 8 > emax) emax = Dtot - accept8[a] / 8;
    time_t t0 = time(NULL);
    load_windows(argv[1]);
    for (int i = 5; i < argc; i++) if (strncmp(argv[i], "--", 2)) load_clusters(argv[i], cmax8, 2 * emax);
    fprintf(logf, "%zu clusters loaded, emax %d, cheap8 %d, %ld s\n", ncl, emax, cheap8, (long)(time(NULL) - t0));
    /* ranks */
    order = malloc(ncl * sizeof(int)); rank_of = malloc(ncl * sizeof(int));
    for (size_t i = 0; i < ncl; i++) order[i] = (int)i;
    qsort(order, ncl, sizeof(int), cmp_rank);
    for (size_t r = 0; r < ncl; r++) rank_of[order[r]] = (int)r;
    /* fixers over the cheap set */
    int *cheap = malloc(ncl * sizeof(int)); int ncheap = 0;
    for (size_t i = 0; i < ncl; i++) if (cl[i].cost8 <= cheap8) cheap[ncheap++] = (int)i;
    qsort(cheap, ncheap, sizeof(int), cmp_cost_idx);
    for (int k = 0; k < NK; k++) for (int s = 0; s < 2; s++) nfix[k][s] = 0;
    int sig[MAXK];
    for (int q = 0; q < ncheap; q++) { cluster_sig(&cl[cheap[q]], sig); for (int k = 0; k < NK; k++) if (sig[k]) nfix[k][sig[k] > 0 ? 0 : 1]++; }
    for (int k = 0; k < NK; k++) for (int s = 0; s < 2; s++) { fix[k][s] = malloc((nfix[k][s] + 1) * sizeof(int)); nfix[k][s] = 0; }
    for (int q = 0; q < ncheap; q++) { cluster_sig(&cl[cheap[q]], sig); for (int k = 0; k < NK; k++) if (sig[k]) { int s = sig[k] > 0 ? 0 : 1; fix[k][s][nfix[k][s]++] = cheap[q]; } }
    fprintf(logf, "%d cheap clusters (cost8 <= %d); fixer lists built, %ld s\n", ncheap, cheap8, (long)(time(NULL) - t0));
    /* phase A */
    for (size_t r = 0; r < ncl; r++) {
        long long s0 = st_states; size_t g0 = ngroups;
        grow_groups(order[r]);
        if (st_states - s0 > 200000) fprintf(logf, "  heavy seed %d (rank %zu, cost8 %d, deficit %d, L1 %d): %lld states, %zu groups\n", order[r], r, cl[order[r]].cost8, cl[order[r]].deficit, cl[order[r]].nsig, st_states - s0, ngroups - g0);
        if ((r + 1) % progress == 0) fprintf(logf, "groups: %zu seeds, %zu groups, states %lld, %ld s\n", r + 1, ngroups, st_states, (long)(time(NULL) - t0));
    }
    fprintf(logf, "phase A done: %zu groups from %zu seeds, states %lld, token states %lld, dead %lld, %ld s\n", ngroups, ncl, st_states, st_tok_states, st_dead, (long)(time(NULL) - t0));
    if (gpath) {
        FILE *gf = fopen(gpath, "w");
        for (size_t g = 0; g < ngroups; g++) { fprintf(gf, "%d %d %d %d %d %d %d |", groups[g].deficit, groups[g].cost8, groups[g].G, groups[g].slots[0], groups[g].slots[1], groups[g].slots[2], groups[g].ntok); for (int q = 0; q < groups[g].n; q++) fprintf(gf, " %d*%d", garena[groups[g].goff + 2 * q + 1], garena[groups[g].goff + 2 * q]); fprintf(gf, "\n"); }
        fclose(gf);
    }
    /* signatures */
    gsorted = malloc(ngroups * sizeof(int32_t));
    for (size_t g = 0; g < ngroups; g++) gsorted[g] = (int32_t)g;
    qsort(gsorted, ngroups, sizeof(int32_t), cmp_group_sig);
    sigs = malloc((ngroups + 1) * sizeof(Sig)); nsigs = 0;
    for (size_t g = 0; g < ngroups; g++) {
        const Group *G_ = &groups[gsorted[g]];
        if (nsigs && sigs[nsigs - 1].deficit == G_->deficit && sigs[nsigs - 1].cost8 == G_->cost8 && sigs[nsigs - 1].G == G_->G && sigs[nsigs - 1].ntok == G_->ntok && !memcmp(sigs[nsigs - 1].slots, G_->slots, 3)) { sigs[nsigs - 1].count++; continue; }
        Sig *S = &sigs[nsigs++]; S->deficit = G_->deficit; S->cost8 = G_->cost8; S->G = G_->G; S->ntok = G_->ntok; memcpy(S->slots, G_->slots, 3); S->first = (int32_t)g; S->count = 1;
    }
    fprintf(logf, "%d coarse signatures, %ld s\n", nsigs, (long)(time(NULL) - t0));
    out = fopen(argv[4], "w"); if (!out) { perror(argv[4]); return 1; }
    hs_init(&outset, 1 << 16);
    int chosen[64], mult[64], slots0[3] = {0, 0, 0};
    knapsack(0, chosen, mult, 0, 0, 0, 0, slots0, 0);
    fclose(out);
    fprintf(logf, "done: coarse nodes %lld fine multisets %lld survive %lld dup %lld skipped_expansions %lld imbalance_fail %lld decl_fail %lld rows_fail %lld hubhub_fail %lld parity_fail %lld, %ld s\n",
            st_coarse, st_fine, st_survive, st_dup, st_skipped_exp, st_fail_imb, st_fail_decl, st_fail_rows, st_fail_hh, st_fail_par, (long)(time(NULL) - t0));
    return st_skipped_exp ? 3 : 0;
}
