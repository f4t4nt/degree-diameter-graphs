/* combine_classes.c: the group combine of combine_groups.c run on *classes* of clusters.
 *
 * Two clusters are in the same class when they agree on everything the shape conditions see:
 * cost, deficit, G, slots, the six rows, hub-to-hub, (2,2)-slot counts, declared slots and the
 * slot imbalance vector sigma. Every shape condition of combine_groups.c (and of pairs108.py /
 * combine108.py) is a function of these class data, so a multiset of classes either passes all
 * conditions for every choice of instances or for none. The combine therefore runs on classes
 * (phase A: forced growth on sigma with cheap fixer classes, rank order (cost desc, deficit
 * desc, class id); phase B: knapsack over group signatures and fine checks), and the number of
 * instance multisets of a surviving class multiset is prod_c C(n_c + m_c - 1, m_c) (n_c
 * instances, multiplicity m_c). Instance multisets (window multisets) are written out while the
 * running total stays below --maxout; beyond that they are only counted (reported as such).
 *
 * Usage: combine_classes WINDOWS.txt N ACCEPT8 OUT.txt CLUSTERS.bin [...] [--log=FILE]
 *        [--cheap8=48] [--maxexp=M] [--maxout=M] [--progress=K] */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "clusters_io.h"

#define MAXP 24
#define MAXTOK 8

static int N, Dtot, nacc, accept8[8], cmax8, emax, cheap8 = 48;
static long long maxexp = 20000000, maxout = 2000000; static size_t progress = 100000;
static FILE *out;

/* classes */
typedef struct { int32_t rep; int32_t first, count; } Class;   /* rep = representative cluster; instances in inst[first .. first+count) */
static Class *cls; static int ncls;
static int32_t *inst;
static int *rank_of, *order;
static int *fix[MAXK][2], nfix[MAXK][2];     /* cheap classes by (coordinate, sign), sorted by cost8 asc */

typedef struct {
    uint8_t n; int32_t id[MAXP]; uint8_t mult[MAXP];
    uint8_t ntok; uint8_t tokk[MAXTOK]; int8_t toks[MAXTOK];
} GState;

typedef struct { int32_t goff; int16_t cost8, deficit; uint8_t n, ntok, G, slots[3]; } Group;
static Group *groups; static size_t ngroups, capgroups;
static int32_t *garena; static size_t ngar, capgar;

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

/* ------------------------------------------------------------------ classes */
static uint64_t *ckey1, *ckey2;
static void class_key(const Cl *c, uint64_t *a, uint64_t *b) {
    uint64_t h1 = 0x3C6EF372FE94F82BULL, h2 = 0xA54FF53A5F1D36F1ULL;
#define FEED(v) do { uint64_t _v = (uint64_t)(uint32_t)(int32_t)(v); h1 = mix64(h1 ^ (_v * 0x9E3779B97F4A7C15ULL + 1)); h2 = mix64(h2 + _v + 0x5BE0CD19137E2179ULL) * 0xD6E8FEB86659FD93ULL ^ (h2 >> 31); } while (0)
    FEED(c->cost8); FEED(c->deficit); FEED(c->G); FEED(c->decl);
    for (int j = 0; j < 3; j++) { FEED(c->slots[j]); FEED(c->hh[j]); FEED(c->sv22[j]); }
    for (int k = 0; k < 6; k++) FEED(c->ex[k]);
    FEED(c->nsig);
    for (int s = 0; s < c->nsig; s++) { FEED((uint8_t)sarena[c->soff + 2 * s]); FEED(sarena[c->soff + 2 * s + 1]); }
#undef FEED
    *a = h1; *b = h2;
}

static int cmp_key(const void *p, const void *q) {
    int32_t a = *(const int32_t *)p, b = *(const int32_t *)q;
    if (ckey1[a] != ckey1[b]) return ckey1[a] < ckey1[b] ? -1 : 1;
    if (ckey2[a] != ckey2[b]) return ckey2[a] < ckey2[b] ? -1 : 1;
    return a - b;
}

static int same_class(const Cl *a, const Cl *b) {
    if (a->cost8 != b->cost8 || a->deficit != b->deficit || a->G != b->G || a->decl != b->decl || a->nsig != b->nsig) return 0;
    if (memcmp(a->slots, b->slots, 3) || memcmp(a->hh, b->hh, sizeof a->hh) || memcmp(a->sv22, b->sv22, 3) || memcmp(a->ex, b->ex, sizeof a->ex)) return 0;
    return !memcmp(sarena + a->soff, sarena + b->soff, 2 * a->nsig);
}

static void build_classes(void) {
    ckey1 = malloc(ncl * 8); ckey2 = malloc(ncl * 8); inst = malloc(ncl * sizeof(int32_t));
    for (size_t i = 0; i < ncl; i++) { class_key(&cl[i], &ckey1[i], &ckey2[i]); inst[i] = (int32_t)i; }
    qsort(inst, ncl, sizeof(int32_t), cmp_key);
    cls = malloc(ncl * sizeof(Class)); ncls = 0;
    for (size_t i = 0; i < ncl; i++) {
        if (ncls && ckey1[inst[i]] == ckey1[cls[ncls - 1].rep] && ckey2[inst[i]] == ckey2[cls[ncls - 1].rep]) {
            if (!same_class(&cl[inst[i]], &cl[cls[ncls - 1].rep])) { fprintf(stderr, "class hash collision\n"); exit(2); }
            cls[ncls - 1].count++; continue;
        }
        cls[ncls].rep = inst[i]; cls[ncls].first = (int32_t)i; cls[ncls].count = 1; ncls++;
    }
    free(ckey1); free(ckey2);
}

static int cmp_rank(const void *p, const void *q) {
    const Cl *a = &cl[cls[*(const int *)p].rep], *b = &cl[cls[*(const int *)q].rep];
    if (a->cost8 != b->cost8) return b->cost8 - a->cost8;
    if (a->deficit != b->deficit) return b->deficit - a->deficit;
    return *(const int *)p - *(const int *)q;
}
static int cmp_cost_idx(const void *p, const void *q) {
    const Cl *a = &cl[cls[*(const int *)p].rep], *b = &cl[cls[*(const int *)q].rep];
    if (a->cost8 != b->cost8) return a->cost8 - b->cost8;
    return *(const int *)p - *(const int *)q;
}

/* ------------------------------------------------------------------ phase A */
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
static void record_group(const GState *s, int cost8, int deficit, int G, const int *slots) {
    if (ngroups == capgroups) { capgroups = capgroups ? capgroups * 2 : 1 << 16; groups = realloc(groups, capgroups * sizeof(Group)); if (!groups) { fprintf(stderr, "oom\n"); exit(2); } }
    if (ngar + 2 * s->n > capgar) { capgar = capgar ? capgar * 2 : 1 << 18; if (ngar + 2 * s->n > capgar) capgar = ngar + 2 * s->n; garena = realloc(garena, capgar * sizeof(int32_t)); if (!garena) { fprintf(stderr, "oom\n"); exit(2); } }
    Group *g = &groups[ngroups++];
    g->goff = (int32_t)ngar; g->n = s->n; g->ntok = s->ntok; g->cost8 = cost8; g->deficit = deficit; g->G = G;
    for (int j = 0; j < 3; j++) g->slots[j] = slots[j];
    for (int i = 0; i < s->n; i++) { garena[ngar++] = s->id[i]; garena[ngar++] = s->mult[i]; }
}

static long long st_states, st_dead, st_tok_states;
static GState *stack; static size_t stack_n, stack_cap;
static void push(const GState *s) { if (stack_n == stack_cap) { stack_cap = stack_cap ? stack_cap * 2 : 1024; stack = realloc(stack, stack_cap * sizeof(GState)); if (!stack) exit(2); } stack[stack_n++] = *s; }

static inline const Cl *CR(int class_id) { return &cl[cls[class_id].rep]; }

static int prefix_upto(int k, int s, int budget) {
    /* number of fixers of (k, s) with cost8 <= budget (lists sorted by cost) */
    int lo = 0, hi = nfix[k][s];
    while (lo < hi) { int mid = (lo + hi) / 2; if (CR(fix[k][s][mid])->cost8 <= budget) lo = mid + 1; else hi = mid; }
    return lo;
}

static inline int addable(const Cl *c, int rseed, int cid, int cost8, int deficit, int G, const int *slots) {
    if (cost8 + c->cost8 > cmax8 || rank_of[cid] < rseed) return 0;
    int nd = deficit + c->deficit;
    if ((c->deficit >= 0 && nd > Dtot + 2) || nd < -2) return 0;
    if (G + c->G > N || slots[0] + c->slots[0] > 32 || slots[1] + c->slots[1] > 32 || slots[2] + c->slots[2] > 32) return 0;
    return 1;
}

static void grow_groups(int seed) {
    HSet visited; hs_init(&visited, 1 << 8);
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
            const Cl *c = CR(s.id[i]); int m = s.mult[i];
            cost8 += m * c->cost8; deficit += m * c->deficit; G += m * c->G;
            for (int j = 0; j < 3; j++) slots[j] += m * c->slots[j];
            cluster_sig_add(c, m, sig);
        }
        for (int i = 0; i < s.ntok; i++) sig[s.tokk[i]] -= s.toks[i];
        int budget = cmax8 - cost8, best = -1, bestn = 1 << 30, tokens_left = 2 * emax - s.ntok;
        for (int k = 0; k < NK && bestn > 0; k++) {
            if (!sig[k]) continue;
            int n = prefix_upto(k, sig[k] > 0 ? 1 : 0, budget);
            if (n < bestn) { bestn = n; best = k; }
        }
        if (best < 0) { record_group(&s, cost8, deficit, G, slots); continue; }
        if (bestn == 0 && tokens_left <= 0) { st_dead++; continue; }
        int sgn = sig[best] > 0 ? 1 : 0;
        int *lst = fix[best][sgn];
        int any = 0;
        for (int q = 0; q < bestn; q++) {
            int cid = lst[q];
            if (!addable(CR(cid), rseed, cid, cost8, deficit, G, slots)) continue;
            GState t = s; gstate_add(&t, cid);
            gstate_hash(&t, &h1, &h2);
            if (!hs_insert(&visited, h1, h2)) continue;
            push(&t); any = 1;
        }
        if (tokens_left > 0) {
            GState t = s; gstate_token(&t, best, sig[best] > 0 ? 1 : -1);
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
static int32_t *gsorted;
static long long st_coarse, st_fine, st_survive, st_skipped_exp, st_fail_imb, st_fail_decl, st_fail_rows, st_fail_hh, st_fail_par, st_dup, st_written;
static long double st_instances = 0;      /* instance multisets of surviving class multisets */
static long long st_inst_by_cost[128], st_inst_by_parts[64];
static HSet outset;

static int cmp_group_sig(const void *p, const void *q) {
    const Group *a = &groups[*(const int32_t *)p], *b = &groups[*(const int32_t *)q];
    if (a->deficit != b->deficit) return b->deficit - a->deficit;
    if (a->cost8 != b->cost8) return a->cost8 - b->cost8;
    if (a->G != b->G) return a->G - b->G;
    for (int j = 0; j < 3; j++) if (a->slots[j] != b->slots[j]) return a->slots[j] - b->slots[j];
    return a->ntok - b->ntok;
}

static long double binom_rep(long long n, int m) { long double r = 1; for (int i = 1; i <= m; i++) r = r * (n + i - 1) / i; return r; }

/* write every instance multiset of the class multiset (clid/clm) as a window multiset */
static int wc_[MAXW];
static void write_instances(const int *clid, const int *clm, int nc, int pos, int *chosen_inst, int u, int G, int cost8, int nparts) {
    if (pos == nc) {
        static int touched[MAXW]; int nt = 0;
        for (int z = 0, w = 0; z < nc; z++) for (int m = 0; m < clm[z]; m++, w++) {
            const Cl *c = &cl[chosen_inst[w]];
            for (int i = 0; i < c->nwin; i++) { int wid = warena[c->woff + 2 * i]; if (!wc_[wid]) touched[nt++] = wid; wc_[wid] += warena[c->woff + 2 * i + 1]; }
        }
        for (int i = 1; i < nt; i++) { int v = touched[i], j = i; while (j > 0 && touched[j - 1] > v) { touched[j] = touched[j - 1]; j--; } touched[j] = v; }
        fprintf(out, "%d %d %d %d |", u, G, cost8, nparts);
        for (int i = 0; i < nt; i++) fprintf(out, " %d:%d", touched[i], wc_[touched[i]]);
        fprintf(out, " | inst");
        for (int z = 0, w = 0; z < nc; z++) for (int m = 0; m < clm[z]; m++, w++) fprintf(out, " %d", chosen_inst[w]);
        fprintf(out, "\n");
        for (int i = 0; i < nt; i++) wc_[touched[i]] = 0;
        st_written++;
        return;
    }
    /* combinations with replacement of clm[pos] instances of class clid[pos] */
    const Class *C = &cls[clid[pos]]; int m = clm[pos];
    int base = 0; for (int z = 0; z < pos; z++) base += clm[z];
    int idx[16]; for (int i = 0; i < m; i++) idx[i] = 0;
    for (;;) {
        for (int i = 0; i < m; i++) chosen_inst[base + i] = inst[C->first + idx[i]];
        write_instances(clid, clm, nc, pos + 1, chosen_inst, u, G, cost8, nparts);
        int i = m - 1;
        while (i >= 0 && idx[i] == C->count - 1) i--;
        if (i < 0) break;
        idx[i]++;
        for (int q = i + 1; q < m; q++) idx[q] = idx[i];
    }
}

static void emit(const int32_t *gl, int ng) {
    static int clid[MAXP * 8], clm[MAXP * 8]; int nc = 0;
    for (int i = 0; i < ng; i++) {
        const Group *g = &groups[gl[i]];
        for (int q = 0; q < g->n; q++) {
            int id = garena[g->goff + 2 * q], m = garena[g->goff + 2 * q + 1];
            int f = -1; for (int z = 0; z < nc; z++) if (clid[z] == id) { f = z; break; }
            if (f < 0) { clid[nc] = id; clm[nc] = m; nc++; } else clm[f] += m;
        }
    }
    st_fine++;
    /* sort by class id for the digest */
    for (int i = 1; i < nc; i++) { int v = clid[i], mv = clm[i], j = i; while (j > 0 && clid[j - 1] > v) { clid[j] = clid[j - 1]; clm[j] = clm[j - 1]; j--; } clid[j] = v; clm[j] = mv; }
    uint64_t h1 = 0x6A09E667F3BCC908ULL, h2 = 0xBB67AE8584CAA73BULL;
    for (int i = 0; i < nc; i++) { uint64_t v = ((uint64_t)(uint32_t)clid[i] << 8) | clm[i]; h1 = mix64(h1 ^ (v + i)); h2 = mix64(h2 + v * 0x9E3779B97F4A7C15ULL); }
    if (!hs_insert(&outset, h1, h2)) { st_dup++; return; }
    int cost8 = 0, deficit = 0, G = 0, decl = 0, ex[6] = {0}, hh[3] = {0, 0, 0}, sv22[3] = {0, 0, 0};
    static int sig[MAXK]; memset(sig, 0, sizeof(int) * NK);
    for (int z = 0; z < nc; z++) {
        const Cl *c = CR(clid[z]); int m = clm[z];
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
    st_survive++;
    long double ninst = 1; int nparts = 0;
    for (int z = 0; z < nc; z++) { ninst *= binom_rep(cls[clid[z]].count, clm[z]); nparts += clm[z]; }
    st_instances += ninst;
    st_inst_by_cost[cost8] += (long long)ninst; st_inst_by_parts[nparts < 63 ? nparts : 63] += (long long)ninst;
    if (st_instances <= (long double)maxout) {
        static int chosen_inst[256];
        write_instances(clid, clm, nc, 0, chosen_inst, u, G, cost8, nparts);
    } else {
        fprintf(out, "# class multiset not expanded (%.0Lf instance multisets): cost8 %d G %d parts %d classes", ninst, cost8, G, nparts);
        for (int z = 0; z < nc; z++) fprintf(out, " %d*%d(n=%d)", clm[z], clid[z], cls[clid[z]].count);
        fprintf(out, "\n");
    }
}

static int32_t exp_gl[64]; static int exp_ng;
static void expand(const int *chosen, const int *mult, int nch, int pos) {
    if (pos == nch) { emit(exp_gl, exp_ng); return; }
    const Sig *S = &sigs[chosen[pos]]; int m = mult[pos];
    int idx[16]; for (int i = 0; i < m; i++) idx[i] = 0;
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

static void knapsack(int start, int *chosen, int *mult, int nch, int deficit, int cost8, int G, const int *slots, int tok) {
    st_coarse++;
    if (nch && deficit == Dtot && (N - G) % 2 == 0) {
        int ok = 0; for (int a = 0; a < nacc; a++) if (accept8[a] == cost8) ok = 1;
        if (ok) {
            int u = (N - G) / 2, g = 32 - 2 * u, excess = Dtot - cost8 / 8;
            if (slots[0] == g && slots[1] == g && slots[2] == g && tok <= 2 * excess) {
                long double total = 1;
                for (int i = 0; i < nch; i++) total *= binom_rep(sigs[chosen[i]].count, mult[i]);
                if (total > (long double)maxexp) { st_skipped_exp++; fprintf(logf, "skipped expansion of %.0Lf group combinations (%d signatures)\n", total, nch); }
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
    if (argc < 6) { fprintf(stderr, "usage: %s WINDOWS.txt N ACCEPT8 OUT.txt CLUSTERS.bin [...] [--log=FILE] [--cheap8=48] [--maxexp=M] [--maxout=M] [--progress=K]\n", argv[0]); return 1; }
    logf = stderr;
    for (int i = 5; i < argc; i++) {
        if (!strncmp(argv[i], "--log=", 6)) { logf = fopen(argv[i] + 6, "a"); if (!logf) { perror("log"); return 1; } }
        else if (!strncmp(argv[i], "--cheap8=", 9)) cheap8 = atoi(argv[i] + 9);
        else if (!strncmp(argv[i], "--maxexp=", 9)) maxexp = atoll(argv[i] + 9);
        else if (!strncmp(argv[i], "--maxout=", 9)) maxout = atoll(argv[i] + 9);
        else if (!strncmp(argv[i], "--progress=", 11)) progress = (size_t)atoll(argv[i] + 11);
    }
    N = atoi(argv[2]); Dtot = 3 * N - 96;
    { char *s = strdup(argv[3]), *p = strtok(s, ","); nacc = 0; cmax8 = 0; while (p) { accept8[nacc++] = atoi(p); if (accept8[nacc - 1] > cmax8) cmax8 = accept8[nacc - 1]; p = strtok(NULL, ","); } }
    emax = 0; for (int a = 0; a < nacc; a++) if (Dtot - accept8[a] / 8 > emax) emax = Dtot - accept8[a] / 8;
    time_t t0 = time(NULL);
    load_windows(argv[1]);
    for (int i = 5; i < argc; i++) if (strncmp(argv[i], "--", 2)) load_clusters(argv[i], cmax8, 2 * emax);
    fprintf(logf, "%zu clusters loaded, emax %d, cheap8 %d, %ld s\n", ncl, emax, cheap8, (long)(time(NULL) - t0));
    build_classes();
    int ncheap_cls = 0; long long ncheap_inst = 0;
    for (int c = 0; c < ncls; c++) if (CR(c)->cost8 <= cheap8) { ncheap_cls++; ncheap_inst += cls[c].count; }
    fprintf(logf, "%d classes (%d cheap classes with %lld instances), %ld s\n", ncls, ncheap_cls, ncheap_inst, (long)(time(NULL) - t0));
    order = malloc(ncls * sizeof(int)); rank_of = malloc(ncls * sizeof(int));
    for (int i = 0; i < ncls; i++) order[i] = i;
    qsort(order, ncls, sizeof(int), cmp_rank);
    for (int r = 0; r < ncls; r++) rank_of[order[r]] = r;
    int *cheap = malloc(ncls * sizeof(int)); int ncheap = 0;
    for (int c = 0; c < ncls; c++) if (CR(c)->cost8 <= cheap8) cheap[ncheap++] = c;
    qsort(cheap, ncheap, sizeof(int), cmp_cost_idx);
    int sig[MAXK];
    for (int k = 0; k < NK; k++) for (int s = 0; s < 2; s++) nfix[k][s] = 0;
    for (int q = 0; q < ncheap; q++) { cluster_sig(CR(cheap[q]), sig); for (int k = 0; k < NK; k++) if (sig[k]) nfix[k][sig[k] > 0 ? 0 : 1]++; }
    for (int k = 0; k < NK; k++) for (int s = 0; s < 2; s++) { fix[k][s] = malloc((nfix[k][s] + 1) * sizeof(int)); nfix[k][s] = 0; }
    for (int q = 0; q < ncheap; q++) { cluster_sig(CR(cheap[q]), sig); for (int k = 0; k < NK; k++) if (sig[k]) { int s = sig[k] > 0 ? 0 : 1; fix[k][s][nfix[k][s]++] = cheap[q]; } }
    fprintf(logf, "fixer lists built, %ld s\n", (long)(time(NULL) - t0));
    for (int r = 0; r < ncls; r++) {
        long long s0 = st_states; size_t g0 = ngroups;
        grow_groups(order[r]);
        if (st_states - s0 > 2000000) fprintf(logf, "  heavy seed class %d (rank %d, cost8 %d, deficit %d, nsig %d, %d instances): %lld states, %zu groups\n", order[r], r, CR(order[r])->cost8, CR(order[r])->deficit, CR(order[r])->nsig, cls[order[r]].count, st_states - s0, ngroups - g0);
        if ((size_t)(r + 1) % progress == 0) { fprintf(logf, "groups: %d seed classes, %zu groups, states %lld, %ld s\n", r + 1, ngroups, st_states, (long)(time(NULL) - t0)); fflush(logf); }
    }
    fprintf(logf, "phase A done: %zu groups from %d seed classes, states %lld, token states %lld, dead %lld, %ld s\n", ngroups, ncls, st_states, st_tok_states, st_dead, (long)(time(NULL) - t0));
    gsorted = malloc((ngroups + 1) * sizeof(int32_t));
    for (size_t g = 0; g < ngroups; g++) gsorted[g] = (int32_t)g;
    qsort(gsorted, ngroups, sizeof(int32_t), cmp_group_sig);
    sigs = malloc((ngroups + 1) * sizeof(Sig)); nsigs = 0;
    for (size_t g = 0; g < ngroups; g++) {
        const Group *G_ = &groups[gsorted[g]];
        if (nsigs && sigs[nsigs - 1].deficit == G_->deficit && sigs[nsigs - 1].cost8 == G_->cost8 && sigs[nsigs - 1].G == G_->G && sigs[nsigs - 1].ntok == G_->ntok && !memcmp(sigs[nsigs - 1].slots, G_->slots, 3)) { sigs[nsigs - 1].count++; continue; }
        Sig *S = &sigs[nsigs++]; S->deficit = G_->deficit; S->cost8 = G_->cost8; S->G = G_->G; S->ntok = G_->ntok; memcpy(S->slots, G_->slots, 3); S->first = (int32_t)g; S->count = 1;
    }
    fprintf(logf, "%d coarse signatures, %ld s\n", nsigs, (long)(time(NULL) - t0)); fflush(logf);
    out = fopen(argv[4], "w"); if (!out) { perror(argv[4]); return 1; }
    hs_init(&outset, 1 << 16);
    int chosen[64], mult[64], slots0[3] = {0, 0, 0};
    knapsack(0, chosen, mult, 0, 0, 0, 0, slots0, 0);
    fclose(out);
    fprintf(logf, "done: coarse nodes %lld fine class multisets %lld distinct %lld survive %lld skipped_expansions %lld imbalance_fail %lld decl_fail %lld rows_fail %lld hubhub_fail %lld parity_fail %lld; instance multisets %.0Lf, written %lld, %ld s\n",
            st_coarse, st_fine, st_fine - st_dup, st_survive, st_skipped_exp, st_fail_imb, st_fail_decl, st_fail_rows, st_fail_hh, st_fail_par, st_instances, st_written, (long)(time(NULL) - t0));
    fprintf(logf, "instance multisets by cost8:"); for (int c = 0; c < 128; c++) if (st_inst_by_cost[c]) fprintf(logf, " %d:%lld", c, st_inst_by_cost[c]); fprintf(logf, "\n");
    fprintf(logf, "instance multisets by parts:"); for (int c = 0; c < 64; c++) if (st_inst_by_parts[c]) fprintf(logf, " %d:%lld", c, st_inst_by_parts[c]); fprintf(logf, "\n");
    return st_skipped_exp ? 3 : 0;
}
