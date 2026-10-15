/* clusters_io.h: window table (export_windows.py format) and cluster file (clusters108.c binary
 * format) loading shared by combine108.c and combine_groups.c, with the slot-imbalance (sigma)
 * coordinates: for block j and foreign vector fv != (2,2), key (j, a) with a = min(fv, 4 - fv) in
 * lexicographic order (Python tuple order, negative entries allowed), count(fv) - count(4 - fv). */
#ifndef CLUSTERS_IO_H
#define CLUSTERS_IO_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NT 20
#define MAXW 4096
#define MAXK 256

typedef struct {
    int id, t, nx, quad, x[3], cost8, deficit, nh, hm[4], hk[4], ns, sj[4], sf0[4], sf1[4], ex[6], hh[3];
} Win;
static Win win[MAXW];
static int NW;
static int coord_of[3][32][32];
static int NK = 0;
static int ka[MAXK][3];            /* coordinate -> (j, a0, a1) */
static int wsig[MAXW][8], wsig_n[MAXW], wsv22[MAXW][3];

typedef struct {
    int32_t woff;        /* offset into warena: nwin pairs (id, cnt) */
    int16_t cost8, deficit, ex[6], hh[3];
    uint8_t nwin, G, quads, decl, slots[3], sv22[3], nsig;
    int32_t soff;        /* offset into sarena: nsig pairs (coord, val) */
} Cl;
static Cl *cl; static size_t ncl, capcl;
static uint16_t *warena; static size_t nwar, capwar;
static int8_t *sarena; static size_t nsar, capsar;
static FILE *logf;

static inline uint64_t mix64(uint64_t z) {
    z ^= z >> 30; z *= 0xbf58476d1ce4e5b9ULL; z ^= z >> 27; z *= 0x94d049bb133111ebULL; z ^= z >> 31;
    return z;
}

static void load_windows(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }
    int nt, dummy[6];
    if (fscanf(f, "NT %d", &nt) != 1 || nt != NT) { fprintf(stderr, "bad NT\n"); exit(1); }
    for (int t = 0; t < NT; t++) if (fscanf(f, "%d %d %d %d %d %d", &dummy[0], &dummy[1], &dummy[2], &dummy[3], &dummy[4], &dummy[5]) != 6) exit(1);
    for (int t = 0; t < NT * NT; t++) if (fscanf(f, "%d", &dummy[0]) != 1) exit(1);
    if (fscanf(f, " NW %d", &NW) != 1 || NW > MAXW) { fprintf(stderr, "bad NW\n"); exit(1); }
    for (int k = 0; k < NW; k++) {
        Win *w = &win[k];
        if (fscanf(f, "%d %d %d %d %d %d %d %d %d %d", &w->id, &w->t, &w->nx, &w->quad, &w->x[0], &w->x[1], &w->x[2], &w->cost8, &w->deficit, &w->nh) != 10 || w->id != k) { fprintf(stderr, "bad window %d\n", k); exit(1); }
        for (int i = 0; i < w->nh; i++) if (fscanf(f, "%d %d", &w->hm[i], &w->hk[i]) != 2) exit(1);
        if (fscanf(f, "%d", &w->ns) != 1) exit(1);
        for (int i = 0; i < w->ns; i++) if (fscanf(f, "%d %d %d", &w->sj[i], &w->sf0[i], &w->sf1[i]) != 3) exit(1);
        for (int i = 0; i < 6; i++) if (fscanf(f, "%d", &w->ex[i]) != 1) exit(1);
        for (int i = 0; i < 3; i++) if (fscanf(f, "%d", &w->hh[i]) != 1) exit(1);
    }
    fclose(f);
    for (int j = 0; j < 3; j++) for (int a = 0; a < 32; a++) for (int b = 0; b < 32; b++) coord_of[j][a][b] = -1;
    for (int k = 0; k < NW; k++) {
        Win *w = &win[k];
        int tmp[MAXK]; memset(tmp, 0, sizeof tmp);
        wsv22[k][0] = wsv22[k][1] = wsv22[k][2] = 0;
        for (int i = 0; i < w->ns; i++) {
            int j = w->sj[i], f0 = w->sf0[i], f1 = w->sf1[i];
            if (f0 == 2 && f1 == 2) { wsv22[k][j]++; continue; }
            int b0 = 4 - f0, b1 = 4 - f1, sign, a0, a1;
            if (f0 < b0 || (f0 == b0 && f1 < b1)) { a0 = f0; a1 = f1; sign = 1; } else { a0 = b0; a1 = b1; sign = -1; }
            if (a0 + 8 < 0 || a1 + 8 < 0 || a0 + 8 >= 32 || a1 + 8 >= 32) { fprintf(stderr, "vector out of range\n"); exit(1); }
            int *c = &coord_of[j][a0 + 8][a1 + 8];
            if (*c < 0) { if (NK >= MAXK) { fprintf(stderr, "too many coordinates\n"); exit(1); } ka[NK][0] = j; ka[NK][1] = a0; ka[NK][2] = a1; *c = NK++; }
            tmp[*c] += sign;
        }
        wsig_n[k] = 0;
        for (int c = 0; c < NK; c++) if (tmp[c]) { if (wsig_n[k] >= 8) { fprintf(stderr, "window sigma overflow\n"); exit(1); } wsig[k][wsig_n[k]++] = c * 256 + (tmp[c] + 128); }
    }
    fprintf(logf, "%d windows, %d sigma coordinates\n", NW, NK);
}

static void add_cluster(const uint16_t *ids, const uint8_t *cnts, int nwin, int decl) {
    if (ncl == capcl) { capcl = capcl ? capcl * 2 : 1 << 16; cl = realloc(cl, capcl * sizeof(Cl)); if (!cl) { fprintf(stderr, "oom\n"); exit(2); } }
    if (nwar + 2 * nwin > capwar) { capwar = capwar ? capwar * 2 : 1 << 20; if (nwar + 2 * nwin > capwar) capwar = nwar + 2 * nwin; warena = realloc(warena, capwar * sizeof(uint16_t)); if (!warena) { fprintf(stderr, "oom\n"); exit(2); } }
    Cl *c = &cl[ncl];
    memset(c, 0, sizeof *c);
    c->woff = (int32_t)nwar; c->nwin = nwin; c->decl = decl;
    int sig[MAXK]; memset(sig, 0, sizeof sig);
    int cost8 = 0, deficit = 0, G = 0, quads = 0, slots[3] = {0, 0, 0}, ex[6] = {0}, hh[3] = {0, 0, 0}, sv22[3] = {0, 0, 0};
    for (int i = 0; i < nwin; i++) {
        warena[nwar++] = ids[i]; warena[nwar++] = cnts[i];
        const Win *w = &win[ids[i]]; int m = cnts[i];
        cost8 += m * w->cost8; deficit += m * w->deficit; G += m; if (w->quad) quads += m;
        for (int j = 0; j < 3; j++) { slots[j] += m * w->x[j]; hh[j] += m * w->hh[j]; sv22[j] += m * wsv22[ids[i]][j]; }
        for (int k = 0; k < 6; k++) ex[k] += m * w->ex[k];
        for (int s = 0; s < wsig_n[ids[i]]; s++) sig[wsig[ids[i]][s] / 256] += m * (wsig[ids[i]][s] % 256 - 128);
    }
    c->cost8 = cost8; c->deficit = deficit; c->G = G; c->quads = quads;
    for (int j = 0; j < 3; j++) { c->slots[j] = slots[j]; c->hh[j] = hh[j]; c->sv22[j] = sv22[j]; }
    for (int k = 0; k < 6; k++) c->ex[k] = ex[k];
    int ns = 0; for (int k = 0; k < NK; k++) if (sig[k]) ns++;
    if (nsar + 2 * ns > capsar) { capsar = capsar ? capsar * 2 : 1 << 20; if (nsar + 2 * ns > capsar) capsar = nsar + 2 * ns; sarena = realloc(sarena, capsar); if (!sarena) { fprintf(stderr, "oom\n"); exit(2); } }
    c->soff = (int32_t)nsar; c->nsig = ns;
    for (int k = 0; k < NK; k++) if (sig[k]) { if (sig[k] > 127 || sig[k] < -127) { fprintf(stderr, "sigma overflow\n"); exit(2); } sarena[nsar++] = (int8_t)k; sarena[nsar++] = (int8_t)sig[k]; }
    ncl++;
}

/* maxdecl: clusters with more declared slots are dropped (they cannot be in a shape of the
   accepted budgets); maxcost8: cost filter */
static void load_clusters(const char *path, int maxcost8, int maxdecl) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    size_t n0 = ncl, dropped = 0;
    for (;;) {
        uint8_t nwin;
        if (fread(&nwin, 1, 1, f) != 1) break;
        uint16_t ids[64]; uint8_t cnts[64];
        for (int i = 0; i < nwin; i++) { if (fread(&ids[i], 2, 1, f) != 1 || fread(&cnts[i], 1, 1, f) != 1) { fprintf(stderr, "truncated\n"); exit(1); } }
        uint8_t nd, nr;
        if (fread(&nd, 1, 1, f) != 1 || fread(&nr, 1, 1, f) != 1) { fprintf(stderr, "truncated\n"); exit(1); }
        uint8_t dv[3 * 16];
        if (nd && fread(dv, 1, 3 * nd, f) != (size_t)(3 * nd)) { fprintf(stderr, "truncated\n"); exit(1); }
        int cost8 = 0; for (int i = 0; i < nwin; i++) cost8 += cnts[i] * win[ids[i]].cost8;
        if (nd > maxdecl || cost8 > maxcost8) { dropped++; continue; }
        add_cluster(ids, cnts, nwin, nd);
    }
    fclose(f);
    fprintf(logf, "%s: %zu clusters loaded, %zu dropped (decl > %d or cost8 > %d)\n", path, ncl - n0, dropped, maxdecl, maxcost8);
}

static inline void cluster_sig(const Cl *c, int *sig) {
    memset(sig, 0, sizeof(int) * NK);
    for (int s = 0; s < c->nsig; s++) sig[(uint8_t)sarena[c->soff + 2 * s]] = sarena[c->soff + 2 * s + 1];
}

static inline void cluster_sig_add(const Cl *c, int m, int *sig) {
    for (int s = 0; s < c->nsig; s++) sig[(uint8_t)sarena[c->soff + 2 * s]] += m * sarena[c->soff + 2 * s + 1];
}

#endif
