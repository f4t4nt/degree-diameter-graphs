/* Voltage-lift search for small-diameter cubic graphs.
 *
 * usage: lift GROUPFILE BASEFILE D exhaust [MAXSOL] [OUTPREFIX]
 *        lift GROUPFILE BASEFILE D sample SECONDS SEED [MAXSOL] [OUTPREFIX]
 *
 * GROUPFILE: m, then m rows of the multiplication table (element 0 = identity),
 *            then naut, then naut rows: automorphisms of G as permutations.
 * BASEFILE:  k nslots ntree, then nslots lines "u v type" (type 0 edge, 1 loop,
 *            2 semi-edge); the first ntree slots are a spanning tree (voltage 0).
 *
 * The lift has vertex set {0..k-1} x G; (u,g) ~ (v, g*a) for a slot (u,v,a):
 * forward dart uses a, backward dart uses a^-1; a loop gives both; a semi-edge
 * needs an involution a.  Left multiplication by G acts as automorphisms, so
 * eccentricity is constant on fibres and k BFS runs (from (u,0)) decide the
 * diameter.  defect = sum over u of #vertices at distance > D from (u,0)
 * (= ordered pairs at distance > D divided by m).  In exhaust mode the BFS
 * stops at the first fibre with a defect, so the reported best_defect is exact
 * only when it is <= full_threshold (near misses are re-evaluated in full);
 * otherwise it is a lower bound.  The best assignment is printed as a BEST
 * line usable with the eval mode.
 *
 * exhaust: enumerates one representative of every Aut(G)-orbit of voltage
 * assignments (iterated stabilisers, lexicographically minimal representative;
 * loop voltages are also identified with their inverses).  Assignments giving
 * non-simple lifts (identity/involution on a loop, equal voltages on parallel
 * edges) are skipped.  sample: hill climbing with random restarts on defect.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXM 256
#define MAXK 32
#define MAXSLOTS 64
#define MAXN (MAXM * MAXK)

static int m, k, n, D;
static int16_t mul[MAXM * MAXM];
static int16_t inv_[MAXM];
static int isinv[MAXM]; /* 1 if non-identity involution */
static int naut;
static int16_t *aut; /* naut x m */

static int nslots, ntree;
static int su[MAXSLOTS], sv[MAXSLOTS], stype[MAXSLOTS];
static int volt[MAXSLOTS];
static int nfree, freeslot[MAXSLOTS];
static int bundle_prev[MAXSLOTS][MAXSLOTS], nbundle_prev[MAXSLOTS]; /* earlier slots parallel to this one */

/* darts: for each base vertex, up to 3 darts: target vertex, slot, direction (0 fwd, 1 bwd) */
static int ndart[MAXK], dv[MAXK][3], dslot[MAXK][3], ddir[MAXK][3];

static int16_t dist_[MAXN];
static int queue_[MAXN];
static int nbr_a[MAXK][3]; /* current group element used by each dart */

static long long evaluated = 0, solutions = 0;
static int best_defect = 1 << 30;
static int best_volt[MAXSLOTS]; /* voltages of the first assignment reaching best_defect */
static int maxsol = 5;
static const char *outprefix = "sol";

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void set_dart_elements(void) {
    for (int u = 0; u < k; u++)
        for (int j = 0; j < ndart[u]; j++) {
            int a = volt[dslot[u][j]];
            nbr_a[u][j] = ddir[u][j] ? inv_[a] : a;
        }
}

/* BFS from (u,0); returns number of vertices at distance > D.  If early, stops
 * at the first unreached vertex. */
static int bfs(int u, int early) {
    memset(dist_, -1, sizeof(int16_t) * n);
    int head = 0, tail = 0;
    int s = u * m;
    dist_[s] = 0;
    queue_[tail++] = s;
    int reached = 1;
    while (head < tail) {
        int x = queue_[head++];
        int dx = dist_[x];
        if (dx >= D) break;
        int bu = x / m, g = x % m;
        const int16_t *row = mul + g * m;
        for (int j = 0; j < ndart[bu]; j++) {
            int y = dv[bu][j] * m + row[nbr_a[bu][j]];
            if (dist_[y] < 0) {
                dist_[y] = dx + 1;
                queue_[tail++] = y;
                reached++;
            }
        }
        (void)early;
    }
    return n - reached;
}

static int defect(int early) {
    int d = 0;
    for (int u = 0; u < k; u++) {
        d += bfs(u, early);
        if (early && d) return d;
    }
    return d;
}

static void write_solution(void) {
    char path[512];
    snprintf(path, sizeof path, "%s_%lld.edges", outprefix, solutions);
    FILE *f = fopen(path, "w");
    if (!f) return;
    for (int x = 0; x < n; x++) {
        int bu = x / m, g = x % m;
        for (int j = 0; j < ndart[bu]; j++) {
            int y = dv[bu][j] * m + mul[g * m + nbr_a[bu][j]];
            if (x < y) fprintf(f, "%d %d\n", x, y);
        }
    }
    fclose(f);
    printf("SOLUTION %s voltages:", path);
    for (int s = 0; s < nslots; s++) printf(" %d", volt[s]);
    printf("\n");
    fflush(stdout);
}

static void note_best(int d) {
    if (d < best_defect) {
        best_defect = d;
        memcpy(best_volt, volt, sizeof(int) * nslots);
    }
}

/* free-slot voltages of the best assignment, in the order eval mode expects */
static void print_best(void) {
    if (best_defect >= (1 << 30)) return;
    printf("BEST defect=%d voltages:", best_defect);
    for (int i = 0; i < nfree; i++) printf(" %d", best_volt[freeslot[i]]);
    printf("\n");
}

/* In exhaust mode the early defect is only the first failing fibre's count (a
 * lower bound on the full defect); near misses up to full_threshold are
 * re-evaluated exactly so that best_defect is exact whenever it is <= full_threshold. */
static int full_threshold = 8;
static int near_threshold = 0; /* env LIFT_NEAR: print every assignment with exact defect <= this */

static void evaluate(int early) {
    set_dart_elements();
    int d = defect(early);
    evaluated++;
    if (early && d > 0 && d <= full_threshold) d = defect(0);
    note_best(d);
    if (d > 0 && d <= near_threshold) {
        printf("NEAR defect=%d voltages:", d);
        for (int i = 0; i < nfree; i++) printf(" %d", volt[freeslot[i]]);
        printf("\n");
    }
    if (d == 0) {
        solutions++;
        if (solutions <= maxsol) write_solution();
    }
}

/* allowed voltage for slot s given earlier assignments */
static int allowed(int s, int g) {
    if (stype[s] == 2) return isinv[g];
    if (stype[s] == 1) {
        if (g == 0 || isinv[g]) return 0;
        return 1;
    }
    for (int i = 0; i < nbundle_prev[s]; i++)
        if (volt[bundle_prev[s][i]] == g) return 0;
    return 1;
}

static int32_t *stab[MAXSLOTS + 1]; /* indices into aut; int32 because naut can exceed 32767 */
static int nstab[MAXSLOTS + 1];
static long long nodes_level[MAXSLOTS + 1];

static void rec(int level) {
    if (level == nfree) {
        evaluate(1);
        return;
    }
    int s = freeslot[level];
    int32_t *S = stab[level];
    int nS = nstab[level];
    int32_t *T = stab[level + 1];
    for (int g = 0; g < m; g++) {
        if (!allowed(s, g)) continue;
        /* minimality under current stabiliser (and inversion for loops) */
        int ok = 1, nT = 0;
        for (int i = 0; i < nS && ok; i++) {
            const int16_t *p = aut + (long)S[i] * m;
            int h = p[g];
            if (h < g) { ok = 0; break; }
            if (stype[s] == 1) {
                int hi = inv_[h];
                if (hi < g) { ok = 0; break; }
                if (h == g || hi == g) T[nT++] = S[i];
            } else if (h == g) T[nT++] = S[i];
        }
        if (!ok) continue;
        volt[s] = g;
        nstab[level + 1] = nT;
        nodes_level[level]++;
        rec(level + 1);
    }
    volt[s] = 0;
}

static uint64_t rng_state = 88172645463325252ULL;
static uint64_t rnd(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static int random_allowed(int s) {
    for (int tries = 0; tries < 1000; tries++) {
        int g = rnd() % m;
        if (allowed(s, g)) return g;
    }
    for (int g = 0; g < m; g++)
        if (allowed(s, g)) return g;
    return -1;
}

static void sample(double seconds) {
    double t0 = now();
    long long restarts = 0;
    int cur[MAXSLOTS];
    while (now() - t0 < seconds) {
        /* random start */
        for (int i = 0; i < nfree; i++) {
            int s = freeslot[i];
            volt[s] = 0;
        }
        for (int i = 0; i < nfree; i++) {
            int s = freeslot[i];
            int g = random_allowed(s);
            if (g < 0) { printf("no allowed voltage for slot %d\n", s); return; }
            volt[s] = g;
        }
        set_dart_elements();
        int d = defect(0);
        evaluated++;
        note_best(d);
        if (d == 0) { solutions++; if (solutions <= maxsol) write_solution(); }
        int stale = 0;
        while (d > 0 && stale < 3 * nfree && now() - t0 < seconds) {
            /* best single-voltage change on a random slot (ties accepted at random) */
            int i = rnd() % nfree;
            int s = freeslot[i];
            int old = volt[s];
            int bestg = old, bestd = d, nties = 0;
            for (int g = 0; g < m; g++) {
                if (g == old || !allowed(s, g)) continue;
                volt[s] = g;
                set_dart_elements();
                int dd = defect(0);
                evaluated++;
                if (dd < bestd) { bestd = dd; bestg = g; nties = 1; }
                else if (dd == bestd && bestd <= d) { nties++; if (rnd() % nties == 0) bestg = g; }
            }
            volt[s] = bestg;
            if (bestd < d) stale = 0; else stale++;
            d = bestd;
            note_best(d);
            if (d == 0) {
                set_dart_elements();
                solutions++;
                if (solutions <= maxsol) write_solution();
                break;
            }
        }
        memcpy(cur, volt, sizeof cur);
        restarts++;
    }
    printf("restarts %lld\n", restarts);
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: %s GROUPFILE BASEFILE D exhaust [MAXSOL] [OUTPREFIX] | sample SECONDS SEED [MAXSOL] [OUTPREFIX]\n", argv[0]);
        return 1;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror(argv[1]); return 1; }
    if (fscanf(f, "%d", &m) != 1 || m > MAXM) { fprintf(stderr, "bad group\n"); return 1; }
    for (int i = 0; i < m * m; i++) { int x; if (fscanf(f, "%d", &x) != 1) return 1; mul[i] = x; }
    if (fscanf(f, "%d", &naut) != 1) return 1;
    aut = malloc(sizeof(int16_t) * (size_t)naut * m);
    for (long i = 0; i < (long)naut * m; i++) { int x; if (fscanf(f, "%d", &x) != 1) return 1; aut[i] = x; }
    fclose(f);
    for (int a = 0; a < m; a++) {
        for (int b = 0; b < m; b++) if (mul[a * m + b] == 0) inv_[a] = b;
        isinv[a] = (a != 0 && mul[a * m + a] == 0);
    }

    f = fopen(argv[2], "r");
    if (!f) { perror(argv[2]); return 1; }
    if (fscanf(f, "%d %d %d", &k, &nslots, &ntree) != 3 || k > MAXK || nslots > MAXSLOTS) { fprintf(stderr, "bad base\n"); return 1; }
    for (int s = 0; s < nslots; s++) {
        if (fscanf(f, "%d %d %d", &su[s], &sv[s], &stype[s]) != 3) return 1;
    }
    fclose(f);
    n = k * m;
    if (n > MAXN) { fprintf(stderr, "too big\n"); return 1; }
    D = atoi(argv[3]);

    /* darts */
    for (int u = 0; u < k; u++) ndart[u] = 0;
    for (int s = 0; s < nslots; s++) {
        int u = su[s], v = sv[s];
        if (stype[s] == 0) {
            dv[u][ndart[u]] = v; dslot[u][ndart[u]] = s; ddir[u][ndart[u]] = 0; ndart[u]++;
            dv[v][ndart[v]] = u; dslot[v][ndart[v]] = s; ddir[v][ndart[v]] = 1; ndart[v]++;
        } else if (stype[s] == 1) {
            dv[u][ndart[u]] = u; dslot[u][ndart[u]] = s; ddir[u][ndart[u]] = 0; ndart[u]++;
            dv[u][ndart[u]] = u; dslot[u][ndart[u]] = s; ddir[u][ndart[u]] = 1; ndart[u]++;
        } else {
            dv[u][ndart[u]] = u; dslot[u][ndart[u]] = s; ddir[u][ndart[u]] = 0; ndart[u]++;
        }
    }
    for (int u = 0; u < k; u++) if (ndart[u] != 3) { fprintf(stderr, "vertex %d has degree %d\n", u, ndart[u]); return 1; }
    nfree = 0;
    for (int s = ntree; s < nslots; s++) freeslot[nfree++] = s;
    for (int s = 0; s < nslots; s++) {
        nbundle_prev[s] = 0;
        if (stype[s] != 0) continue;
        for (int t = 0; t < s; t++)
            if (stype[t] == 0 && ((su[t] == su[s] && sv[t] == sv[s]) || (su[t] == sv[s] && sv[t] == su[s])))
                bundle_prev[s][nbundle_prev[s]++] = t;
    }
    for (int s = 0; s < nslots; s++) volt[s] = 0;

    int ninv = 0;
    for (int a = 0; a < m; a++) ninv += isinv[a];
    if (getenv("LIFT_NEAR")) near_threshold = atoi(getenv("LIFT_NEAR"));
    if (near_threshold > full_threshold) full_threshold = near_threshold;
    double t0 = now();
    if (strcmp(argv[4], "exhaust") == 0) {
        if (argc > 5) maxsol = atoi(argv[5]);
        if (argc > 6) outprefix = argv[6];
        for (int l = 0; l <= nfree; l++) stab[l] = malloc(sizeof(int32_t) * (naut > 0 ? naut : 1));
        for (int i = 0; i < naut; i++) stab[0][i] = i;
        nstab[0] = naut;
        rec(0);
        double dt = now() - t0;
        printf("RESULT mode=exhaust m=%d k=%d n=%d D=%d naut=%d ninv=%d nfree=%d evaluated=%lld solutions=%lld best_defect=%d seconds=%.2f\n",
               m, k, n, D, naut, ninv, nfree, evaluated, solutions, best_defect, dt);
        print_best();
    } else if (strcmp(argv[4], "sample") == 0) {
        if (argc < 7) { fprintf(stderr, "sample needs SECONDS SEED\n"); return 1; }
        double seconds = atof(argv[5]);
        rng_state ^= (uint64_t)atoll(argv[6]) * 0x9E3779B97F4A7C15ULL;
        for (int i = 0; i < 10; i++) rnd();
        if (argc > 7) maxsol = atoi(argv[7]);
        if (argc > 8) outprefix = argv[8];
        sample(seconds);
        double dt = now() - t0;
        printf("RESULT mode=sample m=%d k=%d n=%d D=%d naut=%d ninv=%d nfree=%d evaluated=%lld solutions=%lld best_defect=%d seconds=%.2f\n",
               m, k, n, D, naut, ninv, nfree, evaluated, solutions, best_defect, dt);
        print_best();
    } else if (strcmp(argv[4], "eval") == 0) {
        /* eval v_ntree ... v_{nslots-1}: evaluate one assignment and write it */
        if (argc != 5 + nfree) { fprintf(stderr, "eval needs %d voltages\n", nfree); return 1; }
        for (int i = 0; i < nfree; i++) volt[freeslot[i]] = atoi(argv[5 + i]);
        set_dart_elements();
        int d = defect(0);
        printf("defect=%d (ordered pairs at distance > D = %d)\n", d, d * m);
        solutions = 0;
        write_solution();
    } else {
        fprintf(stderr, "unknown mode\n");
        return 1;
    }
    return 0;
}
