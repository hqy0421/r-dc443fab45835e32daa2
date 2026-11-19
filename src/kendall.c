/*
 * Incremental sliding Kendall's tau-b for a single (already sorted) series.
 *
 * For every time point t the window is the half-open interval (t - width, t]
 * and only records whose two variables are both observed take part in the
 * window.  Records enter and leave the window in time order, so the pair
 * counts are maintained incrementally instead of being recomputed per window.
 *
 * The pair counts are kept in a sparse two-dimensional Fenwick tree over
 * compressed (x rank, y rank) coordinates: a Fenwick tree over x whose nodes
 * each hold a Fenwick tree over the y ranks that actually occur in them.  Each
 * update / query costs O(log^2 n) and the structure occupies O(n log n)
 * memory, giving O(n log^2 n) time and O(n log n) space for n records.
 *
 * Ties are handled exactly: a pair tied in x, in y or in both is neither
 * concordant nor discordant, and the tau-b denominator subtracts the tied
 * pairs on each side.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int lowbit(int i) {
  return i & -i;
}

static int cmp_int(const void *a, const void *b) {
  int x = *(const int *) a;
  int y = *(const int *) b;
  return (x > y) - (x < y);
}

/* number of elements in the sorted array a[0..n) that are <= v */
static int upper_bound_int(const int *a, int n, int v) {
  int lo = 0, hi = n;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (a[mid] <= v) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

typedef struct {
  int mx, my;
  int *ys_off;   /* start of each x node's sorted y ranks in `ys_val` */
  int *ys_len;   /* number of distinct y ranks in each x node */
  int *ys_val;   /* concatenated sorted y ranks */
  int *bit_off;  /* start of each x node's Fenwick storage in `bit` */
  int *bit;      /* concatenated Fenwick values (1..ys_len used) */
  int *bx;       /* 1-D Fenwick over x ranks */
  int *by;       /* 1-D Fenwick over y ranks */
  int *cx;       /* multiplicity of each x rank in the window */
  int *cy;       /* multiplicity of each y rank in the window */
  int nact;      /* number of records currently in the window */
  long long n_conc;  /* concordant pairs in the window */
  long long n_disc;  /* discordant pairs in the window */
  long long n_tie_x; /* pairs tied in x */
  long long n_tie_y; /* pairs tied in y */
} kstate;

static void bit1_add(int *bit, int n, int i, int d) {
  for (; i <= n; i += lowbit(i)) bit[i] += d;
}

static long long bit1_prefix(const int *bit, int i) {
  long long s = 0;
  for (; i > 0; i -= lowbit(i)) s += bit[i];
  return s;
}

static void bit2_add(kstate *s, int xr, int yr, int d) {
  int i;
  for (i = xr; i <= s->mx; i += lowbit(i)) {
    const int *v = s->ys_val + s->ys_off[i];
    int p = upper_bound_int(v, s->ys_len[i], yr);
    int base = s->bit_off[i];
    int j;
    for (j = p; j <= s->ys_len[i]; j += lowbit(j)) s->bit[base + j] += d;
  }
}

/* number of window records with x rank <= xr and y rank <= yr */
static long long bit2_prefix(const kstate *s, int xr, int yr) {
  long long r = 0;
  int i;
  if (xr <= 0 || yr <= 0) return 0;
  for (i = xr; i > 0; i -= lowbit(i)) {
    const int *v = s->ys_val + s->ys_off[i];
    int p = upper_bound_int(v, s->ys_len[i], yr);
    int base = s->bit_off[i];
    int j;
    for (j = p; j > 0; j -= lowbit(j)) r += s->bit[base + j];
  }
  return r;
}

/*
 * Pairs formed by the record (xr, yr) with the records currently in the
 * window.  The record itself may or may not be in the window: it contributes
 * zero to every term, so the same routine serves insertion and removal.
 */
static void ks_contrib(kstate *s, int xr, int yr,
                       long long *pc, long long *pd) {
  long long ax = bit1_prefix(s->bx, xr - 1);
  long long ay = bit1_prefix(s->by, yr - 1);
  long long qll = bit2_prefix(s, xr - 1, yr - 1);
  long long qxl = bit2_prefix(s, xr - 1, yr);
  long long qly = bit2_prefix(s, xr, yr - 1);
  long long qle = bit2_prefix(s, xr, yr);
  *pc = qll + ((long long) s->nact - (ax + s->cx[xr]) -
               (ay + s->cy[yr]) + qle);
  *pd = (ax - qxl) + (ay - qly);
}

static void ks_insert(kstate *s, int xr, int yr) {
  long long c, d;
  ks_contrib(s, xr, yr, &c, &d);
  s->n_conc += c;
  s->n_disc += d;
  s->n_tie_x += s->cx[xr];
  s->n_tie_y += s->cy[yr];
  s->cx[xr]++;
  s->cy[yr]++;
  bit2_add(s, xr, yr, 1);
  bit1_add(s->bx, s->mx, xr, 1);
  bit1_add(s->by, s->my, yr, 1);
  s->nact++;
}

static void ks_remove(kstate *s, int xr, int yr) {
  long long c, d;
  ks_contrib(s, xr, yr, &c, &d);
  s->n_conc -= c;
  s->n_disc -= d;
  s->cx[xr]--;
  s->cy[yr]--;
  s->n_tie_x -= s->cx[xr];
  s->n_tie_y -= s->cy[yr];
  bit2_add(s, xr, yr, -1);
  bit1_add(s->bx, s->mx, xr, -1);
  bit1_add(s->by, s->my, yr, -1);
  s->nact--;
}

/* dense ranks (1..*nuniq) of v[0..m), ties share a rank */
static void rankify(const double *v, int m, int *rank, int *nuniq) {
  double *s = (double *) R_alloc((size_t) (m > 0 ? m : 1), sizeof(double));
  int i, k = 0;
  if (m > 0) memcpy(s, v, (size_t) m * sizeof(double));
  R_rsort(s, m);
  for (i = 0; i < m; i++) {
    if (k == 0 || s[i] != s[k - 1]) s[k++] = s[i];
  }
  *nuniq = k;
  for (i = 0; i < m; i++) {
    int lo = 0, hi = k;
    while (lo < hi) {
      int mid = lo + (hi - lo) / 2;
      if (s[mid] < v[i]) lo = mid + 1;
      else hi = mid;
    }
    rank[i] = lo + 1;
  }
}

SEXP C_slide_kendall(SEXP time_s, SEXP x_s, SEXP y_s, SEXP width_s) {
  const double *time = REAL(time_s);
  const double *x = REAL(x_s);
  const double *y = REAL(y_s);
  const double width = REAL(width_s)[0];
  R_xlen_t n = XLENGTH(time_s);
  R_xlen_t k;
  int m = 0;
  int i, r;

  SEXP n_valid_s = PROTECT(allocVector(INTSXP, n));
  SEXP n_conc_s = PROTECT(allocVector(REALSXP, n));
  SEXP n_disc_s = PROTECT(allocVector(REALSXP, n));
  SEXP tau_s = PROTECT(allocVector(REALSXP, n));
  SEXP out = PROTECT(allocVector(VECSXP, 4));
  SEXP nm = PROTECT(allocVector(STRSXP, 4));
  int *n_valid = INTEGER(n_valid_s);
  double *n_conc = REAL(n_conc_s);
  double *n_disc = REAL(n_disc_s);
  double *tau = REAL(tau_s);

  for (k = 0; k < n; k++) {
    n_valid[k] = 0;
    n_conc[k] = 0;
    n_disc[k] = 0;
    tau[k] = NA_REAL;
  }

  SET_VECTOR_ELT(out, 0, n_valid_s);
  SET_VECTOR_ELT(out, 1, n_conc_s);
  SET_VECTOR_ELT(out, 2, n_disc_s);
  SET_VECTOR_ELT(out, 3, tau_s);
  SET_STRING_ELT(nm, 0, mkChar("n_valid"));
  SET_STRING_ELT(nm, 1, mkChar("n_concordant"));
  SET_STRING_ELT(nm, 2, mkChar("n_discordant"));
  SET_STRING_ELT(nm, 3, mkChar("tau_b"));
  setAttrib(out, R_NamesSymbol, nm);

  /* only records observed in both variables can enter a window */
  int *vidx = (int *) R_alloc((size_t) (n > 0 ? n : 1), sizeof(int));
  for (k = 0; k < n; k++) {
    if (!ISNAN(x[k]) && !ISNAN(y[k])) vidx[m++] = (int) k;
  }

  if (m == 0) {
    UNPROTECT(6);
    return out;
  }

  double *tv = (double *) R_alloc((size_t) m, sizeof(double));
  double *xv = (double *) R_alloc((size_t) m, sizeof(double));
  double *yv = (double *) R_alloc((size_t) m, sizeof(double));
  for (i = 0; i < m; i++) {
    tv[i] = time[vidx[i]];
    xv[i] = x[vidx[i]];
    yv[i] = y[vidx[i]];
  }

  int *rx = (int *) R_alloc((size_t) m, sizeof(int));
  int *ry = (int *) R_alloc((size_t) m, sizeof(int));
  int mx = 0, my = 0;
  rankify(xv, m, rx, &mx);
  rankify(yv, m, ry, &my);

  kstate s;
  memset(&s, 0, sizeof(s));
  s.mx = mx;
  s.my = my;

  /* count how many y ranks fall into each x node */
  int *cnt = (int *) R_alloc((size_t) mx + 1, sizeof(int));
  memset(cnt, 0, ((size_t) mx + 1) * sizeof(int));
  for (i = 0; i < m; i++) {
    int j;
    for (j = rx[i]; j <= mx; j += lowbit(j)) cnt[j]++;
  }

  int total = 0;
  int *ys_off = (int *) R_alloc((size_t) mx + 1, sizeof(int));
  for (i = 1; i <= mx; i++) {
    ys_off[i] = total;
    total += cnt[i];
  }
  int *ys_val = (int *) R_alloc((size_t) (total > 0 ? total : 1), sizeof(int));
  int *pos = (int *) R_alloc((size_t) mx + 1, sizeof(int));
  memcpy(pos, ys_off, ((size_t) mx + 1) * sizeof(int));
  for (i = 0; i < m; i++) {
    int j;
    for (j = rx[i]; j <= mx; j += lowbit(j)) ys_val[pos[j]++] = ry[i];
  }

  /* sort and de-duplicate each node, then lay out its Fenwick storage */
  int *ys_len = (int *) R_alloc((size_t) mx + 1, sizeof(int));
  memset(ys_len, 0, ((size_t) mx + 1) * sizeof(int));
  int btot = 0;
  int *bit_off = (int *) R_alloc((size_t) mx + 1, sizeof(int));
  for (i = 1; i <= mx; i++) {
    int a = ys_off[i];
    int len = cnt[i];
    int p, u = 0;
    if (len > 1) qsort(ys_val + a, (size_t) len, sizeof(int), cmp_int);
    for (p = 0; p < len; p++) {
      if (u == 0 || ys_val[a + p] != ys_val[a + u - 1]) {
        ys_val[a + u++] = ys_val[a + p];
      }
    }
    ys_len[i] = u;
    bit_off[i] = btot;
    btot += u + 1;
  }

  s.ys_off = ys_off;
  s.ys_len = ys_len;
  s.ys_val = ys_val;
  s.bit_off = bit_off;
  s.bit = (int *) R_alloc((size_t) (btot > 0 ? btot : 1), sizeof(int));
  memset(s.bit, 0, (size_t) (btot > 0 ? btot : 1) * sizeof(int));
  s.bx = (int *) R_alloc((size_t) mx + 1, sizeof(int));
  s.by = (int *) R_alloc((size_t) my + 1, sizeof(int));
  memset(s.bx, 0, ((size_t) mx + 1) * sizeof(int));
  memset(s.by, 0, ((size_t) my + 1) * sizeof(int));
  s.cx = (int *) R_alloc((size_t) mx + 1, sizeof(int));
  s.cy = (int *) R_alloc((size_t) my + 1, sizeof(int));
  memset(s.cx, 0, ((size_t) mx + 1) * sizeof(int));
  memset(s.cy, 0, ((size_t) my + 1) * sizeof(int));

  /* sweep the distinct time points of the full series */
  int lo_rec = 0, hi_rec = -1;
  r = 0;
  while (r < (int) n) {
    double t = time[r];
    int j = r;
    long long n0, d1, d2;
    double tv_window = NA_REAL;
    while (j + 1 < (int) n && time[j + 1] == t) j++;

    while (lo_rec <= hi_rec && tv[lo_rec] <= t - width) {
      ks_remove(&s, rx[lo_rec], ry[lo_rec]);
      lo_rec++;
    }
    while (hi_rec + 1 < m && tv[hi_rec + 1] <= t) {
      hi_rec++;
      ks_insert(&s, rx[hi_rec], ry[hi_rec]);
    }

    n0 = (long long) s.nact * (s.nact - 1) / 2;
    d1 = n0 - s.n_tie_x;
    d2 = n0 - s.n_tie_y;
    if (d1 > 0 && d2 > 0) {
      tv_window = (double) (s.n_conc - s.n_disc) /
        sqrt((double) d1 * (double) d2);
    }
    for (k = r; k <= j; k++) {
      n_valid[k] = s.nact;
      n_conc[k] = (double) s.n_conc;
      n_disc[k] = (double) s.n_disc;
      tau[k] = tv_window;
    }
    r = j + 1;
  }

  UNPROTECT(6);
  return out;
}
