#include "func.h"

/* library routines */
double agexp(double r, double x, double y)
{
   double w;

   if (r == 0.0)
      return (exp(2 * x));
   else {
      x = 1 + r * x;
      y = r * y;
      w = x * x + y * y;
      if (r < 0.0)
         return (pow(1 / w, -1 / r));
      else
         return (pow(w, 1 / r));
   }
}

int cholesky(double *c, double *a, double *b, const int n, double eps)
{
   int i, j, k;
   static double *d = NULL, *y, *v, *vp;
   static int size;

   if (d == NULL) {
      d = dgetmem(n * (n + 2));
      y = d + n;
      v = y + n;
      size = n;
   }

   if (n > size) {
      free(d);
      d = dgetmem(n * (n + 2));
      y = d + n;
      v = y + n;
      size = n;
   }

   if (eps < 0.0)
      eps = 1.0e-6;

   for (j = 0; j < n; j++, c += n) {
      d[j] = c[j];
      vp = v + j * n;
      for (k = 0; k < j; k++)
         d[j] -= vp[k] * vp[k] * d[k];

      if (fabs(d[j]) <= eps)
         return (-1);

      for (i = j + 1; i < n; i++) {
         vp = v + i * n;
         vp[j] = c[i];
         for (k = 0; k < j; k++)
            vp[j] -= vp[k] * v[j * n + k] * d[k];
         vp[j] /= d[j];
      }
   }

   for (i = 0; i < n; i++) {
      y[i] = b[i];
      vp = v + i * n;
      for (k = 0; k < i; k++)
         y[i] -= vp[k] * y[k];
   }

   for (i = n - 1; i >= 0; i--) {
      a[i] = y[i] / d[i];
      for (k = i + 1; k < n; k++)
         a[i] -= v[n * k + i] * a[k];
   }
   return (0);
}

/* freada: read ascii */
int freada(double *p, const int bl, FILE * fp)
{
   int c;
   char buf[LINEBUFSIZE];

#if defined(WIN32)
   _setmode(_fileno(fp), _O_TEXT);
#endif

   c = 0;
   while (c < bl) {
      if (fgets(buf, LINEBUFSIZE, fp) == NULL)
         break;
      p[c] = atof(buf);
      c++;
   }
   return (c);
}

/* fritex: wrapper function for fwrite */
int fwritex(void *ptr, const size_t size, const int nitems, FILE * fp)
{
#if defined(WIN32)
   _setmode(_fileno(fp), _O_BINARY);
#endif
   return (fwrite(ptr, size, nitems, fp));
}

/* freadx: wrapper function for fread */
int freadx(void *ptr, const size_t size, const int nitems, FILE * fp)
{
#if defined(WIN32)
   _setmode(_fileno(fp), _O_BINARY);
#endif
   return (fread(ptr, size, nitems, fp));
}

/* --------------- double I/O compile --------------- */
#ifdef DOUBLE
/* fwritef : write double type data */
int fwritef(double *ptr, const size_t size, const int nitems, FILE * fp)
{
   return (fwritex(ptr, size, nitems, fp));
}

/* freadf : read double type data */
int freadf(double *ptr, const size_t size, const int nitems, FILE * fp)
{
   return (freadx(ptr, size, nitems, fp));
}

#else                           /* DOUBLE */
/* --------------- float I/O compile --------------- */

static float *f;
static int items;

/* fwritef : convert double type data to float type and write */
int fwritef(double *ptr, const size_t size, const int nitems, FILE * fp)
{
   int i;
   if (items < nitems) {
      if (f != NULL)
         free(f);
      items = nitems;
      f = fgetmem(items);
   }
   for (i = 0; i < nitems; i++)
      f[i] = ptr[i];

#if defined(WIN32)
   _setmode(_fileno(fp), _O_BINARY);
#endif

   return fwrite(f, sizeof(float), nitems, fp);
}

/* freadf : read float type data and convert to double type */
int freadf(double *ptr, const size_t size, const int nitems, FILE * fp)
{
   int i, n;
   if (items < nitems) {
      if (f != NULL)
         free(f);
      items = nitems;
      f = fgetmem(items);
   }
#if defined(WIN32)
   _setmode(_fileno(fp), _O_BINARY);
#endif

   n = fread(f, sizeof(float), nitems, fp);
   for (i = 0; i < n; i++)
      ptr[i] = f[i];

   return n;
}
#endif                          /* DOUBLE */

void SPTK_byte_swap(void *p, size_t size, size_t num)
{
   char *q, tmp;
   size_t i, j;

   q = (char *) p;

   for (i = 0; i < num; i++) {
      for (j = 0; j < (size / 2); j++) {
         tmp = *(q + j);
         *(q + j) = *(q + (size - 1 - j));
         *(q + (size - 1 - j)) = tmp;
      }
      q += size;
   }
}

int fwrite_little_endian(void *buf, const size_t size,
                         const size_t n, FILE * fp)
{
#ifdef WORDS_BIGENDIAN
   SPTK_byte_swap(buf, size, n);
#endif
   return fwrite(buf, size, n, fp);
}

void fillz(void *ptr, const size_t size, const int nitem)
{
   long n;
   char *p = (char *)ptr;

   n = size * nitem;
   while (n--)
      *p++ = '\0';
}

FILE *getfp(char *name, char *opt)
{
   FILE *fp;

   if ((fp = fopen(name, opt)) == NULL) {
      fprintf(stderr, "Cannot open file %s!\n", name);
      exit(2);
   }

   return (fp);
}

char *getmem(const size_t leng, const size_t size)
{
   char *p = NULL;

   if ((p = (char *) calloc(leng, size)) == NULL) {
      fprintf(stderr, "Cannot allocate memory!\n");
      exit(3);
   }
   return (p);
}

short *sgetmem(const int leng)
{
   return ((short *) getmem((size_t) leng, sizeof(short)));
}

long *lgetmem(const int leng)
{
   return ((long *) getmem((size_t) leng, sizeof(long)));
}

double *dgetmem(const int leng)
{
   return ((double *) getmem((size_t) leng, sizeof(double)));
}

float *fgetmem(const int leng)
{
   return ((float *) getmem((size_t) leng, sizeof(float)));
}

real *rgetmem(const int leng)
{
   return ((real *) getmem((size_t) leng, sizeof(real)));
}

float **ffgetmem(const int leng)
{
   return ((float **) getmem((size_t) leng, sizeof(float *)));
}

double **ddgetmem(const int leng1, const int leng2)
{
   int i, j;
   double **tmp, *tmp2;
   tmp = (double **) getmem((size_t) leng1, sizeof(double *));
   tmp2 = dgetmem(leng1 * leng2);

   for (i = 0, j = 0; i < leng1; i++, j += leng2) {
      tmp[i] = tmp2 + j;
   }

   return (tmp);
}

double gexp(const double r, const double x)
{
   if (r == 0.0) {
      return (exp(x));
   } else {
      if (r < 0.0)
         return (pow(1 / (1 + r * x), -1 / r));
      else
         return (pow(1 + r * x, 1 / r));
   }
}

double glog(const double r, const double x)
{
   if (r == 0.0) {
      return (log(x));
   } else {
      if (r < 0.0)
         return ((pow(1 / x, -r) - 1.0) / r);
      else
         return ((pow(x, r) - 1.0) / r);
   }
}

int ifftr(double *x, double *y, const int l)
{
   int i;
   double *xp, *yp;

   fftr(x, y, l);

   xp = x;
   yp = y;
   i = l;
   while (i--) {
      *xp++ /= l;
      *yp++ /= -l;
   }

   return (0);
}

double invert(double **mat, double **inv, const int n)
{
   int i, j, k, *swap, ii, ik;
   double **copy_mat, *tmpmat, d, u, det, *work;

   copy_mat = (double **) malloc(sizeof(double *) * n);
   tmpmat = dgetmem(n * n);

   for (i = 0, j = 0; i < n; i++, j += n) {
      copy_mat[i] = tmpmat + j;
   }
   for (i = 0; i < n; i++) {
      for (k = 0; k < n; k++) {
         copy_mat[i][k] = mat[i][k];
      }
   }

   swap = (int *) malloc(sizeof(int) * n);
   work = dgetmem(n);

   for (k = 0; k < n; k++) {
      swap[k] = k;
      u = 0.0;
      for (j = 0; j < n; j++) {
         d = fabs(copy_mat[k][j]);
         if (d > u) {
            u = d;
         }
      }
      if (u == 0.0) {
         fprintf(stderr, "Can't calculate inverse matrix!\n");
         exit(1);
      }
      work[k] = 1.0 / u;
   }

   det = 1;
   for (k = 0; k < n; k++) {
      u = -1;
      for (i = k; i < n; i++) {
         ii = swap[i];
         d = fabs(copy_mat[ii][k]) * work[ii];
         if (d > u) {
            u = d;
            j = i;
         }
      }

      ik = swap[j];
      if (j != k) {
         swap[j] = swap[k];
         swap[k] = ik;
         det = -det;
      }

      u = copy_mat[ik][k];
      det *= u;
      if (u == 0.0) {
         fprintf(stderr, "Can't calculate inverse matrix!\n");
         exit(1);
      }
      for (i = k + 1; i < n; i++) {
         ii = swap[i];
         d = (copy_mat[ii][k] /= u);
         for (j = k + 1; j < n; j++) {
            copy_mat[ii][j] -= d * copy_mat[ik][j];
         }
      }
   }

   if (det != 0.0) {
      for (k = 0; k < n; k++) {
         for (i = 0; i < n; i++) {
            ii = swap[i];
            d = (ii == k);
            for (j = 0; j < i; j++) {
               d -= copy_mat[ii][j] * inv[j][k];
            }
            inv[i][k] = d;
         }
         for (i = n - 1; i >= 0; i--) {
            d = inv[i][k];
            ii = swap[i];
            for (j = i + 1; j < n; j++) {
               d -= copy_mat[ii][j] * inv[j][k];
            }
            inv[i][k] = d / copy_mat[ii][i];
         }
      }
   } else {
      fprintf(stderr, "Can't calculate inverse matrix!\n");
      exit(1);
   }

   free(copy_mat[0]);
   free(copy_mat);
   free(swap);
   free(work);

   return (det);
}

static double *tmp;
static int tmpsize = 0;

static void mm(double x[], const int xx, const int xy, double y[], const int yx,
               const int yy, double a[])
{
   int i, j, k;
   double *wx, *wy;

   if (xx == 1 && xy == 1) {
      for (i = yx * yy - 1; i >= 0; i--)
         a[i] = x[0] * y[i];
      return;
   }

   if (xx != yy) {
      fprintf(stderr, "Invalid matrix size x= %d*%d,y= %d*%d\n", xx, xy, yx,
              yy);
      exit(1);
   }

   wx = x;
   for (i = 0; i < xy; i++) {
      for (j = 0; j < yx; j++) {
         wy = &y[j];
         *a = 0;
         for (k = 0; k < xx; k++) {
            *a += *wx * *wy;
            wx++;
            wy += yx;
         }
         wx -= xx;
         a++;
      }
      wx += xx;
   }

   return;
}

void multim(double x[], const int xx, const int xy, double y[], const int yx,
            const int yy, double a[])
{
   int i;

   if (x == a) {
      if (((xy > yy) ? xy : yy) * yx > tmpsize) {
         if (tmp != NULL)
            free(tmp);
         tmpsize = ((xy > yy) ? xy : yy) * yx;
         tmp = (double *) getmem(tmpsize, sizeof(*tmp));
      }
      mm(x, xx, xy, y, yx, yy, tmp);
      if (xx == xy)
         for (i = yx * yy - 1; i >= 0; i--)
            a[i] = tmp[i];
      else
         for (i = xy * yx - 1; i >= 0; i--)
            a[i] = tmp[i];
   } else {
      mm(x, xx, xy, y, yx, yy, a);
   }

   return;
}

static void am(double x[], double y[], const int xx, const int yy, double a[])
{
   int i, j;

   for (i = 0; i < yy; i++)
      for (j = 0; j < xx; j++)
         a[j + i * xx] = x[j + i * xx] + y[j + i * xx];
}

void addm(double x[], double y[], const int xx, const int yy, double a[])
{
   int i;

   if (x == a) {
      if (xx * yy > tmpsize) {
         if (tmp != NULL)
            free(tmp);
         tmpsize = xx * yy;
         tmp = (double *) getmem(tmpsize, sizeof(*tmp));
      }
      am(x, y, xx, yy, tmp);
      for (i = xx * yy - 1; i >= 0; i--)
         a[i] = tmp[i];
   } else {
      am(x, y, xx, yy, a);
   }

   return;
}

void movem(void *a, void *b, const size_t size, const int nitem)
{
   long i;
   char *c = (char *)a;
   char *d = (char *)b;

   i = size * nitem;
   if (c > d)
      while (i--)
         *d++ = *c++;
   else {
      c += i;
      d += i;
      while (i--)
         *--d = *--c;
   }
}

int mseq(void)
{
   static int x = 0x55555555;
   int x0, x28;

   x >>= 1;

   if (x & B0)
      x0 = 1;
   else
      x0 = -1;

   if (x & B28)
      x28 = 1;
   else
      x28 = -1;

   if (x0 + x28)
      x &= B31_;
   else
      x |= B31;

   return (x0);
}

static void mv_mul(double *t, double *x, double *y)
{
   t[0] = x[0] * y[0] + x[1] * y[1];
   t[1] = x[2] * y[0] + x[3] * y[1];

   return;
}

static void mm_mul(double *t, double *x, double *y)
{
   t[0] = x[0] * y[0] + x[1] * y[2];
   t[1] = x[0] * y[1] + x[1] * y[3];
   t[2] = x[2] * y[0] + x[3] * y[2];
   t[3] = x[2] * y[1] + x[3] * y[3];

   return;
}

static int inverse(double *x, double *y, const double eps)
{
   double det;

   det = y[0] * y[3] - y[1] * y[2];

#ifdef WIN32
   if ((fabs(det) < eps) || _isnan(det)) {
#else
   if ((fabs(det) < eps) || isnan(det)) {
#endif
      fprintf(stderr,
              "theq() : determinant of the normal matrix is too small!\n");
      return (-1);
   }

   x[0] = y[3] / det;
   x[1] = -y[1] / det;
   x[2] = -y[2] / det;
   x[3] = y[0] / det;

   return (0);
}

static void crstrns(double *x, double *y)
{
   x[0] = y[3];
   x[1] = y[2];
   x[2] = y[1];
   x[3] = y[0];

   return;
}

static double **mtrx2(const int a, const int b)
{
   int i;
   double **x;

   if (!(x = (double **) calloc((size_t) a, sizeof(*x)))) {
      fprintf(stderr, "mtrx2() in theq() : Cannot allocate memory!\n");
      exit(3);
   }
   for (i = 0; i < a; i++)
      if (!(x[i] = (double *) calloc((size_t) b, sizeof(**x)))) {
         fprintf(stderr, "mtrx2() in theq() : Cannot allocate memory!\n");
         exit(3);
      }

   return (x);
}

static int cal_p0(double **p, double **r, double *b, const int n,
                  const double eps)
{
   double t[4], s[2];

   if (inverse(t, r[0], eps) == -1)
      return (-1);
   s[0] = b[0];
   s[1] = b[n - 1];
   mv_mul(p[0], t, s);

   return (0);
}

static void cal_ex(double *ex, double **r, double **x, const int i)
{
   int j;
   double t[4], s[4];

   s[0] = s[1] = s[2] = s[3] = 0.;

   for (j = 0; j < i; j++) {
      mm_mul(t, r[i - j], x[j]);
      s[0] += t[0];
      s[1] += t[1];
      s[2] += t[2];
      s[3] += t[3];
   }

   ex[0] = s[0];
   ex[1] = s[1];
   ex[2] = s[2];
   ex[3] = s[3];

   return;
}

static void cal_ep(double *ep, double **r, double **p, const int i)
{
   int j;
   double t[2], s[2];

   s[0] = s[1] = 0.;

   for (j = 0; j < i; j++) {
      mv_mul(t, r[i - j], p[j]);
      s[0] += t[0];
      s[1] += t[1];
   }
   ep[0] = s[0];
   ep[1] = s[1];

   return;
}

static int cal_bx(double *bx, double *vx, double *ex, const double eps)
{
   double t[4], s[4];

   crstrns(t, vx);
   if (inverse(s, t, eps) == -1)
      return (-1);
   mm_mul(bx, s, ex);

   return (0);
}

static void cal_x(double **x, double **xx, double *bx, const int i)
{
   int j;
   double t[4], s[4];

   for (j = 1; j < i; j++) {
      crstrns(t, xx[i - j]);
      mm_mul(s, t, bx);
      x[j][0] -= s[0];
      x[j][1] -= s[1];
      x[j][2] -= s[2];
      x[j][3] -= s[3];
   }

   for (j = 1; j < i; j++) {
      xx[j][0] = x[j][0];
      xx[j][1] = x[j][1];
      xx[j][2] = x[j][2];
      xx[j][3] = x[j][3];
   }

   x[i][0] = xx[i][0] = -bx[0];
   x[i][1] = xx[i][1] = -bx[1];
   x[i][2] = xx[i][2] = -bx[2];
   x[i][3] = xx[i][3] = -bx[3];

   return;
}

static void cal_vx(double *vx, double *ex, double *bx)
{
   double t[4], s[4];

   crstrns(t, ex);
   mm_mul(s, t, bx);
   vx[0] -= s[0];
   vx[1] -= s[1];
   vx[2] -= s[2];
   vx[3] -= s[3];

   return;
}

static int cal_g(double *g, double *vx, double *b, double *ep,
                 const int i, const int n, const double eps)
{
   double t[2], s[4], u[4];

   t[0] = b[i] - ep[0];
   t[1] = b[n - 1 - i] - ep[1];
   crstrns(s, vx);

   if (inverse(u, s, eps) == -1)
      return (-1);
   mv_mul(g, u, t);

   return (0);
}

static void cal_p(double **p, double **x, double *g, const int i)
{
   double t[4], s[2];
   int j;

   for (j = 0; j < i; j++) {
      crstrns(t, x[i - j]);
      mv_mul(s, t, g);
      p[j][0] += s[0];
      p[j][1] += s[1];
   }

   p[i][0] = g[0];
   p[i][1] = g[1];

   return;
}

int theq(double *t, double *h, double *a, double *b, const int n, double eps)
{
   static double **r = NULL, **x, **xx, **p;
   static int size;
   double ex[4], ep[2], vx[4], bx[4], g[2];
   int i;

   if (r == NULL) {
      r = mtrx2(n, 4);
      x = mtrx2(n, 4);
      xx = mtrx2(n, 4);
      p = mtrx2(n, 2);
      size = n;
   }
   if (n > size) {
      for (i = 0; i < size; i++) {
         free((char *) r[i]);
         free((char *) x[i]);
         free((char *) xx[i]);
         free((char *) p[i]);
      }
      free((char *) r);
      free((char *) x);
      free((char *) xx);
      free((char *) p);

      r = mtrx2(n, 4);
      x = mtrx2(n, 4);
      xx = mtrx2(n, 4);
      p = mtrx2(n, 2);
      size = n;
   }

   if (eps < 0.0)
      eps = 1.0e-6;

   /* make r */
   for (i = 0; i < n; i++) {
      r[i][0] = r[i][3] = t[i];
      r[i][1] = h[n - 1 + i];
      r[i][2] = h[n - 1 - i];
   }

   /* step 1 */
   x[0][0] = x[0][3] = 1.0;
   if (cal_p0(p, r, b, n, eps) == -1)
      return (-1);

   vx[0] = r[0][0];
   vx[1] = r[0][1];
   vx[2] = r[0][2];
   vx[3] = r[0][3];

   /* step 2 */
   for (i = 1; i < n; i++) {
      cal_ex(ex, r, x, i);
      cal_ep(ep, r, p, i);
      if (cal_bx(bx, vx, ex, eps) == -1)
         return (-1);
      cal_x(x, xx, bx, i);
      cal_vx(vx, ex, bx);
      if (cal_g(g, vx, b, ep, i, n, eps) == -1)
         return (-1);
      cal_p(p, x, g, i);
   }

   /* step 3 */
   for (i = 0; i < n; i++)
      a[i] = p[i][0];

   return (0);
}

int toeplitz(double *t, double *a, double *b, const int n, double eps)
{
   int l, k;
   static double *c = NULL, *cc;
   static int size;
   double rmd, mue, mue2;

   if (c == NULL) {
      c = dgetmem(n + n + 2);
      cc = c + n;
      size = n;
   }
   if (n > size) {
      free(c);
      c = dgetmem(n + n + 2);
      cc = c + n;
      size = n;
   }

   if (eps < 0.0)
      eps = 1.0e-6;

   fillz(c, sizeof(*c), n + 1);

   rmd = t[0];
   if (((rmd < 0.0) ? -rmd : rmd) <= eps)
      return (-1);

   a[0] = b[0] / rmd;

   for (l = 1; l < n; l++) {
      mue = -t[l];
      for (k = 1; k < l; k++)
         mue -= c[k] * t[l - k];
      mue /= rmd;

      for (k = 1; k < l; k++)
         cc[k] = c[k] + mue * c[l - k];
      cc[l] = mue;

      rmd = (1.0 - mue * mue) * rmd;
      if (((rmd < 0.0) ? -rmd : rmd) <= eps)
         return (-1);

      for (k = 1; k <= l; k++)
         c[k] = cc[k];

      mue2 = b[l];
      for (k = 0; k <= l - 1; k++)
         mue2 += c[l - k] * b[k];
      mue2 /= rmd;

      for (k = 0; k < l; k++)
         a[k] += mue2 * c[l - k];
      a[l] = mue2;
   }

   return (0);
}

/* tool routines */

double acep(double x, double *c, const int m, const double lambda,
            const double step, const double tau, const int pd, const double eps)
{
   int i, memory_size;
   static double *cc = NULL, *e, *ep, *d, gg = 1.0;
   static int size;
   double mu, tx;

   memory_size = m + m + m + 3 + (m + 1) * pd * 2;
   if (cc == NULL) {
      cc = dgetmem(memory_size);
      e = cc + m + 1;
      ep = e + m + 1;
      d = ep + m + 1;
      size = memory_size;
   }

   if (memory_size > size) {
      free(cc);
      cc = dgetmem(memory_size);
      e = cc + m + 1;
      ep = e + m + 1;
      d = ep + m + 1;
      size = memory_size;
   }

   for (i = 1; i <= m; i++)
      cc[i] = -c[i];

   x = lmadf(x, cc, m, pd, d);

   for (i = m; i >= 1; i--)
      e[i] = e[i - 1];
   e[0] = x;

   gg = gg * lambda + (1.0 - lambda) * e[0] * e[0];
   c[0] = 0.5 * log(gg);

   gg = (gg < eps) ? eps : gg;
   mu = step / (double) m / gg;
   tx = 2 * (1.0 - tau) * x;

   for (i = 1; i <= m; i++) {
      ep[i] = tau * ep[i] - tx * e[i];
      c[i] -= mu * ep[i];
   }

   return (x);
}

void acorr(double *x, int l, double *r, const int np)
{
   double d;
   int k, i;

   for (k = 0; k <= np; k++) {
      for (d = 0.0, i = 0; i < l - k; i++)
         d += x[i] * x[i + k];
      r[k] = d;
   }

   return;
}

void quicksort(double *x, int left, int right)
{
   int i, j;
   double pivot;
   double tmp;

   i = left;
   j = right;

   pivot = x[(left + right) / 2];

   while (1) {

      while (x[i] > pivot)
         i++;

      while (pivot > x[j])
         j--;
      if (i >= j)
         break;

      tmp = x[i];
      x[i] = x[j];
      x[j] = tmp;

      i++;
      j--;
   }

   if (left < i - 1)
      quicksort(x, left, i - 1);
   if (j + 1 < right)
      quicksort(x, j + 1, right);
}

int vander(double *c, double *a, double *b, const int n, double eps) {
   int i, j;
   double e, numer, denom;
   static double *d = NULL;
   static int size;

   if (d == NULL) {
      d = dgetmem(n);
      size = n;
   }

   if (n > size) {
      free(d);
      d = dgetmem(n);
      size = n;
   }

   if (eps < 0.0)
      eps = 1.0e-6;

   fillz(d, sizeof(*d), n);
   for (j = 0; j < n; j++) {
      for (i = n - j - 1; i < n - 1; i++) {
         d[i] += -c[j] * d[i + 1];
      }
      d[n - 1] += -c[j];
   }

   for (j = 0; j < n; j++) {
      e = 1.0;
      numer = b[n - 1];
      denom = 1.0;
      for (i = n - 2; i >= 0; i--) {
         e = d[i + 1] + c[j] * e;
         numer = numer + b[i] * e;
         denom = denom * c[j] + e;
      }

      if (fabs(denom) <= eps)
         return -1;

      a[j] = numer / denom;
   }

   return 0;
}

unsigned long long nck(int n, int k)
{
   int i;
   unsigned long long p = 1;

   if (2 * k > n)
      k = n - k;
   for (i = 1; i <= k; i++) {
      p *= n--;
      p /= i;
   }

   return p;
}

int acr2csm(double *r, double *csm, const int m)
{
   int i, k, l, n = (m + 1) / 2;
   double sum, tmp;
   static complex *z = NULL;
   static double *u = NULL, *h, *p, *x, *q;
   static int size;

   if (u == NULL || x == NULL) {
      u = dgetmem(m + 1 + n * n + n + 1 + n + n);
      h = u + m + 1;
      p = h + n * n;
      x = p + n + 1;
      q = x + n;
      z = cplx_getmem(n + 1);
      size = m;
   }

   if (m > size) {
      free(u);
      free(z);
      u = dgetmem(m + 1 + n * n + n + 1 + n + n);
      h = u + m + 1;
      p = h + n * n;
      x = p + n + 1;
      q = x + n;
      z = cplx_getmem(n + 1);
      size = m;
   }

   for (l = 0; l <= m; l++) {
      sum = 0.0;
      for (k = 0; k <= l; k++) {
         sum += nck(l, k) * r[abs(2 * k - l)];
      }
      u[l] = sum / pow(2.0, l);
   }

   for (i = 0; i < n * n; i++)
      h[i] = u[i % n + i / n];

   /* solve a Hankel system */
   if (cholesky(h, p, &(u[n]), n, 1.0e-6) == -1)
      return -1;

   for (i = 0; i < n; i++) {
      p[i] *= -1.0;
   }
   p[n] = 1.0;
   for (i = 0; i < (n + 1) / 2; i++) {
      tmp = p[i];
      p[i] = p[n - i];
      p[n - i] = tmp;
   }

   /* solve roots of a polynomial equation */
   root_pol(p, n, z, 1, 1e-12, 1000);

   for (i = 1; i <= n; i++)
      x[i - 1] = z[i].re;
   quicksort(x, 0, n - 1);

   /* save CSM frequencies */
   for (i = 0; i < n; i++) {
      if (fabs(x[i]) > 1.0)
         return -1;

      csm[i] = acos(x[i]);
   }

   /* solve a Van der Monde system */
   if (vander(x, q, &(u[0]), n, 0.0) == -1)
      return -1;

   /* save CSM intensities */
   for (i = 0; i < n; i++)
      csm[i + n] = q[i];

   return 0;
}

double agcep(double x, double *c, const int m, const int stage,
             const double lambda, const double step, const double tau,
             const double eps)
{
   int i;
   static double *eg = NULL, *ep, *d, gg = 1.0, ee = 1.0, tx;
   static int size;
   double mu, ll;

   if (eg == NULL) {
      eg = dgetmem(2 * (m + 1) + m * stage);
      ep = eg + m + 1;
      d = ep + m + 1;
      size = m * stage;
   }
   if (m * stage > size) {
      free(eg);
      eg = dgetmem(2 * (m + 1) + m * stage);
      ep = eg + m + 1;
      d = ep + m + 1;
      size = m * stage;
   }

   ll = 1.0 - lambda;

   eg[m] = d[stage * m - 1];
   x = iglsadf1(x, c, m, stage, d);

   movem(d + (stage - 1) * m, eg, sizeof(*d), m);

   gg = lambda * gg + ll * eg[0] * eg[0];
   gg = (gg < eps) ? eps : gg;
   mu = step / (double) m / gg;
   tx = 2 * (1.0 - tau) * x;

   for (i = 1; i <= m; i++) {
      ep[i] = tau * ep[i] - tx * eg[i];
      c[i] -= mu * ep[i];
   }

   ee = lambda * ee + ll * x * x;
   c[0] = sqrt(ee);

   return (x);
}

double amcep(double x, double *b, const int m, const double a,
             const double lambda, const double step, const double tau,
             const int pd, const double eps)
{
   int i, memory_size;
   static double *bb = NULL, *d, *ep, *e, xx, gg = 1.0;
   static int size;
   double mu, tx;

   memory_size = 3 * (m + 1) + 3 * (pd + 1) + pd * (m + 2);
   if (bb == NULL) {
      bb = dgetmem(memory_size);
      e = bb + m + 1;
      ep = e + m + 1;
      d = ep + m + 1;
      size = memory_size;
   }
   if (memory_size > size) {
      free(bb);
      bb = dgetmem(memory_size);
      e = bb + m + 1;
      ep = e + m + 1;
      d = ep + m + 1;
      size = memory_size;
   }

   for (i = 1; i <= m; i++)
      bb[i] = -b[i];

   x = mlsadf(x, bb, m, a, pd, d);
   phidf(xx, m, a, e);
   xx = x;

   gg = gg * lambda + (1.0 - lambda) * x * x;
   gg = (gg < eps) ? eps : gg;
   b[0] = 0.5 * log(gg);

   mu = step / (double) m / gg;
   tx = 2 * (1.0 - tau) * x;

   for (i = 1; i <= m; i++) {
      ep[i] = tau * ep[i] - tx * e[i];
      b[i] -= mu * ep[i];
   }

   return (x);
}

void phidf(const double x, const int m, double a, double *d)
{
   int i;

   d[0] = a * d[0] + (1.0 - a * a) * x;
   for (i = 1; i < m; i++)
      d[i] += a * (d[i + 1] - d[i - 1]);

   for (i = m; i >= 1; i--)
      d[i] = d[i - 1];

   return;
}

double average(double *x, const int n)
{
   int i;
   double sum = 0.0;

   for (i = 0; i < n; i++)
      sum += x[i];

   return (sum / n);
}

void vaverage(double *x, const int l, const int num, double *ave)
{
   int i, j;

   fillz(ave, sizeof(*ave), l);
   for (i = 0; i < num; i++)
      for (j = 0; j < l; j++)
         ave[j] += *x++;

   for (j = 0; j < l; j++)
      ave[j] /= (double) num;

   return;
}

void b2mc(double *b, double *mc, int m, const double a)
{
   double d, o;

   d = mc[m] = b[m];
   for (m--; m >= 0; m--) {
      o = b[m] + a * d;
      d = b[m];
      mc[m] = o;
   }

   return;
}

void c2acr(double *c, const int m1, double *r, const int m2, const int flng)
{
   int i;
   static double *x = NULL, *y;
   static int size;

   if (x == NULL) {
      x = dgetmem(flng + flng);
      y = x + flng;
      size = flng;
   }
   if (flng > size) {
      free(x);
      x = dgetmem(flng + flng);
      y = x + flng;
      size = flng;
   }

   movem(c, x, sizeof(*c), m1 + 1);
   fillz(&x[m1 + 1], sizeof(*x), flng - m1 - 1);

   fftr(x, y, flng);

   for (i = 0; i < flng; i++)
      x[i] = exp(2.0 * x[i]);

   fftr(x, y, flng);

   for (i = 0; i <= m2; i++)
      r[i] = x[i] / flng;

   return;
}

void c2ir(double *c, const int nc, double *h, const int leng)
{
   int n, k, upl;
   double d;

   h[0] = exp(c[0]);
   for (n = 1; n < leng; n++) {
      d = 0;
      upl = (n >= nc) ? nc - 1 : n;
      for (k = 1; k <= upl; k++)
         d += k * c[k] * h[n - k];
      h[n] = d / n;
   }

   return;
}

void c2ndps(double *c, const int m, double *n, const int l)
{
   int i;
   double *tmp;

   fillz(n, sizeof(*n), l);
   tmp = dgetmem(l);

   // generate mc(m)
   for (i = 1; i < m + 1; i++) {
      n[i] = c[i] * i / 2.0;
      n[l - i] = n[i];
   }
   if (m == l / 2)
      n[m] = n[m] * 2.0;

   fftr(n, tmp, l);

   free(tmp);
}

void ic2ir(double *h, const int leng, double *c, const int nc)
{
   int n, k, upl;
   double d;

   c[0] = log(h[0]);
   for (n = 1; n < nc; n++) {
      d = (n >= leng) ? 0 : n * h[n];
      upl = (n > leng) ? n - leng + 1 : 1;
      for (k = upl; k < n; k++)
         d -= k * c[k] * h[n - k];
      c[n] = d / (n * h[0]);
   }

   return;
}

void c2sp(double *c, const int m, double *x, double *y, const int l)
{
   int m1;

   m1 = m + 1;

   movem(c, x, sizeof(*c), m1);
   fillz(x + m1, sizeof(*x), l - m1);

   fftr(x, y, l);
}

void clip(double *x, const int l, const double min, const double max, double *y)
{
   int i;

   for (i = 0; i < l; i++)
      y[i] = (x[i] < min) ? min : ((x[i] > max) ? max : x[i]);
}

void csm2acr(double *csm, double *r, const int m)
{
   int i, l, n;
   double *frequencies, *intensities;

   n = 2 * m - 1;

   frequencies = csm;
   intensities = csm + m;

   for (l = 0; l <= n; l++) {
      r[l] = 0.0;
      for (i = 0; i < m; i++) {
         r[l] += intensities[i] * cos(l * frequencies[i]);
      }
   }

   return;
}

/* workspace */
static int dct_table_size = 0;
static double *dct_workspace = NULL;
static double *pLocalReal = NULL;
static double *pLocalImag = NULL;
static double *pWeightReal = NULL;
static double *pWeightImag = NULL;

static int dct_table_size_fft = 0;
static double *dct_workspace2 = NULL;
static double *pLocalReal2 = NULL;
static double *pLocalImag2 = NULL;
static double *pWeightReal2 = NULL;
static double *pWeightImag2 = NULL;


int dft(double *pReal, double *pImag, const int nDFTLength)
{
   int k, n;
   double *pTempReal, *pTempImag, TempReal, TempImag;

   pTempReal = dgetmem(nDFTLength);
   pTempImag = dgetmem(nDFTLength);

   memcpy(pTempReal, pReal, sizeof(double) * nDFTLength);
   memcpy(pTempImag, pImag, sizeof(double) * nDFTLength);

   for (k = 0; k < nDFTLength; k++) {
      TempReal = 0;
      TempImag = 0;
      for (n = 0; n < nDFTLength; n++) {
         TempReal += pTempReal[n] * cos(2.0 * PI * n * k / (double) nDFTLength)
             + pTempImag[n] * sin(2.0 * PI * n * k / (double) nDFTLength);
         TempImag += -pTempReal[n] * sin(2.0 * PI * n * k / (double) nDFTLength)
             + pTempImag[n] * cos(2.0 * PI * n * k / (double) nDFTLength);
      }
      pReal[k] = TempReal;
      pImag[k] = TempImag;
   }
   free(pTempReal);
   free(pTempImag);

   return (0);
}

int dct_create_table_fft(const int nSize)
{
   register int k;

   if (nSize == dct_table_size_fft) {
      /* no needs to resize workspace. */
      return (0);
   } else {
      /* release resources to resize workspace. */
      if (dct_workspace2 != NULL) {
         free(dct_workspace2);
         dct_workspace2 = NULL;
      }
      pLocalReal2 = NULL;
      pLocalImag2 = NULL;
      pWeightReal2 = NULL;
      pWeightImag2 = NULL;
   }

   /* getting resources. */
   if (nSize <= 0) {
      dct_table_size_fft = 0;
      return (0);
   } else {
      dct_table_size_fft = nSize;
      dct_workspace2 = dgetmem(dct_table_size_fft * 6);
      pWeightReal2 = dct_workspace2;
      pWeightImag2 = dct_workspace2 + dct_table_size_fft;
      pLocalReal2 = dct_workspace2 + (2 * dct_table_size_fft);
      pLocalImag2 = dct_workspace2 + (4 * dct_table_size_fft);

      for (k = 0; k < dct_table_size_fft; k++) {
         pWeightReal2[k] =
             cos(k * PI / (2.0 * dct_table_size_fft)) /
             sqrt(2.0 * dct_table_size_fft);
         pWeightImag2[k] =
             -sin(k * PI / (2.0 * dct_table_size_fft)) /
             sqrt(2.0 * dct_table_size_fft);
      }
      pWeightReal2[0] /= sqrt(2.0);
      pWeightImag2[0] /= sqrt(2.0);
   }

   return (0);
}

int dct_create_table(const int nSize)
{
   register int k;

   if (nSize == dct_table_size) {
      /* no needs to resize workspace. */
      return (0);
   } else {
      /* release resources to resize workspace. */
      if (dct_workspace != NULL) {
         free(dct_workspace);
         dct_workspace = NULL;
      }
      pLocalReal = NULL;
      pLocalImag = NULL;
      pWeightReal = NULL;
      pWeightImag = NULL;
   }

   /* getting resources. */
   if (nSize <= 0) {
      dct_table_size = 0;
      return (0);
   } else {
      dct_table_size = nSize;
      dct_workspace = dgetmem(dct_table_size * 6);
      pWeightReal = dct_workspace;
      pWeightImag = dct_workspace + dct_table_size;
      pLocalReal = dct_workspace + (2 * dct_table_size);
      pLocalImag = dct_workspace + (4 * dct_table_size);

      for (k = 0; k < dct_table_size; k++) {
         pWeightReal[k] =
             cos(k * PI / (2.0 * dct_table_size)) / sqrt(2.0 * dct_table_size);
         pWeightImag[k] =
             -sin(k * PI / (2.0 * dct_table_size)) / sqrt(2.0 * dct_table_size);
      }
      pWeightReal[0] /= sqrt(2.0);
      pWeightImag[0] /= sqrt(2.0);
   }

   return (0);
}

int dct_based_on_fft(double *pReal, double *pImag, const double *pInReal,
                     const double *pInImag)
{
   register int n, k;


   for (n = 0; n < dct_table_size_fft; n++) {
      pLocalReal2[n] = pInReal[n];
      pLocalImag2[n] = pInImag[n];
      pLocalReal2[dct_table_size_fft + n] = pInReal[dct_table_size_fft - 1 - n];
      pLocalImag2[dct_table_size_fft + n] = pInImag[dct_table_size_fft - 1 - n];
   }


   fft(pLocalReal2, pLocalImag2, dct_table_size_fft * 2);       /* double input */


   for (k = 0; k < dct_table_size_fft; k++) {
      pReal[k] =
          pLocalReal2[k] * pWeightReal2[k] - pLocalImag2[k] * pWeightImag2[k];
      pImag[k] =
          pLocalReal2[k] * pWeightImag2[k] + pLocalImag2[k] * pWeightReal2[k];
   }

   return (0);
}

int dct_based_on_dft(double *pReal, double *pImag, const double *pInReal,
                     const double *pInImag)
{
   register int n, k;

   for (n = 0; n < dct_table_size; n++) {
      pLocalReal[n] = pInReal[n];
      pLocalImag[n] = pInImag[n];
      pLocalReal[dct_table_size + n] = pInReal[dct_table_size - 1 - n];
      pLocalImag[dct_table_size + n] = pInImag[dct_table_size - 1 - n];
   }

   dft(pLocalReal, pLocalImag, dct_table_size * 2);

   for (k = 0; k < dct_table_size; k++) {
      pReal[k] =
          pLocalReal[k] * pWeightReal[k] - pLocalImag[k] * pWeightImag[k];
      pImag[k] =
          pLocalReal[k] * pWeightImag[k] + pLocalImag[k] * pWeightReal[k];
   }

   return (0);
}

void dct(double *in, double *out, const int size, const int m,
         const Boolean dftmode, const Boolean compmode)
{
   int k, i, j, iter;
   double *pReal, *pImag;
   double *x, *y;
   double *x2, *y2;

   x = dgetmem(2 * size);
   y = x + size;
   pReal = dgetmem(2 * size);
   pImag = pReal + size;
   x2 = dgetmem(2 * size);
   y2 = x2 + size;

   for (k = 0; k < size; k++) {
      x[k] = in[k];
      y[k] = in[k + size];
      x2[k] = x[k];
      y2[k] = y[k];
   }

   iter = 0;
   i = size;
   while ((i /= 2) != 0) {
      iter++;
   }
   j = 1;
   for (i = 1; i <= iter; i++) {
      j *= 2;
   }
   if (size != j || dftmode) {
      dct_create_table(size);
      dct_based_on_dft(pReal, pImag, x2, y2);
   } else {
      dct_create_table_fft(size);
      dct_based_on_fft(pReal, pImag, x2, y2);
   }

   for (k = 0; k < m; k++) {
      out[k] = pReal[k];
      if (compmode == TR) {
         out[k + size] = pImag[k];
      }
   }

   free(x);
   free(x2);
   free(pReal);
}

double df2(const double x, const double sf, const double f0p, const double wbp,
           const double f0z, const double wbz, const int fp, const int fz,
           double *buf, int *bufp)
{
   double a[3], b[3];
   double p, e;

   p = 4 * atan(1.0) / sf;
   e = exp(-p * wbz);

   a[0] = 1.0;
   if (fz) {
      a[1] = -2 * e * cos(2 * p * f0z);
      a[2] = e * e;
   } else {
      a[1] = 0;
      a[2] = 0;
   }

   e = exp(-p * wbp);
   b[0] = 1.0;
   if (fp) {
      b[1] = -2 * e * cos(2 * p * f0p);
      b[2] = e * e;
   } else {
      b[1] = 0;
      b[2] = 0;
   }

   return (dfs(x, b, 2, a, 2, buf, bufp));
}

double dfs(double x, double *a, int m, double *b, int n, double *buf, int *bufp)
{
   double y = 0.0;
   int i, p;
   int max;

   n++;
   m++;

   (m < n) ? (max = n) : (max = m);

   x = x * a[0];
   for (i = 1; i < m; i++) {
      if ((p = *bufp + i) >= max)
         p -= max;
      x -= buf[p] * a[i];
   }
   buf[*bufp] = x;
   for (i = 0; i < n; i++) {
      if ((p = *bufp + i) >= max)
         p -= max;
      y += buf[p] * b[i];
   }

   if (--*bufp < 0)
      *bufp += max;

   return (y);
}

double *_sintbl = 0;
int maxfftsize = 0;

static int checkm(const int m)
{
   int k;

   for (k = 4; k <= m; k <<= 1) {
      if (k == m)
         return (0);
   }
   fprintf(stderr, "fft : m must be a integer of power of 2!\n");

   return (-1);
}

int fft(double *x, double *y, const int m)
{
   int j, lmx, li;
   double *xp, *yp;
   double *sinp, *cosp;
   int lf, lix, tblsize;
   int mv2, mm1;
   double t1, t2;
   double arg;
   int checkm(const int);

   /**************
   * RADIX-2 FFT *
   **************/

   if (checkm(m))
      return (-1);

   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }

   lf = maxfftsize / m;
   lmx = m;

   for (;;) {
      lix = lmx;
      lmx /= 2;
      if (lmx <= 1)
         break;
      sinp = _sintbl;
      cosp = _sintbl + maxfftsize / 4;
      for (j = 0; j < lmx; j++) {
         xp = &x[j];
         yp = &y[j];
         for (li = lix; li <= m; li += lix) {
            t1 = *(xp) - *(xp + lmx);
            t2 = *(yp) - *(yp + lmx);
            *(xp) += *(xp + lmx);
            *(yp) += *(yp + lmx);
            *(xp + lmx) = *cosp * t1 + *sinp * t2;
            *(yp + lmx) = *cosp * t2 - *sinp * t1;
            xp += lix;
            yp += lix;
         }
         sinp += lf;
         cosp += lf;
      }
      lf += lf;
   }

   xp = x;
   yp = y;
   for (li = m / 2; li--; xp += 2, yp += 2) {
      t1 = *(xp) - *(xp + 1);
      t2 = *(yp) - *(yp + 1);
      *(xp) += *(xp + 1);
      *(yp) += *(yp + 1);
      *(xp + 1) = t1;
      *(yp + 1) = t2;
   }

   /***************
   * bit reversal *
   ***************/
   j = 0;
   xp = x;
   yp = y;
   mv2 = m / 2;
   mm1 = m - 1;
   for (lmx = 0; lmx < mm1; lmx++) {
      if ((li = lmx - j) < 0) {
         t1 = *(xp);
         t2 = *(yp);
         *(xp) = *(xp + li);
         *(yp) = *(yp + li);
         *(xp + li) = t1;
         *(yp + li) = t2;
      }
      li = mv2;
      while (li <= j) {
         j -= li;
         li /= 2;
      }
      j += li;
      xp = x + j;
      yp = y + j;
   }

   return (0);
}

int fft2(double x[], double y[], const int n)
{
   double *xq, *yq;
   static double *xb = NULL, *yb;
   double *xp, *yp;
   int i, j;
   static int size_f;

   if (xb == NULL) {
      size_f = 2 * n;
      xb = dgetmem(size_f);
      yb = xb + n;
   }
   if (2 * n > size_f) {
      free(xb);
      size_f = 2 * n;
      xb = dgetmem(size_f);
      yb = xb + n;
   }

   for (i = 0; i < n; i++) {
      xp = xb;
      xq = x + i;
      yp = yb;
      yq = y + i;
      for (j = n; --j >= 0; xq += n, yq += n) {
         *xp++ = *xq;
         *yp++ = *yq;
      }

      if (fft(xb, yb, n) < 0)
         return (-1);

      xp = xb;
      xq = x + i;
      yp = yb;
      yq = y + i;
      for (j = n; --j >= 0; xq += n, yq += n) {
         *xq = *xp++;
         *yq = *yp++;
      }
   }

   for (i = n, xp = x, yp = y; --i >= 0; xp += n, yp += n) {
      if (fft(xp, yp, n) < 0)
         return (-1);
   }

   return 0;
}

void fftcep(double *sp, const int flng, double *c, const int m, int itr,
            double ac)
{
   double temp;
   static double *x = NULL, *y;
   static int size;
   int k;

   if (x == NULL) {
      x = dgetmem(flng + flng);
      y = x + flng;
      size = flng;
   }
   if (flng > size) {
      free(x);
      x = dgetmem(flng + flng);
      y = x + flng;
      size = flng;
   }

   movem(sp, x, sizeof(*sp), flng);

   fftr(x, y, flng);
   for (k = 0; k < flng; k++)
      x[k] /= flng;
   for (k = 0; k <= m; k++) {
      c[k] = x[k];
      x[k] = 0;
   }

   ac += 1.0;
   while (--itr > 0) {
      for (k = 1; k <= m; k++)
         x[flng - k] = x[k];

      fftr(x, y, flng);

      for (k = 0; k < flng; k++)
         if (x[k] < 0.0)
            x[k] = 0.0;
         else
            x[k] /= flng;

      fftr(x, y, flng);

      for (k = 0; k <= m; k++) {
         temp = x[k] * ac;
         c[k] += temp;
         x[k] -= temp;
      }
   }
   c[0] *= 0.5;

   if (m == flng / 2)
      c[m] *= 0.5;
}

int fftr(double *x, double *y, const int m)
{
   int i, j;
   double *xp, *yp, *xq;
   double *yq;
   int mv2, n, tblsize;
   double xt, yt, *sinp, *cosp;
   double arg;

   mv2 = m / 2;

   /* separate even and odd  */
   xq = xp = x;
   yp = y;
   for (i = mv2; --i >= 0;) {
      *xp++ = *xq++;
      *yp++ = *xq++;
   }

   if (fft(x, y, mv2) == -1)    /* m / 2 point fft */
      return (-1);


   /***********************
   * SIN table generation *
   ***********************/

   if ((_sintbl == 0) || (maxfftsize < m)) {
      tblsize = m - m / 4 + 1;
      arg = PI / m * 2;
      if (_sintbl != 0)
         free(_sintbl);
      _sintbl = sinp = dgetmem(tblsize);
      *sinp++ = 0;
      for (j = 1; j < tblsize; j++)
         *sinp++ = sin(arg * (double) j);
      _sintbl[m / 2] = 0;
      maxfftsize = m;
   }

   n = maxfftsize / m;
   sinp = _sintbl;
   cosp = _sintbl + maxfftsize / 4;

   xp = x;
   yp = y;
   xq = xp + m;
   yq = yp + m;
   *(xp + mv2) = *xp - *yp;
   *xp = *xp + *yp;
   *(yp + mv2) = *yp = 0;

   for (i = mv2, j = mv2 - 2; --i; j -= 2) {
      ++xp;
      ++yp;
      sinp += n;
      cosp += n;
      yt = *yp + *(yp + j);
      xt = *xp - *(xp + j);
      *(--xq) = (*xp + *(xp + j) + *cosp * yt - *sinp * xt) * 0.5;
      *(--yq) = (*(yp + j) - *yp + *sinp * yt + *cosp * xt) * 0.5;
   }

   xp = x + 1;
   yp = y + 1;
   xq = x + m;
   yq = y + m;

   for (i = mv2; --i;) {
      *xp++ = *(--xq);
      *yp++ = -(*(--yq));
   }

   return (0);
}

int fftr2(double x[], double y[], const int n)
{
   double *xq, *yq;
   static double *xb = NULL, *yb;
   double *xp, *yp;
   int i, j;
   static int size_f;

   if (xb == NULL) {
      size_f = 2 * n;
      xb = dgetmem(size_f);
      yb = xb + n;
   }
   if (2 * n > size_f) {
      free(xb);
      size_f = 2 * n;
      xb = dgetmem(size_f);
      yb = xb + n;
   }

   for (i = 0; i < n; i++) {
      xp = xb;
      xq = x + i;
      for (j = n; --j >= 0; xq += n) {
         *xp++ = *xq;
      }

      if (fftr(xb, yb, n) < 0)
         return (-1);

      xp = xb;
      xq = x + i;
      yp = yb;
      yq = y + i;
      for (j = n; --j >= 0; xq += n, yq += n) {
         *xq = *xp++;
         *yq = *yp++;
      }
   }

   for (i = n, xp = x, yp = y; --i >= 0; xp += n, yp += n) {
      if (fft(xp, yp, n) < 0)
         return (-1);
   }

   return (0);
}

void freqt(double *c1, const int m1, double *c2, const int m2, const double a)
{
   int i, j;
   double b;
   static double *d = NULL, *g;
   static int size;

   if (d == NULL) {
      size = m2;
      d = dgetmem(size + size + 2);
      g = d + size + 1;
   }

   if (m2 > size) {
      free(d);
      size = m2;
      d = dgetmem(size + size + 2);
      g = d + size + 1;
   }

   b = 1 - a * a;
   fillz(g, sizeof(*g), m2 + 1);

   for (i = -m1; i <= 0; i++) {
      if (0 <= m2)
         g[0] = c1[-i] + a * (d[0] = g[0]);
      if (1 <= m2)
         g[1] = b * d[0] + a * (d[1] = g[1]);
      for (j = 2; j <= m2; j++)
         g[j] = d[j - 1] + a * ((d[j] = g[j]) - g[j - 1]);
   }

   movem(g, c2, sizeof(*g), m2 + 1);

   return;
}

void gc2gc(double *c1, const int m1, const double g1, double *c2, const int m2,
           const double g2)
{
   int i, min, k, mk;
   double ss1, ss2, cc;
   static double *ca = NULL;
   static int size;

   if (ca == NULL) {
      ca = dgetmem(m1 + 1);
      size = m1;
   }
   if (m1 > size) {
      free(ca);
      ca = dgetmem(m1 + 1);
      size = m1;
   }

   movem(c1, ca, sizeof(*c1), m1 + 1);

   c2[0] = ca[0];
   for (i = 1; i <= m2; i++) {
      ss1 = ss2 = 0.0;
      min = (m1 < i) ? m1 : i - 1;
      for (k = 1; k <= min; k++) {
         mk = i - k;
         cc = ca[k] * c2[mk];
         ss2 += k * cc;
         ss1 += mk * cc;
      }

      if (i <= m1)
         c2[i] = ca[i] + (g2 * ss2 - g1 * ss1) / i;
      else
         c2[i] = (g2 * ss2 - g1 * ss1) / i;
   }

   return;
}

int gcep(double *xw, const int flng, double *gc, const int m, const double g,
         const int itr1, const int itr2, const double d, const int etype,
         const double e, const double f, const int itype)
{
   int i, j, flag = 0;
   double t, s, eps = 0.0, min, max, dd = 0.0;
   static double *x = NULL, *y, *cr, *ci, *rr, *hr, *hi, *er, *ei;
   static int size;

   if (etype == 1 && e < 0.0) {
      fprintf(stderr, "gcep : value of e must be e>=0!\n");
      exit(1);
   }

   if (etype == 2 && e >= 0.0) {
      fprintf(stderr, "gcep : value of E must be E<0!\n");
      exit(1);
   }

   if (etype == 1) {
      eps = e;
   }

   if (x == NULL) {
      x = dgetmem(9 * flng);
      size = flng;

      y = x + flng;
      cr = y + flng;
      ci = cr + flng;
      rr = ci + flng;
      hr = rr + flng;
      hi = hr + flng;
      er = hi + flng;
      ei = er + flng;
   }

   if (flng > size) {
      free(x);
      x = dgetmem(9 * flng);
      size = flng;

      y = x + flng;
      cr = y + flng;
      ci = cr + flng;
      rr = ci + flng;
      hr = rr + flng;
      hi = hr + flng;
      er = hi + flng;
      ei = er + flng;
   }

   movem(xw, x, sizeof(*x), flng);

   switch (itype) {
   case 0:                     /* windowed data sequence */
      fftr(x, y, flng);
      for (i = 0; i < flng; i++) {
         x[i] = x[i] * x[i] + y[i] * y[i] + eps;        /*  periodogram  */
      }
      break;
   case 1:                     /* dB */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp((x[i] / 20.0) * log(10.0)); /* dB -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 2:                     /* log */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp(x[i]);      /* log -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 3:                     /* amplitude */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 4:                     /* periodogram */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] + eps;
      }
      break;
   default:
      fprintf(stderr, "gcep : Input type %d is not supported!\n", itype);
      exit(1);
   }
   if (itype > 0) {
      for (i = 1; i < flng / 2; i++)
         x[flng - i] = x[i];
   }

   if (etype == 2 && e < 0.0) {
      max = x[0];
      for (i = 1; i < flng; i++) {
         if (max < x[i])
            max = x[i];
      }
      max = sqrt(max);
      min = max * pow(10.0, e / 20.0);  /* floor is 20*log10(min/max) */
      min = min * min;
      for (i = 0; i < flng; i++) {
         if (x[i] < min)
            x[i] = min;
      }
   }

   for (i = 0; i < flng; i++)
      cr[i] = log(x[i]);

   /*  initial value of generalized cepstrum  */
   ifftr(cr, y, flng);          /*  x : IFFT[x]  */
   cr[0] = exp(cr[0] / 2);
   gc2gc(cr, m, 0.0, gc, m, g); /*  gc : generalized cepstrum  */

   /*  Newton-Raphson method  */
   for (j = 1; j <= itr2; j++) {
      fillz(cr, sizeof(*cr), flng);
      movem(&gc[1], &cr[1], sizeof(*cr), m);
      fftr(cr, ci, flng);       /*  cr+jci : FFT[gc]  */

      for (i = 0; i < flng; i++) {
         t = x[i] / agexp(g, cr[i], ci[i]);
         cr[i] = 1 + g * cr[i];
         ci[i] = g * ci[i];
         s = cr[i] * cr[i] + ci[i] * ci[i];
         rr[i] = t / s;
         hr[i] = (cr[i] * cr[i] - ci[i] * ci[i]) * t / (s * s);
         hi[i] = 2 * cr[i] * ci[i] * t / (s * s);
         er[i] = cr[i] * t / s;
         ei[i] = ci[i] * t / s;
      }

      ifftr(rr, y, flng);       /*  rr : r(k)  */
      ifft(hr, hi, flng);       /*  hr : h(k)  */
      ifft(er, ei, flng);       /*  er : e(k)  */
      s = gc[0];                /*  gc[0] : gain  */

      for (i = 1, t = 0.0; i <= m; i++)
         t += er[i] * gc[i];

      t = er[0] + g * t;
      t = sqrt(fabs(t));

      if (j >= itr1) {
         if (fabs((t - dd) / t) < d) {
            flag = 1;
            break;
         }
         dd = t;
      }

      for (i = 2; i <= m + m; i++)
         hr[i] *= 1 + g;

      if (theq(rr, &hr[2], &y[1], &er[1], m, f)) {
         fprintf(stderr, "gcep : Error in theq() at %dth iteration!\n", j);
         exit(1);
      }

      gc[0] = t;

      for (i = 1; i <= m; i++)
         gc[i] += y[i];
   }

   if (flag)
      return (0);
   else
      return (-1);
}

static double gpoledf(double x, double *c, int m, const double g, double *d)
{
   double y = 0.0;

   for (m--; m > 0; m--) {
      y -= c[m + 1] * d[m];
      d[m] = d[m - 1];
   }
   y -= c[1] * d[0];
   y *= g;
   d[0] = (x += y);

   return (x);
}

double glsadf(double x, double *c, const int m, const int n, double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = poledf(x, c, m, &d[m * i]);

   return (x);
}

double glsadf1(double x, double *c, const int m, const int n, double *d)
{
   int i;
   double gamma;

   gamma = -1 / (double) n;

   for (i = 0; i < n; i++)
      x = gpoledf(x, c, m, gamma, &d[m * i]);

   return (x);
}

/* inverse option */

static double gzerodf(double x, double *c, int m, const double g, double *d)
{
   double y = 0.0;

   for (m--; m > 0; m--) {
      y += c[m + 1] * d[m];
      d[m] = d[m - 1];
   }
   y += c[1] * d[0];
   d[0] = x;

   x += y * g;

   return (x);
}

double iglsadf(double x, double *c, const int m, const int n, double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = zerodf1(x, c, m, &d[m * i]);

   return (x);
}

double iglsadf1(double x, double *c, const int m, const int n, double *d)
{
   int i;
   double gamma;

   gamma = -1 / (double) n;

   for (i = 0; i < n; i++)
      x = gzerodf(x, c, m, gamma, &d[m * i]);

   return (x);
}

/* transpose option */

static double glsadfft(double x, double *c, const int m, double *d)
{
   int i;

   x -= d[0];
   d[m] = c[m] * x;
   for (i = m - 1; i >= 1; i--)
      d[i] += c[i] * x;

   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (x);
}

double glsadft(double x, double *c, const int m, const int n, double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = glsadfft(x, c, m, &d[(m + 1) * i]);
   return (x);
}

static double glsadff1t(double x, double *c, const int m, const double g,
                        double *d)
{
   int i;

   x -= d[0] * g;

   d[m] = c[m] * x;
   for (i = m - 1; i >= 1; i--)
      d[i] += c[i] * x;

   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (x);
}

double glsadf1t(double x, double *c, const int m, const int n, double *d)
{
   int i;
   double gamma;

   gamma = -1 / (double) n;

   for (i = 0; i < n; i++)
      x = glsadff1t(x, c, m, gamma, &d[(m + 1) * i]);

   return (x);
}

/* inverse and transpose */

static double iglsadfft(double x, double *c, const int m, double *d)
{
   int i;
   double y;

   y = x + d[0];

   d[m] = c[m] * x;
   for (i = m - 1; i >= 1; i--)
      d[i] += c[i] * x;

   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (y);
}

double iglsadft(double x, double *c, const int m, const int n, double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = iglsadfft(x, c, m, &d[i * (m + 1)]);

   return (x);
}

static double iglsadff1t(double x, double *c, const int m, const double g,
                         double *d)
{
   int i;
   double y;

   y = x + g * d[0];

   d[m] = c[m] * x;
   for (i = m - 1; i >= 1; i--)
      d[i] += c[i] * x;

   for (i = 0; i < m; i++)
      d[i] = c[i + 1];

   return (y);
}

double iglsadf1t(double x, double *c, const int m, const int n, double *d)
{
   int i;
   double g;

   g = -1.0 / (double) n;

   for (i = 0; i < n; i++)
      x = iglsadff1t(x, c, m, g, &d[i * (m + 1)]);

   return (x);
}

int choleski(double **cov, double **S, const int L)
{
   int i, j, k;
   double tmp;

   for (i = 0; i < L; i++) {
      for (j = 0; j < i; j++) {
         tmp = cov[i][j];
         for (k = 0; k < j; k++)
            tmp -= S[i][k] * S[j][k];
         S[i][j] = tmp / S[j][j];
      }
      tmp = cov[i][i];
      for (k = 0; k < i; k++)
         tmp -= S[i][k] * S[i][k];
      if (tmp <= 0) {
         return 0;
      }
      S[i][i] = sqrt(tmp);
   }
   return 1;
}

double cal_ldet(double **var, const int D)
{
   int i, j, l;
   double ldet = 0.0, **tri;

   tri = (double **) malloc(sizeof(double *) * D);
   for (l = 0; l < D; l++)
      tri[l] = dgetmem(D);

   for (i = 0; i < D; i++)
      for (j = 0; j < D; j++)
         tri[i][j] = 0.0;

   if (choleski(var, tri, D)) {
      for (l = 0; l < D; l++)
         ldet += log(tri[l][l]);

      for (l = 0; l < D; l++) {
         free(tri[l]);
      }
      free(tri);

      return (2.0 * ldet);
   } else {
      for (l = 0; l < D; l++) {
         free(tri[l]);
      }
      free(tri);

      return LZERO;
   }
}

double cal_gconst(double *var, const int D)
{
   int d;
   double gconst;

   gconst = D * log(M_2PI);
   for (d = 0; d < D; d++)
      gconst += log(var[d]);

   return (gconst);
}

double cal_gconstf(double **var, const int D)
{
   double gconst, tmp;

   tmp = cal_ldet(var, D);
   if (tmp == LZERO) {
      fprintf(stderr, "WARNING : det is 0!\n");
      return LZERO;
   }
   gconst = D * log(M_2PI);
   gconst += tmp;

   return (gconst);
}

void cal_tri_inv(double **S, double **S_inv, const int L)
{
   int i, j, k;

   for (i = 0; i < L; i++) {
      S_inv[i][i] = 1.0 / S[i][i];
   }
   for (i = 1; i < L; i++)
      for (j = i - 1; j >= 0; j--)
         for (k = j; k < i; k++) {
            S_inv[i][j] = S_inv[i][j] - S[i][k] * S_inv[k][j] / S[i][i];
         }
}

void cal_inv(double **cov, double **inv, const int L)
{
   int i, j, k;
   double **S, **S_inv;

   S = (double **) malloc(sizeof(double *) * L);
   S_inv = (double **) malloc(sizeof(double *) * L);

   for (i = 0; i < L; i++) {
      S[i] = dgetmem(L);
      S_inv[i] = dgetmem(L);
   }

   for (i = 0; i < L; i++) {
      for (j = 0; j < L; j++) {
         S[i][j] = 0.0;
         S_inv[i][j] = 0.0;
         inv[i][j] = 0.0;
      }
   }

   if (choleski(cov, S, L) == 0)
      return;
   cal_tri_inv(S, S_inv, L);

   for (i = 0; i < L; i++)
      for (j = 0; j < L; j++) {
         if (i > j)
            for (k = i; k < L; k++)
               inv[i][j] = inv[i][j] + S_inv[k][i] * S_inv[k][j];
         else
            for (k = j; k < L; k++)
               inv[i][j] = inv[i][j] + S_inv[k][i] * S_inv[k][j];
      }

   for (i = 0; i < L; i++) {
      free(S[i]);
      free(S_inv[i]);
   }
   free(S);
   free(S_inv);
}

void fillz_GMM(GMM * gmm)
{
   int m, l, ll;

   for (m = 0; m < gmm->nmix; m++) {
      gmm->weight[m] = 0.;
      if (gmm->full != TR) {
         for (l = 0; l < gmm->dim; l++) {
            gmm->gauss[m].mean[l] = 0.0;
            gmm->gauss[m].var[l] = 0.0;
         }
      } else {
         for (l = 0; l < gmm->dim; l++) {
            gmm->gauss[m].mean[l] = 0.0;
            for (ll = 0; ll < gmm->dim; ll++) {
               gmm->gauss[m].cov[l][ll] = 0.0;
               gmm->gauss[m].inv[l][ll] = 0.0;
            }
         }
      }
   }
}

void maskCov_GMM(GMM * gmm, const int *dim_list, const int cov_dim,
                 const Boolean block_full, const Boolean block_corr)
{

   int row, col, i, k, l, m, *offset;

   offset = (int *) malloc(sizeof(int) * cov_dim + 1);

   offset[0] = 0;
   for (i = 1; i < cov_dim + 1; i++) {
      offset[i] = offset[i - 1] + dim_list[i - 1];
   }

   for (m = 0; m < gmm->nmix; m++) {
      if (block_full == FA && block_corr == FA) {       /* without -c1 and -c2 */
         for (k = 0; k < gmm->dim; k++) {
            for (l = 0; l < gmm->dim; l++) {
               if (k != l) {
                  gmm->gauss[m].cov[k][l] = 0.0;
               }
            }
         }
      } else if (block_full == FA && block_corr == TR) {        /* with -c1 */
         for (row = 0; row < cov_dim; row++) {
            for (col = 0; col < cov_dim; col++) {
               for (k = offset[row]; k < offset[row] + dim_list[row]; k++) {
                  for (l = offset[col]; l < offset[col] + dim_list[col]; l++) {
                     if (dim_list[row] != dim_list[col]) {
                        gmm->gauss[m].cov[k][l] = 0.0;
                     } else {
                        if (offset[row + 1] - k != offset[col + 1] - l) {
                           gmm->gauss[m].cov[k][l] = 0.0;
                        }
                     }
                  }
               }
            }
         }
      } else if (block_full == TR && block_corr == FA) {        /* with -c2 */
         for (row = 0; row < cov_dim; row++) {
            for (col = 0; col < cov_dim; col++) {
               if (row != col) {
                  for (k = offset[row]; k < offset[row] + dim_list[row]; k++) {
                     for (l = offset[col]; l < offset[col] + dim_list[col]; l++) {
                        gmm->gauss[m].cov[k][l] = 0.0;
                     }
                  }
               }
            }
         }
      } else {                  /* with -c1 and -c2 */
         for (row = 0; row < cov_dim; row++) {
            for (col = 0; col < cov_dim; col++) {
               if (dim_list[row] != dim_list[col]) {
                  for (k = offset[row]; k < offset[row] + dim_list[row]; k++) {
                     for (l = offset[col]; l < offset[col] + dim_list[col]; l++) {
                        gmm->gauss[m].cov[k][l] = 0.0;
                     }
                  }
               }
            }
         }
      }
   }

   free(offset);

}


double log_wgd(const GMM * gmm, const int m, const int l1, const int l2,
               const double *dat)
{
   int l, ll;
   double sum, *diff = NULL, tmp, lwgd;

   sum = gmm->gauss[m].gconst;

   if (gmm->full != TR) {
      for (l = l1; l < l2; l++) {
         tmp = dat[l] - gmm->gauss[m].mean[l];
         sum += (tmp * tmp) / gmm->gauss[m].var[l];
      }
   } else {
      diff = dgetmem(l2);
      for (l = l1; l < l2; l++) {
         diff[l] = dat[l] - gmm->gauss[m].mean[l];
      }
      for (l = l1; l < l2; l++) {
         for (ll = l1, tmp = 0.0; ll < l2; ll++) {
            tmp += diff[ll] * gmm->gauss[m].inv[ll][l];
         }
         sum += tmp * diff[l];
      }
      free(diff);
   }

   lwgd = log(gmm->weight[m]) - 0.5 * sum;

   return (lwgd);
}

double log_add(double logx, double logy)
{
   double swap, diff, minLogExp, z;

   if (logx < logy) {
      swap = logx;
      logx = logy;
      logy = swap;
   }

   diff = logy - logx;
   minLogExp = -log(-LZERO);

   if (diff < minLogExp)
      return ((logx < LSMALL) ? LZERO : logx);
   else {
      z = exp(diff);
      return (logx + log(1.0 + z));
   }
}

double log_outp(const GMM * gmm, const int l1, const int l2, const double *dat)
{
   int m;
   double logwgd, logb;

   for (m = 0, logb = LZERO; m < gmm->nmix; m++) {
      logwgd = log_wgd(gmm, m, l1, l2, dat);
      logb = log_add(logb, logwgd);
   }
   return (logb);
}

int alloc_GMM(GMM * gmm, const int M, const int L, const Boolean full)
{
   int m;
   gmm->nmix = M;
   gmm->dim = L;
   gmm->full = full;
   gmm->weight = dgetmem(M);
   gmm->gauss = (Gauss *) getmem(sizeof(Gauss), M);
   for (m = 0; m < M; m++) {
      gmm->gauss[m].mean = dgetmem(L);

      if (full != TR) {
         gmm->gauss[m].var = dgetmem(L);
      } else {
         gmm->gauss[m].cov = ddgetmem(L, L);
         gmm->gauss[m].inv = ddgetmem(L, L);
      }
   }

   return (0);
}

int load_GMM(GMM * gmm, FILE * fp)
{
   int m, l;

   freadf(gmm->weight, sizeof(*(gmm->weight)), gmm->nmix, fp);
   for (m = 0; m < gmm->nmix; m++) {
      freadf(gmm->gauss[m].mean, sizeof(*(gmm->gauss[m].mean)), gmm->dim, fp);

      if (gmm->full != TR) {
         freadf(gmm->gauss[m].var, sizeof(*(gmm->gauss[m].var)), gmm->dim, fp);
         gmm->gauss[m].gconst = cal_gconst(gmm->gauss[m].var, gmm->dim);
      } else {
         for (l = 0; l < gmm->dim; l++) {
            freadf(gmm->gauss[m].cov[l],
                   sizeof(*(gmm->gauss[m].cov[l])), gmm->dim, fp);
         }
      }
   }

   return (0);
}

int save_GMM(const GMM * gmm, FILE * fp)
{
   int m, i, j;

   fwritef(gmm->weight, sizeof(*(gmm->weight)), gmm->nmix, fp);
   for (m = 0; m < gmm->nmix; m++) {
      if (gmm->full != TR) {
         fwritef(gmm->gauss[m].mean, sizeof(*(gmm->gauss[m].mean)), gmm->dim,
                 fp);
         fwritef(gmm->gauss[m].var, sizeof(*(gmm->gauss[m].var)), gmm->dim, fp);
      } else {
         fwritef(gmm->gauss[m].mean, sizeof(*(gmm->gauss[m].mean)), gmm->dim,
                 fp);
         for (i = 0; i < gmm->dim; i++) {
            for (j = 0; j < i; j++) {
               gmm->gauss[m].cov[j][i] = gmm->gauss[m].cov[i][j];
            }
         }
         for (i = 0; i < gmm->dim; i++) {
            fwritef(gmm->gauss[m].cov[i],
                    sizeof(*(gmm->gauss[m].cov[i])), gmm->dim, fp);
         }
      }
   }

   return (0);
}

int prepareCovInv_GMM(GMM * gmm)
{
   int m;
   for (m = 0; m < gmm->nmix; m++) {
      cal_inv(gmm->gauss[m].cov, gmm->gauss[m].inv, gmm->dim);
   }

   return (0);
}

int prepareGconst_GMM(GMM * gmm)
{
   int m;

   for (m = 0; m < gmm->nmix; m++) {
      if (gmm->full == FA) {
         gmm->gauss[m].gconst = cal_gconst(gmm->gauss[m].var, gmm->dim);
      } else {
         gmm->gauss[m].gconst = cal_gconstf(gmm->gauss[m].cov, gmm->dim);
      }
      if (gmm->gauss[m].gconst == LZERO) {
         return -1;
      }
   }

   return (0);
}

int floorWeight_GMM(GMM * gmm, double floor)
{
   int m;
   double sum_w = 0.0, sum_floor = floor * gmm->nmix;

   for (m = 0; m < gmm->nmix; m++) {
      if (gmm->weight[m] < floor) {
         gmm->weight[m] = floor;
      }
      sum_w += gmm->weight[m];
   }
   if (sum_w != 1.0) {
      for (m = 0; m < gmm->nmix; m++) {
         gmm->weight[m] =
             (1.0 - sum_floor) / (sum_w - sum_floor) * (gmm->weight[m] -
                                                        floor) + floor;
      }
   }

   return (0);
}

int floorVar_GMM(GMM * gmm, double floor)
{
   int m, l;
   if (gmm->full == FA) {
      for (m = 0; m < gmm->nmix; m++) {
         for (l = 0; l < gmm->dim; l++) {
            if (gmm->gauss[m].var[l] < floor) {
               gmm->gauss[m].var[l] = floor;
            }
         }
      }
   } else {
      for (m = 0; m < gmm->nmix; m++) {
         for (l = 0; l < gmm->dim; l++) {
            if (gmm->gauss[m].cov[l][l] < floor) {
               gmm->gauss[m].cov[l][l] = floor;
            }
         }
      }
   }

   return (0);
}

int free_GMM(GMM * gmm)
{
   int m;

   for (m = 0; m < gmm->nmix; m++) {
      free(gmm->gauss[m].mean);

      if (gmm->full != TR) {
         free(gmm->gauss[m].var);
      } else {
         free(gmm->gauss[m].cov[0]);
         free(gmm->gauss[m].inv[0]);
         free(gmm->gauss[m].cov);
         free(gmm->gauss[m].inv);
      }
   }
   free(gmm->gauss);
   free(gmm->weight);
   gmm->nmix = 0;
   gmm->dim = 0;
   gmm->full = FA;
   gmm->weight = NULL;
   gmm->gauss = NULL;

   return (0);
}

void gnorm(double *c1, double *c2, int m, const double g)
{
   double k;

   if (g != 0.0) {
      k = 1.0 + g * c1[0];
      for (; m >= 1; m--)
         c2[m] = c1[m] / k;
      c2[0] = pow(k, 1.0 / g);
   } else {
      movem(&c1[1], &c2[1], sizeof(*c1), m);
      c2[0] = exp(c1[0]);
   }

   return;
}

void grpdelay(double *x, double *gd, const int size, const int is_arma)
{
   static double *y;
   static int fsize;

   double *u, *v;
   int k, size_2;

   if (fsize < size) {
      if (y != NULL)
         free(y);
      fsize = size;
      y = dgetmem(3 * size);
   }
   movem(x, gd, sizeof(*x), size);
   u = y + size;
   v = u + size;

   size_2 = size / 2;

   if (is_arma)
      gd[0] = 1;
   for (k = 0; k < size; ++k)
      u[k] = gd[k] * k;

   fftr(gd, y, size);
   fftr(u, v, size);

   for (k = 0; k <= size_2; k++) {
      gd[k] = (gd[k] * u[k] + y[k] * v[k]) / (gd[k] * gd[k] + y[k] * y[k]);
      if (is_arma)
         gd[k] *= -1;
   }

   return;
}

int histogram(double *x, const int size, const double min, const double max,
              const double step, double *h)
{
   int k, ii, flg = 0;
   int jj;

   k = (int) ((max - min) / step + 1.0);

   fillz(h, sizeof(*h), k);

   for (ii = 0; ii < size; ii++) {
      if ((x[ii] < min) || (x[ii] > max)) {
         flg = 1;
      } else {
         for (jj = 0; jj < k; jj++) {
            if (x[ii] < min + (jj + 1) * step) {
               h[jj] += 1.0;
               break;
            }
         }
      }
   }

   return (flg);
}

int ifft(double *x, double *y, const int m)
{
   int i;

   if (fft(y, x, m) == -1)
      return (-1);

   for (i = m; --i >= 0; ++x, ++y) {
      *x /= m;
      *y /= m;
   }

   return (0);
}

int ifft2(double x[], double y[], const int n)
{
   double *xq, *yq;
   static double *xb = NULL, *yb;
   double *xp, *yp;
   int i, j;
   static int size_f;

   if (xb == NULL) {
      size_f = 2 * n;
      xb = dgetmem(size_f);
      yb = xb + n;
   }
   if (2 * n > size_f) {
      free(xb);
      size_f = 2 * n;
      xb = dgetmem(size_f);
      yb = xb + n;
   }

   for (i = 0; i < n; i++) {
      xp = xb;
      xq = x + i;
      yp = yb;
      yq = y + i;

      for (j = n; --j >= 0; xq += n, yq += n) {
         *xp++ = *xq;
         *yp++ = *yq;
      }

      if (ifft(xb, yb, n) < 0)
         return (-1);

      xp = xb;
      xq = x + i;
      yp = yb;
      yq = y + i;

      for (j = n; --j >= 0; xq += n, yq += n) {
         *xq = *xp++;
         *yq = *yp++;
      }
   }

   for (i = n, xp = x, yp = y; --i >= 0; xp += n, yp += n) {
      if (ifft(xp, yp, n) < 0)
         return (-1);
   }

   return (0);
}

void ignorm(double *c1, double *c2, int m, const double g)
{
   double k;

   k = pow(c1[0], g);
   if (g != 0.0) {
      for (; m >= 1; m--)
         c2[m] = k * c1[m];
      c2[0] = (k - 1.0) / g;
   } else {
      movem(&c1[1], &c2[1], sizeof(*c1), m);
      c2[0] = log(c1[0]);
   }

   return;
}

void imsvq(int *index, double *cb, const int l, int *cbsize, const int stage,
           double *x)
{
   int i, j;
   static double *xx = NULL;
   static int size;

   if (xx == NULL) {
      xx = dgetmem(l);
      size = l;
   }
   if (size > l) {
      free(xx);
      xx = dgetmem(l);
      size = l;
   }

   fillz(x, sizeof(*x), l);

   for (i = 0; i < stage; i++) {
      ivq(index[i], cb, l, xx);

      for (j = 0; j < l; j++)
         x[j] += xx[j];

      cb += cbsize[i] * l;
   }

   return;
}

void ivq(const int index, double *cb, const int l, double *x)
{
   movem((cb + index * l), x, sizeof(*cb), l);

   return;
}

void lbg(double *x, const int l, const int tnum, double *icb, int icbsize,
         double *cb, const int ecbsize, const int iter, const int mintnum,
         const int seed, const int centup, const double delta, const double end)
{
   int i, j, k, it, maxindex, tnum1, tnum2;
   static int *cntcb, *tindex, size, sizex, sizecb;
   unsigned long next = SEED;
   double d0, d1, dl, err, tmp, rand;
   static double *cb1 = NULL;
   double *p, *q, *r;

   if (cb1 == NULL) {
      cb1 = dgetmem(ecbsize * l);
      tindex = (int *) dgetmem(tnum);
      cntcb = (int *) dgetmem(ecbsize);
      size = l;
      sizex = tnum;
      sizecb = ecbsize;
   }
   if (l > size) {
      free(cb1);
      cb1 = dgetmem(ecbsize * l);
      size = l;
   }
   if (tnum > sizex) {
      free(tindex);
      tindex = (int *) dgetmem(tnum);
      sizex = tnum;
   }
   if (sizecb > ecbsize) {
      free(cb1);
      free(cntcb);
      cb1 = dgetmem(ecbsize * l);
      cntcb = (int *) dgetmem(ecbsize);
   }

   movem(icb, cb, sizeof(*icb), icbsize * l);

   if (seed != 1)
      next = srnd((unsigned int) seed);

   for (; icbsize * 2 <= ecbsize;) {
      q = cb;
      r = cb + icbsize * l;
      for (i = 0; i < icbsize; i++) {
         for (j = 0; j < l; j++) {
            dl = delta * nrandom(&next);
            *r = *q - dl;
            r++;
            *q = *q + dl;
            q++;
         }
      }
      icbsize *= 2;

      d0 = MAXVALUE;
      for (it = 1; it <= iter; it++) {
         fillz((double *) cntcb, sizeof(*cntcb), icbsize);
         d1 = 0.0;
         p = x;
         for (i = 0; i < tnum; i++, p += l) {
            tindex[i] = vq(p, cb, l, icbsize);
            cntcb[tindex[i]]++;

            q = cb + tindex[i] * l;
            d1 += edist(p, q, l);
         }


         d1 /= tnum;
         err = abs((d0 - d1) / d1);

         if (err < end)
            break;

         d0 = d1;
         fillz(cb1, sizeof(*cb), icbsize * l);

         p = x;
         for (i = 0; i < tnum; i++) {
            q = cb1 + tindex[i] * l;
            for (j = 0; j < l; j++)
               *q++ += *p++;
         }

         k = maxindex = 0;
         for (i = 0; i < icbsize; i++)
            if (cntcb[i] > k) {
               k = cntcb[i];
               maxindex = i;
            }


         q = cb;
         r = cb1;
         for (i = 0; i < icbsize; i++, r += l, q += l)
            if (cntcb[i] >= mintnum)
               for (j = 0; j < l; j++)
                  q[j] = r[j] / (double) cntcb[i];
            else {
               if (centup == 1) {
                  p = cb + maxindex * l;
                  for (j = 0; j < l; j++) {
                     rand = nrandom(&next);
                     q[j] = p[j] + delta * rand;
                     p[j] = p[j] - delta * rand;
                  }
               } else if (centup == 2) {
                  if (i < icbsize / 2) {
                     p = q + icbsize / 2 * l;
                     tnum1 = cntcb[i];
                     tnum2 = cntcb[i + icbsize / 2];
                     for (j = 0; j < l; j++) {
                        tmp = (tnum2 * q[j] + tnum1 * p[j]) / (tnum1 + tnum2);
                        rand = nrandom(&next);
                        q[j] = tmp + delta * rand;
                        p[j] = tmp - delta * rand;
                     }
                  } else {
                     p = q - icbsize / 2 * l;
                     tnum1 = cntcb[i];
                     tnum2 = cntcb[i - icbsize / 2];
                     for (j = 0; j < l; j++) {
                        tmp = (tnum2 * q[j] + tnum1 * p[j]) / (tnum1 + tnum2);
                        rand = nrandom(&next);
                        q[j] = tmp + delta * rand;
                        p[j] = tmp - delta * rand;
                     }
                  }
               }
            }
      }
      if (icbsize == ecbsize)
         break;
   }

   return;
}

int levdur(double *r, double *a, const int m, double eps)
{
   int l, k, flag = 0;
   double rmd, mue;
   static double *c = NULL;
   static int size;

   if (c == NULL) {
      c = dgetmem(m + 1);
      size = m;
   }

   if (m > size) {
      free(c);
      c = dgetmem(m + 1);
      size = m;
   }

   if (eps < 0.0)
      eps = 0.0;
   rmd = r[0];
#ifdef WIN32
   if ((((rmd < 0.0) ? -rmd : rmd) <= eps) || _isnan(rmd))
      return (-1);
#else
   if ((((rmd < 0.0) ? -rmd : rmd) <= eps) || isnan(rmd))
      return (-1);
#endif
   a[0] = 0.0;

   for (l = 1; l <= m; l++) {
      mue = -r[l];
      for (k = 1; k < l; k++)
         mue -= c[k] * r[l - k];
      mue = mue / rmd;

      for (k = 1; k < l; k++)
         a[k] = c[k] + mue * c[l - k];
      a[l] = mue;

      rmd = (1.0 - mue * mue) * rmd;
#ifdef WIN32
      if ((((rmd < 0.0) ? -rmd : rmd) <= eps) || _isnan(rmd))
         return (-1);
#else
      if ((((rmd < 0.0) ? -rmd : rmd) <= eps) || isnan(rmd))
         return (-1);
#endif
      if (((mue < 0.0) ? -mue : mue) >= 1.0)
         flag = -2;

      for (k = 0; k <= l; k++)
         c[k] = a[k];
   }
   a[0] = sqrt(rmd);

   return (flag);
}

static double pade[] = { 1.0,
   1.0, 0.0,
   1.0, 0.0, 0.0,
   1.0, 0.0, 0.0, 0.0,
   1.0, 0.4999273, 0.1067005, 0.01170221, 0.0005656279,
   1.0, 0.4999391, 0.1107098, 0.01369984, 0.0009564853, 0.00003041721,
   1.0, 0.500000157834843, 0.113600112846183, 0.015133367945131, 0.001258740956606, 0.000062701416552, 0.000001481891776,
   1.0, 0.499802889651314, 0.115274789205577, 0.015997611632083, 0.001452640362652, 0.000087007832645, 0.000003213962732, 0.000000057148619
};

double *ppade;

/****************************************************************

    double lmafir(x, c, d, m, m1, m2)

    double x  : input
    double *c : cepstrum
    int    m  : order of cepstrum
    double *d : delay
    int    m1 : start order
    int    m2 : end order

*****************************************************************/

static double lmafir(double x, double *c, const int m, double *d, const int m1,
                     const int m2)
{
   int i;

   for (i = m - 1; i >= 1; i--)
      d[i] = d[i - 1];
   d[0] = x;
   for (x = 0.0, i = m1; i <= m2; i++)
      x += c[i] * d[i - 1];

   return (x);
}

double lmadf(double x, double *c, const int m, const int pd, double *d)
{
   ppade = &pade[pd * (pd + 1) / 2];

   x = lmadf1(x, c, m, d, pd, 1, 1);    /* D1(z) */
   x = lmadf1(x, c, m, &d[(m + 1) * pd], pd, 2, m);     /* D2(z) */

   return (x);
}

double cascade_lmadf(double x, double *c, const int m, const int pd, double *d,
                     const int block_num, int *block_size)
{
   int i, block_start = 1, block_end = 0;
   ppade = &pade[pd * (pd + 1) / 2];

   for (i = 0; i < block_num; i++) {
      block_end += abs(block_size[i]);
      if (block_size[i] > 0) {
         x = lmadf1(x, c, m, d, pd, block_start, block_end);
      }
      d += (m + 1) * pd;
      block_start += abs(block_size[i]);
   }

   return (x);
}

/****************************************************************

    double lmadf1(x, c, m, d, m1, m2, pd)

    double x  : input
    double *c : cepstrum
    int    m  : order of cepstrum
    double *d : delay
    int    m1 : start order
    int    m2 : end order
    int    pd : order of Pade approximation

*****************************************************************/

double lmadf1(double x, double *c, const int m, double *d, const int pd,
              const int m1, const int m2)
{
   double y, t, *pt;
   int i;

   pt = &d[pd * m];
   t = lmafir(pt[pd - 1], c, m, &d[(pd - 1) * m], m1, m2);
   y = (t *= ppade[pd]);
   x += (1 & pd) ? t : -t;
   for (i = pd - 1; i >= 1; i--) {
      pt[i] = t = lmafir(pt[i - 1], c, m, &d[(i - 1) * m], m1, m2);
      y += (t *= ppade[i]);
      x += (1 & i) ? t : -t;
   }
   y += (pt[0] = x);

   return (y);
}

double lmadf1t(double x, double *b, const int pd, double *d)
{
   double v, out = 0.0, *pt;
   int i;

   pt = &d[pd + 1];

   for (i = pd; i >= 1; i--) {
      d[i] = pt[i - 1];
      pt[i] = d[i] * b[1];
      v = pt[i] * ppade[i];

      x += (1 & i) ? v : -v;
      out += v;
   }

   pt[0] = x;
   out += x;

   return (out);
}

static double lmafirt(double x, double *b, const int m, double *d, const int m1,
                      const int m2)
{
   int i;
   double y = 0.0;

   y = d[1];

   d[m2] = b[m2] * x;
   for (i = m2 - 1; i >= m1; i--)
      d[i] += b[i] * x;
   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (y);
}

double lmadf2t(double x, double *b, const int m, const int pd, double *d,
               const int m1, const int m2)
{
   double v, out = 0.0, *pt;
   int i;

   pt = &d[pd * (m + 1)];

   for (i = pd; i >= 1; i--) {
      pt[i] = lmafirt(pt[i - 1], b, m, &d[(i - 1) * (m + 1)], m1, m2);
      v = pt[i] * ppade[i];

      x += (1 & i) ? v : -v;
      out += v;
   }

   pt[0] = x;
   out += x;

   return (out);
}

double lmadft(double x, double *c, const int m, const int pd, double *d,
              int block_num, int *block_size)
{
   int i, block_start = 2, block_end = 1;
   ppade = &pade[pd * (pd + 1) / 2];

   x = lmadf1t(x, c, pd, d);

   for (i = 1; i < block_num; i++) {
      block_end += abs(block_size[i]);
      if (block_size[i] > 0) {
         x = lmadf2t(x, c, m, pd, &d[2 * (pd + 1)], block_start, block_end);
      }
      d += pd * (m + 1) + (pd + 1);
      block_start += abs(block_size[i]);
   }

   return (x);
}

int lpc(double *x, const int flng, double *a, const int m, const double f)
{
   int flag;
   static double *r = NULL;
   static int size;

   if (r == NULL) {
      r = dgetmem(m + 1);
      size = m;
   }
   if (m > size) {
      free(r);
      r = dgetmem(m + 1);
      size = m;
   }

   acorr(x, flng, r, m);
   flag = levdur(r, a, m, f);

   return (flag);
}

void lpc2c(double *a, int m1, double *c, const int m2)
{
   int i, k, upl;
   double d;

   c[0] = log(a[0]);
   c[1] = -a[1];
   for (k = 2; k <= m2; ++k) {
      upl = (k > m2) ? m2 + 1 : k;

      for (d = 0.0, i = (k > m1) ? k - m1 : 1; i < upl; i++)
         d += i * c[i] * a[k - i];
      c[k] = -d / k;

      if (k <= m1)
         c[k] -= a[k];
   }

   return;
}

static double chebpoly(const double x, double *c, const int mh)
{
   int i;
   double b[3];

   b[1] = b[2] = 0.0;
   for (i = mh; i > 0; i--) {
      b[0] = 2.0 * x * b[1] - b[2] + c[i];
      b[2] = b[1];
      b[1] = b[0];
   }
   b[0] = x * b[1] - b[2] + c[0];

   return (b[0]);
}

int lpc2lsp(double *lpc, double *lsp, const int order, const int numsp,
            const int maxitr, const double eps)
{
   int i;
   double *p1, *p2;
   int mh1, mh2, mh, mm, itr, flag_odd;
   double delta, x0, x1, g0, g1, x, y;
   static double *c1 = NULL, *c2;
   static int size_order;

   delta = 1.0 / (double) numsp;

   flag_odd = 0;
   if (order % 2 == 0)
      mh1 = mh2 = order / 2;
   else {
      mh1 = (order + 1) / 2;
      mh2 = (order - 1) / 2;
      flag_odd = 1;
   }

   if (c1 == NULL) {
      c1 = dgetmem(2 * (mh1 + 1));
      c2 = c1 + (mh1 + 1);
      size_order = order;
   }
   if (order > size_order) {
      free(c1);
      c1 = dgetmem(2 * (mh1 + 1));
      c2 = c1 + (mh1 + 1);
      size_order = order;
   }

   /* calculate symmetric and antisymmetrica polynomials */
   p1 = lpc + 1;
   p2 = lpc + order;
   c1[mh1] = c2[mh2] = 1.0;
   if (flag_odd) {
      c2[mh2 + 1] = 0.0;
      for (i = mh2 - 1; i >= 0; i--) {
         c1[i + 1] = *p1 + *p2;
         c2[i] = *p1++ - *p2-- + c2[i + 2];
      }
      c1[0] = *p1 + *p2;
   } else {
      for (i = mh1 - 1; i >= 0; i--) {
         c1[i] = *p1 + *p2 - c1[i + 1];
         c2[i] = *p1++ - *p2-- + c2[i + 1];
      }
   }
   c1[0] *= 0.5;
   c2[0] *= 0.5;

   /* root search */
   p1 = c1;
   mh = mh1;
   g0 = chebpoly(1.0, p1, mh);

   mm = 0;
   for (x = 1.0 - delta; x > -delta - 1.0; x -= delta) {
      g1 = chebpoly(x, p1, mh);

      if (g0 * g1 <= 0.0) {
         x0 = x + delta;
         x1 = x;

         itr = 0;
         do {
            x = (x0 + x1) / 2.0;
            y = chebpoly(x, p1, mh);

            if (y * g0 < 0.0) {
               x1 = x;
               g1 = y;
            } else {
               x0 = x;
               g0 = y;
            }

            itr++;
         }
         while ((fabs(y) > eps) && (itr < maxitr));

         x = (g1 * x0 - g0 * x1) / (g1 - g0);
         lsp[mm] = acos(x) / PI2;

         mm++;
         if (mm == order)
            return (0);

         if (p1 == c1) {
            p1 = c2;
            mh = mh2;
         } else {
            p1 = c1;
            mh = mh1;
         }

         g1 = chebpoly(x, p1, mh);
      }
      g0 = g1;
   }
   return (-1);
}

int lpc2par(double *a, double *k, const int m)
{
   int i, n, flg = 0;
   double s;
   static double *kk = NULL, *aa;
   static int size;

   if (kk == NULL) {
      kk = dgetmem(m + m + 2);
      aa = kk + m + 1;
      size = m;
   }

   if (m > size) {
      free(kk);
      kk = dgetmem(m + m + 2);
      aa = kk + m + 1;
      size = m;
   }

   movem(a, aa, sizeof(*aa), m + 1);

   kk[0] = aa[0];
   for (n = m; n >= 1; n--) {
      movem(&aa[1], &kk[1], sizeof(*aa), n);

      if (kk[n] >= 1.0 || kk[n] <= -1.0)
         flg = -1;

      s = 1.0 - kk[n] * kk[n];
      for (i = 1; i < n; i++)
         aa[i] = (kk[i] - kk[n] * kk[n - i]) / s;
   }
   movem(kk, k, sizeof(*kk), m + 1);

   return (flg);
}

void lsp2lpc(double *lsp, double *a, const int m)
{
   int i, k, mh1, mh2, flag_odd;
   double xx, xf, xff;
   static double *f = NULL, *p, *q, *a0, *a1, *a2, *b0, *b1, *b2;
   static int size;

   flag_odd = 0;
   if (m % 2 == 0)
      mh1 = mh2 = m / 2;
   else {
      mh1 = (m + 1) / 2;
      mh2 = (m - 1) / 2;
      flag_odd = 1;
   }

   if (f == NULL) {
      f = dgetmem(5 * m + 6);
      p = f + m;
      q = p + mh1;
      a0 = q + mh2;
      a1 = a0 + (mh1 + 1);
      a2 = a1 + (mh1 + 1);
      b0 = a2 + (mh1 + 1);
      b1 = b0 + (mh2 + 1);
      b2 = b1 + (mh2 + 1);
      size = m;
   }
   if (m > size) {
      free(f);
      f = dgetmem(5 * m + 6);
      p = f + m;
      q = p + mh1;
      a0 = q + mh2;
      a1 = a0 + (mh1 + 1);
      a2 = a1 + (mh1 + 1);
      b0 = a2 + (mh1 + 1);
      b1 = b0 + (mh2 + 1);
      b2 = b1 + (mh2 + 1);
      size = m;
   }

   movem(lsp, f, sizeof(*lsp), m);

   fillz(a0, sizeof(*a0), mh1 + 1);
   fillz(b0, sizeof(*b0), mh2 + 1);
   fillz(a1, sizeof(*a1), mh1 + 1);
   fillz(b1, sizeof(*b1), mh2 + 1);
   fillz(a2, sizeof(*a2), mh1 + 1);
   fillz(b2, sizeof(*b2), mh2 + 1);

   /* lsp filter parameters */
   for (i = k = 0; i < mh1; i++, k += 2)
      p[i] = -2.0 * cos(PI2 * f[k]);
   for (i = k = 0; i < mh2; i++, k += 2)
      q[i] = -2.0 * cos(PI2 * f[k + 1]);

   /* impulse response of analysis filter */
   xx = 1.0;
   xf = xff = 0.0;
   for (k = 0; k <= m; k++) {
      if (flag_odd) {
         a0[0] = xx;
         b0[0] = xx - xff;
         xff = xf;
         xf = xx;
      } else {
         a0[0] = xx + xf;
         b0[0] = xx - xf;
         xf = xx;
      }

      for (i = 0; i < mh1; i++) {
         a0[i + 1] = a0[i] + p[i] * a1[i] + a2[i];
         a2[i] = a1[i];
         a1[i] = a0[i];
      }
      for (i = 0; i < mh2; i++) {
         b0[i + 1] = b0[i] + q[i] * b1[i] + b2[i];
         b2[i] = b1[i];
         b1[i] = b0[i];
      }

      if (k != 0)
         a[k - 1] = -0.5 * (a0[mh1] + b0[mh2]);

      xx = 0.0;
   }

   for (i = m - 1; i >= 0; i--)
      a[i + 1] = -a[i];
   a[0] = 1.0;

   return;
}

double log_conv(double x)
{
   double temp;

   temp = log(fabs(x));
   if (temp < LSMALL)
      return LZERO;
   else
      return temp;
}

double log_add2(double x, double y)
{
   double lmin, lmax;

   if (x == y)
      return x + LOG2;

   lmin = (x < y) ? x : y;
   lmax = (x < y) ? y : x;
   if (lmax > lmin + 50)
      return lmax;
   else
      return lmax + log_conv(exp(lmin - lmax) + 1.0);
}

void lsp2sp(double *lsp, const int m, double *x, const int l, const int gain)
{
   int i, p;
   double w, eq1, eq2, ap = 0.0;

   for (p = 0; p < l; p++) {
      eq1 = 0.0;
      eq2 = 0.0;
      w = p * (M_PI / (l - 1));

      if (m % 2 == 0) {
         for (i = 0; i < m / 2; i++) {
            eq1 += 2 * log_conv(cos(w) - cos(lsp[2 * i + gain]));
            eq2 += 2 * log_conv(cos(w) - cos(lsp[2 * i + 1 + gain]));
         }
         eq1 += 2 * log_conv(cos(w / 2));
         eq2 += 2 * log_conv(sin(w / 2));

         ap = m * log(2) + log_add2(eq1, eq2);
      } else {
         for (i = 0; i < (m + 1) / 2; i++)
            eq1 += 2 * log_conv(cos(w) - cos(lsp[2 * i + gain]));
         for (i = 0; i < (m - 1) / 2; i++)
            eq2 += 2 * log_conv(cos(w) - cos(lsp[2 * i + 1 + gain]));
         eq2 += 2 * log_conv(sin(w));

         ap = (m - 1) * log(2) + log_add2(eq1, eq2);
      }

      x[p] = -0.5 * ap;
      if (gain == 1)
         x[p] += lsp[0];
   }
}

int lspcheck(double *lsp, const int ord)
{
   int i;

   for (i = 1; i < ord; i++) {
      if (lsp[i] <= lsp[i - 1])
         return (-1);
   }
   if ((lsp[0] <= 0.0) || (lsp[ord - 1] >= 0.5))
      return (-1);

   return (0);
}

double lspdf_even(double x, double *f, const int m, double *d)
{
   double *d1, *d2, *lsp, x1, x2;
   int i;

   d1 = d + 1;
   d2 = d1 + m;
   lsp = f + 1;
   x1 = x2 = d[0];

   for (i = 0; i < m; i += 2) {
      d1[i] -= 2.0 * x1 * cos(lsp[i]);
      d2[i] -= 2.0 * x2 * cos(lsp[i + 1]);
      d1[i + 1] += x1;
      d2[i + 1] += x2;
      x += d1[i] + d2[i];
      x1 = d1[i + 1];
      x2 = d2[i + 1];
   }

   x -= d2[m - 1] - d1[m - 1];

   for (i = m - 1; i > 0; i--) {
      d1[i] = d1[i - 1];
      d2[i] = d2[i - 1];
   }
   d1[0] = d2[0] = d[0];
   d[0] = -0.5 * x;

   return (x);
}

double lspdf_odd(double x, double *f, const int m, double *d)
{
   int i;
   int mh1;
   double *d1, *d2, *lsp, x1, x2;

   mh1 = (m + 1) / 2;

   d1 = d + 1;
   d2 = d1 + (mh1 + mh1 - 1);
   lsp = f + 1;
   x1 = x2 = d[0];

   for (i = 0; i < m - 1; i += 2) {
      d1[i] -= 2.0 * x1 * cos(lsp[i]);
      d2[i] -= 2.0 * x2 * cos(lsp[i + 1]);
      d1[i + 1] += x1;
      d2[i + 1] += x2;
      x += d1[i] + d2[i];
      x1 = d1[i + 1];
      x2 = d2[i + 1];
   }
   d1[i] -= 2.0 * x1 * cos(lsp[i]);
   x += d1[i] - d2[i];

   for (i = m - 1; i > 0; i--) {
      d1[i] = d1[i - 1];
      d2[i] = d2[i - 1];
   }
   d1[0] = d2[0] = d[0];
   d[0] = -0.5 * x;

   return (x);
}

double ltcdf(double x, double *k, int m, double *d)
{
   x -= k[m] * d[m - 1];
   for (m--; m >= 1; m--) {
      x -= k[m] * d[m - 1];
      d[m] = d[m - 1] + k[m] * x;
   }
   d[0] = x;

   return (x);
}

void mc2b(double *mc, double *b, int m, const double a)
{
   b[m] = mc[m];

   for (m--; m >= 0; m--)
      b[m] = mc[m] - a * b[m + 1];

   return;
}

int mcep(double *xw, const int flng, double *mc, const int m, const double a,
         const int itr1, const int itr2, const double dd, const int etype,
         const double e, const double f, const int itype)
{
   int i, j;
   int flag = 0, f2, m2;
   double t, s, eps = 0.0, min, max;
   static double *x = NULL, *y, *c, *d, *al, *b;
   static int size_x, size_d;

   if (etype == 1 && e < 0.0) {
      fprintf(stderr, "mcep : value of e must be e>=0!\n");
      exit(1);
   }

   if (etype == 2 && e >= 0.0) {
      fprintf(stderr, "mcep : value of E must be E<0!\n");
      exit(1);
   }

   if (etype == 1) {
      eps = e;
   }

   if (x == NULL) {
      x = dgetmem(3 * flng);
      y = x + flng;
      c = y + flng;
      size_x = flng;

      d = dgetmem(3 * m + 3);
      al = d + (m + 1);
      b = al + (m + 1);
      size_d = m;
   }
   if (flng > size_x) {
      free(x);
      x = dgetmem(3 * flng);
      y = x + flng;
      c = y + flng;
      size_x = flng;
   }
   if (m > size_d) {
      free(d);
      d = dgetmem(3 * m + 3);
      al = d + (m + 1);
      b = al + (m + 1);
      size_d = m;
   }

   f2 = flng / 2;
   m2 = m + m;

   movem(xw, x, sizeof(*x), flng);

   switch (itype) {
   case 0:                     /* windowed data sequence */
      fftr(x, y, flng);
      for (i = 0; i < flng; i++) {
         x[i] = x[i] * x[i] + y[i] * y[i] + eps;        /*  periodogram  */
      }
      break;
   case 1:                     /* dB */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp((x[i] / 20.0) * log(10.0)); /* dB -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 2:                     /* log */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp(x[i]);      /* log -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 3:                     /* amplitude */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 4:                     /* periodogram */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] + eps;
      }
      break;
   default:
      fprintf(stderr, "mcep : input type %d is not supported!\n", itype);
      exit(1);
   }
   if (itype > 0) {
      for (i = 1; i < flng / 2; i++)
         x[flng - i] = x[i];
   }

   if (etype == 2 && e < 0.0) {
      max = x[0];
      for (i = 1; i < flng; i++) {
         if (max < x[i])
            max = x[i];
      }
      max = sqrt(max);
      min = max * pow(10.0, e / 20.0);  /* floor is 20*log10(min/max) */
      min = min * min;
      for (i = 0; i < flng; i++) {
         if (x[i] < min)
            x[i] = min;
      }
   }

   for (i = 0; i < flng; i++) {
      if (x[i] <= 0.0) {
         fprintf(stderr,
                 "mcep : periodogram has '0', use '-e' option to floor it!\n");
         exit(1);
      }
      c[i] = log(x[i]);
   }

   /*  1, (-a), (-a)^2, ..., (-a)^M  */
   al[0] = 1.0;
   for (i = 1; i <= m; i++)
      al[i] = -a * al[i - 1];

   /*  initial value of cepstrum  */
   ifftr(c, y, flng);           /*  c : IFFT[x]  */

   c[0] /= 2.0;
   c[f2] /= 2.0;
   freqt(c, f2, mc, m, a);      /*  mc : mel cep.  */
   s = c[0];

   /*  Newton Raphson method  */
   for (j = 1; j <= itr2; j++) {
      fillz(c, sizeof(*c), flng);
      freqt(mc, m, c, f2, -a);  /*  mc : mel cep.  */
      fftr(c, y, flng);         /*  c, y : FFT[mc]  */
      for (i = 0; i < flng; i++)
         c[i] = x[i] / exp(c[i] + c[i]);
      ifftr(c, y, flng);
      frqtr(c, f2, c, m2, a);   /*  c : r(k)  */

      t = c[0];
      if (j >= itr1) {
         if (fabs((t - s) / t) < dd) {
            flag = 1;
            break;
         }
         s = t;
      }

      for (i = 0; i <= m; i++)
         b[i] = c[i] - al[i];
      for (i = 0; i <= m2; i++)
         y[i] = c[i];
      for (i = 0; i <= m2; i += 2)
         y[i] -= c[0];
      for (i = 2; i <= m; i += 2)
         c[i] += c[0];
      c[0] += c[0];

      if (theq(c, y, d, b, m + 1, f)) {
         fprintf(stderr, "mcep : Error in theq() at %dth iteration !\n", j);
         exit(1);
      }

      for (i = 0; i <= m; i++)
         mc[i] += d[i];
   }

   if (flag)
      return (0);
   else
      return (-1);

}

/***************************************************************

    Frequency Transformation for Calculating Coefficients

        void frqtr(c1, m1, c2, m2, a)

        double *c1   : minimum phase sequence
        int m1       : order of minimum phase sequence
        double *c2   : warped sequence
        int m2       : order of warped sequence
        double a     : all-pass constant

***************************************************************/

void frqtr(double *c1, int m1, double *c2, int m2, const double a)
{
   int i, j;
   static double *d = NULL, *g;
   static int size;

   if (d == NULL) {
      size = m2;
      d = dgetmem(size + size + 2);
      g = d + size + 1;
   }

   if (m2 > size) {
      free(d);
      size = m2;
      d = dgetmem(size + size + 2);
      g = d + size + 1;
   }

   fillz(g, sizeof(*g), m2 + 1);

   for (i = -m1; i <= 0; i++) {
      if (0 <= m2) {
         d[0] = g[0];
         g[0] = c1[-i];
      }
      for (j = 1; j <= m2; j++)
         g[j] = d[j - 1] + a * ((d[j] = g[j]) - g[j - 1]);
   }

   movem(g, c2, sizeof(*g), m2 + 1);

   return;
}

double freq_mel(double freq)
{
   return MEL * log(freq / 700.0 + 1.0);
}


double sample_mel(int sample, int num, double fs)
{
   double freq;
   freq = (double) (sample + 1) / (double) (num) * (fs / 2.0);

   return freq_mel(freq);
}

double cal_energy(double *x, const int leng)
{
   int k;
   double energy = 0.0;
   for (k = 0; k < leng; k++)
      energy += x[k] * x[k];

   return ((energy <= 0) ? EZERO : log(energy));
}

void pre_emph(double *x, double *y, const double alpha, const int leng)
{
   int k;
   y[0] = x[0] * (1.0 - alpha);
   for (k = 1; k < leng; k++)
      y[k] = x[k] - x[k - 1] * alpha;
}

void spec(double *x, double *sp, const int leng)
{
   int k, no;
   double *y, *mag;

   no = leng / 2;

   y = dgetmem(leng + no);
   mag = y + leng;

   fftr(x, y, leng);
   for (k = 1; k < no; k++) {
      mag[k] = x[k] * x[k] + y[k] * y[k];
      sp[k] = sqrt(mag[k]);
   }
   free(y);
}

void fbank(double *x, double *fb, const double eps, const double fs,
           const int leng, const int n)
{
   int k, fnum, no, chanNum = 0;
   int *noMel;
   double *w, *countMel;
   double maxMel, kMel;

   no = leng / 2;
   noMel = (int *) getmem((size_t) no, sizeof(int));
   countMel = dgetmem(n + 1 + no);
   w = countMel + n + 1;
   maxMel = freq_mel(fs / 2.0);

   for (k = 0; k <= n; k++)
      countMel[k] = (double) (k + 1) / (double) (n + 1) * maxMel;
   for (k = 1; k < no; k++) {
      kMel = sample_mel(k - 1, no, fs);
      while (countMel[chanNum] < kMel && chanNum <= n)
         chanNum++;
      noMel[k] = chanNum;
   }

   for (k = 1; k < no; k++) {
      chanNum = noMel[k];
      kMel = sample_mel(k - 1, no, fs);
      w[k] = (countMel[chanNum] - kMel) / (countMel[0]);
   }

   for (k = 1; k < no; k++) {
      fnum = noMel[k];
      if (fnum > 0)
         fb[fnum] += x[k] * w[k];
      if (fnum <= n)
         fb[fnum + 1] += (1 - w[k]) * x[k];
   }

   free(noMel);
   free(countMel);

   for (k = 1; k <= n; k++) {
      if (fb[k] < eps)
         fb[k] = eps;
      fb[k] = log(fb[k]);
   }
}



void lifter(double *x, double *y, const int m, const int leng)
{
   int k;
   double theta;
   for (k = 0; k < m; k++) {
      theta = PI * (double) k / (double) leng;
      y[k] = (1.0 + (double) leng / 2.0 * sin(theta)) * x[k];
   }
}

void mfcc(double *in, double *mc, const double sampleFreq, const double alpha,
          const double eps, const int wlng, const int flng, const int m,
          const int n, const int ceplift, const Boolean dftmode,
          const Boolean usehamming)
{
   static double *x = NULL, *px, *wx, *sp, *fb, *dc;
   double energy = 0.0, c0 = 0.0;
   int k;

   if (x == NULL) {
      x = dgetmem(wlng + wlng + flng + flng + 2 * n + 1 + m);
      px = x + wlng;
      wx = px + wlng;
      sp = wx + flng;
      fb = sp + flng;
      dc = fb + 2 * n + 1;
   } else {
      free(x);
      x = dgetmem(wlng + wlng + flng + flng + 2 * n + 1 + m);
      px = x + wlng;
      wx = px + wlng;
      sp = wx + flng;
      fb = sp + flng;
      dc = fb + 2 * n + 1;
   }

   movem(in, x, sizeof(*in), wlng);
   /* calculate energy */
   energy = cal_energy(x, wlng);
   pre_emph(x, px, alpha, wlng);
   /* apply hamming window */
   if (usehamming)
      window(HAMMING, px, wlng, 0);
   for (k = 0; k < wlng; k++)
      wx[k] = px[k];
   spec(wx, sp, flng);
   fillz(fb + 1, 2 * n, sizeof(*fb));
   fbank(sp, fb, eps, sampleFreq, flng, n);
   /* calculate 0'th coefficient */
   for (k = 1; k <= n; k++)
      c0 += fb[k];
   c0 *= sqrt(2.0 / (double) n);
   dct(fb + 1, dc, n, m, dftmode, FA);

   /* liftering */
   if (ceplift > 0)
      lifter(dc, mc, m, ceplift);
   else
      movem(dc, mc, sizeof(*dc), m);

   for (k = 0; k < m - 1; k++)
      mc[k] = mc[k + 1];
   mc[m - 1] = c0;
   mc[m] = energy;

}

void mgc2mgc(double *c1, const int m1, const double a1, const double g1,
             double *c2, const int m2, const double a2, const double g2)
{
   double a;
   static double *ca = NULL;
   static int size_a;

   if (ca == NULL) {
      ca = dgetmem(m1 + 1);
      size_a = m1;
   }
   if (m1 > size_a) {
      free(ca);
      ca = dgetmem(m1 + 1);
      size_a = m1;
   }

   a = (a2 - a1) / (1 - a1 * a2);

   if (a == 0) {
      movem(c1, ca, sizeof(*c1), m1 + 1);
      gnorm(ca, ca, m1, g1);
      gc2gc(ca, m1, g1, c2, m2, g2);
      ignorm(c2, c2, m2, g2);
   } else {
      freqt(c1, m1, c2, m2, a);
      gnorm(c2, c2, m2, g1);
      gc2gc(c2, m2, g1, c2, m2, g2);
      ignorm(c2, c2, m2, g2);
   }

   return;
}

void mgc2sp(double *mgc, const int m, const double a, const double g, double *x,
            double *y, const int flng)
{
   static double *c = NULL;
   static int size;

   if (c == NULL) {
      c = dgetmem(flng / 2 + 1);
      size = flng;
   }
   if (flng > size) {
      free(c);
      c = dgetmem(flng / 2 + 1);
      size = flng;
   }

   mgc2mgc(mgc, m, a, g, c, flng / 2, 0.0, 0.0);
   c2sp(c, flng / 2, x, y, flng);

   return;
}

double mel_conv(double a, double w)
{
   return w + 2.0 * atan(a * sin(w) / (1.0 - a * cos(w)));
}

void mgclsp2sp(double a, double g, double *lsp, const int m, double *x,
               const int l, const int gain)
{
   int i, p;
   double w, eq1, eq2, ap = 0.0;

   for (p = 0; p < l; p++) {
      eq1 = 0.0;
      eq2 = 0.0;
      w = mel_conv(a, p * (M_PI / (l - 1)));

      if (m % 2 == 0) {
         for (i = 0; i < m / 2; i++) {
            eq1 += 2.0 * log_conv(cos(w) - cos(lsp[2 * i + gain]));
            eq2 += 2.0 * log_conv(cos(w) - cos(lsp[2 * i + 1 + gain]));
         }
         eq1 += 2.0 * log_conv(cos(w / 2.0));
         eq2 += 2.0 * log_conv(sin(w / 2.0));

         ap = m * log(2.0) + log_add2(eq1, eq2);
      } else {
         for (i = 0; i < (m + 1) / 2; i++)
            eq1 += 2.0 * log_conv(cos(w) - cos(lsp[2 * i + gain]));
         for (i = 0; i < (m - 1) / 2; i++)
            eq2 += 2.0 * log_conv(cos(w) - cos(lsp[2 * i + 1 + gain]));
         eq2 += 2.0 * log_conv(sin(w));

         ap = (m - 1.0) * log(2.0) + log_add2(eq1, eq2);
      }

      x[p] = -0.5 * ap;
      x[p] *= -(1.0 / g);
      if (gain == 1)
         x[p] += lsp[0];
   }
}

/*  gain(epsilon) calculation  */
static double gain(double *er, double *c, int m, double g)
{
   int i;
   double t;

   if (g != 0.0) {
      for (t = 0.0, i = 1; i <= m; i++)
         t += er[i] * c[i];
      return (er[0] + g * t);
   } else
      return (er[0]);
}

/*  b'(m) to c(m)  */
static void b2c(double *b, int m1, double *c, int m2, double a)
{
   int i, j;
   static double *d = NULL, *g;
   static int size;
   double k;

   if (d == NULL) {
      size = m2;
      d = dgetmem(size + size + 2);
      g = d + size + 1;
   }
   if (m2 > size) {
      free(d);
      size = m2;
      d = dgetmem(size + size + 2);
      g = d + size + 1;
   }

   k = 1 - a * a;

   fillz(g, sizeof(*g), m2 + 1);

   for (i = -m1; i <= 0; i++) {
      d[0] = g[0];
      g[0] = b[-i];

      if (1 <= m2)
         g[1] = k * d[0] + a * (d[1] = g[1]);

      for (j = 2; j <= m2; j++)
         g[j] = d[j - 1] + a * ((d[j] = g[j]) - g[j - 1]);
   }
   movem(g, c, sizeof(*g), m2 + 1);

   return;
}

/*  recursion for p(m)  */
static void ptrans(double *p, int m, double a)
{
   double d, o;

   d = p[m];
   for (m--; m > 0; m--) {
      o = p[m] + a * d;
      d = p[m];
      p[m] = o;
   }
   o = a * d;
   p[m] = (1. - a * a) * p[m] + o + o;

   return;
}

/*  recursion for q(m)  */
static void qtrans(double *q, int m, double a)
{
   int i;
   double d, o;

   m += m;
   i = 1;
   d = q[i];
   for (i++; i <= m; i++) {
      o = q[i] + a * d;
      d = q[i];
      q[i] = o;
   }

   return;
}

int mgcep(double *xw, int flng, double *b, const int m, const double a,
          const double g, const int n, const int itr1, const int itr2,
          const double dd, const int etype, const double e, const double f,
          const int itype)
{
   int i, j, flag = 0;
   static double *x = NULL, *y, *d;
   static int size_x, size_c;
   double ep, epo, eps = 0.0, min, max;

   if (etype == 1 && e < 0.0) {
      fprintf(stderr, "mgcep : value of e must be e>=0!\n");
      exit(1);
   }

   if (etype == 2 && e >= 0.0) {
      fprintf(stderr, "mgcep : value of E must be E<0!\n");
      exit(1);
   }

   if (etype == 1) {
      eps = e;
   }

   if (x == NULL) {
      x = dgetmem(flng + flng);
      y = x + flng;
      size_x = flng;
      d = dgetmem(m + 1);
      size_c = m;
   }
   if (flng > size_x) {
      free(x);
      x = dgetmem(flng + flng);
      y = x + flng;
      size_x = flng;
   }
   if (m > size_c) {
      free(d);
      d = dgetmem(m + 1);
      size_c = m;
   }

   movem(xw, x, sizeof(*x), flng);

   switch (itype) {
   case 0:                     /* windowed data sequence */
      fftr(x, y, flng);
      for (i = 0; i < flng; i++) {
         x[i] = x[i] * x[i] + y[i] * y[i] + eps;        /*  periodogram  */
      }
      break;
   case 1:                     /* dB */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp((x[i] / 20.0) * log(10.0)); /* dB -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 2:                     /* log */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp(x[i]);      /* log -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 3:                     /* amplitude */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 4:                     /* periodogram */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] + eps;
      }
      break;
   default:
      fprintf(stderr, "mgcep : Input type %d is not supported!\n", itype);
      exit(1);
   }
   if (itype > 0) {
      for (i = 1; i < flng / 2; i++)
         x[flng - i] = x[i];
   }

   if (etype == 2 && e < 0.0) {
      max = x[0];
      for (i = 1; i < flng; i++) {
         if (max < x[i])
            max = x[i];
      }
      max = sqrt(max);
      min = max * pow(10.0, e / 20.0);  /* floor is 20*log10(min/max) */
      min = min * min;
      for (i = 0; i < flng; i++) {
         if (x[i] < min)
            x[i] = min;
      }
   }

   /* initial value */
   fillz(b, sizeof(*b), m + 1);
   ep = newton(x, flng, b, m, a, -1.0, n, 0, f);

   if (g != -1.0) {
      if (a != 0.0) {
         ignorm(b, b, m, -1.0); /*  K, b'r(m)    -> br(m)         */
         b2mc(b, b, m, a);      /*  br(m)        -> c~r(m)        */
         gnorm(b, d, m, -1.0);  /*  c~r(m)       -> K~, c~'r(m)   */
      } else
         movem(b, d, sizeof(*b), m + 1);

      gc2gc(d, m, -1.0, b, m, g);       /*  K~, c~'r(m)  -> K~, c~'r'(m)  */

      if (a != 0.0) {
         ignorm(b, b, m, g);    /*  K~, c~'r'(m) -> c~r(m)        */
         mc2b(b, b, m, a);      /*  c~r(m)       -> br(m)         */
         gnorm(b, b, m, g);     /*  br(m)        -> K, b'r'(m)    */
      }
   }

   /*  Newton-Raphson method  */
   if (g != -1.0) {
      for (j = 1; j <= itr2; j++) {
         epo = ep;
         ep = newton(x, flng, b, m, a, g, n, j, f);

         if (j >= itr1)
            if (fabs((epo - ep) / ep) < dd) {
               flag = 1;
               break;
            }
      }
   }

   if (flag)
      return (0);
   else
      return (-1);
}

double newton(double *x, const int flng, double *c, const int m, const double a,
              const double g, const int n, const int j, const double f)
{
   int i, m2;
   double t = 0, s, tr, ti, trr, tii;
   static double *cr = NULL, *ci, *pr, *qr, *qi, *rr, *ri, *b;
   static int size_cr, size_b;

   if (cr == NULL) {
      cr = dgetmem(7 * flng);
      ci = cr + flng;
      pr = ci + flng;
      qr = pr + flng;
      qi = qr + flng;
      rr = qi + flng;
      ri = rr + flng;
      size_cr = flng;

      b = dgetmem(m + 1);
      size_b = m;
   }
   if (flng > size_cr) {
      free(cr);
      cr = dgetmem(7 * flng);
      ci = cr + flng;
      pr = ci + flng;
      qr = pr + flng;
      qi = qr + flng;
      rr = qi + flng;
      ri = rr + flng;
      size_cr = flng;
   }
   if (m > size_b) {
      free(b);
      b = dgetmem(m + 1);
      size_b = m;
   }

   m2 = m + m;

   fillz(cr, sizeof(*cr), flng);
   movem(&c[1], &cr[1], sizeof(*c), m);

   if (a != 0.0)
      b2c(cr, m, cr, n, -a);

   fftr(cr, ci, flng);          /* cr +j ci : FFT[c]  */

   if (g == -1.0)
      movem(x, pr, sizeof(*x), flng);
   else if (g == 0.0)
      for (i = 0; i < flng; i++)
         pr[i] = x[i] / exp(cr[i] + cr[i]);
   else
      for (i = 0; i < flng; i++) {
         tr = 1 + g * cr[i];
         ti = g * ci[i];
         s = (trr = tr * tr) + (tii = ti * ti);
         t = x[i] * pow(s, -1.0 / g);
         pr[i] = (t /= s);
         rr[i] = tr * t;
         ri[i] = ti * t;
         t /= s;
         qr[i] = (trr - tii) * t;
         s = tr * ti * t;
         qi[i] = s + s;
      }

   ifftr(pr, ci, flng);

   if (a != 0.0)
      b2c(pr, n, pr, m2, a);

   if (g == 0.0 || g == -1.0) {
      movem(pr, qr, sizeof(*pr), m2 + 1);
      movem(pr, rr, sizeof(*pr), m + 1);
   } else {
      ifft(qr, qi, flng);
      ifft(rr, ri, flng);

      if (a != 0.0) {
         b2c(qr, n, qr, n, a);
         b2c(rr, n, rr, m, a);
      }
   }

   if (a != 0.0) {
      ptrans(pr, m, a);
      qtrans(qr, m, a);
   }

   /*  c[0] : gain, t : epsilon  */
   if (g != -1.0)
      c[0] = sqrt(t = gain(rr, c, m, g));

   if (g == -1.0)
      fillz(qr, sizeof(*qr), m2 + 1);
   else if (g != 0.0)
      for (i = 2; i <= m2; i++)
         qr[i] *= 1.0 + g;

   if (theq(pr, &qr[2], &b[1], &rr[1], m, f)) {
      fprintf(stderr, "mgcep : Error in theq() at %dth iteration!\n", j);
      exit(1);
   }

   for (i = 1; i <= m; i++)
      c[i] += b[i];

   /*  c[0] : gain, t : epsilon  */
   if (g == -1.0)
      c[0] = sqrt(t = gain(rr, c, m, g));

   return (log(t));
}

static double mglsadff(double x, double *b, const int m, const double a,
                       double *d)
{
   int i;
   double y, aa;

   aa = 1 - a * a;

   y = d[0] * b[1];
   for (i = 1; i < m; i++) {
      d[i] += a * (d[i + 1] - d[i - 1]);
      y += d[i] * b[i + 1];
   }
   x -= y;

   for (i = m; i > 0; i--)
      d[i] = d[i - 1];
   d[0] = a * d[0] + aa * x;

   return (x);
}

double mglsadf(double x, double *b, const int m, const double a, const int n,
               double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = mglsadff(x, b, m, a, &d[i * (m + 1)]);

   return (x);
}

static double mglsadff1(double x, double *b, const int m, const double a,
                        const double g, double *d)
{
   int i;
   double y, aa;

   aa = 1 - a * a;

   y = d[0] * b[1];
   for (i = 1; i < m; i++) {
      d[i] += a * (d[i + 1] - d[i - 1]);
      y += d[i] * b[i + 1];
   }
   x -= g * y;

   for (i = m; i > 0; i--)
      d[i] = d[i - 1];

   d[0] = a * d[0] + aa * x;

   return (x);
}

double mglsadf1(double x, double *b, const int m, const double a, const int n,
                double *d)
{
   int i;
   double g;

   g = -1.0 / (double) n;

   for (i = 0; i < n; i++)
      x = mglsadff1(x, b, m, a, g, &d[i * (m + 1)]);

   return (x);
}

static double mglsadfft(double x, double *b, const int m, const double a,
                        double *d)
{
   int i;

   x -= d[0] * (1.0 - a * a);

   d[m] = b[m] * x + a * d[m - 1];
   for (i = m - 1; i >= 1; i--)
      d[i] += b[i] * x + a * (d[i - 1] - d[i + 1]);

   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (x);
}

double mglsadft(double x, double *b, const int m, const double a, const int n,
                double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = mglsadfft(x, b, m, a, &d[i * (m + 1)]);

   return (x);
}

static double mglsadff1t(double x, double *b, const int m, const double a,
                         const double g, double *d)
{
   int i;

   x -= d[0] * (1.0 - a * a) * g;

   d[m] = b[m] * x + a * d[m - 1];
   for (i = m - 1; i >= 1; i--)
      d[i] += b[i] * x + a * (d[i - 1] - d[i + 1]);

   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (x);
}

double mglsadf1t(double x, double *b, const int m, const double a, const int n,
                 double *d)
{
   int i;
   double g;

   g = -1.0 / (double) n;

   for (i = 0; i < n; i++)
      x = mglsadff1t(x, b, m, a, g, &d[i * (m + 1)]);

   return (x);
}

static double imglsadff(double x, double *b, const int m, const double a,
                        double *d)
{
   int i;
   double y, aa;

   aa = 1 - a * a;

   y = d[0] * b[1];
   for (i = 1; i < m; i++) {
      d[i] += a * (d[i + 1] - d[i - 1]);
      y += d[i] * b[i + 1];
   }
   y += x;

   for (i = m; i > 0; i--)
      d[i] = d[i - 1];

   d[0] = a * d[0] + aa * x;

   return (y);
}

double imglsadf(double x, double *b, const int m, const double a, const int n,
                double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = imglsadff(x, b, m, a, &d[i * (m + 1)]);

   return (x);
}

static double imglsadff1(double x, double *b, const int m, const double a,
                         const double g, double *d)
{
   int i;
   double y, aa;

   aa = 1 - a * a;

   y = d[0] * b[1];
   for (i = 1; i < m; i++) {
      d[i] += a * (d[i + 1] - d[i - 1]);
      y += d[i] * b[i + 1];
   }
   y = g * y + x;

   for (i = m; i > 0; i--)
      d[i] = d[i - 1];

   d[0] = a * d[0] + aa * x;

   return (y);
}

double imglsadf1(double x, double *b, const int m, const double a, const int n,
                 double *d)
{
   int i;
   double g;

   g = -1.0 / (double) n;

   for (i = 0; i < n; i++)
      x = imglsadff1(x, b, m, a, g, &d[i * (m + 1)]);

   return (x);
}

static double imglsadfft(double x, double *b, const int m, const double a,
                         double *d)
{
   int i;
   double y;

   y = x + (1.0 - a * a) * d[0];

   d[m] = b[m] * x + a * d[m - 1];
   for (i = m - 1; i >= 1; i--)
      d[i] += b[i] * x + a * (d[i - 1] - d[i + 1]);

   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (y);
}

double imglsadft(double x, double *b, const int m, const double a, const int n,
                 double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = imglsadfft(x, b, m, a, &d[i * (m + 1)]);

   return (x);
}

static double imglsadff1t(double x, double *b, const int m, const double a,
                          const double g, double *d)
{
   int i;
   double y;

   y = x + g * (1.0 - a * a) * d[0];

   d[m] = b[m] * x + a * d[m - 1];
   for (i = m - 1; i >= 1; i--)
      d[i] += b[i] * x + a * (d[i - 1] - d[i + 1]);

   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (y);
}

double imglsadf1t(double x, double *b, const int m, const double a, const int n,
                  double *d)
{
   int i;
   double g;

   g = -1.0 / (double) n;

   for (i = 0; i < n; i++)
      x = imglsadff1t(x, b, m, a, g, &d[i * (m + 1)]);

   return (x);
}

int str2darray(char *c, double **x)
{
   int i, size, sp;
   char *p, *buf;

   while (isspace(*c))
      c++;
   if (*c == '\0') {
      *x = NULL;
      return (0);
   }

   size = 1;
   sp = 0;
   for (p = c; *p != '\0'; p++) {
      if (!isspace(*p)) {
         if (sp == 1) {
            size++;
            sp = 0;
         }
      } else
         sp = 1;
   }
   buf = getmem(strlen(c), sizeof(*buf));
   *x = dgetmem(size);
   for (i = 0; i < size; i++)
      (*x)[i] = strtod(c, &c);
   return (size);
}

int isfloat(char *c)
{
   int isnum = 0, wfe = 1;
   int i = 0;

   if (strlen(c) == 0)
      return (0);

   if ((c[i] == '+') || (c[i] == '-'))
      i++;
   while ((c[i] >= '0') && (c[i] <= '9')) {
      isnum = 1;
      i++;
   }
   if (c[i] == '.') {
      i++;
      while ((c[i] >= '0') && (c[i] <= '9')) {
         isnum = 1;
         i++;
      }
   }
   if ((c[i] == 'e') || (c[i] == 'E')) {
      wfe = 0;
      i++;
      if ((c[i] == '+') || (c[i] == '-'))
         i++;
      while ((c[i] >= '0') && (c[i] <= '9')) {
         wfe = 1;
         i++;
      }
   }
   if ((c[i] == 'f') || (c[i] == 'F') || (c[i] == 'l') || (c[i] == 'L'))
      i++;

   if ((c[i] == '\0') && isnum && wfe)
      return (1);
   else
      return (0);
}

static double mlsafir(double x, double *b, const int m, const double a,
                      double *d)
{
   double y = 0.0, aa;
   int i;

   aa = 1 - a * a;

   d[0] = x;
   d[1] = aa * d[0] + a * d[1];

   for (i = 2; i <= m; i++) {
      d[i] = d[i] + a * (d[i + 1] - d[i - 1]);
      y += d[i] * b[i];
   }

   for (i = m + 1; i > 1; i--)
      d[i] = d[i - 1];

   return (y);
}

static double mlsadf1(double x, double *b, const double a,
                      const int pd, double *d)
{
   double v, out = 0.0, *pt, aa;
   int i;

   aa = 1 - a * a;
   pt = &d[pd + 1];

   for (i = pd; i >= 1; i--) {
      d[i] = aa * pt[i - 1] + a * d[i];
      pt[i] = d[i] * b[1];
      v = pt[i] * ppade[i];

      x += (1 & i) ? v : -v;
      out += v;
   }

   pt[0] = x;
   out += x;

   return (out);
}

static double mlsadf2(double x, double *b, const int m, const double a,
                      const int pd, double *d)
{
   double v, out = 0.0, *pt;
   int i;

   pt = &d[pd * (m + 2)];

   for (i = pd; i >= 1; i--) {
      pt[i] = mlsafir(pt[i - 1], b, m, a, &d[(i - 1) * (m + 2)]);
      v = pt[i] * ppade[i];

      x += (1 & i) ? v : -v;
      out += v;
   }

   pt[0] = x;
   out += x;

   return (out);
}

double mlsadf(double x, double *b, const int m, const double a, const int pd,
              double *d)
{
   ppade = &pade[pd * (pd + 1) / 2];

   x = mlsadf1(x, b, a, pd, d);
   x = mlsadf2(x, b, m, a, pd, &d[2 * (pd + 1)]);

   return (x);
}


static double mlsafirt(double x, double *b, const int m, const double a,
                       double *d)
{
   int i;
   double y = 0.0;

   y = (1.0 - a * a) * d[0];

   d[m] = b[m] * x + a * d[m - 1];
   for (i = m - 1; i > 1; i--)
      d[i] += b[i] * x + a * (d[i - 1] - d[i + 1]);
   d[1] += a * (d[0] - d[2]);

   for (i = 0; i < m; i++)
      d[i] = d[i + 1];

   return (y);
}

static double mlsadf2t(double x, double *b, const int m, const double a,
                       const int pd, double *d)
{
   double v, out = 0.0, *pt;
   int i;

   pt = &d[pd * (m + 2)];

   for (i = pd; i >= 1; i--) {
      pt[i] = mlsafirt(pt[i - 1], b, m, a, &d[(i - 1) * (m + 2)]);
      v = pt[i] * ppade[i];

      x += (1 & i) ? v : -v;
      out += v;
   }

   pt[0] = x;
   out += x;

   return (out);
}

double mlsadft(double x, double *b, const int m, const double a, const int pd,
               double *d)
{
   ppade = &pade[pd * (pd + 1) / 2];

   x = mlsadf1(x, b, a, pd, d);
   x = mlsadf2t(x, b, m, a, pd, &d[2 * (pd + 1)]);

   return (x);
}

void msvq(double *x, double *cb, const int l, int *cbsize, const int stage,
          int *index)
{
   int i, j;
   double *p;
   static double *xx = NULL;
   static int size;

   if (xx == NULL) {
      xx = dgetmem(l);
      size = l;
   }
   if (size > l) {
      free(xx);
      xx = dgetmem(l);
      size = l;
   }

   movem(x, xx, sizeof(*x), l);

   for (i = 0; i < stage; i++) {
      index[i] = vq(xx, cb, l, cbsize[i]);

      p = cb + index[i] * l;
      for (j = 0; j < l; j++)
         xx[j] -= p[j];

      cb += cbsize[i] * l;
   }
}

void ndps2c(double *n, const int l, double *c, const int m)
{
   int i, no;
   double *nx, *ny;

   no = l / 2;
   nx = dgetmem(l);
   ny = dgetmem(l);

   for (i = 0; i <= no; i++) {
      nx[i] = n[i];
   }

   for (i = 1; i < no; i++) {
      nx[l - i] = nx[i];
   }

   fftr(nx, ny, l);

   c[0] = 0.0;
   for (i = 1; i <= m; i++)
      c[i] = nx[i] / (i * l / 2.0);
   if (m == l / 2)
      c[m] = c[m] / 2.0;

   free(nx);
   free(ny);
}

void norm0(double *x, double *y, int m)
{
   y[0] = 1 / x[0];
   for (; m >= 1; m--)
      y[m] = x[m] * y[0];

   return;
}

static double rnd(unsigned long *next)
{
   double r;

   *next = *next * 1103515245L + 12345;
   r = (*next / 65536L) % 32768L;

   return (r / _RAND_MAX);
}

int nrand(double *p, const int leng, const int seed)
{
   int i;
   unsigned long next;

   if (seed != 1)
      next = srnd((unsigned int) seed);
   for (i = 0; i < leng; i++)
      p[i] = (double) nrandom(&next);

   return (0);
}

double nrandom(unsigned long *next)
{
   static int sw = 0;
   static double r1, r2, s;

   if (sw == 0) {
      sw = 1;
      do {
         r1 = 2 * rnd(next) - 1;
         r2 = 2 * rnd(next) - 1;
         s = r1 * r1 + r2 * r2;
      }
      while (s > 1 || s == 0);
      s = sqrt(-2 * log(s) / s);
      return (r1 * s);
   } else {
      sw = 0;
      return (r2 * s);
   }
}

unsigned long srnd(const unsigned int seed)
{
   return (seed);
}

void par2lpc(double *k, double *a, const int m)
{
   int i, n;

   a[0] = k[0];
   for (n = 1; n <= m; n++) {
      for (i = 1; i < n; i++)
         a[i] = k[i] + k[n] * k[n - i];
      movem(&a[1], &k[1], sizeof(*a), n - 1);
   }
   a[m] = k[m];

   return;
}

void phase(double *p, const int mp, double *z, const int mz, double *ph,
           const int flng, const int unlap)
{
   static double *x;
   static int fsize = 0;
   double *y, *xx, *yy, *py;
   int no, i, offset;
   double pi;

   pi = atan(1.) * 4.;

   no = flng / 2 + 1;

   if (flng > fsize) {
      if (x != NULL)
         free(x);
      fsize = flng;
      x = dgetmem(4 * flng + no);
   }
   y = &x[flng];
   xx = &y[flng];
   yy = &xx[flng];
   py = &yy[flng];

   fillz(x, sizeof(*x), flng);
   fillz(xx, sizeof(*xx), flng);
   movem(z, x, mz + 1, sizeof(*z));
   movem(p, xx, mp + 1, sizeof(*p));

   fftr(x, y, flng);
   xx[0] = 1;
   fftr(xx, yy, flng);
   for (i = 0; i < no; i++) {
      ph[i] = x[i] * xx[i] + y[i] * yy[i];
      py[i] = y[i] * xx[i] - x[i] * yy[i];
   }
   offset = 0;
   i = 0;
   ph[i] = atan2(py[i], ph[i]) / pi;
   i++;
   for (; i < no; i++) {
      ph[i] = atan2(py[i], ph[i]) / pi;
      if (unlap) {
         if (ph[i - 1] - ph[i] - offset > 1)
            offset += 2;
         else if (ph[i] + offset - ph[i - 1] > 1)
            offset -= 2;
         ph[i] += offset;
      }
   }

   return;
}

double poledf(double x, double *a, int m, double *d)
{
   for (m--; m > 0; m--) {
      x -= a[m + 1] * d[m];
      d[m] = d[m - 1];
   }
   x -= a[1] * d[0];
   d[0] = x;

   return (x);
}

double poledft(double x, double *a, int m, double *d)
{
   int i;

   x -= d[0];
   for (i = 1; i < m; i++)
      d[i - 1] = d[i] + a[i] * x;
   d[m - 1] = a[m] * x;

   return (x);
}

void reverse(double *x, const int l)
{
   int i = 0;
   double d;

   while (i < l - i - 1) {
      d = x[i];
      x[i] = x[l - i - 1];
      x[l - i - 1] = d;
      i++;
   }

   return;
}

int rlevdur(double *a, double *r, const int m, double eps)
{
   int i, j;
   double rmd, sum;
   static double **u = NULL, *e = NULL;
   static int size;

   if (u == NULL && e == NULL) {
      u = ddgetmem(m + 1, m + 1);
      e = dgetmem(m + 1);
      size = m;
   }

   if (m > size) {
      free(u);
      free(e);
      u = ddgetmem(m + 1, m + 1);
      e = dgetmem(m + 1);
      size = m;
   }

   for (j = 0; j <= m; j++) {
      u[j][j] = 1.0;
   }

   for (j = 0; j < m; j++) {
      u[m][j] = a[m - j];
   }
   e[m] = a[0] * a[0];

   for (i = m - 1; i > 0; i--) {
      rmd = (1.0 - u[i + 1][0] * u[i + 1][0]);
      if ((rmd < 0.0) ? -rmd : rmd <= eps)
         return (-1);
      for (j = 0; j < i; j++) {
         u[i][i - j - 1] =
             (u[i + 1][i - j] - u[i + 1][0] * u[i + 1][j + 1]) / rmd;
      }
      e[i] = e[i + 1] / rmd;
   }
   e[0] = e[1] / (1.0 - u[1][0] * u[1][0]);

   r[0] = e[0];
   for (i = 1; i <= m; i++) {
      sum = 0.0;
      for (j = 1; j < i; j++) {
         sum -= u[i - 1][i - j - 1] * r[i - j];
      }
      r[i] = sum - u[i][0] * e[i - 1];
   }

   return (0);
}

double rmse(double *x, double *y, const int n)
{
   int i;
   double sub, sum;

   sum = 0.0;
   for (i = 0; i < n; i++) {
      sub = x[i] - y[i];
      sum += sub * sub;
   }

   return (sqrt(sum / n));
}

typedef enum { plus, minus, multiply, divide } opt;

static double rad_root(const double x, const int i)
{
   if (x == 0.0)
      return -1.0;
   else
      return exp(log(x) / i);
}

static complex c_math(complex c1, opt op, complex c2)
{
   double p;
   complex t;

   switch (op) {
   case plus:
      t.re = c1.re + c2.re;
      t.im = c1.im + c2.im;
      break;
   case minus:
      t.re = c1.re - c2.re;
      t.im = c1.im - c2.im;
      break;
   case multiply:
      t.re = c1.re * c2.re - c1.im * c2.im;
      t.im = c1.re * c2.im + c1.im * c2.re;
      break;
   case divide:
      p = c2.re * c2.re + c2.im * c2.im;
      t.re = (c1.re * c2.re + c1.im * c2.im) / p;
      t.im = (c1.im * c2.re - c1.re * c2.im) / p;
      break;
   default:
      t.re = c1.re;
      t.im = c1.im;
      break;
   }
   return t;
}

static double c_mag(complex x)
{
   return sqrt(x.re * x.re + x.im * x.im);
}

static double c_arg(complex x)
{
   return atan2(x.im, x.re);
}

void output_root_pol(complex * x, int odr, int form)
{
   int i, k;
   double mag, arg, *a;

   a = dgetmem(2 * odr);

   switch (form) {
   case 1:
      for (k = i = 0; i < odr; i++) {
         a[k++] = c_mag(x[i + 1]);
         a[k++] = c_arg(x[i + 1]);
      }
      break;
   case 2:
   case 3:
      for (k = i = 0; i < odr; i++) {
         mag = 1 / c_mag(x[i + 1]);
         arg = -c_arg(x[i + 1]);
         if (form == 3) {
            a[k++] = mag;
            a[k++] = arg;
         } else {
            a[k++] = mag * cos(arg);
            a[k++] = mag * sin(arg);
         }
      }
      break;
   case 0:
   default:
      for (k = i = 0; i < odr; i++) {
         a[k++] = x[i + 1].re;
         a[k++] = x[i + 1].im;
      }
      break;
   }

   fwritef(a, sizeof(*a), odr * 2, stdout);

   return;
}

complex *cplx_getmem(const int leng)
{
   int i;
   complex *p = NULL;

   if ((p = (complex *) malloc(sizeof(complex) * leng)) == NULL) {
      fprintf(stderr, "root_pol : Cannot allocate memory!\n");
      exit(3);
   }

   for (i = 0; i < leng; i++)
      p[i].re = p[i].im = 0;

   return p;
}

void root_pol(double *a, const int odr, complex * x, const int a_zero,
              const double eps, const int itrat)
{
   int i, j, k, l;
   double th, th1, th2, cm, cmax;
   complex cden, cnum, c1, *deltx;

   deltx = cplx_getmem(odr + 1);

   if (!a_zero)
      for (i = 1; i <= odr; i++)
         a[i] /= a[0];

   cmax = 0;
   for (i = 2; i <= odr; i++) {
      cm = odr * rad_root(fabs(a[i]), i);
      if (cm > cmax)
         cmax = cm;
   }

   th1 = PI * 2.0 / odr;
   th2 = th1 / 4.0;
   for (i = 1; i <= odr; i++) {
      th = th1 * (i - 1) + th2;
      x[i].re = cmax * cos(th);
      x[i].im = cmax * sin(th);
   }

   l = 1;
   do {
      for (i = 1; i <= odr; i++) {
         cden.re = 1.0;
         cden.im = 0.0;
         cnum.re = 1.0;
         cnum.im = 0.0;
         c1 = x[i];
         for (j = 1; j <= odr; j++) {
            cnum = c_math(cnum, multiply, c1);
            cnum.re += a[j];
            if (j != i)
               cden = c_math(cden, multiply, c_math(c1, minus, x[j]));
         }
         deltx[i] = c_math(cnum, divide, cden);
         x[i] = c_math(c1, minus, deltx[i]);
      }
      k = 1;
      while ((k <= odr) && (c_mag(deltx[k++]) <= eps));
      l++;
   }
   while ((l <= itrat) && (k <= odr));

   if (l > itrat) {
      fprintf(stderr, "root_pol : No convergence!\n");
      exit(1);
   }

   return;
}


static double warp(const double w, const double a, const double t)
{
   double ww, x, y;

   x = w - t;
   y = w + t;

   ww = w + atan2((a * sin(x)), (1.0 - a * cos(x)))
       + atan2((a * sin(y)), (1.0 - a * cos(y)));

   return (ww);
}


/*============================================================*/

static double derivw(const double w, const double a, const double t)
{
   double dw, x, y, a2, aa;

   x = w - t;
   y = w + t;

   a2 = a + a;
   aa = a * a;

   dw = 1.0 + (a * cos(x) - aa) / (1.0 - a2 * cos(x) + aa)
       + (a * cos(y) - aa) / (1.0 - a2 * cos(y) + aa);

   return (dw);
}

/***************************************************************

  No.1  frqt_a    static : *l, size1

  Frequency Transformation of "al" (second term of dE/dc)

      void frqt_a(al, m, fftsz, a, t)

      double *al   : sequence which will be warped
      int m        : order of warped sequence
      int fftsz    : ifft size
      double a     : all-pass constant
      double t     : emphasized frequency (t * pi)

***************************************************************/

static void frqt_a(double *al, const int m, const int fftsz, const double a,
                   const double t)
{
   int i, j;
   double w, b, *ww, *f, *re, *im, *pf, *pl, *next;
   int size_l, size_f, fftsz2;
   static double *l = NULL;
   static int size1, flag_l = 1;

   b = M_2PI / (double) fftsz;

   size_l = m + 1;

   if (l == NULL) {
      flag_l = 0;
      size1 = size_l;
      l = dgetmem(size1);
   } else if (size_l != size1) {
      free(l);
      flag_l = 0;
      size1 = size_l;
      l = dgetmem(size1);
   }

   /*-------  if "l" is not defined  ----------*/

   if (flag_l == 0) {

      ww = dgetmem(fftsz);

      for (j = 0, w = 0.0; j < fftsz; j++, w += b)
         ww[j] = warp(w, a, t);

      fftsz2 = fftsz + fftsz;   /* size of (re + im) */
      size_f = (m + 1) * fftsz2;        /* size of array "f" */
      f = dgetmem(size_f);

      for (i = 0, re = f, im = f + fftsz; i <= m; i++) {

         for (j = 0; j < fftsz; j++)
            *(re++) = cos(ww[j] * i);
         for (j = 0; j < fftsz; j++)
            *(im++) = -sin(ww[j] * i);

         re -= fftsz;
         im -= fftsz;

         ifft(re, im, fftsz);

         re += fftsz2;
         im += fftsz2;
      }

      free(ww);


      /*-------  copy "f" to "l" ----------*/

      for (i = 0, next = f, pf = f, pl = l; i <= m; i++) {
         *(pl++) = *pf;
         next += fftsz2;
         pf = next;
      }

      free(f);
      flag_l = 1;
   }

   movem(l, al, sizeof(*al), m + 1);

   return;
}

/***************************************************************

  No.2  freqt2    static : *g, size2

  Frequency Transformation

      void freqt2(c1, m1, c2, m2, fftsz, a, t)

      double *c1   : minimum phase sequence
      int    m1    : order of minimum phase sequence
      double *c2   : warped sequence
      int    m2    : order of warped sequence
      int    fftsz : ifft size
      double a     : all-pass constant
      double t     : emphasized frequency (t * pi)

***************************************************************/

static void freqt2(double *c1, const int m1, double *c2, const int m2,
                   const int fftsz, const double a, const double t)
{
   int i, j;
   double w, b, *ww, *dw, *f, *re, *im, *pf, *pg, *next;
   int size_g, size_f, fftsz2;
   static double *g = NULL;
   static int size2, flag_g = 1;

   b = M_2PI / (double) fftsz;

   size_g = (m2 + 1) * (m1 + 1);

   if (g == NULL) {
      flag_g = 0;
      size2 = size_g;
      g = dgetmem(size2);
   } else if (size_g != size2) {
      free(g);
      flag_g = 0;
      size2 = size_g;
      g = dgetmem(size2);
   }

   /*-------  if "g" is not defined  ----------*/

   if (flag_g == 0) {
      ww = dgetmem(fftsz);
      dw = dgetmem(fftsz);

      for (j = 0, w = 0.0; j < fftsz; j++, w += b)
         ww[j] = warp(w, a, t);

      for (j = 0, w = 0.0; j < fftsz; j++, w += b)
         dw[j] = derivw(w, a, t);


      fftsz2 = fftsz + fftsz;   /* size of (re + im) */
      size_f = (m2 + 1) * fftsz2;       /* size of array "f" */
      f = dgetmem(size_f);

      for (i = 0, re = f, im = f + fftsz; i <= m2; i++) {

         for (j = 0; j < fftsz; j++)
            *(re++) = cos(ww[j] * i) * dw[j];
         for (j = 0; j < fftsz; j++)
            *(im++) = -sin(ww[j] * i) * dw[j];

         re -= fftsz;
         im -= fftsz;

         ifft(re, im, fftsz);

         for (j = 1; j <= m1; j++)
            re[j] += re[fftsz - j];

         re += fftsz2;
         im += fftsz2;
      }

      free(ww);
      free(dw);


      /*-------  copy "f" to "g" ----------*/

      for (i = 0, next = f, pf = f, pg = g; i <= m2; i++) {
         for (j = 0; j <= m1; j++)
            *(pg++) = *(pf++);
         next += fftsz2;
         pf = next;
      }
      free(f);
      flag_g = 1;

      for (j = 1; j <= m1; j++)
         g[j] *= 0.5;

      for (i = 1; i <= m2; i++)
         g[i * (m1 + 1)] *= 2.0;
   }

   for (i = 0, pg = g; i <= m2; i++)
      for (j = 0, c2[i] = 0.0; j <= m1; j++)
         c2[i] += *(pg++) * c1[j];

   return;
}

/***************************************************************

  No.3  ifreqt2    static : *h, size3

  Inverse Frequency Transformation

      void ifreqt2(c1, m1, c2, m2, fftsz, a, t)

      double *c1   : minimum phase sequence
      int    m1    : order of minimum phase sequence
      double *c2   : warped sequence
      int    m2    : order of warped sequence
      int    fftsz : ifft size
      double a     : all-pass constant
      double t     : emphasized frequency t * pi(rad)

***************************************************************/

static void ifreqt2(double *c1, int m1, double *c2, int m2, int fftsz, double a,
                    double t)
{
   int i, j;
   double w, b, *ww, *f, *re, *im, *pl, *pr, *plnxt, *prnxt, *pf, *ph, *next;
   int size_h, size_f, fftsz2, m12, m11;
   static double *h = NULL;
   static int size3, flag_h = 1;

   b = M_2PI / (double) fftsz;

   size_h = (m2 + 1) * (m1 + 1);

   if (h == NULL) {
      flag_h = 0;
      size3 = size_h;
      h = dgetmem(size3);
   } else if (size_h != size3) {
      free(h);
      flag_h = 0;
      size3 = size_h;
      h = dgetmem(size3);
   }

   /*-------  if "h" is not defined  ----------*/

   if (flag_h == 0) {
      ww = dgetmem(fftsz);

      for (j = 0, w = 0.0; j < fftsz; j++, w += b)
         ww[j] = warp(w, a, t);

      fftsz2 = fftsz + fftsz;   /* size of (re + im) */

      m12 = m1 + m1 + 1;
      size_f = m12 * fftsz2;    /* size of array "f" */
      f = dgetmem(size_f);

      for (i = -m1, re = f, im = f + fftsz; i <= m1; i++) {

         for (j = 0; j < fftsz; j++)
            *(re++) = cos(ww[j] * i);

         for (j = 0; j < fftsz; j++)
            *(im++) = -sin(ww[j] * i);

         re -= fftsz;
         im -= fftsz;

         ifft(re, im, fftsz);

         re += fftsz2;
         im += fftsz2;
      }

      free(ww);

      /*------- b'(n,m)=b(n,m)+b(n,-m) ----------*/

      pl = f;
      pr = f + (m12 - 1) * fftsz2;

      for (i = 0, plnxt = pl, prnxt = pr; i < m1; i++) {
         plnxt += fftsz2;
         prnxt -= fftsz2;

         for (j = 0; j <= m2; j++)
            *(pr++) += *(pl++);

         pl = plnxt;
         pr = prnxt;
      }

      /*-------  copy "f" to "h" ----------*/

      m11 = m1 + 1;
      pf = f + m1 * fftsz2;

      for (j = 0, next = pf; j <= m1; j++) {

         next += fftsz2;

         for (i = 0; i <= m2; i++)
            h[m11 * i + j] = *(pf++);

         pf = next;
      }
      free(f);
      flag_h = 1;

      for (j = 1; j <= m1; j++)
         h[j] *= 0.5;

      for (i = 1; i <= m2; i++)
         h[i * m11] *= 2.0;
   }

   for (i = 0, ph = h; i <= m2; i++)
      for (j = 0, c2[i] = 0.0; j <= m1; j++)
         c2[i] += *(ph++) * c1[j];

   return;
}

/***************************************************************

  No.4  frqtr2    static : *k, size4

  Frequency Transformation for Calculating Coefficients

      void frqtr2(c1, m1, c2, m2, fftsz, a, t)

      double *c1   : minimum phase sequence
      int    m1    : order of minimum phase sequence
      double *c2   : warped sequence
      int    m2    : order of warped sequence
      int    fftsz : frame length (fft size)
      double a     : all-pass constant
      double t     : emphasized frequency

***************************************************************/

static void frqtr2(double *c1, int m1, double *c2, int m2, int fftsz, double a,
                   double t)
{
   int i, j;
   double w, b, *ww, *f, *tc2, *re, *im, *pf, *pk, *next;
   int size_k, size_f, fftsz2;
   static double *k = NULL;
   static int size4, flag_k = 1;

   b = M_2PI / (double) fftsz;

   size_k = (m2 + 1) * (m1 + 1);

   if (k == NULL) {
      flag_k = 0;
      size4 = size_k;
      k = dgetmem(size4);
   } else if (size_k != size4) {
      free(k);
      flag_k = 0;
      size4 = size_k;
      k = dgetmem(size4);
   }

   /*-------  if "k" is not defined  ----------*/

   if (flag_k == 0) {

      ww = dgetmem(fftsz);

      for (j = 0, w = 0.0; j < fftsz; j++, w += b)
         ww[j] = warp(w, a, t);

      fftsz2 = fftsz + fftsz;   /* size of (re + im) */
      size_f = (m2 + 1) * fftsz2;       /* size of array "f" */
      f = dgetmem(size_f);

      for (i = 0, re = f, im = f + fftsz; i <= m2; i++) {

         for (j = 0; j < fftsz; j++)
            *(re++) = cos(ww[j] * i);
         for (j = 0; j < fftsz; j++)
            *(im++) = -sin(ww[j] * i);

         re -= fftsz;
         im -= fftsz;

         ifft(re, im, fftsz);

         for (j = 1; j <= m1; j++)
            re[j] += re[fftsz - j];

         re += fftsz2;
         im += fftsz2;
      }

      free(ww);


      /*-------  copy "f" to "k" ----------*/

      for (i = 0, next = f, pf = f, pk = k; i <= m2; i++) {
         for (j = 0; j <= m1; j++)
            *(pk++) = *(pf++);
         next += fftsz2;
         pf = next;
      }
      free(f);
      flag_k = 1;
   }

   tc2 = dgetmem(m2 + 1);       /*  tmp of c2  */

   for (i = 0, pk = k; i <= m2; i++)
      for (j = 0, tc2[i] = 0.0; j <= m1; j++)
         tc2[i] += *(pk++) * c1[j];

   movem(tc2, c2, sizeof(*c2), m2 + 1);

   free(tc2);

   return;
}

int smcep(double *xw, const int flng, double *mc, const int m, const int fftsz,
          const double a, const double t, const int itr1, const int itr2,
          const double dd, const int etype, const double e, const double f,
          const int itype)
{
   int i, j;
   int flag = 0, f2, m2;
   double u, s, eps = 0.0, min, max;
   static double *x = NULL, *y, *c, *d, *al, *b;
   static int size_x, size_d;

   if (etype == 1 && e < 0.0) {
      fprintf(stderr, "smcep : value of e must be e>=0!\n");
      exit(1);
   }

   if (etype == 2 && e >= 0.0) {
      fprintf(stderr, "smcep : value of E must be E<0!\n");
      exit(1);
   }

   if (etype == 1) {
      eps = e;
   }


   if (x == NULL) {
      x = dgetmem(3 * flng);
      y = x + flng;
      c = y + flng;
      size_x = flng;

      d = dgetmem(3 * m + 3);
      al = d + (m + 1);
      b = al + (m + 1);
      size_d = m;
   }
   if (flng > size_x) {
      free(x);
      x = dgetmem(3 * flng);
      y = x + flng;
      c = y + flng;
      size_x = flng;
   }
   if (m > size_d) {
      free(d);
      d = dgetmem(3 * m + 3);
      al = d + (m + 1);
      b = al + (m + 1);
      size_d = m;
   }

   f2 = flng / 2.;
   m2 = m + m;

   movem(xw, x, sizeof(*x), flng);

   switch (itype) {
   case 0:                     /* windowed data sequence */
      fftr(x, y, flng);
      for (i = 0; i < flng; i++) {
         x[i] = x[i] * x[i] + y[i] * y[i] + eps;        /*  periodogram  */
      }
      break;
   case 1:                     /* dB */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp((x[i] / 20.0) * log(10.0)); /* dB -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 2:                     /* log */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp(x[i]);      /* log -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 3:                     /* amplitude */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 4:                     /* periodogram */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] + eps;
      }
      break;
   default:
      fprintf(stderr, "smcep : Input type %d is not supported!\n", itype);
      exit(1);
   }
   if (itype > 0) {
      for (i = 1; i < flng / 2; i++)
         x[flng - i] = x[i];
   }

   if (etype == 2 && e < 0.0) {
      max = x[0];
      for (i = 1; i < flng; i++) {
         if (max < x[i])
            max = x[i];
      }
      max = sqrt(max);
      min = max * pow(10.0, e / 20.0);  /* floor is 20*log10(min/max) */
      min = min * min;
      for (i = 0; i < flng; i++) {
         if (x[i] < min)
            x[i] = min;
      }
   }

   for (i = 0; i < flng; i++)
      c[i] = log(x[i]);

   /*  1, (-a), (-a)^2, ..., (-a)^M  */

   al[0] = 1.0;
   for (i = 1; i <= m; i++)
      al[i] = 0.0;

   frqt_a(al, m, fftsz, a, t);


   /*  initial value of cepstrum  */
   ifftr(c, y, flng);           /*  c : IFFT[x]  */

   c[0] /= 2.0;
   c[flng / 2] /= 2.0;
   freqt2(c, f2, mc, m, fftsz, a, t);   /*  mc : mel cep.  */

   s = c[0];

   /*  Newton Raphson method  */
   for (j = 1; j <= itr2; j++) {
      fillz(c, sizeof(*c), flng);
      ifreqt2(mc, m, c, f2, fftsz, a, t);       /*  mc : mel cep.  */

      fftr(c, y, flng);         /*  c, y : FFT[mc]  */
      for (i = 0; i < flng; i++)
         c[i] = x[i] / exp(c[i] + c[i]);
      ifftr(c, y, flng);
      frqtr2(c, f2, c, m2, fftsz, a, t);        /*  c : r(k)  */

      u = c[0];
      if (j >= itr1) {
         if (fabs((u - s) / u) < dd) {
            flag = 1;
            break;
         }
         s = u;
      }

      for (i = 0; i <= m; i++)
         b[i] = c[i] - al[i];
      for (i = 0; i <= m2; i++)
         y[i] = c[i];
      for (i = 0; i <= m2; i += 2)
         y[i] -= c[0];
      for (i = 2; i <= m; i += 2)
         c[i] += c[0];
      c[0] += c[0];

      if (theq(c, y, d, b, m + 1, f)) {
         fprintf(stderr, "smcep : Error in theq() at %dth iteration!\n", j);
         exit(1);
      }

      for (i = 0; i <= m; i++)
         mc[i] += d[i];
   }

   if (flag)
      return (0);
   else
      return (-1);
}

/* Fast Algorithm for Linear Prediction with Linear Phase */
static void lplp(double *r, double *c, const int m)
{
   int k, n;
   double pn, alpha, beta, gamma, tz = r[0] / 2, rtz = 1 / tz, to = r[1], rttz =
       0, tto = 1;
   static double *p = NULL, *pp;
   static int size;

   if (p == NULL) {
      p = dgetmem(m + m + 4);
      pp = p + m + 2;
      size = m;
   }
   if (m > size) {
      free(p);
      p = dgetmem(m + m + 2);
      pp = p + m + 1;
      size = m;
   }

   c[0] = 1.0 / r[0];
   p[0] = 1.0;
   p[1] = 0.0;
   pp[0] = 0.0;

   for (k = 1; k <= m; k++) {
      p[k + 1] = 0.0;

      pp[k] = 0.0;
      beta = -tz * rttz;
      alpha = tto * rttz;
      alpha -= (tto = to) * (rttz = rtz);
      pn = p[1] + p[1] + alpha * p[0] + beta * pp[0];
      pp[0] = p[0];
      p[0] = pn;

      for (n = 1; n <= k; n++) {
         pn = p[n + 1] + pp[n - 1] + alpha * p[n] + beta * pp[n];
         pp[n] = p[n];
         p[n] = pn;
      }
      for (n = 1, tz = p[0] * r[k]; n <= k; n++)
         tz += p[n] * (r[k - n] + r[k + n]);

      for (n = 1, to = p[0] * r[1 + k]; n <= k; n++)
         to += p[n] * (r[1 + k - n] + r[1 + k + n]);

      gamma = 0.5 * p[0] * (rtz = 1 / tz);

      for (n = 0; n < k; n++)
         c[n] = c[n] + gamma * p[n];

      c[k] = gamma * p[k];
   }
   for (c[0] = 1.0 / c[0], n = 1; n <= m; n++)
      c[n] *= c[0];

   return;
}

int uels(double *xw, const int flng, double *c, const int m, const int itr1,
         const int itr2, const double dd, const int etype, const double e,
         const int itype)
{
   int i, j, flag = 0;
   double k, eps = 0.0, min, max;
   static double *x = NULL, *r, *cr, *y, *a;
   static int size_x, size_a;

   if (etype == 1 && e < 0.0) {
      fprintf(stderr, "uels : value of e must be e>=0!\n");
      exit(1);
   }

   if (etype == 2 && e >= 0.0) {
      fprintf(stderr, "uels : value of E must be E<0!\n");
      exit(1);
   }

   if (etype == 1) {
      eps = e;
   }

   if (x == NULL) {
      x = dgetmem(4 * flng);
      a = dgetmem(m + 1);
      r = x + flng;
      cr = r + flng;
      y = cr + flng;
      size_x = flng;
      size_a = m;
   }
   if (flng > size_x) {
      free(x);
      x = dgetmem(4 * flng);
      r = x + flng;
      cr = r + flng;
      y = cr + flng;
      size_x = flng;
   }
   if (m > size_a) {
      free(a);
      a = dgetmem(m + 1);
      size_a = m;
   }

   movem(xw, x, sizeof(*xw), flng);

   switch (itype) {
   case 0:                     /* windowed data sequence */
      fftr(x, y, flng);
      for (i = 0; i < flng; i++) {
         x[i] = x[i] * x[i] + y[i] * y[i] + eps;        /*  periodogram  */
      }
      break;
   case 1:                     /* dB */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp((x[i] / 20.0) * log(10.0)); /* dB -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 2:                     /* log */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = exp(x[i]);      /* log -> amplitude spectrum */
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 3:                     /* amplitude */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] * x[i] + eps;      /* amplitude -> periodogram */
      }
      break;
   case 4:                     /* periodogram */
      for (i = 0; i <= flng / 2; i++) {
         x[i] = x[i] + eps;
      }
      break;
   default:
      fprintf(stderr, "uels : Input type %d is not supported!\n", itype);
      exit(1);
   }
   if (itype > 0) {
      for (i = 1; i < flng / 2; i++)
         x[flng - i] = x[i];
   }

   if (etype == 2 && e < 0.0) {
      max = x[0];
      for (i = 1; i < flng; i++) {
         if (max < x[i])
            max = x[i];
      }
      max = sqrt(max);
      min = max * pow(10.0, e / 20.0);  /* floor is 20*log10(min/max) */
      min = min * min;
      for (i = 0; i < flng; i++) {
         if (x[i] < min)
            x[i] = min;
      }
   }

   for (i = 0; i < flng; i++) {
      if (x[i] <= 0) {
         fprintf(stderr,
                 "uels : The log periodogram has '0', use '-e' option!\n");
         exit(1);
      }
      x[i] = cr[i] = log(x[i]);
   }
   ifftr(cr, y, flng);          /*  cr : c(m)  */

   /*  initial value  */
   k = exp(cr[0]);
   for (i = 1; i <= m; i++)
      c[i] = cr[i];

   for (j = 1; j <= itr2; j++) {
      cr[0] = 0.0;

      for (i = 1; i <= m; i++)
         cr[i] = c[i];
      for (; i < flng; i++)
         cr[i] = 0.0;

      fftr(cr, y, flng);        /*  cr+jy : log D(z)  */
      for (i = 0; i < flng; i++)
         r[i] = exp(x[i] - cr[i] - cr[i]);
      ifftr(r, y, flng);        /*  r : autocorr  */

      c[0] = k;
      k = r[0];

      if (j >= itr1) {
         if (fabs((k - c[0]) / c[0]) < dd) {
            flag = 1;
            break;
         }
         k = c[0];
      }

      lplp(r, a, m);
      for (i = 1; i <= m; i++)
         c[i] -= a[i];
   }

   c[0] = 0.5 * log(k);
   if (flag)
      return (0);
   else
      return (-1);
}

double ulaw_c(const double x, const double max, const double mu)
{
   double y;

   y = sign(x) * max * log(1 + mu * abs(x) / max) / log(1 + mu);
   return (y);
}

double ulaw_d(const double x, const double max, const double mu)
{
   double y;

   y = sign(x) * max * (pow(1 + mu, abs(x) / max) - 1) / mu;
   return (y);
}

int vq(double *x, double *cb, const int l, const int cbsize)
{
   int i, index = 0;
   double min = 1e23, dist;

   for (i = 0; i < cbsize; i++) {
      dist = edist(x, cb, l);
      if (dist < min) {
         index = i;
         min = dist;
      }
      cb += l;
   }
   return (index);
}

double edist(double *x, double *y, const int m)
{
   int i;
   double sub, dist = 0.0;

   for (i = 0; i < m; i++) {
      sub = x[i] - y[i];
      dist += sub * sub;
   }

   return (dist);
}

/************************************************
   Blackman window

       double  *blackman(w, leng)

       double  *w   : window values
       int     leng : window length
************************************************/

static double *blackman(double *w, const int leng)
{
   int i;
   double arg, x;
   double *p;

   arg = M_2PI / (leng - 1);
   for (p = w, i = 0; i < leng; i++) {
      x = arg * i;
      *p++ = 0.42 - 0.50 * cos(x) + 0.08 * cos(x + x);
   }
   return (w);
}


/************************************************
   Hamming window

       double  *hamming(w, leng)
       double  *w   : window values
       int     leng : window length
************************************************/

static double *hamming(double *w, const int leng)
{
   int i;
   double arg;
   double *p;

   arg = M_2PI / (leng - 1);
   for (p = w, i = 0; i < leng; i++)
      *p++ = 0.54 - 0.46 * cos(i * arg);

   return (w);
}


/************************************************
   Hanning window

       double  *hanning(w, leng)
       double  *w   : window values
       int     leng : window length
************************************************/

static double *hanning(double *w, const int leng)
{
   int i;
   double arg;
   double *p;

   arg = M_2PI / (leng - 1);
   for (p = w, i = 0; i < leng; i++)
      *p++ = 0.5 * (1 - cos(i * arg));

   return (w);
}


/************************************************
   Bartlett window

       double  *bartlett(w, leng)
       double  *w   : window values
       int     leng : window length
************************************************/

static double *bartlett(double *w, const int leng)
{
   int k, m;
   double *p, slope;

   m = leng / 2;
   slope = 2.0 / (double) (leng - 1);

   for (k = 0, p = w; k < m; k++)
      *p++ = slope * k;
   for (; k < leng; k++)
      *p++ = 2.0 - slope * k;

   return (w);
}


/************************************************
   trapezoid window

       double  *trapezoid(w, leng)
       double  *w   : window values
       int     leng : window length
************************************************/

static double *trapezoid(double *w, const int leng)
{
   int k, m1, m2;
   double *p, slope;

   m1 = leng / 4;
   m2 = (leng * 3) / 4;
   slope = 4.0 / (double) (leng - 1);

   for (k = 0, p = w; k < m1; k++)
      *p++ = slope * k;
   for (; k < m2; k++)
      *p++ = 1.0;
   for (; k < leng; k++)
      *p++ = 4.0 - slope * k;

   return (w);
}


/************************************************
   rectangular window

       double  *rectangular(w, leng)
       double  *w   : window values
       int     leng : window length
************************************************/

static double *rectangular(double *w, const int leng)
{
   int k;
   double *p;

   for (k = 0, p = w; k < leng; k++)
      *p++ = 1.0;

   return (w);
}

double window(Window type, double *x, const int size, const int nflg)
{
   int i;
   static double g;
   static double *w = NULL;
   static Window ptype = (Window) - 1;
   static int psize = -1, pnflg = -1;

   if ((type != ptype) || (size != psize) || (nflg != pnflg)) {
      if (size > psize) {
         if (w != NULL)
            free(w);
         w = dgetmem(size);
      }

      switch (type) {
      case BLACKMAN:
         blackman(w, size);
         break;
      case HAMMING:
         hamming(w, size);
         break;
      case HANNING:
         hanning(w, size);
         break;
      case BARTLETT:
         bartlett(w, size);
         break;
      case TRAPEZOID:
         trapezoid(w, size);
         break;
      case RECTANGULAR:
         rectangular(w, size);
         break;
      default:
         fprintf(stderr, "window : Unknown window type %d!\n", (int) type);
         exit(1);
      }

      switch (nflg) {
      case 1:
         for (i = 0, g = 0.0; i < size; i++)
            g += w[i] * w[i];
         g = sqrt(g);
         for (i = 0; i < size; i++)
            w[i] /= g;
         break;
      case 2:
         for (i = 0, g = 0.0; i < size; i++)
            g += w[i];
         for (i = 0; i < size; i++)
            w[i] /= g;
         break;
      case 0:
      default:
         g = 1.0;
      }

      ptype = type;
      psize = size;
      pnflg = nflg;
   }

   for (i = 0; i < size; i++)
      x[i] = x[i] * w[i];

   return (g);
}

static double sgn(const double x)
{
   if (x >= 0)
      return (0.5);
   else
      return (-0.5);
}

double zcross(double *x, const int fl, const int n)
{
   int i;
   double z = 0;

   for (i = 0; i < fl; i++)
      z += fabs(sgn(x[i + 1]) - sgn(x[i]));
   if (n)
      z /= fl;

   return (z);
}

double zerodf(double x, double *b, int m, double *d)
{
   double out;

   out = b[0] * x;

   for (m--; m > 0; m--) {
      out += b[m + 1] * d[m];
      d[m] = d[m - 1];
   }
   out += b[1] * d[0];
   d[0] = x;

   return (out);
}

double zerodft(double x, double *b, const int m, double *d)
{
   int i;
   double out;

   out = b[0] * x + d[0];

   for (i = 1; i < m; i++)
      d[i - 1] = b[i] * x + d[i];

   d[m - 1] = b[m] * x;

   return (out);
}

double zerodf1(double x, double *b, int m, double *d)
{
   double out;

   out = x;
   for (m--; m > 0; m--) {
      out += b[m + 1] * d[m];
      d[m] = d[m - 1];
   }
   out += b[1] * d[0];
   d[0] = x;

   return (out);
}

double zerodf1t(double x, double *b, const int m, double *d)
{
   int i;
   double out;

   out = x + d[0];

   for (i = 1; i < m; i++)
      d[i - 1] = b[i] * x + d[i];

   d[m - 1] = b[m] * x;

   return (out);
}