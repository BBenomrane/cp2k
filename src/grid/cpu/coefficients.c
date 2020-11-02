#include <assert.h>
#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__LIBXSMM)
#include <libxsmm.h>
#endif
#include "../common/grid_common.h"
#include "tensor_local.h"
#include "coefficients.h"
#include "private_header.h"

void transform_xyz_to_triangular(const tensor *const coef,
                                 double *const coef_xyz) {
  assert(coef != NULL);
  assert(coef_xyz != NULL);

  int lxyz = 0;
  const int lp = (coef->size[0] - 1);
  for (int lzp = 0; lzp <= lp; lzp++) {
    for (int lyp = 0; lyp <= lp - lzp; lyp++) {
      for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++, lxyz++) {
        coef_xyz[lxyz] = idx3(coef[0], lzp, lyp, lxp);
      }
    }
  }
}

void transform_yxz_to_triangular(const tensor *const coef,
                                 double *const coef_xyz) {
  assert(coef != NULL);
  assert(coef_xyz != NULL);
  int lxyz = 0;
  const int lp = (coef->size[0] - 1);
  for (int lzp = 0; lzp <= lp; lzp++) {
    for (int lyp = 0; lyp <= lp - lzp; lyp++) {
      for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++, lxyz++) {
        coef_xyz[lxyz] = idx3(coef[0], lyp, lxp, lzp);
      }
    }
  }
}

void transform_triangular_to_xyz(const double *const coef_xyz,
                                 tensor *const coef) {
  assert(coef != NULL);
  assert(coef_xyz != NULL);
  int lxyz = 0;
  const int lp = coef->size[0] - 1;
  for (int lzp = 0; lzp <= lp; lzp++) {
    for (int lyp = 0; lyp <= lp - lzp; lyp++) {
      for (int lxp = 0; lxp <= lp - lzp - lyp; lxp++, lxyz++) {
        idx3(coef[0], lzp, lyp, lxp) = coef_xyz[lxyz];
      }
      /* initialize the remaining coefficients to zero */
      for (int lxp = lp - lzp - lyp + 1; lxp <= lp; lxp++)
        idx3(coef[0], lzp, lyp, lxp) = 0.0;
    }
  }
}

// *****************************************************************************
void grid_prepare_coef_dgemm(
    const int *lmin, const int *lmax, const int lp, const double prefactor,
    const tensor *const alpha,  // [3][lb_max+1][la_max+1][lp+1]
    const tensor *const pab,
    tensor *coef_xyz)  //[lp+1][lp+1][lp+1]
{
  /* can be done with dgemms as well, since it is a change of basis from (x -
   * x1) (x - x2) to (x - x12)^alpha */

  assert(alpha != NULL);
  assert(coef_xyz != NULL);
  assert(coef_xyz->data != NULL);
  memset(coef_xyz->data, 0, coef_xyz->alloc_size_ * sizeof(double));
  // we need a proper fix for that. We can use the tensor structure for this

  double *coef_xyt =
      grid_allocate_scratch(sizeof(double) * (lp + 1) * (lp + 1));
  double *coef_xtt = grid_allocate_scratch(sizeof(double) * (lp + 1));

  for (int lzb = 0; lzb <= lmax[1]; lzb++) {
    for (int lza = 0; lza <= lmax[0]; lza++) {
      memset(coef_xyt, 0, sizeof(double) * (lp + 1) * (lp + 1));
      for (int lyb = 0; lyb <= lmax[1] - lzb; lyb++) {
        for (int lya = 0; lya <= lmax[0] - lza; lya++) {
          const int lxpm = (lmax[1] - lzb - lyb) + (lmax[0] - lza - lya);
          for (int i = 0; i <= lxpm; i++) {
            coef_xtt[i] = 0.0;
          }
          for (int lxb = imax(lmin[1] - lzb - lyb, 0);
               lxb <= lmax[1] - lzb - lyb; lxb++) {
            for (int lxa = imax(lmin[0] - lza - lya, 0);
                 lxa <= lmax[0] - lza - lya; lxa++) {
              const int ico = coset(lxa, lya, lza);
              const int jco = coset(lxb, lyb, lzb);
              const double p_ele = prefactor * idx2(pab[0], jco, ico);
              const double *__restrict__ src = &idx4(alpha[0], 0, lxb, lxa, 0);
#pragma GCC ivdep
              for (int lxp = 0; lxp <= lxa + lxb; lxp++) {
                coef_xtt[lxp] += p_ele * src[lxp];
              }
            }
          }
          for (int lyp = 0; lyp <= lya + lyb; lyp++) {
            const double pe = idx4(alpha[0], 1, lyb, lya, lyp);
            double *__restrict__ dst = &coef_xyt[lyp * (lp + 1)];
#pragma GCC ivdep
            for (int lxp = 0; lxp <= lp - lza - lzb - lya - lyb; lxp++) {
              dst[lxp] += pe * coef_xtt[lxp];
            }
          }
        }
      }
      /* I need to permute two of the indices for the orthogonal case */
      for (int lzp = 0; lzp <= lza + lzb; lzp++) {
        for (int lyp = 0; lyp <= lp - lza - lzb; lyp++) {
          const double pe = idx4(alpha[0], 2, lzb, lza, lzp);
          const double *__restrict__ src = &coef_xyt[lyp * (lp + 1)];

#pragma GCC ivdep
          for (int lxp = 0; lxp <= lp - lza - lzb - lyp; lxp++) {
            idx3(coef_xyz[0], lxp, lzp, lyp) += pe * src[lxp];
          }
        }
      }
    }
  }
  grid_free_scratch(coef_xyt);
  grid_free_scratch(coef_xtt);
}

// *****************************************************************************
void grid_prepare_alpha_dgemm(const double ra[3], const double rb[3],
                              const double rp[3], const int *lmax,
                              tensor *alpha) {
  assert(alpha != NULL);
  // Initialize with zeros.
  memset(alpha->data, 0, alpha->alloc_size_ * sizeof(double));

  //
  //   compute polynomial expansion coefs -> (x-a)**lxa (x-b)**lxb -> sum_{ls}
  //   alpha(ls,lxa,lxb,1)*(x-p)**ls
  //

  for (int iaxis = 0; iaxis < 3; iaxis++) {
    const double drpa = rp[iaxis] - ra[iaxis];
    const double drpb = rp[iaxis] - rb[iaxis];
    for (int lxa = 0; lxa <= lmax[0]; lxa++) {
      for (int lxb = 0; lxb <= lmax[1]; lxb++) {
        double binomial_k_lxa = 1.0;
        double a = 1.0;
        for (int k = 0; k <= lxa; k++) {
          double binomial_l_lxb = 1.0;
          double b = 1.0;
          for (int l = 0; l <= lxb; l++) {
            idx4(alpha[0], iaxis, lxb, lxa, lxa - l + lxb - k) +=
                binomial_k_lxa * binomial_l_lxb * a * b;
            binomial_l_lxb *= ((double)(lxb - l)) / ((double)(l + 1));
            b *= drpb;
          }
          binomial_k_lxa *= ((double)(lxa - k)) / ((double)(k + 1));
          a *= drpa;
        }
      }
    }
  }
}

/* this function computes the coefficients initially expressed in the cartesian
 * space to the grid space. It is inplane and can also be done with
 * matrix-matrix multiplication. It is in fact a tensor reduction. */

void grid_transform_coef_xzy_to_ikj(const double dh[3][3],
                                    const tensor *coef_xyz) {
  assert(coef_xyz != NULL);
  const int lp = coef_xyz->size[0] - 1;
  tensor coef_ijk;

  /* this tensor corresponds to the term
   * $v_{11}^{k_{11}}v_{12}^{k_{12}}v_{13}^{k_{13}}
   * v_{21}^{k_{21}}v_{22}^{k_{22}}v_{23}^{k_{23}}
   * v_{31}^{k_{31}}v_{32}^{k_{32}} v_{33}^{k_{33}}$ in Eq.26 found section
   * III.A of the notes */
  tensor hmatgridp;

  initialize_tensor_3(&coef_ijk, coef_xyz->size[0], coef_xyz->size[1],
                      coef_xyz->size[2]);

  coef_ijk.data = grid_allocate_scratch(sizeof(double) * coef_ijk.alloc_size_);

  if (coef_ijk.data == NULL) abort();

  memset(coef_ijk.data, 0, sizeof(double) * coef_ijk.alloc_size_);
  initialize_tensor_3(&hmatgridp, coef_xyz->size[0], 3, 3);

  hmatgridp.data =
      grid_allocate_scratch(sizeof(double) * hmatgridp.alloc_size_);

  // transform using multinomials
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      idx3(hmatgridp, 0, j, i) = 1.0;
      for (int k = 1; k <= lp; k++) {
        idx3(hmatgridp, k, j, i) = idx3(hmatgridp, k - 1, j, i) * dh[j][i];
      }
    }
  }

  const int lpx = lp;
  for (int klx = 0; klx <= lpx; klx++) {
    for (int jlx = 0; jlx <= lpx - klx; jlx++) {
      for (int ilx = 0; ilx <= lpx - klx - jlx; ilx++) {
        const int lx = ilx + jlx + klx;
        const int lpy = lp - lx;
        const double tx = idx3(hmatgridp, ilx, 0, 0) *
                          idx3(hmatgridp, jlx, 1, 0) *
                          idx3(hmatgridp, klx, 2, 0) * fac[lx] * inv_fac[klx] *
                          inv_fac[jlx] * inv_fac[ilx];

        for (int kly = 0; kly <= lpy; kly++) {
          for (int jly = 0; jly <= lpy - kly; jly++) {
            for (int ily = 0; ily <= lpy - kly - jly; ily++) {
              const int ly = ily + jly + kly;
              const int lpz = lp - lx - ly;
              const double ty = tx * idx3(hmatgridp, ily, 0, 1) *
                                idx3(hmatgridp, jly, 1, 1) *
                                idx3(hmatgridp, kly, 2, 1) * fac[ly] *
                                inv_fac[kly] * inv_fac[jly] * inv_fac[ily];
              for (int klz = 0; klz <= lpz; klz++) {
                for (int jlz = 0; jlz <= lpz - klz; jlz++) {
                  for (int ilz = 0; ilz <= lpz - klz - jlz; ilz++) {
                    const int lz = ilz + jlz + klz;
                    const int il = ilx + ily + ilz;
                    const int jl = jlx + jly + jlz;
                    const int kl = klx + kly + klz;
                    // const int lijk= coef_map[kl][jl][il];
                    /* the fac table is the factorial. It
                     * would be better to use the
                     * multinomials. */
                    idx3(coef_ijk, il, kl, jl) +=
                        idx3(coef_xyz[0], lx, lz, ly) * ty *
                        idx3(hmatgridp, ilz, 0, 2) *
                        idx3(hmatgridp, jlz, 1, 2) *
                        idx3(hmatgridp, klz, 2, 2) * fac[lz] * inv_fac[klz] *
                        inv_fac[jlz] * inv_fac[ilz];
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  memcpy(coef_xyz->data, coef_ijk.data, sizeof(double) * coef_ijk.alloc_size_);
  grid_free_scratch(coef_ijk.data);
  grid_free_scratch(hmatgridp.data);
}

/* Rotate the coefficients computed in the local grid coordinates to the
 * cartesians coorinates. The order of the indices indicates how the
 * coefficients are stored */
void grid_transform_coef_jik_to_yxz(const double dh[3][3],
                                    const tensor *coef_xyz) {
  assert(coef_xyz);
  const int lp = coef_xyz->size[0] - 1;
  tensor coef_ijk;

  /* this tensor corresponds to the term
   * $v_{11}^{k_{11}}v_{12}^{k_{12}}v_{13}^{k_{13}}
   * v_{21}^{k_{21}}v_{22}^{k_{22}}v_{23}^{k_{23}}
   * v_{31}^{k_{31}}v_{32}^{k_{32}} v_{33}^{k_{33}}$ in Eq.26 found section
   * III.A of the notes */
  tensor hmatgridp;

  initialize_tensor_3(&coef_ijk, coef_xyz->size[0], coef_xyz->size[1],
                      coef_xyz->size[2]);

  coef_ijk.data = grid_allocate_scratch(sizeof(double) * coef_ijk.alloc_size_);
  if (coef_ijk.data == NULL) abort();

  memset(coef_ijk.data, 0, sizeof(double) * coef_ijk.alloc_size_);
  initialize_tensor_3(&hmatgridp, coef_xyz->size[0], 3, 3);

  hmatgridp.data =
      grid_allocate_scratch(sizeof(double) * hmatgridp.alloc_size_);

  // transform using multinomials
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      idx3(hmatgridp, 0, j, i) = 1.0;
      for (int k = 1; k <= lp; k++) {
        idx3(hmatgridp, k, j, i) = idx3(hmatgridp, k - 1, j, i) * dh[j][i];
      }
    }
  }

  const int lpx = lp;
  for (int klx = 0; klx <= lpx; klx++) {
    for (int jlx = 0; jlx <= lpx - klx; jlx++) {
      for (int ilx = 0; ilx <= lpx - klx - jlx; ilx++) {
        const int lx = ilx + jlx + klx;
        const int lpy = lp - lx;
        for (int kly = 0; kly <= lpy; kly++) {
          for (int jly = 0; jly <= lpy - kly; jly++) {
            for (int ily = 0; ily <= lpy - kly - jly; ily++) {
              const int ly = ily + jly + kly;
              const int lpz = lp - lx - ly;
              for (int klz = 0; klz <= lpz; klz++) {
                for (int jlz = 0; jlz <= lpz - klz; jlz++) {
                  for (int ilz = 0; ilz <= lpz - klz - jlz; ilz++) {
                    const int lz = ilz + jlz + klz;
                    const int il = ilx + ily + ilz;
                    const int jl = jlx + jly + jlz;
                    const int kl = klx + kly + klz;
                    // const int lijk= coef_map[kl][jl][il];
                    /* the fac table is the factorial. It
                     * would be better to use the
                     * multinomials. */
                    idx3(coef_ijk, ly, lx, lz) +=
                        idx3(coef_xyz[0], jl, il, kl) *
                        idx3(hmatgridp, ilx, 0, 0) *
                        idx3(hmatgridp, jlx, 1, 0) *
                        idx3(hmatgridp, klx, 2, 0) *
                        idx3(hmatgridp, ily, 0, 1) *
                        idx3(hmatgridp, jly, 1, 1) *
                        idx3(hmatgridp, kly, 2, 1) *
                        idx3(hmatgridp, ilz, 0, 2) *
                        idx3(hmatgridp, jlz, 1, 2) *
                        idx3(hmatgridp, klz, 2, 2) * fac[lx] * fac[ly] *
                        fac[lz] /
                        (fac[ilx] * fac[ily] * fac[ilz] * fac[jlx] * fac[jly] *
                         fac[jlz] * fac[klx] * fac[kly] * fac[klz]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  memcpy(coef_xyz->data, coef_ijk.data, sizeof(double) * coef_ijk.alloc_size_);
  grid_free_scratch(coef_ijk.data);
  grid_free_scratch(hmatgridp.data);
}
