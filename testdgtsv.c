/* testlapack.c - test program for the LAPACK dgtsv function
 *
 * Copyright (C) 2004  Jochen Voss.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * $Id: testdgtsv.c 5702 2004-05-12 21:13:41Z voss $
 */

#include <stdio.h>

double l[] = {
  -1, -2, -1, -1
};

double d[] = {
  2, 2, 3, 3, 1
};

static long
dstevx(long N, long NRHS, double *DL, double *D, double *DU, double *B,
      long LDB)
{
  extern void dstevx_(const long *Np, const long *NRHSp, double *DL,
		     double *D, double *DU, double *B, const long *LDBp,
		     long *INFOp);
  long info;
  dstevx_(&N, &NRHS, DL, D, DU, B, &LDB, &info);
  return info;
}

int
main()
{
  int i, info;

  info = dstevx(5, 1, l, d, u, x, 5);
  if (info != 0) fprintf(stderr, "failure with error %d\n", info);

  for (i=0; i<5; ++i) printf("%5.1f\n", x[i]);

  return 0;
}
