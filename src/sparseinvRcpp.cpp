/* sparseinv: An R Software package for computing the sparse inverse subset with the Takahashi equations with large datasets.
 Copyright (c) 2018 Andrew Zammit-Mangion
 Author: Andrew Zammit-Mangion, azm (at) uow.edu.au

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.*/

/* This program was adapted from sparseinv.C. The copyright note there was the following.*/

/* sparsinv: computes the sparse inverse subset, using Takahashi's equations.

 On input, the pattern of Z must be equal to the symbolic Cholesky
factorization of A+A', where A=(L+I)*(U+I).  The pattern of L+U must be a
subset of Z.  Z must have zero-free diagonal.  These conditions are
difficult to check, so they are assumed to hold.  Results will be completely
wrong if the conditions do not hold.

This function performs the same amount of work as the initial LU
factorization, assuming that the pattern of P*A*Q is symmetric.  For large
matrices, this function can take a lot more time than LU in MATLAB, even if
P*A*Q is symmetric.  This is because LU is a multifrontal method, whereas
this sparseinv function is based on gather/scatter operations.

The basic integer type is an Int, or ptrdiff_t, which is 32 bits on a 32
bits and 64 bits on a 64 bit system.  The function returns the flop count as
an Int.  This will not overflow on a 64 bit system but might on a 32 bit.
The total work is flops + O(n + nnz(Z)).  Since flops > n and flops > nnz(Z),
this is O(flops).

Copyright 2011, Timothy A. Davis, University of Florida
*/
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sparseinv2(int n, IntegerVector Lp, IntegerVector Li, NumericVector Lx, NumericVector d, IntegerVector Up, IntegerVector Uj, NumericVector Ux, IntegerVector Zp,  IntegerVector Zi) {

    double ljk, zkj ;
    int j, i, k, p, znz, pdiag, up, zp, flops = n ;
    double *z = (double *) calloc(n,sizeof(double));
    int *Zdiagp = (int *) malloc((n)*sizeof(int));
    int *Lmunch = (int *) malloc((n)*sizeof(int));


    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* clear the numerical values of Z */
    int lengthZx = Zp[n] ;
    NumericVector Zx(lengthZx);
    znz = Zp[n] ;
    for (p = 0 ; p < znz ; p++)
    {
        Zx[p] = 0 ;
    }

    /* find the diagonal of Z and initialize it */
    for (j = 0 ; j < n ; j++)
    {
        pdiag = -1 ;
        for (p = Zp [j] ; p < Zp [j+1] && pdiag == -1 ; p++)
        {
            if (Zi [p] == j)
            {
                pdiag = p ;
                Zx [p] = 1 / d [j] ;
            }
        }
        Zdiagp [j] = pdiag ;
        if (pdiag == -1) return (-1) ;  /* Z must have a zero-free diagonal */
    }

    /* Lmunch [k] points to the last entry in column k of L */
    for (k = 0 ; k < n ; k++)
    {
        Lmunch [k] = Lp [k+1] - 1 ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute the sparse inverse subset */
    /* ---------------------------------------------------------------------- */

    for (j = (n)-1 ; j >= 0 ; j--)
    {

        /* ------------------------------------------------------------------ */
        /* scatter Z (:,j) into z workspace */
        /* ------------------------------------------------------------------ */

        /* only the lower triangular part is needed, since the upper triangular
         part is all zero */
        for (p = Zdiagp [j] ; p < Zp [j+1] ; p++)
        {
            z [Zi [p]] = Zx [p] ;
        }

        /* ------------------------------------------------------------------ */
        /* compute the strictly upper triangular part of Z (:,j) */
        /* ------------------------------------------------------------------ */

        /* for k = (j-1):-1:1 but only for the entries Z(k,j) */
        for (p = Zdiagp [j]-1 ; p >= Zp [j] ; p--)
        {
            /* Z (k,j) = - U (k,k+1:n) * Z (k+1:n,j) */
            k = Zi [p] ;
            zkj = 0 ;
            flops += (Up [k+1] - Up [k]) ;
            for (up = Up [k] ; up < Up [k+1] ; up++)
            {
                /* skip the diagonal of U, if present */
                i = Uj [up] ;
                if (i > k)
                {
                    zkj -= Ux [up] * z [i] ;
                }
            }
            z [k] = zkj ;
        }

        /* ------------------------------------------------------------------ */
        /* left-looking update to lower triangular part of Z */
        /* ------------------------------------------------------------------ */

        /* for k = (j-1):-1:1 but only for the entries Z(k,j) */
        for (p = Zdiagp [j]-1 ; p >= Zp [j] ; p--)
        {
            k = Zi [p] ;

            /* ljk = L (j,k) */
            if (Lmunch [k] < Lp [k] || Li [Lmunch [k]] != j)
            {
                /* L (j,k) is zero, so there is no work to do */
                continue ;
            }
            ljk = Lx [Lmunch [k]--] ;

            /* Z (k+1:n,k) = Z (k+1:n,k) - Z (k+1:n,j) * L (j,k) */
            flops += (Zp [k+1] - Zdiagp [k]) ;
            for (zp = Zdiagp [k] ; zp < Zp [k+1] ; zp++)
            {
                Zx [zp] -= z [Zi [zp]] * ljk ;
            }
        }

        /* ------------------------------------------------------------------ */
        /* gather Z (:,j) back from z workspace */
        /* ------------------------------------------------------------------ */

        for (p = Zp [j] ; p < Zp [j+1] ; p++)
        {
            i = Zi [p] ;
            Zx [p] = z [i] ;
            z [i] = 0 ;
        }
    }

    free(z);
    free(Zdiagp);
    free(Lmunch);

    return Zx;
}
