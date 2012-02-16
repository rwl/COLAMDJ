/**
 * Copyright (c) 1998-2007, Timothy A. Davis, All Rights Reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
 */

package edu.ufl.cise.colamd.tdouble.test;

import junit.framework.TestCase;

import static edu.ufl.cise.colamd.tdouble.Dcolamd.COLAMD_STATS;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.colamd;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.COLAMD_report;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.symamd;
import static edu.ufl.cise.colamd.tdouble.Dcolamd.SYMAMD_report;

/**
 * COLAMD / SYMAMD example
 *
 * colamd example of use, to order the columns of a 5-by-4 matrix with
 * 11 nonzero entries in the following nonzero pattern, with default knobs.
 *
 *   x 0 x 0
 *   x 0 x x
 *   0 x x 0
 *   0 0 x x
 *   x x 0 0
 *
 * symamd example of use, to order the rows and columns of a 5-by-5
 * matrix with 13 nonzero entries in the following nonzero pattern,
 * with default knobs.
 *
 *   x x 0 0 0
 *   x x x x 0
 *   0 x x 0 0
 *   0 x 0 x x
 *   0 0 0 x x
 *
 * (where x denotes a nonzero value).
*/
public class Dcolamd_example extends TestCase {

	int A_NNZ = 11 ;
	int A_NROW = 5 ;
	int A_NCOL = 4 ;
	int ALEN = 150 ;

	int B_NNZ = 4 ;
	int B_N = 5 ;

	public void test_colamd() {

		/* ====================================================================== */
		/* input matrix A definition */
		/* ====================================================================== */

		int[] A = new int [ALEN] ;
		int[] AA = new int [] {

				0, 1, 4,	/* row indices of nonzeros in column 0 */
				2, 4,		/* row indices of nonzeros in column 1 */
				0, 1, 2, 3,	/* row indices of nonzeros in column 2 */
				1, 3} ;		/* row indices of nonzeros in column 3 */
		System.arraycopy(AA, 0, A, 0, AA.length) ;

		int[] p = new int [] {

				0,		/* column 0 is in A [0..2] */
				3,		/* column 1 is in A [3..4] */
				5,		/* column 2 is in A [5..8] */
				9,		/* column 3 is in A [9..10] */
				A_NNZ} ;	/* number of nonzeros in A */

		/* ====================================================================== */
		/* input matrix B definition */
		/* ====================================================================== */

		int[] B = new int [] {		/* Note: only strictly lower triangular part */
		    				/* is included, since symamd ignores the */
						/* diagonal and upper triangular part of B. */

				1,		/* row indices of nonzeros in column 0 */
				2, 3,		/* row indices of nonzeros in column 1 */
						/* row indices of nonzeros in column 2 (none) */
				4		/* row indices of nonzeros in column 3 */
		} ;				/* row indices of nonzeros in column 4 (none) */

		int[] q = new int[] {

				0,		/* column 0 is in B [0] */
				1,		/* column 1 is in B [1..2] */
				3,		/* column 2 is empty */
				3,		/* column 3 is in B [3] */
				4,		/* column 4 is empty */
				B_NNZ} ;	/* number of nonzeros in strictly lower B */

		/* ====================================================================== */
		/* other variable definitions */
		/* ====================================================================== */

		int[] perm = new int [B_N+1] ;		/* note the size is N+1 */
		int[] stats = new int [COLAMD_STATS] ;	/* for colamd and symamd output statistics */

		int row, col, pp, length, ok ;

		/* ====================================================================== */
		/* dump the input matrix A */
		/* ====================================================================== */

		System.out.printf ("colamd %d-by-%d input matrix:\n", A_NROW, A_NCOL) ;
		for (col = 0 ; col < A_NCOL ; col++)
		{
			length = p [col+1] - p [col] ;
		    	System.out.printf ("Column %d, with %d entries:\n", col, length) ;
			for (pp = p [col] ; pp < p [col+1] ; pp++)
			{
				row = A [pp] ;
				System.out.printf ("    row %d\n", row) ;
			}
		}

		/* ====================================================================== */
		/* order the matrix.  Note that this destroys A and overwrites p */
		/* ====================================================================== */

		ok = colamd (A_NROW, A_NCOL, ALEN, A, p, null, stats) ;
		COLAMD_report (stats) ;

		if (ok == 0)
		{
			System.out.printf ("colamd error!\n") ;
			fail () ;
		}

		/* ====================================================================== */
		/* print the column ordering */
		/* ====================================================================== */

		System.out.printf ("colamd column ordering:\n") ;
		System.out.printf ("1st column: %d\n", p [0]) ;
		System.out.printf ("2nd column: %d\n", p [1]) ;
		System.out.printf ("3rd column: %d\n", p [2]) ;
		System.out.printf ("4th column: %d\n", p [3]) ;

		assertEquals(1, p [0]) ;
		assertEquals(0, p [1]) ;
		assertEquals(2, p [2]) ;
		assertEquals(3, p [3]) ;

		/* ====================================================================== */
		/* dump the strictly lower triangular part of symmetric input matrix B */
		/* ====================================================================== */

		System.out.printf ("\n\nsymamd %d-by-%d input matrix:\n", B_N, B_N) ;
		System.out.printf ("Entries in strictly lower triangular part:\n") ;
		for (col = 0 ; col < B_N ; col++)
		{
			length = q [col+1] - q [col] ;
		    	System.out.printf ("Column %d, with %d entries:\n", col, length) ;
			for (pp = q [col] ; pp < q [col+1] ; pp++)
			{
				row = B [pp] ;
				System.out.printf ("    row %d\n", row) ;
			}
		}

		/* ====================================================================== */
		/* order the matrix B.  Note that this does not modify B or q. */
		/* ====================================================================== */

		ok = symamd (B_N, B, q, perm, null, stats) ;
		SYMAMD_report (stats) ;

		if (ok == 0)
		{
			System.out.printf ("symamd error!\n") ;
			fail () ;
		}

		/* ====================================================================== */
		/* print the symmetric ordering */
		/* ====================================================================== */

		System.out.printf ("symamd column ordering:\n") ;
		System.out.printf ("1st row/column: %d\n", perm [0]) ;
		System.out.printf ("2nd row/column: %d\n", perm [1]) ;
		System.out.printf ("3rd row/column: %d\n", perm [2]) ;
		System.out.printf ("4th row/column: %d\n", perm [3]) ;
		System.out.printf ("5th row/column: %d\n", perm [4]) ;

		assertEquals(0, perm [0]) ;
		assertEquals(2, perm [1]) ;
		assertEquals(1, perm [2]) ;
		assertEquals(3, perm [3]) ;
		assertEquals(4, perm [4]) ;
	}

}
