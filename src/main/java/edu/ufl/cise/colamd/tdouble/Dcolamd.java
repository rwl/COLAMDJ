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

package edu.ufl.cise.colamd.tdouble;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * COLAMD / SYMAMD
 *
 * colamd:  an approximate minimum degree column ordering algorithm,
 * for LU factorization of symmetric or unsymmetric matrices,
 * QR factorization, least squares, interior point methods for
 * linear programming problems, and other related problems.
 *
 * symamd:  an approximate minimum degree ordering algorithm for Cholesky
 * factorization of symmetric matrices.
 *
 * Purpose:
 *
 * Colamd computes a permutation Q such that the Cholesky factorization of
 * (AQ)'(AQ) has less fill-in and requires fewer floating point operations
 * than A'A.  This also provides a good ordering for sparse partial
 * pivoting methods, P(AQ) = LU, where Q is computed prior to numerical
 * factorization, and P is computed during numerical factorization via
 * conventional partial pivoting with row interchanges.  Colamd is the
 * column ordering method used in SuperLU, part of the ScaLAPACK library.
 * It is also available as built-in function in MATLAB Version 6,
 * available from MathWorks, Inc. (http://www.mathworks.com).  This
 * routine can be used in place of colmmd in MATLAB.
 *
 * Symamd computes a permutation P of a symmetric matrix A such that the
 * Cholesky factorization of PAP' has less fill-in and requires fewer
 * floating point operations than A.  Symamd constructs a matrix M such
 * that M'M has the same nonzero pattern of A, and then orders the columns
 * of M using colmmd.  The column ordering of M is then returned as the
 * row and column ordering P of A.
 *
 * Authors:
 *
 * The authors of the code itself are Stefan I. Larimore and Timothy A.
 * Davis (davis at cise.ufl.edu), University of Florida.  The algorithm was
 * developed in collaboration with John Gilbert, Xerox PARC, and Esmond
 * Ng, Oak Ridge National Laboratory.
 *
 * Acknowledgements:
 *
 * This work was supported by the National Science Foundation, under
 * grants DMS-9504974 and DMS-9803599.
 */
public class Dcolamd {

/* ========================================================================== */
/* === COLAMD version ======================================================= */
/* ========================================================================== */

	/* COLAMD Version 2.4 and later will include the following definitions.
	 * As an example, to test if the version you are using is 2.4 or later:
	 *
	 *	if (COLAMD_VERSION >= COLAMD_VERSION_CODE (2,4)) ...
	 */

	public static String COLAMD_DATE = "Jan 25, 2011" ;
	public static int COLAMD_VERSION_CODE(int main, int sub)
	{
		return ((main) * 1000 + (sub)) ;
	}
	public static int COLAMD_MAIN_VERSION = 2 ;
	public static int COLAMD_SUB_VERSION = 7 ;
	public static int COLAMD_SUBSUB_VERSION = 3 ;
	public static int COLAMD_VERSION = COLAMD_VERSION_CODE(COLAMD_MAIN_VERSION,
			COLAMD_SUB_VERSION) ;

/* ========================================================================== */
/* === Knob and statistics definitions ====================================== */
/* ========================================================================== */

	/* size of the knobs [ ] array.  Only knobs [0..1] are currently used. */
	public static int COLAMD_KNOBS = 20 ;

	/* number of output statistics.  Only stats [0..6] are currently used. */
	public static int COLAMD_STATS = 20 ;

	/* knobs [0] and stats [0]: dense row knob and output statistic. */
	public static int COLAMD_DENSE_ROW = 0 ;

	/* knobs [1] and stats [1]: dense column knob and output statistic. */
	public static int COLAMD_DENSE_COL = 1 ;

	/* knobs [2]: aggressive absorption */
	public static int COLAMD_AGGRESSIVE = 2 ;

	/* stats [2]: memory defragmentation count output statistic */
	public static int COLAMD_DEFRAG_COUNT = 2 ;

	/* stats [3]: colamd status:  zero OK, > 0 warning or notice, < 0 error */
	public static int COLAMD_STATUS = 3 ;

	/* stats [4..6]: error info, or info on jumbled columns */
	public static final int COLAMD_INFO1 = 4 ;
	public static final int COLAMD_INFO2 = 5 ;
	public static final int COLAMD_INFO3 = 6 ;

	/* error codes returned in stats [3]: */
	public static final int COLAMD_OK				= (0) ;
	public static final int COLAMD_OK_BUT_JUMBLED			= (1) ;
	public static final int COLAMD_ERROR_A_not_present		= (-1) ;
	public static final int COLAMD_ERROR_p_not_present		= (-2) ;
	public static final int COLAMD_ERROR_nrow_negative		= (-3) ;
	public static final int COLAMD_ERROR_ncol_negative		= (-4) ;
	public static final int COLAMD_ERROR_nnz_negative		= (-5) ;
	public static final int COLAMD_ERROR_p0_nonzero			= (-6) ;
	public static final int COLAMD_ERROR_A_too_small		= (-7) ;
	public static final int COLAMD_ERROR_col_length_negative	= (-8) ;
	public static final int COLAMD_ERROR_row_index_out_of_bounds	= (-9) ;
	public static final int COLAMD_ERROR_out_of_memory		= (-10) ;
	public static final int COLAMD_ERROR_internal_error		= (-999) ;

/* ========================================================================== */
/* === Scaffolding code definitions  ======================================== */
/* ========================================================================== */

	public static boolean NDEBUG = true;

	public static boolean NPRINT = true;

	/*
	   Our "scaffolding code" philosophy:  In our opinion, well-written library
	   code should keep its "debugging" code, and just normally have it turned off
	   by the compiler so as not to interfere with performance.  This serves
	   several purposes:

	   (1) assertions act as comments to the reader, telling you what the code
		expects at that point.  All assertions will always be true (unless
		there really is a bug, of course).

	   (2) leaving in the scaffolding code assists anyone who would like to modify
		the code, or understand the algorithm (by reading the debugging output,
		one can get a glimpse into what the code is doing).

	   (3) (gasp!) for actually finding bugs.  This code has been heavily tested
		and "should" be fully functional and bug-free ... but you never know...

	    The code will become outrageously slow when debugging is
	    enabled.  To control the level of debugging output, set an environment
	    variable D to 0 (little), 1 (some), 2, 3, or 4 (lots).  When debugging,
	    you should see the following message on the standard output:

	    	colamd: debug version, D = 1 (THIS WILL BE SLOW!)

	    or a similar message for symamd.  If you don't, then debugging has not
	    been enabled.

	*/

	private static int Int_MAX = Integer.MAX_VALUE;

/* ========================================================================== */
/* === Definitions ========================================================== */
/* ========================================================================== */

	protected static double sqrt (double a)
	{
		return Math.sqrt (a) ;
	}

	private static final int MAX (int a, int b)
	{
		return (((a) > (b)) ? (a) : (b)) ;

	}

	private static final int MIN (int a, int b)
	{
		return (((a) < (b)) ? (a) : (b)) ;
	}

	private static final double MAX (double a, double b)
	{
		return (((a) > (b)) ? (a) : (b)) ;

	}

	private static int DENSE_DEGREE (double alpha, int n)
	{
		return ((int) MAX (16.0, (alpha) * sqrt ((double) (n)))) ;
	}

	private static int ONES_COMPLEMENT (int r)
	{
		return (-(r)-1) ;
	}

	private static final int TRUE  = (1) ;
	private static final int FALSE = (0) ;

	private static final int EMPTY = (-1) ;

	/* Row and column status */
	private static final int ALIVE = (0) ;
	private static final int DEAD = (-1) ;

	/* Column status */
	private static final int DEAD_PRINCIPAL = (-1) ;
	private static final int DEAD_NON_PRINCIPAL	= (-2) ;

	/* Row and column status update and checking. */
	private static boolean ROW_IS_DEAD(Colamd_Row[] Row, int r)
	{
		return ROW_IS_MARKED_DEAD (Row [r].mark()) ;
	}
	private static boolean ROW_IS_MARKED_DEAD(int row_mark)
	{
		return (row_mark < ALIVE) ;
	}
	private static boolean ROW_IS_ALIVE(Colamd_Row[] Row, int r)
	{
		return (Row [r].mark() >= ALIVE) ;
	}
	private static boolean COL_IS_DEAD(Colamd_Col[] Col, int c)
	{
		return (Col [c].start < ALIVE) ;
	}
	private static boolean COL_IS_ALIVE(Colamd_Col[] Col, int c)
	{
		return (Col [c].start >= ALIVE) ;
	}
	private static boolean COL_IS_DEAD_PRINCIPAL(Colamd_Col[] Col, int c)
	{
		return (Col [c].start == DEAD_PRINCIPAL) ;
	}
	private static void KILL_ROW(Colamd_Row[] Row, int r)
	{
		Row [r].mark(DEAD) ;
	}
	private static void KILL_PRINCIPAL_COL(Colamd_Col[] Col, int c)
	{
		Col [c].start = DEAD_PRINCIPAL ;
	}
	private static void KILL_NON_PRINCIPAL_COL(Colamd_Col[] Col, int c)
	{
		Col [c].start = DEAD_NON_PRINCIPAL ;
	}

/* ========================================================================== */
/* === Colamd reporting mechanism =========================================== */
/* ========================================================================== */

	/**
	 * In Java, matrices are 0-based and indices are reported as such
	 * in *_report.
	 */
	private static int INDEX (int i)
	{
		return (i) ;
	}

	/**
	 * All output goes through PRINTF.
	 */
	private static void PRINTF (String format, Object... args)
	{
		if (!NPRINT)
		{
			System.out.printf (format, args) ;
		}
	}

	/** debug print level */
	public static int colamd_debug = 0 ;

	protected static void DEBUG0 (String format, Object... args)
	{
		if (!NDEBUG)
		{
			PRINTF (format, args) ;
		}
	}

	protected static void DEBUG1(String format, Object... args)
	{
		if (!NDEBUG)
		{
			if (colamd_debug >= 1) PRINTF (format, args) ;
		}
	}

	protected static void DEBUG2(String format, Object... args)
	{
		if (!NDEBUG)
		{
			if (colamd_debug >= 2) PRINTF (format, args) ;
		}
	}

	protected static void DEBUG3(String format, Object... args)
	{
		if (!NDEBUG)
		{
			if (colamd_debug >= 3) PRINTF (format, args) ;
		}
	}

	protected static void DEBUG4(String format, Object... args)
	{
		if (!NDEBUG)
		{
			if (colamd_debug >= 4) PRINTF (format, args) ;
		}
	}

	protected static void ASSERT (boolean a)
	{
		if (!NDEBUG)
		{
			assert a ;
		}
	}

	protected static void ASSERT (int a)
	{
		ASSERT (a != 0) ;
	}

/* ========================================================================== */
/* === USER-CALLABLE ROUTINES: ============================================== */
/* ========================================================================== */

/* ========================================================================== */
/* === colamd_recommended =================================================== */
/* ========================================================================== */

	/*
	    The colamd_recommended routine returns the suggested size for Alen.  This
	    value has been determined to provide good balance between the number of
	    garbage collections and the memory requirements for colamd.  If any
	    argument is negative, or if integer overflow occurs, a 0 is returned as an
	    error condition.  2*nnz space is required for the row and column
	    indices of the matrix. COLAMD_C (n_col) + COLAMD_R (n_row) space is
	    required for the Col and Row arrays, respectively, which are internal to
	    colamd (roughly 6*n_col + 4*n_row).  An additional n_col space is the
	    minimal amount of "elbow room", and nnz/5 more space is recommended for
	    run time efficiency.

	    Alen is approximately 2.2*nnz + 7*n_col + 4*n_row + 10.

	    This function is not needed when using symamd.
	*/

	/**
	 * add two values of type int, and check for integer overflow
	 */
	private static int t_add (int a, int b, int[] ok)
	{
	    (ok[0]) = (ok[0] != 0) && ((a + b) >= MAX (a,b)) ? 1 : 0;
	    return ((ok[0] != 0) ? (a + b) : 0) ;
	}

	/**
	 * compute a*k where k is a small integer, and check for integer overflow
	 */
	static int t_mult (int a, int k, int[] ok)
	{
	    int i, s = 0 ;
	    for (i = 0 ; i < k ; i++)
	    {
		s = t_add (s, a, ok) ;
	    }
	    return (s) ;
	}

	/**
	 * size of the Col and Row structures
	 */
	private static int COLAMD_C(int n_col, int[] ok)
	{
//		return ((t_mult (t_add (n_col, 1, ok), sizeof (Colamd_Col), ok) / sizeof (int))) ;
		return t_add (n_col, 1, ok) ;
	}

	private static int COLAMD_R(int n_row, int[] ok)
	{
//	    return ((t_mult (t_add (n_row, 1, ok), sizeof (Colamd_Row), ok) / sizeof (int))) ;
	    return t_add (n_row, 1, ok) ;
	}

	/**
	 *
	 * @param nnz number of nonzeros in A. This must
	 * be the same value as p [n_col] in the call to
	 * colamd - otherwise you will get a wrong value
	 * of the recommended memory to use.
	 * @param n_row number of rows in A
	 * @param n_col number of columns in A
	 * @return recommended value of Alen. Returns 0
	 * if any input argument is negative.  The use of this routine
	 * is optional.  Not needed for symamd, which dynamically allocates
	 * its own memory.
	 */
	public static int COLAMD_recommended (int nnz, int n_row, int n_col)
	{
		int s ;
		int[] ok = new int [] { TRUE } ;
		if (nnz < 0 || n_row < 0 || n_col < 0)
		{
			return (0) ;
		}
		s = t_mult (nnz, 2, ok) ;	    /* 2*nnz */
//		c = COLAMD_C (n_col, ok) ;	    /* size of column structures */
//		r = COLAMD_R (n_row, ok) ;	    /* size of row structures */
//		s = t_add (s, c, ok) ;
//		s = t_add (s, r, ok) ;
		s = t_add (s, n_col, ok) ;	    /* elbow room */
		s = t_add (s, nnz/5, ok) ;	    /* elbow room */
		ok[0] = (s < Int_MAX) ? 1 : 0;
		return (ok[0] != 0 ? s : 0) ;
	}


/* ========================================================================== */
/* === colamd_set_defaults ================================================== */
/* ========================================================================== */

	/*
	    The colamd_set_defaults routine sets the default values of the user-
	    controllable parameters for colamd and symamd:

		Colamd: rows with more than max (16, knobs [0] * sqrt (n_col))
		entries are removed prior to ordering.  Columns with more than
		max (16, knobs [1] * sqrt (MIN (n_row,n_col))) entries are removed
		prior to ordering, and placed last in the output column ordering.

		Symamd: Rows and columns with more than max (16, knobs [0] * sqrt (n))
		entries are removed prior to ordering, and placed last in the
		output ordering.

		knobs [0]	dense row control

		knobs [1]	dense column control

		knobs [2]	if nonzero, do aggresive absorption

		knobs [3..19]	unused, but future versions might use this

	*/

	/**
	 * knobs [0] and knobs [1] control dense row and col detection:
	 *
	 * Colamd: rows with more than
	 * max (16, knobs [COLAMD_DENSE_ROW] * sqrt (n_col))
	 * entries are removed prior to ordering.  Columns with more than
	 * max (16, knobs [COLAMD_DENSE_COL] * sqrt (MIN (n_row,n_col)))
	 * entries are removed prior to
	 * ordering, and placed last in the output column ordering.
	 *
	 * Symamd: uses only knobs [COLAMD_DENSE_ROW], which is knobs [0].
	 * Rows and columns with more than
	 * max (16, knobs [COLAMD_DENSE_ROW] * sqrt (n))
	 * entries are removed prior to ordering, and placed last in the
	 * output ordering.
	 *
	 * COLAMD_DENSE_ROW and COLAMD_DENSE_COL are defined as 0 and 1,
	 * respectively, in colamd.h.  Default values of these two knobs
	 * are both 10.  Currently, only knobs [0] and knobs [1] are
	 * used, but future versions may use more knobs.  If so, they will
	 * be properly set to their defaults by the future version of
	 * colamd_set_defaults, so that the code that calls colamd will
	 * not need to change, assuming that you either use
	 * colamd_set_defaults, or pass a (double *) NULL pointer as the
	 * knobs array to colamd or symamd.
	 *
	 * knobs [2]: aggressive absorption
	 *
	 * knobs [COLAMD_AGGRESSIVE] controls whether or not to do
	 * aggressive absorption during the ordering.  Default is TRUE.
	 *
	 * @param knobs knob array
	 */
	public static void COLAMD_set_defaults (double[] knobs)
	{
		/* === Local variables ============================================== */

		int i ;

		if (knobs == null || knobs.length == 0)
		{
			return ;			/* no knobs to initialize */
		}
		for (i = 0 ; i < COLAMD_KNOBS ; i++)
		{
			knobs [i] = 0 ;
		}
		knobs [COLAMD_DENSE_ROW] = 10 ;
		knobs [COLAMD_DENSE_COL] = 10 ;
		knobs [COLAMD_AGGRESSIVE] = TRUE ;	/* default: do aggressive absorption*/
	}

/* ========================================================================== */
/* === symamd =============================================================== */
/* ========================================================================== */

	/**
	 * The symamd routine computes an ordering P of a symmetric sparse
	 * matrix A such that the Cholesky factorization PAP' = LL' remains
	 * sparse.  It is based on a column ordering of a matrix M constructed
	 * so that the nonzero pattern of M'M is the same as A.  The matrix A
	 * is assumed to be symmetric; only the strictly lower triangular part
	 * is accessed.  You must pass your selected memory allocator (usually
	 * calloc/free or mxCalloc/mxFree) to symamd, for it to allocate
	 * memory for the temporary matrix M.
	 *
	 * @param n number of rows and columns of A. Restriction:  n >= 0.
	 * Symamd returns FALSE if n is negative.
	 * @param A an integer array of size nnz, where nnz = p [n].
	 *
	 * The row indices of the entries in column c of the matrix are
	 * held in A [(p [c]) ... (p [c+1]-1)].  The row indices in a
	 * given column c need not be in ascending order, and duplicate
	 * row indices may be present.  However, symamd will run faster
	 * if the columns are in sorted order with no duplicate entries.
	 *
	 * The matrix is 0-based.  That is, rows are in the range 0 to
	 * n-1, and columns are in the range 0 to n-1.  Symamd
	 * returns FALSE if any row index is out of range.
	 *
	 * The contents of A are not modified.
	 * @param an integer array of size n+1.  On input, it holds the
	 * "pointers" for the column form of the matrix A.  Column c of
	 * the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
	 * entry, p [0], must be zero, and p [c] <= p [c+1] must hold
	 * for all c in the range 0 to n-1.  The value p [n] is
	 * thus the total number of entries in the pattern of the matrix A.
	 * Symamd returns FALSE if these conditions are not met.
	 *
	 * The contents of p are not modified.
	 * @param perm On output, if symamd returns TRUE, the array perm holds the
	 * permutation P, where perm [0] is the first index in the new
	 * ordering, and perm [n-1] is the last.  That is, perm [k] = j
	 * means that row and column j of A is the kth column in PAP',
	 * where k is in the range 0 to n-1 (perm [0] = j means
	 * that row and column j of A are the first row and column in
	 * PAP').  The array is used as a workspace during the ordering,
	 * which is why it must be of length n+1, not just n.
	 * @param knobs parameters (uses defaults if NULL)
	 * @param stats Statistics on the ordering, and error status.
	 * Symamd returns FALSE if stats is not present.
	 *
	 * stats [0]:  number of dense or empty row and columns ignored
	 * 		(and ordered last in the output permutation
	 * 		perm).  Note that a row/column can become
	 * 		"empty" if it contains only "dense" and/or
	 * 		"empty" columns/rows.
	 *
	 * stats [1]:  (same as stats [0])
	 *
	 * stats [2]:  number of garbage collections performed.
	 *
	 * stats [3]:  status code.  < 0 is an error code.
	 * 	    > 1 is a warning or notice.
	 *
	 * 		0	OK.  Each column of the input matrix contained
	 * 			row indices in increasing order, with no
	 * 			duplicates.
	 *
	 * 		1	OK, but columns of input matrix were jumbled
	 * 			(unsorted columns or duplicate entries).  Symamd
	 * 			had to do some extra work to sort the matrix
	 * 			first and remove duplicate entries, but it
	 * 			still was able to return a valid permutation
	 * 			(return value of symamd was TRUE).
	 *
	 * 				stats [4]: highest numbered column that
	 * 					is unsorted or has duplicate
	 * 					entries.
	 * 				stats [5]: last seen duplicate or
	 * 					unsorted row index.
	 * 				stats [6]: number of duplicate or
	 * 					unsorted row indices.
	 *
	 * 		-1	A is a null pointer
	 *
	 * 		-2	p is a null pointer
	 *
	 * 		-3	(unused, see colamd.c)
	 *
	 * 		-4 	n is negative
	 *
	 * 				stats [4]: n
	 *
	 * 		-5	number of nonzeros in matrix is negative
	 *
	 * 				stats [4]: # of nonzeros (p [n]).
	 *
	 * 		-6	p [0] is nonzero
	 *
	 * 				stats [4]: p [0]
	 *
	 * 		-7	(unused)
	 *
	 * 		-8	a column has a negative number of entries
	 *
	 * 				stats [4]: column with < 0 entries
	 * 				stats [5]: number of entries in col
	 *
	 * 		-9	a row index is out of bounds
	 *
	 * 				stats [4]: column with bad row index
	 * 				stats [5]: bad row index
	 * 				stats [6]: n_row, # of rows of matrx
	 *
	 * 		-10	out of memory (unable to allocate temporary
	 * 			workspace for M or count arrays using the
	 * 			"allocate" routine passed into symamd).
	 *
	 * Future versions may return more statistics in the stats array.
	 * @param allocate pointer to calloc
	 * @param release pointer to free
	 * @return TRUE if OK, FALSE otherwise
	 */
	public static int symamd (int n, int[] A, int[] p, int[] perm,
			double[] knobs, int[] stats)
	{
		/* === Local variables ============================================== */

		int[] count ;   /* length of each column of M, and col pointer*/
		int[] mark ;    /* mark array for finding duplicate entries */
		int[] M ;       /* row indices of matrix M */
		int Mlen ;      /* length of M */
		int n_row ;     /* number of rows in M */
		int nnz ;       /* number of entries in A */
		int i ;         /* row index of A */
		int j ;         /* column index of A */
		int k ;         /* row index of M */
		int mnz ;       /* number of nonzeros in M */
		int pp ;        /* index into a column of A */
		int last_row ;  /* last row seen in the current column */
		int length ;    /* number of nonzeros in a column */

		double[] cknobs = new double[COLAMD_KNOBS] ;         /* knobs for colamd */
		double[] default_knobs = new double[COLAMD_KNOBS] ;  /* default knobs for colamd */

		if (!NDEBUG)
		{
//			colamd_get_debug ("symamd") ;
		}

		/* === Check the input arguments ==================================== */

		if (stats == null)
		{
			DEBUG0 ("symamd: stats not present\n") ;
			return (FALSE) ;
		}
		for (i = 0 ; i < COLAMD_STATS ; i++)
		{
			stats [i] = 0 ;
		}
		stats [COLAMD_STATUS] = COLAMD_OK ;
		stats [COLAMD_INFO1] = -1 ;
		stats [COLAMD_INFO2] = -1 ;

		if (A == null)
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
			DEBUG0 ("symamd: A not present\n") ;
			return (FALSE) ;
		}

		if (p == null)		/* p is not present */
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
			DEBUG0 ("symamd: p not present\n") ;
			return (FALSE) ;
		}

		if (n < 0)		/* n must be >= 0 */
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
			stats [COLAMD_INFO1] = n ;
			DEBUG0 ("symamd: n negative %d\n", n) ;
			return (FALSE) ;
		}

		nnz = p [n] ;
		if (nnz < 0)	/* nnz must be >= 0 */
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
			stats [COLAMD_INFO1] = nnz ;
			DEBUG0 ("symamd: number of entries negative %d\n", nnz) ;
			return (FALSE) ;
		}

		if (p [0] != 0)
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero ;
			stats [COLAMD_INFO1] = p [0] ;
			DEBUG0 ("symamd: p[0] not zero %d\n", p [0]) ;
			return (FALSE) ;
		}

		/* === If no knobs, set default knobs =============================== */

		if (knobs == null)
		{
			COLAMD_set_defaults (default_knobs) ;
			knobs = default_knobs ;
		}

		/* === Allocate count and mark ====================================== */

		try
		{
			count = new int [n+1] ;
		} catch (OutOfMemoryError e) {
			stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
			DEBUG0 ("symamd: allocate count (size %d) failed\n", n+1) ;
			return (FALSE) ;
		}

		try
		{
			mark = new int[n+1] ;
		} catch (OutOfMemoryError e) {
			stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
			count = null ;
			DEBUG0 ("symamd: allocate mark (size %d) failed\n", n+1) ;
			return (FALSE) ;
		}

		/* === Compute column counts of M, check if A is valid ============== */

		stats [COLAMD_INFO3] = 0 ;  /* number of duplicate or unsorted row indices*/

		for (i = 0 ; i < n ; i++)
		{
			mark [i] = -1 ;
		}

		for (j = 0 ; j < n ; j++)
		{
			last_row = -1 ;

			length = p [j+1] - p [j] ;
			if (length < 0)
			{
				/* column pointers must be non-decreasing */
				stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
				stats [COLAMD_INFO1] = j ;
				stats [COLAMD_INFO2] = length ;
				count = null ;
				mark = null ;
				DEBUG0 ("symamd: col %d negative length %d\n", j, length) ;
				return (FALSE) ;
			}

			for (pp = p [j] ; pp < p [j+1] ; pp++)
			{
				i = A [pp] ;
				if (i < 0 || i >= n)
				{
					/* row index i, in column j, is out of bounds */
					stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
					stats [COLAMD_INFO1] = j ;
					stats [COLAMD_INFO2] = i ;
					stats [COLAMD_INFO3] = n ;
					count = null ;
					mark = null ;
					DEBUG0 ("symamd: row %d col %d out of bounds\n", i, j) ;
					return (FALSE) ;
				}

				if (i <= last_row || mark [i] == j)
				{
					/* row index is unsorted or repeated (or both), thus col */
					/* is jumbled.  This is a notice, not an error condition. */
					stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
					stats [COLAMD_INFO1] = j ;
					stats [COLAMD_INFO2] = i ;
					(stats [COLAMD_INFO3]) ++ ;
					DEBUG1 ("symamd: row %d col %d unsorted/duplicate\n", i, j) ;
				}

				if (i > j && mark [i] != j)
				{
					/* row k of M will contain column indices i and j */
					count [i]++ ;
					count [j]++ ;
				}

				/* mark the row as having been seen in this column */
				mark [i] = j ;

				last_row = i ;
			}
		}

		/* === Compute column pointers of M ================================= */

		/* use output permutation, perm, for column pointers of M */
		perm [0] = 0 ;
		for (j = 1 ; j <= n ; j++)
		{
			perm [j] = perm [j-1] + count [j-1] ;
		}
		for (j = 0 ; j < n ; j++)
		{
			count [j] = perm [j] ;
		}

		/* === Construct M ================================================== */

		mnz = perm [n] ;
		n_row = mnz / 2 ;
		Mlen = COLAMD_recommended (mnz, n_row, n) ;
		try
		{
			M = new int [Mlen] ;
			DEBUG0 ("symamd: M is %d-by-%d with %d entries, Mlen = %d\n",
					n_row, n, mnz, Mlen) ;
		} catch (OutOfMemoryError e) {
			stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
			count = null ;
			mark = null;
			DEBUG0 ("symamd: allocate M (size %g) failed\n", (double) Mlen) ;
			return (FALSE) ;
		}

		k = 0 ;

		if (stats [COLAMD_STATUS] == COLAMD_OK)
		{
			/* Matrix is OK */
			for (j = 0 ; j < n ; j++)
			{
				ASSERT (p [j+1] - p [j] >= 0) ;
				for (pp = p [j] ; pp < p [j+1] ; pp++)
				{
					i = A [pp] ;
					ASSERT (i >= 0 && i < n) ;
					if (i > j)
					{
						/* row k of M contains column indices i and j */
						M [count [i]++] = k ;
						M [count [j]++] = k ;
						k++ ;
					}
				}
			}
		}
		else
		{
			/* Matrix is jumbled.  Do not add duplicates to M.  Unsorted cols OK. */
			DEBUG0 ("symamd: Duplicates in A.\n") ;
			for (i = 0 ; i < n ; i++)
			{
				mark [i] = -1 ;
			}
			for (j = 0 ; j < n ; j++)
			{
				ASSERT (p [j+1] - p [j] >= 0) ;
				for (pp = p [j] ; pp < p [j+1] ; pp++)
				{
					i = A [pp] ;
					ASSERT (i >= 0 && i < n) ;
					if (i > j && mark [i] != j)
					{
						/* row k of M contains column indices i and j */
						M [count [i]++] = k ;
						M [count [j]++] = k ;
						k++ ;
						mark [i] = j ;
					}
				}
			}
		}

		/* count and mark no longer needed */
		count = null ;
		mark = null ;
		ASSERT (k == n_row) ;

		/* === Adjust the knobs for M ======================================= */

		for (i = 0 ; i < COLAMD_KNOBS ; i++)
		{
			cknobs [i] = knobs [i] ;
		}

		/* there are no dense rows in M */
		cknobs [COLAMD_DENSE_ROW] = -1 ;
		cknobs [COLAMD_DENSE_COL] = knobs [COLAMD_DENSE_ROW] ;

		/* === Order the columns of M ======================================= */

		colamd (n_row, n, Mlen, M, perm, cknobs, stats) ;

		/* Note that the output permutation is now in perm */

		/* === get the statistics for symamd from colamd ==================== */

		/* a dense column in colamd means a dense row and col in symamd */
		stats [COLAMD_DENSE_ROW] = stats [COLAMD_DENSE_COL] ;

		/* === Free M ======================================================= */

		M = null ;
		DEBUG0 ("symamd: done.\n") ;
		return (TRUE) ;
	}

/* ========================================================================== */
/* === colamd =============================================================== */
/* ========================================================================== */

	/*
	    The colamd routine computes a column ordering Q of a sparse matrix
	    A such that the LU factorization P(AQ) = LU remains sparse, where P is
	    selected via partial pivoting.   The routine can also be viewed as
	    providing a permutation Q such that the Cholesky factorization
	    (AQ)'(AQ) = LL' remains sparse.
	*/

	/**
	 * Computes a column ordering (Q) of A such that P(AQ)=LU or
	 * (AQ)'AQ=LL' have less fill-in and require fewer floating point
	 * operations than factorizing the unpermuted matrix A or A'A,
	 * respectively.
	 *
	 * @param n_row number of rows in A. Restriction:  n_row >= 0.
	 * Colamd returns FALSE if n_row is negative.
	 * @param n_col number of columns in A. Restriction:  n_col >= 0.
	 * Colamd returns FALSE if n_col is negative.
	 * @param Alen length of A. Restriction (see note):
	 * Alen >= 2*nnz + 6*(n_col+1) + 4*(n_row+1) + n_col
	 * Colamd returns FALSE if these conditions are not met.
	 *
	 * Note:  this restriction makes an modest assumption regarding
	 * the size of the two typedef's structures in colamd.h.
	 * We do, however, guarantee that
	 *
	 *     Alen >= colamd_recommended (nnz, n_row, n_col)
	 *
	 * will be sufficient.  Note: the macro version does not check
	 * for integer overflow, and thus is not recommended.  Use
	 * the colamd_recommended routine instead.
	 * @param A row indices of A.
	 *
	 * A is an integer array of size Alen.  Alen must be at least as
	 * large as the bare minimum value given above, but this is very
	 * low, and can result in excessive run time.  For best
	 * performance, we recommend that Alen be greater than or equal to
	 * colamd_recommended (nnz, n_row, n_col), which adds
	 * nnz/5 to the bare minimum value given above.
	 *
	 * On input, the row indices of the entries in column c of the
	 * matrix are held in A [(p [c]) ... (p [c+1]-1)].  The row indices
	 * in a given column c need not be in ascending order, and
	 * duplicate row indices may be be present.  However, colamd will
	 * work a little faster if both of these conditions are met
	 * (Colamd puts the matrix into this format, if it finds that the
	 * the conditions are not met).
	 *
	 * The matrix is 0-based.  That is, rows are in the range 0 to
	 * n_row-1, and columns are in the range 0 to n_col-1.  Colamd
	 * returns FALSE if any row index is out of range.
	 *
	 * The contents of A are modified during ordering, and are
	 * undefined on output.
	 * @param p pointers to columns in A.
	 *
	 * p is an integer array of size n_col+1.  On input, it holds the
	 * "pointers" for the column form of the matrix A.  Column c of
	 * the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
	 * entry, p [0], must be zero, and p [c] <= p [c+1] must hold
	 * for all c in the range 0 to n_col-1.  The value p [n_col] is
	 * thus the total number of entries in the pattern of the matrix A.
	 * Colamd returns FALSE if these conditions are not met.
	 *
	 * On output, if colamd returns TRUE, the array p holds the column
	 * permutation (Q, for P(AQ)=LU or (AQ)'(AQ)=LL'), where p [0] is
	 * the first column index in the new ordering, and p [n_col-1] is
	 * the last.  That is, p [k] = j means that column j of A is the
	 * kth pivot column, in AQ, where k is in the range 0 to n_col-1
	 * (p [0] = j means that column j of A is the first column in AQ).
	 *
	 * If colamd returns FALSE, then no permutation is returned, and
	 * p is undefined on output.
	 * @param knobs parameters (uses defaults if NULL)
	 * @param stats output statistics and error codes.
	 *
	 * Statistics on the ordering, and error status.
	 * Colamd returns FALSE if stats is not present.
	 *
	 * stats [0]:  number of dense or empty rows ignored.
	 *
	 * stats [1]:  number of dense or empty columns ignored (and
	 * 		ordered last in the output permutation p)
	 * 		Note that a row can become "empty" if it
	 * 		contains only "dense" and/or "empty" columns,
	 * 		and similarly a column can become "empty" if it
	 * 		only contains "dense" and/or "empty" rows.
	 *
	 * stats [2]:  number of garbage collections performed.
	 * 		This can be excessively high if Alen is close
	 * 		to the minimum required value.
	 *
	 * stats [3]:  status code.  < 0 is an error code.
	 * 	    > 1 is a warning or notice.
	 *
	 * 		0	OK.  Each column of the input matrix contained
	 * 			row indices in increasing order, with no
	 * 			duplicates.
	 *
	 * 		1	OK, but columns of input matrix were jumbled
	 * 			(unsorted columns or duplicate entries).  Colamd
	 * 			had to do some extra work to sort the matrix
	 * 			first and remove duplicate entries, but it
	 * 			still was able to return a valid permutation
	 * 			(return value of colamd was TRUE).
	 *
	 * 				stats [4]: highest numbered column that
	 * 					is unsorted or has duplicate
	 * 					entries.
	 * 				stats [5]: last seen duplicate or
	 * 					unsorted row index.
	 * 				stats [6]: number of duplicate or
	 * 					unsorted row indices.
	 *
	 * 		-1	A is a null pointer
	 *
	 * 		-2	p is a null pointer
	 *
	 * 		-3 	n_row is negative
	 *
	 * 				stats [4]: n_row
	 *
	 * 		-4	n_col is negative
	 *
	 * 				stats [4]: n_col
	 *
	 * 		-5	number of nonzeros in matrix is negative
	 *
	 * 				stats [4]: number of nonzeros, p [n_col]
	 *
	 * 		-6	p [0] is nonzero
	 *
	 * 				stats [4]: p [0]
	 *
	 * 		-7	A is too small
	 *
	 * 				stats [4]: required size
	 * 				stats [5]: actual size (Alen)
	 *
	 * 		-8	a column has a negative number of entries
	 *
	 * 				stats [4]: column with < 0 entries
	 * 				stats [5]: number of entries in col
	 *
	 * 		-9	a row index is out of bounds
	 *
	 * 				stats [4]: column with bad row index
	 * 				stats [5]: bad row index
	 * 				stats [6]: n_row, # of rows of matrx
	 *
	 * 		-10	(unused; see symamd.c)
	 *
	 * 		-999	(unused; see symamd.c)
	 *
	 * Future versions may return more statistics in the stats array.
	 * @return TRUE if successful, FALSE otherwise
	 */
	public static int colamd (int n_row, int n_col, int Alen, int[] A,
			int[] p, double[] knobs, int[] stats)
	{
		/* === Local variables ============================================== */

		int i ;			/* loop index */
		int nnz ;			/* nonzeros in A */
		int Row_size ;		/* size of Row [], in integers */
		int Col_size ;		/* size of Col [], in integers */
		int need ;		/* minimum required length of A */
		Colamd_Row[] Row ;		/* pointer into A of Row [0..n_row] array */
		Colamd_Col[] Col ;		/* pointer into A of Col [0..n_col] array */
		int[] n_col2 = new int [1] ;		/* number of non-dense, non-empty columns */
		int[] n_row2 = new int [1] ;		/* number of non-dense, non-empty rows */
		int ngarbage ;		/* number of garbage collections performed */
		int[] max_deg = new int [1] ;		/* maximum row degree */
		double[] default_knobs = new double[COLAMD_KNOBS] ;	/* default knobs array */
		int aggressive ;		/* do aggressive absorption */
		int[] ok ;

		if (!NDEBUG)
		{
//			colamd_get_debug ("colamd") ;
		}

		/* === Check the input arguments ==================================== */

		if (stats == null)
		{
			DEBUG0 ("colamd: stats not present\n") ;
			return (FALSE) ;
		}
		for (i = 0 ; i < COLAMD_STATS ; i++)
		{
			stats [i] = 0 ;
		}
		stats [COLAMD_STATUS] = COLAMD_OK ;
		stats [COLAMD_INFO1] = -1 ;
		stats [COLAMD_INFO2] = -1 ;

		if (A == null)		/* A is not present */
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
			DEBUG0 ("colamd: A not present\n") ;
			return (FALSE) ;
		}

		if (p == null)		/* p is not present */
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
			DEBUG0 ("colamd: p not present\n") ;
			return (FALSE) ;
		}

		if (n_row < 0)	/* n_row must be >= 0 */
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_nrow_negative ;
			stats [COLAMD_INFO1] = n_row ;
			DEBUG0 ("colamd: nrow negative %d\n", n_row) ;
			return (FALSE) ;
		}

		if (n_col < 0)	/* n_col must be >= 0 */
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
			stats [COLAMD_INFO1] = n_col ;
			DEBUG0 ("colamd: ncol negative %d\n", n_col) ;
			return (FALSE) ;
		}

		nnz = p [n_col] ;
		if (nnz < 0)	/* nnz must be >= 0 */
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
			stats [COLAMD_INFO1] = nnz ;
			DEBUG0 ("colamd: number of entries negative %d\n", nnz) ;
			return (FALSE) ;
		}

		if (p [0] != 0)
		{
			stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero	;
			stats [COLAMD_INFO1] = p [0] ;
			DEBUG0 ("colamd: p[0] not zero %d\n", p [0]) ;
			return (FALSE) ;
		}

		/* === If no knobs, set default knobs =============================== */

		if (knobs == null)
		{
			COLAMD_set_defaults (default_knobs) ;
			knobs = default_knobs ;
		}

		aggressive = (knobs [COLAMD_AGGRESSIVE] != FALSE) ? 1 : 0;

		/* === Allocate the Row and Col arrays from array A ================= */

		ok = new int [] { TRUE } ;
		Col_size = COLAMD_C (n_col, ok) ;	    /* size of Col array of structs */
		Row_size = COLAMD_R (n_row, ok) ;	    /* size of Row array of structs */

		/* need = 2*nnz + n_col + Col_size + Row_size ; */
		need = t_mult (nnz, 2, ok) ;
		need = t_add (need, n_col, ok) ;
//		need = t_add (need, Col_size, ok) ;
//		need = t_add (need, Row_size, ok) ;

		if ((ok[0] == 0) || need > (int) Alen || need > Int_MAX)
		{
			/* not enough space in array A to perform the ordering */
			stats [COLAMD_STATUS] = COLAMD_ERROR_A_too_small ;
			stats [COLAMD_INFO1] = need ;
			stats [COLAMD_INFO2] = Alen ;
			DEBUG0 ("colamd: Need Alen >= %d, given only Alen = %d\n", need,Alen);
			return (FALSE) ;
		}

//		Alen -= Col_size + Row_size ;
		Col = new Colamd_Col [Col_size] ;  //A [Alen] ;
		Row = new Colamd_Row [Row_size] ;  //A [Alen + Col_size] ;

		for (i = 0; i < Col_size; i++)
			Col [i] = new Colamd_Col() ;
		for (i = 0; i < Row_size; i++)
			Row [i] = new Colamd_Row() ;

		/* === Construct the row and column data structures ================= */

		if (init_rows_cols (n_row, n_col, Row, Col, A, p, stats) == 0)
		{
			/* input matrix is invalid */
			DEBUG0 ("colamd: Matrix invalid\n") ;
			return (FALSE) ;
		}

		/* === Initialize scores, kill dense rows/columns =================== */

		init_scoring (n_row, n_col, Row, Col, A, p, knobs,
				n_row2, n_col2, max_deg) ;

		/* === Order the supercolumns ======================================= */

		ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
				n_col2[0], max_deg[0], 2*nnz, aggressive) ;

		/* === Order the non-principal columns ============================== */

		order_children (n_col, Col, p) ;

		/* === Return statistics in stats =================================== */

		stats [COLAMD_DENSE_ROW] = n_row - n_row2[0] ;
		stats [COLAMD_DENSE_COL] = n_col - n_col2[0] ;
		stats [COLAMD_DEFRAG_COUNT] = ngarbage ;
		DEBUG0 ("colamd: done.\n") ;
		return (TRUE) ;
	}


/* ========================================================================== */
/* === colamd_report ======================================================== */
/* ========================================================================== */

	public static void COLAMD_report (int[] stats)
	{
		print_report ("colamd", stats) ;
	}


/* ========================================================================== */
/* === symamd_report ======================================================== */
/* ========================================================================== */

	public static void SYMAMD_report (int[] stats)
	{
		print_report ("symamd", stats) ;
	}



/* ========================================================================== */
/* === NON-USER-CALLABLE ROUTINES: ========================================== */
/* ========================================================================== */

/* There are no user-callable routines beyond this point in the file */


/* ========================================================================== */
/* === init_rows_cols ======================================================= */
/* ========================================================================== */

	/**
	 * Takes the column form of the matrix in A and creates the row form of the
	 * matrix.  Also, row and column attributes are stored in the Col and Row
	 * structs.  If the columns are un-sorted or contain duplicate row indices,
	 * this routine will also sort and remove duplicate row indices from the
	 * column form of the matrix.  Returns FALSE if the matrix is invalid,
	 * TRUE otherwise.  Not user-callable.
	 *
	 * @param n_row number of rows of A
	 * @param n_col number of columns of A
	 * @param Row of size n_row+1
	 * @param Col of size n_col+1
	 * @param A row indices of A, of size Alen
	 * @param p pointers to columns in A, of size n_col+1
	 * @param stats colamd statistics
	 * @return TRUE if OK, or FALSE otherwise
	 */
	private static int init_rows_cols (int n_row, int n_col, Colamd_Row[] Row,
			Colamd_Col[] Col, int[] A, int[] p, int[] stats)
	{
		/* === Local variables ============================================== */

		int col ;		/* a column index */
		int row ;		/* a row index */
		int cp ;		/* a column pointer */
		int cp_end ;		/* a pointer to the end of a column */
		int rp ;		/* a row pointer */
		int rp_end ;		/* a pointer to the end of a row */
		int last_row ;		/* previous row */

		/* === Initialize columns, and check column pointers ================ */

		for (col = 0 ; col < n_col ; col++)
		{
			Col [col].start = p [col] ;
			Col [col].length = p [col+1] - p [col] ;

			if (Col [col].length < 0)
			{
				/* column pointers must be non-decreasing */
				stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
				stats [COLAMD_INFO1] = col ;
				stats [COLAMD_INFO2] = Col [col].length ;
				DEBUG0 ("colamd: col %d length %d < 0\n", col, Col [col].length) ;
				return (FALSE) ;
			}

			Col [col].thickness(1) ;
			Col [col].score(0) ;
			Col [col].prev(EMPTY) ;
			Col [col].degree_next(EMPTY) ;
		}

		/* p [0..n_col] no longer needed, used as "head" in subsequent routines */

		/* === Scan columns, compute row degrees, and check row indices ===== */

		stats [COLAMD_INFO3] = 0 ;	/* number of duplicate or unsorted row indices*/

		for (row = 0 ; row < n_row ; row++)
		{
			Row [row].length = 0 ;
			Row [row].mark(-1) ;
		}

		for (col = 0 ; col < n_col ; col++)
		{
			last_row = -1 ;

			cp = p [col] ;
			cp_end = p [col+1] ;

			while (cp < cp_end)
			{
				row = A [cp++] ;

				/* make sure row indices within range */
				if (row < 0 || row >= n_row)
				{
					stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
					stats [COLAMD_INFO1] = col ;
					stats [COLAMD_INFO2] = row ;
					stats [COLAMD_INFO3] = n_row ;
					DEBUG0 ("colamd: row %d col %d out of bounds\n", row, col) ;
					return (FALSE) ;
				}

				if (row <= last_row || Row [row].mark() == col)
				{
					/* row index are unsorted or repeated (or both), thus col */
					/* is jumbled.  This is a notice, not an error condition. */
					stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
					stats [COLAMD_INFO1] = col ;
					stats [COLAMD_INFO2] = row ;
					(stats [COLAMD_INFO3]) ++ ;
					DEBUG1 ("colamd: row %d col %d unsorted/duplicate\n",row,col);
				}

				if (Row [row].mark() != col)
				{
					Row [row].length++ ;
				}
				else
				{
					/* this is a repeated entry in the column, */
					/* it will be removed */
					Col [col].length-- ;
				}

				/* mark the row as having been seen in this column */
				Row [row].mark(col) ;

				last_row = row ;
			}
		}

		/* === Compute row pointers ========================================= */

		/* row form of the matrix starts directly after the column */
		/* form of matrix in A */
		Row [0].start = p [n_col] ;
		Row [0].p(Row [0].start) ;
		Row [0].mark(-1) ;
		for (row = 1 ; row < n_row ; row++)
		{
			Row [row].start = Row [row-1].start + Row [row-1].length ;
			Row [row].p(Row [row].start) ;
			Row [row].mark(-1) ;
		}

		/* === Create row form ============================================== */

		if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
		{
			/* if cols jumbled, watch for repeated row indices */
			for (col = 0 ; col < n_col ; col++)
			{
				cp = p [col] ;
				cp_end = p [col+1] ;

				while (cp < cp_end)
				{
					row = A [cp++] ;

					if (Row [row].mark() != col)
					{
						A [Row [row].p()] = col ;
						Row [row].p( Row [row].p() + 1 ) ;
						Row [row].mark(col) ;
					}
				}
			}
		}
		else
		{
			/* if cols not jumbled, we don't need the mark (this is faster) */
			for (col = 0 ; col < n_col ; col++)
			{
				cp = p [col] ;
				cp_end = p [col+1] ;
				while (cp < cp_end)
				{
					A [Row [A [cp]].p()] = col ;
					Row [A [cp]].p( Row [A [cp]].p() + 1 ) ;
					cp++ ;
				}
			}
		}

		/* === Clear the row marks and set row degrees ====================== */

		for (row = 0 ; row < n_row ; row++)
		{
			Row [row].mark(0) ;
			Row [row].degree(Row [row].length) ;
		}

		/* === See if we need to re-create columns ========================== */

		if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
		{
			DEBUG0 ("colamd: reconstructing column form, matrix jumbled\n") ;

			if (!NDEBUG)
			{
				/* make sure column lengths are correct */
				for (col = 0 ; col < n_col ; col++)
				{
					p [col] = Col [col].length ;
				}
				for (row = 0 ; row < n_row ; row++)
				{
					rp = Row [row].start ;
					rp_end = rp + Row [row].length ;
					while (rp < rp_end)
					{
						p [A [rp++]]-- ;
					}
				}
				for (col = 0 ; col < n_col ; col++)
				{
					ASSERT (p [col] == 0) ;
				}
				/* now p is all zero (different than when debugging is turned off) */
			} /* NDEBUG */

			/* === Compute col pointers ========================================= */

			/* col form of the matrix starts at A [0]. */
			/* Note, we may have a gap between the col form and the row */
			/* form if there were duplicate entries, if so, it will be */
			/* removed upon the first garbage collection */
			Col [0].start = 0 ;
			p [0] = Col [0].start ;
			for (col = 1 ; col < n_col ; col++)
			{
				/* note that the lengths here are for pruned columns, i.e. */
				/* no duplicate row indices will exist for these columns */
				Col [col].start = Col [col-1].start + Col [col-1].length ;
				p [col] = Col [col].start ;
			}

			/* === Re-create col form =========================================== */

			for (row = 0 ; row < n_row ; row++)
			{
				rp = Row [row].start ;
				rp_end = rp + Row [row].length ;
				while (rp < rp_end)
				{
					A [(p [A [rp++]])++] = row ;
				}
			}
		}

		/* === Done.  Matrix is not (or no longer) jumbled ================== */

		return (TRUE) ;
	}


/* ========================================================================== */
/* === init_scoring ========================================================= */
/* ========================================================================== */

	/**
	 * Kills dense or empty columns and rows, calculates an initial score for
	 * each column, and places all columns in the degree lists.init_rows_cols
	 *
	 * @param n_row number of rows of A
	 * @param n_col number of columns of A
	 * @param Row of size n_row+1
	 * @param Col of size n_col+1
	 * @param A column form and row form of A
	 * @param head of size n_col+1
	 * @param knobs parameters
	 * @param p_n_row2 size 1, number of non-dense, non-empty rows
	 * @param p_n_col2 size 1, number of non-dense, non-empty columns
	 * @param p_max_deg size 1, maximum row degree
	 */
	private static void init_scoring (int n_row, int n_col, Colamd_Row[] Row,
			Colamd_Col[] Col, int[] A, int[] head, double[] knobs,
			int[] p_n_row2, int[] p_n_col2, int[] p_max_deg)
	{
		/* === Local variables ============================================== */

		int c ;                /* a column index */
		int r, row ;           /* a row index */
		int cp ;               /* a column pointer */
		int deg ;              /* degree of a row or column */
		int cp_end ;           /* a pointer to the end of a column */
		int new_cp ;           /* new column pointer */
		int col_length ;       /* length of pruned column */
		int score ;            /* current column score */
		int n_col2 ;           /* number of non-dense, non-empty columns */
		int n_row2 ;           /* number of non-dense, non-empty rows */
		int dense_row_count ;  /* remove rows with more entries than this */
		int dense_col_count ;  /* remove cols with more entries than this */
		int min_score ;        /* smallest column score */
		int max_deg ;          /* maximum row degree */
		int next_col ;         /* Used to add to degree list.*/

		int debug_count = 0 ;  /* debug only. */

		/* === Extract knobs ================================================ */

		/* Note: if knobs contains a NaN, this is undefined: */
		if (knobs [COLAMD_DENSE_ROW] < 0)
		{
			/* only remove completely dense rows */
			dense_row_count = n_col-1 ;
		}
		else
		{
			dense_row_count = DENSE_DEGREE (knobs [COLAMD_DENSE_ROW], n_col) ;
		}
		if (knobs [COLAMD_DENSE_COL] < 0)
		{
			/* only remove completely dense columns */
			dense_col_count = n_row-1 ;
		}
		else
		{
			dense_col_count = DENSE_DEGREE (knobs [COLAMD_DENSE_COL],
					MIN (n_row, n_col)) ;
		}

		DEBUG1 ("colamd: densecount: %d %d\n", dense_row_count, dense_col_count) ;
		max_deg = 0 ;
		n_col2 = n_col ;
		n_row2 = n_row ;

		/* === Kill empty columns =========================================== */

		/* Put the empty columns at the end in their natural order, so that LU */
		/* factorization can proceed as far as possible. */
		for (c = n_col-1 ; c >= 0 ; c--)
		{
			deg = Col [c].length ;
			if (deg == 0)
			{
				/* this is a empty column, kill and order it last */
				Col [c].order(--n_col2) ;
				KILL_PRINCIPAL_COL (Col, c) ;
			}
		}
		DEBUG1 ("colamd: null columns killed: %d\n", n_col - n_col2) ;

		/* === Kill dense columns =========================================== */

		/* Put the dense columns at the end, in their natural order */
		for (c = n_col-1 ; c >= 0 ; c--)
		{
			/* skip any dead columns */
			if (COL_IS_DEAD (Col, c))
			{
				continue ;
			}
			deg = Col [c].length ;
			if (deg > dense_col_count)
			{
				/* this is a dense column, kill and order it last */
				Col [c].order(--n_col2) ;
				/* decrement the row degrees */
				cp = Col [c].start ;
				cp_end = cp + Col [c].length ;
				while (cp < cp_end)
				{
					Row [A [cp]].degree( Row [A [cp]].degree() - 1 ) ;
					cp++ ;
				}
				KILL_PRINCIPAL_COL (Col, c) ;
			}
		}
		DEBUG1 ("colamd: Dense and null columns killed: %d\n", n_col - n_col2) ;

		/* === Kill dense and empty rows ==================================== */

		for (r = 0 ; r < n_row ; r++)
		{
			deg = Row [r].degree() ;
			ASSERT (deg >= 0 && deg <= n_col) ;
			if (deg > dense_row_count || deg == 0)
			{
				/* kill a dense or empty row */
				KILL_ROW (Row, r) ;
				--n_row2 ;
			}
			else
			{
				/* keep track of max degree of remaining rows */
				max_deg = MAX (max_deg, deg) ;
			}
		}
		DEBUG1 ("colamd: Dense and null rows killed: %d\n", n_row - n_row2) ;

		/* === Compute initial column scores ================================ */

		/* At this point the row degrees are accurate.  They reflect the number */
		/* of "live" (non-dense) columns in each row.  No empty rows exist. */
		/* Some "live" columns may contain only dead rows, however.  These are */
		/* pruned in the code below. */

		/* now find the initial matlab score for each column */
		for (c = n_col-1 ; c >= 0 ; c--)
		{
			/* skip dead column */
			if (COL_IS_DEAD (Col, c))
			{
				continue ;
			}
			score = 0 ;
			cp = Col [c].start ;
			new_cp = cp ;
			cp_end = cp + Col [c].length ;
			while (cp < cp_end)
			{
				/* get a row */
				row = A [cp++] ;
				/* skip if dead */
				if (ROW_IS_DEAD (Row, row))
				{
					continue ;
				}
				/* compact the column */
				A [new_cp++] = row ;
				/* add row's external degree */
				score += Row [row].degree() - 1 ;
				/* guard against integer overflow */
				score = MIN (score, n_col) ;
			}
			/* determine pruned column length */
			col_length = (new_cp - Col [c].start) ;
			if (col_length == 0)
			{
				/* a newly-made null column (all rows in this col are "dense" */
				/* and have already been killed) */
				DEBUG2 ("Newly null killed: %d\n", c) ;
				Col [c].order(--n_col2) ;
				KILL_PRINCIPAL_COL (Col, c) ;
			}
			else
			{
				/* set column length and set score */
				ASSERT (score >= 0) ;
				ASSERT (score <= n_col) ;
				Col [c].length = col_length ;
				Col [c].score(score) ;
			}
		}
		DEBUG1 ("colamd: Dense, null, and newly-null columns killed: %d\n",
				n_col-n_col2) ;

		/* At this point, all empty rows and columns are dead.  All live columns */
		/* are "clean" (containing no dead rows) and simplicial (no supercolumns */
		/* yet).  Rows may contain dead columns, but all live rows contain at */
		/* least one live column. */

		if (!NDEBUG)
		{
			debug_structures (n_row, n_col, Row, Col, A, n_col2) ;
		}

		/* === Initialize degree lists ========================================== */

		if (!NDEBUG)
		{
			debug_count = 0 ;
		}

		/* clear the hash buckets */
		for (c = 0 ; c <= n_col ; c++)
		{
			head [c] = EMPTY ;
		}
		min_score = n_col ;
		/* place in reverse order, so low column indices are at the front */
		/* of the lists.  This is to encourage natural tie-breaking */
		for (c = n_col-1 ; c >= 0 ; c--)
		{
			/* only add principal columns to degree lists */
			if (COL_IS_ALIVE (Col, c))
			{
				DEBUG4 ("place %d score %d minscore %d ncol %d\n",
						c, Col [c].score(), min_score, n_col) ;

				/* === Add columns score to DList =============================== */

				score = Col [c].score() ;

				ASSERT (min_score >= 0) ;
				ASSERT (min_score <= n_col) ;
				ASSERT (score >= 0) ;
				ASSERT (score <= n_col) ;
				ASSERT (head [score] >= EMPTY) ;

				/* now add this column to dList at proper score location */
				next_col = head [score] ;
				Col [c].prev(EMPTY) ;
				Col [c].degree_next(next_col) ;

				/* if there already was a column with the same score, set its */
				/* previous pointer to this new column */
				if (next_col != EMPTY)
				{
					Col [next_col].prev(c) ;
				}
				head [score] = c ;

				/* see if this score is less than current min */
				min_score = MIN (min_score, score) ;

				if (!NDEBUG)
				{
					debug_count++ ;
				}

			}
		}

		if (!NDEBUG)
		{
			DEBUG1 ("colamd: Live cols %d out of %d, non-princ: %d\n",
					debug_count, n_col, n_col-debug_count) ;
			ASSERT (debug_count == n_col2) ;
			debug_deg_lists (n_row, n_col, Row, Col, head, min_score, n_col2, max_deg) ;
		} /* NDEBUG */

		/* === Return number of remaining columns, and max row degree ======= */

		p_n_col2[0] = n_col2 ;
		p_n_row2[0] = n_row2 ;
		p_max_deg[0] = max_deg ;
	}


/* ========================================================================== */
/* === find_ordering ======================================================== */
/* ========================================================================== */

	/**
	 * Order the principal columns of the supercolumn form of the matrix
	 * (no supercolumns on input).  Uses a minimum approximate column minimum
	 * degree ordering method.  Not user-callable.
	 *
	 * @param n_row number of rows of A
	 * @param n_col number of columns of A
	 * @param Alen size of A, 2*nnz + n_col or larger
	 * @param Row of size n_row+1
	 * @param Col of size n_col+1
	 * @param A column form and row form of A
	 * @param head of size n_col+1
	 * @param n_col2 Remaining columns to order
	 * @param max_deg Maximum row degree
	 * @param pfree index of first free slot (2*nnz on entry)
	 * @param aggressive
	 * @return the number of garbage collections
	 */
	private static int find_ordering (int n_row, int n_col, int Alen,
			Colamd_Row[] Row, Colamd_Col[] Col, int[] A, int[] head, int n_col2,
			int max_deg, int pfree, int aggressive)
	{
		/* === Local variables ============================================== */

		int k ;                 /* current pivot ordering step */
		int pivot_col ;         /* current pivot column */
		int cp ;                /* a column pointer */
		int rp ;                /* a row pointer */
		int pivot_row ;         /* current pivot row */
		int new_cp ;            /* modified column pointer */
		int new_rp ;            /* modified row pointer */
		int pivot_row_start ;   /* pointer to start of pivot row */
		int pivot_row_degree ;  /* number of columns in pivot row */
		int pivot_row_length ;  /* number of supercolumns in pivot row */
		int pivot_col_score ;   /* score of pivot column */
		int needed_memory ;     /* free space needed for pivot row */
		int cp_end ;            /* pointer to the end of a column */
		int rp_end ;            /* pointer to the end of a row */
		int row ;               /* a row index */
		int col ;               /* a column index */
		int max_score ;         /* maximum possible score */
		int cur_score ;         /* score of current column */
		/* FIXME unsigned */ int hash ; /* hash value for supernode detection */
		int head_column ;       /* head of hash bucket */
		int first_col ;         /* first column in hash bucket */
		int tag_mark ;          /* marker value for mark array */
		int row_mark ;          /* Row [row].shared2.mark */
		int set_difference ;    /* set difference size of row with pivot row */
		int min_score ;         /* smallest column score */
		int col_thickness ;     /* "thickness" (no. of columns in a supercol) */
		int max_mark ;          /* maximum value of tag_mark */
		int pivot_col_thickness ; /* number of columns represented by pivot col */
		int prev_col ;          /* Used by Dlist operations. */
		int next_col ;          /* Used by Dlist operations. */
		int ngarbage ;          /* number of garbage collections performed */

		int debug_d ;           /* debug loop counter */
		int debug_step = 0 ;    /* debug loop counter */

		/* === Initialization and clear mark ================================ */

		max_mark = Int_MAX - n_col ;
		tag_mark = clear_mark (0, max_mark, n_row, Row) ;
		min_score = 0 ;
		ngarbage = 0 ;
		DEBUG1 ("colamd: Ordering, n_col2=%d\n", n_col2) ;

		/* === Order the columns ============================================ */

		for (k = 0 ; k < n_col2 ; /* 'k' is incremented below */)
		{

			if (!NDEBUG)
			{
				if (debug_step % 100 == 0)
				{
					DEBUG2 ("\n...       Step k: %d out of n_col2: %d\n", k, n_col2) ;
				}
				else
				{
					DEBUG3 ("\n----------Step k: %d out of n_col2: %d\n", k, n_col2) ;
				}
				debug_step++ ;
				debug_deg_lists (n_row, n_col, Row, Col, head,
						min_score, n_col2-k, max_deg) ;
				debug_matrix (n_row, n_col, Row, Col, A) ;
			} /* NDEBUG */

			/* === Select pivot column, and order it ============================ */

			/* make sure degree list isn't empty */
			ASSERT (min_score >= 0) ;
			ASSERT (min_score <= n_col) ;
			ASSERT (head [min_score] >= EMPTY) ;

			if (!NDEBUG)
			{
				for (debug_d = 0 ; debug_d < min_score ; debug_d++)
				{
					ASSERT (head [debug_d] == EMPTY) ;
				}
			} /* NDEBUG */

			/* get pivot column from head of minimum degree list */
			while (head [min_score] == EMPTY && min_score < n_col)
			{
				min_score++ ;
			}
			pivot_col = head [min_score] ;
			ASSERT (pivot_col >= 0 && pivot_col <= n_col) ;
			next_col = Col [pivot_col].degree_next() ;
			head [min_score] = next_col ;
			if (next_col != EMPTY)
			{
				Col [next_col].prev(EMPTY) ;
			}

			ASSERT (COL_IS_ALIVE (Col, pivot_col)) ;

			/* remember score for defrag check */
			pivot_col_score = Col [pivot_col].score() ;

			/* the pivot column is the kth column in the pivot order */
			Col [pivot_col].order(k) ;

			/* increment order count by column thickness */
			pivot_col_thickness = Col [pivot_col].thickness() ;
			k += pivot_col_thickness ;
			ASSERT (pivot_col_thickness > 0) ;
			DEBUG3 ("Pivot col: %d thick %d\n", pivot_col, pivot_col_thickness) ;

			/* === Garbage_collection, if necessary ============================= */

			needed_memory = MIN (pivot_col_score, n_col - k) ;
			if (pfree + needed_memory >= Alen)
			{
				pfree = garbage_collection (n_row, n_col, Row, Col, A, pfree) ;
				ngarbage++ ;
				/* after garbage collection we will have enough */
				ASSERT (pfree + needed_memory < Alen) ;
				/* garbage collection has wiped out the Row[].shared2.mark array */
				tag_mark = clear_mark (0, max_mark, n_row, Row) ;

				if (!NDEBUG)
				{
					debug_matrix (n_row, n_col, Row, Col, A) ;
				} /* NDEBUG */
			}

			/* === Compute pivot row pattern ==================================== */

			/* get starting location for this new merged row */
			pivot_row_start = pfree ;

			/* initialize new row counts to zero */
			pivot_row_degree = 0 ;

			/* tag pivot column as having been visited so it isn't included */
			/* in merged pivot row */
			Col [pivot_col].thickness(-pivot_col_thickness) ;

			/* pivot row is the union of all rows in the pivot column pattern */
			cp = Col [pivot_col].start ;
			cp_end = cp + Col [pivot_col].length ;
			while (cp < cp_end)
			{
				/* get a row */
				row = A [cp++] ;
				DEBUG4 ("Pivot col pattern %d %d\n", ROW_IS_ALIVE (Row, row) ? 1 : 0, row) ;
				/* skip if row is dead */
				if (ROW_IS_ALIVE (Row, row))
				{
					rp = Row [row].start ;
					rp_end = rp + Row [row].length ;
					while (rp < rp_end)
					{
						/* get a column */
						col = A [rp++] ;
						/* add the column, if alive and untagged */
						col_thickness = Col [col].thickness() ;
						if (col_thickness > 0 && COL_IS_ALIVE (Col, col))
						{
							/* tag column in pivot row */
							Col [col].thickness(-col_thickness) ;
							ASSERT (pfree < Alen) ;
							/* place column in pivot row */
							A [pfree++] = col ;
							pivot_row_degree += col_thickness ;
						}
					}
				}
			}

			/* clear tag on pivot column */
			Col [pivot_col].thickness(pivot_col_thickness) ;
			max_deg = MAX (max_deg, pivot_row_degree) ;

			if (!NDEBUG)
			{
				DEBUG3 ("check2\n") ;
				debug_mark (n_row, Row, tag_mark, max_mark) ;
			} /* NDEBUG */

			/* === Kill all rows used to construct pivot row ==================== */

			/* also kill pivot row, temporarily */
			cp = Col [pivot_col].start ;
			cp_end = cp + Col [pivot_col].length ;
			while (cp < cp_end)
			{
				/* may be killing an already dead row */
				row = A [cp++] ;
				DEBUG3 ("Kill row in pivot col: %d\n", row) ;
				KILL_ROW (Row, row) ;
			}

			/* === Select a row index to use as the new pivot row =============== */

			pivot_row_length = pfree - pivot_row_start ;
			if (pivot_row_length > 0)
			{
				/* pick the "pivot" row arbitrarily (first row in col) */
				pivot_row = A [Col [pivot_col].start] ;
				DEBUG3 ("Pivotal row is %d\n", pivot_row) ;
			}
			else
			{
				/* there is no pivot row, since it is of zero length */
				pivot_row = EMPTY ;
				ASSERT (pivot_row_length == 0) ;
			}
			ASSERT (Col [pivot_col].length > 0 || pivot_row_length == 0) ;

			/* === Approximate degree computation =============================== */

			/* Here begins the computation of the approximate degree.  The column */
			/* score is the sum of the pivot row "length", plus the size of the */
			/* set differences of each row in the column minus the pattern of the */
			/* pivot row itself.  The column ("thickness") itself is also */
			/* excluded from the column score (we thus use an approximate */
			/* external degree). */

			/* The time taken by the following code (compute set differences, and */
			/* add them up) is proportional to the size of the data structure */
			/* being scanned - that is, the sum of the sizes of each column in */
			/* the pivot row.  Thus, the amortized time to compute a column score */
			/* is proportional to the size of that column (where size, in this */
			/* context, is the column "length", or the number of row indices */
			/* in that column).  The number of row indices in a column is */
			/* monotonically non-decreasing, from the length of the original */
			/* column on input to colamd. */

			/* === Compute set differences ====================================== */

			DEBUG3 ("** Computing set differences phase. **\n") ;

			/* pivot row is currently dead - it will be revived later. */

			DEBUG3 ("Pivot row: ") ;
			/* for each column in pivot row */
			rp = pivot_row_start ;
			rp_end = rp + pivot_row_length ;
			while (rp < rp_end)
			{
				col = A [rp++] ;
				ASSERT (COL_IS_ALIVE (Col, col) && col != pivot_col) ;
				DEBUG3 ("Col: %d\n", col) ;

				/* clear tags used to construct pivot row pattern */
				col_thickness = -Col [col].thickness() ;
				ASSERT (col_thickness > 0) ;
				Col [col].thickness(col_thickness) ;

				/* === Remove column from degree list =========================== */

				cur_score = Col [col].score() ;
				prev_col = Col [col].prev() ;
				next_col = Col [col].degree_next() ;
				ASSERT (cur_score >= 0) ;
				ASSERT (cur_score <= n_col) ;
				ASSERT (cur_score >= EMPTY) ;
				if (prev_col == EMPTY)
				{
					head [cur_score] = next_col ;
				}
				else
				{
					Col [prev_col].degree_next(next_col) ;
				}
				if (next_col != EMPTY)
				{
					Col [next_col].prev(prev_col) ;
				}

				/* === Scan the column ========================================== */

				cp = Col [col].start ;
				cp_end = cp + Col [col].length ;
				while (cp < cp_end)
				{
					/* get a row */
					row = A [cp++] ;
					row_mark = Row [row].mark() ;
					/* skip if dead */
					if (ROW_IS_MARKED_DEAD (row_mark))
					{
						continue ;
					}
					ASSERT (row != pivot_row) ;
					set_difference = row_mark - tag_mark ;
					/* check if the row has been seen yet */
					if (set_difference < 0)
					{
						ASSERT (Row [row].degree() <= max_deg) ;
						set_difference = Row [row].degree() ;
					}
					/* subtract column thickness from this row's set difference */
					set_difference -= col_thickness ;
					ASSERT (set_difference >= 0) ;
					/* absorb this row if the set difference becomes zero */
					if (set_difference == 0 && aggressive != 0)
					{
						DEBUG3 ("aggressive absorption. Row: %d\n", row) ;
						KILL_ROW (Row, row) ;
					}
					else
					{
						/* save the new mark */
						Row [row].mark(set_difference + tag_mark) ;
					}
				}
			}

			if (!NDEBUG)
			{
				debug_deg_lists (n_row, n_col, Row, Col, head,
						min_score, n_col2-k-pivot_row_degree, max_deg) ;
			} /* NDEBUG */

			/* === Add up set differences for each column ======================= */

			DEBUG3 ("** Adding set differences phase. **\n") ;

			/* for each column in pivot row */
			rp = pivot_row_start ;
			rp_end = rp + pivot_row_length ;
			while (rp < rp_end)
			{
				/* get a column */
				col = A [rp++] ;
				ASSERT (COL_IS_ALIVE (Col, col) && col != pivot_col) ;
				hash = 0 ;
				cur_score = 0 ;
				cp = Col [col].start ;
				/* compact the column */
				new_cp = cp ;
				cp_end = cp + Col [col].length ;

				DEBUG4 ("Adding set diffs for Col: %d.\n", col) ;

				while (cp < cp_end)
				{
					/* get a row */
					row = A [cp++] ;
					ASSERT(row >= 0 && row < n_row) ;
					row_mark = Row [row].mark() ;
					/* skip if dead */
					if (ROW_IS_MARKED_DEAD (row_mark))
					{
						DEBUG4 (" Row %d, dead\n", row) ;
						continue ;
					}
					DEBUG4 (" Row %d, set diff %d\n", row, row_mark-tag_mark);
					ASSERT (row_mark >= tag_mark) ;
					/* compact the column */
					A [new_cp++] = row ;
					/* compute hash function */
					hash += row ;
					/* add set difference */
					cur_score += row_mark - tag_mark ;
					/* integer overflow... */
					cur_score = MIN (cur_score, n_col) ;
				}

				/* recompute the column's length */
				Col [col].length = (int) (new_cp - Col [col].start) ;

				/* === Further mass elimination ================================= */

				if (Col [col].length == 0)
				{
					DEBUG4 ("further mass elimination. Col: %d\n", col) ;
					/* nothing left but the pivot row in this column */
					KILL_PRINCIPAL_COL (Col, col) ;
					pivot_row_degree -= Col [col].thickness() ;
					ASSERT (pivot_row_degree >= 0) ;
					/* order it */
					Col [col].order(k) ;
					/* increment order count by column thickness */
					k += Col [col].thickness() ;
				}
				else
				{
					/* === Prepare for supercolumn detection ==================== */

					DEBUG4 ("Preparing supercol detection for Col: %d.\n", col) ;

					/* save score so far */
					Col [col].score(cur_score) ;

					/* add column to hash table, for supercolumn detection */
					hash %= n_col + 1 ;

					DEBUG4 (" Hash = %d, n_col = %d.\n", hash, n_col) ;
					ASSERT (((int) hash) <= n_col) ;

					head_column = head [hash] ;
					if (head_column > EMPTY)
					{
						/* degree list "hash" is non-empty, use prev (shared3) of */
						/* first column in degree list as head of hash bucket */
						first_col = Col [head_column].headhash() ;
						Col [head_column].headhash(col) ;
					}
					else
					{
						/* degree list "hash" is empty, use head as hash bucket */
						first_col = - (head_column + 2) ;
						head [hash] = - (col + 2) ;
					}
					Col [col].hash_next(first_col) ;

					/* save hash function in Col [col].shared3.hash */
					Col [col].hash( (int) hash ) ;
					ASSERT (COL_IS_ALIVE (Col, col)) ;
				}
			}

			/* The approximate external column degree is now computed.  */

			/* === Supercolumn detection ======================================== */

			DEBUG3 ("** Supercolumn detection phase. **\n") ;

			if (!NDEBUG)
			{
				detect_super_cols (n_col, Row,
						Col, A, head, pivot_row_start, pivot_row_length) ;
			} else {
				detect_super_cols (Col, A, head, pivot_row_start, pivot_row_length) ;
			}

			/* === Kill the pivotal column ====================================== */

			KILL_PRINCIPAL_COL (Col, pivot_col) ;

			/* === Clear mark =================================================== */

			tag_mark = clear_mark (tag_mark+max_deg+1, max_mark, n_row, Row) ;

			if (!NDEBUG)
			{
				DEBUG3 ("check3\n") ;
				debug_mark (n_row, Row, tag_mark, max_mark) ;
			} /* NDEBUG */

			/* === Finalize the new pivot row, and column scores ================ */

			DEBUG3 ("** Finalize scores phase. **\n") ;

			/* for each column in pivot row */
			rp = pivot_row_start ;
			/* compact the pivot row */
			new_rp = rp ;
			rp_end = rp + pivot_row_length ;
			while (rp < rp_end)
			{
				col = A [rp++] ;
				/* skip dead columns */
				if (COL_IS_DEAD (Col, col))
				{
					continue ;
				}
				A [new_rp++] = col ;
				/* add new pivot row to column */
				A [Col [col].start + (Col [col].length++)] = pivot_row ;

				/* retrieve score so far and add on pivot row's degree. */
				/* (we wait until here for this in case the pivot */
				/* row's degree was reduced due to mass elimination). */
				cur_score = Col [col].score() + pivot_row_degree ;

				/* calculate the max possible score as the number of */
				/* external columns minus the 'k' value minus the */
				/* columns thickness */
				max_score = n_col - k - Col [col].thickness() ;

				/* make the score the external degree of the union-of-rows */
				cur_score -= Col [col].thickness() ;

				/* make sure score is less or equal than the max score */
				cur_score = MIN (cur_score, max_score) ;
				ASSERT (cur_score >= 0) ;

				/* store updated score */
				Col [col].score(cur_score) ;

				/* === Place column back in degree list ========================= */

				ASSERT (min_score >= 0) ;
				ASSERT (min_score <= n_col) ;
				ASSERT (cur_score >= 0) ;
				ASSERT (cur_score <= n_col) ;
				ASSERT (head [cur_score] >= EMPTY) ;
				next_col = head [cur_score] ;
				Col [col].degree_next(next_col) ;
				Col [col].prev(EMPTY) ;
				if (next_col != EMPTY)
				{
					Col [next_col].prev(col) ;
				}
				head [cur_score] = col ;

				/* see if this score is less than current min */
				min_score = MIN (min_score, cur_score) ;

			}

			if (!NDEBUG)
			{
				debug_deg_lists (n_row, n_col, Row, Col, head,
						min_score, n_col2-k, max_deg) ;
			} /* NDEBUG */

			/* === Resurrect the new pivot row ================================== */

			if (pivot_row_degree > 0)
			{
				/* update pivot row length to reflect any cols that were killed */
				/* during super-col detection and mass elimination */
				Row [pivot_row].start  = pivot_row_start ;
				Row [pivot_row].length = (int) (new_rp - pivot_row_start) ;
				ASSERT (Row [pivot_row].length > 0) ;
				Row [pivot_row].degree(pivot_row_degree) ;
				Row [pivot_row].mark(0) ;
				/* pivot row is no longer dead */

				DEBUG1 ("Resurrect Pivot_row %d deg: %d\n",
						pivot_row, pivot_row_degree) ;
			}
		}

		/* === All principal columns have now been ordered ================== */

		return (ngarbage) ;
	}


/* ========================================================================== */
/* === order_children ======================================================= */
/* ========================================================================== */

	/**
	 * The find_ordering routine has ordered all of the principal columns (the
	 * representatives of the supercolumns).  The non-principal columns have not
	 * yet been ordered.  This routine orders those columns by walking up the
	 * parent tree (a column is a child of the column which absorbed it).  The
	 * final permutation vector is then placed in p [0 ... n_col-1], with p [0]
	 * being the first column, and p [n_col-1] being the last.  It doesn't look
	 * like it at first glance, but be assured that this routine takes time linear
	 * in the number of columns.  Although not immediately obvious, the time
	 * taken by this routine is O (n_col), that is, linear in the number of
	 * columns.  Not user-callable.
	 *
	 * @param n_col number of columns of A
	 * @param Col of size n_col+1
	 * @param p p [0 ... n_col-1] is the column permutation
	 */
	private static void order_children (int n_col, Colamd_Col[] Col, int[] p)
	{
		/* === Local variables ============================================== */

		int i ;          /* loop counter for all columns */
		int c ;          /* column index */
		int parent ;     /* index of column's parent */
		int order ;      /* column's order */

		/* === Order each non-principal column ============================== */

		for (i = 0 ; i < n_col ; i++)
		{
			/* find an un-ordered non-principal column */
			ASSERT (COL_IS_DEAD (Col, i)) ;
			if (!COL_IS_DEAD_PRINCIPAL (Col, i) && Col [i].order() == EMPTY)
			{
				parent = i ;
				/* once found, find its principal parent */
				do
				{
					parent = Col [parent].parent() ;
				} while (!COL_IS_DEAD_PRINCIPAL (Col, parent)) ;

				/* now, order all un-ordered non-principal columns along path */
				/* to this parent.  collapse tree at the same time */
				c = i ;
				/* get order of parent */
				order = Col [parent].order() ;

				do
				{
					ASSERT (Col [c].order() == EMPTY) ;

					/* order this column */
					Col [c].order(order++) ;
					/* collaps tree */
					Col [c].parent(parent) ;

					/* get immediate parent of this column */
					c = Col [c].parent() ;

					/* continue until we hit an ordered column.  There are */
					/* guarranteed not to be anymore unordered columns */
					/* above an ordered column */
				} while (Col [c].order() == EMPTY) ;

				/* re-order the super_col parent to largest order for this group */
				Col [parent].order(order) ;
			}
		}

		/* === Generate the permutation ===================================== */

		for (c = 0 ; c < n_col ; c++)
		{
			p [Col [c].order()] = c ;
		}
	}


/* ========================================================================== */
/* === detect_super_cols ==================================================== */
/* ========================================================================== */

	/**
	 * Detects supercolumns by finding matches between columns in the hash buckets.
	 * Check amongst columns in the set A [row_start ... row_start + row_length-1].
	 * The columns under consideration are currently *not* in the degree lists,
	 * and have already been placed in the hash buckets.
	 *
	 * The hash bucket for columns whose hash function is equal to h is stored
	 * as follows:
	 *
	 * if head [h] is >= 0, then head [h] contains a degree list, so:
	 *
	 * 		head [h] is the first column in degree bucket h.
	 * 		Col [head [h]].headhash gives the first column in hash bucket h.
	 *
	 * otherwise, the degree list is empty, and:
	 *
	 * 		-(head [h] + 2) is the first column in hash bucket h.
	 *
	 * For a column c in a hash bucket, Col [c].shared3.prev is NOT a "previous
	 * column" pointer.  Col [c].shared3.hash is used instead as the hash number
	 * for that column.  The value of Col [c].shared4.hash_next is the next column
	 * in the same hash bucket.
	 *
	 * Assuming no, or "few" hash collisions, the time taken by this routine is
	 * linear in the sum of the sizes (lengths) of each column whose score has
	 * just been computed in the approximate degree computation.
	 * Not user-callable.
	 *
	 * @param Col of size n_col+1
	 * @param A row indices of A
	 * @param head head of degree lists and hash buckets
	 * @param row_start pointer to set of columns to check
	 * @param row_length number of columns to check
	 */
	private static void detect_super_cols (Colamd_Col[] Col, int[] A,
			int[] head, int row_start, int row_length)
	{
		detect_super_cols (0, null, Col, A, head, row_start, row_length) ;
	}

	/**
	 * debug version
	 *
	 * @param n_col number of columns of A
	 * @param Row of size n_row+1
	 * @param Col of size n_col+1
	 * @param A row indices of A
	 * @param head head of degree lists and hash buckets
	 * @param row_start pointer to set of columns to check
	 * @param row_length number of columns to check
	 */
	private static void detect_super_cols (int n_col, Colamd_Row[] Row,
		Colamd_Col[] Col, int[] A, int[] head, int row_start, int row_length)
	{
		/* === Local variables ============================================== */

		int hash ;         /* hash value for a column */
		int rp ;           /* pointer to a row */
		int c ;            /* a column index */
		int super_c ;      /* column index of the column to absorb into */
		int cp1 ;          /* column pointer for column super_c */
		int cp2 ;          /* column pointer for column c */
		int length ;       /* length of column super_c */
		int prev_c ;       /* column preceding c in hash bucket */
		int i ;            /* loop counter */
		int rp_end ;       /* pointer to the end of the row */
		int col ;          /* a column index in the row to check */
		int head_column ;  /* first column in hash bucket or degree list */
		int first_col ;    /* first column in hash bucket */

		/* === Consider each column in the row ============================== */

		rp = row_start ;
		rp_end = rp + row_length ;
		while (rp < rp_end)
		{
			col = A [rp++] ;
			if (COL_IS_DEAD (Col, col))
			{
				continue ;
			}

			/* get hash number for this column */
			hash = Col [col].hash() ;
			ASSERT (hash <= n_col) ;

			/* === Get the first column in this hash bucket ===================== */

			head_column = head [hash] ;
			if (head_column > EMPTY)
			{
				first_col = Col [head_column].headhash() ;
			}
			else
			{
				first_col = - (head_column + 2) ;
			}

			/* === Consider each column in the hash bucket ====================== */

			for (super_c = first_col ; super_c != EMPTY ;
					super_c = Col [super_c].hash_next())
			{
				ASSERT (COL_IS_ALIVE (Col, super_c)) ;
				ASSERT (Col [super_c].hash() == hash) ;
				length = Col [super_c].length ;

				/* prev_c is the column preceding column c in the hash bucket */
				prev_c = super_c ;

				/* === Compare super_c with all columns after it ================ */

				for (c = Col [super_c].hash_next() ;
						c != EMPTY ; c = Col [c].hash_next())
				{
					ASSERT (c != super_c) ;
					ASSERT (COL_IS_ALIVE (Col, c)) ;
					ASSERT (Col [c].hash() == hash) ;

					/* not identical if lengths or scores are different */
					if (Col [c].length != length ||
							Col [c].score() != Col [super_c].score())
					{
						prev_c = c ;
						continue ;
					}

					/* compare the two columns */
					cp1 = Col [super_c].start ;
					cp2 = Col [c].start ;

					for (i = 0 ; i < length ; i++)
					{
						/* the columns are "clean" (no dead rows) */
						ASSERT (ROW_IS_ALIVE (Row, A [cp1]))  ;
						ASSERT (ROW_IS_ALIVE (Row, A [cp2]))  ;
						/* row indices will same order for both supercols, */
						/* no gather scatter nessasary */
						if (A [cp1++] != A [cp2++])
						{
							break ;
						}
					}

					/* the two columns are different if the for-loop "broke" */
					if (i != length)
					{
						prev_c = c ;
						continue ;
					}

					/* === Got it!  two columns are identical =================== */

					ASSERT (Col [c].score() == Col [super_c].score()) ;

					Col [super_c].thickness( Col [super_c].thickness() + Col [c].thickness() ) ;
					Col [c].parent(super_c) ;
					KILL_NON_PRINCIPAL_COL (Col, c) ;
					/* order c later, in order_children() */
					Col [c].order(EMPTY) ;
					/* remove c from hash bucket */
					Col [prev_c].hash_next(Col [c].hash_next()) ;
				}
			}

			/* === Empty this hash bucket ======================================= */

			if (head_column > EMPTY)
			{
				/* corresponding degree list "hash" is not empty */
				Col [head_column].headhash(EMPTY) ;
			}
			else
			{
				/* corresponding degree list "hash" is empty */
				head [hash] = EMPTY ;
			}
		}
	}


/* ========================================================================== */
/* === garbage_collection =================================================== */
/* ========================================================================== */

	/**
	 * Defragments and compacts columns and rows in the workspace A.  Used when
	 * all avaliable memory has been used while performing row merging.  Returns
	 * the index of the first free position in A, after garbage collection.  The
	 * time taken by this routine is linear in the size of the array A, which is
	 * itself linear in the number of nonzeros in the input matrix.
	 * Not user-callable.
	 *
	 * @param n_row number of rows
	 * @param n_col number of columns
	 * @param Row row info
	 * @param Col column info
	 * @param A A [0 ... Alen-1] holds the matrix
	 * @param pfree
	 * @return the new value of pfree
	 */
	private static int garbage_collection (int n_row, int n_col,
			Colamd_Row[] Row, Colamd_Col[] Col, int[] A, int pfree)
	{
		/* === Local variables ============================================== */

		int psrc ;    /* source pointer */
		int pdest ;   /* destination pointer */
		int j ;       /* counter */
		int r ;       /* a row index */
		int c ;       /* a column index */
		int length ;  /* length of a row or column */

		int debug_rows = 0 ;
		if (!NDEBUG)
		{
			DEBUG2 ("Defrag..\n") ;
			for (psrc = 0 ; psrc < pfree ; psrc++) ASSERT (A [psrc] >= 0) ;
			debug_rows = 0 ;
		}

		/* === Defragment the columns ======================================= */

		pdest = 0 ;
		for (c = 0 ; c < n_col ; c++)
		{
			if (COL_IS_ALIVE (Col, c))
			{
				psrc = Col [c].start ;

				/* move and compact the column */
				ASSERT (pdest <= psrc) ;
				Col [c].start = (int) (pdest - 0) ;
				length = Col [c].length ;
				for (j = 0 ; j < length ; j++)
				{
					r = A [psrc++] ;
					if (ROW_IS_ALIVE (Row, r))
					{
						A [pdest] = r ;
					}
				}
				Col [c].length = (int) (pdest - Col [c].start) ;
			}
		}

		/* === Prepare to defragment the rows =============================== */

		for (r = 0 ; r < n_row ; r++)
		{
			if (ROW_IS_DEAD (Row, r) || (Row [r].length == 0))
			{
				/* This row is already dead, or is of zero length.  Cannot compact
				 * a row of zero length, so kill it.  NOTE: in the current version,
				 * there are no zero-length live rows.  Kill the row (for the first
				 * time, or again) just to be safe. */
				KILL_ROW (Row, r) ;
			}
			else
			{
				/* save first column index in Row [r].shared2.first_column */
				psrc = Row [r].start ;
				Row [r].first_column(A [psrc]) ;
				ASSERT (ROW_IS_ALIVE (Row, r)) ;
				/* flag the start of the row with the one's complement of row */
				A [psrc] = ONES_COMPLEMENT (r) ;
				if (!NDEBUG)
				{
					debug_rows++ ;
				} /* NDEBUG */
			}
		}

		/* === Defragment the rows ========================================== */

		psrc = pdest ;
		while (psrc < pfree)
		{
			/* find a negative number ... the start of a row */
			if (A [psrc++] < 0)
			{
				psrc-- ;
				/* get the row index */
				r = ONES_COMPLEMENT (A [psrc]) ;
				ASSERT (r >= 0 && r < n_row) ;
				/* restore first column index */
				A [psrc] = Row [r].first_column() ;
				ASSERT (ROW_IS_ALIVE (Row, r)) ;
				ASSERT (Row [r].length > 0) ;
				/* move and compact the row */
				ASSERT (pdest <= psrc) ;
				Row [r].start = (int) (pdest - 0) ;
				length = Row [r].length ;
				for (j = 0 ; j < length ; j++)
				{
					c = A [psrc++] ;
					if (COL_IS_ALIVE (Col, c))
					{
						A [pdest++] = c ;
					}
				}
				Row [r].length = (int) (pdest - Row [r].start) ;
				ASSERT (Row [r].length > 0) ;
				if (!NDEBUG)
				{
					debug_rows-- ;
				} /* NDEBUG */
			}
		}
		/* ensure we found all the rows */
		ASSERT (debug_rows == 0) ;

		/* === Return the new value of pfree ================================ */

		return ((int) (pdest - 0)) ;
	}


/* ========================================================================== */
/* === clear_mark =========================================================== */
/* ========================================================================== */

	/**
	 * Clears the Row [].shared2.mark array, and returns the new tag_mark.
	 *
	 * @param tag_mark new value of tag_mark
	 * @param max_mark max allowed value of tag_mark
	 * @param n_row number of rows in A
	 * @param Row Row [0 ... n_row-1].shared2.mark is set to zero
	 * @return the new value for tag_mark
	 */
	private static int clear_mark (int tag_mark, int max_mark, int n_row,
			Colamd_Row[] Row)
	{
		int r ;

		if (tag_mark <= 0 || tag_mark >= max_mark)
		{
			for (r = 0 ; r < n_row ; r++)
			{
				if (ROW_IS_ALIVE (Row, r))
				{
					Row [r].mark(0) ;
				}
			}
			tag_mark = 1 ;
		}

		return (tag_mark) ;
	}


/* ========================================================================== */
/* === print_report ========================================================= */
/* ========================================================================== */

	/**
	 *
	 * @param method
	 * @param stats
	 */
	private static void print_report (String method, int[] stats)
	{
		int i1, i2, i3 ;

		PRINTF ("\n%s version %d.%d, %s: ", method,
				COLAMD_MAIN_VERSION, COLAMD_SUB_VERSION, COLAMD_DATE) ;

		if (stats == null)
		{
			PRINTF ("No statistics available.\n") ;
			return ;
		}

		i1 = stats [COLAMD_INFO1] ;
		i2 = stats [COLAMD_INFO2] ;
		i3 = stats [COLAMD_INFO3] ;

		if (stats [COLAMD_STATUS] >= 0)
		{
			PRINTF ("OK.  ") ;
		}
		else
		{
			PRINTF ("ERROR.  ") ;
		}

		switch (stats [COLAMD_STATUS])
		{

			case COLAMD_OK_BUT_JUMBLED:

				PRINTF("Matrix has unsorted or duplicate row indices.\n") ;

				PRINTF("%s: number of duplicate or out-of-order row indices: %d\n",
						method, i3) ;

				PRINTF("%s: last seen duplicate or out-of-order row index:   %d\n",
						method, INDEX (i2)) ;

				PRINTF("%s: last seen in column:                             %d",
						method, INDEX (i1)) ;

				/* no break - fall through to next case instead */

			case COLAMD_OK:

				PRINTF("\n") ;

				PRINTF("%s: number of dense or empty rows ignored:           %d\n",
						method, stats [COLAMD_DENSE_ROW]) ;

				PRINTF("%s: number of dense or empty columns ignored:        %d\n",
						method, stats [COLAMD_DENSE_COL]) ;

				PRINTF("%s: number of garbage collections performed:         %d\n",
						method, stats [COLAMD_DEFRAG_COUNT]) ;
				break ;

			case COLAMD_ERROR_A_not_present:

				PRINTF("Array A (row indices of matrix) not present.\n") ;
				break ;

			case COLAMD_ERROR_p_not_present:

				PRINTF("Array p (column pointers for matrix) not present.\n") ;
				break ;

			case COLAMD_ERROR_nrow_negative:

				PRINTF("Invalid number of rows (%d).\n", i1) ;
				break ;

			case COLAMD_ERROR_ncol_negative:

				PRINTF("Invalid number of columns (%d).\n", i1) ;
				break ;

			case COLAMD_ERROR_nnz_negative:

				PRINTF("Invalid number of nonzero entries (%d).\n", i1) ;
				break ;

			case COLAMD_ERROR_p0_nonzero:

				PRINTF("Invalid column pointer, p [0] = %d, must be zero.\n", i1);
				break ;

			case COLAMD_ERROR_A_too_small:

				PRINTF("Array A too small.\n") ;
				PRINTF("        Need Alen >= %d, but given only Alen = %d.\n",
						i1, i2) ;
				break ;

			case COLAMD_ERROR_col_length_negative:

				PRINTF
				("Column %d has a negative number of nonzero entries (%d).\n",
						INDEX (i1), i2) ;
				break ;

			case COLAMD_ERROR_row_index_out_of_bounds:

				PRINTF
				("Row index (row %d) out of bounds (%d to %d) in column %d.\n",
						INDEX (i2), INDEX (0), INDEX (i3-1), INDEX (i1)) ;
				break ;

			case COLAMD_ERROR_out_of_memory:

				PRINTF("Out of memory.\n") ;
				break ;

			/* v2.4: internal-error case deleted */
		}
	}


/* ========================================================================== */
/* === colamd debugging routines ============================================ */
/* ========================================================================== */

	/**
	 * At this point, all empty rows and columns are dead.  All live columns
	 * are "clean" (containing no dead rows) and simplicial (no supercolumns
	 * yet).  Rows may contain dead columns, but all live rows contain at
	 * least one live column.
	 *
	 * @param n_row
	 * @param n_col
	 * @param Row
	 * @param Col
	 * @param A
	 * @param n_col2
	 */
	private static void debug_structures (int n_row, int n_col,
			Colamd_Row[] Row, Colamd_Col[] Col, int[] A, int n_col2)
	{
		if (!NDEBUG)
		{
			/* === Local variables ============================================== */

			int i ;
			int c ;
			int cp ;
			int cp_end ;
			int len ;
			int score ;
			int r ;
			int rp ;
			int rp_end ;
			int deg ;

			/* === Check A, Row, and Col ======================================== */

			for (c = 0 ; c < n_col ; c++)
			{
				if (COL_IS_ALIVE (Col, c))
				{
					len = Col [c].length ;
					score = Col [c].score() ;
					DEBUG4 ("initial live col %5d %5d %5d\n", c, len, score) ;
					ASSERT (len > 0) ;
					ASSERT (score >= 0) ;
					ASSERT (Col [c].thickness() == 1) ;
					cp = Col [c].start ;
					cp_end = cp + len ;
					while (cp < cp_end)
					{
						r = A [cp++] ;
						ASSERT (ROW_IS_ALIVE (Row, r)) ;
					}
				}
				else
				{
					i = Col [c].order() ;
					ASSERT (i >= n_col2 && i < n_col) ;
				}
			}

			for (r = 0 ; r < n_row ; r++)
			{
				if (ROW_IS_ALIVE (Row, r))
				{
					i = 0 ;
					len = Row [r].length ;
					deg = Row [r].degree() ;
					ASSERT (len > 0) ;
					ASSERT (deg > 0) ;
					rp = Row [r].start ;
					rp_end = rp + len ;
					while (rp < rp_end)
					{
						c = A [rp++] ;
						if (COL_IS_ALIVE (Col, c))
						{
							i++ ;
						}
					}
					ASSERT (i > 0) ;
				}
			}
		}
	}


/* ========================================================================== */
/* === debug_deg_lists ====================================================== */
/* ========================================================================== */

	/**
	 * Prints the contents of the degree lists.  Counts the number of columns
	 * in the degree list and compares it to the total it should have.  Also
	 * checks the row degrees.
	 *
	 * @param n_row
	 * @param n_col
	 * @param Row
	 * @param Col
	 * @param head
	 * @param min_score
	 * @param should
	 * @param max_deg
	 */
	private static void debug_deg_lists (int n_row, int n_col,
			Colamd_Row[] Row, Colamd_Col[] Col, int[] head, int min_score,
			int should, int max_deg)
	{
		if (!NDEBUG)
		{
			/* === Local variables ============================================== */

			int deg ;
			int col ;
			int have ;
			int row ;

			/* === Check the degree lists ======================================= */

			if (n_col > 10000 && colamd_debug <= 0)
			{
				return ;
			}
			have = 0 ;
			DEBUG4 ("Degree lists: %d\n", min_score) ;
			for (deg = 0 ; deg <= n_col ; deg++)
			{
				col = head [deg] ;
				if (col == EMPTY)
				{
					continue ;
				}
				DEBUG4 ("%d:", deg) ;
				while (col != EMPTY)
				{
					DEBUG4 (" %d", col) ;
					have += Col [col].thickness() ;
					ASSERT (COL_IS_ALIVE (Col, col)) ;
					col = Col [col].degree_next() ;
				}
				DEBUG4 ("\n") ;
			}
			DEBUG4 ("should %d have %d\n", should, have) ;
			ASSERT (should == have) ;

			/* === Check the row degrees ======================================== */

			if (n_row > 10000 && colamd_debug <= 0)
			{
				return ;
			}
			for (row = 0 ; row < n_row ; row++)
			{
				if (ROW_IS_ALIVE (Row, row))
				{
					ASSERT (Row [row].degree() <= max_deg) ;
				}
			}
		}
	}


/* ========================================================================== */
/* === debug_mark =========================================================== */
/* ========================================================================== */

	/**
	 * Ensures that the tag_mark is less that the maximum and also ensures that
	 * each entry in the mark array is less than the tag mark.
	 *
	 * @param n_row
	 * @param Row
	 * @param tag_mark
	 * @param max_mark
	 */
	private static void debug_mark (int n_row, Colamd_Row[] Row, int tag_mark,
			int max_mark)
	{
		if (!NDEBUG)
		{
			/* === Local variables ============================================== */

			int r ;

			/* === Check the Row marks ========================================== */

			ASSERT (tag_mark > 0 && tag_mark <= max_mark) ;
			if (n_row > 10000 && colamd_debug <= 0)
			{
				return ;
			}
			for (r = 0 ; r < n_row ; r++)
			{
				ASSERT (Row [r].mark() < tag_mark) ;
			}
		}
	}


/* ========================================================================== */
/* === debug_matrix ========================================================= */
/* ========================================================================== */

	/**
	 * Prints out the contents of the columns and the rows.
	 *
	 * @param n_row
	 * @param n_col
	 * @param Row
	 * @param Col
	 * @param A
	 */
	private static void debug_matrix (int n_row, int n_col,
			Colamd_Row[] Row, Colamd_Col[] Col, int[] A)
	{
		if (!NDEBUG)
		{
			/* === Local variables ============================================== */

			int r ;
			int c ;
			int rp ;
			int rp_end ;
			int cp ;
			int cp_end ;

			/* === Dump the rows and columns of the matrix ====================== */

			if (colamd_debug < 3)
			{
				return ;
			}
			DEBUG3 ("DUMP MATRIX:\n") ;
			for (r = 0 ; r < n_row ; r++)
			{
				DEBUG3 ("Row %d alive? %d\n", r, ROW_IS_ALIVE (Row, r) ? 1 : 0) ;
				if (ROW_IS_DEAD (Row, r))
				{
					continue ;
				}
				DEBUG3 ("start %d length %d degree %d\n",
						Row [r].start, Row [r].length, Row [r].degree()) ;
				rp = Row [r].start ;
				rp_end = rp + Row [r].length ;
				while (rp < rp_end)
				{
					c = A [rp++] ;
					DEBUG4 ("	%d col %d\n", COL_IS_ALIVE (Col, c) ? 1 : 0, c) ;
				}
			}

			for (c = 0 ; c < n_col ; c++)
			{
				DEBUG3 ("Col %d alive? %d\n", c, COL_IS_ALIVE (Col, c) ? 1 : 0) ;
				if (COL_IS_DEAD (Col, c))
				{
					continue ;
				}
				DEBUG3 ("start %d length %d shared1 %d shared2 %d\n",
						Col [c].start, Col [c].length,
						Col [c].thickness(), Col [c].score()) ;
				cp = Col [c].start ;
				cp_end = cp + Col [c].length ;
				while (cp < cp_end)
				{
					r = A [cp++] ;
					DEBUG4 ("	%d row %d\n", ROW_IS_ALIVE (Row, r) ? 1 : 0, r) ;
				}
			}
		}
	}

	@SuppressWarnings("unused")
	private static void colamd_get_debug (String method)
	{
		if (!NDEBUG)
		{
			File f ;
			colamd_debug = 0 ;  /* no debug printing */
			f = new File("debug") ;
			if (!f.exists())
			{
				colamd_debug = 0 ;
			}
			else
			{
				try {
					FileReader fr ;
					fr = new FileReader(f) ;
					BufferedReader br ;
					br = new BufferedReader(fr) ;
					colamd_debug = Integer.valueOf( br.readLine() ) ;
					br.close() ;
					fr.close() ;
				} catch (IOException e) {
					PRINTF ("%s: AMD_debug_init, " +
							"error reading debug.amd file", method) ;
				}
			}
			DEBUG0 ("%s: debug version, D = %d (THIS WILL BE SLOW!)\n",
					method, colamd_debug) ;
		}
	}

}
