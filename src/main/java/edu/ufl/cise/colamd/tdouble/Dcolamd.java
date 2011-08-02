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
	public static int COLAMD_INFO1 = 4 ;
	public static int COLAMD_INFO2 = 5 ;
	public static int COLAMD_INFO3 = 6 ;

	/* error codes returned in stats [3]: */
	public static int COLAMD_OK				= (0) ;
	public static int COLAMD_OK_BUT_JUMBLED			= (1) ;
	public static int COLAMD_ERROR_A_not_present		= (-1) ;
	public static int COLAMD_ERROR_p_not_present		= (-2) ;
	public static int COLAMD_ERROR_nrow_negative		= (-3) ;
	public static int COLAMD_ERROR_ncol_negative		= (-4) ;
	public static int COLAMD_ERROR_nnz_negative		= (-5) ;
	public static int COLAMD_ERROR_p0_nonzero			= (-6) ;
	public static int COLAMD_ERROR_A_too_small		= (-7) ;
	public static int COLAMD_ERROR_col_length_negative	= (-8) ;
	public static int COLAMD_ERROR_row_index_out_of_bounds	= (-9) ;
	public static int COLAMD_ERROR_out_of_memory		= (-10) ;
	public static int COLAMD_ERROR_internal_error		= (-999) ;

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

	private static String ID = "%d" ;

	private static int Int_MAX = Integer.MAX_VALUE;

	/* ========================================================================== */
	/* === Row and Column structures ============================================ */
	/* ========================================================================== */

	/* User code that makes use of the colamd/symamd routines need not directly */
	/* reference these structures.  They are used only for colamd_recommended. */

	public class Colamd_Col {

		/* index for A of first row in this column, or DEAD */
		/* if column is dead */
	    public int start ;
	    /* number of rows in this column */
	    public int length ;

	    /* number of original columns represented by this */
		/* col, if the column is alive */
	    public int thickness ;
	    /* parent in parent tree super-column structure, if */
		/* the column is dead */
	    public int parent ;

	    /* the score used to maintain heap, if col is alive */
	    public int score ;
	    /* pivot ordering of this column, if col is dead */
	    public int order ;

	    /* head of a hash bucket, if col is at the head of */
		/* a degree list */
	    public int headhash ;
	    /* hash value, if col is not in a degree list */
	    public int hash ;
	    /* previous column in degree list, if col is in a */
		/* degree list (but not at the head of a degree list) */
	    public int prev ;

	    /* next column, if col is in a degree list */
	    public int degree_next ;
	    /* next column, if col is in a hash list */
	    public int hash_next ;

	}

	public class Colamd_Row {

		/* index for A of first col in this row */
		public int start ;
	    /* number of principal columns in this row */
		public int length ;

	    /* number of principal & non-principal columns in row */
		public int degree ;
		/* used as a row pointer in init_rows_cols () */
		public int p ;

		/* for computing set differences and marking dead rows*/
		public int mark ;
		/* first column in row (used in garbage collection) */
		public int first_column ;

	}

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

	private static final double MIN (double a, double b)
	{
		return (((a) < (b)) ? (a) : (b)) ;
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
//	private static ROW_IS_DEAD(r)			ROW_IS_MARKED_DEAD (Row[r].shared2.mark)
//	private static ROW_IS_MARKED_DEAD(row_mark)	(row_mark < ALIVE)
//	private static ROW_IS_ALIVE(r)			(Row [r].shared2.mark >= ALIVE)
//	private static COL_IS_DEAD(c)			(Col [c].start < ALIVE)
//	private static COL_IS_ALIVE(c)			(Col [c].start >= ALIVE)
//	private static COL_IS_DEAD_PRINCIPAL(c)	(Col [c].start == DEAD_PRINCIPAL)
//	private static KILL_ROW(r)			{ Row [r].shared2.mark = DEAD ; }
//	private static KILL_PRINCIPAL_COL(c)		{ Col [c].start = DEAD_PRINCIPAL ; }
//	private static KILL_NON_PRINCIPAL_COL(c)	{ Col [c].start = DEAD_NON_PRINCIPAL ; }

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
	private static int colamd_debug = 0 ;

	protected static void AMD_DEBUG0 (String format, Object... args)
	{
		if (!NDEBUG)
		{
			PRINTF (format, args) ;
		}
	}

	protected static void AMD_DEBUG1(String format, Object... args)
	{
		if (!NDEBUG)
		{
			if (colamd_debug >= 1) PRINTF (format, args) ;
		}
	}

	protected static void AMD_DEBUG2(String format, Object... args)
	{
		if (!NDEBUG)
		{
			if (colamd_debug >= 2) PRINTF (format, args) ;
		}
	}

	protected static void AMD_DEBUG3(String format, Object... args)
	{
		if (!NDEBUG)
		{
			if (colamd_debug >= 3) PRINTF (format, args) ;
		}
	}

	protected static void AMD_DEBUG4(String format, Object... args)
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

}
