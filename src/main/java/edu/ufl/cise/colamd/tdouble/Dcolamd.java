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
	private static int t_add (int a, int b, int ok)
	{
	    (ok) = (ok != 0) && ((a + b) >= MAX (a,b)) ? 1 : 0;
	    return ((ok != 0) ? (a + b) : 0) ;
	}

	/**
	 * compute a*k where k is a small integer, and check for integer overflow
	 */
	static int t_mult (int a, int k, int ok)
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
	private static int COLAMD_C(int n_col, int ok)
	{
		return ((t_mult (t_add (n_col, 1, ok), sizeof (Colamd_Col), ok) / sizeof (int))) ;
	}

	private static int COLAMD_R(int n_row, int ok)
	{
	    return ((t_mult (t_add (n_row, 1, ok), sizeof (Colamd_Row), ok) / sizeof (int))) ;
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
	    int s, c, r ;
	    int ok = TRUE ;
	    if (nnz < 0 || n_row < 0 || n_col < 0)
	    {
		return (0) ;
	    }
	    s = t_mult (nnz, 2, ok) ;	    /* 2*nnz */
	    c = COLAMD_C (n_col, ok) ;	    /* size of column structures */
	    r = COLAMD_R (n_row, ok) ;	    /* size of row structures */
	    s = t_add (s, c, ok) ;
	    s = t_add (s, r, ok) ;
	    s = t_add (s, n_col, ok) ;	    /* elbow room */
	    s = t_add (s, nnz/5, ok) ;	    /* elbow room */
	    ok = ok != 0 && (s < Int_MAX) ? 1 : 0;
	    return (ok != 0 ? s : 0) ;
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
	    /* === Local variables ================================================== */

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
	    double[] knobs, int[] stats, Object allocate, Object release)
	{
		/* === Local variables ================================================== */

	    int[] count ;		/* length of each column of M, and col pointer*/
	    int[] mark ;			/* mark array for finding duplicate entries */
	    int[] M ;			/* row indices of matrix M */
	    int Mlen ;		/* length of M */
	    int n_row ;			/* number of rows in M */
	    int nnz ;			/* number of entries in A */
	    int i ;			/* row index of A */
	    int j ;			/* column index of A */
	    int k ;			/* row index of M */
	    int mnz ;			/* number of nonzeros in M */
	    int pp ;			/* index into a column of A */
	    int last_row ;		/* last row seen in the current column */
	    int length ;		/* number of nonzeros in a column */

	    double[] cknobs = new double[COLAMD_KNOBS] ;		/* knobs for colamd */
	    double[] default_knobs = new double[COLAMD_KNOBS] ;	/* default knobs for colamd */

	    if (!NDEBUG)
	    {
	    	colamd_get_debug ("symamd") ;
	    }

	    /* === Check the input arguments ======================================== */

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

	    /* === If no knobs, set default knobs =================================== */

	    if (knobs == null)
	    {
		COLAMD_set_defaults (default_knobs) ;
		knobs = default_knobs ;
	    }

	    /* === Allocate count and mark ========================================== */

	    count = (int[]) ((allocate) (n+1, sizeof (int))) ;
	    if (count == null)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
		DEBUG0 ("symamd: allocate count (size %d) failed\n", n+1) ;
		return (FALSE) ;
	    }

	    mark = (int[]) ((allocate) (n+1, sizeof (int))) ;
	    if (!mark)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
		(release) ((void[]) count) ;
		DEBUG0 ("symamd: allocate mark (size %d) failed\n", n+1) ;
		return (FALSE) ;
	    }

	    /* === Compute column counts of M, check if A is valid ================== */

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
		    (*release) ((void *) count) ;
		    (*release) ((void *) mark) ;
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
			(*release) ((void *) count) ;
			(*release) ((void *) mark) ;
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

	    /* v2.4: removed free(mark) */

	    /* === Compute column pointers of M ===================================== */

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

	    /* === Construct M ====================================================== */

	    mnz = perm [n] ;
	    n_row = mnz / 2 ;
	    Mlen = COLAMD_recommended (mnz, n_row, n) ;
	    M = (int[]) ((allocate) (Mlen, sizeof (int))) ;
	    DEBUG0 ("symamd: M is %d-by-%d with %d entries, Mlen = %g\n",
	    	n_row, n, mnz, (double) Mlen) ;

	    if (!M)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
		(*release) ((void *) count) ;
		(*release) ((void *) mark) ;
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
		/* v2.4: free(mark) moved below */
	    }

	    /* count and mark no longer needed */
	    (*release) ((void *) count) ;
	    (*release) ((void *) mark) ;	/* v2.4: free (mark) moved here */
	    ASSERT (k == n_row) ;

	    /* === Adjust the knobs for M =========================================== */

	    for (i = 0 ; i < COLAMD_KNOBS ; i++)
	    {
		cknobs [i] = knobs [i] ;
	    }

	    /* there are no dense rows in M */
	    cknobs [COLAMD_DENSE_ROW] = -1 ;
	    cknobs [COLAMD_DENSE_COL] = knobs [COLAMD_DENSE_ROW] ;

	    /* === Order the columns of M =========================================== */

	    /* v2.4: colamd cannot fail here, so the error check is removed */
	    (void) COLAMD_MAIN (n_row, n, (int) Mlen, M, perm, cknobs, stats) ;

	    /* Note that the output permutation is now in perm */

	    /* === get the statistics for symamd from colamd ======================== */

	    /* a dense column in colamd means a dense row and col in symamd */
	    stats [COLAMD_DENSE_ROW] = stats [COLAMD_DENSE_COL] ;

	    /* === Free M =========================================================== */

	    (*release) ((void *) M) ;
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
		/* === Local variables ================================================== */

	    int i ;			/* loop index */
	    int nnz ;			/* nonzeros in A */
	    int Row_size ;		/* size of Row [], in integers */
	    int Col_size ;		/* size of Col [], in integers */
	    int need ;		/* minimum required length of A */
	    Colamd_Row Row ;		/* pointer into A of Row [0..n_row] array */
	    Colamd_Col Col ;		/* pointer into A of Col [0..n_col] array */
	    int n_col2 ;		/* number of non-dense, non-empty columns */
	    int n_row2 ;		/* number of non-dense, non-empty rows */
	    int ngarbage ;		/* number of garbage collections performed */
	    int max_deg ;		/* maximum row degree */
	    double[] default_knobs = new double[COLAMD_KNOBS] ;	/* default knobs array */
	    int aggressive ;		/* do aggressive absorption */
	    int ok ;

	    if (!NDEBUG)
	    {
	    	colamd_get_debug ("colamd") ;
	    }

	    /* === Check the input arguments ======================================== */

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

	    /* === If no knobs, set default knobs =================================== */

	    if (knobs == null)
	    {
		COLAMD_set_defaults (default_knobs) ;
		knobs = default_knobs ;
	    }

	    aggressive = (knobs [COLAMD_AGGRESSIVE] != FALSE) ? 1 : 0;

	    /* === Allocate the Row and Col arrays from array A ===================== */

	    ok = TRUE ;
	    Col_size = COLAMD_C (n_col, ok) ;	    /* size of Col array of structs */
	    Row_size = COLAMD_R (n_row, ok) ;	    /* size of Row array of structs */

	    /* need = 2*nnz + n_col + Col_size + Row_size ; */
	    need = t_mult (nnz, 2, ok) ;
	    need = t_add (need, n_col, ok) ;
	    need = t_add (need, Col_size, ok) ;
	    need = t_add (need, Row_size, ok) ;

	    if (!(ok != 0) || need > (int) Alen || need > Int_MAX)
	    {
		/* not enough space in array A to perform the ordering */
		stats [COLAMD_STATUS] = COLAMD_ERROR_A_too_small ;
		stats [COLAMD_INFO1] = need ;
		stats [COLAMD_INFO2] = Alen ;
		DEBUG0 ("colamd: Need Alen >= %d, given only Alen = %d\n", need,Alen);
		return (FALSE) ;
	    }

	    Alen -= Col_size + Row_size ;
	    Col = (Colamd_Col) A [Alen] ;
	    Row = (Colamd_Row) A [Alen + Col_size] ;

	    /* === Construct the row and column data structures ===================== */

	    if (!init_rows_cols (n_row, n_col, Row, Col, A, p, stats))
	    {
		/* input matrix is invalid */
		DEBUG0 ("colamd: Matrix invalid\n") ;
		return (FALSE) ;
	    }

	    /* === Initialize scores, kill dense rows/columns ======================= */

	    init_scoring (n_row, n_col, Row, Col, A, p, knobs,
		n_row2, n_col2, max_deg) ;

	    /* === Order the supercolumns =========================================== */

	    ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
		n_col2, max_deg, 2*nnz, aggressive) ;

	    /* === Order the non-principal columns ================================== */

	    order_children (n_col, Col, p) ;

	    /* === Return statistics in stats ======================================= */

	    stats [COLAMD_DENSE_ROW] = n_row - n_row2 ;
	    stats [COLAMD_DENSE_COL] = n_col - n_col2 ;
	    stats [COLAMD_DEFRAG_COUNT] = ngarbage ;
	    DEBUG0 (("colamd: done.\n")) ;
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

}
