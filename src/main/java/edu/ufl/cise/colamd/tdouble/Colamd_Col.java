package edu.ufl.cise.colamd.tdouble;

public class Colamd_Col {

	public Colamd_Col() {

	}

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
