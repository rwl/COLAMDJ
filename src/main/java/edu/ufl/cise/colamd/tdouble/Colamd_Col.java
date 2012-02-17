package edu.ufl.cise.colamd.tdouble;

public class Colamd_Col {

	/* index for A of first row in this column, or DEAD */
	/* if column is dead */
	public int start ;

	/* number of rows in this column */
	public int length ;

	private int shared1 ;
	private int shared2 ;
	private int shared3 ;
	private int shared4 ;

	public Colamd_Col () {

	}

	/* number of original columns represented by this */
	/* col, if the column is alive */
	public int thickness () {
		return shared1 ;
	}
	public void thickness (int thickness) {
		shared1 = thickness ;
	}

	/* parent in parent tree super-column structure, if */
	/* the column is dead */
	public int parent () {
		return shared1 ;
	}

	public void parent (int parent) {
		shared1 = parent ;
	}


	/* the score used to maintain heap, if col is alive */
	public int score () {
		return shared2 ;
	}

	public void score (int score) {
		shared2 = score ;
	}

	/* pivot ordering of this column, if col is dead */
	public int order () {
		return shared2 ;
	}

	public void order (int order) {
		shared2 = order ;
	}


	/* head of a hash bucket, if col is at the head of */
	/* a degree list */
	public int headhash () {
		return shared3 ;
	}

	public void headhash (int headhash) {
		shared3 = headhash ;
	}

	/* hash value, if col is not in a degree list */
	public int hash () {
		return shared3 ;
	}

	public void hash (int hash) {
		shared3 = hash ;
	}

	/* previous column in degree list, if col is in a */
	/* degree list (but not at the head of a degree list) */
	public int prev () {
		return shared3 ;
	}

	public void prev (int prev) {
		shared3 = prev ;
	}


	/* next column, if col is in a degree list */
	public int degree_next () {
		return shared4 ;
	}

	public void degree_next (int degree_next) {
		shared4 = degree_next ;
	}

	/* next column, if col is in a hash list */
	public int hash_next () {
		return shared4 ;
	}

	public void hash_next (int hash_next) {
		shared4 = hash_next ;
	}

}
