package edu.ufl.cise.colamd.tdouble;

public class Colamd_Row {

	/* index for A of first col in this row */
	public int start ;
	/* number of principal columns in this row */
	public int length ;

	private int shared1 ;
	private int shared2 ;

	public Colamd_Row() {

	}

	/* number of principal & non-principal columns in row */
	public int degree () {
		return shared1 ;
	}

	public void degree (int degree) {
		shared1 = degree ;
	}

	/* used as a row pointer in init_rows_cols () */
	public int p () {
		return shared1 ;
	}

	public void p (int p) {
		shared1 = p ;
	}


	/* for computing set differences and marking dead rows*/
	public int mark () {
		return shared2 ;
	}

	public void mark (int mark) {
		shared2 = mark ;
	}

	/* first column in row (used in garbage collection) */
	public int first_column () {
		return shared2 ;
	}

	public void first_column (int first_column) {
		shared2 = first_column ;
	}

}
