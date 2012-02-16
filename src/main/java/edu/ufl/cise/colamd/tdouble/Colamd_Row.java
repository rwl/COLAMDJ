package edu.ufl.cise.colamd.tdouble;

public class Colamd_Row {

	public Colamd_Row() {

	}

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
