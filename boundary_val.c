#include "boundary_val.h"
#include "init.h"


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
) {
	int i, j;
	for (i = 1; i <= imax; ++i) {
		/* (15) */
		V[i][0] = 0.0;
		V[i][jmax] = 0.0;
		/* (16) */
		U[i][0] = -U[i][1];
		/* Apply condition for the moving lid */
		U[i][jmax + 1] = 2.0 - U[i][jmax];
	}

	for (j = 1; j <= jmax; ++j) {
		/* (15) */
		U[0][j] = 0.0;
		U[imax][j] = 0.0;
		/* (16) */
		V[0][j] = -V[1][j];
		V[imax + 1][j] = -V[imax][j];
	}
}
