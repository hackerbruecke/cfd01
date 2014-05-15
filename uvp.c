#include "uvp.h"
#include "math.h"

/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
) {
	double cond1, umax, vmax, cfl1, cfl2, mincond;
	int i, j;

	cond1 = Re/2.0 * dx*dx*dy*dy / (dx*dx + dy*dy);

	umax = fabs(U[0][0]);
	vmax = fabs(V[0][0]);
	/* Find dominant values in the array */
	for (i = 0; i <= imax+1; ++i) {
		for (j = 0; j <= jmax+1; ++j) {
			umax = fmax(umax, fabs(U[i][j]));
			vmax = fmax(vmax, fabs(U[i][j]));
		}
	}

	cfl1 = dx/umax;
	cfl2 = dy/vmax;
	mincond = fmin(cond1, fmin(cfl1, cfl2));

	*dt = tau * mincond;
}

/*
 * Calculate approximation for second derivative in dx
 */
static double second_derivative_dx(double **F, int i, int j, double dx)
{
	return (F[i+1][j] - 2.0*F[i][j] + F[i-1][j]) / (dx*dx);
}

/*
 * Calculate approximation for second derivative in dy
 */
static double second_derivative_dy(double **F, int i, int j, double dy)
{
	return (F[i][j+1] - 2.0*F[i][j] + F[i][j-1]) / (dy*dy);
}

/*
 * Calculate derivative of u^2 in dx
 */
static double derivative_u_square_dx(double **U, double dx, double alpha, int i, int j)
{
	double a = (U[i][j] + U[i+1][j]) / 2.0;
	double b = (U[i-1][j] + U[i][j]) / 2.0;
	double c = (U[i][j] - U[i+1][j]) / 2.0;
	double d = (U[i-1][j] - U[i][j]) / 2.0;
	return (a*a - b*b)/dx + alpha/dx * (fabs(a)*c - fabs(b)*d);
}

/*
 * Calculate derivative of v^2 in dy
 */
static double derivative_v_square_dy(double **V, double dy, double alpha, int i, int j)
{
	double a = (V[i][j] + V[i][j+1]) / 2.0;
	double b = (V[i][j-1] + V[i][j]) / 2.0;
	double c = (V[i][j] - V[i][j+1]) / 2.0;
	double d = (V[i][j-1] - V[i][j]) / 2.0;
	return (a*a - b*b)/dy + alpha/dy * (fabs(a)*c - fabs(b)*d);
}

/*
 * Calculate derivative of u*v in dx
 */
static double derivative_uv_dx(double **U, double **V, double dx, double alpha, int i, int j)
{
	double au = (U[i][j] + U[i][j+1]) / 2.0;
	double av = (V[i][j] + V[i+1][j]) / 2.0;
	double bu = (U[i-1][j] + U[i-1][j+1]) / 2.0;
	double bv = (V[i-1][j] + V[i][j]) / 2.0;
	double avm = (V[i][j] - V[i+1][j]) / 2.0;
	double bvm = (V[i-1][j] - V[i][j]) / 2.0;
	return (au*av - bu*bv)/dx + alpha/dx*(fabs(au)*avm - fabs(bu)*bvm);
}

/*
 * Calculate derivative of u*v in dy
 */
static double derivative_uv_dy(double **U, double **V, double dy, double alpha, int i, int j)
{
	double av = (V[i][j] + V[i+1][j]) / 2.0;
	double au = (U[i][j] + U[i][j+1]) / 2.0;
	double bv = (V[i][j-1] + V[i+1][j-1]) / 2.0;
	double bu = (U[i][j-1] + U[i][j]) / 2.0;
	double aum = (U[i][j] - U[i][j+1]) / 2.0;
	double bum = (U[i][j-1] - U[i][j]) / 2.0;
	return (av*au - bv*bu)/dy + alpha/dy*(fabs(av)*aum - fabs(bv)*bum);
}

/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */
void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
) {
	int i, j;
	/* Required for computation of F */
	double d2ux, d2uy, du2x, duvy;
	/* Required for computation of G */
	double duvx, dv2y, d2vy, d2vx;

	/* Computation of F */
	for (i = 1; i < imax; ++i) {
		for (j = 1; j <= jmax; ++j) {
			du2x = derivative_u_square_dx(U, dx, alpha, i, j);
			duvy = derivative_uv_dy(U, V, dy, alpha, i, j);
			/* d^2u/dx^2 */
			d2ux = second_derivative_dx(U, i, j, dx);
			/* d^2u/dy^2 */
			d2uy = second_derivative_dy(U, i, j, dy);
			F[i][j] = U[i][j] + dt * ((d2ux + d2uy)/Re - du2x - duvy + GX);
		}
	}
	/* Computation of G */
	for (i = 1; i <= imax; ++i) {
		for (j = 1; j < jmax; ++j) {
			duvx = derivative_uv_dx(U, V, dx, alpha, i, j);
			dv2y = derivative_v_square_dy(V, dy, alpha, i, j);
			/* d^2v/dx^2 */
			d2vx = second_derivative_dx(V, i, j, dx);
			/* d^2v/dy^2 */
			d2vy = second_derivative_dy(V, i, j, dy);
			G[i][j] = V[i][j] + dt * ((d2vx + d2vy)/Re - duvx - dv2y + GY);
		}
	}

	/* Apply Neumann boundary conditions */
	for (j = 1; j <= jmax; ++j) {
		F[0][j] = U[0][j];
	}
	for (i = 1; i <= imax; ++i) {
		G[i][0] = V[i][0];
	}
}

/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
) {
	int i, j;
	for (i = 1; i <= imax; ++i) {
		for (j = 1; j <= jmax; ++j) {
			RS[i][j] = ((F[i][j] - F[i-1][j])/dx + (G[i][j] - G[i][j-1])/dy) / dt;
		}
	}
}

/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
) {
	int i, j;
	for (i = 1; i < imax; ++i) {
		for (j = 1; j <= jmax; ++j) {
			U[i][j] = F[i][j] - dt/dx*(P[i+1][j] - P[i][j]);
		}
	}

	for (i = 1; i <= imax; ++i) {
		for (j = 1; j < jmax; ++j) {
			V[i][j] = G[i][j] - dt/dy*(P[i][j+1] - P[i][j]);
		}
	}
}

