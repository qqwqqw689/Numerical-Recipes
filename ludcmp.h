#include "nr3.h"

struct LUdcmp
// object for solving linear equations A*x=b using LU decomposition, and related functions.
{
	Int n;
	MatDoub lu; // stores the decomposition
	VecInt indx; // stores the row permutation effected by the parital pivoting.
	Doub d; // used by det
	// d is output as +/-1 depending on whether the number of row interchanges
	// was even or odd;
	LUdcmp(MatDoub_I &a); // constructor.argument is the matrix A.
	void solve(VecDoub_I &b, VecDoub_O &x); // slove for a single right-hand side.
	void solve(MatDoub_I &b, MatDoub_O &x); // solve for multiple right-hand sides
	void inverse(MatDoub_O &ainv); // calculate matrix inverse
	Doub det(); // return the determinant of A
	void mprove(VecDoub_I &b, VecDoub_IO &x); // 
	MatDoub_I &aref; // used only by mprove
};
LUdcmp::LUdcmp(MatDoub_I &a) : n(a.nrows()), lu(a), aref(a), indx(n) {
	// a : a matrix a[0..n-1][0..n-1];
	// this routine replaces it by the LU decomposition of a rowwise permutation of itself;
	const Doub TINY=1.0e-40; // a small number
	Int i,imax,j,k;
	Doub big,temp;
	VecDoub vv(n); // vv stores the implicit scaling of each row.
	d=1.0; // no row interchanges yet.
	// one additional wrinkle:
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=abs(lu[i][j])) > big) big=temp; // the largest element in each row
		if (big == 0.0) throw("Singular matrix in LUdcmp");
		vv[i]=1.0/big;
	}
	for (k=0;k<n;k++) {
		big=0.0;
		imax=k;
		for (i=k;i<n;i++) {
			temp=vv[i]*abs(lu[i][k]);
			if (temp > big) {
				big=temp;
				imax=i;
			}
		}
		if (k != imax) {
			for (j=0;j<n;j++) {
				temp=lu[imax][j];
				lu[imax][j]=lu[k][j];
				lu[k][j]=temp;
			}
			d = -d; // change the parity of d
			vv[imax]=vv[k]; // 
		}
		indx[k]=imax;
		if (lu[k][k] == 0.0) lu[k][k]=TINY;
		for (i=k+1;i<n;i++) {
			temp=lu[i][k] /= lu[k][k];
			for (j=k+1;j<n;j++)
				lu[i][j] -= temp*lu[k][j];
		}
	}
}

void LUdcmp::solve(VecDoub_I &b, VecDoub_O &x)
{
	// sloves the set of n linear equation A*x = b.
	// b[0..n-1] is input as the right-hand side vector b.
	// x returns the solution vector x.
	Int i,ii=0,ip,j;
	Doub sum;
	if (b.size() != n || x.size() != n)
		throw("LUdcmp::solve bad sizes");
	for (i=0;i<n;i++) x[i] = b[i];
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0) // when ii is set to a positive value, it will become the index of the first nonvanishing element of b.
			for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
		else if (sum != 0.0) // a nonzero element was encountered,so from now on we will
			ii=i+1;
		x[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=x[i];
		for (j=i+1;j<n;j++) sum -= lu[i][j]*x[j];
		x[i]=sum/lu[i][i]; // store a component of the solution vector X.
	}
}

void LUdcmp::solve(MatDoub_I &b, MatDoub_O &x)
{
	// solves m sets of n linear equations A*X=B.
	// the matrix b[0..n-1][0..m-1] inputs the right-hand sides.
	// x[0..n-1][0..m-1] returns the solution.
	int i,j,m=b.ncols();
	if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
		throw("LUdcmp::solve bad sizes");
	VecDoub xx(n);
	for (j=0;j<m;j++) {
		for (i=0;i<n;i++) xx[i] = b[i][j];
		solve(xx,xx);
		for (i=0;i<n;i++) x[i][j] = xx[i];
	}
}
void LUdcmp::inverse(MatDoub_O &ainv)
{
	// using the stored LU decomposition, return in ainv the matrix inverse.
	Int i,j;
	ainv.resize(n,n);
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) ainv[i][j] = 0.;
		ainv[i][i] = 1.;
	}
	solve(ainv,ainv);
}
Doub LUdcmp::det()
{
	Doub dd = d;
	for (Int i=0;i<n;i++) dd *= lu[i][i];
	return dd;
}
void LUdcmp::mprove(VecDoub_I &b, VecDoub_IO &x)
{
	Int i,j;
	VecDoub r(n);
	for (i=0;i<n;i++) {
		Ldoub sdp = -b[i];
		for (j=0;j<n;j++)
			sdp += (Ldoub)aref[i][j] * (Ldoub)x[j];
		r[i]=sdp;
	}
	solve(r,r);
	for (i=0;i<n;i++) x[i] -= r[i];
}