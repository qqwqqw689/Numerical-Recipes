#include "nr3.h"

void gaussj(MatDoub_IO &a, MatDoub_IO &b)
{
	// the input matrix is a[0..n-1][0..n-1], b[0..n-1][0..m-1].
	// simply picking the largest(in magnitude) available element as pivot
	// is a very good choice.
	// The matrix inverse of A is gradually built up in A
	// as the original A is destoryed.
	Int i,icol,irow,j,k,l,ll,n=a.nrows(),m=b.ncols();
	Doub big,dum,pivinv;
	VecInt indxc(n),indxr(n),ipiv(n);
	// These integer arrays are used for bookkeeping on the pivoting.
	for (j=0;j<n;j++) ipiv[j]=0; 
	for (i=0;i<n;i++) { // // This is the main loop over the columns to be reduced.
		big=0.0;
		for (j=0;j<n;j++) // this is the outer loop of the search for a pivot element.
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (abs(a[j][k]) >= big) {
							big=abs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		// We  now have the pivot element, so we interchange rows, if indeed, to put the pivot
		// element on the diagonal.
		if (irow != icol) {
			for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
		}
		indxr[i]=irow;
		// the row in which that pivot element was originally located.
		indxc[i]=icol;
		// the column of the (i+1)th pivot element, is the (i+1)th the column
		// that is reduced.
		if (a[icol][icol] == 0.0) throw("gaussj: Singular Matrix");
		// We are now ready to divide the pivot row by the pivot element.
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		// Next, we reduce the rows expect for the pivot one, of course.
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	// interchanging pairs of column in the reverse order that the
	// permutation was built up.
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
}

// Overloaded version with no right-hand sides.
// Replaces a by its inverse. 
void gaussj(MatDoub_IO &a)
{
	MatDoub b(a.nrows(),0); // Dummy vector with zero columns.
	gaussj(a,b);
}
