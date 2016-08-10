
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "invert_matrix.h"

 /* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
int invert_matrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse)
{
    using namespace boost::numeric::ublas;

	typedef permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix<double> A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<double> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}
