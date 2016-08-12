
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "invert_matrix.h"

int invert_matrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse)
{
    using namespace boost::numeric::ublas;
	//copy of the input
	matrix<double> A(input);
	permutation_matrix<std::size_t> pm(A.size1());

	//LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0) return false;

	inverse.assign(identity_matrix<double> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);
	return 0;
}
