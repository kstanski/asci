
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "solver.h"

int invert_matrix(const bnu::matrix<double>& input, bnu::matrix<double>& inverse)
{
    using namespace boost::numeric::ublas;
	//copy of the input
	matrix<double> A(input);
	permutation_matrix<std::size_t> pm(A.size1());

	//LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0) return 1;

	inverse.assign(identity_matrix<double> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);
	return 0;
}

int solve_linear_system(const bnu::matrix<double>& A, bnu::vector<double>& x, bnu::vector<double>& y)
{
    using namespace boost::numeric::ublas;
	matrix<double> A_cp(A);
	vector<double> y_cp (y);
	permutation_matrix<std::size_t> pm(A.size1());
	int res = lu_factorize(A_cp, pm);
	if (res != 0) return 1;
	lu_substitute(A_cp, pm, y_cp);
    x = y_cp;

	return 0;
}
