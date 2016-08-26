#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include "solver.h"

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
