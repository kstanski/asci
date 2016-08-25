#ifndef INVERT_MATRIX_H_INCLUDED
#define INVERT_MATRIX_H_INCLUDED

#include <boost/numeric/ublas/matrix.hpp>
namespace bnu = boost::numeric::ublas;
int invert_matrix(const bnu::matrix<double>& input, bnu::matrix<double>& inverse);
int solve_linear_system(const bnu::matrix<double>& A, bnu::vector<double>& x, bnu::vector<double>& y);

#endif // INVERT_MATRIX_H_INCLUDED
