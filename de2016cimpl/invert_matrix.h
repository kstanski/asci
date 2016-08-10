#ifndef INVERT_MATRIX_H_INCLUDED
#define INVERT_MATRIX_H_INCLUDED

#include <boost/numeric/ublas/matrix.hpp>
int invert_matrix(const boost::numeric::ublas::matrix<double>& input, boost::numeric::ublas::matrix<double>& inverse);

#endif // INVERT_MATRIX_H_INCLUDED
