#ifndef STATS_H_INCLUDED
#define STATS_H_INCLUDED

#include <boost/numeric/ublas/vector.hpp>

namespace bnu = boost::numeric::ublas;

typedef struct s
{
    double re_max;
    double mre;
    double mae;
    double rmse;
} Stats;

Stats produce_stats(bnu::vector<double> vals, bnu::vector<double> prediction);
int output_plot_data(bnu::vector<double> vals, bnu::vector<double> prediction);

#endif // STATS_H_INCLUDED
