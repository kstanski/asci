
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "stats.h"

Stats produce_stats(bnu::vector<double> vals, bnu::vector<double> prediction)
{
    Stats s;
    int len = vals.size();
    for (int p_idx=0; p_idx<len; p_idx++)
    {
        double err = vals(p_idx) - prediction(p_idx);
        double re = std::abs(err/vals(p_idx));
        s.mre += re;
        if (re > s.re_max) s.re_max = re;
        s.mae += std::abs(err);
        s.rmse += pow(err,2);
    }
    s.mre /= len;
    s.mae /= len;
    s.rmse /= len;
    s.rmse = sqrt(s.rmse);

    std::cout << "RE_MAX: " << s.re_max*100 << "%" << std::endl;
    std::cout << "MRE: " << s.mre*100 << "%" << std::endl;
    std::cout << "MAE: " << s.mae << std::endl;
    std::cout << "RMSE: " << s.rmse << std::endl;
    return s;
}

int output_plot_data(bnu::vector<double> vals, bnu::vector<double> prediction)
{
    std::ofstream file;
    file.open ("data.txt");
    int len = vals.size();
    for (int idx=0; idx<len; idx++)
        file << vals(idx) << " " << prediction(idx) << std::endl;
    file.close();
    return 0;
}
