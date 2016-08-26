#include <math.h>
#include "stratify.h"

int stratify(Molecule **sorted_arr, Molecule **train_arr, int train_no, Molecule **validate_arr, int validate_no)
{
    int len = train_no + validate_no;
    int small_last = 0;
    int large_last = 0;
    int low_step, high_step, missing, step, large_no, small_no;
    Molecule **small_arr, **large_arr;

    if (train_no < validate_no)
    {
        small_arr = train_arr;
        small_no = train_no;
        large_arr = validate_arr;
        large_no = validate_no;
        low_step = floor(len/train_no);
        high_step = low_step + 1;
        missing = train_no - floor(len/high_step);
    } else {
        small_arr = validate_arr;
        small_no = validate_no;
        large_arr = train_arr;
        large_no = train_no;
        low_step = floor(len/validate_no);
        high_step = low_step + 1;
        missing = validate_no - floor(len/high_step);
    }

    for (int i=0; i<len; i++)
    {
        if (i < low_step*high_step*missing) step = low_step;
        else step = high_step;

        if (small_last < small_no && (i%step == 0 || large_last+1 == large_no))
        {
            small_arr[small_last] = sorted_arr[i];
            small_last++;
        } else {
            large_arr[large_last] = sorted_arr[i];
            large_last++;
        }
    }
    return 0;
}
