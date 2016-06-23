
%The lengths of set_1 and set_2 have to add up to the length of sorted_set.

function [set_1,set_2] = stratify(sorted_set,set_1,set_2)

    l_sorted = length(sorted_set);
    l1 = length(set_1);
    set_1_last = 0;
    set_2_last = 0;

    low_step = floor(l_sorted/l1);
    high_step = low_step+1;
    missing = l1-floor(l_sorted/high_step);

    for i = 1:l_sorted
        if i < low_step*high_step*missing    %make up for missing ones
            step = low_step;
        else
            step = high_step;
        end
        if rem(i,step) == 0
            set_1_last = set_1_last+1;
            set_1(set_1_last) = sorted_set(i);
        else
            set_2_last = set_2_last+1;
            set_2(set_2_last) = sorted_set(i);
        end
    end

end