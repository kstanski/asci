#include <fstream>
#include <string>
#include <iostream>

#include "read_data.h"

int main()
{
    char filename[] = "test_data.xyz";
    read_data(filename,2);

    return 0 ;
}
