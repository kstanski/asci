#ifndef SETUP_H_INCLUDED
#define SETUP_H_INCLUDED

#include <boost/numeric/ublas/vector.hpp>

#include "setup.h"
#include "stratify.h"



#include <fstream>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <time.h>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "molecule.h"
#include "neighbourhood.h"
#include "power_spectrum.h"
#include "solver.h"
#include "stratify.h"
#include "descriptor.h"
#include "local_similarity.h"
#include "structural_similarity.h"
#include "stats.h"

namespace bnu = boost::numeric::ublas;

typedef struct dataset
{
    int train_no;
    Descriptor *train_desc;
    double *train_val;
    int validate_no;
    Descriptor *validate_desc;
    double *validate_val;
} Dataset;

Dataset *setup(const char *filename, int molecules_no, int train_no, int validate_no);
Dataset *create_dataset(int train_no, int validate_no);

#endif // SETUP_H_INCLUDED
