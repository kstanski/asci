#include <complex>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/geometry.hpp>

#include <stdlib.h>

#include "power_spectrum.h"
#include "neighbourhood.h"

Power_spectrum coords2power_spectrum(double **coords, int coords_no)
{
    /* create a power spectrum (double array) */
    int ps_length = (L_MAX+1)*pow(N_MAX,2);
    Power_spectrum ps(ps_length);

    /* convert cartesian to spherical */
    double phi[coords_no], theta[coords_no], r[coords_no];
    cart2sph(coords, coords_no, phi, theta, r);

    /* radial basis functions */
    double rbf[coords_no][N_MAX];
    for (int idx=0; idx<coords_no; idx++)
    {
        for (int n=0; n<N_MAX; n++)
        {
            rbf[idx][n] = radial_basis_function(r[idx],CUTOFF,n,N_MAX);
        }
    }

    /* initialise 2D complex array for spherical harmonics terms */
    /*
    std::complex<double> *sh[coords_no];
    for (int atom_idx=0; atom_idx<coords_no; atom_idx++)
    {
        sh[atom_idx] = (std::complex<double> *) malloc((2*L_MAX+1)*sizeof(std::complex<double>));
        if (sh[atom_idx] == NULL)
        {
            fprintf(stderr,"Not enough memory: spherical harmonics array");
            exit (1);
        }
    }
    */

    std::complex<double> sh[coords_no][2*L_MAX+1];

    int temp_idx;
    for (int l=0; l<=L_MAX ; l++)
    {
        /* fill the spherical harmonics array */
        for (int atom_idx=0; atom_idx<coords_no; atom_idx++)
        {
            for (int m_idx=0; m_idx<2*l+1; m_idx++)
            {
                sh[atom_idx][m_idx] = boost::math::spherical_harmonic(l,m_idx-l,theta[atom_idx],phi[atom_idx]);
            }
        }

        /* compute power spectrum terms */
        for (int n1=0; n1<N_MAX; n1++)
        {
            for (int n2=0; n2<=n1; n2++)
            {
                for (int m_idx=0; m_idx<2*l+1; m_idx++)
                {
                    std::complex<double> c1 = 0;
                    std::complex<double> c2 = 0;
                    for (int atom_idx=0; atom_idx<coords_no; atom_idx++)
                    {
                        c1 += rbf[atom_idx][n1]*sh[atom_idx][m_idx];
                        c2 += rbf[atom_idx][n2]*sh[atom_idx][m_idx];
                    }

                    temp_idx = l*N_MAX*N_MAX + n1*N_MAX + n2;
                    ps(temp_idx) += std::conj(c1)*c2;
                    if (n1 != n2)
                    {
                        temp_idx = l*N_MAX*N_MAX + n2*N_MAX + n1;
                        ps(temp_idx) += std::conj(c2)*c1;
                    }
                }

                temp_idx = l*N_MAX*N_MAX + n1*N_MAX + n2;
                ps(temp_idx) *= sqrt(8/(2*l+1));
                if (n1 != n2)
                {
                    temp_idx = l*N_MAX*N_MAX + n2*N_MAX + n1;
                    ps(temp_idx) *= sqrt(8/(2*l+1));
                }
            }
        }
    }

    /* normalise */
    std::complex<double> norm = norm_2(ps);
    if (norm != std::complex<double>(0,0)) ps /= norm;

    return ps;
}

int cart2sph(double **coords, int coords_no, double *phi, double *theta, double *r)
{
    namespace bg = boost::geometry;
    bg::model::point<double, 3, bg::cs::cartesian> cart;
    bg::model::point<double, 3, bg::cs::spherical<bg::radian> > sph;

    for (int idx=0; idx<coords_no; idx++)
    {
        bg::set<0>(cart, coords[idx][0]);
        bg::set<1>(cart, coords[idx][1]);
        bg::set<2>(cart, coords[idx][2]);
        bg::transform(cart, sph);
        phi[idx] = bg::get<0>(sph);
        theta[idx] = bg::get<1>(sph);
        r[idx] = bg::get<2>(sph);
    }
    return 0;
}
/*
double get_element(Power_spectrum ps,int l, int n1, int n2)
{
    if (l > L_MAX || n1 > N_MAX-1 || n2 > N_MAX-1)
    {
        fprintf(stderr,"Index out of bound: power spectrum (get)");
        exit (1);
    }
    int idx = l*pow(N_MAX,2) + n1*N_MAX + n2;
    return ps[idx];
}

int set_element(Power_spectrum ps,int l, int n1, int n2, double element)
{
    if (l > L_MAX || n1 > N_MAX-1 || n2 > N_MAX-1)
    {
        fprintf(stderr,"Index out of bound: power spectrum (set)");
        exit (1);
    }
    int idx = l*pow(N_MAX,2) + n1*N_MAX + n2;
    ps[idx] = element;
    return 0;
}

int spherical_harmonics(std::complex<double> **sh, int coords_no, int l, double *theta, double *phi)
{
    for (int idx=0; idx<coords_no; idx++)
    {
        for (int m=-l; m<=l; m++)
        {
            sh[idx][m+l] = boost::math::spherical_harmonic(l,m,theta[idx],phi[idx]);
        }
    }
    return 0;
}
*/
double radial_basis_function(double r,double cutoff,int n,int n_max)
{
    double S[n_max][n_max];
    for (int i=0; i<n_max; i++)
        for (int j=0; j<n_max; j++)
            S[i][j] = sqrt((5+2*i)*(5+2*j))/(5+i+j);

    //W = S^-0.5

    double g = 0;
    double n_alpha;
    for (int alpha=0; alpha<n_max; alpha++)
    {
        n_alpha = sqrt(pow(cutoff,(2*alpha+5))/(2*alpha+5));
        g += S[n][alpha] * pow(cutoff-r,alpha+2)/n_alpha;    //change to W
    }

    return g;
}
