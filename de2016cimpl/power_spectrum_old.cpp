#include <complex>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/geometry.hpp>

#include "power_spectrum.h"
#include "neighbourhood.h"

Power_spectrum coords2power_spectrum(double **coords, int coords_no)
{
    /* create a power spectrum (double array) */
    int ps_length = (L_MAX+1)*pow(N_MAX,2);
    Power_spectrum ps = (Power_spectrum) malloc(ps_length*sizeof(double));
    if (ps == NULL)
    {
        fprintf(stderr,"Not enough memory: Power_spectrum");
        exit (1);
    }
    memset(ps,0,ps_length);

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

    /* compute power spectrum terms */
    std::complex<double> sh[coords_no][2*L_MAX+1];
    double temp;
    for (int l=0; l<=L_MAX ; l++)
    {
        //spherical_harmonics(sh,coords_no,l,theta,phi);
        /* compute spherical harmonics */
        for (int idx=0; idx<coords_no; idx++)
        {
            for (int m=-l; m<=l; m++)
            {
                sh[idx][m+l] = boost::math::spherical_harmonic(l,m,theta[idx],phi[idx]);
            }
        }

        /* compute power spectrum */
        for (int n1=0; n1<N_MAX; n1++)
        {
            for (int n2=0; n2<=n1; n2++)
            {
                for (int m=0; m<=2*l; m++)
                {
                    std::complex<double> c1 = 0;
                    std::complex<double> c2 = 0;
                    for (int atom_idx=0; atom_idx<coords_no; atom_idx++)
                    {
                        c1 += rbf[atom_idx][n1]*sh[atom_idx][m];
                        c2 += rbf[atom_idx][n2]*sh[atom_idx][m];
                    }
                    temp = get_element(ps,l,n1,n2);
                    double bla = std::real(std::conj(c1)*c2);
                    set_element(ps,l,n1,n2,temp + std::real(std::conj(c1)*c2));
                    if (n1 != n2)
                    {
                        temp = get_element(ps,l,n2,n1);
                        set_element(ps,l,n2,n1,temp + std::real(std::conj(c2)*c1));
                    }
                }
                temp = get_element(ps,l,n1,n2);
                set_element(ps,l,n1,n2,temp * sqrt(8/(2*l+1)));
                if (n1 != n2)
                {
                    temp = get_element(ps,l,n2,n1);
                    set_element(ps,l,n2,n1,temp * sqrt(8/(2*l+1)));
                }
            }
        }
    }

    /* normalise */
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

double radial_basis_function(double r,double cutoff,int n,int n_max)
{
    return 1.0;
}

int free_power_spectrum(Power_spectrum ps)
{
    free(ps);
    return 0;
}
