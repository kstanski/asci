#include <complex>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/geometry.hpp>

#include <stdlib.h>

#include "power_spectrum.h"
#include "neighbourhood.h"

Power_spectrum coords2power_spectrum(Position *coords, int coords_no)
{
    /* create a power spectrum (vector of doubles) */
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

    std::complex<double> sh[coords_no][2*L_MAX+1];
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

                    ps(get_ps_idx(l,n1,n2)) += std::conj(c1)*c2;
                    if (n1 != n2)
                    {
                        ps(get_ps_idx(l,n2,n1)) += std::conj(c2)*c1;
                    }
                }

                ps(get_ps_idx(l,n1,n2)) *= sqrt(8/(2*l+1));
                if (n1 != n2)
                {
                    ps(get_ps_idx(l,n2,n1)) *= sqrt(8/(2*l+1));
                }
            }
        }
    }

    /* normalise */
    std::complex<double> norm = norm_2(ps);
    if (norm != std::complex<double>(0,0)) ps /= norm;

    return ps;
}

int cart2sph(Position *coords, int coords_no, double *phi, double *theta, double *r)
{
    namespace bg = boost::geometry;
    bg::model::point<double, DIMENSIONS, bg::cs::spherical<bg::radian> > sph;

    for (int idx=0; idx<coords_no; idx++)
    {
        bg::transform(coords[idx], sph);
        phi[idx] = bg::get<0>(sph);
        theta[idx] = bg::get<1>(sph);
        r[idx] = bg::get<2>(sph);
    }
    return 0;
}

int get_ps_idx(int l, int n1, int n2)
{
    return l*N_MAX*N_MAX + n1*N_MAX + n2;
}

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
