#include <complex>
#include <math.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/geometry.hpp>
#include <boost/simd/pack.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/sum.hpp>

#include <stdlib.h>

#include "power_spectrum.h"
#include "neighbourhood.h"
#include "invert_matrix.h"

Power_spectrum coords2power_spectrum(Position *coords, int coords_no)
{
    /* create a power spectrum (vector of doubles) */
    //int ps_length = (L_MAX+1)*pow(N_MAX,2);
    Power_spectrum ps = (Power_spectrum) malloc(PS_LEN*sizeof(ps_element_type));

    /* convert cartesian to spherical */
    double phi[coords_no], theta[coords_no], r[coords_no];
    cart2sph(coords, coords_no, phi, theta, r);

    /* radial basis functions */
    ps_element_type rbf[coords_no][N_MAX];
    for (int idx=0; idx<coords_no; idx++)
    {
        for (int n=0; n<N_MAX; n++)
        {
            rbf[idx][n] = radial_basis_function(r[idx],CUTOFF,n,N_MAX);
        }
    }

    ps_element_type sh[coords_no][2*L_MAX+1];
    for (int l=0; l<=L_MAX ; l++)
    {
        /* fill the spherical harmonics array */
        for (int atom_idx=0; atom_idx<coords_no; atom_idx++)
        {
            for (int m_idx=0; m_idx<2*l+1; m_idx++)
            {
                sh[atom_idx][m_idx] = sh_real_form(l,m_idx-l,theta[atom_idx],phi[atom_idx]);
            }
        }

        /* compute power spectrum terms */
        for (int n1=0; n1<N_MAX; n1++)
        {
            for (int n2=0; n2<=n1; n2++)
            {
                ps_element_type ps_elem1 = 0;
                ps_element_type ps_elem2 = 0;
                for (int m_idx=0; m_idx<2*l+1; m_idx++)
                {
                    ps_element_type c1 = 0;
                    ps_element_type c2 = 0;
                    for (int atom_idx=0; atom_idx<coords_no; atom_idx++)
                    {
                        c1 += rbf[atom_idx][n1]*sh[atom_idx][m_idx];
                        c2 += rbf[atom_idx][n2]*sh[atom_idx][m_idx];
                    }

                    ps_elem1 += c1*c2;
                    if (n1 != n2)
                    {
                        ps_elem2 += c2*c1;
                    }
                }
                ps[get_ps_idx(l,n1,n2)] = ps_elem1 * sqrt(8/(2*l+1));
                if (n1 != n2)
                {
                    ps[get_ps_idx(l,n2,n1)] = ps_elem2 * sqrt(8/(2*l+1));
                }
            }
        }
    }

    normalise(ps);

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

ps_element_type sh_real_form(int l, int m, double theta, double phi)
{
    std::complex<ps_element_type> sh = (std::complex<ps_element_type>) boost::math::spherical_harmonic(l,std::abs(m),theta,phi);
    if (m == 0)
        return std::real(sh);
    else if (m < 0)
        return sqrt(2) * pow(-1,m) * std::imag(sh);
    else
        return sqrt(2) * pow(-1,m) * std::real(sh);
}

int get_ps_idx(int l, int n1, int n2)
{
    return l*N_MAX*N_MAX + n1*N_MAX + n2;
}

ps_element_type radial_basis_function(double r,double cutoff,int n,int n_max)
{
    namespace bnu = boost::numeric::ublas;
    bnu::matrix<double> S(n_max,n_max), W(n_max,n_max);
    for (int i=0; i<n_max; i++)
    {
        S(i,i) = 1;
        for (int j=0; j<i; j++)
        {
             double val = sqrt((5+2*i)*(5+2*j))/(5+i+j);
             S(i,j) = val;
             S(j,i) = val;
        }
    }

    //can also try W = S^-0.5
    //invert_matrix(S,W);
    W = S;

    ps_element_type g = 0;
    double n_alpha;
    for (int alpha=0; alpha<n_max; alpha++)
    {
        n_alpha = sqrt(pow(cutoff,(2*alpha+5))/(2*alpha+5));
        g += W(n,alpha) * pow(cutoff-r,alpha+2)/n_alpha;
    }

    return g;
}

double dot_prod(Power_spectrum A, Power_spectrum B)
{
    namespace bs = boost::simd;
    using pack = bs::pack<ps_element_type>;
    size_t pack_card = bs::cardinal_of<pack>();

    pack p_dot{0}, p_A, p_B;
    int iter_no = PS_LEN/pack_card;

    for (int i=0; i<iter_no; i++)
    {
        int idx = i*pack_card;
        p_A = bs::load<pack>(A + idx);
        p_B = bs::load<pack>(B + idx);
        p_dot += p_A * p_B;
    }

    double dot = bs::sum(p_dot);
    for (int i=iter_no*pack_card; i<PS_LEN; i++)
        dot += A[i] * B[i];

    return dot;
}

int normalise(Power_spectrum ps)
{
    ps_element_type norm = sqrt(dot_prod(ps,ps));
    if (norm != 0)
    {
        for (int i=0; i<PS_LEN; i++) ps[i] /= norm;
    }
    return 0;
}
