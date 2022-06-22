#ifndef MEDO_HPP
#define MEDO_HPP

#include "eigen_short.hpp"

#include <vector>

namespace exoukf {

    struct est_t {
        vec<7> xm;
        mat<7,7> Pxx;
    };

    struct meas_t {
        vec<2> z;
        mat<2,2> Pww;
        double t;
        ALNEW
    };

    struct coe_t {
        double sma;
        double ecc;
        double inc;
        double lan;
        double aop;
        double mae;
        double per;
    };

    struct sigpts_t {
        int npts;
        mat<7> X;
        vec<>  w;
        sigpts_t(int npts_in);
    };

    est_t update(const est_t& prior, const meas_t& meas,
            const sigpts_t& sigpts);

    vec<2> get_z(double t, cvec<7> x);

    double get_phi(double t, double n, cvec<2> e);

    vec<7> get_x(const coe_t& coe, double R);

    std::vector<meas_t> gen_meas(cvec<7> xtru, int nmeas, double ti,
            double tf, double stdw, int seed = 0);

    est_t get_prior(double stdx, double stdy, double stdn, double nm);

    mat<7,7> mat_sqrt(cmat<7,7> A);

    sigpts_t stroud_5_3();
    sigpts_t stroud_5_4();
    sigpts_t stroud_7_1();
    sigpts_t stroud_7_2();
    sigpts_t unscented(double k);

}

#endif
