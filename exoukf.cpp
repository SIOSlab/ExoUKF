#include "exoukf.hpp"

#include "angles.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

using namespace exoukf;

sigpts_t::sigpts_t(int npts_in) :
    npts(npts_in),
    X(7, npts),
    w(npts)
{}

est_t exoukf::update(const est_t& prior, const meas_t& meas,
        const sigpts_t& sigpts) {

    mat<7,7> L = prior.Pxx.llt().matrixL();

    mat<7> Xc = L * sigpts.X;

    mat<7> X = Xc.colwise() + prior.xm;

    mat<2> Z(2, sigpts.npts);
    for (int j = 0; j < sigpts.npts; j++)
        Z.col(j) = get_z(meas.t, X.col(j));

    vec<2> zm = Z * sigpts.w;

    mat<2> Zc = Z.colwise() - zm;

    mat<2,2> Pzz = Zc * sigpts.w.asDiagonal() * Zc.transpose() + meas.Pww;
    mat<7,2> Pxz = Xc * sigpts.w.asDiagonal() * Zc.transpose();

    mat<7,2> K = Pxz * Pzz.inverse();

    est_t post = prior;

    post.xm += K * (meas.z - zm);

    post.Pxx -= K * Pxz.transpose();

    return post;

}

vec<2> exoukf::get_z(double t, cvec<7> x) {

    mat<2,2> T;
    T(0, 0) = x(0);
    T(1, 0) = x(1);
    T(0, 1) = x(2);
    T(1, 1) = x(3);

    vec<2> eta = x.segment<2>(4);

    vec<2> e = eta / sqrt(1 + eta.squaredNorm());

    double f = sqrt(1 - e.squaredNorm());

    //double n = x(6);
    double n = 2 * pi / exp(x(6));

    double phi = get_phi(t, n, e);

    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    double s = (e.x() * sin_phi - e.y() * cos_phi) / (1 + f);

    vec<2> ro;
    ro.x() = cos_phi - e.x() + e.y() * s;
    ro.y() = sin_phi - e.y() - e.x() * s;

    return T * ro;

}

double exoukf::get_phi(double t, double n, cvec<2> e) {

    constexpr double tol = 1E-6;

    double dM = fmod(n * t, 2*pi);

    double phi = dM;

    double dphi;

    int k = 0;

    do {

        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        
        dphi = (phi - dM - e.x()*sin_phi + e.y()*cos_phi) / 
            (1 - e.x() * cos_phi - e.y() * sin_phi);

        phi -= dphi;

        k++;

        if (k > 100) {
            std::cout << "ITERATION FAILED!" << std::endl;
            std::cout << "n = " << n << std::endl;
            std::cout << "t = " << t << std::endl;
            std::cout << "e = " << e.x() << ", " << e.y() << std::endl;
            exit(EXIT_FAILURE);
        }

    } while (fabs(dphi) > tol);

    return phi;

}

vec<7> exoukf::get_x(const coe_t& coe, double R) {

    double a = coe.sma;
    double e = coe.ecc;

    double i = deg2rad(coe.inc);
    double W = deg2rad(coe.lan);
    double w = deg2rad(coe.aop + coe.mae);

    double M0 = deg2rad(coe.mae);

//    double n = 2 * pi / coe.per;
   
    double s = a / R;

    double A = s * ( cos(w)*cos(W) - sin(w)*sin(W)*cos(i));
    double B = s * ( cos(w)*sin(W) + sin(w)*cos(W)*cos(i));
    double F = s * (-sin(w)*cos(W) - cos(w)*sin(W)*cos(i));
    double G = s * (-sin(w)*sin(W) + cos(w)*cos(W)*cos(i));

    double etax =  e * cos(M0) / sqrt(1 - e*e);
    double etay = -e * sin(M0) / sqrt(1 - e*e);
        
    vec<7> x;

    x(0) = A;
    x(1) = B;
    x(2) = F;
    x(3) = G;
    x(4) = etax;
    x(5) = etay;
//    x(6) = n;
    x(6) = log(coe.per);

    return x;

}

std::vector<meas_t> exoukf::gen_meas(cvec<7> xtru, int nmeas, double ti,
            double tf, double stdw, int seed) {

    std::mt19937_64 mt(seed);

    std::normal_distribution<double> nd;

    vec<> t = vec<>::LinSpaced(nmeas, ti, tf); 

    mat<2,2> Pww = mat<2,2>::Identity() * stdw * stdw;

    std::vector<meas_t> meas(nmeas);

    for (int k = 0; k < nmeas; k++) {
       
        vec<2> w;
        w.x() = stdw * nd(mt);
        w.y() = stdw * nd(mt);

        meas[k].z = get_z(t(k), xtru) + w;
        meas[k].Pww = Pww;
        meas[k].t = t(k);

    }

    return meas;

}

est_t exoukf::get_prior(double stdx, double stdy, double stdn, double nm) {

    est_t prior;
    
    prior.xm.setZero();
    prior.xm(6) = nm;

    prior.Pxx.setIdentity();
    prior.Pxx.diagonal().head<4>() *= stdx * stdx;
    prior.Pxx.diagonal().segment<2>(4) *= stdy * stdy;
    prior.Pxx(6,6) = stdn * stdn;

    return prior;

}

mat<7,7> exoukf::mat_sqrt(cmat<7,7> A) {
    using namespace Eigen;
    //BDCSVD<mat<>> svd(A, ComputeFullU);
    JacobiSVD<mat<7,7>> svd(A, ComputeFullU);
    return svd.matrixU() * svd.singularValues().cwiseSqrt().asDiagonal();    
}

using namespace std;

#include "stroud.hpp"

sigpts_t exoukf::stroud_5_3() {

    int npts = en_r2_05_3_size(7);

    sigpts_t sigpts(npts);

    en_r2_05_3(7, npts, sigpts.X.data(), sigpts.w.data());

    double wsum = sigpts.w.sum();

    sigpts.w /= wsum;

    sigpts.X *= sqrt(2);

    return sigpts;

}

sigpts_t exoukf::stroud_5_4() {

    int npts = en_r2_05_4_size(7);

    sigpts_t sigpts(npts);

    en_r2_05_4(7, npts, sigpts.X.data(), sigpts.w.data());

    double wsum = sigpts.w.sum();

    sigpts.w /= wsum;

    sigpts.X *= sqrt(2);

    return sigpts;

}

sigpts_t exoukf::stroud_7_1() {

    int npts = en_r2_07_1_size(7);

    sigpts_t sigpts(npts);

    en_r2_07_1(7, 1, npts, sigpts.X.data(), sigpts.w.data());

    double wsum = sigpts.w.sum();

    sigpts.w /= wsum;

    sigpts.X *= sqrt(2);

    return sigpts;

}

sigpts_t exoukf::stroud_7_2() {

    int npts = en_r2_07_2_size(7);

    sigpts_t sigpts(npts);

    en_r2_07_2(7, npts, sigpts.X.data(), sigpts.w.data());

    double wsum = sigpts.w.sum();

    sigpts.w /= wsum;

    sigpts.X *= sqrt(2);

    return sigpts;

}

sigpts_t exoukf::unscented(double k) {

    sigpts_t sigpts(15);

    double s = sqrt(7 + k);

    double ws = 0.5 / (7 + k);
    double wc =   k / (7 + k);

    sigpts.w.setConstant(ws);
    sigpts.w(7) = wc;

    sigpts.X.col(7).setZero();
    sigpts.X.leftCols<7>()  =  s * mat<7,7>::Identity();
    sigpts.X.rightCols<7>() = -s * mat<7,7>::Identity();

    return sigpts;

}
