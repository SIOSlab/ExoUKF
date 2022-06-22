#include "exoukf.hpp"

#include "eigen_csv.hpp"

#include <cmath>
#include <iostream>

int main() {

    using namespace exoukf;

    coe_t coe;
    coe.sma = 3;
    coe.ecc = 0.1;
    coe.inc = 30;
    coe.lan = 120;
    coe.aop = 45;
    coe.mae = 75;
    coe.per = 5;

    double R = 2.5; 

    vec<7> xtru = get_x(coe, R);

    int nmeas = 5;

    double ti = 0;
    double tf = 1;
    
    double stdw = 0.005;
   
    vec<9> par;
    par << coe.sma, coe.ecc, coe.inc, coe.lan, coe.aop, coe.mae, coe.per, R, stdw;
    eigen_csv::write(par, "par.csv");
    
    eigen_csv::write(xtru, "xtru.csv");

    double stdx = 1;
    double stdy = 1;
    double stdn = 1;
    double nm = 0;

    est_t prior = get_prior(stdx, stdy, stdn, nm);

    sigpts_t sigpts = stroud_7_2();

    int npass = 20;
    
    int ntrials = 100;

    mat<> xi_sqerr(npass+1, ntrials);
    mat<> eta_sqerr(npass+1, ntrials);
    mat<> lambda_sqerr(npass+1, ntrials);

    std::vector<est_t> ests(npass+1);
        
    est_t est = prior;

    for (int m = 0; m < ntrials; m++) {

        std::vector<meas_t> meas = gen_meas(xtru, nmeas, ti, tf, stdw, m);

        ests[0] = prior;

        for (int n = 1; n <= npass; n++) {

            for (const meas_t& m : meas)
                est = update(est, m, sigpts); 

            ests[n] = est;

        }

        for (int n = 0; n <= npass; n++) {

            vec<7> dx = ests[n].xm - xtru;

            xi_sqerr(n, m) = dx.head<4>().squaredNorm();

            eta_sqerr(n, m) = dx.segment<2>(4).squaredNorm();

            lambda_sqerr(n, m) = dx(6) * dx(6);

        }

    }

    vec<> xi_rmse = xi_sqerr.rowwise().mean().cwiseSqrt(); 
    vec<> eta_rmse = eta_sqerr.rowwise().mean().cwiseSqrt(); 
    vec<> lambda_rmse = lambda_sqerr.rowwise().mean().cwiseSqrt(); 

    vec<> xi_rmse_rel = xi_rmse / xtru.head<4>().norm();
    vec<> eta_rmse_rel = eta_rmse / xtru.segment<2>(4).norm();
    vec<> lambda_rmse_rel = lambda_rmse / abs(xtru(6));

    vec<> pass = vec<>::LinSpaced(npass+1, 0, npass);

    mat<> table(npass+1, 7);
    table << pass, xi_rmse, xi_rmse_rel, eta_rmse, eta_rmse_rel, 
          lambda_rmse, lambda_rmse_rel;

    eigen_csv::write(table, "results.csv");

}
