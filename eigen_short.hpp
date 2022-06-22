#ifndef EIGEN_SHORT_HPP
#define EIGEN_SHORT_HPP

#include <Eigen/Dense>

template <class T>
    using map = Eigen::Map<T>;

template <int M = Eigen::Dynamic, int N = Eigen::Dynamic>
    using mat = Eigen::Matrix<double, M, N>;

template <int M = Eigen::Dynamic, int N = Eigen::Dynamic>
    using rmat = Eigen::Ref<mat<M,N>>;

template <int M = Eigen::Dynamic, int N = Eigen::Dynamic>
    using cmat = const Eigen::Ref<const mat<M,N>>;

template <int M = Eigen::Dynamic>
    using vec = mat<M,1>;

template <int M = Eigen::Dynamic>
    using rvec = rmat<M,1>;

template <int M = Eigen::Dynamic>
    using cvec = cmat<M,1>;

using quat = Eigen::Quaterniond;

using angax = Eigen::AngleAxisd;

#define ALNEW EIGEN_MAKE_ALIGNED_OPERATOR_NEW

#endif
