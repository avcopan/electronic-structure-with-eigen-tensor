#ifndef HELPERS_HPP_
#define HELPERS_HPP_

#include <iostream>
#include <Eigen/Dense>                    // Eigen::Matrix
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor
#include <utility>                        // std::pair
#include "type.hpp"

// overload % for matrix multiplication in the global namespace
Eigen::Tensor<double, 2> operator%(const Eigen::Tensor<double, 2>& L, const Eigen::Tensor<double, 2>& R);

namespace matricks {

// aliases
using IndexPair = Eigen::IndexPair<int>;
template<int naxes>
using Contraction = std::array<IndexPair, naxes>;
template<int naxes>
using Permutation = std::array<int, naxes>;

Eigen::MatrixXd          matrix(Eigen::Tensor<double, 2>);
Eigen::Tensor<double, 2> tensor(Eigen::MatrixXd);
Eigen::Tensor<double, 1> tensor(Eigen::VectorXd);
double trace(Eigen::Tensor<double, 2>);
Eigen::Tensor<double, 2> transpose(Eigen::Tensor<double, 2>);
Eigen::Tensor<double, 2> inverse(Eigen::Tensor<double, 2>);
Eigen::Tensor<double, 2> sqrth(Eigen::Tensor<double, 2>);
std::pair<Eigen::Tensor<double, 1>, Eigen::Tensor<double, 2> > eigh(Eigen::Tensor<double, 2>);

}

#endif // HELPERS_HPP_
