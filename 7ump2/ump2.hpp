#ifndef UMP2_HPP_
#define UMP2_HPP_

#include <iostream> // std::cout
#include <cstdio>   // std::printf
#include <liboptions/liboptions.h>        // psi::Options
#include <unsupported/Eigen/CXX11/Tensor> // Eigen::Tensor
#include "../6uhf/uhf.hpp"

namespace ump2 {

using IndexPair = Eigen::IndexPair<int>;
template<class C> using Shared      = boost::shared_ptr<C>;
template<int n>   using Contraction = std::array<IndexPair, n>;
template<int n>   using Tuple       = std::array<int      , n>;

class UMP2 {
private:
  int nbf_;
  int naocc_;
  int nbocc_;
  int navir_;
  int nbvir_;
  double hfenergy_;
  double cenergy_;
  psi::Options options_;
  Shared<uhf::UHF> uhf_;
  Eigen::Tensor<double, 1> ae_;
  Eigen::Tensor<double, 1> be_;
  Eigen::Tensor<double, 4> aag_;
  Eigen::Tensor<double, 4> abg_;
  Eigen::Tensor<double, 4> bbg_;
public:
  UMP2(Shared<uhf::UHF>, psi::Options&);
  double compute_energy();
};

} // ump2

#endif // UMP2_HPP_
