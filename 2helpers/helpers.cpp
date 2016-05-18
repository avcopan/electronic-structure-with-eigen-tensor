#include "helpers.hpp"

Eigen::Tensor<double, 2> operator%(const Eigen::Tensor<double, 2>& L, const Eigen::Tensor<double, 2>& R) {
  std::array<Eigen::IndexPair<int>, 1> matmul({Eigen::IndexPair<int>(1, 0)});
  return L.contract(R, matmul);
}

namespace matricks {

Eigen::MatrixXd matrix(Eigen::Tensor<double, 2> m) {
  Eigen::Map<Eigen::MatrixXd> M(m.data(), m.dimension(0), m.dimension(1));
  return M;
}

Eigen::Tensor<double, 2> tensor(Eigen::MatrixXd M) {
  Eigen::TensorMap<Eigen::Tensor<double, 2> > m(M.data(), M.rows(), M.cols());
  return m;
}

Eigen::Tensor<double, 1> tensor(Eigen::VectorXd V) {
  Eigen::TensorMap<Eigen::Tensor<double, 1> > v(V.data(), V.size());
  return v;
}

double trace(Eigen::Tensor<double, 2> m) {
  return matrix(m).trace();
}

Eigen::Tensor<double, 2> transpose(Eigen::Tensor<double, 2> m) {
  Permutation<2> transpose{{1, 0}};
  return m.shuffle(transpose);
}

Eigen::Tensor<double, 2> inverse(Eigen::Tensor<double, 2> m) {
  Eigen::MatrixXd minv = matrix(m).inverse();
  return tensor(minv);
}

Eigen::Tensor<double, 2> sqrth(Eigen::Tensor<double, 2> m) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(matrix(m));
  Eigen::VectorXd x = eigensolver.eigenvalues().array().pow(1./2);
  Eigen::MatrixXd u = eigensolver.eigenvectors();
  Eigen::MatrixXd msqrt = u * x.asDiagonal() * u.transpose();
  return tensor(msqrt);
}

std::pair<Eigen::Tensor<double, 1>, Eigen::Tensor<double, 2> > eigh(Eigen::Tensor<double, 2> m) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(matrix(m));
  Eigen::Tensor<double, 1> eigenvalues  = tensor( eigensolver.eigenvalues()  );
  Eigen::Tensor<double, 2> eigenvectors = tensor( eigensolver.eigenvectors() );
  return std::make_pair(eigenvalues, eigenvectors);  
}

}
