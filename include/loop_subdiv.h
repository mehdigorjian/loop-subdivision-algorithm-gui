#ifndef LOOP_SUBDIVISION_H
#define LOOP_SUBDIVISION_H
#include <Eigen/Core>

void loop_subdivision(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const int num_iters,
    Eigen::MatrixXd &SV,
    Eigen::MatrixXi &SF);
#endif
