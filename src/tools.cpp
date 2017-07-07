#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse = VectorXd(4);
  rmse << 0,0,0,0;
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    std::cout << "CalculateRMSE:invalid arguments" << std::endl;
    return rmse;
  }
  for (int i=0; i<estimations.size(); i++) {
    VectorXd d = ground_truth[i] - estimations[i];
    d = d.array().square();
    rmse = rmse + d;
  }
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;

}
