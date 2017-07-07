#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2*n_aug_+1;
  lambda_ = 3 - n_aug_;
  Xsig_pred_ = MatrixXd::Zero(n_aug_, n_sig_);
  weights_.fill(0.5/(lambda_+n_aug_));
  weights_(0) = lambda_/(lambda_+n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  VectorXd z = meas_package.raw_measurements_;
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float px, py, cz = cos(z[1]), sz = sin(z[1]);
      px = cz * z[0];
      py = sz * z[0];
      x_ << px,py,z[2],z[1],0.;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << z[0],z[1],0,0,0;
    }
    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  double dt = (meas_package.timestamp_ - time_us_)/1000000.;
  time_us_ = meas_package.timestamp_;
  Prediction(dt);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    UpdateRadar(meas_package);
  else
    UpdateLidar(meas_package);
}

MatrixXd UKF::generateSigmaPoints(void) {
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sig_);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  //create augmented covariance matrix
  P_aug.block(0,0,n_x_,n_x_) = P_;
  P_aug.block(n_x_,n_x_,2,2) << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_;
  //create square root matrix
  MatrixXd P_sqrt = P_aug.llt().matrixL();
  //create augmented sigma points
  MatrixXd X_expand = MatrixXd(n_aug_,n_aug_);
  X_expand.colwise() = x_aug;
  P_sqrt *= sqrt(lambda_+n_aug_);
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0,1,n_aug_,n_aug_) = X_expand + P_sqrt; 
  Xsig_aug.block(0,1+n_aug_,n_aug_,n_aug_) = X_expand - P_sqrt; 

  return Xsig_aug;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_aug = generateSigmaPoints();

  //predict sigma points
  for (int i=0; i<n_sig_; i++) {
      double px = Xsig_aug.col(i)[0];
      double py = Xsig_aug.col(i)[1];
      double v = Xsig_aug.col(i)[2];
      double yaw = Xsig_aug.col(i)[3];
      double yawd = Xsig_aug.col(i)[4];
      double u_a = Xsig_aug.col(i)[5];
      double u_yawdd = Xsig_aug.col(i)[6];
      
      VectorXd d = VectorXd::Zero(n_x_), u = VectorXd(n_x_);
      if (fabs(yawd) < 0.00001) {
          d[0] = v*cos(yaw)*delta_t;
          d[1] = v*sin(yaw)*delta_t;
      } else {
          d[0] = v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw));
          d[1] = v/yawd*(-cos(yaw+yawd*delta_t)+cos(yaw));
          d[3] = yawd*delta_t;
      }
      u[0] = 0.5*delta_t*delta_t*cos(yaw)*u_a;
      u[1] = 0.5*delta_t*delta_t*sin(yaw)*u_a;
      u[2] = u_a * delta_t;
      u[3] = 0.5*delta_t*delta_t*u_yawdd;
      u[4] = u_yawdd * delta_t;
      
      Xsig_pred_.col(i) = Xsig_aug.col(i).head(5) + d + u;
  }

  VectorXd x_ = Xsig_pred_ * weights_;
  //predict state covariance matrix
  MatrixXd d = Xsig_pred_.colwise() - x_;
  MatrixXd P_ = MatrixXd::Zero(n_x_,n_x_);
  for (int i=0; i<n_sig_; i++) {
      d(3,i) = atan2(sin(d(3,i)),cos(d(3,i)));
      P_ += weights_(i) * d.col(i) * d.col(i).transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;

  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  for (int i=0; i<n_sig_; i++) {
      double px = Xsig_pred_.col(i)[0];
      double py = Xsig_pred_.col(i)[1];
      double v = Xsig_pred_.col(i)[2];
      double yaw = Xsig_pred_.col(i)[3];
      double rho = sqrt(px*px + py*py);
      double phi = atan2(py,px);
      double rhod = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;
      Zsig.col(i) << rho, phi, rhod;
  }
  //calculate mean predicted measurement
  VectorXd z_pred = Zsig * weights_;
  //calculate measurement covariance matrix S
  MatrixXd d = Zsig.colwise() - z_pred;
  MatrixXd S = MatrixXd::Zero(n_z,n_z);
  for (int i=0; i<n_sig_; i++) {
    d(1,i) = atan2(sin(d(1,i)),cos(d(1,i)));
    S += weights_(i)*d.col(i)*d.col(i).transpose();
  }
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  //calculate cross correlation matrix
  for (int i=0; i<n_sig_; i++) {
      VectorXd xdiff = Xsig_pred_.col(i) - x_;
      VectorXd zdiff = Zsig.col(i) - z_pred;
      xdiff(3) = atan2(sin(xdiff(3)), cos(xdiff(3)));
      zdiff(1) = atan2(sin(zdiff(1)), cos(zdiff(1)));
      Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //update state mean and covariance matrix
  VectorXd zdiff = meas_package.raw_measurements_ - z_pred;
  zdiff(1) = atan2(sin(zdiff(1)), cos(zdiff(1)));
  x_ = x_ + K*zdiff;
  P_ -= K*S*K.transpose();
}
