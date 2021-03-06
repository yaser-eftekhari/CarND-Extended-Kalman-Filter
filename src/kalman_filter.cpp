#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init() {
  // State vector
  x_ = VectorXd(4);
  // Transition matrix
  F_ = MatrixXd::Identity(4, 4);
  // Measurement matrix (used in laser case)
  H_ = MatrixXd(2, 4);
  // Process covariance matrix
  Q_ = MatrixXd(4, 4);
  // State covariance matrix P
  P_ = MatrixXd(4, 4);
  // Measurement covariance matrix - radar
  R_R_ = MatrixXd(3, 3);
  // Measurement covariance matrix - laser
  R_L_ = MatrixXd(2, 2);

  H_ << 1, 0, 0, 0,
			  0, 1, 0, 0;

  P_ <<  1, 0, 0, 0,
	       0, 1, 0, 0,
	       0, 0, 1000, 0,
	       0, 0, 0, 1000;

  R_L_ << 0.0225, 0,
          0, 0.0225;

  R_R_ << 0.09, 0, 0,
          0, 0.0009, 0,
          0, 0, 0.09;

  x_size = x_.size();
  I = MatrixXd::Identity(x_size, x_size);
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_L_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];

  float phi = atan2(py,px);
  float rho = sqrt(px*px + py*py);
  float rhoDot;

  if(rho < 0.0001) {
    cout << "Error: rho too small." << endl;
  } else {
    rhoDot = (px*vx + py*vy)/rho;
  }

  VectorXd z_pred = VectorXd(3);
  z_pred << rho, phi, rhoDot;
  MatrixXd Hj = tools_.CalculateJacobian(x_);

  VectorXd y = z - z_pred;
  /* Normalize phi */
  if(y[1] > 0)
    y[1] = fmod(y[1], M_PI);
  else
    y[1] = fmod(y[1], -M_PI);

	MatrixXd Hjt = Hj.transpose();
	MatrixXd S = Hj * P_ * Hjt + R_R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Hjt;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	P_ = (I - K * Hj) * P_;
}
