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
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

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
  
  lambda_ = 3 - n_aug_;
    
  weights_ = VectorXd(2*n_aug_+1);
  
  //define spreading parameter
  weights_(0) = lambda_ / (lambda_+n_aug_);
  
  for(int i=0; i<2*n_aug_; i++) 
	  weights_(i+1)=1/(2*(lambda_+n_aug_));
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
  
  if(is_initialized_){
  
    //how much time has elapsed
//    cout << "New T:" << meas_package.timestamp_ << endl;
//    cout << "Delta T long long:" << meas_package.timestamp_ - time_us_ << endl;
    double delta_t = double(meas_package.timestamp_ - time_us_) / 1000000;
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
//      cout << "LIDAR" << endl;	
      //predict
      Prediction(delta_t);    
      //constrat with measurement and update
      UpdateLidar(meas_package);   
      //update the last timestamp
      time_us_ = meas_package.timestamp_;
    }
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
//        cout << "RADAR" << endl;
        //predict
        Prediction(delta_t);    
        //constrat with measurement and update
        UpdateRadar(meas_package);
        //update the last timestamp
        time_us_ = meas_package.timestamp_;
    }
    
  }else{
    //initialise
    time_us_ = meas_package.timestamp_;
    //ekf_.x_ << 1, 1, 1, 1;
    x_ = VectorXd(n_x_);
       
    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
	// cout << "Initialize LIDAR" << endl;
    }else{
      /**
      Initialize state.
      */
      double rho = meas_package.raw_measurements_[0]; 
      double theta = meas_package.raw_measurements_[1];
      double rhod = meas_package.raw_measurements_[2];
      x_ << rho * cos(theta), rho * sin(theta), rhod, 0, 0;
	//cout << "Initialize RADAR" << endl;
    }
    P_ = MatrixXd(5,5);
    P_.fill(0.0);
    P_(0,0)=1;
    P_(1,1)=1;
    P_(2,2)=10;
    P_(3,3)=10;
    P_(4,4)=10;

	//cout << "Initial x: " << x_ << endl;
	//cout << "Initial P: " << P_ << endl;
	//cout << "Initial t: " << time_us_ << endl;
    
    is_initialized_ = true;
  } 
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
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //  cout << "Dt" << endl << delta_t << endl;
  //  cout << "x_aug" << endl << x_aug << endl;
  
  P_aug.fill(0.0);
  //create augmented covariance matrix
  P_aug.topLeftCorner( n_x_, n_x_ ) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd root =  sqrt(lambda_ + n_aug_) * A;
    
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++){
      Xsig_aug.col(1 + i) = x_aug + root.col(i);
      Xsig_aug.col(1 + n_aug_ + i) = x_aug - root.col(i);
//    Xsig_aug(3,1 + i) = Normalize(Xsig_aug(3,1 + i));
//    Xsig_aug(3,1 + n_aug_ + i) = Normalize(Xsig_aug(3,1 + i));
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
}
