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

  is_initialized_ = false;

  //previous_timestamp_ = 0;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
	
  n_x_ = 5;

  
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  
  Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;//30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;//30;

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
  



  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);


  
  time_us_ = 0.0;


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

  if (!is_initialized_) {
    cout << "UKF: " << endl;
    x_ << 1, 1, 1, 1, 1;

    P_ << .1, 0, 0, 0, 0,
          0, .1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0 , 1;
    
	time_us_ = meas_package.timestamp_;

	
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double ro = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      //double ro_dot = meas_package.raw_measurements_(2);

      x_ << ro*cos(phi), ro*sin(phi) , 1, 1, 0.1;
    
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 1, 1, 0.1;
    }


    
    is_initialized_ = true;
    return;
}

	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
	time_us_ = meas_package.timestamp_;
	
	Prediction(delta_t);
	
	
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLidar (meas_package);
     }
	
	cout << "x_ = " << x_ << endl;
	cout << "P_ = " << P_ << endl;
	

}


  
  


void UKF::Prediction(double delta_t) {
	
	
	VectorXd x_aug = VectorXd(n_aug_);
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;
	
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P_;
	P_aug(5,5) = std_a_*std_a_;
	P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
	Xsig_aug.col(0)  = x_aug;
	for (int i = 0; i< n_aug_; i++){
		Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
	}


	
	for (int i=0; i<2*n_aug_ + 1; i++){
       
       double p_x = Xsig_aug(0,i);
       double p_y = Xsig_aug(1,i);
       double v = Xsig_aug(2,i);
       double yaw = Xsig_aug(3,i);
       double yaw_d = Xsig_aug(4,i);
       double nu_a = Xsig_aug(5,i);
       double nu_yawdd = Xsig_aug(6,i);
       
       double px_p, py_p, pv, pyaw, pyaw_d;
       
       if (fabs(yaw_d) > 0.001){
           px_p = p_x + (v/yaw_d)*(sin(yaw+yaw_d*delta_t)-sin(yaw));
           py_p = p_y + (v/yaw_d)*(cos(yaw)-cos(yaw+yaw_d*delta_t));
       }
       
       else { 
           px_p = p_x + v*cos(yaw)*delta_t;
           py_p = p_y + v*sin(yaw)*delta_t;
       }
        
        pv = v;
        pyaw = yaw + yaw_d*delta_t;
        pyaw_d = yaw_d;
       
       
       px_p = px_p + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
       py_p = py_p + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
       pv = pv + delta_t*nu_a;
       pyaw = pyaw + 0.5*delta_t*delta_t*nu_yawdd;
       pyaw_d = pyaw_d + delta_t * nu_yawdd;
       
       
	
       Xsig_pred(0,i) = px_p;
       Xsig_pred(1,i) = py_p;
       Xsig_pred(2,i) = pv;
       Xsig_pred(3,i) = pyaw;
       Xsig_pred(4,i) = pyaw_d;  

   }
	
	
	double weight_0 = lambda_/(lambda_+n_aug_);
	weights_(0) = weight_0;
	for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
		double weight = 0.5/(n_aug_+lambda_);
		weights_(i) = weight;
	}  


     
	x_.fill(0.0);
	for (int i=0; i < 2*n_aug_+1; i++){
      x_ = x_ + weights_(i) * Xsig_pred.col(i);
	}
  
	P_.fill(0.0);
  
	for (int i=0; i < 2*n_aug_+1; i++){
      VectorXd x_diff = Xsig_pred.col(i) - x_;
      while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
      while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
	}
    
  }

  



void UKF::UpdateLidar(MeasurementPackage meas_package) {

	
	
	
	
	//VectorXd z = VectorXd(n_x_);
	VectorXd z = meas_package.raw_measurements_;
	
	MatrixXd H_laser_ = MatrixXd(2, 5);
	H_laser_<< 1, 0, 0, 0, 0,
				0, 1, 0, 0, 0;
	
	MatrixXd R_laser_ = MatrixXd(2, 2);
	R_laser_ << std_laspx_*std_laspx_, 0,
				0, std_laspy_*std_laspy_;
	
  	VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;
}



 
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	

  int n_z = 3;
	VectorXd z = VectorXd(n_z);
	z = meas_package.raw_measurements_;
	

	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
	VectorXd z_pred = VectorXd(n_z);
	MatrixXd S = MatrixXd(n_z,n_z);
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	
	for (int i=0;i<2*n_aug_ + 1;i++){
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        double v = Xsig_pred(2,i);
        double yaw = Xsig_pred(3,i);
        double yaw_d = Xsig_pred(4,i);
        
        
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);
        Zsig(1,i) = atan2(p_y,p_x);
        Zsig(2,i) = (p_x*v*cos(yaw) + p_y*v*sin(yaw))/sqrt(p_x*p_x + p_y*p_y);
     }
     
     z_pred.fill(0.0);
     
     for (int i=0;i<2*n_aug_ +1;i++){
         z_pred = z_pred + weights_(i)*Zsig.col(i);
     }
     
     S.fill(0.0);
     
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
    
  }
  
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_*std_radrd_;
    S = S + R;
	
	
	Tc.fill(0.0);
    
    for (int i=0;i<2*n_aug_+1; i++){
        VectorXd x_diff = Xsig_pred.col(i)-x_;
        
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        
        VectorXd z_diff = Zsig.col(i)-z_pred;
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        
        Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
    }
    
    
    MatrixXd Kc = Tc * S.inverse();
    
    VectorXd z_diff2 = z-z_pred;
    while (z_diff2(1)> M_PI) z_diff2(1)-=2.*M_PI;
    while (z_diff2(1)<-M_PI) z_diff2(1)+=2.*M_PI;
    
    x_ = x_ + Kc*z_diff2;
    
    P_ = P_ - Kc*S*Kc.transpose();
	



}

