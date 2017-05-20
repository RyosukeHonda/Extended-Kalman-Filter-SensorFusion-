#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout <<"Kalman Filter Initialization" << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;

    //the initial transition matrix F_
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;

    //state covariance matrix P
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;

    H_laser_ << 1 , 0 , 0 , 0,
               0 , 1 , 0 , 0;

     //measurement matrix
    Hj_ << 1 , 1 , 0 , 0,
           1 , 1 , 0 , 0,
           1 , 1 , 1 , 1;

     //set the acceleration noise components
     noise_ax = 9;
     noise_ay = 9;

     R_laser_ <<  0.0225, 0,
                  0, 0.0225;

      //measurement covariance
     R_radar_ << 0.09 , 0 , 0 ,
                 0 , 0.0009 , 0 ,
                 0 , 0 , 0.09 ;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //x = rho * cos(phi), y = rho *sin(phi)
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];
      float x_cart =  rho * cos(phi);
      float y_cart =  rho * sin(phi);
      float vx_cart = rho_dot * cos(phi);
      float vy_cart = rho_dot * sin(phi);

      //avoid 0 division when initializing
      if (x_cart == 0 or y_cart ==0){

        return;
      }

      ekf_.x_ << x_cart, y_cart,vx_cart,vy_cart;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      //avoid 0 division when initializing
      if (measurement_pack.raw_measurements_[0] == 0 or measurement_pack.raw_measurements_[1] == 0){
        return;
      }

      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1] , 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds

  previous_timestamp_ = measurement_pack.timestamp_;
  //cout << dt << endl;

  //Modify the F matrix so that the time is integrated
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    double range = sqrt(pow(ekf_.x_[0],2) + pow(ekf_.x_[1],2));
    double bearing = atan2(ekf_.x_[1],ekf_.x_[0]);
    double range_rate = ((ekf_.x_[0]*ekf_.x_[2]+ekf_.x_[1]*ekf_.x_[3])/(sqrt(pow(ekf_.x_[0],2) + pow(ekf_.x_[1],2))));
    MatrixXd zpred(3, 1);
    zpred << range, bearing, range_rate;

    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, zpred);
    //ekf_.Update(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
