#include <RcppArmadillo.h>
#include <stdexcept>
#include <math.h>
#include <stdlib.h> /* srand, rand */
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]

List log_IC_2var_Case3(arma::vec sharedtype, List HB, arma::vec Gompertz, arma::vec Weibull,
                       arma::vec nb_points_integral, arma::vec alpha_inter_intra,
                       arma::vec alpha_y_slope, List alpha_z, List gamma_B, arma::vec beta, arma::vec beta_slope,
                       arma::mat b_y, arma::mat b_y_slope, arma::vec wk, arma::vec rep_wk, arma::vec sigma_inter, arma::vec sigma_intra,
                       int delta2_i, arma::rowvec Z_01_i, arma::rowvec Z_02_i, arma::rowvec Z_12_i, arma::rowvec X_T_i, arma::vec U_T_i,
                       arma::rowvec Xslope_T_i, arma::vec Uslope_T_i, arma::mat X_GK_T_i, arma::mat U_GK_T_i, arma::mat Xslope_GK_T_i,
                       arma::mat Uslope_GK_T_i, arma::mat X_GK_L_T_i, arma::mat U_GK_L_T_i, arma::mat Xslope_GK_L_T_i, arma::mat Uslope_GK_L_T_i,
                       arma::mat X_GK_0_LT_i, arma::mat U_GK_0_LT_i, arma::mat Xslope_GK_0_LT_i, arma::mat Uslope_GK_0_LT_i,
                       arma::mat X_GK_T0_i, arma::mat U_GK_T0_i, arma::mat Xslope_GK_T0_i, arma::mat Uslope_GK_T0_i,
                       double Time_T_i, double Time_L_T_i, double Time_T0_i,arma::vec st_T_i, arma::mat st_0_LT_i, arma::vec st_L_T_i, arma::vec st_T0_i,
                       arma::vec ck,
                       arma::vec B_T_i_01, arma::vec B_T_i_02, arma::vec B_T_i_12,
                       arma::mat Bs_T_i_01, arma::mat Bs_T_i_02, arma::mat Bs_T_i_12,
                       arma::mat Bs_0_LT_i_01, arma::mat Bs_0_LT_i_02, arma::mat Bs_0_LT_i_12,
                       arma::mat Bs_L_T_i_01,  arma::mat Bs_L_T_i_02, arma::mat Bs_L_T_i_12,
                       arma::mat Bs_T0_i_01, arma::mat Bs_T0_i_02, bool left_trunc
){

  // parameters
  bool dep_cv_01 = sharedtype[0];
  bool dep_slope_01 = sharedtype[1];
  bool dep_var_inter_01 = sharedtype[2];
  bool dep_var_intra_01 = sharedtype[3];
  bool dep_cv_02 = sharedtype[4];
  bool dep_slope_02 = sharedtype[5];
  bool dep_var_inter_02 = sharedtype[6];
  bool dep_var_intra_02 = sharedtype[7];
  bool dep_cv_12 = sharedtype[8];
  bool dep_slope_12 = sharedtype[9];
  bool dep_var_inter_12= sharedtype[10];
  bool dep_var_intra_12 = sharedtype[11];
  const std::string& hazard_baseline_01 = HB[0];
  const std::string& hazard_baseline_02 = HB[1];
  const std::string& hazard_baseline_12 = HB[2];
  double Gompertz_1_01 = Gompertz[0];
  double Gompertz_2_01 = Gompertz[1];
  double Gompertz_1_02 = Gompertz[2];
  double Gompertz_2_02 = Gompertz[3];
  double Gompertz_1_12 = Gompertz[4];
  double Gompertz_2_12 = Gompertz[5];
  double shape_01 = Weibull[0];
  double shape_02 = Weibull[1];
  double shape_12 = Weibull[2];
  int S = nb_points_integral[0];
  int nb_pointsGK = nb_points_integral[1];
  double alpha_inter_01 = alpha_inter_intra[0];
  double alpha_inter_02 = alpha_inter_intra[1];
  double alpha_inter_12 = alpha_inter_intra[2];
  double alpha_intra_01 = alpha_inter_intra[3];
  double alpha_intra_02 = alpha_inter_intra[4];
  double alpha_intra_12 = alpha_inter_intra[5];
  double alpha_y_01 = alpha_y_slope[0];
  double alpha_y_02 = alpha_y_slope[1];
  double alpha_y_12 = alpha_y_slope[2];
  double alpha_slope_01 = alpha_y_slope[3];
  double alpha_slope_02 = alpha_y_slope[4];
  double alpha_slope_12 = alpha_y_slope[5];
  arma::vec alpha_z_01 = alpha_z[0];
  arma::vec alpha_z_02 = alpha_z[1];
  arma::vec alpha_z_12 = alpha_z[2];
  //Rcout << "The value of v : \n" << 3 << "\n";
  arma::vec gamma_01 = gamma_B[0];
  arma::vec gamma_02 = gamma_B[1];
  arma::vec gamma_12 = gamma_B[2];


  // Survival part
  ///// h
  arma::vec h_02_T_i(S,fill::ones);
  arma::vec etaBaseline_01_T_i(S,fill::zeros);
  arma::mat survLong_01_T_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_02_T_i(S,fill::zeros);
  arma::mat survLong_02_T_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_02_T0_i(S,fill::zeros);
  arma::mat survLong_02_T0_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_01_T0_i(S,fill::zeros);
  arma::mat survLong_01_T0_i(S,nb_pointsGK,fill::zeros);

  arma::vec h_12_T_i(S,fill::ones);
  arma::vec etaBaseline_12_T_i(S,fill::zeros);
  arma::mat survLong_12_T_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_01_L_T_i(S,fill::zeros);
  arma::mat survLong_01_L_T_i(S,nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_01_0_LT_i(S,fill::zeros);
  arma::mat survLong_01_0_LT_i(S,nb_pointsGK*nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_02_0_LT_i(S,fill::zeros);
  arma::mat survLong_02_0_LT_i(S,nb_pointsGK*nb_pointsGK,fill::zeros);
  arma::vec etaBaseline_12_0_LT_i(S,fill::zeros);
  arma::mat survLong_12_0_LT_i(S,nb_pointsGK*nb_pointsGK,fill::zeros);

  arma::mat CV_T;
  arma::mat current_GK_T;
  arma::mat slope_T;
  arma::mat slope_GK_T;
  arma::mat current_GK_T0;
  arma::mat slope_GK_T0;

  arma::mat current_GK_L_T;
  arma::mat current_GK_0_LT;
  arma::mat slope_GK_L_T;
  arma::mat slope_GK_0_LT;

  if(dep_var_inter_12){
    h_12_T_i = h_12_T_i%exp(alpha_inter_12*sigma_inter);
    etaBaseline_12_T_i = etaBaseline_12_T_i + alpha_inter_12*sigma_inter;
    etaBaseline_12_0_LT_i = etaBaseline_12_0_LT_i + alpha_inter_12*sigma_inter;
  }
  if(dep_var_inter_02){
    h_02_T_i = h_02_T_i%exp(alpha_inter_02*sigma_inter);
    etaBaseline_02_T_i = etaBaseline_02_T_i + alpha_inter_02*sigma_inter;
    etaBaseline_02_0_LT_i = etaBaseline_02_0_LT_i + alpha_inter_02*sigma_inter;
    if(left_trunc){
      etaBaseline_02_T0_i = etaBaseline_02_T0_i + alpha_inter_02*sigma_inter;
    }
  }
  if(dep_var_inter_01){
    etaBaseline_01_T_i = etaBaseline_01_T_i + alpha_inter_01*sigma_inter;
    etaBaseline_01_L_T_i = etaBaseline_01_L_T_i + alpha_inter_01*sigma_inter;
    etaBaseline_01_0_LT_i = etaBaseline_01_0_LT_i + alpha_inter_01*sigma_inter;
    if(left_trunc){
      etaBaseline_01_T0_i = etaBaseline_01_T0_i + alpha_inter_01*sigma_inter;
    }
  }
  if(dep_var_intra_12){
    h_12_T_i = h_12_T_i%exp(alpha_intra_12*sigma_intra);
    etaBaseline_12_T_i = etaBaseline_12_T_i + alpha_intra_12*sigma_intra;
    etaBaseline_12_0_LT_i = etaBaseline_12_0_LT_i + alpha_intra_12*sigma_intra;
  }
  if(dep_var_intra_02){
    h_02_T_i = h_02_T_i%exp(alpha_intra_02*sigma_intra);
    etaBaseline_02_T_i = etaBaseline_02_T_i + alpha_intra_02*sigma_intra;
    etaBaseline_02_0_LT_i = etaBaseline_02_0_LT_i + alpha_intra_02*sigma_intra;
    if(left_trunc){
      etaBaseline_02_T0_i = etaBaseline_02_T0_i + alpha_intra_02*sigma_intra;
    }
  }
  if(dep_var_intra_01){
    etaBaseline_01_T_i = etaBaseline_01_T_i + alpha_intra_01*sigma_intra;
    etaBaseline_01_L_T_i = etaBaseline_01_L_T_i + alpha_intra_01*sigma_intra;
    etaBaseline_01_0_LT_i = etaBaseline_01_0_LT_i + alpha_intra_01*sigma_intra;
    if(left_trunc){
      etaBaseline_01_T0_i = etaBaseline_01_T0_i + alpha_intra_01*sigma_intra;
    }
  }
  if(dep_cv_01 || dep_cv_02 || dep_cv_12){
    CV_T = arma::dot(beta, X_T_i) + b_y*U_T_i;
    current_GK_T = arma::repmat(beta.t()*X_GK_T_i.t(),S,1)+b_y*U_GK_T_i.t();
    current_GK_L_T = arma::repmat(beta.t()*X_GK_L_T_i.t(),S,1)+b_y*U_GK_L_T_i.t();
    current_GK_0_LT = arma::repmat(beta.t()*X_GK_0_LT_i.t(),S,1)+b_y*U_GK_0_LT_i.t();
    if(left_trunc){
      current_GK_T0 = arma::repmat(beta.t()*X_GK_T0_i.t(),S,1)+b_y*U_GK_T0_i.t();
    }
    if(dep_cv_01){
      survLong_01_T_i = survLong_01_T_i + alpha_y_01*current_GK_T;
      survLong_01_L_T_i = survLong_01_L_T_i + alpha_y_01*current_GK_L_T;
      survLong_01_0_LT_i = survLong_01_0_LT_i + alpha_y_01*current_GK_0_LT;
      if(left_trunc){
        survLong_01_T0_i = survLong_01_T0_i + alpha_y_01*current_GK_T0;
      }
    }
    if(dep_cv_02){
      h_02_T_i = h_02_T_i%exp(alpha_y_02*CV_T);
      survLong_02_T_i = survLong_02_T_i + alpha_y_02*current_GK_T;
      survLong_02_0_LT_i = survLong_02_0_LT_i + alpha_y_02*current_GK_0_LT;
      if(left_trunc){
        survLong_02_T0_i = survLong_02_T0_i + alpha_y_02*current_GK_T0;
      }
    }
    if(dep_cv_12){
      survLong_12_0_LT_i = survLong_12_0_LT_i + alpha_y_12*current_GK_0_LT;
      survLong_12_T_i = survLong_12_T_i + alpha_y_12*current_GK_T;
      h_12_T_i = h_12_T_i%exp(alpha_y_12*CV_T);
    }
  }
  if(dep_slope_01 || dep_slope_02 || dep_slope_12){
    slope_T = arma::dot(beta_slope, Xslope_T_i)+b_y_slope*Uslope_T_i;
    slope_GK_T = arma::repmat(beta_slope.t()*Xslope_GK_T_i.t(),S,1) + b_y_slope*Uslope_GK_T_i.t();
    slope_GK_L_T = arma::repmat(beta_slope.t()*Xslope_GK_L_T_i.t(),S,1) + b_y_slope*Uslope_GK_L_T_i.t();
    slope_GK_0_LT = arma::repmat(beta_slope.t()*Xslope_GK_0_LT_i.t(),S,1) + b_y_slope*Uslope_GK_0_LT_i.t();
    if(left_trunc){
      slope_GK_T0 = arma::repmat(beta_slope.t()*Xslope_GK_T0_i,S,1)+b_y_slope*Uslope_GK_T0_i.t();
    }
    if(dep_slope_01){
      survLong_01_T_i = survLong_01_T_i + alpha_slope_01*slope_GK_T;
      survLong_01_L_T_i = survLong_01_L_T_i + alpha_slope_01*slope_GK_L_T;
      survLong_01_0_LT_i = survLong_01_0_LT_i + alpha_slope_01*slope_GK_0_LT;
      if(left_trunc){
        survLong_01_T0_i = survLong_01_T0_i + alpha_slope_01*slope_GK_T0;
      }
    }
    if(dep_slope_02){
      h_02_T_i = h_02_T_i%exp(alpha_slope_02*slope_T);
      survLong_02_T_i = survLong_02_T_i + alpha_slope_02*slope_GK_T;
      survLong_02_0_LT_i = survLong_02_0_LT_i + alpha_slope_02*slope_GK_0_LT;
      if(left_trunc){
        survLong_02_T0_i = survLong_02_T0_i + alpha_slope_02*slope_GK_T0;
      }
    }
    if(dep_slope_12){
      survLong_12_0_LT_i = survLong_12_0_LT_i + alpha_slope_12*slope_GK_0_LT;
      survLong_12_T_i = survLong_12_T_i + alpha_slope_12*slope_GK_T;
      h_12_T_i = h_12_T_i%exp(alpha_slope_12*slope_T);
    }
  }
  ///// h0
  ///////// 1-2
  double h_0_12_T_i;
  arma::vec h_0_GK_12_T_i;
  arma::vec h_0_GK_12_0_LT_i;
  if(hazard_baseline_12 == "Exponential"){
    h_0_12_T_i = 1;
    h_0_GK_12_T_i = wk;
    h_0_GK_12_0_LT_i = rep_wk;
  }
  if(hazard_baseline_12 == "Weibull"){
    h_0_12_T_i = shape_12*(pow(Time_T_i,(shape_12-1)));
    h_0_GK_12_T_i = shape_12*(pow(st_T_i,shape_12-1))%wk;
    h_0_GK_12_0_LT_i = shape_12*(pow(st_0_LT_i,shape_12-1));
    h_0_GK_12_0_LT_i = h_0_GK_12_0_LT_i.t()%rep_wk;
  }
  if(hazard_baseline_12 == "Gompertz"){
    h_0_12_T_i = Gompertz_1_12*exp(Gompertz_2_12*Time_T_i);
    h_0_GK_12_T_i = Gompertz_1_12*exp(Gompertz_2_12*st_T_i)%wk;
    h_0_GK_12_0_LT_i = Gompertz_1_12*exp(Gompertz_2_12*st_0_LT_i);
    h_0_GK_12_0_LT_i = h_0_GK_12_0_LT_i.t()%rep_wk;
  }
  if(hazard_baseline_12 == "Splines"){
    h_0_12_T_i = exp(arma::dot(gamma_12,B_T_i_12));
    h_0_GK_12_T_i = wk%exp(Bs_T_i_12*gamma_12);
    h_0_GK_12_0_LT_i = exp(Bs_0_LT_i_12*gamma_12)%rep_wk;
  }
  double predsurv_12;
  if(Z_12_i.is_empty()){
    predsurv_12 = 0;
  }
  else{
    predsurv_12 = arma::dot(alpha_z_12, Z_12_i);
  }

  h_12_T_i = h_0_12_T_i*exp(predsurv_12)*h_12_T_i;

  etaBaseline_12_T_i = etaBaseline_12_T_i + predsurv_12;
  survLong_12_T_i = exp(survLong_12_T_i)*h_0_GK_12_T_i;
  arma::vec A_12_T_i;
  A_12_T_i = (exp(etaBaseline_12_T_i)%survLong_12_T_i*(Time_T_i/2)); //% : multiplication element par element

  etaBaseline_12_0_LT_i = exp(etaBaseline_12_0_LT_i + predsurv_12);
  survLong_12_0_LT_i = exp(survLong_12_0_LT_i);
  survLong_12_0_LT_i = survLong_12_0_LT_i%arma::repelem(h_0_GK_12_0_LT_i.t(),S,1);
  arma::mat A_12_0_LT_i;
  A_12_0_LT_i = arma::repelem(etaBaseline_12_0_LT_i,1,nb_pointsGK*nb_pointsGK)%survLong_12_0_LT_i;

  ///////// 0-1
  arma::vec h_0_GK_01_L_T_i;
  arma::vec h_0_GK_01_0_LT_i;
  arma::vec h_0_GK_01_T0_i;
  arma::vec h_0_GK_01_T_i;
  if(hazard_baseline_01 == "Exponential"){
    h_0_GK_01_L_T_i = wk;
    h_0_GK_01_T_i = wk;
    h_0_GK_01_0_LT_i = rep_wk;
    if(left_trunc){
      h_0_GK_01_T0_i = wk;
    }
  }
  if(hazard_baseline_01 == "Weibull"){
    h_0_GK_01_T_i = shape_01*(pow(st_T_i,shape_01-1))%wk;
    h_0_GK_01_L_T_i = shape_01*(pow(st_L_T_i,shape_01-1))%wk;
    h_0_GK_01_0_LT_i = shape_01*(pow(st_0_LT_i,shape_01-1));
    h_0_GK_01_0_LT_i = h_0_GK_01_0_LT_i.t()%rep_wk;
    if(left_trunc){
      h_0_GK_01_T0_i = shape_01*(pow(st_T0_i,shape_01-1));
    }
  }
  if(hazard_baseline_01 == "Gompertz"){
    h_0_GK_01_T_i = Gompertz_1_01*exp(Gompertz_2_01*st_T_i)%wk;
    h_0_GK_01_L_T_i = Gompertz_1_01*exp(Gompertz_2_01*st_L_T_i)%wk;
    h_0_GK_01_0_LT_i = Gompertz_1_01*exp(Gompertz_2_01*st_0_LT_i);
    h_0_GK_01_0_LT_i = h_0_GK_01_0_LT_i.t()%rep_wk;
    if(left_trunc){
      h_0_GK_01_T0_i = Gompertz_1_01*exp(st_T0_i*Gompertz_2_01)%wk;
    }
  }
  if(hazard_baseline_01 == "Splines"){
    h_0_GK_01_T_i = wk%exp(Bs_T_i_01*gamma_01);
    h_0_GK_01_L_T_i = wk%exp(Bs_L_T_i_01*gamma_01);
    h_0_GK_01_0_LT_i = exp(Bs_0_LT_i_01*gamma_01)%rep_wk;
    if(left_trunc){
      h_0_GK_01_T0_i = wk%exp(Bs_T0_i_01*gamma_01);
    }
  }
  double predsurv_01;
  if(Z_01_i.is_empty()){
    predsurv_01 = 0;
  }
  else{
    predsurv_01 = arma::dot(alpha_z_01, Z_01_i);
  }
  etaBaseline_01_L_T_i = exp(etaBaseline_01_L_T_i + predsurv_01);
  survLong_01_L_T_i = exp(survLong_01_L_T_i)%arma::repelem(h_0_GK_01_L_T_i.t(),S,1);
  arma::mat A_01_L_T_i;
  A_01_L_T_i = arma::repelem(etaBaseline_01_L_T_i,1,nb_pointsGK)%survLong_01_L_T_i*(Time_L_T_i/2);

  etaBaseline_01_T_i = exp(etaBaseline_01_T_i + predsurv_01);
  survLong_01_T_i = exp(survLong_01_T_i)*h_0_GK_01_T_i;
  arma::mat A_01_T_i;
  A_01_T_i = etaBaseline_01_T_i%survLong_01_T_i*(Time_T_i/2);

  etaBaseline_01_0_LT_i = exp(etaBaseline_01_0_LT_i + predsurv_01);
  survLong_01_0_LT_i = exp(survLong_01_0_LT_i);
  survLong_01_0_LT_i = survLong_01_0_LT_i%arma::repelem(h_0_GK_01_0_LT_i.t(),S,1);
  arma::mat A_01_0_LT_i;
  A_01_0_LT_i = arma::repelem(etaBaseline_01_0_LT_i,1,nb_pointsGK*nb_pointsGK)%survLong_01_0_LT_i;

  arma::vec A_01_T0_i;
  if(left_trunc){
    etaBaseline_01_T0_i = etaBaseline_01_T0_i + predsurv_01;
    survLong_01_T0_i = exp(survLong_01_T0_i)*h_0_GK_01_T0_i;
    A_01_T0_i = (exp(etaBaseline_01_T0_i)%survLong_01_T0_i*(Time_T0_i/2));
  }

  ///////// 0-2
  arma::vec h_0_GK_02_0_LT_i;
  arma::vec h_0_GK_02_T0_i;
  arma::vec h_0_GK_02_T_i;
  double h_0_02_T_i;
  if(hazard_baseline_02 == "Exponential"){
    h_0_02_T_i = 1;
    h_0_GK_02_T_i = wk;
    h_0_GK_02_0_LT_i = rep_wk;
    if(left_trunc){
      h_0_GK_02_T0_i = wk;
    }
  }
  if(hazard_baseline_02 == "Weibull"){
    h_0_02_T_i = shape_02*(pow(Time_T_i,(shape_02-1)));
    h_0_GK_02_T_i = shape_02*(pow(st_T_i,shape_02-1))%wk;
    h_0_GK_02_0_LT_i = shape_02*(pow(st_0_LT_i,shape_02-1));
    h_0_GK_02_0_LT_i = h_0_GK_02_0_LT_i.t()%rep_wk;
    if(left_trunc){
      h_0_GK_02_T0_i = shape_02*(pow(st_T0_i,shape_02-1));
    }
  }
  if(hazard_baseline_02 == "Gompertz"){
    h_0_02_T_i = Gompertz_1_02*exp(Gompertz_2_02*Time_T_i);
    h_0_GK_02_T_i = Gompertz_1_02*exp(Gompertz_2_02*st_T_i)%wk;
    h_0_GK_02_0_LT_i = Gompertz_1_02*exp(Gompertz_2_02*st_0_LT_i);
    h_0_GK_02_0_LT_i = h_0_GK_02_0_LT_i.t()%rep_wk;
    if(left_trunc){
      h_0_GK_02_T0_i = Gompertz_1_02*exp(st_T0_i*Gompertz_2_02)%wk;
    }
  }
  if(hazard_baseline_02 == "Splines"){
    h_0_02_T_i = exp(arma::dot(gamma_02,B_T_i_02));
    h_0_GK_02_T_i = wk%exp(Bs_T_i_02*gamma_02);
    h_0_GK_02_0_LT_i = exp(Bs_0_LT_i_02*gamma_02)%rep_wk;
    if(left_trunc){
      h_0_GK_02_T0_i = wk%exp(Bs_T0_i_02*gamma_02);
    }
  }
  double predsurv_02;
  if(Z_02_i.is_empty()){
    predsurv_02 = 0;
  }
  else{
    predsurv_02 = arma::dot(alpha_z_02, Z_02_i);
  }


  etaBaseline_02_T_i = exp(etaBaseline_02_T_i + predsurv_02);
  survLong_02_T_i = exp(survLong_02_T_i)*h_0_GK_02_T_i;
  arma::mat A_02_T_i;
  A_02_T_i = etaBaseline_02_T_i%survLong_02_T_i*(Time_T_i/2);

  etaBaseline_02_0_LT_i = exp(etaBaseline_02_0_LT_i + predsurv_02);
  survLong_02_0_LT_i = exp(survLong_02_0_LT_i);
  survLong_02_0_LT_i = survLong_02_0_LT_i%arma::repelem(h_0_GK_02_0_LT_i.t(),S,1);
  arma::mat A_02_0_LT_i;
  A_02_0_LT_i = arma::repelem(etaBaseline_02_0_LT_i,1,nb_pointsGK*nb_pointsGK)%survLong_02_0_LT_i;

  h_02_T_i = h_0_02_T_i*exp(predsurv_02)*h_02_T_i;

  arma::vec A_02_T0_i;
  if(left_trunc){
    etaBaseline_02_T0_i = etaBaseline_02_T0_i + predsurv_02;
    survLong_02_T0_i = exp(survLong_02_T0_i)*h_0_GK_02_T0_i;
    A_02_T0_i = (exp(etaBaseline_02_T0_i)%survLong_02_T0_i*(Time_T0_i/2));
  }

  arma::mat A_0_LT = -A_01_0_LT_i - A_02_0_LT_i + A_12_0_LT_i;
  arma::mat A_0_LT_red(S,nb_pointsGK,fill::zeros);
  for(int column = 0; column < nb_pointsGK; ++column){
    A_0_LT_red.col(column) = sum(A_0_LT.cols(column*nb_pointsGK, (column+1)*nb_pointsGK-1),  1);
  }
  A_0_LT_red = A_0_LT_red%arma::repelem(ck.t(),S,1);
  arma::mat A_12_T_i_rep = arma::repelem(A_12_T_i,1,nb_pointsGK);
  arma::vec SurvTotCase3 = log(exp(-A_01_T_i - A_02_T_i)%(pow(h_02_T_i,delta2_i))+(pow(h_12_T_i,delta2_i))%sum(A_01_L_T_i%exp(A_0_LT_red-A_12_T_i_rep),1));

  double den = 0;
  if(left_trunc){
    den = log(sum(exp(-A_01_T0_i - A_02_T0_i)))-log(S);
  }




  List ret;
  ret["SurvTotCase3"] = SurvTotCase3;

  ret["den"] = den;
  return ret;









}
