#' log_llh_2var_IC_rcpp
#'
#' @param param
#' @param hazard_baseline_01
#' @param sharedtype_01
#' @param hazard_baseline_02
#' @param sharedtype_02
#' @param hazard_baseline_12
#' @param sharedtype_12
#' @param ord.splines
#' @param nb.beta
#' @param Zq
#' @param nb_pointsGK
#' @param nb.e.a
#' @param S
#' @param wk
#' @param sk_GK
#' @param nb.alpha
#' @param variability_inter_visit
#' @param variability_intra_visit
#' @param correlated_re
#' @param Case1
#' @param Case1bis
#' @param Case2
#' @param Case3
#' @param nbCase1
#' @param nbCase1bis
#' @param nbCase2
#' @param nbCase3
#' @param left_trunc
#' @param knots.hazard_baseline.splines_01
#' @param knots.hazard_baseline.splines_02
#' @param knots.hazard_baseline.splines_12
#' @param index_beta_slope
#' @param control
#'
#' @return
#' @export
#'
#' @examples
log_llh_2var_IC_rcpp <- function(param,hazard_baseline_01, sharedtype_01,
                                 hazard_baseline_02, sharedtype_02,
                                 hazard_baseline_12, sharedtype_12,
                                 ord.splines, nb.beta, Zq, nb_pointsGK,
                                 nb.e.a, S, wk, rep_wk, sk_GK, nb.alpha,
                                 variability_inter_visit, variability_intra_visit,
                                 correlated_re, Case1, Case1bis, Case2, Case3,
                                 nbCase1, nbCase1bis, nbCase2, nbCase3, left_trunc,
                                 knots.hazard_baseline.splines_01,
                                 knots.hazard_baseline.splines_02,
                                 knots.hazard_baseline.splines_12,
                                 index_beta_slope,
                                 control){

  # Initialiser certains paramètres
  Gompertz.1_01 <- 0; Gompertz.2_01 <- 0; Gompertz.1_02 <- 0; Gompertz.2_02 <- 0; Gompertz.1_12 <- 0; Gompertz.2_12 <- 0
  shape_01 <- 0; shape_02 <- 0; shape_12 <- 0
  alpha.inter_01 <- 0; alpha.inter_02 <- 0; alpha.inter_12 <- 0
  alpha.intra_01 <- 0; alpha.intra_02 <- 0; alpha.intra_12 <- 0
  alpha.current_01 <- 0; alpha.current_02 <- 0; alpha.current_12 <- 0
  alpha.slope_01 <- 0; alpha.slope_02 <- 0; alpha.slope_12 <- 0
  alpha_01 <- c(0); alpha_02 <- c(0); alpha_12 <- c(0);
  gamma_01 <- c(0,3); gamma_02 <- c(0); gamma_12 <- c(0);
  beta_slope <- c(0); b_y_slope <- as.matrix(1); sigma_inter <- c(0); sigma_intra <- c(0);
  Xslope_T_i <- c(0); Uslope_T_i <- c(0); Xslope_GK_T_i <- as.matrix(1); Uslope_GK_T_i <- as.matrix(1);
  Xslope_L_i <- c(0); Uslope_L_i <- c(0); Xslope_GK_L_i <- as.matrix(1); Uslope_GK_L_i <- as.matrix(1);
  Xslope_GK_L_R_i <- as.matrix(1); Uslope_GK_L_R_i <- as.matrix(1);
  Xslope_GK_L_T_i <- as.matrix(1); Uslope_GK_L_T_i <- as.matrix(1);
  Xslope_GK_0_LR_i <- as.matrix(1); Uslope_GK_0_LR_i <- as.matrix(1); X_GK_T0_i<- as.matrix(1); U_GK_T0_i<- as.matrix(1); Xslope_GK_T0_i<- as.matrix(1); Uslope_GK_T0_i<- as.matrix(1);
  Xslope_GK_0_LT_i <- as.matrix(1); Uslope_GK_0_LT_i <- as.matrix(1);
  st_T0_i <- c(0); B_T_i_01 <- c(0); B_T_i_02 <- c(0); B_T_i_12 <- c(0); B_L_i_02 <- c(0); B_L_i_12 <- c(0);
  Bs_T_i_01 <- as.matrix(1);Bs_T_i_02 <- as.matrix(1); Bs_T_i_12 <- as.matrix(1); Bs_L_i_01 <- as.matrix(1);Bs_L_i_02 <- as.matrix(1); Bs_L_i_12 <- as.matrix(1);
  Bs_0_LR_i_01 <- as.matrix(1); Bs_0_LR_i_02 <- as.matrix(1); Bs_0_LR_i_12 <- as.matrix(1);
  Bs_0_LT_i_01 <- as.matrix(1); Bs_0_LT_i_02 <- as.matrix(1); Bs_0_LT_i_12 <- as.matrix(1);
  Bs_L_R_i_01 <- as.matrix(1); Bs_L_R_i_02 <- as.matrix(1); Bs_L_R_i_12<- as.matrix(1);
  Bs_L_T_i_01 <- as.matrix(1); Bs_L_T_i_02 <- as.matrix(1); Bs_L_T_i_12<- as.matrix(1);
  Bs_T0_i_01 <- as.matrix(1); Bs_T0_i_02 <- as.matrix(1); Bs_T0_i_12 <- as.matrix(1); Time_T0_i <- 0;
  st_T_i <- c(0); st_0_LR_i <- as.matrix(1); st_L_R_i <- c(0); st_T0_i <- c(0); st_0_LT_i <- as.matrix(1); st_L_T_i <- c(0);
  #Manage parameter
  curseur <- 1
  ## Risque 01
  ### Hazard baseline
  if(hazard_baseline_01 == "Weibull"){
    shape_01 <- param[curseur]**2
    curseur <- curseur + 1
  }

  if(hazard_baseline_01 == "Gompertz"){
    Gompertz.1_01 <- param[curseur]**2
    Gompertz.2_01 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(hazard_baseline_01 == "Splines"){
    gamma_01 <- param[(curseur):(curseur+ord.splines+1)]
    curseur <- curseur + ord.splines + 2
  }
  ### Covariables :
  nb.alpha_01 <- nb.alpha[1]
  if(nb.alpha_01 >=1){
    alpha_01 <-  param[(curseur):(curseur+nb.alpha_01-1)]
    curseur <- curseur+nb.alpha_01
  }
  ### Association
  if("current value" %in% sharedtype_01){
    alpha.current_01 <-  param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% sharedtype_01){
    alpha.slope_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("inter visit variability" %in% sharedtype_01){
    alpha.inter_01 <- param[curseur]
    curseur <- curseur + 1
  }
  if("intra visit variability" %in% sharedtype_01){
    alpha.intra_01 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Risque 02
  if(hazard_baseline_02 == "Weibull"){
    shape_02 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(hazard_baseline_02 == "Gompertz"){
    Gompertz.1_02 <- param[curseur]**2
    Gompertz.2_02 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(hazard_baseline_02 == "Splines"){
    gamma_02 <- param[(curseur):(curseur+ord.splines+1)]
    curseur <- curseur + ord.splines + 2
  }
  ### Covariables :
  nb.alpha_02 <- nb.alpha[2]
  if(nb.alpha_02 >=1){
    alpha_02 <-  param[(curseur):(curseur+nb.alpha_02-1)]
    curseur <- curseur+nb.alpha_02
  }
  ### Association
  if("current value" %in% sharedtype_02){
    alpha.current_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% sharedtype_02){
    alpha.slope_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("inter visit variability" %in% sharedtype_02){
    alpha.inter_02 <- param[curseur]
    curseur <- curseur + 1
  }
  if("intra visit variability" %in% sharedtype_02){
    alpha.intra_02 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Risque 12
  if(hazard_baseline_12 == "Weibull"){
    shape_12 <- param[curseur]**2
    curseur <- curseur + 1
  }
  if(hazard_baseline_12 == "Gompertz"){
    Gompertz.1_12 <- param[curseur]**2
    Gompertz.2_12 <- param[curseur+1]
    curseur <- curseur + 2
  }
  if(hazard_baseline_12 == "Splines"){
    gamma_12 <- param[(curseur):(curseur+ord.splines+1)]
    curseur <- curseur + ord.splines + 2
  }
  ### Covariables :
  nb.alpha_12 <- nb.alpha[3]
  if(nb.alpha_12 >=1){
    alpha_12 <- param[(curseur):(curseur+nb.alpha_12-1)]
    curseur <- curseur+nb.alpha_12
  }
  ### Association
  if("current value" %in% sharedtype_12){
    alpha.current_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("slope" %in% sharedtype_12){
    alpha.slope_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("inter visit variability" %in% sharedtype_12){
    alpha.inter_12 <- param[curseur]
    curseur <- curseur + 1
  }
  if("intra visit variability" %in% sharedtype_12){
    alpha.intra_12 <- param[curseur]
    curseur <- curseur + 1
  }
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+nb.beta-1)]
  if( "slope" %in% sharedtype_01 || "slope" %in% sharedtype_02 || "slope" %in% sharedtype_12){
    beta_slope <- beta[index_beta_slope]
  }
  curseur <- curseur+nb.beta
  if(variability_inter_visit){
    mu.inter <- param[curseur]
    curseur <- curseur + 1
  }
  if(variability_intra_visit){
    mu.intra <- param[curseur]
    curseur <- curseur + 1
  }
  if(!variability_inter_visit && !variability_intra_visit){
    sigma.epsilon <- param[curseur]
    curseur <- curseur +1
  }
  ### Cholesky matrix for random effects
  if(variability_inter_visit && variability_intra_visit){
    if(correlated_re){
      C1 <- matrix(rep(0,(length(param)-curseur+1)**2),nrow=length(param)-curseur+1,ncol=length(param)-curseur+1)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }
    else{
      borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
      C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
      C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
      C2 <-matrix(c(param[(borne1+1)], 0,param[borne1+2], param[borne1+3]),nrow=2,ncol=2, byrow = TRUE)
      C3 <- matrix(rep(0,2*nb.e.a), ncol = nb.e.a)
      C4 <- matrix(rep(0,2*nb.e.a), nrow = nb.e.a)
      Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
      Cholesky <- as.matrix(Cholesky)
    }
  }
  else{
    if(variability_inter_visit){
      if(correlated_re){
        C1 <- matrix(rep(0,(length(param)-curseur+1)**2),nrow=length(param)-curseur+1,ncol=length(param)-curseur+1)
        C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
        Cholesky <- C1
      }
      else{
        borne1 <- curseur + choose(n = nb.e.a, k = 2) + nb.e.a - 1
        C1 <- matrix(rep(0,(nb.e.a)**2),nrow=nb.e.a,ncol=nb.e.a)
        C1[lower.tri(C1, diag=T)] <- param[curseur:borne1]
        C2 <-matrix(c(param[(borne1+1)]),nrow=1,ncol=1, byrow = TRUE)
        C3 <- matrix(rep(0,nb.e.a), ncol = 1)
        C4 <- matrix(rep(0,nb.e.a), nrow = 1)
        Cholesky <- rbind(cbind(C1,C4),cbind(C3,C2))
        Cholesky <- as.matrix(Cholesky)
      }
    }
    if(!variability_inter_visit && !variability_intra_visit){
      C1 <- matrix(rep(0,(length(param)-curseur)**2),nrow=length(param)-curseur,ncol=length(param)-curseur)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }
    else{
      stop("not implemented yet")
    }
  }

  # Manage random effects
  random.effects <- Zq%*%t(Cholesky)
  b_y <- random.effects[,1:nb.e.a]
  b_y <- matrix(b_y, ncol = nb.e.a)
  if("slope" %in% sharedtype_01 || "slope" %in% sharedtype_02 || "slope" %in% sharedtype_12 ){
    b_y_slope <- as.matrix(b_y[,-1])
  }
  if(variability_inter_visit){
    b_inter <- random.effects[,nb.e.a+1]
    sigma_inter <- exp(mu.inter + b_inter)
    var.inter <- sigma_inter**2
  }
  if(variability_intra_visit){
    b_intra <- random.effects[,nb.e.a+2]
    sigma_intra <- exp(mu.intra + b_intra)
    var.intra <- sigma_intra**2
  }
  ll_glob <- 0

  # Creations entrees rcpp
  sharedtype <- c("current value" %in% sharedtype_01, "slope" %in% sharedtype_01, "inter visit variability" %in% sharedtype_01, "intra visit variability" %in% sharedtype_01,
                  "current value" %in% sharedtype_02, "slope" %in% sharedtype_02, "inter visit variability" %in% sharedtype_02, "intra visit variability" %in% sharedtype_02,
                  "current value" %in% sharedtype_12, "slope" %in% sharedtype_12, "inter visit variability" %in% sharedtype_12, "intra visit variability" %in% sharedtype_12)
  HB <- list(hazard_baseline_01, hazard_baseline_02, hazard_baseline_12)
  Weibull <- c(shape_01, shape_02, shape_12)
  Gompertz <- c(Gompertz.1_01, Gompertz.2_01, Gompertz.1_02, Gompertz.2_02, Gompertz.1_12, Gompertz.2_12)
  alpha_inter_intra <- c(alpha.inter_01,alpha.inter_02,alpha.inter_12, alpha.intra_01,alpha.intra_02,alpha.intra_12)
  alpha_y_slope <- c(alpha.current_01,alpha.current_02,alpha.current_12, alpha.slope_01,alpha.slope_02,alpha.slope_12)
  alpha_z <- list(alpha_01, alpha_02, alpha_12)
  nb_points_integral <- c(S, nb_pointsGK)
  gamma_z0 <- list(gamma_01, gamma_02, gamma_12)

  ### Case 1
  if(nbCase1 != 0){
    print("Case1 go")
    delta2 <- Case1[["delta2"]]; Z_12 <- Case1[["Z_12"]]; Time_T <- Case1[["Time_T"]]; st_T <- Case1[["st_T"]]
    X_GK_T <- Case1[["X_GK_T"]]; U_GK_T <- Case1[["U_GK_T"]]
    Xslope_GK_T <- Case1[["Xslope_GK_T"]];Uslope_GK_T <- Case1[["Uslope_GK_T"]]
    X_T <- Case1[["X_T"]]; U_T <- Case1[["U_T"]]; Xslope_T <- Case1[["Xslope_T"]]; Uslope_T <- Case1[["Uslope_T"]]
    X_base <- Case1[["X_base"]]; U_base <- Case1[["U_base"]]; y.new <- Case1[["y.new"]];
    ID.visit <- Case1[["ID.visit"]]; offset <- Case1[["offset"]];
    B_T_01 <- Case1[["B_T_01"]]; Bs_T_01 <- Case1[["Bs_T_01"]]
    B_T_02 <- Case1[["B_T_02"]]; Bs_T_02 <- Case1[["Bs_T_02"]]
    B_T_12 <- Case1[["B_T_12"]]; Bs_T_12 <- Case1[["Bs_T_12"]]
    if(left_trunc){
      Time_T0 <- Case1[["Time_T0"]]; st_T0 <- Case1[["st_T0"]]; Bs_T0_01 <- Case1[["Bs_T0_01"]]
      Bs_T0_02 <- Case1[["Bs_T0_02"]];
      X_GK_T0 <- Case1[["X_GK_T0"]];U_GK_T0 <- Case1[["U_GK_T0"]]
      Xslope_GK_T0 <- Case1[["Xslope_GK_T0"]];Uslope_GK_T0 <- Case1[["Uslope_GK_T0"]]
    }
    ### Intégrale L_i to R_i :
    Z_01 <- Case1[["Z_01"]]; Z_02 <- Case1[["Z_02"]];
    X_GK_L_R <- Case1[["X_GK_L_R"]]; U_GK_L_R <- Case1[["U_GK_L_R"]]
    Xslope_GK_L_R <- Case1[["Xslope_GK_L_R"]];Uslope_GK_L_R <- Case1[["Uslope_GK_L_R"]]
    Bs_L_R_01 <- Case1[["Bs_L_R_01"]]; Bs_L_R_02 <- Case1[["Bs_L_R_02"]]; Bs_L_R_12 <- Case1[["Bs_L_R_12"]];
    Time_L_R <- Case1[["Time_L_R"]]; st_L_R <- Case1[["st_L_R"]] #Time_L_R = R-L # st_L_R = (xk+1)/2*(R-L)+L
    ### Integrale de 0 à u (u dans L_i-R_i)
    st_0_LR <- Case1[["st_0_LR"]]; X_0_LR <- Case1[["X_0_LR"]]; U_0_LR <- Case1[["U_0_LR"]];
    Xslope_0_LR <- Case1[["Xslope_0_LR"]]; Uslope_0_LR <- Case1[["Uslope_0_LR"]]; Time_L <- Case1[["Time_L"]];
    Bs_0_LR_01 <- Case1[["Bs_0_LR_01"]]; Bs_0_LR_02 <- Case1[["Bs_0_LR_02"]]; Bs_0_LR_12 <- Case1[["Bs_0_LR_12"]];

    for(i in 1:nbCase1){
      if("current value" %in% sharedtype_12 || "current value" %in% sharedtype_01 || "current value" %in% sharedtype_02){
        X_T_i <- X_T[i,]; U_T_i <- U_T[i,]
        X_GK_T_i <- X_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]; U_GK_T_i <- U_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        X_GK_L_R_i <- X_GK_L_R[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]; U_GK_L_R_i <- U_GK_L_R[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        X_GK_0_LR_i <- X_0_LR[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        U_GK_0_LR_i <- U_0_LR[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        if(left_trunc && ("current value" %in% sharedtype_01 || "current value" %in% sharedtype_02)){
          X_GK_T0_i <- X_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_T0_i <- U_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }

      }
      if("slope" %in% sharedtype_12 || "slope" %in% sharedtype_01 || "slope" %in% sharedtype_02){
        Xslope_T_i <- Xslope_T[i,];Uslope_T_i <- Uslope_T[i,]
        Xslope_GK_T_i <- Xslope_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_T_i <- Uslope_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Xslope_GK_L_R_i <- Xslope_GK_L_R[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_L_R_i <- Uslope_GK_L_R[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Xslope_GK_0_LR_i <- Xslope_0_LR[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        Uslope_GK_0_LR_i <- Uslope_0_LR[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        if(left_trunc && ("slope" %in% sharedtype_01 || "slope" %in% sharedtype_02)){
          Xslope_GK_T0_i <- Xslope_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_T0_i <- Uslope_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Weibull" %in% c(hazard_baseline_01,hazard_baseline_02, hazard_baseline_12) ||"Gompertz" %in% c(hazard_baseline_01,hazard_baseline_02, hazard_baseline_12)){
        st_T_i <- st_T[i,]
        st_0_LR_i <- st_0_LR[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        st_L_R_i <- st_L_R[i,]
        if(left_trunc){
          st_T0_i <- st_T0[i,]
        }
      }
      if("Splines" %in% hazard_baseline_01){
        Bs_L_R_i_01 <- Bs_L_R_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Bs_0_LR_i_01 <- Bs_0_LR_01[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        if(left_trunc){
          Bs_T0_i_01 <- Bs_T0_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Splines" %in% hazard_baseline_02){
        Bs_0_LR_i_02 <- Bs_0_LR_02[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        if(left_trunc){
          Bs_T0_i_02 <- Bs_T0_02[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Splines" %in% hazard_baseline_12){
        B_T_i_12 <- B_T_12[i,]
        Bs_T_i_12 <- Bs_T_12[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Bs_0_LR_i_12 <- Bs_0_LR_12[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
      }
      Z_01_i <- Z_01[i,]
      Z_02_i <- Z_02[i,]
      Z_12_i <- Z_12[i,]
      Time_T_i <- Time_T[i]
      ck <- ((sk_GK+1)/4)*Time_L_R[i]+Time_L[i]/2
      Time_L_R_i <- Time_L_R[i]
      if(left_trunc){
        Time_T0_i <- Time_T0[i]}
      delta2_i <- delta2[i]
      log_ind_surv <- log_IC_2var_Case1(sharedtype, HB, Gompertz, Weibull,
                                        nb_points_integral, alpha_inter_intra,
                                        alpha_y_slope, alpha_z,  gamma_z0,  beta,  beta_slope,b_y,
                                        b_y_slope,  wk,  rep_wk,  sigma_inter,  sigma_intra,
                                       delta2_i, Z_01_i, Z_02_i, Z_12_i, X_T_i,  U_T_i,
                                       Xslope_T_i,  Uslope_T_i,  X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                       Uslope_GK_T_i,  X_GK_L_R_i,  U_GK_L_R_i,  Xslope_GK_L_R_i,  Uslope_GK_L_R_i,
                                       X_GK_0_LR_i,  U_GK_0_LR_i,  Xslope_GK_0_LR_i,  Uslope_GK_0_LR_i,
                                       X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i,
                                       Time_T_i,  Time_L_R_i,  Time_T0_i, st_T_i,  st_0_LR_i,  st_L_R_i,  st_T0_i,
                                       ck,
                                       B_T_i_01,  B_T_i_02,  B_T_i_12,
                                       Bs_T_i_01,  Bs_T_i_02,  Bs_T_i_12,
                                       Bs_0_LR_i_01,  Bs_0_LR_i_02,  Bs_0_LR_i_12,
                                       Bs_L_R_i_01,   Bs_L_R_i_02,  Bs_L_R_i_12,
                                       Bs_T0_i_01,  Bs_T0_i_02,  Bs_T0_i_12,left_trunc
      )

      ## Longitudinal part
      X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[i+1]-offset[i])
      U_i <- U_base[offset[i]:(offset[i+1]-1),]
      U_i <- matrix(U_i, nrow = offset[i+1]-offset[i])
      y_i <- y.new[offset[i]:(offset[i+1]-1)]
      ID.visit_i <- ID.visit[offset[i]:(offset[i+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      f_Y_b_sigma <- rep(0,S)
      for(id.visit in 1:length(unique(ID.visit_i))){

        X_base_i.id.visit <- X_base_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1),]
        X_base_i.id.visit <- matrix(X_base_i.id.visit, nrow = offset_ID_i[id.visit+1]-offset_ID_i[id.visit])
        X_base_i.id.visit <- matrix(X_base_i.id.visit[1,], nrow = 1)

        U_i.id.visit <- U_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1),]
        U_i.id.visit <- matrix(U_i.id.visit, nrow = offset_ID_i[id.visit+1]-offset_ID_i[id.visit])
        U_i.id.visit <- matrix(U_i.id.visit[1,],nrow=1)

        y_i.id.visit <- y_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1)]

        if(is.null(nrow(X_base_i.id.visit))){
          stop("There is a something wrong")
        }
        else{
          CV <- (X_base_i.id.visit%*%beta)[1,1] + b_y%*%t(U_i.id.visit)
          n_ij <- length(y_i.id.visit)
          if(n_ij == 1){
            f_Y_b_sigma <- f_Y_b_sigma + log(dnorm(y_i.id.visit,CV,var.inter+var.intra))
          }
          else{
            if(variability_inter_visit && variability_intra_visit){
              if(n_ij == 2){
                f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*sqrt(var.intra*(2*var.inter+var.intra)))) -
                  (1/(2*var.intra*(var.intra+2*var.inter)))*((((rep(y_i.id.visit[1],S)-CV)**2)*(var.intra+var.inter))-2*(var.inter*(rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[2],S)-CV)) +
                                                               (((rep(y_i.id.visit[2],S)-CV)**2)*(var.intra+var.inter)))
              }
              else{
                if(n_ij == 3){
                  f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*var.intra*sqrt((3*var.inter+var.intra)))) -
                    (1/(2*(var.intra**2)*(var.intra+3*var.inter)))*((var.intra*(var.intra+2*var.inter)*((rep(y_i.id.visit[1],S)-CV)**2 + (rep(y_i.id.visit[2],S)-CV)**2 + (rep(y_i.id.visit[3],S)-CV)**2))-
                                                                      2*var.inter*var.intra*((rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[2],S)-CV) + (rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[3],S)-CV) + (rep(y_i.id.visit[2],S)-CV)*(rep(y_i.id.visit[3],S)-CV)))
                }
                else{
                  somme1 <- 0
                  somme2 <- 0
                  for(k in 1:n_ij){
                    somme1 <- somme1 + (rep(y_i.id.visit[k],S)-CV)**2
                    if(k != n_ij){
                      for(l in (k+1):n_ij){
                        somme2 <- somme2 + (rep(y_i.id.visit[k],S)-CV)*(rep(y_i.id.visit[l],S)-CV)
                      }
                    }
                  }
                  f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*sqrt((var.intra**(n_ij-1))*(var.intra+n_ij*var.inter)))) -
                    (1/(2*(var.intra**(n_ij-1))*(var.intra+n_ij*var.inter)))*((var.intra**(n_ij-2))*(var.intra+(n_ij-1)*var.inter)*somme1 -
                                                                                2*var.inter*(var.intra**(n_ij-2))*somme2)
                }
              }
            }
            else{
              stop("Not implemented in this program.")
            }
          }


        }
      }

      log_dens_int <- f_Y_b_sigma + log_ind_surv$SurvTotCase1
      Clogexp <- max(log_dens_int) - 500
      log_dens_int <- log_dens_int - Clogexp
      log_dens <- Clogexp + log(sum(exp(log_dens_int))) - log(S)
      ## Left truncation
      if(left_trunc){
        log_dens <- log_dens - log_ind_surv$den
      }
      ll_glob <- ll_glob + log_dens
    }
    print("end Case 1")
  }

  if(nbCase1bis != 0){
    print("Case1bis go")
    delta2 <- Case1bis[["delta2"]]; Z_12 <- Case1bis[["Z_12"]]; Z_01 <- Case1bis[["Z_01"]]; Z_02 <- Case1bis[["Z_02"]]
    Time_T <- Case1bis[["Time_T"]];Time_L <- Case1bis[["Time_L"]]
    st_T <- Case1bis[["st_T"]];st_L <- Case1bis[["st_L"]]
    X_GK_L <- Case1bis[["X_GK_L"]];U_GK_L <- Case1bis[["U_GK_L"]]
    Xslope_GK_L <- Case1bis[["Xslope_GK_L"]];Uslope_GK_L <- Case1bis[["Uslope_GK_L"]]
    X_GK_T <- Case1bis[["X_GK_T"]];U_GK_T <- Case1bis[["U_GK_T"]]
    Xslope_GK_T <- Case1bis[["Xslope_GK_T"]];Uslope_GK_T <- Case1bis[["Uslope_GK_T"]]
    X_T <- Case1bis[["X_T"]]; U_T <- Case1bis[["U_T"]]; Xslope_T <- Case1bis[["Xslope_T"]]; Uslope_T <- Case1bis[["Uslope_T"]]
    X_L <- Case1bis[["X_L"]]; U_L <- Case1bis[["U_L"]]; Xslope_L <- Case1bis[["Xslope_L"]]; Uslope_L <- Case1bis[["Uslope_L"]]
    X_base <- Case1bis[["X_base"]]; U_base <- Case1bis[["U_base"]]; y.new <- Case1bis[["y.new"]]
    ID.visit <- Case1bis[["ID.visit"]]; offset <- Case1bis[["offset"]];
    B_T_01 <- Case1bis[["B_T_01"]]; Bs_T_01 <- Case1bis[["Bs_T_01"]]; B_L_01 <- Case1bis[["B_L_01"]]; Bs_L_01 <- Case1bis[["Bs_L_01"]]
    B_T_02 <- Case1bis[["B_T_02"]]; Bs_T_02 <- Case1bis[["Bs_T_02"]]; B_L_02 <- Case1bis[["B_L_02"]]; Bs_L_02 <- Case1bis[["Bs_L_02"]]
    B_T_12 <- Case1bis[["B_T_12"]]; Bs_T_12 <- Case1bis[["Bs_T_12"]]; B_L_12 <- Case1bis[["B_L_12"]]; Bs_L_12 <- Case1bis[["Bs_L_12"]]
    if(left_trunc){
      Time_T0 <- Case1bis[["Time_T0"]]; st_T0 <- Case1bis[["st_T0"]]; Bs_T0_01 <- Case1bis[["Bs_T0_01"]]; Bs_T0_02 <- Case1bis[["Bs_T0_02"]]
      X_GK_T0 <- Case1bis[["X_GK_T0"]];U_GK_T0 <- Case1bis[["U_GK_T0"]]
      Xslope_GK_T0 <- Case1bis[["Xslope_GK_T0"]];Uslope_GK_T0 <- Case1bis[["Uslope_GK_T0"]]
    }
    for(i in 1:nbCase1bis){
      if("current value" %in% sharedtype_12 || "current value" %in% sharedtype_01 || "current value" %in% sharedtype_02){
        X_T_i <- X_T[i,];U_T_i <- U_T[i,]
        X_L_i <- X_L[i,];U_L_i <- U_L[i,]
        X_GK_T_i <- X_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_T_i <- U_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        X_GK_L_i <- X_GK_L[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_L_i <- U_GK_L[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          X_GK_T0_i <- X_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_T0_i <- U_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("slope" %in% sharedtype_12 || "slope" %in% sharedtype_01 || "slope" %in% sharedtype_02){
        Xslope_T_i <- Xslope_T[i,];Uslope_T_i <- Uslope_T[i,]
        Xslope_L_i <- Xslope_L[i,];Uslope_L_i <- Uslope_L[i,]
        Xslope_GK_T_i <- Xslope_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_T_i <- Uslope_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Xslope_GK_L_i <- Xslope_GK_L[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_L_i <- Uslope_GK_L[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          Xslope_GK_T0_i <- Xslope_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_T0_i <- Uslope_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }

      if("Weibull" %in% c(hazard_baseline_01,hazard_baseline_02, hazard_baseline_12) ||"Gompertz" %in% c(hazard_baseline_01,hazard_baseline_02, hazard_baseline_12)){
        st_T_i <- st_T[i,]
        st_L_i <- st_L[i,]
        if(left_trunc){
          st_T0_i <- st_T0[i,]
        }
      }
      if("Splines" %in% hazard_baseline_01){
        B_L_i_01 <- B_L_01[i,]
        Bs_L_i_01 <- Bs_L_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          Bs_T0_i_01 <- Bs_T0_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Splines" %in% hazard_baseline_02){
        Bs_L_i_02 <- Bs_L_02[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          Bs_T0_i_02 <- Bs_T0_02[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Splines" %in% hazard_baseline_12){
        B_T_i_12 <- B_T_12[i,]
        Bs_T_i_12 <- Bs_T_12[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Bs_L_i_12 <- Bs_L_12[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
      }
      Z_01_i <- Z_01[i,]
      Z_02_i <- Z_02[i,]
      Z_12_i <- Z_12[i,]
      Time_T_i <- Time_T[i]
      Time_L_i <- Time_L[i]
      if(left_trunc){
        Time_T0_i <- Time_T0[i]
      }
      delta2_i <- delta2[i]

      log_ind_surv <- log_IC_2var_Case1bis(sharedtype,  HB,  Gompertz,  Weibull,
                                           nb_points_integral, alpha_inter_intra,
                                           alpha_y_slope, alpha_z, gamma, beta, beta_slope,
                                           b_y, b_y_slope, wk, sigma_inter, sigma_intra,
                                           delta2_i, Z_01_i,  Z_02_i,  Z_12_i,  X_T_i, U_T_i,
                                           Xslope_T_i,Uslope_T_i, X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                           Uslope_GK_T_i,  X_L_i,  U_L_i,
                                            Xslope_L_i, Uslope_L_i, X_GK_L_i, U_GK_L_i, Xslope_GK_L_i,
                                            Uslope_GK_L_i,
                                            X_GK_T0_i, U_GK_T0_i, Xslope_GK_T0_i, Uslope_GK_T0_i,
                                           Time_T_i, Time_L_i, Time_T0_i, st_T_i, st_L_i, st_T0_i,
                                           B_T_i_01, B_T_i_02, B_T_i_12,
                                           Bs_T_i_01, Bs_T_i_02,  Bs_T_i_12,
                                           B_L_i_01, B_L_i_02,  B_L_i_12,
                                           Bs_L_i_01, Bs_L_i_02, Bs_L_i_12,
                                           Bs_T0_i_01, Bs_T0_i_02, Bs_T0_i_12, left_trunc
      )

      ## Longitudinal part
      X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[i+1]-offset[i])
      U_i <- U_base[offset[i]:(offset[i+1]-1),]
      U_i <- matrix(U_i, nrow = offset[i+1]-offset[i])
      y_i <- y.new[offset[i]:(offset[i+1]-1)]
      ID.visit_i <- ID.visit[offset[i]:(offset[i+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      f_Y_b_sigma <- rep(0,S)
      for(id.visit in 1:length(unique(ID.visit_i))){

        X_base_i.id.visit <- X_base_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1),]
        X_base_i.id.visit <- matrix(X_base_i.id.visit, nrow = offset_ID_i[id.visit+1]-offset_ID_i[id.visit])
        X_base_i.id.visit <- matrix(X_base_i.id.visit[1,], nrow = 1)

        U_i.id.visit <- U_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1),]
        U_i.id.visit <- matrix(U_i.id.visit, nrow = offset_ID_i[id.visit+1]-offset_ID_i[id.visit])
        U_i.id.visit <- matrix(U_i.id.visit[1,],nrow=1)

        y_i.id.visit <- y_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1)]

        if(is.null(nrow(X_base_i.id.visit))){
          stop("There is a something wrong")
        }
        else{
          CV <- (X_base_i.id.visit%*%beta)[1,1] + b_y%*%t(U_i.id.visit)
          n_ij <- length(y_i.id.visit)
          if(n_ij == 1){
            f_Y_b_sigma <- f_Y_b_sigma + log(dnorm(y_i.id.visit,CV,var.inter+var.intra))
          }
          else{
            if(variability_inter_visit && variability_intra_visit){
              if(n_ij == 2){
                f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*sqrt(var.intra*(2*var.inter+var.intra)))) -
                  (1/(2*var.intra*(var.intra+2*var.inter)))*((((rep(y_i.id.visit[1],S)-CV)**2)*(var.intra+var.inter))-2*(var.inter*(rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[2],S)-CV)) +
                                                               (((rep(y_i.id.visit[2],S)-CV)**2)*(var.intra+var.inter)))
              }
              else{
                if(n_ij == 3){
                  f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*var.intra*sqrt((3*var.inter+var.intra)))) -
                    (1/(2*(var.intra**2)*(var.intra+3*var.inter)))*((var.intra*(var.intra+2*var.inter)*((rep(y_i.id.visit[1],S)-CV)**2 + (rep(y_i.id.visit[2],S)-CV)**2 + (rep(y_i.id.visit[3],S)-CV)**2))-
                                                                      2*var.inter*var.intra*((rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[2],S)-CV) + (rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[3],S)-CV) + (rep(y_i.id.visit[2],S)-CV)*(rep(y_i.id.visit[3],S)-CV)))
                }
                else{
                  somme1 <- 0
                  somme2 <- 0
                  for(k in 1:n_ij){
                    somme1 <- somme1 + (rep(y_i.id.visit[k],S)-CV)**2
                    if(k != n_ij){
                      for(l in (k+1):n_ij){
                        somme2 <- somme2 + (rep(y_i.id.visit[k],S)-CV)*(rep(y_i.id.visit[l],S)-CV)
                      }
                    }
                  }
                  f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*sqrt((var.intra**(n_ij-1))*(var.intra+n_ij*var.inter)))) -
                    (1/(2*(var.intra**(n_ij-1))*(var.intra+n_ij*var.inter)))*((var.intra**(n_ij-2))*(var.intra+(n_ij-1)*var.inter)*somme1 -
                                                                                2*var.inter*(var.intra**(n_ij-2))*somme2)
                }
              }
            }
            else{
              stop("Not implemented in this program.")
            }
          }


        }
      }

      log_dens_int <- f_Y_b_sigma + log_ind_surv$SurvTotCase1bis
      Clogexp <- max(log_dens_int) - 500
      log_dens_int <- log_dens_int - Clogexp
      log_dens <- Clogexp + log(sum(exp(log_dens_int))) - log(S)
      ## Left truncation
      if(left_trunc){
        log_dens <- log_dens - log_ind_surv$den
      }
    }
  }

  if(nbCase2 != 0){
    print("Case2 go")
    delta2 <- Case2[["delta2"]]; Z_01 <- Case2[["Z_01"]]; Z_02 <- Case2[["Z_02"]]
    Time_T <- Case2[["Time_T"]]; st_T <- Case2[["st_T"]];
    X_GK_T <- Case2[["X_GK_T"]];U_GK_T <- Case2[["U_GK_T"]]
    Xslope_GK_T <- Case2[["Xslope_GK_T"]];Uslope_GK_T <- Case2[["Uslope_GK_T"]]
    X_T <- Case2[["X_T"]]; U_T <- Case2[["U_T"]]; Xslope_T <- Case2[["Xslope_T"]]; Uslope_T <- Case2[["Uslope_T"]]
    X_base <- Case2[["X_base"]]; U_base <- Case2[["U_base"]]; y.new <- Case2[["y.new"]]
    ID.visit <- Case2[["ID.visit"]]; offset <- Case2[["offset"]];
    B_T_01 <- Case2[["B_T_01"]]; Bs_T_01 <- Case2[["Bs_T_01"]]
    B_T_02 <- Case2[["B_T_02"]]; Bs_T_02 <- Case2[["Bs_T_02"]]
    if(left_trunc){
      Time_T0 <- Case2[["Time_T0"]]; st_T0 <- Case2[["st_T0"]]; Bs_T0_01 <- Case2[["Bs_T0_01"]]; Bs_T0_02 <- Case2[["Bs_T0_02"]]
      X_GK_T0 <- Case2[["X_GK_T0"]];U_GK_T0 <- Case2[["U_GK_T0"]]
      Xslope_GK_T0 <- Case2[["Xslope_GK_T0"]];Uslope_GK_T0 <- Case2[["Uslope_GK_T0"]]
    }

    for(i in 1:nbCase2){

      if("current value" %in% sharedtype_01 || "current value" %in% sharedtype_02){
        X_T_i <- X_T[i,];U_T_i <- U_T[i,]
        X_GK_T_i <- X_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_T_i <- U_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          X_GK_T0_i <- X_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_T0_i <- U_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("slope" %in% sharedtype_01 || "slope" %in% sharedtype_02){
        Xslope_T_i <- Xslope_T[i,];Uslope_T_i <- Uslope_T[i,]
        Xslope_GK_T_i <- Xslope_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_T_i <- Uslope_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          Xslope_GK_T0_i <- Xslope_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_T0_i <- Uslope_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }

      if("Weibull" %in% c(hazard_baseline_01,hazard_baseline_02) ||"Gompertz" %in% c(hazard_baseline_01,hazard_baseline_02)){
        st_T_i <- st_T[i,]
        if(left_trunc){
          st_T0_i <- st_T0[i,]
        }
      }
      if("Splines" %in% hazard_baseline_01){
        Bs_T_i_01 <- Bs_T_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          Bs_T0_i_01 <- Bs_T0_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Splines" %in% hazard_baseline_02){
        B_T_i_02 <- B_T_02[i,]
        Bs_T_i_02 <- Bs_T_02[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        if(left_trunc){
          Bs_T0_i_02 <- Bs_T0_02[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      Z_01_i <- Z_01[i,]
      Z_02_i <- Z_02[i,]
      Time_T_i <- Time_T[i]
      if(left_trunc){
        Time_T0_i <- Time_T0[i]
      }
      delta2_i <- delta2[i]
      log_ind_surv <- log_IC_2var_Case2(sharedtype, HB, Gompertz, Weibull,
                                        nb_points_integral, alpha_inter_intra,
                                        alpha_y_slope, alpha_z, gamma,  beta,  beta_slope,
                                        b_y,  b_y_slope,  wk,  sigma_inter,  sigma_intra,
                                        delta2_i, Z_01_i, Z_02_i, X_T_i,  U_T_i,
                                        Xslope_T_i,  Uslope_T_i,  X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                         Uslope_GK_T_i,
                                         X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i,
                                         Time_T_i,   Time_T0_i, st_T_i,   st_T0_i,
                                         B_T_i_02,
                                         Bs_T_i_01,  Bs_T_i_02,
                                         Bs_T0_i_01,  Bs_T0_i_02,  left_trunc
      )
      ## Longitudinal part
      X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[i+1]-offset[i])
      U_i <- U_base[offset[i]:(offset[i+1]-1),]
      U_i <- matrix(U_i, nrow = offset[i+1]-offset[i])
      y_i <- y.new[offset[i]:(offset[i+1]-1)]
      ID.visit_i <- ID.visit[offset[i]:(offset[i+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      f_Y_b_sigma <- rep(0,S)
      for(id.visit in 1:length(unique(ID.visit_i))){

        X_base_i.id.visit <- X_base_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1),]
        X_base_i.id.visit <- matrix(X_base_i.id.visit, nrow = offset_ID_i[id.visit+1]-offset_ID_i[id.visit])
        X_base_i.id.visit <- matrix(X_base_i.id.visit[1,], nrow = 1)

        U_i.id.visit <- U_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1),]
        U_i.id.visit <- matrix(U_i.id.visit, nrow = offset_ID_i[id.visit+1]-offset_ID_i[id.visit])
        U_i.id.visit <- matrix(U_i.id.visit[1,],nrow=1)

        y_i.id.visit <- y_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1)]

        if(is.null(nrow(X_base_i.id.visit))){
          stop("There is a something wrong")
        }
        else{
          CV <- (X_base_i.id.visit%*%beta)[1,1] + b_y%*%t(U_i.id.visit)
          n_ij <- length(y_i.id.visit)
          if(n_ij == 1){
            f_Y_b_sigma <- f_Y_b_sigma + log(dnorm(y_i.id.visit,CV,var.inter+var.intra))
          }
          else{
            if(variability_inter_visit && variability_intra_visit){
              if(n_ij == 2){
                f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*sqrt(var.intra*(2*var.inter+var.intra)))) -
                  (1/(2*var.intra*(var.intra+2*var.inter)))*((((rep(y_i.id.visit[1],S)-CV)**2)*(var.intra+var.inter))-2*(var.inter*(rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[2],S)-CV)) +
                                                               (((rep(y_i.id.visit[2],S)-CV)**2)*(var.intra+var.inter)))
              }
              else{
                if(n_ij == 3){
                  f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*var.intra*sqrt((3*var.inter+var.intra)))) -
                    (1/(2*(var.intra**2)*(var.intra+3*var.inter)))*((var.intra*(var.intra+2*var.inter)*((rep(y_i.id.visit[1],S)-CV)**2 + (rep(y_i.id.visit[2],S)-CV)**2 + (rep(y_i.id.visit[3],S)-CV)**2))-
                                                                      2*var.inter*var.intra*((rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[2],S)-CV) + (rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[3],S)-CV) + (rep(y_i.id.visit[2],S)-CV)*(rep(y_i.id.visit[3],S)-CV)))
                }
                else{
                  somme1 <- 0
                  somme2 <- 0
                  for(k in 1:n_ij){
                    somme1 <- somme1 + (rep(y_i.id.visit[k],S)-CV)**2
                    if(k != n_ij){
                      for(l in (k+1):n_ij){
                        somme2 <- somme2 + (rep(y_i.id.visit[k],S)-CV)*(rep(y_i.id.visit[l],S)-CV)
                      }
                    }
                  }
                  f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*sqrt((var.intra**(n_ij-1))*(var.intra+n_ij*var.inter)))) -
                    (1/(2*(var.intra**(n_ij-1))*(var.intra+n_ij*var.inter)))*((var.intra**(n_ij-2))*(var.intra+(n_ij-1)*var.inter)*somme1 -
                                                                                2*var.inter*(var.intra**(n_ij-2))*somme2)
                }
              }
            }
            else{
              stop("Not implemented in this program.")
            }
          }


        }
      }

      log_dens_int <- f_Y_b_sigma + log_ind_surv$SurvTotCase2
      Clogexp <- max(log_dens_int) - 500
      log_dens_int <- log_dens_int - Clogexp
      log_dens <- Clogexp + log(sum(exp(log_dens_int))) - log(S)
      ## Left truncation
      if(left_trunc){
        log_dens <- log_dens - log_ind_surv$den
      }
      ll_glob <- ll_glob + log_dens

    }
  }

  if(nbCase3 != 0){
    print("Case3 go")
    delta2 <- Case3[["delta2"]]; Z_01 <- Case3[["Z_01"]]; Z_02 <- Case3[["Z_02"]]; Z_12 <- Case3[["Z_12"]]
    Time_T <- Case3[["Time_T"]]; st_T <- Case3[["st_T"]]; Time_L <- Case3[["Time_L"]]
    X_GK_T <- Case3[["X_GK_T"]];U_GK_T <- Case3[["U_GK_T"]]
    Xslope_GK_T <- Case3[["Xslope_GK_T"]];Uslope_GK_T <- Case3[["Uslope_GK_T"]]
    X_T <- Case3[["X_T"]]; U_T <- Case3[["U_T"]]; Xslope_T <- Case3[["Xslope_T"]]; Uslope_T <- Case3[["Uslope_T"]]
    X_base <- Case3[["X_base"]]; U_base <- Case3[["U_base"]]; y.new <- Case3[["y.new"]]
    ID.visit <- Case3[["ID.visit"]]; offset <- Case3[["offset"]];
    B_T_01 <- Case3[["B_T_01"]]; Bs_T_01 <- Case3[["Bs_T_01"]]
    B_T_02 <- Case3[["B_T_02"]]; Bs_T_02 <- Case3[["Bs_T_02"]]
    B_T_12 <- Case3[["B_T_12"]]; Bs_T_12 <- Case3[["Bs_T_12"]]
    if(left_trunc){
      Time_T0 <- Case3[["Time_T0"]]; st_T0 <- Case3[["st_T0"]]; Bs_T0_01 <- Case3[["Bs_T0_01"]]; Bs_T0_02 <- Case3[["Bs_T0_02"]]
      X_GK_T0 <- Case3[["X_GK_T0"]];U_GK_T0 <- Case3[["U_GK_T0"]]
      Xslope_GK_T0 <- Case3[["Xslope_GK_T0"]];Uslope_GK_T0 <- Case3[["Uslope_GK_T0"]]
    }
    ### Partie dépendante de l'intégrale :
    X_GK_L_T <- Case3[["X_GK_L_T"]]; U_GK_L_T <- Case3[["U_GK_L_T"]]
    Xslope_GK_L_T <- Case3[["Xslope_GK_L_T"]];Uslope_GK_L_T<- Case3[["Uslope_GK_L_T"]]
    Bs_L_T_01 <- Case3[["Bs_L_T_01"]]; Time_L_T <- Case3[["Time_L_T"]]; st_L_T <- Case3[["st_L_T"]] #Time_L_R = R-L # st_L_R = (xk+1)/2*(R-L)+L
    ### Integrale de 0 à u (u dans L_i-T_i)
    st_0_LT <- Case3[["st_0_LT"]]; X_0_LT <- Case3[["X_0_LT"]]; U_0_LT <- Case3[["U_0_LT"]];
    Xslope_0_LT <- Case3[["Xslope_0_LT"]]; Uslope_0_LT <- Case3[["Uslope_0_LT"]];
    Bs_0_LT_01 <- Case3[["Bs_0_LT_01"]]; Bs_0_LT_02 <- Case3[["Bs_0_LT_02"]]; Bs_0_LT_12 <- Case3[["Bs_0_LT_12"]];

    for(i in 1:nbCase3){
      if("current value" %in% sharedtype_01 || "current value" %in% sharedtype_02 || "current value" %in% sharedtype_12){
        X_T_i <- X_T[i,];U_T_i <- U_T[i,]
        X_GK_T_i <- X_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_T_i <- U_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        X_GK_L_T_i <- X_GK_L_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_L_T_i <- U_GK_L_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        #####()
        X_GK_0_LT_i <- X_0_LT[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        U_GK_0_LT_i <- U_0_LT[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        if(left_trunc){
          X_GK_T0_i <- X_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];U_GK_T0_i <- U_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("slope" %in% sharedtype_01 || "slope" %in% sharedtype_02 || "slope" %in% sharedtype_12){
        Xslope_T_i <- Xslope_T[i,];Uslope_T_i <- Uslope_T[i,]
        Xslope_GK_T_i <- Xslope_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_T_i <- Uslope_GK_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Xslope_GK_L_T_i <- Xslope_GK_L_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_L_T_i <- Uslope_GK_L_T[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        #####()
        Xslope_0_LT_i <- Xslope_0_LT[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        Uslope_0_LT_i <- Uslope_0_LT[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        if(left_trunc){
          Xslope_GK_T0_i <- Xslope_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),];Uslope_GK_T0_i <- Uslope_GK_T0[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Weibull" %in% c(hazard_baseline_01,hazard_baseline_02, hazard_baseline_12) ||"Gompertz" %in% c(hazard_baseline_01,hazard_baseline_02, hazard_baseline_12)){
        st_T_i <- st_T[i,]
        st_0_LT_i <- st_0_LT[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        st_L_T_i <- st_L_T[i,]
        if(left_trunc){
          st_T0_i <- st_T0[i,]
        }
      }
      if("Splines" %in% hazard_baseline_01){
        Bs_T_i_01 <- Bs_T_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Bs_L_T_i_01 <- Bs_L_T_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Bs_0_LT_i_01 <- Bs_0_LT_01[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        if(left_trunc){
          Bs_T0_i_01 <- Bs_T0_01[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Splines" %in% hazard_baseline_02){
        B_T_i_02 <- B_T_02[i,]
        Bs_T_i_02 <- Bs_T_02[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Bs_0_LT_i_02 <- Bs_0_LT_02[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
        if(left_trunc){
          Bs_T0_i_02 <- Bs_T0_02[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        }
      }
      if("Splines" %in% hazard_baseline_12){
        B_T_i_12 <- B_T_12[i,]
        Bs_T_i_12 <- Bs_T_12[(nb_pointsGK*(i-1)+1):(nb_pointsGK*i),]
        Bs_0_LT_i_12 <- Bs_0_LT_12[(nb_pointsGK*nb_pointsGK*(i-1)+1):(nb_pointsGK*nb_pointsGK*i),]
      }
      Z_01_i <- Z_01[i,]
      Z_02_i <- Z_02[i,]
      Z_12_i <- Z_12[i,]
      Time_T_i <- Time_T[i]
      ck <- ((sk_GK+1)/4)*Time_L_T[i]+Time_L[i]/2
      Time_L_T_i <- Time_L_T[i]
      if(left_trunc){
        Time_T0_i <- Time_T0[i]
      }
      delta2_i <- delta2[i]

      log_ind_surv <- log_IC_2var_Case3( sharedtype,  HB,  Gompertz,  Weibull,
                                         nb_points_integral,  alpha_inter_intra,
                                         alpha_y_slope,  alpha_z,  gamma,  beta,  beta_slope,
                                         b_y,  b_y_slope,  wk,  rep_wk,  sigma_inter,  sigma_intra,
                                        delta2_i,  Z_01_i,  Z_02_i,  Z_12_i,  X_T_i,  U_T_i,
                                         Xslope_T_i,  Uslope_T_i,  X_GK_T_i,  U_GK_T_i,  Xslope_GK_T_i,
                                         Uslope_GK_T_i,  X_GK_L_T_i,  U_GK_L_T_i,  Xslope_GK_L_T_i,  Uslope_GK_L_T_i,
                                         X_GK_0_LT_i,  U_GK_0_LT_i,  Xslope_GK_0_LT_i,  Uslope_GK_0_LT_i,
                                         X_GK_T0_i,  U_GK_T0_i,  Xslope_GK_T0_i,  Uslope_GK_T0_i,
                                         Time_T_i,  Time_L_T_i,  Time_T0_i, st_T_i,  st_0_LT_i,  st_L_T_i,  st_T0_i,
                                         ck,
                                         B_T_i_01,  B_T_i_02,  B_T_i_12,
                                         Bs_T_i_01,  Bs_T_i_02,  Bs_T_i_12,
                                         Bs_0_LT_i_01,  Bs_0_LT_i_02,  Bs_0_LT_i_12,
                                         Bs_L_T_i_01,   Bs_L_T_i_02,  Bs_L_T_i_12,
                                         Bs_T0_i_01,  Bs_T0_i_02, left_trunc
      )

      ## Longitudinal part
      X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
      X_base_i <- matrix(X_base_i, nrow = offset[i+1]-offset[i])
      U_i <- U_base[offset[i]:(offset[i+1]-1),]
      U_i <- matrix(U_i, nrow = offset[i+1]-offset[i])
      y_i <- y.new[offset[i]:(offset[i+1]-1)]
      ID.visit_i <- ID.visit[offset[i]:(offset[i+1]-1)]
      offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
      f_Y_b_sigma <- rep(0,S)
      for(id.visit in 1:length(unique(ID.visit_i))){

        X_base_i.id.visit <- X_base_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1),]
        X_base_i.id.visit <- matrix(X_base_i.id.visit, nrow = offset_ID_i[id.visit+1]-offset_ID_i[id.visit])
        X_base_i.id.visit <- matrix(X_base_i.id.visit[1,], nrow = 1)

        U_i.id.visit <- U_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1),]
        U_i.id.visit <- matrix(U_i.id.visit, nrow = offset_ID_i[id.visit+1]-offset_ID_i[id.visit])
        U_i.id.visit <- matrix(U_i.id.visit[1,],nrow=1)

        y_i.id.visit <- y_i[offset_ID_i[id.visit]:(offset_ID_i[id.visit+1]-1)]

        if(is.null(nrow(X_base_i.id.visit))){
          stop("There is a something wrong")
        }
        else{
          CV <- (X_base_i.id.visit%*%beta)[1,1] + b_y%*%t(U_i.id.visit)
          n_ij <- length(y_i.id.visit)
          if(n_ij == 1){
            f_Y_b_sigma <- f_Y_b_sigma + log(dnorm(y_i.id.visit,CV,var.inter+var.intra))
          }
          else{
            if(variability_inter_visit && variability_intra_visit){
              if(n_ij == 2){
                f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*sqrt(var.intra*(2*var.inter+var.intra)))) -
                  (1/(2*var.intra*(var.intra+2*var.inter)))*((((rep(y_i.id.visit[1],S)-CV)**2)*(var.intra+var.inter))-2*(var.inter*(rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[2],S)-CV)) +
                                                               (((rep(y_i.id.visit[2],S)-CV)**2)*(var.intra+var.inter)))
              }
              else{
                if(n_ij == 3){
                  f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*var.intra*sqrt((3*var.inter+var.intra)))) -
                    (1/(2*(var.intra**2)*(var.intra+3*var.inter)))*((var.intra*(var.intra+2*var.inter)*((rep(y_i.id.visit[1],S)-CV)**2 + (rep(y_i.id.visit[2],S)-CV)**2 + (rep(y_i.id.visit[3],S)-CV)**2))-
                                                                      2*var.inter*var.intra*((rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[2],S)-CV) + (rep(y_i.id.visit[1],S)-CV)*(rep(y_i.id.visit[3],S)-CV) + (rep(y_i.id.visit[2],S)-CV)*(rep(y_i.id.visit[3],S)-CV)))
                }
                else{
                  somme1 <- 0
                  somme2 <- 0
                  for(k in 1:n_ij){
                    somme1 <- somme1 + (rep(y_i.id.visit[k],S)-CV)**2
                    if(k != n_ij){
                      for(l in (k+1):n_ij){
                        somme2 <- somme2 + (rep(y_i.id.visit[k],S)-CV)*(rep(y_i.id.visit[l],S)-CV)
                      }
                    }
                  }
                  f_Y_b_sigma <- f_Y_b_sigma + log(1/(((2*pi)**(n_ij/2))*sqrt((var.intra**(n_ij-1))*(var.intra+n_ij*var.inter)))) -
                    (1/(2*(var.intra**(n_ij-1))*(var.intra+n_ij*var.inter)))*((var.intra**(n_ij-2))*(var.intra+(n_ij-1)*var.inter)*somme1 -
                                                                                2*var.inter*(var.intra**(n_ij-2))*somme2)
                }
              }
            }
            else{
              stop("Not implemented in this program.")
            }
          }


        }
      }

      log_dens_int <- f_Y_b_sigma + log_ind_surv$SurvTotCase3
      Clogexp <- max(log_dens_int) - 500
      log_dens_int <- log_dens_int - Clogexp
      log_dens <- Clogexp + log(sum(exp(log_dens_int))) - log(S)
      ## Left truncation
      if(left_trunc){
        log_dens <- log_dens - log_ind_surv$den
      }
      ll_glob <- ll_glob + log_dens


    }
  }

  if(is.na(ll_glob)){
    print(param)
    ll_glob <- -1E09
  }
  print(ll_glob)
  ll_glob

}
