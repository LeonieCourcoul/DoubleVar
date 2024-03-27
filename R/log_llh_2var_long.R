#' Title
#'
#' @param param
#' @param nb.beta
#' @param Zq
#' @param nb.e.a
#' @param S
#' @param variability_inter_visit
#' @param variability_intra_visit
#' @param correlated_re
#' @param X_base
#' @param U_base
#' @param y.new
#' @param ID.visit
#' @param Ind
#' @param offset
#'
#' @return
#' @export
#'
#' @examples
log_llh_2var_long <- function(param,nb.beta, Zq,
                              nb.e.a, S,
                              variability_inter_visit, variability_intra_visit,
                              correlated_re, X_base, U_base, y.new, ID.visit, Ind, offset
){#31


  #browser()
  #Manage parameter
  curseur <- 1
  ## Marker
  ### Fiexd effects
  beta <- param[curseur:(curseur+nb.beta-1)]
  curseur <- curseur+nb.beta
  if(variability_inter_visit){
    mu.inter <- param[curseur]
    curseur <- curseur + 1
  }
  if(variability_intra_visit){
    mu.intra <- param[curseur]
    curseur <- curseur + 1
  }
  ### Cholesky matrix for random effects
  if(variability_inter_visit && variability_intra_visit){
    if(correlated_re){

      C1 <- matrix(rep(0,(nb.e.a+2)**2),nrow=nb.e.a+2,ncol=nb.e.a+2)
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
        C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
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
        Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
        Cholesky <- as.matrix(Cholesky)
      }
    }
    else{
      if(variability_intra_visit){
        if(correlated_re){
          C1 <- matrix(rep(0,(nb.e.a+1)**2),nrow=nb.e.a+1,ncol=nb.e.a+1)
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
          Cholesky <- rbind(cbind(C1,C3),cbind(C4,C2))
          Cholesky <- as.matrix(Cholesky)
        }
      }
    }
    if(!variability_inter_visit && !variability_intra_visit){
      C1 <- matrix(rep(0,(length(param)-curseur)**2),nrow=length(param)-curseur,ncol=length(param)-curseur)
      C1[lower.tri(C1, diag=T)] <- param[curseur:length(param)]
      Cholesky <- C1
    }

  }
  # Manage random effects
  random.effects <- Zq%*%t(Cholesky)
  b_y <- random.effects[,1:nb.e.a]
  b_y <- matrix(b_y, ncol = nb.e.a)
  if(variability_inter_visit){
    b_inter <- random.effects[,nb.e.a+1]
    sigma.inter <- exp(mu.inter + b_inter)
    var.inter <- sigma.inter**2
  }
  if(variability_intra_visit){
    b_intra <- random.effects[,nb.e.a+2]
    sigma.intra <- exp(mu.intra + b_intra)
    var.intra <- sigma.intra**2
  }
  ll_glob <- 0
  for(i in 1:Ind){

    # Manage random effects
    random.effects <- Zq%*%t(Cholesky)
    b_y <- random.effects[,1:nb.e.a]
    b_y <- matrix(b_y, ncol = nb.e.a)
    if(variability_inter_visit){
      b_inter <- random.effects[,nb.e.a+1]
      sigma_inter <- exp(mu.inter + b_inter)
      var.inter <- sigma_inter**2
    }
    else{
      sigma_inter <- rep(sigma.epsilon.inter,S)
      var.inter <- rep(sigma.epsilon.inter**2,S)
    }
    if(variability_intra_visit){
      b_intra <- random.effects[,ncol(random.effects)]
      sigma_intra <- exp(mu.intra + b_intra)
      var.intra <- sigma_intra**2
    }
    else{
      sigma_intra <- rep(sigma.epsilon.intra,S)
      var.intra <- rep(sigma.epsilon.intra**2,S)
    }
    sigma_inter_intra <- list(sigma_inter, sigma_intra, var.inter+var.intra, var.inter, var.intra, var.intra*(2*var.inter+var.intra))

    X_base_i <- X_base[offset[i]:(offset[i+1]-1),]
    X_base_i <- matrix(X_base_i, nrow = offset[i+1]-offset[i])
    X_base_i <- unique(X_base_i)
    U_i <- U_base[offset[i]:(offset[i+1]-1),]
    U_i <- matrix(U_i, nrow = offset[i+1]-offset[i])
    U_base_i <- unique(U_i)
    y_i <- y.new[offset[i]:(offset[i+1]-1)]
    ID.visit_i <- ID.visit[offset[i]:(offset[i+1]-1)]
    offset_ID_i <- as.vector(c(1, 1 + cumsum(tapply(ID.visit_i, ID.visit_i, length))))
    len_visit_i <- length(unique(ID.visit_i))
    log_ind <- log_llh_2var_longInd( S,  sigma_inter_intra, X_base_i, U_base_i,
                         y_i, len_visit_i,offset_ID_i, beta,
                         b_y)

    log_dens_int <- log_ind
    Clogexp <- max(log_dens_int) - 500
    log_dens_int <- log_dens_int - Clogexp
    log_dens <- Clogexp + log(sum(exp(log_dens_int))) - log(S)
    ll_glob <- ll_glob + log_dens

  }
  print(ll_glob)
  ll_glob

}
