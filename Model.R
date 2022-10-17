####SEIRpred
## Deterministic SEIR model
#### By Suyi 2022/09/24

#Model
SEIRpred <- function(pars, 
                     init_settings) {
  tmp_ret=init_settings$var_trans_fun(pars)
  beta_vec=init_settings$beta
  PAAT_vec=tmp_ret
  stage_intervals=init_settings$stage_intervals
  n_stage=length(stage_intervals)
  
  ##
  Di <- init_settings$Di
  Dp <- init_settings$Dp
  De <- init_settings$De
  Dq <- init_settings$Dq
  alpha <- init_settings$alpha
  alpha_2 <- init_settings$alpha_2
  Dh <- init_settings$Dh
  N <- init_settings$N
  ContactMatrix_list <- init_settings$ContactMatrix_list
  r = init_settings$r
  p_c = init_settings$p_c
  p_v2 = init_settings$p_v2
  p_v3 = init_settings$p_v3
  Ve2 = init_settings$Ve2
  Ve3 = init_settings$Ve3
  init_states <- init_settings$init_states
  days_to_fit <- init_settings$days_to_fit
  beta_percentate <- init_settings$beta_percentate
  
  ## ODE function based on deterministic SEIR model
  update_func <- function(stage_pars, states_old) {
    ## stage pars
    beta <- stage_pars[[1]]
    PAAT <- stage_pars[[2]]
    ContactMatrix_use <- stage_pars[[3]]
    
    ## old states number: c(S, E, P, I, A, R, Eq, Pq, Aq, Iq)
    S <- states_old[(1:num_of_age_group) + num_of_age_group * 0]
    E <- states_old[(1:num_of_age_group) + num_of_age_group * 1]
    P <- states_old[(1:num_of_age_group) + num_of_age_group * 2]
    I <- states_old[(1:num_of_age_group) + num_of_age_group * 3]
    A <- states_old[(1:num_of_age_group) + num_of_age_group * 4]
    R <- states_old[(1:num_of_age_group) + num_of_age_group * 5]
    Eq <- states_old[(1:num_of_age_group) + num_of_age_group * 6]
    Pq <- states_old[(1:num_of_age_group) + num_of_age_group * 7]
    Aq <- states_old[(1:num_of_age_group) + num_of_age_group * 8]
    Iq <- states_old[(1:num_of_age_group) + num_of_age_group * 9]
    
    CM <- ContactMatrix_use %*% (alpha * A + alpha * P + I)
    ## new values
    S_new <- S + -beta * CM * S * (1 - (p_v2-p_v3) * Ve2 - p_v3 * Ve3) / N
    E_new <- E + (1 - p_c) * beta * CM * S * (1 - (p_v2-p_v3) * Ve2 - p_v3 * Ve3) / N - E * De
    
    P_new <- P + E * De - (1 - PAAT) * P * Dp - PAAT * P * alpha_2
    A_new <- A + (1 - r) * (1 - PAAT) * P * Dp - (1 - PAAT) * A * Di - PAAT * A * alpha_2 
    I_new <- I + r * (1 - PAAT) * P * Dp - Dq * I
    
    Eq_new <- Eq + p_c * beta * CM * S * (1 - (p_v2-p_v3) * Ve2 - p_v3 * Ve3) / N - Eq * De
    Pq_new <- Pq + Eq * De + PAAT * P * alpha_2 - Pq * Dp 
    Aq_new <- Aq + (1 - r) * Pq * Dp + PAAT * A * alpha_2 - Aq * Di
    Iq_new <- Iq + r * Pq * Dp + Dq * I - Dh * I
    
    R_new <- R + (1 - PAAT) * A * Di + Aq * Di + Dh * I
    
    Onset_expect <- PAAT*P*Dp
    ##
    return(c(S_new, E_new, P_new, I_new, A_new, Eq_new, Pq_new, Iq_new, Aq_new, R_new, Onset_expect))
  }
  ## matrix for results
  states_mat <- matrix(0, length(days_to_fit), length(init_states) + 11)
  states_mat[, 1] <- days_to_fit
  colnames(states_mat) <- c("time", 
                            "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
                            "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", 
                            "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", 
                            "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9", "I10", 
                            "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", 
                            "Eq1", "Eq2", "Eq3", "Eq4", "Eq5", "Eq6", "Eq7", "Eq8", "Eq9", "Eq10", 
                            "Pq1", "Pq2", "Pq3", "Pq4", "Pq5", "Pq6", "Pq7", "Pq8", "Pq9", "Pq10", 
                            "Iq1", "Iq2", "Iq3", "Iq4", "Iq5", "Iq6", "Iq7", "Iq8", "Iq9", "Iq10", 
                            "Aq1", "Aq2", "Aq3", "Aq4", "Aq5", "Aq6", "Aq7", "Aq8", "Aq9", "Aq10", 
                            "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10",
                            "Onset_expect1", "Onset_expect2", "Onset_expect3", "Onset_expect4", "Onset_expect5", "Onset_expect6", "Onset_expect7", "Onset_expect8", "Onset_expect9", "Onset_expect10")
  
  myold_states <- init_states
  
  for (i_stage in 1:n_stage) {
    #print(222)
    stage_pars_setings <- list(
      beta = beta_vec,
      PAAT = PAAT_vec[i_stage],
      ContactMatrix = ContactMatrix_list[[i_stage]]
    )
    for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
      myold_states <- states_mat[d, -1]
    }
  }
  return(states_mat)
}
