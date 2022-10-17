####params & initial conditions####
#### by Suyi ####
#### 2022.09.24####

#load Shanghai CM data
load("C:/Users/lenovo/Desktop/data/SH_CM.RData")
#age group & population settings
num_of_age_group <- 10
pop <- c(36011.6,129846.5,75153.8,166697.1,181417,147494.3,143708.97,82008.3,41638.46,27417.96)


#init settings
generate_init_condi <- function(Di = 1/5,
                                Dp = 1/(3.2 - 2),
                                De = 1/2,
                                Dq = 1/2,
                                Dh = 1/14,
                                N = as.numeric(pop),
                                alpha = 0.6,
                                alpha_2 = 1,
                                r = 0.5,
                                p_c = 0.05,
                                p_v2=c(0,0.9588,1,1,1,1,1,0.9301,	0.7445,	0.3617),
                                p_v3 =c(0,0,0.0006,	0.5450,	0.6662,	0.7780,	0.8569,	0.8818,	0.6746,	0.2129),
                                Ve2=0.04371984795,
                                Ve3=0.2640376479,
                                beta=c(0.04966975, 0.04966975, 0.04966975, 0.15695642, 0.19867902, 0.19867902, 0.22450729, 0.24338180, 0.26225630, 0.26225630)
) {
  
  stopifnot(r>=0 & r<=1 & Di>=0 & Dp>=0 & De>=0 & Dp>=0 & alpha>=0 & alpha<=1 & Dh>=0 & alpha_2>=0 & alpha_2<=1 & p_c>0 & all(p_v2>=0) & all(p_v3>=0) & Ve2>0 &Ve3>0)
  
  #Real data from CDC
  realData_all <- readxl::read_excel("C:/Users/lenovo/Desktop/data/obs(Sanya)0914.xlsx",sheet = 1) 
  realData <- realData_all
  
  # % of beta
  beta_percentage <- beta/sum(beta)
  
  #cases from Aug 1 to Sep 3
  daily_new_case <- realData[1:43, 8]
  daily_new_case_all <- realData$total
  
  #CM
  ContactMatrix_t <- get("ContactMatrix", envir = globalenv())
  COntactMatrix_list <- list(matrix(ContactMatrix, ncol = 10, nrow = 10),
                             matrix(ContactMatrix3, ncol = 10, nrow = 10),
                             matrix(ContactMatrix6, ncol = 10, nrow = 10),
                             matrix(ContactMatrix8, ncol = 10, nrow = 10)) 
  
  ##Initial condition
  aug1_idx = 1
  ##
  E0_all <- sum(daily_new_case_all[(aug1_idx+round(1/Dp)):(aug1_idx+round(1/Dp)+round(1/De)-1)]) / r 
  E0 <- E0_all * beta_percentage
  # E0 <- (40 + 23 + 47) / r                                              
  P0_all <- 10
  P0 <- P0_all * beta_percentage
  # P0 <- (41 + 34) / r
  #I0_all <- sum(daily_new_case_all[(aug1_idx-round(1/Di)):(aug1_idx-1)])   
  #I0 <- round(I0_all * beta_percentage)
  # I0 <- 11 + 13 + 10       
  I0=c(rep(0,5),1,rep(0,4))
  A0 <- I0 * (1 - r) / r
  S0 <- pop - E0 - P0 - I0 - A0
  
  #S0=c(36011.6,129846.5,75153.8,166697.1,181417,147494.3,143708.97,82008.3,41638.46,27417.96)
  #E0=rep(0,10)
  #P0=rep(0,10)
  #A0=rep(0,10)

  R0=rep(0,10) 
  Eq0=rep(0,10)
  Pq0=rep(0,10)
  Aq0=rep(0,10)
  Iq0=rep(0,10)
  
  init_states <- round(c(S = S0, E = E0, P = P0, I = I0, A = A0, R = R0, Eq = Eq0, Pq = Pq0, Aq = Aq0, Iq = Iq0), 0)
  
  transform_var=function(pars) {
    ##
    PAAT1 <- pars[1]
    PAAT2 <- pars[2]
    PAAT3 <- pars[3]
    PAAT4 <- pars[4]
    PAAT_vec <- c(PAAT1,PAAT2,PAAT3,PAAT4)
    return(PAAT_vec)
  }
  
  return(list(Di=Di,
              Dp=Dp,
              De=De,
              Dq=Dq,
              alpha=alpha,
              alpha_2 = alpha_2,
              Dh=Dh,
              N = N,
              r = r,
              p_c = p_c,
              p_v2 = p_v2,
              p_v3 = p_v3,
              Ve2 = Ve2,
              Ve3 = Ve3,
              daily_new_case = daily_new_case, 
              daily_new_case_all = daily_new_case_all, 
              init_states = init_states,
              beta = beta,
              beta_percentage = beta_percentage,
              days_to_fit=1:43,
              stage_intervals=list(
                c(start=1, end=3),
                c(start=4, end=9),
                c(start=10, end=15),
                c(start=16, end=43)
              ),
              var_trans_fun=transform_var,
              ContactMatrix_list = COntactMatrix_list,
              par_lower = c(PAAT1 = 1e-5, PAAT2 = 1e-5, PAAT3 = 1e-5, PAAT4 = 1e-5),
              par_upper = c(PAAT1 = 1, PAAT2 = 1, PAAT3 = 1, PAAT4 = 1)))
}

# get_init_sets_list is an alias of generate_init_condi in order not to break existing code
init_settings = generate_init_condi(Di = 1/5,
                                    Dp = 1/(3.2 - 2),
                                    De = 1/2,
                                    Dq = 1/2,
                                    Dh = 1/14,
                                    N = as.numeric(pop),
                                    alpha = 0.6,
                                    alpha_2 = 1,
                                    r = 0.5,
                                    p_c = 0.05,
                                    p_v2=c(0,0.9588,1,1,1,1,1,0.9301,	0.7445,	0.3617),
                                    p_v3 =c(0,0,0.0006,	0.5450,	0.6662,	0.7780,	0.8569,	0.8818,	0.6746,	0.2129),
                                    Ve2=0.04371984795,
                                    Ve3=0.2640376479)


PAAT_mean1 <- 0.5
PAAT_sd1 <- 0.5
PAAT_mean2 <-0.5
PAAT_sd2 <- 0.5
PAAT_mean3 <- 0.5
PAAT_sd3 <- 0.5
PAAT_mean4 <- 0.5
PAAT_sd4 <-0.5
