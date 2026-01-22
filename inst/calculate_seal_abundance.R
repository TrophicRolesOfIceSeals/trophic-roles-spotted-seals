## calculate variance given a lognormal CI
#provide lower CI, higher CI, estimate, and the nominal coverage (e.g., 0.95 for 95%)

get_var_from_lognormal_CI = function(low,high,est,coverage){
  C = high/est
  z = qnorm(1-0.5*(1-coverage)) #for 2 sided interval
  est^2 * (exp((log(C)/z)^2)-1)  #variance
}

#wBS 2012
var_wBS_2012 = get_var_from_lognormal_CI(151092,338052,226003,0.95)

#wBS 2013
var_wBS_2013 = get_var_from_lognormal_CI(62681,130670,90501,0.95)

#eBS 2012
var_eBS_2012 = get_var_from_lognormal_CI(233280,360417,289962,0.95)

#eBS 2013
var_eBS_2013 = get_var_from_lognormal_CI(113638,157268,133685,0.95)


N_2012 = 226003 + 289962
N_2013 = 90501 + 133685
var_2012 = var_wBS_2012 + var_eBS_2012
var_2013 = var_wBS_2013 + var_eBS_2013

N_wt = (var_2013*N_2012+var_2012*N_2013)/(var_2012+var_2013)
var_N = 1/(var_2012+var_2013)^2 * (var_2013^2*var_2012 + var_2012^2*var_2013)
var_log_N = 1/N_wt^2 * var_N
mu_log_N = log(N_wt)-0.5*var_log_N

#check
mean(exp(rnorm(1000000,mu_log_N,sqrt(var_log_N))))

mu_log_N   #12.45542
var_log_N  #0.00556712

#plot histogram
N = data.frame("N"=exp(rnorm(1000000,mu_log_N,sqrt(var_log_N))))

library(ggplot2)
png("N_hist.png",width = 8, height = 8, units='in', res = 600)
ggplot(N)+geom_density(aes(x=N))+ylab('Probability density')+
  theme(text = element_text(size=14))+xlab('Abundance')
#ggplot(N)+geom_histogram(aes(x=N,y=..density..))+ylab('Probability density')+
#  theme(text = element_text(size=14))+xlab('Abundance')
dev.off()

