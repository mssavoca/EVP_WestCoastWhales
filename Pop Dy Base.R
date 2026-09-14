#EVP: h(t) = h0(t) + h_anth(t)
#h(t) = total mortality at time t
#h0(t) = natural mortality 
#h_anth(t) = anthropogenic mortality 
#h_anth(t) = h_bycatch(t) + h_strike(t) + h_climate_t

bycatch_min = 3
bycatch_max = 53
strike_min = 8
strike_max = 53
climate_min = -500
climate_max = 500

N0_mean = 4973 #estimate from Calambokidis and Barlow 2020
N0_min = 4776 #lower and upper 20th percentiles 
N0_max = 5178

M_avg = 1-0.982 #from Gabriele et al. 2022 - non-calf survival before MHW (in AK, will want to replace )
M_sd = 0.003 #actually the SE but they're using it in the error bounds

#in Barlow et al. 2011, used by Curtis et al. 2022: 0.96 annual adult survival, 11% birth rate, 0.85 semi-annual survival rate for new calves   
birth_rate = 0.11
birth_surv = 0.85

n_years = 2050-2018 + 1
n_sims = 10



N = as.data.frame(matrix(nrow = n_years, ncol = n_sims))

set.seed(123)

for(s in 1:ncol(N)){
  N[1,s] = round(runif(n = 1, min = N0_min, max = N0_max)) #initial population size at year 1 (2018)
  
  for(y in 2:nrow(N)){
    #randomly select mortality 
    bycatch = round(runif(n = 1, min = bycatch_min, max = bycatch_max))
    strike = round(runif(n = 1, min = strike_min, max = strike_max))
    #climate = round(runif(n = 1, min = climate_min, max = climate_max))
    
    #convert mortality numbers to rates
    h_by = bycatch/N[y-1,s]
    h_strike = strike/N[y-1,s]
    #h_clim = climate/N[y-1,s]
    
    #natural mortality rate 
    M = rnorm(n = 1, mean = M_avg, sd = M_sd)
    
    #total mortality
    Z = M + h_by + h_strike #+ h_clim
    
    #survival
    surv_rate = exp(-Z)
    #the potential adjustment for climate is so large, so it can lead survival rate to be greater than 1
    #surv_rate = ifelse(surv_rate > 1, 1, surv_rate)
    
    #number that survive
    N_surv = round(rbinom(n = 1, size = N[y-1,s], prob = surv_rate))
    
    #births 
    N_birth = round(rbinom(n = 1, size = round(0.35*N[y-1,s]), prob = birth_rate))
    
    #and then can also decrement the number of calves by the calf survival rate 
    N_birth_surv = round(rbinom(n = 1, size = N_birth, prob = birth_surv))
    
    #the number in the next year
    N[y,s] = N_surv + N_birth_surv
  }
}

plot(x = 2018:2050, y = N[,1], type = "l")
for(i in 2:ncol(N)){
  lines(2018:2050, y = N[,i])
}

########################
# modifying to closely follow Stewart et al.
#######################


#############
#logistic growth
###########
#wedekin et al. 2017 estimates intrinsic growth for humpbacks in Brazil as 12%
#felix and haase 2025 estimaet intrinsic growth for humpbacks in Ecuador as 11.9%
#N(t+1) = N(t) + rN(t)*(1 - Nt/K) - H

N0 = 4973
K = 8000
r = 0.119
bycatch_min = 10 #3
bycatch_max = 70 #53
strike_min = 8
strike_max = 53
#anth_inc = seq(0, 5, length = 100) 
anth_inc = c(rep(3, 10), rep(0, 90))
#anth_inc = rep(0, length = 100)
growth_dec = rep(0, length = 100) #chronic stress over time on the growth rates 
#growth_dec = seq(0, 0.4, length = 100)

N = as.data.frame(matrix(nrow = 100, ncol = 100))
for(s in 1:ncol(N)){
  N[1,s] = N0
  
  for(y in 2:nrow(N)){
    # randomly select mortality 
   H = round(runif(n = 1, min = bycatch_min + (anth_inc[y]*bycatch_min), 
                   max = bycatch_max + (anth_inc[y]*bycatch_max))) +
     round(runif(n = 1, min = strike_min + (anth_inc[y]*strike_min), 
                 max = strike_max + (strike_max*anth_inc[y])))
    
    N[y,s] = N[y-1,s] + (r - r*growth_dec[y])*N[y-1,s]*(1 - N[y-1,s]/K) - H
  }
}

plot(x = 1:100, y = N[,1], type = "l")
for(i in 2:ncol(N)){
  lines(1:100, y = N[,i])
}

