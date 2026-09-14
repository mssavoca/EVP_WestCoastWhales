##########################################################################
# Simulation Code for Density-Dependent and -Independent Effects of 
# Environmental Stochasticity on Population Dynamics
# From:
# The future of baleen whales: recoveries, environmental constraints, and climate change
# Joshua D. Stewart, M. Tim Tinker, Robert L. Brownell Jr., Andrew J. Read
#
# Simulation code developed by M. Tim Tinker and Joshua D. Stewart
#


# Required packages:

library(tidyverse)
library(boot)
library(patchwork)

# Note: vital rates are calculated using an instantaneous mortality (hazards) formulation
# to facilitate comparison of density-dependent vs. density-independent env. stochasticity  

##########################################################################
# Environmental Variability Only----------------------------------------
# Version a) effects of environmental variation are primarily density-dependent 
##########################################################################

# Important to set seed so density-dependent and -independent simulations have the 
# same set of environmental effects
set.seed(123)

n.sims <- 1000
n.years <- 100
sig = 0.3     # variation in survival rates caused by environmental conditions (SD in log hazard ratio)
sigB = sig/3   # variation in birth rates caused by environmental conditions (SD in log hazard ratio)
ObsCV <- 0.10  # CV for observation error of abundance estimates
rho = 0.25     # temporal autocorrelation in env. conditions

MeanK<-25000
theta = 2   # logistic shape parameter

LHz0 = log(.031)  # base log hazard: corresponds to survival rate of 0.97 as N -> 0, and lambda = 1.066
B0 = .1            # base birth rate: realized birth rate as N -> 0 
BHz0 = log(-log(B0))

# Density dependent parameters alpha (affects realized birth rate) and gamma (affects survival)
alph = .19 
gamm = .66

# Anthropogenic hazards (for depleted pop)
Hz_anth = 0.05

# Convenience functions for calculating birth and survival rates
B_calc = function(BHz0,alph,epB,sgB,N,K,thta){ 
  B = exp(-exp(BHz0 + (alph + (epB - (sgB^2)/2)) * (N/K)^thta))
  return(B) 
}

S_calc = function(LHz0,gamm,epS,sgS,N,K,thta,Hz_anth){ 
  S = exp(-(exp(LHz0 + (gamm + (epS - (sgS^2)/2)) * (N/K)^thta) + Hz_anth))
  return(S) 
}

# calculate max growth rate at low density (this value should be 1.066)
br0 = B_calc(BHz0,alph,((sigB^2)/2),sigB,100,MeanK,theta)
sv0 = S_calc(LHz0,gamm,((sig^2)/2),sig,100,MeanK,theta,0)
lambda_0 = (1 + br0) * sv0; round(lambda_0,3)
# calculate max growth rate at K (this value should be 1)
brK = B_calc(BHz0,alph,((sigB^2)/2),sigB,MeanK,MeanK,theta)
svK = S_calc(LHz0,gamm,((sig^2)/2),sig,MeanK,MeanK,theta,0)
lambda_K = (1 + brK) * svK ; round(lambda_K,3)

Scenario_names = c("Depleted","Recovered")

# Function to find adjustment factor for realized K associated with set of env. conditions
findK_ddstoch = function(x,eb,es){ 
  K_adj = x
  br = exp(-exp(BHz0 + ( alph + eb ) * K_adj^theta )) # - (sgB^2)/2
  sv = exp(-exp(LHz0 + ( gamm + es ) * K_adj^theta )) #  - (sg^2)/2 
  dif = (1 - (1 + br) * sv)^2
  return(dif) 
}
print("This value should be equal to 1:")
adj = optimize(findK_ddstoch, c(0.1, 5), eb = 0, es = 0, tol = 0.00000001, maximum = FALSE); format(adj$minimum,digits = 3)

# Create simulation objects
K <- eps <- epsB <- br_dep <- br_rec <- sv_dep <- sv_rec <- N_dep <- N_rec<- N_dep_obs <- N_rec_obs <- array(dim=c(n.sims,n.years))

N_dep[,1] <- 3000 # Year 1 abundance
N_rec[,1] <- 3000 # Year 1 abundance

for(j in 1:n.sims){
  
  K[j,1] <- MeanK # Year 1 K
  eps[j,1] <- 0
  epsB[j,1] <- 0
  
  set.seed(123+j)
  zrnd = rnorm(n.years-1,0,1)
  
  for(i in 1:(n.years-1)){
    
    # Environmental conditions, with autocorrelation (AR[1])
    eps[j,i+1] <- rho * eps[j,i] + (1 - rho) * (zrnd[i]*sig)
    epsB[j,i+1] <- rho * epsB[j,i] + (1 - rho) * (zrnd[i]*sigB)
    # Find realized K value that corresponds to these conditions
    adj = optimize(findK_ddstoch, c(0.1, 5), eb = epsB[j,i+1], es = eps[j,i+1], tol = 0.00000001, maximum = FALSE)
    K[j,i+1] <- MeanK * adj$minimum
    
    # Calculate annual birth and survival rates for high anthropogenic mortality scenario (depleted)
    br_dep[j,i] = B_calc(BHz0,alph,epsB[j,i+1],sigB,N_dep[j,i],MeanK,theta)
    sv_dep[j,i] = S_calc(LHz0,gamm,eps[j,i+1],sig,N_dep[j,i],MeanK,theta,Hz_anth)  
    
    # Calculate next year's abundance
    N_dep[j,i+1] <- N_dep[j,i] * (1 + br_dep[j,i]) * sv_dep[j,i]
    
    # Calculate annual birth and natural mortality rates for low anthropogenic mortality scenario (recovered)
    br_rec[j,i] = B_calc(BHz0,alph,epsB[j,i+1],sigB,N_rec[j,i],MeanK,theta)
    sv_rec[j,i] = S_calc(LHz0,gamm,eps[j,i+1],sig,N_rec[j,i],MeanK,theta,0)  
    
    # Calculate next year's abundance
    N_rec[j,i+1] <- N_rec[j,i] * (1 + br_rec[j,i]) * sv_rec[j,i]
    
  }#i
  
  # Observation error
  for(i in 1:n.years){
    N_dep_obs[j,i] <- rnorm(1,mean=N_dep[j,i],sd=N_dep[j,i]*ObsCV)
    N_rec_obs[j,i]  <- rnorm(1,mean=N_rec[j,i],sd=N_rec[j,i]*ObsCV)
  }#i
  
}#j

# Process simulation outputs for plotting:
KLong  <- as.data.frame(K)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="K") %>% mutate(Year=as.numeric(substr(Year,2,4)))

N_dep_obs_Lng <- as.data.frame(N_dep_obs)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="Obs") %>% mutate(Year=as.numeric(substr(Year,2,4)))
N_rec_obs_Lng <- as.data.frame(N_rec_obs)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="Obs") %>% mutate(Year=as.numeric(substr(Year,2,4)))

N_dep_Lng <- as.data.frame(N_dep)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="N") %>% mutate(Year=as.numeric(substr(Year,2,4))) %>% left_join(N_dep_obs_Lng)
N_rec_Lng <- as.data.frame(N_rec)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="N") %>% mutate(Year=as.numeric(substr(Year,2,4))) %>% left_join(N_rec_obs_Lng)

B_dep_Lng <- as.data.frame(br_dep)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="B") %>% mutate(Year=as.numeric(substr(Year,2,4)))
B_rec_Lng <- as.data.frame(br_rec)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="B") %>% mutate(Year=as.numeric(substr(Year,2,4))) 

S_dep_Lng <- as.data.frame(sv_dep)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="S") %>% mutate(Year=as.numeric(substr(Year,2,4)))
S_rec_Lng <- as.data.frame(sv_rec)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="S") %>% mutate(Year=as.numeric(substr(Year,2,4))) 

# Select one simulation to highlight in figures:
SimSelect <- 20

# If you want to plot a subset of simulations with observation error included:
ObsYears <- sample(1:100,25)

# Figure plotting:
ClimVar_DD <- ggplot(KLong[KLong$Simulation <= 100,],aes(x=Year,y=K,group=Simulation))+
  geom_line(color='gray50',alpha=0.1)+
  geom_line(data=N_rec_Lng[N_rec_Lng$Simulation <= 100,],aes(x=Year,y=N,group=Simulation),inherit.aes=F,color='navy',alpha=0.05)+
  geom_line(data=N_dep_Lng[N_dep_Lng$Simulation <= 100,],aes(x=Year,y=N,group=Simulation),inherit.aes=F,color='orange',alpha=0.03)+
  geom_line(data=filter(KLong,Simulation==SimSelect),aes(x=Year,y=K),inherit.aes=F,color='black',linewidth=0.5,alpha=0.75)+
  geom_line(data=filter(N_rec_Lng,Simulation==SimSelect),aes(x=Year,y=N),inherit.aes=F,color='navy',linewidth=1)+
  geom_line(data=filter(N_dep_Lng,Simulation==SimSelect),aes(x=Year,y=N),inherit.aes=F,color='orange',linewidth=1)+
  coord_cartesian(ylim=c(0,min(max(KLong$K),2.5*MeanK)))+
  theme_minimal() +
  ylab("Abundance")+
  ggtitle("A) Density-Dependent Effects (Applied to K)",subtitle="Population Response to Environmental Stochasticity")#+


Depleted_Summary_DD <- N_dep_Lng %>% group_by(Simulation) %>% filter(Year>50) %>% summarize(NMean = mean(N),NSD = sd(N), 
                                                           MeanLambda = mean(N[2:n()]/N[1:(n()-1)]),SDLambda = sd(N[2:n()]/N[1:(n()-1)]),                              
                                                           MaxDecl=max((N[2:n()]-N[1:(n()-1)])/N[1:(n()-1)])) %>% mutate(NCV = NSD/NMean) %>% mutate(Scenario="Depleted")
Recovered_Summary_DD <- N_rec_Lng %>% group_by(Simulation) %>% filter(Year>50) %>% summarize(NMean = mean(N),NSD = sd(N), 
                                                           MeanLambda = mean(N[2:n()]/N[1:(n()-1)]),SDLambda = sd(N[2:n()]/N[1:(n()-1)]),                               
                                                           MaxDecl=max((N[2:n()]-N[1:(n()-1)])/N[1:(n()-1)])) %>% mutate(NCV = NSD/NMean) %>% mutate(Scenario="Recovered")

CVSummaries_DD <- bind_rows(Depleted_Summary_DD,Recovered_Summary_DD)

AbundVar_DD <- ggplot(CVSummaries_DD,aes(x=NCV,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlab("Abundance Coefficient of Variation")

LambdaVar_DD <- ggplot(CVSummaries_DD,aes(x=SDLambda,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(0,NA)+
  xlab(expression(sigma))+
  ggtitle(label="",subtitle="B) Variation in annual growth rates")

MaxDecl_DD<-ggplot(CVSummaries_DD,aes(x=-MaxDecl,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(NA,0)+
  xlab("Proportional Change")+
  ggtitle(label="",subtitle="C) Maximum Proportional Decline")

print(ClimVar_DD / (LambdaVar_DD+MaxDecl_DD) +plot_layout(guides = 'collect',heights=c(3,1,1)))


BRate_DD <- ggplot()+
  geom_line(data=B_rec_Lng[B_rec_Lng$Simulation <= 100,],aes(x=Year,y=B,group=Simulation),inherit.aes=F,color='navy',alpha=0.03)+
  geom_line(data=B_dep_Lng[B_dep_Lng$Simulation <= 100,],aes(x=Year,y=B,group=Simulation),inherit.aes=F,color='orange',alpha=0.03)+
  geom_line(data=filter(B_rec_Lng,Simulation==SimSelect),aes(x=Year,y=B),inherit.aes=F,color='navy',linewidth=1)+
  geom_line(data=filter(B_dep_Lng,Simulation==SimSelect),aes(x=Year,y=B),inherit.aes=F,color='orange',linewidth=1)+
  theme_minimal() +
  ylim(0,0.15)+
  ggtitle(label="",subtitle="D) Birth Rate")


SRate_DD <- ggplot()+
  geom_line(data=S_rec_Lng[S_rec_Lng$Simulation <= 100,],aes(x=Year,y=S,group=Simulation),inherit.aes=F,color='navy',alpha=0.03)+
  geom_line(data=S_dep_Lng[S_dep_Lng$Simulation <= 100,],aes(x=Year,y=S+Hz_anth,group=Simulation),inherit.aes=F,color='orange',alpha=0.03)+
  geom_line(data=filter(S_rec_Lng,Simulation==SimSelect),aes(x=Year,y=S),inherit.aes=F,color='navy',linewidth=1)+
  geom_line(data=filter(S_dep_Lng,Simulation==SimSelect),aes(x=Year,y=S+Hz_anth),inherit.aes=F,color='orange',linewidth=1)+
  theme_minimal() +
  ylim(0.86,1)+
  ggtitle(label="",subtitle="E) Survival Rate (Natural)")

print(BRate_DD / (SRate_DD) +plot_layout(guides = 'collect',heights=c(1,1)))

print(ClimVar_DD / (LambdaVar_DD+MaxDecl_DD) / (BRate_DD + SRate_DD) +plot_layout(guides = 'collect',heights=c(3,1,1,1,1)))


##########################################################################
# Environmental Variability Only----------------------------------------
# Version b) effects of environmental variation are primarily density-INDEPENDENT
##########################################################################


# Set same seed as above simulation for consistent environmental effects
set.seed(123)

n.sims <- 1000
n.years <- 100
sig = 0.3     # variation in survival rates caused by environmental conditions (SD in log hazard ratio)
sigB = sig/3   # variation in birth rates caused by environmental conditions (SD in log hazard ratio)
ObsCV <- 0.10  # CV for observation error of abundance estimates
rho = 0.25     # temporal autocorrelation in env. conditions

MeanK<-25000
theta = 2

LHz0 = log(.031)  # base log hazard: corresponds to survival rate of 0.97 as N -> 0, and lambda = 1.066
B0 = .1            # base birth rate: realized birth rate as N -> 0 
BHz0 = log(-log(B0))

# Density dependent parameters alpha (affects realized birth rate) and gamma (affects survival)
alph = .19 
gamm = .66

# Anthropogenic hazards (for depleted pop)
Hz_anth = 0.05

# calculate max growth rate at low density (this value should be 1.066)
br0 = B_calc(BHz0,alph,((sigB^2)/2),sigB,100,MeanK,theta)
sv0 = S_calc(LHz0,gamm,((sig^2)/2),sig,100,MeanK,theta,0)
lambda_0 = (1 + br0) * sv0; round(lambda_0,3)
# calculate max growth rate at K (this value should be 1)
brK = B_calc(BHz0,alph,((sigB^2)/2),sigB,MeanK,MeanK,theta)
svK = S_calc(LHz0,gamm,((sig^2)/2),sig,MeanK,MeanK,theta,0)
lambda_K = (1 + brK) * svK ; round(lambda_K,3)

Scenario_names = c("Depleted","Recovered")

eps <- epsB <- br_dep <- br_rec <- sv_dep <- sv_rec <- N_dep <- N_rec<- N_dep_obs <- N_rec_obs <- array(dim=c(n.sims,n.years))

N_dep[,1] <- 3000 # Year 1 abundance
N_rec[,1] <- 3000 # Year 1 abundance

for(j in 1:n.sims){
  
  eps[j,1] <- 0
  epsB[j,1] <- 0
  
  set.seed(123+j)
  zrnd = rnorm(n.years-1,0,1)
  
  for(i in 1:(n.years-1)){
    
    # Environmental conditions, with autocorrelation (AR[1])
    eps[j,i+1] <- rho * eps[j,i] + (1 - rho) * (zrnd[i]*sig)
    epsB[j,i+1] <- rho * epsB[j,i] + (1 - rho) * (zrnd[i]*sigB)
    
    # Calculate annual birth and survival rates for high anthropogenic mortality scenario (depleted)
    br_dep[j,i] = B_calc(BHz0,alph,epsB[j,i+1],sigB,N_dep[j,i],MeanK,theta)
    sv_dep[j,i] = S_calc(LHz0,gamm,eps[j,i+1],sig,N_dep[j,i],MeanK,theta,Hz_anth)  
    
    # Calculate next year's abundance
    N_dep[j,i+1] <- N_dep[j,i] * (1 + br_dep[j,i]) * sv_dep[j,i]
    
    # Calculate annual birth and natural mortality rates for low anthropogenic mortality scenario (recovered)
    br_rec[j,i] = B_calc(BHz0,alph,epsB[j,i+1],sigB,N_rec[j,i],MeanK,theta)
    sv_rec[j,i] = S_calc(LHz0,gamm,eps[j,i+1],sig,N_rec[j,i],MeanK,theta,0)  
    
    # Calculate next year's abundance
    N_rec[j,i+1] <- N_rec[j,i] * (1 + br_rec[j,i]) * sv_rec[j,i]
    
  }#i
  
  # Observation error
  for(i in 1:n.years){
    N_dep_obs[j,i] <- rnorm(1,mean=N_dep[j,i],sd=N_dep[j,i]*ObsCV)
    N_rec_obs[j,i]  <- rnorm(1,mean=N_rec[j,i],sd=N_rec[j,i]*ObsCV)
  }#i
  
}#j

N_dep_obs_Lng <- as.data.frame(N_dep_obs)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="Obs") %>% mutate(Year=as.numeric(substr(Year,2,4)))
N_rec_obs_Lng <- as.data.frame(N_rec_obs)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="Obs") %>% mutate(Year=as.numeric(substr(Year,2,4)))

N_dep_Lng <- as.data.frame(N_dep)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="N") %>% mutate(Year=as.numeric(substr(Year,2,4))) %>% left_join(N_dep_obs_Lng)
N_rec_Lng <- as.data.frame(N_rec)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="N") %>% mutate(Year=as.numeric(substr(Year,2,4))) %>% left_join(N_rec_obs_Lng)

B_dep_Lng <- as.data.frame(br_dep)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="B") %>% mutate(Year=as.numeric(substr(Year,2,4)))
B_rec_Lng <- as.data.frame(br_rec)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="B") %>% mutate(Year=as.numeric(substr(Year,2,4))) 

S_dep_Lng <- as.data.frame(sv_dep)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="S") %>% mutate(Year=as.numeric(substr(Year,2,4)))
S_rec_Lng <- as.data.frame(sv_rec)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="S") %>% mutate(Year=as.numeric(substr(Year,2,4))) 


SimSelect <- 20
ObsYears <- sample(1:100,25)

# Figure plotting:

ClimVar_DI <- ggplot()+
  geom_line(data=N_rec_Lng[N_rec_Lng$Simulation <= 100,],aes(x=Year,y=N,group=Simulation),inherit.aes=F,color='navy',alpha=0.03)+
  geom_line(data=N_dep_Lng[N_dep_Lng$Simulation <= 100,],aes(x=Year,y=N,group=Simulation),inherit.aes=F,color='orange',alpha=0.03)+
  geom_line(data=filter(N_rec_Lng,Simulation==SimSelect),aes(x=Year,y=N),inherit.aes=F,color='navy',linewidth=1)+
  geom_line(data=filter(N_dep_Lng,Simulation==SimSelect),aes(x=Year,y=N),inherit.aes=F,color='orange',linewidth=1)+
  geom_hline(yintercept = MeanK, color='gray50',linetype="dashed")+
  coord_cartesian(ylim=c(0,min(max(KLong$K),2.5*MeanK)))+
  theme_minimal() +
  ylab("Abundance")+
  ggtitle("F) Density-Independent Effects (Applied to r)")#,subtitle="Population Response to Environmental Stochasticity")#+

Depleted_Summary_DI <- N_dep_Lng %>% group_by(Simulation) %>% filter(Year>50) %>% summarize(NMean = mean(N),NSD = sd(N), 
                                                                                            MeanLambda = mean(N[2:n()]/N[1:(n()-1)]),SDLambda = sd(N[2:n()]/N[1:(n()-1)]),                              
                                                                                            MaxDecl=max((N[2:n()]-N[1:(n()-1)])/N[1:(n()-1)])) %>% mutate(NCV = NSD/NMean) %>% mutate(Scenario="Depleted")
Recovered_Summary_DI <- N_rec_Lng %>% group_by(Simulation) %>% filter(Year>50) %>% summarize(NMean = mean(N),NSD = sd(N), 
                                                                                             MeanLambda = mean(N[2:n()]/N[1:(n()-1)]),SDLambda = sd(N[2:n()]/N[1:(n()-1)]),                               
                                                                                             MaxDecl=max((N[2:n()]-N[1:(n()-1)])/N[1:(n()-1)])) %>% mutate(NCV = NSD/NMean) %>% mutate(Scenario="Recovered")
CVSummaries_DI <- bind_rows(Depleted_Summary_DI,Recovered_Summary_DI)

AbundVar_DI <- ggplot(CVSummaries_DI,aes(x=NCV,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlab("Abundance Coefficient of Variation")

LambdaVar_DI <- ggplot(CVSummaries_DI,aes(x=SDLambda,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(0,NA)+
  xlab(expression(sigma))+
  ggtitle(label="",subtitle="G) Variation in annual growth rates")

MaxDecl_DI<-ggplot(CVSummaries_DI,aes(x=-MaxDecl,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(NA,0)+
  xlab("Proportional Change")+
  ggtitle(label="",subtitle="H) Maximum Proportional Decline")

print(ClimVar_DI / (LambdaVar_DI+MaxDecl_DI) +plot_layout(guides = 'collect',heights=c(3,1,1)))

BRate_DI <- ggplot()+
  geom_line(data=B_rec_Lng[B_rec_Lng$Simulation <= 100,],aes(x=Year,y=B,group=Simulation),inherit.aes=F,color='navy',alpha=0.03)+
  geom_line(data=B_dep_Lng[B_dep_Lng$Simulation <= 100,],aes(x=Year,y=B,group=Simulation),inherit.aes=F,color='orange',alpha=0.03)+
  geom_line(data=filter(B_rec_Lng,Simulation==SimSelect),aes(x=Year,y=B),inherit.aes=F,color='navy',linewidth=1)+
  geom_line(data=filter(B_dep_Lng,Simulation==SimSelect),aes(x=Year,y=B),inherit.aes=F,color='orange',linewidth=1)+
  theme_minimal() +
  ylim(0,0.15)+
  ggtitle(label="",subtitle="I) Birth Rate")


SRate_DI <- ggplot()+
  geom_line(data=S_rec_Lng[S_rec_Lng$Simulation <= 100,],aes(x=Year,y=S,group=Simulation),inherit.aes=F,color='navy',alpha=0.03)+
  geom_line(data=S_dep_Lng[S_dep_Lng$Simulation <= 100,],aes(x=Year,y=S+Hz_anth,group=Simulation),inherit.aes=F,color='orange',alpha=0.03)+
  geom_line(data=filter(S_rec_Lng,Simulation==SimSelect),aes(x=Year,y=S),inherit.aes=F,color='navy',linewidth=1)+
  geom_line(data=filter(S_dep_Lng,Simulation==SimSelect),aes(x=Year,y=S+Hz_anth),inherit.aes=F,color='orange',linewidth=1)+
  theme_minimal() +
  ylim(0.86,1)+
  ggtitle(label="",subtitle="J) Survival Rate (Natural)")

print(BRate_DI / (SRate_DI) +plot_layout(guides = 'collect',heights=c(1,1)))

print(ClimVar_DI / (LambdaVar_DI+MaxDecl_DI) / (BRate_DI + SRate_DI) +plot_layout(guides = 'collect',heights=c(3,1,1,1,1)))

# Manuscript Figure 2 layout and plotting:

design <- "
AAFF
AAFF
BCGH
DEIJ
DEIJ
"
ClimVar_DD +LambdaVar_DD +MaxDecl_DD +BRate_DD +SRate_DD +ClimVar_DI +LambdaVar_DI +MaxDecl_DI +BRate_DI +SRate_DI +plot_layout(design=design,guides = 'collect')




##########################################################################
## Environmental Variability PLUS Climate Reduction in K-----------------
# Version a) effects of environmental variation are primarily density-dependent 
##########################################################################

# Set new seed (needs to be the same as next climate simulation)
set.seed(345)


K <- MeanKAdj <- eps <- epsB <- br_dep <- br_rec <- sv_dep <- sv_rec <- N_dep <- N_rec<- N_dep_obs <- N_rec_obs <- array(dim=c(n.sims,n.years))

N_dep[,1] <- 3000 # Year 1 abundance
N_rec[,1] <- 3000 # Year 1 abundance

ClimateStart <- 50 # Year in which climate-induced decline in K begins
ClimateBeta <- 75 # Beta distrib parm: average annual decline in env. conditions due to climate

for(j in 1:n.sims){
  
  K[j,1] <- MeanK # Year 1 K
  eps[j,1] <- 0
  epsB[j,1] <- 0
  
  set.seed(123+j)
  zrnd = rnorm(n.years-1,0,1)
  
  for(i in 1:(n.years-1)){
    
    # Reduced mean K from climate impacts
    if(i<ClimateStart){
      MeanKAdj[j,i] <- MeanK
    }else{
      MeanKAdj[j,i] <- MeanKAdj[j,i-1]*rbeta(1,ClimateBeta,1)
    } 
    
    # Environmental conditions, with autocorrelation (AR[1])
    eps[j,i+1] <- rho * eps[j,i] + (1 - rho) * (zrnd[i]*sig)
    epsB[j,i+1] <- rho * epsB[j,i] + (1 - rho) * (zrnd[i]*sigB)
    # Find realized K value that corresponds to these conditions
    adj = optimize(findK_ddstoch, c(0.1, 5), eb = epsB[j,i+1], es = eps[j,i+1], tol = 0.00000001, maximum = FALSE)
    K[j,i+1] <- MeanKAdj[j,i] * adj$minimum
    
    # Calculate annual birth and survival rates for high anthropogenic mortality scenario (depleted)
    
    br_dep[j,i] = B_calc(BHz0,alph,epsB[j,i+1],sigB,N_dep[j,i],MeanKAdj[j,i],theta)
    sv_dep[j,i] = S_calc(LHz0,gamm,eps[j,i+1],sig,N_dep[j,i],MeanKAdj[j,i],theta,Hz_anth)  
    
    # Calculate next year's abundance
    N_dep[j,i+1] <- N_dep[j,i] * (1 + br_dep[j,i]) * sv_dep[j,i]
    
    # Calculate annual birth and natural mortality rates for low anthropogenic mortality scenario (recovered)
    br_rec[j,i] = B_calc(BHz0,alph,epsB[j,i+1],sigB,N_rec[j,i],MeanKAdj[j,i],theta)
    sv_rec[j,i] = S_calc(LHz0,gamm,eps[j,i+1],sig,N_rec[j,i],MeanKAdj[j,i],theta,0)  
    
    # Calculate next year's abundance
    N_rec[j,i+1] <- N_rec[j,i] * (1 + br_rec[j,i]) * sv_rec[j,i]
    
  }#i
  
  # Observation error
  for(i in 1:n.years){
    N_dep_obs[j,i] <- rnorm(1,mean=N_dep[j,i],sd=N_dep[j,i]*ObsCV)
    N_rec_obs[j,i]  <- rnorm(1,mean=N_rec[j,i],sd=N_rec[j,i]*ObsCV)
  }#i
  
}#j

# Simulation summaries for figure plotting:
KLong  <- as.data.frame(K)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="K") %>% mutate(Year=as.numeric(substr(Year,2,4)))

N_dep_obs_Lng <- as.data.frame(N_dep_obs)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="Obs") %>% mutate(Year=as.numeric(substr(Year,2,4)))
N_rec_obs_Lng <- as.data.frame(N_rec_obs)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="Obs") %>% mutate(Year=as.numeric(substr(Year,2,4)))

N_dep_Lng <- as.data.frame(N_dep)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="N") %>% mutate(Year=as.numeric(substr(Year,2,4))) %>% left_join(N_dep_obs_Lng)
N_rec_Lng <- as.data.frame(N_rec)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="N") %>% mutate(Year=as.numeric(substr(Year,2,4))) %>% left_join(N_rec_obs_Lng)

SimSelect <- 20
ObsYears <- sample(1:100,25)

# Figure plotting:

ClimDecl_DD <- ggplot(KLong[KLong$Simulation <= 100,],aes(x=Year,y=K,group=Simulation))+
  geom_line(color='gray50',alpha=0.1)+
  geom_line(data=N_rec_Lng[N_rec_Lng$Simulation <= 100,],aes(x=Year,y=N,group=Simulation),inherit.aes=F,color='navy',alpha=0.05)+
  geom_line(data=N_dep_Lng[N_dep_Lng$Simulation <= 100,],aes(x=Year,y=N,group=Simulation),inherit.aes=F,color='orange',alpha=0.03)+
  geom_line(data=filter(KLong,Simulation==SimSelect),aes(x=Year,y=K),inherit.aes=F,color='black',linewidth=0.5,alpha=0.75)+
  geom_line(data=filter(N_rec_Lng,Simulation==SimSelect),aes(x=Year,y=N),inherit.aes=F,color='navy',linewidth=1)+
  geom_line(data=filter(N_dep_Lng,Simulation==SimSelect),aes(x=Year,y=N),inherit.aes=F,color='orange',linewidth=1)+
  ylim(0,min(max(KLong$K),2.5*MeanK)) +
  ylab("Abundance")+
  theme_minimal() +
  ggtitle("A) Density-Dependent Effects (Applied to K)",subtitle="Population Response to Climate Decline")#+


Depleted_Summary_DD <- N_dep_Lng %>% group_by(Simulation) %>% filter(Year>50) %>% summarize(NMean = mean(N),NSD = sd(N), 
                                                                                         MeanLambda = mean(N[2:n()]/N[1:(n()-1)]),SDLambda = sd(N[2:n()]/N[1:(n()-1)]),                              
                                                                                         MaxDecl=max((N[2:n()]-N[1:(n()-1)])/N[1:(n()-1)])) %>% mutate(NCV = NSD/NMean) %>% mutate(Scenario="Depleted")
Recovered_Summary_DD <- N_rec_Lng %>% group_by(Simulation) %>% filter(Year>50) %>% summarize(NMean = mean(N),NSD = sd(N), 
                                                                                          MeanLambda = mean(N[2:n()]/N[1:(n()-1)]),SDLambda = sd(N[2:n()]/N[1:(n()-1)]),                               
                                                                                          MaxDecl=max((N[2:n()]-N[1:(n()-1)])/N[1:(n()-1)])) %>% mutate(NCV = NSD/NMean) %>% mutate(Scenario="Recovered")
CVSummaries_DD <- bind_rows(Depleted_Summary_DD,Recovered_Summary_DD)

LowMortDecSummary_DD <- N_rec_Lng %>% group_by(Simulation) %>% mutate(EndN = N[n()]) %>% filter(N==max(N)) %>% mutate(PropDecline = (EndN-N)/N) %>% mutate(Scenario="Recovered")
HighMortDecSummary_DD <- N_dep_Lng %>% group_by(Simulation) %>% mutate(EndN = N[n()]) %>% filter(N==max(N)) %>% mutate(PropDecline = (EndN-N)/N) %>% mutate(Scenario="Depleted")

DecSummaries_DD <- bind_rows(LowMortDecSummary_DD,HighMortDecSummary_DD)

AbundVar_DD <- ggplot(CVSummaries_DD,aes(x=NCV,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlab("Abundance Coefficient of Variation")

LambdaVar_DD <- ggplot(CVSummaries_DD,aes(x=SDLambda,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(0,NA)+
  xlab("Variation in annual growth rates")

MaxDecl_DD<-ggplot(CVSummaries_DD,aes(x=-MaxDecl,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(NA,0)+
  xlab("Maximum Proportional Decline")

PropDecl_DD<-ggplot(DecSummaries_DD,aes(x=PropDecline,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(NA,0.1)+
  xlab("Proportional Change")+
  ggtitle(label="",subtitle="C) Decline from Peak Abundance")

Inflect_DD<-ggplot(DecSummaries_DD,aes(x=Year-ClimateStart,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlab("Years Between Climate Inflection\nand Peak Abundance")+
  ggtitle(label="",subtitle="B) Responsiveness")

print(ClimDecl_DD / (Inflect_DD+PropDecl_DD) + plot_layout(guides = 'collect',heights=c(3,1,1)))



##########################################################################
## Environmental Variability PLUS Climate Reduction in Rmax-----------------
# Version b) effects of environmental variation are primarily density-INDEPENDENT
##########################################################################

# Set same seed as previous simulation:
set.seed(345)


M_HazAdj <- B_HazAdj <- eps <- epsB <- br_dep <- br_rec <- sv_dep <- sv_rec <- N_dep <- N_rec<- N_dep_obs <- N_rec_obs <- array(dim=c(n.sims,n.years))

N_dep[,1] <- 3000 # Year 1 abundance
N_rec[,1] <- 3000 # Year 1 abundance

ClimateStart <- 50 # Year in which climate-induced decline in K begins
ClimateBeta <- 150 # Beta distrib parm: average annual decline in env. conditions due to climate

for(j in 1:n.sims){
  
  eps[j,1] <- 0
  epsB[j,1] <- 0
  
  set.seed(123+j)
  zrnd = rnorm(n.years-1,0,1)
  
  for(i in 1:(n.years-1)){
    
    # Reduced base vital rates from climate impacts
    if(i<ClimateStart){
      B_HazAdj[j,i] <- BHz0
      M_HazAdj[j,i] <- LHz0
    }else{
      B_HazAdj[j,i] <- B_HazAdj[j,i-1] + rbeta(1,1,ClimateBeta)
      M_HazAdj[j,i] <- M_HazAdj[j,i-1] + rbeta(1,1,ClimateBeta)
    } 
    
    # Environmental conditions, with autocorrelation (AR[1])
    eps[j,i+1] <- rho * eps[j,i] + (1 - rho) * (zrnd[i]*sig)
    epsB[j,i+1] <- rho * epsB[j,i] + (1 - rho) * (zrnd[i]*sigB)
    
    # Calculate annual birth and survival rates for high anthropogenic mortality scenario (depleted)
    br_dep[j,i] = B_calc(B_HazAdj[j,i],alph,epsB[j,i+1],sigB,N_dep[j,i],MeanK,theta)
    sv_dep[j,i] = S_calc(M_HazAdj[j,i],gamm,eps[j,i+1],sig,N_dep[j,i],MeanK,theta,Hz_anth)  
    
    # Calculate next year's abundance
    N_dep[j,i+1] <- N_dep[j,i] * (1 + br_dep[j,i]) * sv_dep[j,i]
    
    # Calculate annual birth and natural mortality rates for low anthropogenic mortality scenario (recovered)
    br_rec[j,i] = B_calc(B_HazAdj[j,i],alph,epsB[j,i+1],sigB,N_rec[j,i],MeanK,theta)
    sv_rec[j,i] = S_calc(M_HazAdj[j,i],gamm,eps[j,i+1],sig,N_rec[j,i],MeanK,theta,0)  
    
    # Calculate next year's abundance
    N_rec[j,i+1] <- N_rec[j,i] * (1 + br_rec[j,i]) * sv_rec[j,i]
    
  }#i
  
  # Observation error
  for(i in 1:n.years){
    N_dep_obs[j,i] <- rnorm(1,mean=N_dep[j,i],sd=N_dep[j,i]*ObsCV)
    N_rec_obs[j,i]  <- rnorm(1,mean=N_rec[j,i],sd=N_rec[j,i]*ObsCV)
  }#i
  
}#j

# Simulation summaries for plotting:
N_dep_obs_Lng <- as.data.frame(N_dep_obs)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="Obs") %>% mutate(Year=as.numeric(substr(Year,2,4)))
N_rec_obs_Lng <- as.data.frame(N_rec_obs)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="Obs") %>% mutate(Year=as.numeric(substr(Year,2,4)))

N_dep_Lng <- as.data.frame(N_dep)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="N") %>% mutate(Year=as.numeric(substr(Year,2,4))) %>% left_join(N_dep_obs_Lng)
N_rec_Lng <- as.data.frame(N_rec)%>%mutate(Simulation=1:n.sims) %>% pivot_longer(cols=-c(Simulation),names_to=c("Year"),values_to="N") %>% mutate(Year=as.numeric(substr(Year,2,4))) %>% left_join(N_rec_obs_Lng)

SimSelect <- 20
ObsYears <- sample(1:100,25)

# Figure plotting:

ClimDecl_DI <- ggplot() + 
  geom_line(data=N_rec_Lng[N_rec_Lng$Simulation <= 100,],aes(x=Year,y=N,group=Simulation),inherit.aes=F,color='navy',alpha=0.03)+
  geom_line(data=N_dep_Lng[N_dep_Lng$Simulation <= 100,],aes(x=Year,y=N,group=Simulation),inherit.aes=F,color='orange',alpha=0.03)+
  geom_line(data=filter(N_rec_Lng,Simulation==SimSelect),aes(x=Year,y=N),inherit.aes=F,color='navy',linewidth=1)+
  geom_line(data=filter(N_dep_Lng,Simulation==SimSelect),aes(x=Year,y=N),inherit.aes=F,color='orange',linewidth=1)+
  ylim(0,min(max(KLong$K),2.5*MeanK)) +
  geom_hline(yintercept = MeanK, color='gray50',linetype="dashed")+
  ylab("Abundance")+
  theme_minimal() +
  ggtitle("D) Density-Independent Effects (Applied to r)",subtitle="")

Depleted_Summary_DI <- N_dep_Lng %>% group_by(Simulation) %>% filter(Year>50) %>% summarize(NMean = mean(N),NSD = sd(N), 
                                                                                            MeanLambda = mean(N[2:n()]/N[1:(n()-1)]),SDLambda = sd(N[2:n()]/N[1:(n()-1)]),                              
                                                                                            MaxDecl=max((N[2:n()]-N[1:(n()-1)])/N[1:(n()-1)])) %>% mutate(NCV = NSD/NMean) %>% mutate(Scenario="Depleted")
Recovered_Summary_DI <- N_rec_Lng %>% group_by(Simulation) %>% filter(Year>50) %>% summarize(NMean = mean(N),NSD = sd(N), 
                                                                                             MeanLambda = mean(N[2:n()]/N[1:(n()-1)]),SDLambda = sd(N[2:n()]/N[1:(n()-1)]),                               
                                                                                             MaxDecl=max((N[2:n()]-N[1:(n()-1)])/N[1:(n()-1)])) %>% mutate(NCV = NSD/NMean) %>% mutate(Scenario="Recovered")
CVSummaries_DI <- bind_rows(Depleted_Summary_DI,Recovered_Summary_DI)

LowMortDecSummary_DI <- N_rec_Lng %>% group_by(Simulation) %>% mutate(EndN = N[n()]) %>% filter(N==max(N)) %>% mutate(PropDecline = (EndN-N)/N) %>% mutate(Scenario="Recovered")
HighMortDecSummary_DI <- N_dep_Lng %>% group_by(Simulation) %>% mutate(EndN = N[n()]) %>% filter(N==max(N)) %>% mutate(PropDecline = (EndN-N)/N) %>% mutate(Scenario="Depleted")

DecSummaries_DI <- bind_rows(LowMortDecSummary_DI,HighMortDecSummary_DI)

AbundVar_DI <- ggplot(CVSummaries_DI,aes(x=NCV,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlab("Abundance Coefficient of Variation")

LambdaVar_DI <- ggplot(CVSummaries_DI,aes(x=SDLambda,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(0,NA)+
  xlab("Variation in annual growth rates")

MaxDecl_DI<-ggplot(CVSummaries_DI,aes(x=-MaxDecl,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(NA,0)+
  xlab("Maximum Proportional Decline")

PropDecl_DI<-ggplot(DecSummaries_DI,aes(x=PropDecline,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlim(NA,0.1)+
  xlab("Proportional Change") +
  ggtitle(label="",subtitle="F) Decline from Peak Abundance")

Inflect_DI<-ggplot(DecSummaries_DI,aes(x=Year-ClimateStart,fill=Scenario)) +
  geom_density(alpha=0.5)+
  scale_fill_manual(values=c('orange','navy'))+
  theme_minimal() +
  xlab("Years Between Climate Inflection\nand Peak Abundance")+
  ggtitle(label="",subtitle="E) Responsiveness")

print(ClimDecl_DI / (Inflect_DI+PropDecl_DI) + plot_layout(guides = 'collect',heights=c(3,1,1)))



# Manuscript Figure 3 layout and plotting:

design2 <- "
AAFF
AAFF
BCGH
"
ClimDecl_DD +Inflect_DD+PropDecl_DD +ClimDecl_DI +Inflect_DI+PropDecl_DI +plot_layout(design=design2,guides = 'collect')
