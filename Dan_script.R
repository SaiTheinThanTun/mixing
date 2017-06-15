
library(deSolve)
library(TSA)

#non-reactive parameters
# define the number of weeks to run the model
#time step
dt<-1/12
startyear<-2010
stopyear<-2025
maxt<-stopyear-startyear
#calcs for time steps
times <- seq(0, maxt, by = dt)
tsteps<-length(times)


runFIRST <- function(initprev, scenario, param)
{
  
  parameters = c(scenario,
                 gamma = 0.1,         #recovery rate
                 alpha= .05,            #mixing between populations
                 delta = 100,            #relationship between transmission in pop1 and pop2
                 epsilonh=0.23,      # per bite probability of an infectious mosquito infecting a human
                 epsilonm=0.5,       # per bite probability of an infectious human infecting a mosquito
                 b=365/3,            # per mosquito rate of biting
                 deltam=365/14,      # rate leaving latent period          
                 gammam=365/10,       # mosquito death rate
                 bh = 1,                # bites per human per year
                 mu = 50,            #life expectancy
                 omega = 1/2,           #duration of immunity in years
                 lossd = 30,          #loss of drug prophylactic effect
                 #                nuTr = 14,            #days of infectiousness after ACT
                 timei = 2018,
                 primon = 1,            #primaquine use
                 nuTrp = 7,           # days of infectiosness after treatment ACT+primaquine [N]
                 ps = 90,            # % of all non-immune new infections that are clinical [N]
                 pr = 20,            # % of all immune new infections that are clinical [N]
                 bh = 7,
                 muC = 1,                #imported cases through inmigration
                 nuC = 14,                     # days of symptoms in the absence of treatment [N]
                 param)
  
  # MODEL INITIAL CONDITIONS
  # population size
  initP1<-10000 
  initP2<-10000 
  
  initS1<-0.5*(1-initprev)*initP1
  initI1<-initprev*initP2
  initR1<-0.5*(1-initprev)*initP2
  initS2<-0.5*(1-initprev)*initP2
  initI2<-initprev*initP2
  initR2<-0.5*(1-initprev)*initP2
  
  state <- c(Y = 0, Cinc_det1 = 0, Cinc_tot1 = 0,
             S1 = initS1, I1 = initI1, R1 = initR1, 
             Cinc_det2 = 0, Cinc_tot2 = 0,
             S2 = initS2, I2 = initI2, R2 = initR2 
  )
  
  sir <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      
      beta1<-b*epsilonh*epsilonm*bh/((bh*epsilonh+deltam)*(gammam/(gammam+deltam)))
      beta2=(beta1 *delta)    
      
      mu = 1/mu              #life expectancy in years
      mu_out = mu + muC          #all the ways you can leave a state
      #    lossd<-365/lossd
      #    lossd<-1/((1/lossd)-(1/nuTr))
      #    nTr<-365/nuTr
      #   nTrp<-365/nuTrp
      #    nuTr<- primon*((Y<timei)*nTr+(Y>timei)*nTrp)+(1-primon)*nTr
      muC<-muC/1000
      nuC<-365/nuC
      ps = ps/100
      pr = pr/100
      covEDATi<-0.9*covEDATi/100
      covEDAT0<-0.9*covEDAT0/100    
      P1 <- (S1+I1+R1)
      P2 <- (S2+I2+R2)
      
      lam1 = (beta1*(I1/(P1+P2)) + alpha*beta1*(I2/(P1+P2)))
      lam2 = (beta2*(I2/(P1+P2)) + alpha*beta2*(I1/(P1+P2)))
      
      timei<-timei-startyear
      wsiEDAT<-(1-(Y<=timei))*(Y<=(timei+EDATscale))*((Y-timei)/EDATscale)+1*(Y>=(timei+EDATscale))
      covEDAT<-(1-wsiEDAT)*covEDAT0+wsiEDAT*covEDATi
      tau <- covEDAT
      
      # rate of change
      dY <- 1
      
      #pop 1
      dCinc_det1 <- ps*tau*lam1*S1+pr*tau*lam1*R1+pr*tau*lam1*I1
      dCinc_tot1 <- ps*lam1*S1+pr*lam1*R1+pr*lam1*I1
      dS1 <- mu*P1 - mu_out*S1 + omega*R1 - lam1*S1
      dI1 <- muC*P1 - mu_out*I1 + lam1*S1 - nuC*I1 + lam1*R1
      dR1 <- - mu_out*R1 - omega*R1 - lam1*R1    + nuC*I1 
      
      #pop 2
      dCinc_det2 <- ps*tau*lam2*S2+pr*tau*lam2*R2+pr*tau*lam2*I2
      dCinc_tot2 <- ps*lam2*S2+pr*lam2*R2+pr*lam2*I2
      dS2 <- mu*P2 - mu_out*S2 + omega * R2 -lam1 * S2
      dI2 <- muC*P2 - mu_out*I2 + lam1*S2 - nuC*I2 + lam1*R2 
      dR2 <- - mu_out*R2 - omega*R2 - lam1*R2    + nuC*I2    
      
      list(c(dY, dCinc_det1, dCinc_tot1, dS1, dI1, dR1, dCinc_det2, dCinc_tot2, dS2, dI2, dR2))  
    })
  }
  
  out = ode(y = state, times = times, fun = sir, parms = parameters)
  
  # MODEL OUTPUTS
  ipop <- c(5,6,7, 10,11,12)
  iinc_det <- c(3, 8)
  iinc_tot <- c(4, 9)
  iprev <- c(6, 11)
  
  
  # population
  times<-out[,1]+startyear
  pop<-rowSums(out[,ipop])
  iinc_det <- rowSums(out[,iinc_det])
  iinc_tot <- rowSums(out[,iinc_tot])
  
  # clinical incidence detected per 1000 per month
  tci_det <- out[,iinc_det]
  clinmonth_det <- tci_det
  clinmonth_det[1] <- 0
  clinmonth_det[2:length(times)] <- 1000*(tci_det[2:length(times)] - tci_det[1:(length(times)-1)])/pop[2:length(times)]
  
  # clinical incidence total per 1000 per month
  tci_tot <- out[,iinc_tot]
  clinmonth_tot <- tci_tot
  clinmonth_tot[1] <- 0
  clinmonth_tot[2:length(times)] <- 1000*(tci_tot[2:length(times)] - tci_tot[1:(length(times)-1)])/pop[2:length(times)]
  
  
  # % prevalence
  prevalence <- 100*rowSums(out[,iprev])/pop
  GMSout<-matrix(NA,nrow=length(times),ncol=4)
  GMSout[,1]<-times
  GMSout[,2]<-clinmonth_det
  GMSout[,3]<-clinmonth_tot
  GMSout[,4]<-prevalence
  
  inc1 = out[,6]
  inc2 = out[,11]  
  detinc1 = out[,3]
  detinc2 = out[,8]
  time = out[,1]
  
  
  return(GMSout)
  
} 


parametersR <- c(
  EDATscale = 1,
  covEDATi = 90,
  covEDAT0 = 30)

API <- 2.5
# initial prevalence
initprevR <- 0.001*API


scenario_0<-c(EDATon = 0,
              ITNon = 0,
              RCDon = 0,
              RCDcoex = 0,
              IRSon = 0,
              MDAon = 0,
              primon = 0,
              MSATon = 0)

scenario_i<-c(EDATon = 1,
              ITNon = 0,
              RCDon = 0,
              RCDcoex = 0,
              IRSon = 0,
              MDAon = 0,
              primon = 0,
              MSATon = 0)

GMSout0 <- runFIRST(initprevR, scenario_0, parametersR)

par(mfrow=c(1,3))
#plot(GMSouti[,1], GMSouti[,2], type="l")
plot(GMSout0[,1], GMSout0[,2], type="l")
plot(GMSout0[,1], GMSout0[,3], type="l")
plot(GMSout0[,1], GMSout0[,4], type="l")


###################
#parameters <- c(gamma = 0.1,         #recovery rate
alpha= .5,             #mixing between populations
delta = 2,            #relationship between transmission in pop1 and pop2
epsilonh=0.23,      # per bite probability of an infectious mosquito infecting a human
epsilonm=0.5,       # per bite probability of an infectious human infecting a mosquito
b=365/3,            # per mosquito rate of biting
deltam=365/14,      # rate leaving latent period          
gammam=365/10,       # mosquito death rate
bh = 1,                # bites per human per year
mu = 50,            #life expectancy
omega = 1/2,           #duration of immunity in years
lossd = 30,          #loss of drug prophylactic effect
nuTr = 14,            #days of infectiousness after ACT
timei = 2018,
primon = 1,            #primaquine use
nuTrp = 7           # days of infectiosness after treatment ACT+primaquine [N]

)

#times = seq(2010, 2025, 1/12)

#out <- as.data.frame(ode(y = initprevR, times = times, func = runFIRST, parms = parametersR))

#out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))

out$time <- NULL

matplot([,1], GMSout0, type = "l", xlab = "Time", ylab = "Ss, Is and Rs", 
        main = "SIR Model", lwd = 2, lty = 1, bty = "l", col = 2:4)
