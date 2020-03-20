library('deSolve')
library('ggplot2')
library('reshape2')

# This model is based on the interactive SEIR model of Covid-19 developed by Dr Alison Lynn Hill
#   Model: https://alhill.shinyapps.io/COVID19seir/
#   Source Code: https://github.com/alsnhll/SEIR_COVID19/blob/master/COVID19seir/server.R

# Parameters

IncubPeriod=5  #Incubation period, days
DurMildInf=6 #Duration of mild infections, days
FracSevere=15/100 #Fraction of infections that are severe
FracCritical=5/100 #Fraction of infections that are critical
FracMild=1-FracSevere-FracCritical  #Fraction of infections that are mild
ProbDeath=0.4  #Probability of dying given critical infection
CFR=ProbDeath*FracCritical/100 #Case fatality rate (fraction of infections resulting in death)
TimeICUDeath=10 #Time from ICU admission to death, days
DurHosp=4 #Duration of hospitalization, days

N=10^(7)

b1=0.33/N
b2=0*b1
b3=0*b1
b=c(b1,b2,b3)

# Model Parameters
a=1/IncubPeriod

g1=(1/DurMildInf)*FracMild
p1=(1/DurMildInf)-g1

p2=(1/DurHosp)*(FracCritical/(FracSevere+FracCritical))
g2=(1/DurHosp)-p2

if(FracCritical==0){
  u=0
}else{
  u=(1/TimeICUDeath)*(CFR/FracCritical)
}

g3=(1/TimeICUDeath)-u
p = c(b1,b2,b3,a,g1,p1,p2,g2,u,g3)


SetODEs_SEIR=function(t,y,p){
  S = y[1]
  E = y[2]
  I1 = y[3]
  I2 = y[4]
  I3 = y[5]
  R = y[6]
  D = y[7]
  
  with(as.list(p),{
    
    # if (t>150 & t<500) b1 = 0.2/N
    # if (t>150) b1 = 0.2/N
    
    dS.dt = -(b1*I1+b2*I2+b3*I3)*S
    dE.dt=(b1*I1+b2*I2+b3*I3)*S-a*E
    dI1.dt=a*E-g1*I1-p1*I1
    dI2.dt=p1*I1-g2*I2-p2*I2
    dI3.dt=p2*I2-g3*I3-u*I3
    dR.dt=g1*I1+g2*I2+g3*I3
    dD.dt=u*I3
    
    return(list(c(dS.dt, dE.dt, dI1.dt, dI2.dt, dI3.dt, dR.dt, dD.dt)))
  })
}

# Initial values
S0 = 1e7 - 1
E0 = 0
I1_0 = 1
I2_0 = 0
I3_0 = 0
R = 0
D = 0
y0 = c(S0,E0,I1_0,I2_0,I3_0,R,D)

t = seq(from = 0, to = 1000, by = 1)
sol = ode(y = y0, times = t, func = SetODEs_SEIR, parms = p)
sol.df = as.data.frame(sol)
names(sol.df) = c('time','S','E','I1','I2','I3','R','D')

sol.df_long = melt(sol.df, id="time")
names(sol.df_long) = c('Day','State','Population')
ggplot(data = sol.df_long,
       aes(x = Day, y = Population, color=State)) + 
  geom_line() + ggtitle('COVID-19')



# Intervention - Social Distancing
SetODEs_SEIR_Intervention=function(t,y,p){
  S = y[1]
  E = y[2]
  I1 = y[3]
  I2 = y[4]
  I3 = y[5]
  R = y[6]
  D = y[7]
  
  with(as.list(p),{
    
    if (t>int_begin & t<int_end) b1 = int_level/N
    
    dS.dt = -(b1*I1+b2*I2+b3*I3)*S
    dE.dt=(b1*I1+b2*I2+b3*I3)*S-a*E
    dI1.dt=a*E-g1*I1-p1*I1
    dI2.dt=p1*I1-g2*I2-p2*I2
    dI3.dt=p2*I2-g3*I3-u*I3
    dR.dt=g1*I1+g2*I2+g3*I3
    dD.dt=u*I3
    
    return(list(c(dS.dt, dE.dt, dI1.dt, dI2.dt, dI3.dt, dR.dt, dD.dt)))
  })
}

# INTERVENTION - SOCIAL DISTANCING
int_level_vals = c(0.3,0.25,0.2,0.15,0.1)
I3.socialdistancing = data.frame(matrix(nrow=1001,ncol=length(int_level_vals)+1))
names(I3.socialdistancing) = c('time','0.3','0.25','0.2','0.15','0.1')
I3.socialdistancing$time = 1:1001
for (i in 1:length(int_level_vals)) {
  int_begin = 150
  int_end = 1000
  int_level = int_level_vals[i]
  p_new = c(p,int_begin,int_end,int_level)
  sol.temp = ode(y = y0, times = t, func = SetODEs_SEIR_Intervention, parms = p_new)
  sol.temp = as.data.frame(sol.temp)
  names(sol.temp) = c('time','S','E','I1','I2','I3','R','D')
  I3.socialdistancing[,i+1] = sol.temp$I3
}
I3.socialdistancing_long = melt(I3.socialdistancing, id="time")
names(I3.socialdistancing_long) = c('Day','Intervention_Level','CriticalCases')
ggplot(data = I3.socialdistancing_long,
       aes(x = Day, y = CriticalCases, color=Intervention_Level)) + 
  geom_line() + ggtitle('Intervention Levels')

# LENGTH OF INTERVENTION
int_end_vals = c(14, # 2 weeks
                 30, # 1 month
                 60, # 2 months
                 120, # 4 months
                 180, # 6 months
                 270, # 9 months
                 365 # 1 year
                 )
I3.interventionlength = data.frame(matrix(nrow=1001,ncol=length(int_end_vals)+1))
names(I3.interventionlength) = c('time','2weeks','1month','2months','4months','6months','9months','12months')
I3.interventionlength$time = 1:1001
for (i in 1:length(int_end_vals)) {
  int_begin = 150
  int_end = int_begin+int_end_vals[i]
  int_level = 0.2
  p_new = c(p,int_begin,int_end,int_level)
  sol.temp = ode(y = y0, times = t, func = SetODEs_SEIR_Intervention, parms = p_new)
  sol.temp = as.data.frame(sol.temp)
  names(sol.temp) = c('time','S','E','I1','I2','I3','R','D')
  I3.interventionlength[,i+1] = sol.temp$I3
}
I3.interventionlength_long = melt(I3.interventionlength, id="time")
names(I3.interventionlength_long) = c('Day','Intervention_Length','CriticalCases')
ggplot(data = I3.interventionlength_long,
       aes(x = Day, y = CriticalCases, color=Intervention_Length)) + 
  geom_line() + ggtitle('Intervention Lengths')

# Intervention - Social Distancing with Loss of Immunity
SetODEs_SEIR_ImmunityLoss=function(t,y,p){
  S = y[1]
  E = y[2]
  I1 = y[3]
  I2 = y[4]
  I3 = y[5]
  R = y[6]
  D = y[7]
  
  with(as.list(p),{
    
    # if (t>150 & t<500) b1 = 0.2/N
    # if (t>150) b1 = 0.2/N
    if (t>int_begin & t<int_end) b1 = int_level/N
    
    dS.dt = -(b1*I1+b2*I2+b3*I3)*S + c*R
    dE.dt=(b1*I1+b2*I2+b3*I3)*S-a*E
    dI1.dt=a*E-g1*I1-p1*I1
    dI2.dt=p1*I1-g2*I2-p2*I2
    dI3.dt=p2*I2-g3*I3-u*I3
    dR.dt=g1*I1+g2*I2+g3*I3 - c*R
    dD.dt=u*I3
    
    return(list(c(dS.dt, dE.dt, dI1.dt, dI2.dt, dI3.dt, dR.dt, dD.dt)))
  })
}

c_vals = c(30, # 1 month
           60, # 2 months
           120, # 4 months
           180, # 6 months
           365, # 1 year
           365000 # 100 years (i.e. permanent immunity)
)
I3.immunityloss = data.frame(matrix(nrow=1001,ncol=length(c_vals)+1))
names(I3.immunityloss) = c('time','1month','2months','4months','6months','12months','Permanent')
I3.immunityloss$time = 1:1001
for (i in 1:length(c_vals)) {
  int_begin = 150
  int_end = 1000
  int_level = 0.2
  c = 1/c_vals[i]
  p_new = c(p,int_begin,int_end,int_level,c)
  sol.temp = ode(y = y0, times = t, func = SetODEs_SEIR_ImmunityLoss, parms = p_new)
  sol.temp = as.data.frame(sol.temp)
  names(sol.temp) = c('time','S','E','I1','I2','I3','R','D')
  I3.immunityloss[,i+1] = sol.temp$I3
}
I3.immunityloss_long = melt(I3.immunityloss, id="time")
names(I3.immunityloss_long) = c('Day','Length_of_Immunity','CriticalCases')
ggplot(data = I3.immunityloss_long,
       aes(x = Day, y = CriticalCases, color=Length_of_Immunity)) + 
  geom_line() + ggtitle('Loss of Immunity')
