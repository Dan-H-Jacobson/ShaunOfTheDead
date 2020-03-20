library('deSolve')
library('ggplot2')
library('reshape2')

sizr <- function(t,y,p) {
  
  S = y[1] # Susceptible
  I = y[2] # Infected
  Z = y[3] # Zombies
  R = y[4] # Removed (Dead Zombies)
  
  with(as.list(p),{
    
    if(t>=6  ) b2= 5/7.44e6 # Initially, the average human killed 5 zombies/day
    if(t>=6.3) b2=10/7.44e6 # Later, the average human would kill 10 zombies/day
    
    dS = -b1*S*Z 
    dI = b1*S*Z - m*I
    dZ = m*I -b2*S*Z
    dR = b2*S*Z
    
    return(list(c(dS,dI,dZ,dR)))
    
  })
  
}

N = 7.44e6 # Population of London (2004)
b1 = 5/N # transmission rate
m = 3 # zombification rate = 1/(time spent infected) = 1/(8 hours)
b2 = 0 # rate of zombie killing
p = c(b1,m,b2)

S0 = N-1
I0 = 1 # begin with a single infected individual
Z0 = 0
R0 = 0
y0 = c(S0,I0,Z0,R0)

t = seq(from = 0, to = 7, by = .1) # One week split into 2.4hour parts
sol = ode(y = y0, times = t, func = sizr, parms = p)
sol.df = as.data.frame(sol)
names(sol.df) = c('time','Susceptible','Infected','Zombie','DeadZombie')

sol.df_long = melt(sol.df, id="time")
names(sol.df_long) = c('Day','State','Population')
ggplot(data = sol.df_long,
       aes(x = Day, y = Population, color=State)) + 
  geom_line() + xlim(4,7) + ggtitle('Shaun of the Dead')

# Final Statistics
print('=== Final Breakdown ===')
print(paste('Susceptible:',round(sol.df$Susceptible[nrow(sol.df)])))
print(paste('Infected:',round(sol.df$Infected[nrow(sol.df)])))
print(paste('Zombies:',round(sol.df$Zombie[nrow(sol.df)])))
print(paste('DeadZombies:',round(sol.df$DeadZombie[nrow(sol.df)])))




