# ShaunOfTheDead
Scripts containing compartmental models describing a) the zombie outbreak in the 2004 film Shaun of the Dead, and b) various scenarios of the 2020 Covid-19 pandemic


# ShaunOfTheDead.R

This script contains a single SIZR model.
Initial values: S0 = population of London in 2004 (~7.4million) - 1; I0 = 1
Parameters: infection rate = 5/zombie.day, zombification rate = 3/day, zombie death rate = 0 (init) then increasing

# Covid_19.R

This script contains three SEIR models:
  SetODEs_SEIR: standard function created by Dr Alison Lynn Hill (https://github.com/alsnhll/SEIR_COVID19/blob/master/COVID19seir/server.R)
  SetODEs_SEIR_Intervention: includes changes to infection rate, with specified start and end days, and ability to set level of intervention
  SetODEs_SEIR_ImmunityLoss: includes loss of immunity dynamic, in which removed individuals re-enter the susceptible group.

Initial values: S0 = 10million
Parameters: default parameters set in the Hill model (see above).
