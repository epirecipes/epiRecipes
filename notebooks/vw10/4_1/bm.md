## Example 4.1 - original Berkeley Madonna code


```
{================================================================================================}
{                            SEIR model describing the transmission dynamics of an immunizing infection in a closed population                    }
{================================================================================================}


METHOD RK4

STARTTIME = 0
STOPTIME = 100 
DT = 0.02

{Differential equations}
d/dt (Susceptible) = -beta*Susceptible*Infectious 
d/dt (Preinfectious) =  beta*Susceptible*Infectious - Preinfectious*infectious_rate 
d/dt (Infectious) = Preinfectious*infectious_rate - Infectious*rec_rate 
d/dt (Immune) =  Infectious*rec_rate 
d/dt (Cum_infectious) = Preinfectious*infectious_rate

{Initial values}
      INIT Susceptible = Sus_0
      INIT Preinfectious = Preinfectious_0
      INIT Infectious = Infectious_0
      INIT Immune = Immune_0
      INIT Cum_infectious = Infectious_0


{ ===========================================================================================}
{                                                         INITIAL VALUES                                                                                                      }
{  These are assigned in INIT=box  for the corresponding compartment in the flowchart window (flowchart version only)                }
{ ===========================================================================================}

Sus_0 = total_popn-Infectious_0 -Immune_0                                                ; initial number of susceptible individuals
Preinfectious_0 = 0                                                                                   ; initial number of pre-infectious individuals
Infectious_0 = 1                                                                                        ; initial number of infectious individuals
Immune_0 = total_popn*prop_immune_0                                                     ; initial number of immune individuals



{ ===========================================================================================}
{                                                         TRANSMISSION AND INFECTION-RELATED PARAMETERS                                     }
{                              Note that these are in daily units, unless otherwise specified                                                                  }
{ Once the parameter has been set up, its value should be altered in the Parameters window or using the sliders,                      }
{ since changes to the parameter values in this window will not necessarily be recognized by Berkeley Madonna                       }
{ ===========================================================================================}

preinfectious_period = 8                                                                                  ; average pre-infectious period
infectious_period = 7                                                                                      ;  average infectious period
R0 = 13                                                                                                         ; R0 
infectious_rate = 1/preinfectious_period                                                            ; rate at which individuals become infectious
rec_rate = 1/infectious_period                                                                          ; rate at which individuals recover to become immune
beta = R0/(total_popn*infectious_period)                                                           ; rate at which two specific individuals come into effective contact per unit time
prop_immune_0 = 0.0                                                                                     ; proportion immune at the start


{ ===========================================================================================}
{                                                         DEMOGRAPHY-RELATED  PARAMETERS                                                              }
{                              Note that these are in daily units, unless otherwise specified                                                                  }
{ Once the parameter has been set up, its value should be altered in the Parameters window or using the sliders,                      }
{ since changes to the parameter values in this window won't necessarily be recognized by Berkeley Madonna                          }
{ ===========================================================================================}

total_popn = 100000                                                                                     ; total population size


{ ===========================================================================================}
{                                                         USEFUL STATISTICS                                                                                              }
{ ===========================================================================================}

new_infns = beta*Susceptible*Infectious                                                        ; number of new infections per unit time
new_infectious = Preinfectious*infectious_rate                                                ; number of new infectious persons per unit time
prop_sus = Susceptible/total_popn                                                                ; proportion of the population that is susceptible
Rn = R0*prop_sus                                                                                        ; the net reproduction number
prop_imm = 1-prop_sus                                                                                ; proportion of the population that is immune
sum_pop = susceptible+preinfectious+infectious+immune                               ; sum of the number of individuals in the population
ln_new_infectious = if (time>0) then logn(new_infectious)  else 0                      ; natural log of the no. of new infectious persons per day
ln_infectious = logn(infectious)                                                                       ; natural log of the no. of infectious persons
ln_cum_infectious = logn(cum_infectious)                                                       ; natural log of the cumulative no. of infectious persons

```