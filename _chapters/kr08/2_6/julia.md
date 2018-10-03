---
interact_link: notebooks/kr08/2_6/julia.ipynb
title: 'Julia'
permalink: 'chapters/kr08/2_6/julia'
previouschapter:
  url: chapters/kr08/2_6/r_desolve
  title: 'R using deSolve'
nextchapter:
  url: chapters/kr08/3_1/intro
  title: 'Program 3.1: SIS with risk groups'
redirect_from:
  - 'chapters/kr08/2-6/julia'
---

### Program 2.6: SEIR model in Julia using DifferentialEquations

Author: Lloyd Chapman @LloydChapman

Date: 2018-10-01


{:.input_area}
```julia
using DifferentialEquations
```


{:.input_area}
```julia
function seir_ode(dY,Y,p,t)
 dY[1] = p[4]-p[1]*Y[1]*Y[3]-p[4]*Y[1]
 dY[2] = p[1]*Y[1]*Y[3]-(p[2]+p[4])*Y[2]
 dY[3] = p[2]*Y[2] - (p[3]+p[4])*Y[3]
end
```




{:.output_data_text}
```
seir_ode (generic function with 1 method)
```




{:.input_area}
```julia
par=[520/365,1/60,1/30,774835/(65640000*365)]
init=[0.8,0.1,0.1]
tspan=(0.0,365.0)
```




{:.output_data_text}
```
(0.0, 365.0)
```




{:.input_area}
```julia
seir_prob = ODEProblem(seir_ode,init,tspan,par)
```




{:.output_data_text}
```
DiffEqBase.ODEProblem with uType Array{Float64,1} and tType Float64. In-place: true
timespan: (0.0, 365.0)
u0: [0.8, 0.1, 0.1]
```




{:.input_area}
```julia
sol=solve(seir_prob);
```


{:.input_area}
```julia
using Plots
```


{:.input_area}
```julia
R=ones(1,size(sol,2))-sum(sol,1);
```


{:.input_area}
```julia
plot(sol.t,[sol',R'],xlabel="Time",ylabel="Proportion")
```




![svg](../../../images/chapters/kr08/2_6/julia_9_0.svg)


