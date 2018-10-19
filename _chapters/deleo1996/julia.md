---
interact_link: notebooks/deleo1996/julia.ipynb
title: 'Julia'
permalink: 'chapters/deleo1996/julia'
previouschapter:
  url: chapters/deleo1996/intro
  title: 'Scaling model'
nextchapter:
  url: chapters/deleo1996/python
  title: 'Python'
redirect_from:
  - 'chapters/deleo1996/julia'
---

### De Leo et al. scaling model in Julia

*Author*: Christopher Davis

*Date*: 2018-10-02


{:.input_area}
```julia
using DifferentialEquations
```


{:.input_area}
```julia
micro_1 = @ode_def Micro1 begin
    dS = r*(1-S/K)*S - β*S*I
    dI = β*S*I-(μ+α)*I
    end β r μ K α
```




{:.output_data_text}
```
(::Micro1) (generic function with 4 methods)
```




{:.input_area}
```julia
w = 1;
m = 10;
β = 0.0247*m*w^0.44;
r = 0.6*w^-0.27;
μ = 0.4*w^-0.26;
K = 16.2*w^-0.7;
α = (m-1)*μ;
```


{:.input_area}
```julia
parms = [β,r,μ,K,α];
init = [K,1.];
tspan = (0.0,10.0);
```


{:.input_area}
```julia
sir_prob = ODEProblem(micro_1,init,tspan,parms)
```




{:.output_data_text}
```
DiffEqBase.ODEProblem with uType Array{Float64,1} and tType Float64. In-place: true
timespan: (0.0, 10.0)
u0: [16.2, 1.0]
```




{:.input_area}
```julia
sir_sol = solve(sir_prob);
```

#### Visualisation


{:.input_area}
```julia
using Plots
```


{:.input_area}
```julia
plot(sir_sol,xlabel="Time",ylabel="Number")
```




![svg](../../images/chapters/deleo1996/julia_10_0.svg)



#### Threshold criterion for transmission rate


{:.input_area}
```julia
m = [5,10,20,40]
ws = 10.^linspace(-3,3,601)
βs = zeros(601,4)
for i = 1:4
    βs[:,i] = 0.0247*m[i]*ws.^0.44
end
plot(ws,βs,xlabel="Weight",ylabel="\\beta_min", xscale=:log10,yscale=:log10, label=["m = 5" "m = 10" "m = 20" "m = 40"],lw=3)
```




![svg](../../images/chapters/deleo1996/julia_12_0.svg)


