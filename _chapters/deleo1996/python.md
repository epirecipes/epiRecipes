---
interact_link: notebooks/deleo1996/python.ipynb
title: 'Python'
permalink: 'chapters/deleo1996/python'
previouschapter:
  url: chapters/deleo1996/julia
  title: 'Julia'
nextchapter:
  url: chapters/deleo1996/r_desolve
  title: 'R using deSolve'
redirect_from:
  - 'chapters/deleo1996/python'
---

### Python using SciPy

*Author*: Christopher Davis

*Date*: 2018-10-02


{:.input_area}
```python
import numpy as np
from scipy.integrate import ode, solve_ivp
```


{:.input_area}
```python
def micro_1(times,init,parms):
    beta, r, mu, K, alpha = parms
    S,I = init
    # ODEs
    dS = r*(1-S/K)*S - beta*S*I
    dI = beta*S*I-(mu + alpha)*I
    return [dS,dI]
```


{:.input_area}
```python
w = 1
m = 10
beta = 0.0247*m*w**0.44
r = 0.6*w**-0.27
mu = 0.4*w**-0.26
K = 16.2*w**-0.7
alpha = (m-1)*mu
```


{:.input_area}
```python
parms = [beta,r,mu,K,alpha]
init = [K,1.]
times = np.linspace(0,10,101)
```


{:.input_area}
```python
sir_sol = solve_ivp(fun=lambda t, y: micro_1(t, y, parms), t_span=[min(times),max(times)], y0=init, t_eval=times)
```

#### Visualisation


{:.input_area}
```python
import matplotlib.pyplot as plt
```


{:.input_area}
```python
plt.plot(sir_sol.t,sir_sol.y[0],color="red",linewidth=2, label = "S(t)")
plt.plot(sir_sol.t,sir_sol.y[1],color="blue",linewidth=2, label = "I(t)")
plt.xlabel("Time")
plt.ylabel("Number")
plt.legend()
```




{:.output_data_text}
```
<matplotlib.legend.Legend at 0x7f776a919198>
```




![png](../../images/chapters/deleo1996/python_9_1.png)


#### Threshold criterion for transmission rate


{:.input_area}
```python
m = [5,10,20,40]
ws = 10**np.linspace(-3,3,601)
betas = np.zeros((601,4))
for i in range(4):
    betas[:,i] = 0.0247*m[i]*ws**0.44
plt.loglog(ws,betas)
plt.xlabel("Weight")
plt.ylabel(r'$\beta_{min}$')
plt.legend(("m = 5", "m = 10", "m = 20", "m = 40"))
```




{:.output_data_text}
```
<matplotlib.legend.Legend at 0x7f7769b707f0>
```




![png](../../images/chapters/deleo1996/python_11_1.png)

