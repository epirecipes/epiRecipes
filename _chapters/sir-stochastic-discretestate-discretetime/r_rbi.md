---
interact_link: notebooks/sir-stochastic-discretestate-discretetime/r_rbi.ipynb
title: 'R using LibBi'
permalink: 'chapters/sir-stochastic-discretestate-discretetime/r_rbi'
previouschapter:
  url: chapters/sir-stochastic-discretestate-discretetime/r_pomp
  title: 'R using POMP'
nextchapter:
  url: chapters/sir-stochastic-discretestate-discretetime/python
  title: 'Python'
redirect_from:
  - 'chapters/sir-stochastic-discretestate-discretetime/r-rbi'
---


{:.input_area}
```R
library(rbi)
```


{:.input_area}
```R
sir.bi <- "
model sir {
  const N = 403
  const timestep = 0.1
  state S, I, R, Y
  param bet, gamm, iota
  noise infection, recovery
  sub parameter {
    bet <- 2.62110617498984
    gamm <- 0.5384615384615384
    iota <- 0.5
  }
  sub initial {
    S <- 60
    I <- 1
    R <- 342
    Y <- 0    
  }
  sub transition (delta=timestep) {
    inline lambd = bet*(I+iota)/N
    inline ifrac = 1.0-exp(-lambd*timestep)
    inline rfrac = 1.0-exp(-gamm*timestep)
    infection ~ binomial(S,ifrac)
    recovery ~ binomial(I,rfrac)
    S <- S - infection
    I <- I + infection - recovery
    R <- R + recovery
    Y <- Y + infection
  }
}
"
```


{:.input_area}
```R
sir.bi.model <- bi_model(lines=sir.bi)
sir.bi.obj <- libbi(model=sir.bi.model)
sir.bi.sample <- sample(sir.bi.obj,
                        target="joint",
                        end_time=54.0,
                        nsamples=1000,
                        seed=1)
sir.bi.results <- bi_read(sir.bi.sample$output_file_name)
```


{:.input_area}
```R
sir.bi.results
```
