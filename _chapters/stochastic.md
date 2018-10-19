---
title: 'Simple stochastic models'
permalink: 'chapters/stochastic'
previouschapter:
  url: chapters/deleo1996/r_desolve
  title: 'R using deSolve'
nextchapter:
  url: chapters/sir-stochastic-discretestate-continuoustime/intro
  title: 'Continuous time SIR'
redirect_from:
  - 'chapters/stochastic'
---
## Simple stochastic models

Stochastic, or chance, effects are important when population sizes are small, and are particularly strong when the underlying rates (such as infection and recovery rates) are high. Most epidemics start with a small number of infected individuals, and oscillatory epidemics may frequently be associated with low numbers of infectives during 'troughs'. We are often interested in small outbreaks, where stochastic effects may play a major role.

There are a number of ways of modeling stochastic epidemics. The most widely used algorithm is the [Doob-Gillespie](https://en.wikipedia.org/wiki/Gillespie_algorithm) algorithm, which generates exact trajectories of the system when the underlying rates are constant. Alternative approaches such as [uniformization](https://en.wikipedia.org/wiki/Uniformization_(probability_theory)) can be used when rates vary over time. It is often computationally advantageous to approximate the system either by considering discrete time steps (see e.g. [tau leaping](https://en.wikipedia.org/wiki/Tau-leaping)), and/or approximating the discrete number of individuals by a continuous value through the use of a [stochastic differential equation](https://en.wikipedia.org/wiki/Euler%E2%80%93Maruyama_method).
