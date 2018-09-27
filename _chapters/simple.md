---
title: 'Simple deterministic models'
permalink: 'chapters/simple'
previouschapter:
  url: 
  title: ''
nextchapter:
  url: chapters/sir/intro
  title: 'SIR model'
redirect_from:
  - 'chapters/simple'
---
## Simple deterministic models

Many epidemiological models are compartmental, meaning that the population is divided up into several discrete subpopulations, with flows between the subpopulations modeled using ordinary differential equations (ODEs). In most cases, these models are nonlinear, and numerical techniques are required to solve the dynamics over time. or a given set of parameter values, the model will generate the same output, hence these models are known as deterministic. These models generate output - e.g. the number of infected individuals over time - that are continuous in time and in number, i.e. it is possible to have 1.1 infected individuals, and so on. We can think of this as approximating the average of many stochastic epidemics over time - more on this in subsequent sections.

Algorithms for solving ODEs are implemented in many computer languages, making it relatively straightforward to simulate these models.
