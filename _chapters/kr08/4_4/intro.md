---
title: 'Program 4.4: Multi-host SEIR'
permalink: 'chapters/kr08/4_4/intro'
previouschapter:
  url: chapters/kr08/3_4/r_desolve
  title: 'R using deSolve'
nextchapter:
  url: chapters/kr08/4_4/r_desolve
  title: 'R using deSolve'
redirect_from:
  - 'chapters/kr08/4-4/intro'
---

## Program 4.4: Multi-host SIR model

This is an example of a multi-host SIR model where there are two hosts (humans and mosquitoes). Transmissions within human or mosquito populations is also possible. It is adapted from a MATLAB code example in Keeling & Rohani.

$\frac{dX_H}{dt}=\nu_H-r(T_{HM}Y_M+T_{HH}Y_H)X_H-\mu_HX_H$

$\frac{dX_M}{dt}=\nu_M-r(T_{MH}Y_H+T_{MM}Y_M)X_M-\mu_MX_M$

$\frac{dY_H}{dt}=r(T_{HM}Y_M+T_{HH}Y_H)X_H-\mu_HY_H-\gamma_HY_H$

$\frac{dY_M}{dt}=r(T_{MH}Y_H+T_{MM}Y_M)X_M-\mu_MY_M-\gamma_MY_M$

where:

$X_H$ & $X_M$ are the human and mosquito susceptible populations and $Y_H$ & $Y_M$ are the infected ones.
$\nu_H$ & $\nu_M$ are the birth rates of humans and mosquitoes.
$r$ is the mosquito biting rate of humans.
$T$ is a 2x2 matrix of transmission probabilities between and within human($H$) and mosquito($M$) populations.
$\mu_H$ & $\mu_M$ are the per capita death rate of humans and mosquitoes.
$\gamma_H$ & $\gamma_M$ are the recovery rates of humans and mosquitoes.


### References

- [Program 4.4 from Keeling & Rohani (2007)](http://homepages.warwick.ac.uk/~masfz/ModelingInfectiousDiseases/Chapter4/Program_4.4/Program_4_4.m)
