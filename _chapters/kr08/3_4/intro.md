---
title: 'Program 3.4: Age structured SEIR'
permalink: 'chapters/kr08/3_4/intro'
previouschapter:
  url: chapters/kr08/3_2/r_desolve
  title: 'R using deSolve'
nextchapter:
  url: chapters/kr08/3_4/r_desolve
  title: 'R using deSolve'
redirect_from:
  - 'chapters/kr08/3-4/intro'
---

## Program 3.4: Age structured SEIR

Program 3.4 implement an SEIR model with four age-classes and yearly aging, closely matching the implications of grouping individuals into school cohorts. The four age-classes modelled are 0-6, 6-10, 10-20 and 20+ years old.
Key to this model are two basic assumptions:

1. Only individuals in the adult class give birth, and only individuals in the adult class die.

2. Births and deaths are continuous, but aging only happens once per year.
