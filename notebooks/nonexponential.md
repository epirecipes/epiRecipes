## Nonexponential passage times

Author: Simon Frost @sdwfrost

Date: 2018-10-19

Simple models (see e.g. [the standard SIR model](http://epirecip.es/epicookbook/chapters/sir/intro)) assume that individuals transition from infectious to recovered at a constant rate, such that the distribution of infectious periods follows [an exponential distribution](https://en.wikipedia.org/wiki/Exponential_distribution). For many infectious, this is an unrealistic assumption, and the distribution of infectious periods has a smaller variance than the exponential.

A number of approaches can be used to incorporate non-exponential distributions in deterministic epidemic models.

- The method of stages
- Partial differential equation models
- Integral equation models

For stochastic models, the method of stages can also be used. Alternatively, discrete-event simulation offers a simpler approach to implement these models, although at greater computational expense.
