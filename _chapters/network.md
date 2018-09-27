---
title: 'Network models'
permalink: 'chapters/network'
previouschapter:
  url: chapters/sirforced/r
  title: 'R'
nextchapter:
  url: chapters/sircn/intro
  title: 'An edge based SIR model on a configuration network'
redirect_from:
  - 'chapters/network'
---
## Models on networks

Most models assume that populations are well-mixed, or consider multiple sub-populations in which individuals within each sub-population are assumed to mix randomly with other members of the same sub-population. While increasing the number of sub-populations can accommodate heterogeneity in contacts to a certain degree, it may be more biologically realistic to assume that individuals contacts are structured as a network, which may be static or may change over time.

One approach to modeling the spread of infection over a network is to use detailed simulations. This approach is used in the R package [`EpiModel`](http://www.epimodel.org/). Alternatively, compartmental models that approximate the dynamics on a network can be used, such as [pair approximations](http://doi.org/10.1098/rspb.1997.0159) or [edge-based models](http://doi.org/10.1098/rsif.2011.0403).
