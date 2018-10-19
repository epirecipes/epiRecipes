---
title: 'Host-vector models'
permalink: 'chapters/hostvector'
previouschapter:
  url: chapters/erlang/r
  title: 'R'
nextchapter:
  url: chapters/1host1vector/intro
  title: 'One host, one vector'
redirect_from:
  - 'chapters/hostvector'
---
## Host-vector models

*Author*: Carl A. B. Pearson @pearsonca
*Date*: 2018-10-02

### Overview

Vector-borne infections like malaria, dengue, yellow fever, and so on result in substantial health burdens, particularly in lower resource settings and associated lower quality infrastructure (e.g., un-screened windows, ad hoc water collection and storage) as well as higher rates of outdoor work.

The simplest ODE models of these systems represent just the natural history of infection in the host and vector species, along with other typical maximally simplifying assumptions (e.g., heterogeneous mixing, constant populations).

This demonstration starts with a single host and single vector, then progresses to multiple hosts and a single vector, and finally to a multiple-host-multiple-vector system.  With each additional step, the code uses more advanced Julia features.