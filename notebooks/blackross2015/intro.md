
## Exact recursive expressions for final size

Author: Sangeeta Bhatia @sangeetabhatia03
Date: 2018-10-02

Here we implement a model for estimating the final size distributions for an outbreak. The estimation is based on a simple SIR model but the method can be generalised to other models.

The details are available in [Black and Ross (2015)](https://doi.org/10.1016/j.jtbi.2014.11.029). Instead of keeping tab of the population counts for each compartment ($S$, $I$, $R$), we keep a record of the number of infection and recovery events. At time t, let the number of infection events be $Z_1$ and the number of recovery events be $Z_t$. Then the state of our model is $(Z_1, Z_2)$. From one state to another, we can have at most one  infection event or a single recovery event. Thus from $(Z_1, Z_2)$, the system can transition to $(Z_1 + 1, Z_2)$, $(Z_1, Z_2 + 1)$. To estimate the distribution of the final size, one needs to calculate the probability of 
all paths leading up to the state from which the system cannot transition $(Z, Z)$. The clever bit in the paper is the way the states are indexed. This indexing allows the probabilities to be calculated easily.


### References

[Black and Ross (2015)](https://doi.org/10.1016/j.jtbi.2014.11.029).
