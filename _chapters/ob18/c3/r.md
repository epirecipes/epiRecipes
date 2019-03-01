---
interact_link: notebooks/ob18/c3/r.ipynb
title: 'Original R code'
permalink: 'chapters/ob18/c3/r'
previouschapter:
  url: chapters/ob18/c3/intro
  title: 'Chapter 3: R0'
nextchapter:
  url: chapters/ob18/c4/intro
  title: 'Chapter 4: FoI and age-dependent incidence'
redirect_from:
  - 'chapters/ob18/c3/r'
---

## Chapter 3 from Bjornstad (2018): $R_0$

Let us explore using weekly measles-data from the 2003 outbreak in Niamey, Niger . The data is available as `niamey` in the `epimdr`-package. The `tot_cases`-column represents the total incidence across the city for each week of the outbreak.


{:.input_area}
```R
library(epimdr)
```


{:.input_area}
```R
data(niamey)
head(niamey[, 1:5])
```


<div markdown="0">
<table>
<thead><tr><th scope=col>absweek</th><th scope=col>week</th><th scope=col>tot_cases</th><th scope=col>tot_mort</th><th scope=col>lethality</th></tr></thead>
<tbody>
	<tr><td>1       </td><td>45      </td><td>11      </td><td>0       </td><td>0.000000</td></tr>
	<tr><td>2       </td><td>46      </td><td>12      </td><td>1       </td><td>8.333333</td></tr>
	<tr><td>3       </td><td>47      </td><td>15      </td><td>0       </td><td>0.000000</td></tr>
	<tr><td>4       </td><td>48      </td><td>14      </td><td>1       </td><td>7.142857</td></tr>
	<tr><td>5       </td><td>49      </td><td>30      </td><td>0       </td><td>0.000000</td></tr>
	<tr><td>6       </td><td>50      </td><td>41      </td><td>1       </td><td>2.439024</td></tr>
</tbody>
</table>

</div>


We can do a visual inspection to identify the initial period of exponential growth:


{:.input_area}
```R
par(mar = c(5,5,2,5))
plot(niamey$absweek, niamey$tot_cases, type = "b", 
     xlab = "Week", ylab = "Incidence")
par(new = TRUE)
plot(niamey$absweek, niamey$cum_cases, type = "l", 
     col = "red", axes = FALSE, xlab = NA, ylab = NA, log = "y")
axis(side = 4)
mtext(side = 4, line = 4, "Cumulative incidence")
legend("topleft", legend = c("Cases", "Cumulative"),
       lty = c(1,1), pch = c(1,NA), col = c("black", "red"))
```


![png](../../../images/chapters/ob18/c3/r_4_0.png)


The cumulative incidence looks pretty log-linear for the first 6 weeks or so (Fig. \[fig:niameyweekly\]). The data is weekly and the serial interval for measles is around 10-12 days, thus V is around 1.5-1.8 weeks; We calculate $R_0$ assuming either $1.5$ or $1.8$:


{:.input_area}
```R
fit=lm(log(cum_cases)~absweek, subset=absweek<7, 
     data=niamey)
r=fit$coef["absweek"]
V=c(1.5, 1.8)
V*r+1
```


<div markdown="0">
<ol class=list-inline>
	<li>1.69423334567417</li>
	<li>1.833080014809</li>
</ol>

</div>


So a fast-and-furious estimate of the reproductive ratio for this outbreak places it in the 1.5-2 range. Measles exhibits recurrent epidemics in the presence of various vaccination campaigns in Niger, so this number represents an estimate of the *effective* reproductive ratio, $R_E$, at the beginning of this epidemic.

In their analysis of the SARS epidemics, showed that for an infection with distinct latent and infectious periods a more refined estimate is given by $R = V r + 1 + f (1-f) (V r)^2$, where f is the ratio of infectious period to serial interval. For measles the infectious period is around 5 days:


{:.input_area}
```R
V = c(1.5, 1.8)
f = (5/7)/V
V * r + 1 + f * (1 - f) * (V * r)^2
```


<div markdown="0">
<ol class=list-inline>
	<li>1.81445011031601</li>
	<li>1.99919772595046</li>
</ol>

</div>


### Maximum likelihood: the chain-binomial model

We employ a standard recipe, for doing a ‘non-standard’ maximum likelihood analysis . The first step is to write a function for the likelihood. Conditional on some parameters, the function returns the negative log-likelihood of observing the data given the model. The likelihood, which is the probability of observing data given a model and some parameter values, is the working-horse of a large part of statistics. R has inbuilt `dxxxx`-functions to calculate the likelihood for any conceivable probability distribution. The function to calculate a binomial likelihood is `dbinom`. We can thus define a likelihood-function for the chain-binomial model [7]:


{:.input_area}
```R
llik.cb = function(S0, beta, I) {
    n = length(I)
    S = floor(S0 - cumsum(I[-n]))
    p = 1 - exp(-beta * (I[-n])/S0)
    L = -sum(dbinom(I[-1], S, p, log = TRUE))
    return(L)
}
```

For the the real statistical analysis (below), the two parameters will be estimated simultaneously. However, in order to ease into the idea of likelihood estimation we will consider the two sequentially and visualize the likelihood by plotting it over a grid of potential values. We illustrate with the data on measles from one of the three different reporting centers in Niamey, Niger from 2003 . We first need to aggregate the data into 2-week intervals which is roughly the serial interval for measles. The epidemic in district 1 lasted for 30 weeks (the 31st week is a zero)[8]:


{:.input_area}
```R
twoweek = rep(1:15, each = 2)
y = sapply(split(niamey$cases_1[1:30], twoweek), sum)
sum(y)
```


<div markdown="0">
5920
</div>


In district 1 there were 5920 cases during the epidemics, so $S_0$ needs to be at least that number. In the above parameterization $R_E \simeq \beta$, lets initially assume a candidate value of 6500 for $S_0$ and calculate the likelihood for each candidate value of $\beta$ between 1 and 10 by 0.1 (fig. \[fig:profbeta\]):


{:.input_area}
```R
S0cand = 6500       
betacand = seq(0,10, by = .1)
ll = rep(NA, length(betacand))
for(i in 1:length(betacand)){
     ll[i] = llik.cb(S0 = S0cand, beta = betacand[i],
     I = y)}
plot(ll ~ betacand, ylab = "Neg log-lik", 
     xlab = expression(beta)) 
betacand[which.min(ll)]
```


<div markdown="0">
2.3
</div>



![png](../../../images/chapters/ob18/c3/r_14_1.png)


We follow the convention of using the negative log-likelihood in the profile. Intuitively, one may think that it would be more natural to consider the *likelihood* itself (the probability of observing the data, given particular parameter values). However, since this would be a product of small numbers (one for each observation), computers are not precise enough to distinguish the joint probability from zero if the data set is large. If our $S_0$ guess is right then $\beta$ should be around 2.3. We can do a similar check for $S_0$ (assuming $\beta$ is 2.3). The grid-value associated with the highest likelihood value is 7084.8 (fig. \[fig:profRE\]), so our original $S_0$ guess was good but not perfect.


{:.input_area}
```R
betacand = 2.3
S0cand = seq(5920,8000, length = 101)
ll = rep(NA, length = 101)
for(i in 1:101){
     ll[i] = llik.cb(S0 = S0cand[i], beta = betacand, 
     I = y)}
plot(ll ~ S0cand, ylab = "Neg log-lik", 
     xlab = expression(S[0]))
S0cand[which.min(ll)]
```


<div markdown="0">
7084.8
</div>



![png](../../../images/chapters/ob18/c3/r_16_1.png)


For a proper analysis we minimize the negative log-likelihood by varying both parameters simultaneously. We can do this using the generic `optim`-function or the `mle2`-function in the `bbmle`-package. The `mle2`-function uses `optim` to find maximum likelihood estimates, but also provides confidence intervals, profile likelihoods and a variety of other useful measures . We summarize the basic pertinent likelihood theory for these other measures in section \[sec:c7lik\].


{:.input_area}
```R
library(bbmle)
fit = mle2(llik.cb, start = list(S0 = 7085, beta = 2.3), 
    method = "Nelder-Mead",data = list(I = y))
summary(fit)
confint(fit)
```


{:.output_data_text}
```
Maximum likelihood estimation

Call:
mle2(minuslogl = llik.cb, start = list(S0 = 7085, beta = 2.3), 
    method = "Nelder-Mead", data = list(I = y))

Coefficients:
       Estimate Std. Error z value     Pr(z)    
S0   7.8158e+03 1.3022e+02  60.019 < 2.2e-16 ***
beta 1.8931e+00 3.6968e-02  51.209 < 2.2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

-2 log L: 841.831 
```



<div markdown="0">
<table>
<thead><tr><th></th><th scope=col>2.5 %</th><th scope=col>97.5 %</th></tr></thead>
<tbody>
	<tr><th scope=row>S0</th><td>7577.967212</td><td>8088.641095</td></tr>
	<tr><th scope=row>beta</th><td>   1.820943</td><td>   1.966336</td></tr>
</tbody>
</table>

</div>


So the joint MLE estimates are $S_0 = 7816$ (CI: <span>7578, 8088</span>) and $\beta = 1.89$ (CI: {1.82, 1.97}).

Applying statistical tools to biological models – like the chain-binomial – can usefully highlight uncertainties due to parametric interdependencies. In the case of a ‘simple epidemic’ like the measles outbreak considered here, for example, it is conceivable that similar epidemic trajectories can arise from having a large number of initial susceptibles and a low transmission rate, or a more moderate number of susceptibles and a higher transmission rate. We can quantify this through considering the correlation matrix among the parameters of our likelihood analysis; `vcov` calculates their variance-covariance matrix from which we can calculate standard errors according to `sqrt(diag(vcov(fit)))` and `cov2cor` converts this to a correlation matrix. As intuition suggested there is a strong negative correlation between the estimates of the $\beta$ and $S_0$ parameters.


{:.input_area}
```R
cov2cor(vcov(fit))
```


<div markdown="0">
<table>
<thead><tr><th></th><th scope=col>S0</th><th scope=col>beta</th></tr></thead>
<tbody>
	<tr><th scope=row>S0</th><td> 1.0000000</td><td>-0.7444261</td></tr>
	<tr><th scope=row>beta</th><td>-0.7444261</td><td> 1.0000000</td></tr>
</tbody>
</table>

</div>


## Stochastic simulation

The chain-binomial is both a statistical model for estimation and a stochastic model for dynamics. We can thus write a function to simulate dynamics using the estimated parameters.


{:.input_area}
```R
sim.cb=function(S0, beta, I0){
   I=I0
   S=S0
   i=1
   while(!any(I==0)){
        i=i+1
        I[i]=rbinom(1, size=S[i-1], prob=1-
           exp(-beta*I[i-1]/S0))
        S[i]=S[i-1]-I[i]
   }
   out=data.frame(S=S, I=I)
   return(out)
}
```

We superimpose 100 stochastic simulations on the observed epidemic. The simulations from the chain-binomial model brackets the observed epidemic nicely (fig. \[fig:cbsim\]), suggesting that the model is a reasonable first approximation to the underlying dynamics. We will revisit on this case study in the context of outbreak-response vaccination in section \[sec:orv\].


{:.input_area}
```R
plot(y, type="n", xlim=c(1,18), 
   ylab="Predicted/observed", xlab="Week")
for(i in 1:100){
     sim=sim.cb(S0=floor(coef(fit)["S0"]), 
     beta=coef(fit)["beta"], I0=11)
     lines(sim$I, col=grey(.5))
}
points(y, type="b", col=2)
```


![png](../../../images/chapters/ob18/c3/r_24_0.png)


Further examples
----------------

### Influenza A/H1N1 1977

The `flu` dataset in the `epimdr`-package represents the number of children confined to bed each day during a 1978 outbreak of the reemerging influenza A/H1N1 strain in a boarding school in North England (fig. \[fig:flu\]). This subtype of influenza had been absent from human circulation after the A/H2N2 pandemic of 1957 but reemerge (presumably from some laboratory freezer) in 1977. The school had 763 boys of which 512 boys were confined to bed sometime during the outbreak. None of the boys would have had previous exposure to A/H1N1.

The typical time of illness was 5-7 days. Since the data is number confined to bed each day, the data is not incidence but (a proxy for) *prevalence*. The data looks pretty log-linear for the first 5 days. Family studies have been used to estimate the serial interval for flu between 2 and 4 days . Volunteer studies shows the mean infectious period around 5 days .


{:.input_area}
```R
data(flu)
plot(flu$day, flu$cases, type="b", xlab="Day", 
     ylab="In bed", log="y")
tail(flu)
```


<div markdown="0">
<table>
<thead><tr><th></th><th scope=col>day</th><th scope=col>cases</th></tr></thead>
<tbody>
	<tr><th scope=row>9</th><td> 9 </td><td>192</td></tr>
	<tr><th scope=row>10</th><td>10 </td><td>126</td></tr>
	<tr><th scope=row>11</th><td>11 </td><td> 70</td></tr>
	<tr><th scope=row>12</th><td>12 </td><td> 28</td></tr>
	<tr><th scope=row>13</th><td>13 </td><td> 12</td></tr>
	<tr><th scope=row>14</th><td>14 </td><td>  5</td></tr>
</tbody>
</table>

</div>



![png](../../../images/chapters/ob18/c3/r_26_1.png)


The “fast-and-furious" estimate of $R_0$ is thus:


{:.input_area}
```R
fit=lm(log(cases)~day, subset=day<=5, 
     data=flu)
r=fit$coef["day"]
V=c(2,3)
V*r+1
```


<div markdown="0">
<ol class=list-inline>
	<li>3.17188439700297</li>
	<li>4.25782659550445</li>
</ol>

</div>


This is higher than most estimates of $R_0$ of pandemic flu (which typically lies in the 1.5-2.5 interval). However, contact rates within a boarding school is likely to be higher than average across human populations as a whole.

### Ebola Sierra Leone 2014-15

The CDC’s record for the 2014-15 ebola outbreak in Sierra Leone is in the `ebola`-dataset. The serial interval for ebola is estimated at around 15 days with an incubation period of 11 days. The mean time to hospitalization is 5 days and mean time to death or dismissal was 5 and 11 days, respectively . The data is the back-calculated incidence as the difference of the cumulative cases reported by the [](http://www.cdc.gov/vhf/ebola/outbreaks/2014-west-africa/cumulative-cases-graphs.html). Because of the complexities of reporting and revisions of case-load through time, this lead to some negative numbers for certain dates. These were set to zero as a crude fix.


{:.input_area}
```R
data(ebola)
par(mar = c(5,5,2,5))
plot(ebola$day, ebola$cases, type="b", xlab="Week", 
     ylab="Incidence")
par(new=T)
plot(ebola$day, ebola$cum_cases, type="l", col="red",
      axes=FALSE, xlab=NA, ylab=NA, log="y")
axis(side = 4)
mtext(side = 4, line = 4, "Cumulative incidence")
legend("right", legend=c("Cases", "Cumulative"),
     lty=c(1,1), pch=c(1,NA), col=c("black", "red"))
tail(ebola)
```


<div markdown="0">
<table>
<thead><tr><th></th><th scope=col>date</th><th scope=col>day</th><th scope=col>cum_cases</th><th scope=col>cases</th></tr></thead>
<tbody>
	<tr><th scope=row>98</th><td>7/8/15 </td><td>468    </td><td>13945  </td><td>34     </td></tr>
	<tr><th scope=row>99</th><td>7/15/15</td><td>475    </td><td>13982  </td><td>37     </td></tr>
	<tr><th scope=row>100</th><td>7/22/15</td><td>482    </td><td>14001  </td><td>19     </td></tr>
	<tr><th scope=row>101</th><td>7/29/15</td><td>489    </td><td>14061  </td><td>60     </td></tr>
	<tr><th scope=row>102</th><td>8/5/15 </td><td>496    </td><td>14089  </td><td>28     </td></tr>
	<tr><th scope=row>103</th><td>8/12/15</td><td>503    </td><td>14122  </td><td>33     </td></tr>
</tbody>
</table>

</div>



![png](../../../images/chapters/ob18/c3/r_30_1.png)


We first use the regression method with Lipsitch’s correction:


{:.input_area}
```R
fit=lm(log(cum_cases)~day, subset=day<100, 
   data=ebola)
r=fit$coef["day"]
V=15
f=.5
V*r+1+f*(1-f)*(V*r)^2
```


<div markdown="0">
<strong>day:</strong> 1.69881135562533
</div>


We next aggregate the data in two-week increments roughly corresponding to the serial interval so we can apply the removal method.


{:.input_area}
```R
#Data aggregation
cases=sapply(split(ebola$cases, 
   floor((ebola$day-.1)/14)), sum)
sum(cases)

#Removal MLE
fit = mle2(llik.cb, start = list(S0 = 20000, beta = 2), 
   method = "Nelder-Mead",data = list(I = cases))
summary(fit)
```


<div markdown="0">
14721
</div>



{:.output_data_text}
```
Maximum likelihood estimation

Call:
mle2(minuslogl = llik.cb, start = list(S0 = 20000, beta = 2), 
    method = "Nelder-Mead", data = list(I = cases))

Coefficients:
       Estimate Std. Error    z value     Pr(z)    
S0   2.7731e+04 2.5949e-07 1.0687e+11 < 2.2e-16 ***
beta 1.4237e+00 1.1783e-02 1.2083e+02 < 2.2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

-2 log L: 5546.683 
```


Because of the difference in magnitude of the estimates of $S_0$ (in the ten thusands) and $R_0$ (around 1.4), the numerical method used to calculate confidence intervals struggles, so we suggest starting standard errors for the `confint`-function.


{:.input_area}
```R
confint(fit, std.err = c(100, 0.1))
```


<div markdown="0">
<table>
<thead><tr><th></th><th scope=col>2.5 %</th><th scope=col>97.5 %</th></tr></thead>
<tbody>
	<tr><th scope=row>S0</th><td>26393.579452</td><td>29287.725327</td></tr>
	<tr><th scope=row>beta</th><td>    1.384683</td><td>    1.463184</td></tr>
</tbody>
</table>

</div>


The removal and Lipsitch methods provides comparable estimates that are somewhat lower than those concluded by more elaborate analyses by the WHO team for the Sierra Leone outbreak .

### Ebola DRC 1995

| Disease     | Serial Interval | Location    | Year |
|:------------|:----------------|:------------|:----:|
| Ebola       | 14d             | DRC         | 1995 |
|             |                 | Uganda      | 2000 |
| SARS        | 5d              | Hong Kong   | 2003 |
|             |                 | Singapore   |      |
| Hog cholera | 7d              | Netherlands | 1997 |
| FMD         | 21d             | UK          | 2000 |

The `ferrari`-dataset holds the incidence data for a number of outbreaks – Ebola DRC ’95, Ebola Uganda ’00, SARS Hong Kong ’03, SARS Singapore ’03, Hog Cholera Netherlands ’97 and Foot-and-mouth UK ’00 – aggregated by disease-specific serial intervals . Let’s explore the 1995 DRC ebola outbreak.


{:.input_area}
```R
names(ferrari)
ferrari$Ebolacases95
sum(ferrari$Ebolacases95, na.rm = TRUE)
y = c(na.omit(ferrari$Ebolacases95))
```


<div markdown="0">
<ol class=list-inline>
	<li>'Eboladeaths00'</li>
	<li>'Ebolacases00'</li>
	<li>'Ebolacases95'</li>
	<li>'FMDfarms'</li>
	<li>'HogCholera'</li>
	<li>'SarsHk'</li>
	<li>'SarsSing'</li>
</ol>

</div>



<div markdown="0">
<ol class=list-inline>
	<li>4</li>
	<li>6</li>
	<li>5</li>
	<li>18</li>
	<li>36</li>
	<li>99</li>
	<li>40</li>
	<li>17</li>
	<li>4</li>
	<li>1</li>
	<li>&lt;NA&gt;</li>
	<li>&lt;NA&gt;</li>
	<li>&lt;NA&gt;</li>
	<li>&lt;NA&gt;</li>
	<li>&lt;NA&gt;</li>
</ol>

</div>



<div markdown="0">
230
</div>


The number of initial susceptibles must be larger than the summed incidence so we make an initial guess of 300.


{:.input_area}
```R
fit = mle2(llik.cb, method = "Nelder-Mead", 
    start = list(S0 = 300, beta = 2), 
    data = list(I = y))
fit
confint(fit, std.err = 2)
```


{:.output_data_text}
```

Call:
mle2(minuslogl = llik.cb, start = list(S0 = 300, beta = 2), method = "Nelder-Mead", 
    data = list(I = y))

Coefficients:
        S0       beta 
241.118108   3.181465 

Log-likelihood: -48.3 
```



<div markdown="0">
<table>
<thead><tr><th></th><th scope=col>2.5 %</th><th scope=col>97.5 %</th></tr></thead>
<tbody>
	<tr><th scope=row>S0</th><td>233.973778</td><td>254.051292</td></tr>
	<tr><th scope=row>beta</th><td>  2.692505</td><td>  3.718357</td></tr>
</tbody>
</table>

</div>


The estimated $R_0$ is 3.2. It thus appears that the Ebola outbreak in DRC in 1995 was more explosive than in Sierra Leone in 2014. This could be due to aggregation across a larger geographic area of the latter and/or the more intensive public health interventions. We will revisit on the DRC outbreak using the ‘next-generation matrix’ method.

### Contact tracing

The directional transmission graph of gonorrhoea among 89 individuals is in the `gonnet`-data set. The initial cluster of 17 cases all frequented the same bar, each infected between 0 and 7 other partners with 2.17 as the average. We can use the `statnet`-package to visualize the chains of transmission:


{:.input_area}
```R
require(statnet)
data(gonnet)
nwt = network(gonnet, directed = TRUE)
plot(nwt, vertex.col = c(0, rep(1, 17), rep(2, 71)))
```


![png](../../../images/chapters/ob18/c3/r_42_0.png)


The subsequent infections, in turn, infected between 0 and 6 partners with an average of 0.62. The drop is (i) due to the sexual network being depleted of suscpetibles, and (ii) because infection across heterogenous networks will differentially infect individuals according to their number of contacts . Epidemics across social networks is the topic of chapter \[chap:c11\] and we will revisit on this network therein.

### Advanced: The next-generation matrix

For epidemics that are not simple linear chains, it is less straight-forward to calculate $R_0$ from parameterized models using the ‘logical method’. The [](https://en.wikipedia.org/wiki/Next-generation_matrix) is the general approach that work for all compartmental models of any complexity . It is done in a sequence of steps:

1.  Identify all $n$ infected compartments,

2.  Construct a $n \times 1$ matrix, $\vec{F}$, that contains expressions for all *completely new* infections entering each infected compartment,

3.  Construct a $n \times 1$ matrix, $\vec{V^-}$, that contains expressions for all losses out of each infected compartment,

4.  Construct a $n \times 1$ matrix, $\vec{V^+}$, that contains expressions for all gains into each infected compartment that does *not* represent *new* infections but transfers among infectious classes,

5.  Construct a $n \times 1$ matrix, $\vec{V}= \vec{V^-} - \vec{V^+}$,

6.  Generate two $n \times n$ Jaccobian matrices, $\vec{f}$ and $\vec{v}$ that are the partial derivatives of $\vec{F}$ and $\vec{V}$ with respect to the $n$ infectious state variables,

7.  Evaluate the matrices at the disease free equilibrium (dfe) finally

8.  $R_0$ is the greatest eigenvalue of $\vec{f} \vec{v}^{-1}|_{dfe}$.

#### SEIR

This is quite an elaborate scheme, so we will try it out first for the SEIR model for which we already know the answer. Unfortunately, `R` can not do vectorized *symbolic* calculations, so we need to do this, one matrix element at a time. In chapter 1, we discussed how to use `expression` to do symbolic calculations in `R`. The `quote`-function is an alternative way to define mathematical expressions; `substitute` allows some simple additional manipulations.

Step 1: Infected classes are $E$ and $I$, let us label them 1 and 2.

Step 2: All new infections: $dE/dt = \beta S I / N$, $dI/dt =0$


{:.input_area}
```R
F1 = quote(beta * S * I/N)
F2 = 0
```

Step 3: All losses $dE/dt= (\mu+\sigma) E$, $dI/dt=(\mu + \alpha + \gamma) I$


{:.input_area}
```R
Vm1 = quote(mu * E + sigma * E)
Vm2 = quote(mu * I + alpha * I + gamma * I)
```

Step 4 : All gained transfers $dE/dt= 0$, $dI/dt=(\sigma) E$


{:.input_area}
```R
Vp1 = 0
Vp2 = quote(sigma * E)
```

Step 5: Subtract Vp from Vm


{:.input_area}
```R
V1 = substitute(a - b, list(a = Vm1, b = Vp1))
V2 = substitute(a - b, list(a = Vm2, b = Vp2))
```

Step 6: Generate the partial derivatives for the two Jacobians


{:.input_area}
```R
f11 = D(F1, "E"); f12 = D(F1, "I")
f21 = D(F2, "E"); f22 = D(F2, "I")

v11 = D(V1, "E"); v12 = D(V1, "I")
v21 = D(V2, "E"); v22 = D(V2, "I")
```

Step 7: Assuming N=1, the disease free equilibrium (dfe) is $S=1, E=0, I=0, R=0$. We also need values for other parameters. Assuming a weekly time-step and something chickenpox-like we may use $\mu=0$, $\alpha=0$, $\beta=5$, $\gamma=.8$, $\sigma=1.2$, and $N=1$.


{:.input_area}
```R
paras = list(S = 1, E = 0, I = 0, R = 0, mu = 0, 
   alpha = 0, beta = 5, gamma = .8, sigma = 1.2, N = 1)
f = with(paras,
   matrix(c(eval(f11),eval(f12),eval(f21),
   eval(f22)), nrow = 2, byrow = TRUE))
v=with(paras,
   matrix(c(eval(v11),eval(v12),eval(v21),
   eval(v22)), nrow=2, byrow=TRUE))
```

Step 8: Calculate the largest eigenvalue of $f \times \mbox{inverse}(v)$. Note that the function for inverting matrices in `R` is `solve`.


{:.input_area}
```R
max(eigen(f %*% solve(v))$values)
```


<div markdown="0">
6.25
</div>


Let us check that the next-generation method and the ‘flow’ method are in agreement recalling that for the SEIR-flow $R_0 =  \frac{\sigma}{\sigma +\mu} \frac{\beta}{\gamma+\mu+\alpha}$.


{:.input_area}
```R
with(paras,
sigma/(sigma+mu)*beta/(gamma+mu+alpha))
```


<div markdown="0">
6.25
</div>


#### SEIHFR

| Parameter    | Meaning                        |  Value|
|:-------------|:-------------------------------|------:|
| $N$          | Population size                |       |
| $1/\sigma$   | Incubation period              |     7d|
| $1/\gamma_h$ | Onset to hospitalization       |     5d|
| $1/\gamma_f$ | Onset to death                 |   9.6d|
| $1/\gamma_r$ | Onset to recovery              |    10d|
| $1/\eta_f$   | Hospitalization to death       |   4.6d|
| $1/\eta_r$   | Hospitalization to recovery    |     5d|
| $1/\chi$     | Death to burial                |     2d|
| $\Theta$     | Proportion hospitalized        |    80%|
| $\Lambda$    | Case fatality ratio            |    81%|
| $\beta_i$    | Transmission rate in community |  0.588|
| $\beta_h$    | Transmission rate in hospital  |  0.794|
| $\beta_f$    | Transmission rate at funeral   |  7.653|

The SEIHFR model forms the foundation for many of the recent Ebola models. The model has 5 compartments corresponding to Susceptible, Exposed, Infectious in community, Infectious in hospital, Dead but not yet buried, and removed (either buried or immune). The model is more complex than previous compartmental models and can not be represented by a simple linear chain. The parameterization used here is motivated by the original formulation of , but the notation conforms to the other sections of this book; Each infectious compartment contributes to the force of infection through their individual $\beta$s. There are two branching-points in the flow: The hospitalization of a fraction $\Theta$ of the infectious cases after an average time of $1/\gamma_h$ days following onset of symptoms, and the death of a fraction $\Lambda$ of the $I$- and $H$-class after an average time of $1/\gamma_f$ days and $1/\eta_f$ days, respectively. For the 1995 DRC outbreak, assumed that hospitalization affected transmission rates but not duration of infection or probability of dying. Model parameters are given above, and the model equations are:

$$\begin{aligned}
\frac{dS}{dt}&=&-(\beta_i I + \beta_h H +\beta_f F) S / N\\
\frac{dE}{dt}&=& (\beta_i I + \beta_h H +\beta_f F) S / N - \sigma E\\
\frac{dI}{dt}&=&\sigma E - \Theta \gamma_h I -(1-\Theta)(1-\Lambda)\gamma_r I - (1-\Theta) \Lambda \gamma_f I\\
\frac{dH}{dt}&=&\Theta \gamma_h I - \Lambda \eta_f H - (1-\Lambda) \eta_r H\\
\frac{dF}{dt}&=& (1-\Theta)(1-\Lambda)\gamma_r I+ \Lambda \eta_f H - \chi F\\
\frac{dR}{dt}&=&(1-\Theta)(1-\Lambda)\gamma_r I + (1-\Lambda) \eta_r H+\chi F
\end{aligned}$$

There are 4 infected compartments ($E$, $I$, $H$ and $F$), thus $\vec{F}$, $\vec{V^-}$, and $\vec{V^+}$ will be $4 \times 1$ matrices, and $\vec{f}$ and $\vec{v}$ will be $4 \times 4$ matrices.

Step 1: Infected classes are $E$, $I$, $H$ and $F$, let us label them F1 - 4.

Step 2: All new infections $dE/dt = \beta S I / N$, $dI/dt =0$


{:.input_area}
```R
F1 = expression(betai * S * I / N + betah* S * H / N +
   betaf * S * F / N)
F2=0
F3=0
F4=0
```

Step 3: All losses


{:.input_area}
```R
Vm1 = quote(sigma * E)
Vm2 = quote(Theta * gammah * I + (1 - Theta) * (1-
   Lambda) * gammar * I + (1 - Theta) * Lambda * 
   gammaf * I)
Vm3 = quote(Lambda * etaf * H + (1 - Lambda) * etar * H)
Vm4 = quote(chi * F)
```

Step 4: All gained transfers


{:.input_area}
```R
Vp1 = 0
Vp2 = quote(sigma * E)
Vp3 = quote(Theta * gammah * I)
Vp4 = quote((1 - Theta) * (1 - Lambda) * gammar * I+ 
   Lambda * etaf * H)
```

Step 5: Subtract Vp from Vm


{:.input_area}
```R
V1 = substitute(a - b, list(a = Vm1, b = Vp1))
V2 = substitute(a - b, list(a = Vm2, b = Vp2))
V3 = substitute(a - b, list(a = Vm3, b = Vp3))
V4 = substitute(a - b, list(a = Vm4, b = Vp4))
```

Step 6: Generate the partial derivatives for the two Jacobians


{:.input_area}
```R
f11 = D(F1, "E"); f12 = D(F1, "I"); f13 = D(F1, "H") 
     f14 = D(F1, "F")
f21 = D(F2, "E"); f22 = D(F2, "I"); f23 = D(F2, "H") 
     f24 = D(F2, "F")
f31 = D(F3, "E"); f32 = D(F3, "I"); f33 = D(F3, "H") 
     f34 = D(F3, "F")
f41 = D(F4, "E"); f42 = D(F4, "I"); f43 = D(F4, "H") 
     f44 = D(F4, "F")

v11 = D(V1, "E"); v12 = D(V1, "I"); v13 = D(V1, "H")
     v14 = D(V1, "F")
v21 = D(V2, "E"); v22 = D(V2, "I"); v23 = D(V2, "H")
     v24 = D(V2, "F")
v31 = D(V3, "E"); v32 = D(V3, "I"); v33 = D(V3, "H")
     v34 = D(V3, "F")
v41 = D(V4, "E"); v42 = D(V4, "I"); v43 = D(V4, "H")
     v44 = D(V4, "F")
```

Step 7: Disease free equilibrium: the dfe is $S=1, E=0, I=0, H=0, F=0, R=0$. We also need values for other parameters. We use the estimates from the DRC 1995 outbreak scaled as weekly rates from tables and appendices of .


{:.input_area}
```R
gammah = 1/5 * 7
gammaf = 1/9.6 * 7
gammar = 1/10 * 7
chi = 1/2 * 7
etaf = 1/4.6 * 7
etar = 1/5 * 7
paras = list(S = 1,E = 0, I = 0, H = 0, F = 0,R = 0,
     sigma = 1/7*7, Theta = 0.81, Lambda = 0.81, betai = 0.588, 
     betah = 0.794, betaf = 7.653, N = 1, gammah = gammah,
     gammaf = gammaf, gammar = gammar, etaf = etaf, 
     etar = etar, chi = chi)

f = with(paras, 
matrix(c(eval(f11), eval(f12), eval(f13), eval(f14),
     eval(f21), eval(f22), eval(f23), eval(f24),
     eval(f31), eval(f32), eval(f33), eval(f34),
     eval(f41), eval(f42), eval(f43), eval(f44)),
     nrow = 4, byrow = T))

v = with(paras, 
matrix(c(eval(v11), eval(v12), eval(v13), eval(v14),
     eval(v21), eval(v22), eval(v23), eval(v24),
     eval(v31), eval(v32), eval(v33), eval(v34),
     eval(v41), eval(v42), eval(v43), eval(v44)),
     nrow = 4, byrow = T))
```

Step 8: Calculate the largest eigenvalue of $f \times \mbox{inverse}(v)$


{:.input_area}
```R
max(eigen(f %*% solve(v))$values)
```


<div markdown="0">
2.58242862879894
</div>

