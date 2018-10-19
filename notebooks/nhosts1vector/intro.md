## One Host, One Vector


*Author*: Carl A. B. Pearson @pearsonca
*Date*: 2018-10-02

For this simple model, we represent a single host, $H$, and single vector, $V$, which each start as susceptible ($S$ compartment for both species) to an infectious agent.  When infected, each species enters an incubation period ($E$) before becoming infectious ($I$).  The host species can clear the infection, after which we assume is subsequently immune to infection (which is typical for viral pathogens, though decidely not accurate for other important infections, like malaria).  We represent population turnover in both the host and vector species.  Lastly, we assume that the probabilty of infection is the same for both host and vector - i.e., a susceptible host being bit by an infectious vector becomes infected with the same probabilty as a susceptible vector biting an infectious host.

The model parameters are:

 - $\sigma_H$, $\sigma_V$: the incubation rates for hosts & vectors (units: per time)
 - $\mu_H$, $\mu_V$: the mortality rates for hosts & vectors (units: per time)
 - $\lambda$: the clearance (or recovery) rate for hosts (units: per time)
 - $\beta$: the infection rate (units: per capita per time)

The infection rate is the combination of two factors: the number of bites per vector per time ($c$), the probability of infection per bite ($p$).

Another way to think about this factor:

host infections per time = new infections (infectious bites) per mosquito per time * # of $I$ mosquitos * fraction striking susceptible hosts ($S_H/N_H$)

mosquito infections per time = new infections (infectious bites) per mosquito per time * # of S mosquitos * fraction striking infectious hosts ($I_H/N_H$)

The state equations are:

$$
N_H = S_H + E_H + I_H + R_H, \dot{N_H}=0
\dot{S_H} = \mu_H N_H - \frac{\beta}{N_H} S_H I_V - \mu_H S_H = \mu_H(E_H + I_H + R_H) - \frac{\beta}{N_H} S_H I_V
\dot{E_H} = \frac{\beta}{N_H} S_H I_V - (\sigma_H + \mu_H) E_H
\dot{I_H} = \sigma_H E_H - (\lambda + \mu_H) I_H
\dot{R_H} = \lambda I_H - \mu_H R_H
\dot{S_V} = \mu_V N_V - \frac{\beta}{N_H} S_V I_H - \mu_V S_V = \mu_V(E_V + I_V) - \frac{\beta}{N_H} S_V I_H
\dot{E_V} = \frac{\beta}{N_H} S_V I_H - (\sigma_V + \mu_V) E_V
\dot{I_V} = \sigma_V E_V - \mu_V I_V
$$

#### $N$ Hosts, One Vector

With $N$ hosts, these equations change to accommodate different properties between hosts, and the fact that multiple hosts contribute to mosquito infection:

$$
N_H^i = S_H^i + E_H^i + I_H^i + R_H^i, \dot{N_H^i}=0
N_H = \sum_{i=1}^N N_H^i
\dot{S_H}^i = \mu_H^i(E_H^i + I_H^i + R_H^i) - \frac{\beta^i}{N_H} S_H^i I_V
\dot{E_H}^i = \frac{\beta^i}{N_H} S_H^i I_V - (\sigma_H^i + \mu_H^i) E_H^i
\dot{I_H}^i = \sigma_H^i E_H^i - (\lambda^i + \mu_H^i) I_H^i
\dot{R_H}^i = \lambda^i I_H^i - \mu_H^i R_H^i
\dot{S_V} = \mu_V(E_V + I_V) - \frac{\sum_{i=1}^N\beta^i I_H^i}{N_H} S_V
\dot{E_V} = \frac{\sum_{i=1}^N\beta^i I_H^i}{N_H} S_V - (\sigma_V + \mu_V) E_V
\dot{I_V} = \sigma_V E_V - \mu_V I_V
$$

#### $N$ Hosts, $M$ Vectors

With $M$ vectors as well, this $N$ hosts-style extension can be repeated: 

$$
N_H^i = S_H^i + E_H^i + I_H^i + R_H^i, \dot{N_H^i}=0
N_H = \sum_{i=1}^N N_H^i
\dot{S_H}^i = \mu_H^i(E_H^i + I_H^i + R_H^i) - \frac{\sum_{j=1}^M\beta^{ij}I_V^j}{N_H} S_H^i
\dot{E_H}^i = \frac{\sum_{j=1}^M\beta^{ij}I_V^j}{N_H} S_H^i - (\sigma_H^i + \mu_H^i) E_H^i
\dot{I_H}^i = \sigma_H^i E_H^i - (\lambda^i + \mu_H^i) I_H^i
\dot{R_H}^i = \lambda^i I_H^i - \mu_H^i R_H^i
\dot{S_V}^j = \mu_V^j(E_V^j + I_V^j) - \frac{\sum_{i=1}^N\beta^{ij} I_H^i}{N_H} S_V^j
\dot{E_V}^j = \frac{\sum_{i=1}^N\beta^{ij} I_H^i}{N_H} S_V^j - (\sigma_V^j + \mu_V^j) E_V^j
\dot{I_V}^j = \sigma_V^j E_V^j - \mu_V^j I_V^j
$$
