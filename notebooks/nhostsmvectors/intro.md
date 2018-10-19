
## $N$ host SEIR, $M$ vector SEI model

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
