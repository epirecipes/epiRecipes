---
interact_link: notebooks/1host1vector/julia.ipynb
title: 'Julia'
permalink: 'chapters/1host1vector/julia'
previouschapter:
  url: chapters/1host1vector/intro
  title: 'One host, one vector'
nextchapter:
  url: chapters/nhosts1vector/intro
  title: 'N hosts, one vector'
redirect_from:
  - 'chapters/1host1vector/julia'
---

## One Host SEIR + One Vector SEI model in Julia

*Author*: Carl A. B. Pearson  
*Date*: 2 OCT 2018

First, create a system of ODEs conforming to the required signature:


{:.input_area}
```julia
function F(du,u,p,t)
    S_H, E_H, I_H, R_H, S_V, E_V, I_V = u
    
    # host dynamics
    host_infection = (p.β*S_H*I_V)/p.N_H
    host_mortality = p.μ_H .* u[1:4] # include S_H, so easier to remove mortality
    host_births = sum(host_mortality)
    host_progression = p.σ_H*E_H
    recovery = p.λ*I_H
    
    du[1] = -host_infection + host_births
    du[2] = host_infection - host_progression
    du[3] = host_progression - recovery
    du[4] = recovery
    du[1:4] -= host_mortality
    
    # vector dynamics
    vec_infection = (p.β*S_V*I_H)/p.N_H
    vec_mortality = p.μ_V .* u[5:7] # include S_V, so easier to remove mortality
    vec_births = sum(vec_mortality)
    vec_progression = p.σ_V*E_V
    
    du[5] = -vec_infection + vec_births
    du[6] = vec_infection - vec_progression
    du[7] = vec_progression
    du[5:7] -= vec_mortality
    
end
```




{:.output_data_text}
```
F (generic function with 1 method)
```



Set initial conditions and dynamic parameters, then apply the ODE solver.  Note this will also display the solver run time; this will be slow for the first run, though quite fast if you re-run it (even after changing parameters or initial conditions).


{:.input_area}
```julia
using DifferentialEquations
using IterableTables, DataFrames
using NamedTuples
# nb: in >= Julia v0.7, can eliminate this import
#  and the @NT syntax
u0 = [
    S_H=100.0,   E_H=0.0, I_H=1.0, R_H=0.0,
    S_V=10000.0, E_V=0.0, I_V=0.0
]
p = @NT(
  μ_H=1/365, μ_V=1/30, σ_H=1/3, σ_V=1/7, λ=1/14,
  β=0.05, N_H = sum(u0[1:4])
)
tspan = (0.0, 365.0)
prob = ODEProblem(F, u0, tspan, p)
sol = @time solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8,saveat=linspace(0,365,365*10+1))

df = DataFrame(sol)
rename!(df, Dict(:timestamp => :t,
  :value1 => :S_H, :value2 => :E_H, :value3 => :I_H, :value4 => :R_H,
  :value5 => :S_V, :value6 => :E_V, :value7 => :I_V
))
mlt = melt(df,:t) # convert results into long format for plotting
mlt[:host] = contains.(string.(mlt[:variable]),"H"); # tag which entries are host vs vector
df
```

{:.output_stream}
```
  3.965556 seconds (3.46 M allocations: 192.139 MiB, 1.64% gc time)

```




<div markdown="0">
<table class="data-frame"><thead><tr><th></th><th>t</th><th>S_H</th><th>E_H</th><th>I_H</th><th>R_H</th><th>S_V</th><th>E_V</th><th>I_V</th></tr></thead><tbody><tr><th>1</th><td>0.0</td><td>100.0</td><td>0.0</td><td>1.0</td><td>0.0</td><td>10000.0</td><td>0.0</td><td>0.0</td></tr><tr><th>2</th><td>0.1</td><td>100.00026814164985</td><td>5.745401288701111e-6</td><td>0.9926106550187419</td><td>0.007115457930112796</td><td>9999.50761607706</td><td>0.4888811777907883</td><td>0.0035027451477435017</td></tr><tr><th>3</th><td>0.2</td><td>100.00050177829016</td><td>4.525792070955373e-5</td><td>0.9852765739814515</td><td>0.014176389807693479</td><td>9999.020532914787</td><td>0.9655880127981552</td><td>0.013879072414963943</td></tr><tr><th>4</th><td>0.3</td><td>100.00066737984612</td><td>0.00015040567741793047</td><td>0.9779990053871515</td><td>0.021183209089312854</td><td>9998.538704686205</td><td>1.4303611143782462</td><td>0.030934199416335188</td></tr><tr><th>5</th><td>0.4</td><td>100.00073237948628</td><td>0.0003510652956395811</td><td>0.9707802135129932</td><td>0.028136341705096376</td><td>9998.062084857795</td><td>1.8834377432379645</td><td>0.05447739896513851</td></tr><tr><th>6</th><td>0.5</td><td>100.00066515369967</td><td>0.0006752070897406001</td><td>0.963623406156655</td><td>0.03503623305394339</td><td>9997.5906257077</td><td>2.3250523506169727</td><td>0.08432194168190127</td></tr><tr><th>7</th><td>0.6</td><td>100.00043500273485</td><td>0.0011489770693011803</td><td>0.9565326657056807</td><td>0.041883354490183024</td><td>9997.124277880717</td><td>2.7554370730797264</td><td>0.12028504620147329</td></tr><tr><th>8</th><td>0.7</td><td>100.00001213141086</td><td>0.0017967756786389124</td><td>0.9495128836036143</td><td>0.04867820930690247</td><td>9996.662989979597</td><td>3.1748221839374784</td><td>0.16218783646455837</td></tr><tr><th>9</th><td>0.8</td><td>99.99936763012019</td><td>0.0026413340720121085</td><td>0.9425696975069388</td><td>0.05542133830087409</td><td>9996.206708186579</td><td>3.5834365085203954</td><td>0.20985530490011298</td></tr><tr><th>10</th><td>0.9</td><td>99.99847345594463</td><td>0.00370378797506614</td><td>0.935709431179207</td><td>0.06211332490111005</td><td>9995.755375916528</td><td>3.9815078018550603</td><td>0.2631162816177185</td></tr><tr><th>11</th><td>1.0</td><td>99.99730241455798</td><td>0.00500374674442986</td><td>0.9289390389534903</td><td>0.06875479974410777</td><td>9995.308933509126</td><td>4.369263081472214</td><td>0.32180340940160024</td></tr><tr><th>12</th><td>1.1</td><td>99.9958281411404</td><td>0.006559364435790144</td><td>0.9222660492620383</td><td>0.07534644516179213</td><td>9994.867317928649</td><td>4.746928949742029</td><td>0.38575312160914144</td></tr><tr><th>13</th><td>1.2</td><td>99.99402508304243</td><td>0.008387402053565977</td><td>0.9156985160745399</td><td>0.08188899882947635</td><td>9994.43046252214</td><td>5.11473185068202</td><td>0.45480562717867207</td></tr><tr><th>14</th><td>1.3</td><td>99.99186848062837</td><td>0.010503295303662025</td><td>0.9092449664472831</td><td>0.08838325762069417</td><td>9993.998296764821</td><td>5.4728983393429</td><td>0.528804895835597</td></tr><tr><th>15</th><td>1.4</td><td>99.98933435043577</td><td>0.012921211210917136</td><td>0.9029143578800672</td><td>0.09483008047325146</td><td>9993.57074607381</td><td>5.821655275564482</td><td>0.6075986506255675</td></tr><tr><th>16</th><td>1.5</td><td>99.98639946576401</td><td>0.01565411301160372</td><td>0.8967160295734984</td><td>0.10123039165090672</td><td>9993.147731596253</td><td>6.161230044422309</td><td>0.6910383593243803</td></tr><tr><th>17</th><td>1.6</td><td>99.98304134053252</td><td>0.01871381087842438</td><td>0.8906596656812472</td><td>0.10758518290781208</td><td>9992.72917007232</td><td>6.491850695317421</td><td>0.7789792323625186</td></tr><tr><th>18</th><td>1.7</td><td>99.97923820920528</td><td>0.022111024888270103</td><td>0.8847552497174251</td><td>0.11389551618904671</td><td>9992.314973663997</td><td>6.813746116346348</td><td>0.8712802196570351</td></tr><tr><th>19</th><td>1.8</td><td>99.9749690114076</td><td>0.025855430503449094</td><td>0.8790130328327947</td><td>0.12016252525617195</td><td>9991.905049855464</td><td>7.127146132897464</td><td>0.9678040116372855</td></tr><tr><th>20</th><td>1.9</td><td>99.97021337156639</td><td>0.029955717708979524</td><td>0.8734434930259793</td><td>0.12638741769866832</td><td>9991.499301331016</td><td>7.432281627876439</td><td>1.0684170411069502</td></tr><tr><th>21</th><td>2.0</td><td>99.9649515828934</td><td>0.03441963633162285</td><td>0.8680573044314721</td><td>0.1325714763435111</td><td>9991.097625891161</td><td>7.7293846231762835</td><td>1.1729894856613798</td></tr><tr><th>22</th><td>2.1</td><td>99.95916458894199</td><td>0.03925404480785851</td><td>0.8628653058725415</td><td>0.13871606037763487</td><td>9990.699916391362</td><td>8.018688335574598</td><td>1.2813952730626423</td></tr><tr><th>23</th><td>2.2</td><td>99.95283396455962</td><td>0.044464961696772165</td><td>0.8578784670388516</td><td>0.1448226067047723</td><td>9990.306060664585</td><td>8.300427250046441</td><td>1.3935120853669343</td></tr><tr><th>24</th><td>2.3</td><td>99.94594190064402</td><td>0.05005760246742591</td><td>0.853107866491243</td><td>0.15089263039732773</td><td>9989.915941503974</td><td>8.57483713113998</td><td>1.5092213648849446</td></tr><tr><th>25</th><td>2.4</td><td>99.9384711825094</td><td>0.05603643341442274</td><td>0.8485646584553629</td><td>0.15692772562082902</td><td>9989.52943661802</td><td>8.842155061154639</td><td>1.6284083208239846</td></tr><tr><th>26</th><td>2.5</td><td>99.93040517498355</td><td>0.06240520682808333</td><td>0.8442600521800281</td><td>0.16292956600835462</td><td>9989.14641861804</td><td>9.102619447381544</td><td>1.7509619345776433</td></tr><tr><th>27</th><td>2.6</td><td>99.92172780354208</td><td>0.06916700268500087</td><td>0.8402052889902647</td><td>0.16889990478267222</td><td>9988.76675502676</td><td>9.356470006726369</td><td>1.8767749665126119</td></tr><tr><th>28</th><td>2.7</td><td>99.91242353405677</td><td>0.07632427518811552</td><td>0.8364116156084753</td><td>0.17484057514665405</td><td>9988.390308269225</td><td>9.603947768150924</td><td>2.00574396262303</td></tr><tr><th>29</th><td>2.8</td><td>99.90247735877468</td><td>0.08387888092361276</td><td>0.8328902702538351</td><td>0.18075349004788635</td><td>9988.016935701668</td><td>9.845295039200272</td><td>2.1377692591325457</td></tr><tr><th>30</th><td>2.9</td><td>99.89187477418837</td><td>0.0918321245595308</td><td>0.8296524592467466</td><td>0.18664064200537545</td><td>9987.646489642992</td><td>10.080755368056721</td><td>2.2727549889506884</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table>
</div>



Now that we have a solution, we want to view what is happening in host vs mosquito population:


{:.input_area}
```julia
using Gadfly
fig1a = plot(mlt[mlt[:host] .== true,:], x=:t, y=:value, color=:variable, Geom.line)
fig1b = plot(mlt[mlt[:host] .!= true,:], x=:t, y=:value, color=:variable, Geom.line)
vstack(fig1a,fig1b)
```




<div markdown="0">
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     xmlns:gadfly="http://www.gadflyjl.org/ns"
     version="1.2"
     width="141.42mm" height="100mm" viewBox="0 0 141.42 100"
     stroke="none"
     fill="#000000"
     stroke-width="0.3"
     font-size="3.88"

     id="img-8311e7c5">
<g class="plotroot xscalable yscalable" id="img-8311e7c5-1">
  <g font-size="3.88" font-family="'PT Sans','Helvetica Neue','Helvetica',sans-serif" fill="#564A55" stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-2">
    <g transform="translate(74.93,88.39)">
      <g class="primitive">
        <text text-anchor="middle" dy="0.6em">t</text>
      </g>
    </g>
  </g>
  <g class="guide xlabels" font-size="2.82" font-family="'PT Sans Caption','Helvetica Neue','Helvetica',sans-serif" fill="#6C606B" id="img-8311e7c5-3">
    <g transform="translate(-86.79,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-500</text>
      </g>
    </g>
    <g transform="translate(-63.69,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-400</text>
      </g>
    </g>
    <g transform="translate(-40.59,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-300</text>
      </g>
    </g>
    <g transform="translate(-17.48,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-200</text>
      </g>
    </g>
    <g transform="translate(5.62,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-100</text>
      </g>
    </g>
    <g transform="translate(28.72,84.39)" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">0</text>
      </g>
    </g>
    <g transform="translate(51.83,84.39)" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">100</text>
      </g>
    </g>
    <g transform="translate(74.93,84.39)" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">200</text>
      </g>
    </g>
    <g transform="translate(98.04,84.39)" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">300</text>
      </g>
    </g>
    <g transform="translate(121.14,84.39)" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">400</text>
      </g>
    </g>
    <g transform="translate(144.24,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">500</text>
      </g>
    </g>
    <g transform="translate(167.35,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">600</text>
      </g>
    </g>
    <g transform="translate(190.45,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">700</text>
      </g>
    </g>
    <g transform="translate(213.55,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">800</text>
      </g>
    </g>
    <g transform="translate(236.66,84.39)" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">900</text>
      </g>
    </g>
    <g transform="translate(-63.69,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-400</text>
      </g>
    </g>
    <g transform="translate(-59.07,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-380</text>
      </g>
    </g>
    <g transform="translate(-54.45,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-360</text>
      </g>
    </g>
    <g transform="translate(-49.83,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-340</text>
      </g>
    </g>
    <g transform="translate(-45.21,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-320</text>
      </g>
    </g>
    <g transform="translate(-40.59,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-300</text>
      </g>
    </g>
    <g transform="translate(-35.96,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-280</text>
      </g>
    </g>
    <g transform="translate(-31.34,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-260</text>
      </g>
    </g>
    <g transform="translate(-26.72,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-240</text>
      </g>
    </g>
    <g transform="translate(-22.1,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-220</text>
      </g>
    </g>
    <g transform="translate(-17.48,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-200</text>
      </g>
    </g>
    <g transform="translate(-12.86,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-180</text>
      </g>
    </g>
    <g transform="translate(-8.24,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-160</text>
      </g>
    </g>
    <g transform="translate(-3.62,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-140</text>
      </g>
    </g>
    <g transform="translate(1,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-120</text>
      </g>
    </g>
    <g transform="translate(5.62,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-100</text>
      </g>
    </g>
    <g transform="translate(10.24,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-80</text>
      </g>
    </g>
    <g transform="translate(14.86,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-60</text>
      </g>
    </g>
    <g transform="translate(19.48,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-40</text>
      </g>
    </g>
    <g transform="translate(24.1,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-20</text>
      </g>
    </g>
    <g transform="translate(28.72,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">0</text>
      </g>
    </g>
    <g transform="translate(33.35,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">20</text>
      </g>
    </g>
    <g transform="translate(37.97,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">40</text>
      </g>
    </g>
    <g transform="translate(42.59,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">60</text>
      </g>
    </g>
    <g transform="translate(47.21,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">80</text>
      </g>
    </g>
    <g transform="translate(51.83,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">100</text>
      </g>
    </g>
    <g transform="translate(56.45,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">120</text>
      </g>
    </g>
    <g transform="translate(61.07,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">140</text>
      </g>
    </g>
    <g transform="translate(65.69,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">160</text>
      </g>
    </g>
    <g transform="translate(70.31,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">180</text>
      </g>
    </g>
    <g transform="translate(74.93,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">200</text>
      </g>
    </g>
    <g transform="translate(79.55,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">220</text>
      </g>
    </g>
    <g transform="translate(84.17,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">240</text>
      </g>
    </g>
    <g transform="translate(88.79,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">260</text>
      </g>
    </g>
    <g transform="translate(93.41,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">280</text>
      </g>
    </g>
    <g transform="translate(98.04,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">300</text>
      </g>
    </g>
    <g transform="translate(102.66,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">320</text>
      </g>
    </g>
    <g transform="translate(107.28,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">340</text>
      </g>
    </g>
    <g transform="translate(111.9,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">360</text>
      </g>
    </g>
    <g transform="translate(116.52,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">380</text>
      </g>
    </g>
    <g transform="translate(121.14,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">400</text>
      </g>
    </g>
    <g transform="translate(125.76,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">420</text>
      </g>
    </g>
    <g transform="translate(130.38,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">440</text>
      </g>
    </g>
    <g transform="translate(135,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">460</text>
      </g>
    </g>
    <g transform="translate(139.62,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">480</text>
      </g>
    </g>
    <g transform="translate(144.24,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">500</text>
      </g>
    </g>
    <g transform="translate(148.86,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">520</text>
      </g>
    </g>
    <g transform="translate(153.48,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">540</text>
      </g>
    </g>
    <g transform="translate(158.1,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">560</text>
      </g>
    </g>
    <g transform="translate(162.73,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">580</text>
      </g>
    </g>
    <g transform="translate(167.35,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">600</text>
      </g>
    </g>
    <g transform="translate(171.97,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">620</text>
      </g>
    </g>
    <g transform="translate(176.59,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">640</text>
      </g>
    </g>
    <g transform="translate(181.21,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">660</text>
      </g>
    </g>
    <g transform="translate(185.83,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">680</text>
      </g>
    </g>
    <g transform="translate(190.45,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">700</text>
      </g>
    </g>
    <g transform="translate(195.07,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">720</text>
      </g>
    </g>
    <g transform="translate(199.69,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">740</text>
      </g>
    </g>
    <g transform="translate(204.31,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">760</text>
      </g>
    </g>
    <g transform="translate(208.93,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">780</text>
      </g>
    </g>
    <g transform="translate(213.55,84.39)" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">800</text>
      </g>
    </g>
    <g transform="translate(-86.79,84.39)" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-500</text>
      </g>
    </g>
    <g transform="translate(28.72,84.39)" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">0</text>
      </g>
    </g>
    <g transform="translate(144.24,84.39)" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">500</text>
      </g>
    </g>
    <g transform="translate(259.76,84.39)" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">1000</text>
      </g>
    </g>
    <g transform="translate(-63.69,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-400</text>
      </g>
    </g>
    <g transform="translate(-52.14,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-350</text>
      </g>
    </g>
    <g transform="translate(-40.59,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-300</text>
      </g>
    </g>
    <g transform="translate(-29.03,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-250</text>
      </g>
    </g>
    <g transform="translate(-17.48,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-200</text>
      </g>
    </g>
    <g transform="translate(-5.93,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-150</text>
      </g>
    </g>
    <g transform="translate(5.62,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-100</text>
      </g>
    </g>
    <g transform="translate(17.17,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-50</text>
      </g>
    </g>
    <g transform="translate(28.72,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">0</text>
      </g>
    </g>
    <g transform="translate(40.28,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">50</text>
      </g>
    </g>
    <g transform="translate(51.83,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">100</text>
      </g>
    </g>
    <g transform="translate(63.38,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">150</text>
      </g>
    </g>
    <g transform="translate(74.93,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">200</text>
      </g>
    </g>
    <g transform="translate(86.48,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">250</text>
      </g>
    </g>
    <g transform="translate(98.04,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">300</text>
      </g>
    </g>
    <g transform="translate(109.59,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">350</text>
      </g>
    </g>
    <g transform="translate(121.14,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">400</text>
      </g>
    </g>
    <g transform="translate(132.69,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">450</text>
      </g>
    </g>
    <g transform="translate(144.24,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">500</text>
      </g>
    </g>
    <g transform="translate(155.79,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">550</text>
      </g>
    </g>
    <g transform="translate(167.35,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">600</text>
      </g>
    </g>
    <g transform="translate(178.9,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">650</text>
      </g>
    </g>
    <g transform="translate(190.45,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">700</text>
      </g>
    </g>
    <g transform="translate(202,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">750</text>
      </g>
    </g>
    <g transform="translate(213.55,84.39)" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">800</text>
      </g>
    </g>
  </g>
  <g class="guide colorkey" id="img-8311e7c5-4">
    <g fill="#4C404B" font-size="2.82" font-family="'PT Sans','Helvetica Neue','Helvetica',sans-serif" id="img-8311e7c5-5">
      <g transform="translate(126.95,66.04)" id="img-8311e7c5-6" class="color_S_V">
        <g class="primitive">
          <text dy="0.35em">S_V</text>
        </g>
      </g>
      <g transform="translate(126.95,69.67)" id="img-8311e7c5-7" class="color_E_V">
        <g class="primitive">
          <text dy="0.35em">E_V</text>
        </g>
      </g>
      <g transform="translate(126.95,73.3)" id="img-8311e7c5-8" class="color_I_V">
        <g class="primitive">
          <text dy="0.35em">I_V</text>
        </g>
      </g>
    </g>
    <g stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-9">
      <g transform="translate(125.05,66.04)" id="img-8311e7c5-10" fill="#00BFFF" class="color_S_V">
        <path d="M-0.91,-0.91 L 0.91 -0.91 0.91 0.91 -0.91 0.91 z" class="primitive"/>
      </g>
      <g transform="translate(125.05,69.67)" id="img-8311e7c5-11" fill="#D4CA3A" class="color_E_V">
        <path d="M-0.91,-0.91 L 0.91 -0.91 0.91 0.91 -0.91 0.91 z" class="primitive"/>
      </g>
      <g transform="translate(125.05,73.3)" id="img-8311e7c5-12" fill="#FF6DAE" class="color_I_V">
        <path d="M-0.91,-0.91 L 0.91 -0.91 0.91 0.91 -0.91 0.91 z" class="primitive"/>
      </g>
    </g>
    <g fill="#362A35" font-size="3.88" font-family="'PT Sans','Helvetica Neue','Helvetica',sans-serif" stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-13">
      <g transform="translate(124.14,62.22)" id="img-8311e7c5-14">
        <g class="primitive">
          <text>variable</text>
        </g>
      </g>
    </g>
  </g>
  <g clip-path="url(#img-8311e7c5-15)">
    <g id="img-8311e7c5-16">
      <g pointer-events="visible" opacity="1" fill="#000000" fill-opacity="0.000" stroke="#000000" stroke-opacity="0.000" class="guide background" id="img-8311e7c5-17">
        <g transform="translate(74.93,67.86)" id="img-8311e7c5-18">
          <path d="M-48.21,-12.86 L 48.21 -12.86 48.21 12.86 -48.21 12.86 z" class="primitive"/>
        </g>
      </g>
      <g class="guide ygridlines xfixed" stroke-dasharray="0.5,0.5" stroke-width="0.2" stroke="#D0D0E0" id="img-8311e7c5-19">
        <g transform="translate(74.93,111.29)" id="img-8311e7c5-20" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,100.43)" id="img-8311e7c5-21" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,89.57)" id="img-8311e7c5-22" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,78.72)" id="img-8311e7c5-23" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,67.86)" id="img-8311e7c5-24" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,57)" id="img-8311e7c5-25" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,46.14)" id="img-8311e7c5-26" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,35.28)" id="img-8311e7c5-27" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,24.43)" id="img-8311e7c5-28" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,100.43)" id="img-8311e7c5-29" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,99.34)" id="img-8311e7c5-30" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,98.26)" id="img-8311e7c5-31" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,97.17)" id="img-8311e7c5-32" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,96.09)" id="img-8311e7c5-33" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,95)" id="img-8311e7c5-34" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,93.92)" id="img-8311e7c5-35" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,92.83)" id="img-8311e7c5-36" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,91.74)" id="img-8311e7c5-37" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,90.66)" id="img-8311e7c5-38" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,89.57)" id="img-8311e7c5-39" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,88.49)" id="img-8311e7c5-40" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,87.4)" id="img-8311e7c5-41" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,86.32)" id="img-8311e7c5-42" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,85.23)" id="img-8311e7c5-43" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,84.14)" id="img-8311e7c5-44" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,83.06)" id="img-8311e7c5-45" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,81.97)" id="img-8311e7c5-46" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,80.89)" id="img-8311e7c5-47" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,79.8)" id="img-8311e7c5-48" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,78.72)" id="img-8311e7c5-49" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,77.63)" id="img-8311e7c5-50" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,76.54)" id="img-8311e7c5-51" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,75.46)" id="img-8311e7c5-52" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,74.37)" id="img-8311e7c5-53" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,73.29)" id="img-8311e7c5-54" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,72.2)" id="img-8311e7c5-55" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,71.11)" id="img-8311e7c5-56" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,70.03)" id="img-8311e7c5-57" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,68.94)" id="img-8311e7c5-58" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,67.86)" id="img-8311e7c5-59" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,66.77)" id="img-8311e7c5-60" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,65.69)" id="img-8311e7c5-61" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,64.6)" id="img-8311e7c5-62" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,63.51)" id="img-8311e7c5-63" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,62.43)" id="img-8311e7c5-64" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,61.34)" id="img-8311e7c5-65" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,60.26)" id="img-8311e7c5-66" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,59.17)" id="img-8311e7c5-67" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,58.09)" id="img-8311e7c5-68" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,57)" id="img-8311e7c5-69" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,55.91)" id="img-8311e7c5-70" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,54.83)" id="img-8311e7c5-71" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,53.74)" id="img-8311e7c5-72" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,52.66)" id="img-8311e7c5-73" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,51.57)" id="img-8311e7c5-74" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,50.49)" id="img-8311e7c5-75" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,49.4)" id="img-8311e7c5-76" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,48.31)" id="img-8311e7c5-77" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,47.23)" id="img-8311e7c5-78" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,46.14)" id="img-8311e7c5-79" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,45.06)" id="img-8311e7c5-80" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,43.97)" id="img-8311e7c5-81" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,42.89)" id="img-8311e7c5-82" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,41.8)" id="img-8311e7c5-83" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,40.71)" id="img-8311e7c5-84" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,39.63)" id="img-8311e7c5-85" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,38.54)" id="img-8311e7c5-86" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,37.46)" id="img-8311e7c5-87" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,36.37)" id="img-8311e7c5-88" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,35.28)" id="img-8311e7c5-89" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,100.43)" id="img-8311e7c5-90" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,78.72)" id="img-8311e7c5-91" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,57)" id="img-8311e7c5-92" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,35.28)" id="img-8311e7c5-93" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,100.43)" id="img-8311e7c5-94" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,98.26)" id="img-8311e7c5-95" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,96.09)" id="img-8311e7c5-96" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,93.92)" id="img-8311e7c5-97" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,91.74)" id="img-8311e7c5-98" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,89.57)" id="img-8311e7c5-99" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,87.4)" id="img-8311e7c5-100" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,85.23)" id="img-8311e7c5-101" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,83.06)" id="img-8311e7c5-102" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,80.89)" id="img-8311e7c5-103" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,78.72)" id="img-8311e7c5-104" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,76.54)" id="img-8311e7c5-105" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,74.37)" id="img-8311e7c5-106" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,72.2)" id="img-8311e7c5-107" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,70.03)" id="img-8311e7c5-108" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,67.86)" id="img-8311e7c5-109" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,65.69)" id="img-8311e7c5-110" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,63.51)" id="img-8311e7c5-111" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,61.34)" id="img-8311e7c5-112" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,59.17)" id="img-8311e7c5-113" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,57)" id="img-8311e7c5-114" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,54.83)" id="img-8311e7c5-115" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,52.66)" id="img-8311e7c5-116" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,50.49)" id="img-8311e7c5-117" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,48.31)" id="img-8311e7c5-118" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,46.14)" id="img-8311e7c5-119" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,43.97)" id="img-8311e7c5-120" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,41.8)" id="img-8311e7c5-121" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,39.63)" id="img-8311e7c5-122" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,37.46)" id="img-8311e7c5-123" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
        <g transform="translate(74.93,35.28)" id="img-8311e7c5-124" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-48.21,0 L 48.21 0" class="primitive"/>
        </g>
      </g>
      <g class="guide xgridlines yfixed" stroke-dasharray="0.5,0.5" stroke-width="0.2" stroke="#D0D0E0" id="img-8311e7c5-125">
        <g transform="translate(-86.79,67.86)" id="img-8311e7c5-126" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-63.69,67.86)" id="img-8311e7c5-127" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-40.59,67.86)" id="img-8311e7c5-128" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-17.48,67.86)" id="img-8311e7c5-129" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(5.62,67.86)" id="img-8311e7c5-130" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(28.72,67.86)" id="img-8311e7c5-131" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(51.83,67.86)" id="img-8311e7c5-132" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(74.93,67.86)" id="img-8311e7c5-133" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(98.04,67.86)" id="img-8311e7c5-134" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(121.14,67.86)" id="img-8311e7c5-135" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(144.24,67.86)" id="img-8311e7c5-136" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(167.35,67.86)" id="img-8311e7c5-137" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(190.45,67.86)" id="img-8311e7c5-138" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(213.55,67.86)" id="img-8311e7c5-139" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(236.66,67.86)" id="img-8311e7c5-140" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-63.69,67.86)" id="img-8311e7c5-141" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-59.07,67.86)" id="img-8311e7c5-142" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-54.45,67.86)" id="img-8311e7c5-143" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-49.83,67.86)" id="img-8311e7c5-144" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-45.21,67.86)" id="img-8311e7c5-145" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-40.59,67.86)" id="img-8311e7c5-146" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-35.96,67.86)" id="img-8311e7c5-147" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-31.34,67.86)" id="img-8311e7c5-148" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-26.72,67.86)" id="img-8311e7c5-149" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-22.1,67.86)" id="img-8311e7c5-150" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-17.48,67.86)" id="img-8311e7c5-151" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-12.86,67.86)" id="img-8311e7c5-152" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-8.24,67.86)" id="img-8311e7c5-153" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-3.62,67.86)" id="img-8311e7c5-154" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(1,67.86)" id="img-8311e7c5-155" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(5.62,67.86)" id="img-8311e7c5-156" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(10.24,67.86)" id="img-8311e7c5-157" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(14.86,67.86)" id="img-8311e7c5-158" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(19.48,67.86)" id="img-8311e7c5-159" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(24.1,67.86)" id="img-8311e7c5-160" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(28.72,67.86)" id="img-8311e7c5-161" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(33.35,67.86)" id="img-8311e7c5-162" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(37.97,67.86)" id="img-8311e7c5-163" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(42.59,67.86)" id="img-8311e7c5-164" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(47.21,67.86)" id="img-8311e7c5-165" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(51.83,67.86)" id="img-8311e7c5-166" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(56.45,67.86)" id="img-8311e7c5-167" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(61.07,67.86)" id="img-8311e7c5-168" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(65.69,67.86)" id="img-8311e7c5-169" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(70.31,67.86)" id="img-8311e7c5-170" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(74.93,67.86)" id="img-8311e7c5-171" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(79.55,67.86)" id="img-8311e7c5-172" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(84.17,67.86)" id="img-8311e7c5-173" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(88.79,67.86)" id="img-8311e7c5-174" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(93.41,67.86)" id="img-8311e7c5-175" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(98.04,67.86)" id="img-8311e7c5-176" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(102.66,67.86)" id="img-8311e7c5-177" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(107.28,67.86)" id="img-8311e7c5-178" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(111.9,67.86)" id="img-8311e7c5-179" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(116.52,67.86)" id="img-8311e7c5-180" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(121.14,67.86)" id="img-8311e7c5-181" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(125.76,67.86)" id="img-8311e7c5-182" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(130.38,67.86)" id="img-8311e7c5-183" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(135,67.86)" id="img-8311e7c5-184" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(139.62,67.86)" id="img-8311e7c5-185" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(144.24,67.86)" id="img-8311e7c5-186" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(148.86,67.86)" id="img-8311e7c5-187" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(153.48,67.86)" id="img-8311e7c5-188" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(158.1,67.86)" id="img-8311e7c5-189" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(162.73,67.86)" id="img-8311e7c5-190" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(167.35,67.86)" id="img-8311e7c5-191" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(171.97,67.86)" id="img-8311e7c5-192" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(176.59,67.86)" id="img-8311e7c5-193" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(181.21,67.86)" id="img-8311e7c5-194" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(185.83,67.86)" id="img-8311e7c5-195" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(190.45,67.86)" id="img-8311e7c5-196" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(195.07,67.86)" id="img-8311e7c5-197" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(199.69,67.86)" id="img-8311e7c5-198" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(204.31,67.86)" id="img-8311e7c5-199" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(208.93,67.86)" id="img-8311e7c5-200" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(213.55,67.86)" id="img-8311e7c5-201" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-86.79,67.86)" id="img-8311e7c5-202" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(28.72,67.86)" id="img-8311e7c5-203" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(144.24,67.86)" id="img-8311e7c5-204" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(259.76,67.86)" id="img-8311e7c5-205" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-63.69,67.86)" id="img-8311e7c5-206" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-52.14,67.86)" id="img-8311e7c5-207" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-40.59,67.86)" id="img-8311e7c5-208" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-29.03,67.86)" id="img-8311e7c5-209" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-17.48,67.86)" id="img-8311e7c5-210" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-5.93,67.86)" id="img-8311e7c5-211" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(5.62,67.86)" id="img-8311e7c5-212" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(17.17,67.86)" id="img-8311e7c5-213" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(28.72,67.86)" id="img-8311e7c5-214" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(40.28,67.86)" id="img-8311e7c5-215" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(51.83,67.86)" id="img-8311e7c5-216" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(63.38,67.86)" id="img-8311e7c5-217" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(74.93,67.86)" id="img-8311e7c5-218" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(86.48,67.86)" id="img-8311e7c5-219" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(98.04,67.86)" id="img-8311e7c5-220" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(109.59,67.86)" id="img-8311e7c5-221" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(121.14,67.86)" id="img-8311e7c5-222" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(132.69,67.86)" id="img-8311e7c5-223" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(144.24,67.86)" id="img-8311e7c5-224" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(155.79,67.86)" id="img-8311e7c5-225" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(167.35,67.86)" id="img-8311e7c5-226" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(178.9,67.86)" id="img-8311e7c5-227" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(190.45,67.86)" id="img-8311e7c5-228" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(202,67.86)" id="img-8311e7c5-229" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(213.55,67.86)" id="img-8311e7c5-230" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
      </g>
      <g class="plotpanel" id="img-8311e7c5-231">
        <metadata>
          <boundingbox value="26.724999999999994mm 55.0mm 96.4138562373095mm 25.715000000000003mm"/>
          <unitbox value="-8.656710503949542 10921.022334791618 417.3134210078991 -11842.044669583238"/>
        </metadata>
        <g stroke-width="0.3" fill="#000000" fill-opacity="0.000" stroke-dasharray="none" id="img-8311e7c5-232">
          <g transform="translate(70.89,58.83)" id="img-8311e7c5-233" class="geometry color_S_V" stroke="#00BFFF">
            <path fill="none" d="M-42.16,-1.83 L -42.14 -1.83 -42.12 -1.82 -42.09 -1.82 -42.07 -1.82 -42.05 -1.82 -42.03 -1.82 -42 -1.82 -41.98 -1.82 -41.96 -1.82 -41.93 -1.82 -41.91 -1.82 -41.89 -1.81 -41.86 -1.81 -41.84 -1.81 -41.82 -1.81 -41.79 -1.81 -41.77 -1.81 -41.75 -1.81 -41.72 -1.81 -41.7 -1.81 -41.68 -1.81 -41.66 -1.81 -41.63 -1.81 -41.61 -1.8 -41.59 -1.8 -41.56 -1.8 -41.54 -1.8 -41.52 -1.8 -41.49 -1.8 -41.47 -1.8 -41.45 -1.8 -41.42 -1.8 -41.4 -1.8 -41.38 -1.8 -41.36 -1.8 -41.33 -1.79 -41.31 -1.79 -41.29 -1.79 -41.26 -1.79 -41.24 -1.79 -41.22 -1.79 -41.19 -1.79 -41.17 -1.79 -41.15 -1.79 -41.12 -1.79 -41.1 -1.79 -41.08 -1.79 -41.05 -1.79 -41.03 -1.78 -41.01 -1.78 -40.99 -1.78 -40.96 -1.78 -40.94 -1.78 -40.92 -1.78 -40.89 -1.78 -40.87 -1.78 -40.85 -1.78 -40.82 -1.78 -40.8 -1.78 -40.78 -1.78 -40.75 -1.78 -40.73 -1.77 -40.71 -1.77 -40.69 -1.77 -40.66 -1.77 -40.64 -1.77 -40.62 -1.77 -40.59 -1.77 -40.57 -1.77 -40.55 -1.77 -40.52 -1.77 -40.5 -1.77 -40.48 -1.76 -40.45 -1.76 -40.43 -1.76 -40.41 -1.76 -40.38 -1.76 -40.36 -1.76 -40.34 -1.76 -40.32 -1.76 -40.29 -1.76 -40.27 -1.75 -40.25 -1.75 -40.22 -1.75 -40.2 -1.75 -40.18 -1.75 -40.15 -1.75 -40.13 -1.75 -40.11 -1.75 -40.08 -1.74 -40.06 -1.74 -40.04 -1.74 -40.02 -1.74 -39.99 -1.74 -39.97 -1.74 -39.95 -1.74 -39.92 -1.73 -39.9 -1.73 -39.88 -1.73 -39.85 -1.73 -39.83 -1.73 -39.81 -1.73 -39.78 -1.72 -39.76 -1.72 -39.74 -1.72 -39.71 -1.72 -39.69 -1.72 -39.67 -1.72 -39.65 -1.71 -39.62 -1.71 -39.6 -1.71 -39.58 -1.71 -39.55 -1.7 -39.53 -1.7 -39.51 -1.7 -39.48 -1.7 -39.46 -1.7 -39.44 -1.69 -39.41 -1.69 -39.39 -1.69 -39.37 -1.69 -39.35 -1.68 -39.32 -1.68 -39.3 -1.68 -39.28 -1.67 -39.25 -1.67 -39.23 -1.67 -39.21 -1.67 -39.18 -1.66 -39.16 -1.66 -39.14 -1.66 -39.11 -1.65 -39.09 -1.65 -39.07 -1.65 -39.04 -1.64 -39.02 -1.64 -39 -1.64 -38.98 -1.63 -38.95 -1.63 -38.93 -1.63 -38.91 -1.62 -38.88 -1.62 -38.86 -1.61 -38.84 -1.61 -38.81 -1.61 -38.79 -1.6 -38.77 -1.6 -38.74 -1.59 -38.72 -1.59 -38.7 -1.58 -38.68 -1.58 -38.65 -1.58 -38.63 -1.57 -38.61 -1.57 -38.58 -1.56 -38.56 -1.56 -38.54 -1.55 -38.51 -1.55 -38.49 -1.54 -38.47 -1.54 -38.44 -1.53 -38.42 -1.52 -38.4 -1.52 -38.37 -1.51 -38.35 -1.51 -38.33 -1.5 -38.31 -1.5 -38.28 -1.49 -38.26 -1.48 -38.24 -1.48 -38.21 -1.47 -38.19 -1.46 -38.17 -1.46 -38.14 -1.45 -38.12 -1.44 -38.1 -1.44 -38.07 -1.43 -38.05 -1.42 -38.03 -1.41 -38.01 -1.41 -37.98 -1.4 -37.96 -1.39 -37.94 -1.38 -37.91 -1.37 -37.89 -1.37 -37.87 -1.36 -37.84 -1.35 -37.82 -1.34 -37.8 -1.33 -37.77 -1.32 -37.75 -1.31 -37.73 -1.3 -37.7 -1.3 -37.68 -1.29 -37.66 -1.28 -37.64 -1.27 -37.61 -1.26 -37.59 -1.25 -37.57 -1.24 -37.54 -1.23 -37.52 -1.21 -37.5 -1.2 -37.47 -1.19 -37.45 -1.18 -37.43 -1.17 -37.4 -1.16 -37.38 -1.15 -37.36 -1.14 -37.34 -1.12 -37.31 -1.11 -37.29 -1.1 -37.27 -1.09 -37.24 -1.07 -37.22 -1.06 -37.2 -1.05 -37.17 -1.03 -37.15 -1.02 -37.13 -1.01 -37.1 -0.99 -37.08 -0.98 -37.06 -0.96 -37.03 -0.95 -37.01 -0.93 -36.99 -0.92 -36.97 -0.9 -36.94 -0.89 -36.92 -0.87 -36.9 -0.86 -36.87 -0.84 -36.85 -0.82 -36.83 -0.81 -36.8 -0.79 -36.78 -0.77 -36.76 -0.76 -36.73 -0.74 -36.71 -0.72 -36.69 -0.7 -36.67 -0.69 -36.64 -0.67 -36.62 -0.65 -36.6 -0.63 -36.57 -0.61 -36.55 -0.59 -36.53 -0.57 -36.5 -0.55 -36.48 -0.53 -36.46 -0.51 -36.43 -0.49 -36.41 -0.47 -36.39 -0.45 -36.36 -0.43 -36.34 -0.41 -36.32 -0.39 -36.3 -0.37 -36.27 -0.34 -36.25 -0.32 -36.23 -0.3 -36.2 -0.28 -36.18 -0.25 -36.16 -0.23 -36.13 -0.21 -36.11 -0.18 -36.09 -0.16 -36.06 -0.14 -36.04 -0.11 -36.02 -0.09 -36 -0.06 -35.97 -0.04 -35.95 -0.01 -35.93 0.01 -35.9 0.04 -35.88 0.06 -35.86 0.09 -35.83 0.12 -35.81 0.14 -35.79 0.17 -35.76 0.2 -35.74 0.22 -35.72 0.25 -35.69 0.28 -35.67 0.31 -35.65 0.33 -35.63 0.36 -35.6 0.39 -35.58 0.42 -35.56 0.45 -35.53 0.47 -35.51 0.5 -35.49 0.53 -35.46 0.56 -35.44 0.59 -35.42 0.62 -35.39 0.65 -35.37 0.68 -35.35 0.71 -35.33 0.74 -35.3 0.77 -35.28 0.8 -35.26 0.83 -35.23 0.86 -35.21 0.89 -35.19 0.92 -35.16 0.95 -35.14 0.98 -35.12 1.01 -35.09 1.04 -35.07 1.07 -35.05 1.11 -35.02 1.14 -35 1.17 -34.98 1.2 -34.96 1.23 -34.93 1.26 -34.91 1.29 -34.89 1.32 -34.86 1.35 -34.84 1.39 -34.82 1.42 -34.79 1.45 -34.77 1.48 -34.75 1.51 -34.72 1.54 -34.7 1.57 -34.68 1.6 -34.66 1.63 -34.63 1.67 -34.61 1.7 -34.59 1.73 -34.56 1.76 -34.54 1.79 -34.52 1.82 -34.49 1.85 -34.47 1.88 -34.45 1.91 -34.42 1.94 -34.4 1.97 -34.38 2 -34.35 2.03 -34.33 2.06 -34.31 2.09 -34.29 2.12 -34.26 2.15 -34.24 2.18 -34.22 2.21 -34.19 2.24 -34.17 2.27 -34.15 2.29 -34.12 2.32 -34.1 2.35 -34.08 2.38 -34.05 2.41 -34.03 2.43 -34.01 2.46 -33.99 2.49 -33.96 2.52 -33.94 2.54 -33.92 2.57 -33.89 2.6 -33.87 2.62 -33.85 2.65 -33.82 2.68 -33.8 2.7 -33.78 2.73 -33.75 2.75 -33.73 2.78 -33.71 2.8 -33.68 2.83 -33.66 2.85 -33.64 2.87 -33.62 2.9 -33.59 2.92 -33.57 2.95 -33.55 2.97 -33.52 2.99 -33.5 3.01 -33.48 3.04 -33.45 3.06 -33.43 3.08 -33.41 3.1 -33.38 3.13 -33.36 3.15 -33.34 3.17 -33.32 3.19 -33.29 3.21 -33.27 3.23 -33.25 3.25 -33.22 3.27 -33.2 3.29 -33.18 3.31 -33.15 3.33 -33.13 3.35 -33.11 3.36 -33.08 3.38 -33.06 3.4 -33.04 3.42 -33.01 3.44 -32.99 3.45 -32.97 3.47 -32.95 3.49 -32.92 3.5 -32.9 3.52 -32.88 3.53 -32.85 3.55 -32.83 3.57 -32.81 3.58 -32.78 3.6 -32.76 3.61 -32.74 3.63 -32.71 3.64 -32.69 3.65 -32.67 3.67 -32.65 3.68 -32.62 3.7 -32.6 3.71 -32.58 3.72 -32.55 3.73 -32.53 3.75 -32.51 3.76 -32.48 3.77 -32.46 3.78 -32.44 3.79 -32.41 3.81 -32.39 3.82 -32.37 3.83 -32.34 3.84 -32.32 3.85 -32.3 3.86 -32.28 3.87 -32.25 3.88 -32.23 3.89 -32.21 3.9 -32.18 3.91 -32.16 3.92 -32.14 3.92 -32.11 3.93 -32.09 3.94 -32.07 3.95 -32.04 3.96 -32.02 3.96 -32 3.97 -31.98 3.98 -31.95 3.99 -31.93 3.99 -31.91 4 -31.88 4.01 -31.86 4.01 -31.84 4.02 -31.81 4.03 -31.79 4.03 -31.77 4.04 -31.74 4.04 -31.72 4.05 -31.7 4.05 -31.67 4.06 -31.65 4.06 -31.63 4.07 -31.61 4.07 -31.58 4.08 -31.56 4.08 -31.54 4.08 -31.51 4.09 -31.49 4.09 -31.47 4.1 -31.44 4.1 -31.42 4.1 -31.4 4.1 -31.37 4.11 -31.35 4.11 -31.33 4.11 -31.31 4.12 -31.28 4.12 -31.26 4.12 -31.24 4.12 -31.21 4.12 -31.19 4.13 -31.17 4.13 -31.14 4.13 -31.12 4.13 -31.1 4.13 -31.07 4.13 -31.05 4.13 -31.03 4.13 -31 4.13 -30.98 4.14 -30.96 4.14 -30.94 4.14 -30.91 4.14 -30.89 4.14 -30.87 4.14 -30.84 4.14 -30.82 4.14 -30.8 4.13 -30.77 4.13 -30.75 4.13 -30.73 4.13 -30.7 4.13 -30.68 4.13 -30.66 4.13 -30.64 4.13 -30.61 4.13 -30.59 4.13 -30.57 4.12 -30.54 4.12 -30.52 4.12 -30.5 4.12 -30.47 4.12 -30.45 4.11 -30.43 4.11 -30.4 4.11 -30.38 4.11 -30.36 4.11 -30.33 4.1 -30.31 4.1 -30.29 4.1 -30.27 4.09 -30.24 4.09 -30.22 4.09 -30.2 4.09 -30.17 4.08 -30.15 4.08 -30.13 4.08 -30.1 4.07 -30.08 4.07 -30.06 4.07 -30.03 4.06 -30.01 4.06 -29.99 4.05 -29.97 4.05 -29.94 4.05 -29.92 4.04 -29.9 4.04 -29.87 4.03 -29.85 4.03 -29.83 4.03 -29.8 4.02 -29.78 4.02 -29.76 4.01 -29.73 4.01 -29.71 4 -29.69 4 -29.66 3.99 -29.64 3.99 -29.62 3.98 -29.6 3.98 -29.57 3.97 -29.55 3.97 -29.53 3.96 -29.5 3.96 -29.48 3.95 -29.46 3.95 -29.43 3.94 -29.41 3.94 -29.39 3.93 -29.36 3.93 -29.34 3.92 -29.32 3.91 -29.3 3.91 -29.27 3.9 -29.25 3.9 -29.23 3.89 -29.2 3.88 -29.18 3.88 -29.16 3.87 -29.13 3.87 -29.11 3.86 -29.09 3.85 -29.06 3.85 -29.04 3.84 -29.02 3.83 -28.99 3.83 -28.97 3.82 -28.95 3.82 -28.93 3.81 -28.9 3.8 -28.88 3.8 -28.86 3.79 -28.83 3.78 -28.81 3.78 -28.79 3.77 -28.76 3.76 -28.74 3.76 -28.72 3.75 -28.69 3.74 -28.67 3.73 -28.65 3.73 -28.63 3.72 -28.6 3.71 -28.58 3.71 -28.56 3.7 -28.53 3.69 -28.51 3.69 -28.49 3.68 -28.46 3.67 -28.44 3.66 -28.42 3.66 -28.39 3.65 -28.37 3.64 -28.35 3.63 -28.32 3.63 -28.3 3.62 -28.28 3.61 -28.26 3.6 -28.23 3.6 -28.21 3.59 -28.19 3.58 -28.16 3.57 -28.14 3.57 -28.12 3.56 -28.09 3.55 -28.07 3.54 -28.05 3.54 -28.02 3.53 -28 3.52 -27.98 3.51 -27.96 3.5 -27.93 3.5 -27.91 3.49 -27.89 3.48 -27.86 3.47 -27.84 3.47 -27.82 3.46 -27.79 3.45 -27.77 3.44 -27.75 3.43 -27.72 3.43 -27.7 3.42 -27.68 3.41 -27.65 3.4 -27.63 3.39 -27.61 3.38 -27.59 3.38 -27.56 3.37 -27.54 3.36 -27.52 3.35 -27.49 3.34 -27.47 3.34 -27.45 3.33 -27.42 3.32 -27.4 3.31 -27.38 3.3 -27.35 3.3 -27.33 3.29 -27.31 3.28 -27.29 3.27 -27.26 3.26 -27.24 3.25 -27.22 3.25 -27.19 3.24 -27.17 3.23 -27.15 3.22 -27.12 3.21 -27.1 3.2 -27.08 3.2 -27.05 3.19 -27.03 3.18 -27.01 3.17 -26.98 3.16 -26.96 3.15 -26.94 3.15 -26.92 3.14 -26.89 3.13 -26.87 3.12 -26.85 3.11 -26.82 3.1 -26.8 3.09 -26.78 3.09 -26.75 3.08 -26.73 3.07 -26.71 3.06 -26.68 3.05 -26.66 3.04 -26.64 3.04 -26.62 3.03 -26.59 3.02 -26.57 3.01 -26.55 3 -26.52 2.99 -26.5 2.99 -26.48 2.98 -26.45 2.97 -26.43 2.96 -26.41 2.95 -26.38 2.94 -26.36 2.93 -26.34 2.93 -26.31 2.92 -26.29 2.91 -26.27 2.9 -26.25 2.89 -26.22 2.88 -26.2 2.87 -26.18 2.87 -26.15 2.86 -26.13 2.85 -26.11 2.84 -26.08 2.83 -26.06 2.82 -26.04 2.82 -26.01 2.81 -25.99 2.8 -25.97 2.79 -25.95 2.78 -25.92 2.77 -25.9 2.77 -25.88 2.76 -25.85 2.75 -25.83 2.74 -25.81 2.73 -25.78 2.72 -25.76 2.71 -25.74 2.71 -25.71 2.7 -25.69 2.69 -25.67 2.68 -25.64 2.67 -25.62 2.66 -25.6 2.66 -25.58 2.65 -25.55 2.64 -25.53 2.63 -25.51 2.62 -25.48 2.61 -25.46 2.61 -25.44 2.6 -25.41 2.59 -25.39 2.58 -25.37 2.57 -25.34 2.56 -25.32 2.56 -25.3 2.55 -25.28 2.54 -25.25 2.53 -25.23 2.52 -25.21 2.51 -25.18 2.51 -25.16 2.5 -25.14 2.49 -25.11 2.48 -25.09 2.47 -25.07 2.46 -25.04 2.46 -25.02 2.45 -25 2.44 -24.97 2.43 -24.95 2.42 -24.93 2.42 -24.91 2.41 -24.88 2.4 -24.86 2.39 -24.84 2.38 -24.81 2.38 -24.79 2.37 -24.77 2.36 -24.74 2.35 -24.72 2.34 -24.7 2.33 -24.67 2.33 -24.65 2.32 -24.63 2.31 -24.61 2.3 -24.58 2.29 -24.56 2.29 -24.54 2.28 -24.51 2.27 -24.49 2.26 -24.47 2.25 -24.44 2.25 -24.42 2.24 -24.4 2.23 -24.37 2.22 -24.35 2.21 -24.33 2.21 -24.3 2.2 -24.28 2.19 -24.26 2.18 -24.24 2.18 -24.21 2.17 -24.19 2.16 -24.17 2.15 -24.14 2.14 -24.12 2.14 -24.1 2.13 -24.07 2.12 -24.05 2.11 -24.03 2.11 -24 2.1 -23.98 2.09 -23.96 2.08 -23.94 2.07 -23.91 2.07 -23.89 2.06 -23.87 2.05 -23.84 2.04 -23.82 2.04 -23.8 2.03 -23.77 2.02 -23.75 2.01 -23.73 2.01 -23.7 2 -23.68 1.99 -23.66 1.98 -23.63 1.98 -23.61 1.97 -23.59 1.96 -23.57 1.95 -23.54 1.95 -23.52 1.94 -23.5 1.93 -23.47 1.92 -23.45 1.92 -23.43 1.91 -23.4 1.9 -23.38 1.89 -23.36 1.89 -23.33 1.88 -23.31 1.87 -23.29 1.86 -23.27 1.86 -23.24 1.85 -23.22 1.84 -23.2 1.83 -23.17 1.83 -23.15 1.82 -23.13 1.81 -23.1 1.81 -23.08 1.8 -23.06 1.79 -23.03 1.78 -23.01 1.78 -22.99 1.77 -22.96 1.76 -22.94 1.75 -22.92 1.75 -22.9 1.74 -22.87 1.73 -22.85 1.73 -22.83 1.72 -22.8 1.71 -22.78 1.71 -22.76 1.7 -22.73 1.69 -22.71 1.68 -22.69 1.68 -22.66 1.67 -22.64 1.66 -22.62 1.66 -22.6 1.65 -22.57 1.64 -22.55 1.64 -22.53 1.63 -22.5 1.62 -22.48 1.62 -22.46 1.61 -22.43 1.6 -22.41 1.59 -22.39 1.59 -22.36 1.58 -22.34 1.57 -22.32 1.57 -22.29 1.56 -22.27 1.55 -22.25 1.55 -22.23 1.54 -22.2 1.53 -22.18 1.53 -22.16 1.52 -22.13 1.51 -22.11 1.51 -22.09 1.5 -22.06 1.49 -22.04 1.49 -22.02 1.48 -21.99 1.47 -21.97 1.47 -21.95 1.46 -21.93 1.46 -21.9 1.45 -21.88 1.44 -21.86 1.44 -21.83 1.43 -21.81 1.42 -21.79 1.42 -21.76 1.41 -21.74 1.4 -21.72 1.4 -21.69 1.39 -21.67 1.38 -21.65 1.38 -21.62 1.37 -21.6 1.37 -21.58 1.36 -21.56 1.35 -21.53 1.35 -21.51 1.34 -21.49 1.33 -21.46 1.33 -21.44 1.32 -21.42 1.32 -21.39 1.31 -21.37 1.3 -21.35 1.3 -21.32 1.29 -21.3 1.28 -21.28 1.28 -21.26 1.27 -21.23 1.27 -21.21 1.26 -21.19 1.25 -21.16 1.25 -21.14 1.24 -21.12 1.24 -21.09 1.23 -21.07 1.22 -21.05 1.22 -21.02 1.21 -21 1.21 -20.98 1.2 -20.95 1.2 -20.93 1.19 -20.91 1.18 -20.89 1.18 -20.86 1.17 -20.84 1.17 -20.82 1.16 -20.79 1.15 -20.77 1.15 -20.75 1.14 -20.72 1.14 -20.7 1.13 -20.68 1.13 -20.65 1.12 -20.63 1.11 -20.61 1.11 -20.59 1.1 -20.56 1.1 -20.54 1.09 -20.52 1.09 -20.49 1.08 -20.47 1.07 -20.45 1.07 -20.42 1.06 -20.4 1.06 -20.38 1.05 -20.35 1.05 -20.33 1.04 -20.31 1.04 -20.28 1.03 -20.26 1.02 -20.24 1.02 -20.22 1.01 -20.19 1.01 -20.17 1 -20.15 1 -20.12 0.99 -20.1 0.99 -20.08 0.98 -20.05 0.98 -20.03 0.97 -20.01 0.97 -19.98 0.96 -19.96 0.96 -19.94 0.95 -19.92 0.94 -19.89 0.94 -19.87 0.93 -19.85 0.93 -19.82 0.92 -19.8 0.92 -19.78 0.91 -19.75 0.91 -19.73 0.9 -19.71 0.9 -19.68 0.89 -19.66 0.89 -19.64 0.88 -19.61 0.88 -19.59 0.87 -19.57 0.87 -19.55 0.86 -19.52 0.86 -19.5 0.85 -19.48 0.85 -19.45 0.84 -19.43 0.84 -19.41 0.83 -19.38 0.83 -19.36 0.82 -19.34 0.82 -19.31 0.81 -19.29 0.81 -19.27 0.8 -19.25 0.8 -19.22 0.79 -19.2 0.79 -19.18 0.78 -19.15 0.78 -19.13 0.77 -19.11 0.77 -19.08 0.76 -19.06 0.76 -19.04 0.76 -19.01 0.75 -18.99 0.75 -18.97 0.74 -18.94 0.74 -18.92 0.73 -18.9 0.73 -18.88 0.72 -18.85 0.72 -18.83 0.71 -18.81 0.71 -18.78 0.7 -18.76 0.7 -18.74 0.69 -18.71 0.69 -18.69 0.69 -18.67 0.68 -18.64 0.68 -18.62 0.67 -18.6 0.67 -18.58 0.66 -18.55 0.66 -18.53 0.65 -18.51 0.65 -18.48 0.65 -18.46 0.64 -18.44 0.64 -18.41 0.63 -18.39 0.63 -18.37 0.62 -18.34 0.62 -18.32 0.61 -18.3 0.61 -18.27 0.61 -18.25 0.6 -18.23 0.6 -18.21 0.59 -18.18 0.59 -18.16 0.58 -18.14 0.58 -18.11 0.58 -18.09 0.57 -18.07 0.57 -18.04 0.56 -18.02 0.56 -18 0.55 -17.97 0.55 -17.95 0.55 -17.93 0.54 -17.91 0.54 -17.88 0.53 -17.86 0.53 -17.84 0.53 -17.81 0.52 -17.79 0.52 -17.77 0.51 -17.74 0.51 -17.72 0.51 -17.7 0.5 -17.67 0.5 -17.65 0.49 -17.63 0.49 -17.6 0.49 -17.58 0.48 -17.56 0.48 -17.54 0.47 -17.51 0.47 -17.49 0.47 -17.47 0.46 -17.44 0.46 -17.42 0.45 -17.4 0.45 -17.37 0.45 -17.35 0.44 -17.33 0.44 -17.3 0.44 -17.28 0.43 -17.26 0.43 -17.24 0.42 -17.21 0.42 -17.19 0.42 -17.17 0.41 -17.14 0.41 -17.12 0.4 -17.1 0.4 -17.07 0.4 -17.05 0.39 -17.03 0.39 -17 0.39 -16.98 0.38 -16.96 0.38 -16.93 0.38 -16.91 0.37 -16.89 0.37 -16.87 0.36 -16.84 0.36 -16.82 0.36 -16.8 0.35 -16.77 0.35 -16.75 0.35 -16.73 0.34 -16.7 0.34 -16.68 0.34 -16.66 0.33 -16.63 0.33 -16.61 0.33 -16.59 0.32 -16.57 0.32 -16.54 0.31 -16.52 0.31 -16.5 0.31 -16.47 0.3 -16.45 0.3 -16.43 0.3 -16.4 0.29 -16.38 0.29 -16.36 0.29 -16.33 0.28 -16.31 0.28 -16.29 0.28 -16.26 0.27 -16.24 0.27 -16.22 0.27 -16.2 0.26 -16.17 0.26 -16.15 0.26 -16.13 0.25 -16.1 0.25 -16.08 0.25 -16.06 0.24 -16.03 0.24 -16.01 0.24 -15.99 0.23 -15.96 0.23 -15.94 0.23 -15.92 0.23 -15.9 0.22 -15.87 0.22 -15.85 0.22 -15.83 0.21 -15.8 0.21 -15.78 0.21 -15.76 0.2 -15.73 0.2 -15.71 0.2 -15.69 0.19 -15.66 0.19 -15.64 0.19 -15.62 0.18 -15.59 0.18 -15.57 0.18 -15.55 0.18 -15.53 0.17 -15.5 0.17 -15.48 0.17 -15.46 0.16 -15.43 0.16 -15.41 0.16 -15.39 0.15 -15.36 0.15 -15.34 0.15 -15.32 0.15 -15.29 0.14 -15.27 0.14 -15.25 0.14 -15.23 0.13 -15.2 0.13 -15.18 0.13 -15.16 0.12 -15.13 0.12 -15.11 0.12 -15.09 0.12 -15.06 0.11 -15.04 0.11 -15.02 0.11 -14.99 0.11 -14.97 0.1 -14.95 0.1 -14.92 0.1 -14.9 0.09 -14.88 0.09 -14.86 0.09 -14.83 0.09 -14.81 0.08 -14.79 0.08 -14.76 0.08 -14.74 0.07 -14.72 0.07 -14.69 0.07 -14.67 0.07 -14.65 0.06 -14.62 0.06 -14.6 0.06 -14.58 0.06 -14.56 0.05 -14.53 0.05 -14.51 0.05 -14.49 0.05 -14.46 0.04 -14.44 0.04 -14.42 0.04 -14.39 0.03 -14.37 0.03 -14.35 0.03 -14.32 0.03 -14.3 0.02 -14.28 0.02 -14.25 0.02 -14.23 0.02 -14.21 0.01 -14.19 0.01 -14.16 0.01 -14.14 0.01 -14.12 0 -14.09 0 -14.07 -0 -14.05 -0 -14.02 -0.01 -14 -0.01 -13.98 -0.01 -13.95 -0.01 -13.93 -0.02 -13.91 -0.02 -13.89 -0.02 -13.86 -0.02 -13.84 -0.03 -13.82 -0.03 -13.79 -0.03 -13.77 -0.03 -13.75 -0.04 -13.72 -0.04 -13.7 -0.04 -13.68 -0.04 -13.65 -0.04 -13.63 -0.05 -13.61 -0.05 -13.58 -0.05 -13.56 -0.05 -13.54 -0.06 -13.52 -0.06 -13.49 -0.06 -13.47 -0.06 -13.45 -0.07 -13.42 -0.07 -13.4 -0.07 -13.38 -0.07 -13.35 -0.07 -13.33 -0.08 -13.31 -0.08 -13.28 -0.08 -13.26 -0.08 -13.24 -0.09 -13.22 -0.09 -13.19 -0.09 -13.17 -0.09 -13.15 -0.09 -13.12 -0.1 -13.1 -0.1 -13.08 -0.1 -13.05 -0.1 -13.03 -0.11 -13.01 -0.11 -12.98 -0.11 -12.96 -0.11 -12.94 -0.11 -12.91 -0.12 -12.89 -0.12 -12.87 -0.12 -12.85 -0.12 -12.82 -0.12 -12.8 -0.13 -12.78 -0.13 -12.75 -0.13 -12.73 -0.13 -12.71 -0.13 -12.68 -0.14 -12.66 -0.14 -12.64 -0.14 -12.61 -0.14 -12.59 -0.15 -12.57 -0.15 -12.55 -0.15 -12.52 -0.15 -12.5 -0.15 -12.48 -0.16 -12.45 -0.16 -12.43 -0.16 -12.41 -0.16 -12.38 -0.16 -12.36 -0.17 -12.34 -0.17 -12.31 -0.17 -12.29 -0.17 -12.27 -0.17 -12.24 -0.17 -12.22 -0.18 -12.2 -0.18 -12.18 -0.18 -12.15 -0.18 -12.13 -0.18 -12.11 -0.19 -12.08 -0.19 -12.06 -0.19 -12.04 -0.19 -12.01 -0.19 -11.99 -0.2 -11.97 -0.2 -11.94 -0.2 -11.92 -0.2 -11.9 -0.2 -11.88 -0.21 -11.85 -0.21 -11.83 -0.21 -11.81 -0.21 -11.78 -0.21 -11.76 -0.21 -11.74 -0.22 -11.71 -0.22 -11.69 -0.22 -11.67 -0.22 -11.64 -0.22 -11.62 -0.22 -11.6 -0.23 -11.57 -0.23 -11.55 -0.23 -11.53 -0.23 -11.51 -0.23 -11.48 -0.24 -11.46 -0.24 -11.44 -0.24 -11.41 -0.24 -11.39 -0.24 -11.37 -0.24 -11.34 -0.25 -11.32 -0.25 -11.3 -0.25 -11.27 -0.25 -11.25 -0.25 -11.23 -0.25 -11.21 -0.26 -11.18 -0.26 -11.16 -0.26 -11.14 -0.26 -11.11 -0.26 -11.09 -0.26 -11.07 -0.27 -11.04 -0.27 -11.02 -0.27 -11 -0.27 -10.97 -0.27 -10.95 -0.27 -10.93 -0.28 -10.9 -0.28 -10.88 -0.28 -10.86 -0.28 -10.84 -0.28 -10.81 -0.28 -10.79 -0.28 -10.77 -0.29 -10.74 -0.29 -10.72 -0.29 -10.7 -0.29 -10.67 -0.29 -10.65 -0.29 -10.63 -0.3 -10.6 -0.3 -10.58 -0.3 -10.56 -0.3 -10.54 -0.3 -10.51 -0.3 -10.49 -0.3 -10.47 -0.31 -10.44 -0.31 -10.42 -0.31 -10.4 -0.31 -10.37 -0.31 -10.35 -0.31 -10.33 -0.32 -10.3 -0.32 -10.28 -0.32 -10.26 -0.32 -10.23 -0.32 -10.21 -0.32 -10.19 -0.32 -10.17 -0.33 -10.14 -0.33 -10.12 -0.33 -10.1 -0.33 -10.07 -0.33 -10.05 -0.33 -10.03 -0.33 -10 -0.34 -9.98 -0.34 -9.96 -0.34 -9.93 -0.34 -9.91 -0.34 -9.89 -0.34 -9.87 -0.34 -9.84 -0.35 -9.82 -0.35 -9.8 -0.35 -9.77 -0.35 -9.75 -0.35 -9.73 -0.35 -9.7 -0.35 -9.68 -0.35 -9.66 -0.36 -9.63 -0.36 -9.61 -0.36 -9.59 -0.36 -9.56 -0.36 -9.54 -0.36 -9.52 -0.36 -9.5 -0.37 -9.47 -0.37 -9.45 -0.37 -9.43 -0.37 -9.4 -0.37 -9.38 -0.37 -9.36 -0.37 -9.33 -0.37 -9.31 -0.38 -9.29 -0.38 -9.26 -0.38 -9.24 -0.38 -9.22 -0.38 -9.2 -0.38 -9.17 -0.38 -9.15 -0.38 -9.13 -0.39 -9.1 -0.39 -9.08 -0.39 -9.06 -0.39 -9.03 -0.39 -9.01 -0.39 -8.99 -0.39 -8.96 -0.39 -8.94 -0.4 -8.92 -0.4 -8.89 -0.4 -8.87 -0.4 -8.85 -0.4 -8.83 -0.4 -8.8 -0.4 -8.78 -0.4 -8.76 -0.4 -8.73 -0.41 -8.71 -0.41 -8.69 -0.41 -8.66 -0.41 -8.64 -0.41 -8.62 -0.41 -8.59 -0.41 -8.57 -0.41 -8.55 -0.42 -8.53 -0.42 -8.5 -0.42 -8.48 -0.42 -8.46 -0.42 -8.43 -0.42 -8.41 -0.42 -8.39 -0.42 -8.36 -0.42 -8.34 -0.43 -8.32 -0.43 -8.29 -0.43 -8.27 -0.43 -8.25 -0.43 -8.22 -0.43 -8.2 -0.43 -8.18 -0.43 -8.16 -0.43 -8.13 -0.43 -8.11 -0.44 -8.09 -0.44 -8.06 -0.44 -8.04 -0.44 -8.02 -0.44 -7.99 -0.44 -7.97 -0.44 -7.95 -0.44 -7.92 -0.44 -7.9 -0.45 -7.88 -0.45 -7.86 -0.45 -7.83 -0.45 -7.81 -0.45 -7.79 -0.45 -7.76 -0.45 -7.74 -0.45 -7.72 -0.45 -7.69 -0.45 -7.67 -0.46 -7.65 -0.46 -7.62 -0.46 -7.6 -0.46 -7.58 -0.46 -7.55 -0.46 -7.53 -0.46 -7.51 -0.46 -7.49 -0.46 -7.46 -0.46 -7.44 -0.47 -7.42 -0.47 -7.39 -0.47 -7.37 -0.47 -7.35 -0.47 -7.32 -0.47 -7.3 -0.47 -7.28 -0.47 -7.25 -0.47 -7.23 -0.47 -7.21 -0.47 -7.19 -0.48 -7.16 -0.48 -7.14 -0.48 -7.12 -0.48 -7.09 -0.48 -7.07 -0.48 -7.05 -0.48 -7.02 -0.48 -7 -0.48 -6.98 -0.48 -6.95 -0.48 -6.93 -0.49 -6.91 -0.49 -6.88 -0.49 -6.86 -0.49 -6.84 -0.49 -6.82 -0.49 -6.79 -0.49 -6.77 -0.49 -6.75 -0.49 -6.72 -0.49 -6.7 -0.49 -6.68 -0.5 -6.65 -0.5 -6.63 -0.5 -6.61 -0.5 -6.58 -0.5 -6.56 -0.5 -6.54 -0.5 -6.52 -0.5 -6.49 -0.5 -6.47 -0.5 -6.45 -0.5 -6.42 -0.5 -6.4 -0.51 -6.38 -0.51 -6.35 -0.51 -6.33 -0.51 -6.31 -0.51 -6.28 -0.51 -6.26 -0.51 -6.24 -0.51 -6.21 -0.51 -6.19 -0.51 -6.17 -0.51 -6.15 -0.51 -6.12 -0.52 -6.1 -0.52 -6.08 -0.52 -6.05 -0.52 -6.03 -0.52 -6.01 -0.52 -5.98 -0.52 -5.96 -0.52 -5.94 -0.52 -5.91 -0.52 -5.89 -0.52 -5.87 -0.52 -5.85 -0.52 -5.82 -0.53 -5.8 -0.53 -5.78 -0.53 -5.75 -0.53 -5.73 -0.53 -5.71 -0.53 -5.68 -0.53 -5.66 -0.53 -5.64 -0.53 -5.61 -0.53 -5.59 -0.53 -5.57 -0.53 -5.54 -0.53 -5.52 -0.53 -5.5 -0.54 -5.48 -0.54 -5.45 -0.54 -5.43 -0.54 -5.41 -0.54 -5.38 -0.54 -5.36 -0.54 -5.34 -0.54 -5.31 -0.54 -5.29 -0.54 -5.27 -0.54 -5.24 -0.54 -5.22 -0.54 -5.2 -0.54 -5.18 -0.55 -5.15 -0.55 -5.13 -0.55 -5.11 -0.55 -5.08 -0.55 -5.06 -0.55 -5.04 -0.55 -5.01 -0.55 -4.99 -0.55 -4.97 -0.55 -4.94 -0.55 -4.92 -0.55 -4.9 -0.55 -4.87 -0.55 -4.85 -0.55 -4.83 -0.56 -4.81 -0.56 -4.78 -0.56 -4.76 -0.56 -4.74 -0.56 -4.71 -0.56 -4.69 -0.56 -4.67 -0.56 -4.64 -0.56 -4.62 -0.56 -4.6 -0.56 -4.57 -0.56 -4.55 -0.56 -4.53 -0.56 -4.51 -0.56 -4.48 -0.56 -4.46 -0.57 -4.44 -0.57 -4.41 -0.57 -4.39 -0.57 -4.37 -0.57 -4.34 -0.57 -4.32 -0.57 -4.3 -0.57 -4.27 -0.57 -4.25 -0.57 -4.23 -0.57 -4.2 -0.57 -4.18 -0.57 -4.16 -0.57 -4.14 -0.57 -4.11 -0.57 -4.09 -0.57 -4.07 -0.58 -4.04 -0.58 -4.02 -0.58 -4 -0.58 -3.97 -0.58 -3.95 -0.58 -3.93 -0.58 -3.9 -0.58 -3.88 -0.58 -3.86 -0.58 -3.84 -0.58 -3.81 -0.58 -3.79 -0.58 -3.77 -0.58 -3.74 -0.58 -3.72 -0.58 -3.7 -0.58 -3.67 -0.58 -3.65 -0.59 -3.63 -0.59 -3.6 -0.59 -3.58 -0.59 -3.56 -0.59 -3.53 -0.59 -3.51 -0.59 -3.49 -0.59 -3.47 -0.59 -3.44 -0.59 -3.42 -0.59 -3.4 -0.59 -3.37 -0.59 -3.35 -0.59 -3.33 -0.59 -3.3 -0.59 -3.28 -0.59 -3.26 -0.59 -3.23 -0.59 -3.21 -0.6 -3.19 -0.6 -3.17 -0.6 -3.14 -0.6 -3.12 -0.6 -3.1 -0.6 -3.07 -0.6 -3.05 -0.6 -3.03 -0.6 -3 -0.6 -2.98 -0.6 -2.96 -0.6 -2.93 -0.6 -2.91 -0.6 -2.89 -0.6 -2.86 -0.6 -2.84 -0.6 -2.82 -0.6 -2.8 -0.6 -2.77 -0.6 -2.75 -0.6 -2.73 -0.61 -2.7 -0.61 -2.68 -0.61 -2.66 -0.61 -2.63 -0.61 -2.61 -0.61 -2.59 -0.61 -2.56 -0.61 -2.54 -0.61 -2.52 -0.61 -2.5 -0.61 -2.47 -0.61 -2.45 -0.61 -2.43 -0.61 -2.4 -0.61 -2.38 -0.61 -2.36 -0.61 -2.33 -0.61 -2.31 -0.61 -2.29 -0.61 -2.26 -0.61 -2.24 -0.61 -2.22 -0.61 -2.19 -0.62 -2.17 -0.62 -2.15 -0.62 -2.13 -0.62 -2.1 -0.62 -2.08 -0.62 -2.06 -0.62 -2.03 -0.62 -2.01 -0.62 -1.99 -0.62 -1.96 -0.62 -1.94 -0.62 -1.92 -0.62 -1.89 -0.62 -1.87 -0.62 -1.85 -0.62 -1.83 -0.62 -1.8 -0.62 -1.78 -0.62 -1.76 -0.62 -1.73 -0.62 -1.71 -0.62 -1.69 -0.62 -1.66 -0.62 -1.64 -0.63 -1.62 -0.63 -1.59 -0.63 -1.57 -0.63 -1.55 -0.63 -1.52 -0.63 -1.5 -0.63 -1.48 -0.63 -1.46 -0.63 -1.43 -0.63 -1.41 -0.63 -1.39 -0.63 -1.36 -0.63 -1.34 -0.63 -1.32 -0.63 -1.29 -0.63 -1.27 -0.63 -1.25 -0.63 -1.22 -0.63 -1.2 -0.63 -1.18 -0.63 -1.16 -0.63 -1.13 -0.63 -1.11 -0.63 -1.09 -0.63 -1.06 -0.63 -1.04 -0.63 -1.02 -0.64 -0.99 -0.64 -0.97 -0.64 -0.95 -0.64 -0.92 -0.64 -0.9 -0.64 -0.88 -0.64 -0.85 -0.64 -0.83 -0.64 -0.81 -0.64 -0.79 -0.64 -0.76 -0.64 -0.74 -0.64 -0.72 -0.64 -0.69 -0.64 -0.67 -0.64 -0.65 -0.64 -0.62 -0.64 -0.6 -0.64 -0.58 -0.64 -0.55 -0.64 -0.53 -0.64 -0.51 -0.64 -0.49 -0.64 -0.46 -0.64 -0.44 -0.64 -0.42 -0.64 -0.39 -0.64 -0.37 -0.64 -0.35 -0.64 -0.32 -0.64 -0.3 -0.65 -0.28 -0.65 -0.25 -0.65 -0.23 -0.65 -0.21 -0.65 -0.18 -0.65 -0.16 -0.65 -0.14 -0.65 -0.12 -0.65 -0.09 -0.65 -0.07 -0.65 -0.05 -0.65 -0.02 -0.65 0 -0.65 0.02 -0.65 0.05 -0.65 0.07 -0.65 0.09 -0.65 0.12 -0.65 0.14 -0.65 0.16 -0.65 0.18 -0.65 0.21 -0.65 0.23 -0.65 0.25 -0.65 0.28 -0.65 0.3 -0.65 0.32 -0.65 0.35 -0.65 0.37 -0.65 0.39 -0.65 0.42 -0.65 0.44 -0.65 0.46 -0.65 0.49 -0.65 0.51 -0.66 0.53 -0.66 0.55 -0.66 0.58 -0.66 0.6 -0.66 0.62 -0.66 0.65 -0.66 0.67 -0.66 0.69 -0.66 0.72 -0.66 0.74 -0.66 0.76 -0.66 0.79 -0.66 0.81 -0.66 0.83 -0.66 0.85 -0.66 0.88 -0.66 0.9 -0.66 0.92 -0.66 0.95 -0.66 0.97 -0.66 0.99 -0.66 1.02 -0.66 1.04 -0.66 1.06 -0.66 1.09 -0.66 1.11 -0.66 1.13 -0.66 1.16 -0.66 1.18 -0.66 1.2 -0.66 1.22 -0.66 1.25 -0.66 1.27 -0.66 1.29 -0.66 1.32 -0.66 1.34 -0.66 1.36 -0.66 1.39 -0.66 1.41 -0.66 1.43 -0.67 1.46 -0.67 1.48 -0.67 1.5 -0.67 1.52 -0.67 1.55 -0.67 1.57 -0.67 1.59 -0.67 1.62 -0.67 1.64 -0.67 1.66 -0.67 1.69 -0.67 1.71 -0.67 1.73 -0.67 1.76 -0.67 1.78 -0.67 1.8 -0.67 1.83 -0.67 1.85 -0.67 1.87 -0.67 1.89 -0.67 1.92 -0.67 1.94 -0.67 1.96 -0.67 1.99 -0.67 2.01 -0.67 2.03 -0.67 2.06 -0.67 2.08 -0.67 2.1 -0.67 2.13 -0.67 2.15 -0.67 2.17 -0.67 2.19 -0.67 2.22 -0.67 2.24 -0.67 2.26 -0.67 2.29 -0.67 2.31 -0.67 2.33 -0.67 2.36 -0.67 2.38 -0.67 2.4 -0.67 2.43 -0.67 2.45 -0.67 2.47 -0.67 2.5 -0.67 2.52 -0.67 2.54 -0.67 2.56 -0.68 2.59 -0.68 2.61 -0.68 2.63 -0.68 2.66 -0.68 2.68 -0.68 2.7 -0.68 2.73 -0.68 2.75 -0.68 2.77 -0.68 2.8 -0.68 2.82 -0.68 2.84 -0.68 2.86 -0.68 2.89 -0.68 2.91 -0.68 2.93 -0.68 2.96 -0.68 2.98 -0.68 3 -0.68 3.03 -0.68 3.05 -0.68 3.07 -0.68 3.1 -0.68 3.12 -0.68 3.14 -0.68 3.17 -0.68 3.19 -0.68 3.21 -0.68 3.23 -0.68 3.26 -0.68 3.28 -0.68 3.3 -0.68 3.33 -0.68 3.35 -0.68 3.37 -0.68 3.4 -0.68 3.42 -0.68 3.44 -0.68 3.47 -0.68 3.49 -0.68 3.51 -0.68 3.53 -0.68 3.56 -0.68 3.58 -0.68 3.6 -0.68 3.63 -0.68 3.65 -0.68 3.67 -0.68 3.7 -0.68 3.72 -0.68 3.74 -0.68 3.77 -0.68 3.79 -0.68 3.81 -0.68 3.84 -0.68 3.86 -0.68 3.88 -0.68 3.9 -0.68 3.93 -0.68 3.95 -0.68 3.97 -0.69 4 -0.69 4.02 -0.69 4.04 -0.69 4.07 -0.69 4.09 -0.69 4.11 -0.69 4.14 -0.69 4.16 -0.69 4.18 -0.69 4.2 -0.69 4.23 -0.69 4.25 -0.69 4.27 -0.69 4.3 -0.69 4.32 -0.69 4.34 -0.69 4.37 -0.69 4.39 -0.69 4.41 -0.69 4.44 -0.69 4.46 -0.69 4.48 -0.69 4.51 -0.69 4.53 -0.69 4.55 -0.69 4.57 -0.69 4.6 -0.69 4.62 -0.69 4.64 -0.69 4.67 -0.69 4.69 -0.69 4.71 -0.69 4.74 -0.69 4.76 -0.69 4.78 -0.69 4.81 -0.69 4.83 -0.69 4.85 -0.69 4.87 -0.69 4.9 -0.69 4.92 -0.69 4.94 -0.69 4.97 -0.69 4.99 -0.69 5.01 -0.69 5.04 -0.69 5.06 -0.69 5.08 -0.69 5.11 -0.69 5.13 -0.69 5.15 -0.69 5.18 -0.69 5.2 -0.69 5.22 -0.69 5.24 -0.69 5.27 -0.69 5.29 -0.69 5.31 -0.69 5.34 -0.69 5.36 -0.69 5.38 -0.69 5.41 -0.69 5.43 -0.69 5.45 -0.69 5.48 -0.69 5.5 -0.69 5.52 -0.69 5.54 -0.69 5.57 -0.69 5.59 -0.69 5.61 -0.69 5.64 -0.69 5.66 -0.69 5.68 -0.69 5.71 -0.69 5.73 -0.69 5.75 -0.69 5.78 -0.69 5.8 -0.69 5.82 -0.69 5.85 -0.69 5.87 -0.69 5.89 -0.69 5.91 -0.69 5.94 -0.69 5.96 -0.7 5.98 -0.7 6.01 -0.7 6.03 -0.7 6.05 -0.7 6.08 -0.7 6.1 -0.7 6.12 -0.7 6.15 -0.7 6.17 -0.7 6.19 -0.7 6.21 -0.7 6.24 -0.7 6.26 -0.7 6.28 -0.7 6.31 -0.7 6.33 -0.7 6.35 -0.7 6.38 -0.7 6.4 -0.7 6.42 -0.7 6.45 -0.7 6.47 -0.7 6.49 -0.7 6.52 -0.7 6.54 -0.7 6.56 -0.7 6.58 -0.7 6.61 -0.7 6.63 -0.7 6.65 -0.7 6.68 -0.7 6.7 -0.7 6.72 -0.7 6.75 -0.7 6.77 -0.7 6.79 -0.7 6.82 -0.7 6.84 -0.7 6.86 -0.7 6.88 -0.7 6.91 -0.7 6.93 -0.7 6.95 -0.7 6.98 -0.7 7 -0.7 7.02 -0.7 7.05 -0.7 7.07 -0.7 7.09 -0.7 7.12 -0.7 7.14 -0.7 7.16 -0.7 7.19 -0.7 7.21 -0.7 7.23 -0.7 7.25 -0.7 7.28 -0.7 7.3 -0.7 7.32 -0.7 7.35 -0.7 7.37 -0.7 7.39 -0.7 7.42 -0.7 7.44 -0.7 7.46 -0.7 7.49 -0.7 7.51 -0.7 7.53 -0.7 7.55 -0.7 7.58 -0.7 7.6 -0.7 7.62 -0.7 7.65 -0.7 7.67 -0.7 7.69 -0.7 7.72 -0.7 7.74 -0.7 7.76 -0.7 7.79 -0.7 7.81 -0.7 7.83 -0.7 7.86 -0.7 7.88 -0.7 7.9 -0.7 7.92 -0.7 7.95 -0.7 7.97 -0.7 7.99 -0.7 8.02 -0.7 8.04 -0.7 8.06 -0.7 8.09 -0.7 8.11 -0.7 8.13 -0.7 8.16 -0.7 8.18 -0.7 8.2 -0.7 8.22 -0.7 8.25 -0.7 8.27 -0.7 8.29 -0.7 8.32 -0.7 8.34 -0.7 8.36 -0.7 8.39 -0.7 8.41 -0.7 8.43 -0.7 8.46 -0.7 8.48 -0.7 8.5 -0.7 8.53 -0.7 8.55 -0.7 8.57 -0.7 8.59 -0.7 8.62 -0.7 8.64 -0.7 8.66 -0.7 8.69 -0.7 8.71 -0.7 8.73 -0.7 8.76 -0.7 8.78 -0.7 8.8 -0.7 8.83 -0.7 8.85 -0.7 8.87 -0.7 8.89 -0.7 8.92 -0.7 8.94 -0.7 8.96 -0.7 8.99 -0.7 9.01 -0.7 9.03 -0.7 9.06 -0.7 9.08 -0.7 9.1 -0.7 9.13 -0.7 9.15 -0.7 9.17 -0.7 9.2 -0.7 9.22 -0.7 9.24 -0.7 9.26 -0.7 9.29 -0.7 9.31 -0.7 9.33 -0.7 9.36 -0.7 9.38 -0.7 9.4 -0.7 9.43 -0.7 9.45 -0.71 9.47 -0.71 9.5 -0.71 9.52 -0.71 9.54 -0.71 9.56 -0.71 9.59 -0.71 9.61 -0.71 9.63 -0.71 9.66 -0.71 9.68 -0.71 9.7 -0.71 9.73 -0.71 9.75 -0.71 9.77 -0.71 9.8 -0.71 9.82 -0.71 9.84 -0.71 9.87 -0.71 9.89 -0.71 9.91 -0.71 9.93 -0.71 9.96 -0.71 9.98 -0.71 10 -0.71 10.03 -0.71 10.05 -0.71 10.07 -0.71 10.1 -0.71 10.12 -0.71 10.14 -0.71 10.17 -0.71 10.19 -0.71 10.21 -0.71 10.23 -0.71 10.26 -0.71 10.28 -0.71 10.3 -0.71 10.33 -0.71 10.35 -0.71 10.37 -0.71 10.4 -0.71 10.42 -0.71 10.44 -0.71 10.47 -0.71 10.49 -0.71 10.51 -0.71 10.54 -0.71 10.56 -0.71 10.58 -0.71 10.6 -0.71 10.63 -0.71 10.65 -0.71 10.67 -0.71 10.7 -0.71 10.72 -0.71 10.74 -0.71 10.77 -0.71 10.79 -0.71 10.81 -0.71 10.84 -0.71 10.86 -0.71 10.88 -0.71 10.9 -0.71 10.93 -0.71 10.95 -0.71 10.97 -0.71 11 -0.71 11.02 -0.71 11.04 -0.71 11.07 -0.71 11.09 -0.71 11.11 -0.71 11.14 -0.71 11.16 -0.71 11.18 -0.71 11.21 -0.71 11.23 -0.71 11.25 -0.71 11.27 -0.71 11.3 -0.71 11.32 -0.71 11.34 -0.71 11.37 -0.71 11.39 -0.71 11.41 -0.71 11.44 -0.71 11.46 -0.71 11.48 -0.71 11.51 -0.71 11.53 -0.71 11.55 -0.71 11.57 -0.71 11.6 -0.71 11.62 -0.71 11.64 -0.71 11.67 -0.71 11.69 -0.71 11.71 -0.71 11.74 -0.71 11.76 -0.71 11.78 -0.71 11.81 -0.71 11.83 -0.71 11.85 -0.71 11.88 -0.71 11.9 -0.71 11.92 -0.71 11.94 -0.71 11.97 -0.71 11.99 -0.71 12.01 -0.71 12.04 -0.71 12.06 -0.71 12.08 -0.71 12.11 -0.71 12.13 -0.71 12.15 -0.71 12.18 -0.71 12.2 -0.71 12.22 -0.71 12.24 -0.71 12.27 -0.71 12.29 -0.71 12.31 -0.71 12.34 -0.71 12.36 -0.71 12.38 -0.71 12.41 -0.71 12.43 -0.71 12.45 -0.71 12.48 -0.71 12.5 -0.71 12.52 -0.71 12.55 -0.71 12.57 -0.71 12.59 -0.71 12.61 -0.71 12.64 -0.71 12.66 -0.71 12.68 -0.71 12.71 -0.71 12.73 -0.71 12.75 -0.71 12.78 -0.71 12.8 -0.71 12.82 -0.71 12.85 -0.71 12.87 -0.71 12.89 -0.71 12.91 -0.71 12.94 -0.71 12.96 -0.71 12.98 -0.71 13.01 -0.71 13.03 -0.71 13.05 -0.71 13.08 -0.71 13.1 -0.71 13.12 -0.71 13.15 -0.71 13.17 -0.71 13.19 -0.71 13.22 -0.71 13.24 -0.71 13.26 -0.71 13.28 -0.71 13.31 -0.71 13.33 -0.71 13.35 -0.71 13.38 -0.71 13.4 -0.71 13.42 -0.71 13.45 -0.71 13.47 -0.71 13.49 -0.71 13.52 -0.71 13.54 -0.71 13.56 -0.71 13.58 -0.71 13.61 -0.71 13.63 -0.71 13.65 -0.71 13.68 -0.71 13.7 -0.71 13.72 -0.71 13.75 -0.71 13.77 -0.71 13.79 -0.71 13.82 -0.71 13.84 -0.71 13.86 -0.71 13.89 -0.71 13.91 -0.71 13.93 -0.71 13.95 -0.71 13.98 -0.71 14 -0.71 14.02 -0.71 14.05 -0.71 14.07 -0.71 14.09 -0.71 14.12 -0.71 14.14 -0.71 14.16 -0.71 14.19 -0.71 14.21 -0.71 14.23 -0.71 14.25 -0.71 14.28 -0.71 14.3 -0.71 14.32 -0.71 14.35 -0.71 14.37 -0.71 14.39 -0.71 14.42 -0.71 14.44 -0.71 14.46 -0.71 14.49 -0.71 14.51 -0.71 14.53 -0.71 14.56 -0.71 14.58 -0.71 14.6 -0.71 14.62 -0.71 14.65 -0.71 14.67 -0.71 14.69 -0.71 14.72 -0.71 14.74 -0.71 14.76 -0.71 14.79 -0.71 14.81 -0.71 14.83 -0.71 14.86 -0.71 14.88 -0.71 14.9 -0.71 14.92 -0.71 14.95 -0.71 14.97 -0.71 14.99 -0.71 15.02 -0.71 15.04 -0.71 15.06 -0.71 15.09 -0.71 15.11 -0.71 15.13 -0.71 15.16 -0.71 15.18 -0.71 15.2 -0.71 15.23 -0.71 15.25 -0.71 15.27 -0.71 15.29 -0.71 15.32 -0.71 15.34 -0.71 15.36 -0.71 15.39 -0.71 15.41 -0.71 15.43 -0.71 15.46 -0.71 15.48 -0.71 15.5 -0.71 15.53 -0.71 15.55 -0.71 15.57 -0.71 15.59 -0.71 15.62 -0.71 15.64 -0.71 15.66 -0.71 15.69 -0.71 15.71 -0.71 15.73 -0.71 15.76 -0.71 15.78 -0.71 15.8 -0.71 15.83 -0.71 15.85 -0.71 15.87 -0.71 15.9 -0.71 15.92 -0.71 15.94 -0.71 15.96 -0.71 15.99 -0.71 16.01 -0.71 16.03 -0.71 16.06 -0.71 16.08 -0.71 16.1 -0.71 16.13 -0.71 16.15 -0.71 16.17 -0.71 16.2 -0.71 16.22 -0.71 16.24 -0.71 16.26 -0.71 16.29 -0.71 16.31 -0.71 16.33 -0.71 16.36 -0.71 16.38 -0.71 16.4 -0.71 16.43 -0.71 16.45 -0.71 16.47 -0.71 16.5 -0.71 16.52 -0.71 16.54 -0.71 16.57 -0.71 16.59 -0.71 16.61 -0.71 16.63 -0.71 16.66 -0.71 16.68 -0.71 16.7 -0.71 16.73 -0.71 16.75 -0.71 16.77 -0.71 16.8 -0.71 16.82 -0.71 16.84 -0.71 16.87 -0.71 16.89 -0.71 16.91 -0.71 16.93 -0.71 16.96 -0.71 16.98 -0.71 17 -0.71 17.03 -0.71 17.05 -0.71 17.07 -0.71 17.1 -0.71 17.12 -0.71 17.14 -0.71 17.17 -0.71 17.19 -0.71 17.21 -0.71 17.24 -0.71 17.26 -0.71 17.28 -0.71 17.3 -0.71 17.33 -0.71 17.35 -0.71 17.37 -0.71 17.4 -0.71 17.42 -0.71 17.44 -0.71 17.47 -0.71 17.49 -0.71 17.51 -0.71 17.54 -0.71 17.56 -0.71 17.58 -0.71 17.6 -0.71 17.63 -0.71 17.65 -0.71 17.67 -0.71 17.7 -0.71 17.72 -0.71 17.74 -0.71 17.77 -0.71 17.79 -0.71 17.81 -0.71 17.84 -0.71 17.86 -0.71 17.88 -0.71 17.91 -0.71 17.93 -0.71 17.95 -0.71 17.97 -0.71 18 -0.71 18.02 -0.71 18.04 -0.71 18.07 -0.71 18.09 -0.71 18.11 -0.71 18.14 -0.71 18.16 -0.71 18.18 -0.71 18.21 -0.71 18.23 -0.71 18.25 -0.71 18.27 -0.71 18.3 -0.71 18.32 -0.71 18.34 -0.71 18.37 -0.71 18.39 -0.71 18.41 -0.71 18.44 -0.71 18.46 -0.71 18.48 -0.71 18.51 -0.71 18.53 -0.71 18.55 -0.71 18.58 -0.71 18.6 -0.71 18.62 -0.71 18.64 -0.71 18.67 -0.71 18.69 -0.71 18.71 -0.71 18.74 -0.71 18.76 -0.71 18.78 -0.71 18.81 -0.71 18.83 -0.71 18.85 -0.71 18.88 -0.71 18.9 -0.71 18.92 -0.71 18.94 -0.71 18.97 -0.71 18.99 -0.71 19.01 -0.71 19.04 -0.71 19.06 -0.71 19.08 -0.71 19.11 -0.71 19.13 -0.71 19.15 -0.71 19.18 -0.71 19.2 -0.71 19.22 -0.71 19.25 -0.71 19.27 -0.71 19.29 -0.71 19.31 -0.71 19.34 -0.71 19.36 -0.71 19.38 -0.71 19.41 -0.71 19.43 -0.71 19.45 -0.71 19.48 -0.71 19.5 -0.71 19.52 -0.71 19.55 -0.71 19.57 -0.71 19.59 -0.71 19.61 -0.71 19.64 -0.71 19.66 -0.71 19.68 -0.71 19.71 -0.71 19.73 -0.71 19.75 -0.71 19.78 -0.71 19.8 -0.71 19.82 -0.71 19.85 -0.71 19.87 -0.71 19.89 -0.71 19.92 -0.71 19.94 -0.71 19.96 -0.71 19.98 -0.71 20.01 -0.71 20.03 -0.71 20.05 -0.71 20.08 -0.71 20.1 -0.71 20.12 -0.71 20.15 -0.71 20.17 -0.71 20.19 -0.71 20.22 -0.71 20.24 -0.71 20.26 -0.71 20.28 -0.71 20.31 -0.71 20.33 -0.71 20.35 -0.71 20.38 -0.71 20.4 -0.71 20.42 -0.71 20.45 -0.71 20.47 -0.71 20.49 -0.71 20.52 -0.71 20.54 -0.71 20.56 -0.71 20.59 -0.71 20.61 -0.71 20.63 -0.71 20.65 -0.71 20.68 -0.71 20.7 -0.71 20.72 -0.71 20.75 -0.71 20.77 -0.71 20.79 -0.71 20.82 -0.71 20.84 -0.71 20.86 -0.71 20.89 -0.71 20.91 -0.71 20.93 -0.71 20.95 -0.71 20.98 -0.71 21 -0.71 21.02 -0.71 21.05 -0.71 21.07 -0.71 21.09 -0.71 21.12 -0.71 21.14 -0.71 21.16 -0.71 21.19 -0.71 21.21 -0.71 21.23 -0.71 21.26 -0.71 21.28 -0.71 21.3 -0.71 21.32 -0.71 21.35 -0.71 21.37 -0.71 21.39 -0.71 21.42 -0.71 21.44 -0.71 21.46 -0.71 21.49 -0.71 21.51 -0.71 21.53 -0.71 21.56 -0.71 21.58 -0.71 21.6 -0.71 21.62 -0.71 21.65 -0.71 21.67 -0.71 21.69 -0.71 21.72 -0.71 21.74 -0.71 21.76 -0.71 21.79 -0.71 21.81 -0.71 21.83 -0.71 21.86 -0.71 21.88 -0.71 21.9 -0.71 21.93 -0.71 21.95 -0.71 21.97 -0.71 21.99 -0.71 22.02 -0.71 22.04 -0.71 22.06 -0.71 22.09 -0.71 22.11 -0.71 22.13 -0.71 22.16 -0.71 22.18 -0.71 22.2 -0.71 22.23 -0.71 22.25 -0.71 22.27 -0.71 22.29 -0.71 22.32 -0.71 22.34 -0.71 22.36 -0.71 22.39 -0.71 22.41 -0.71 22.43 -0.71 22.46 -0.71 22.48 -0.71 22.5 -0.71 22.53 -0.71 22.55 -0.71 22.57 -0.71 22.6 -0.71 22.62 -0.71 22.64 -0.71 22.66 -0.71 22.69 -0.71 22.71 -0.71 22.73 -0.71 22.76 -0.71 22.78 -0.71 22.8 -0.71 22.83 -0.71 22.85 -0.71 22.87 -0.71 22.9 -0.71 22.92 -0.71 22.94 -0.71 22.96 -0.71 22.99 -0.71 23.01 -0.71 23.03 -0.71 23.06 -0.71 23.08 -0.71 23.1 -0.71 23.13 -0.71 23.15 -0.71 23.17 -0.71 23.2 -0.71 23.22 -0.71 23.24 -0.71 23.27 -0.71 23.29 -0.71 23.31 -0.71 23.33 -0.71 23.36 -0.71 23.38 -0.71 23.4 -0.71 23.43 -0.71 23.45 -0.71 23.47 -0.71 23.5 -0.71 23.52 -0.71 23.54 -0.71 23.57 -0.71 23.59 -0.71 23.61 -0.71 23.63 -0.71 23.66 -0.71 23.68 -0.71 23.7 -0.71 23.73 -0.71 23.75 -0.71 23.77 -0.71 23.8 -0.71 23.82 -0.71 23.84 -0.71 23.87 -0.71 23.89 -0.71 23.91 -0.71 23.94 -0.71 23.96 -0.71 23.98 -0.71 24 -0.71 24.03 -0.71 24.05 -0.71 24.07 -0.71 24.1 -0.71 24.12 -0.71 24.14 -0.71 24.17 -0.71 24.19 -0.71 24.21 -0.71 24.24 -0.71 24.26 -0.71 24.28 -0.71 24.3 -0.71 24.33 -0.71 24.35 -0.71 24.37 -0.71 24.4 -0.71 24.42 -0.71 24.44 -0.71 24.47 -0.71 24.49 -0.71 24.51 -0.71 24.54 -0.71 24.56 -0.71 24.58 -0.71 24.61 -0.71 24.63 -0.71 24.65 -0.71 24.67 -0.71 24.7 -0.71 24.72 -0.71 24.74 -0.71 24.77 -0.71 24.79 -0.71 24.81 -0.71 24.84 -0.71 24.86 -0.71 24.88 -0.71 24.91 -0.71 24.93 -0.71 24.95 -0.71 24.97 -0.71 25 -0.71 25.02 -0.71 25.04 -0.71 25.07 -0.71 25.09 -0.71 25.11 -0.71 25.14 -0.71 25.16 -0.71 25.18 -0.71 25.21 -0.71 25.23 -0.71 25.25 -0.71 25.28 -0.71 25.3 -0.71 25.32 -0.71 25.34 -0.71 25.37 -0.71 25.39 -0.71 25.41 -0.71 25.44 -0.71 25.46 -0.71 25.48 -0.71 25.51 -0.71 25.53 -0.71 25.55 -0.71 25.58 -0.71 25.6 -0.71 25.62 -0.71 25.64 -0.71 25.67 -0.71 25.69 -0.71 25.71 -0.71 25.74 -0.71 25.76 -0.71 25.78 -0.71 25.81 -0.71 25.83 -0.71 25.85 -0.71 25.88 -0.71 25.9 -0.71 25.92 -0.71 25.95 -0.71 25.97 -0.71 25.99 -0.71 26.01 -0.71 26.04 -0.71 26.06 -0.71 26.08 -0.71 26.11 -0.71 26.13 -0.71 26.15 -0.71 26.18 -0.71 26.2 -0.71 26.22 -0.71 26.25 -0.71 26.27 -0.71 26.29 -0.71 26.31 -0.71 26.34 -0.71 26.36 -0.71 26.38 -0.71 26.41 -0.71 26.43 -0.71 26.45 -0.71 26.48 -0.71 26.5 -0.71 26.52 -0.71 26.55 -0.71 26.57 -0.71 26.59 -0.71 26.62 -0.71 26.64 -0.71 26.66 -0.71 26.68 -0.71 26.71 -0.71 26.73 -0.71 26.75 -0.71 26.78 -0.71 26.8 -0.71 26.82 -0.71 26.85 -0.71 26.87 -0.71 26.89 -0.71 26.92 -0.71 26.94 -0.71 26.96 -0.71 26.98 -0.71 27.01 -0.71 27.03 -0.71 27.05 -0.71 27.08 -0.71 27.1 -0.71 27.12 -0.71 27.15 -0.71 27.17 -0.71 27.19 -0.71 27.22 -0.71 27.24 -0.71 27.26 -0.71 27.29 -0.71 27.31 -0.71 27.33 -0.71 27.35 -0.71 27.38 -0.71 27.4 -0.71 27.42 -0.71 27.45 -0.71 27.47 -0.71 27.49 -0.71 27.52 -0.71 27.54 -0.71 27.56 -0.71 27.59 -0.71 27.61 -0.71 27.63 -0.71 27.65 -0.71 27.68 -0.71 27.7 -0.71 27.72 -0.71 27.75 -0.71 27.77 -0.71 27.79 -0.71 27.82 -0.71 27.84 -0.71 27.86 -0.71 27.89 -0.71 27.91 -0.71 27.93 -0.71 27.96 -0.71 27.98 -0.71 28 -0.71 28.02 -0.71 28.05 -0.71 28.07 -0.71 28.09 -0.71 28.12 -0.71 28.14 -0.71 28.16 -0.71 28.19 -0.71 28.21 -0.71 28.23 -0.71 28.26 -0.71 28.28 -0.71 28.3 -0.71 28.32 -0.71 28.35 -0.71 28.37 -0.71 28.39 -0.71 28.42 -0.71 28.44 -0.71 28.46 -0.71 28.49 -0.71 28.51 -0.71 28.53 -0.71 28.56 -0.71 28.58 -0.71 28.6 -0.71 28.63 -0.71 28.65 -0.71 28.67 -0.71 28.69 -0.71 28.72 -0.71 28.74 -0.71 28.76 -0.71 28.79 -0.71 28.81 -0.71 28.83 -0.71 28.86 -0.71 28.88 -0.71 28.9 -0.71 28.93 -0.71 28.95 -0.71 28.97 -0.71 28.99 -0.71 29.02 -0.71 29.04 -0.71 29.06 -0.71 29.09 -0.71 29.11 -0.71 29.13 -0.71 29.16 -0.71 29.18 -0.71 29.2 -0.71 29.23 -0.71 29.25 -0.71 29.27 -0.71 29.3 -0.71 29.32 -0.71 29.34 -0.71 29.36 -0.71 29.39 -0.71 29.41 -0.71 29.43 -0.71 29.46 -0.71 29.48 -0.71 29.5 -0.71 29.53 -0.71 29.55 -0.71 29.57 -0.71 29.6 -0.71 29.62 -0.71 29.64 -0.71 29.66 -0.71 29.69 -0.71 29.71 -0.71 29.73 -0.71 29.76 -0.71 29.78 -0.71 29.8 -0.71 29.83 -0.71 29.85 -0.71 29.87 -0.71 29.9 -0.71 29.92 -0.71 29.94 -0.71 29.97 -0.71 29.99 -0.71 30.01 -0.71 30.03 -0.71 30.06 -0.71 30.08 -0.71 30.1 -0.71 30.13 -0.71 30.15 -0.71 30.17 -0.71 30.2 -0.71 30.22 -0.71 30.24 -0.71 30.27 -0.71 30.29 -0.71 30.31 -0.71 30.33 -0.71 30.36 -0.71 30.38 -0.71 30.4 -0.71 30.43 -0.71 30.45 -0.71 30.47 -0.71 30.5 -0.71 30.52 -0.71 30.54 -0.71 30.57 -0.71 30.59 -0.71 30.61 -0.71 30.64 -0.71 30.66 -0.71 30.68 -0.71 30.7 -0.71 30.73 -0.71 30.75 -0.71 30.77 -0.71 30.8 -0.71 30.82 -0.71 30.84 -0.71 30.87 -0.71 30.89 -0.71 30.91 -0.71 30.94 -0.71 30.96 -0.71 30.98 -0.71 31 -0.71 31.03 -0.71 31.05 -0.71 31.07 -0.71 31.1 -0.71 31.12 -0.71 31.14 -0.71 31.17 -0.71 31.19 -0.71 31.21 -0.71 31.24 -0.71 31.26 -0.71 31.28 -0.71 31.31 -0.71 31.33 -0.71 31.35 -0.71 31.37 -0.71 31.4 -0.71 31.42 -0.71 31.44 -0.71 31.47 -0.71 31.49 -0.71 31.51 -0.71 31.54 -0.71 31.56 -0.71 31.58 -0.71 31.61 -0.71 31.63 -0.71 31.65 -0.71 31.67 -0.71 31.7 -0.71 31.72 -0.71 31.74 -0.71 31.77 -0.71 31.79 -0.71 31.81 -0.71 31.84 -0.71 31.86 -0.71 31.88 -0.71 31.91 -0.71 31.93 -0.71 31.95 -0.71 31.98 -0.71 32 -0.71 32.02 -0.71 32.04 -0.71 32.07 -0.71 32.09 -0.71 32.11 -0.71 32.14 -0.71 32.16 -0.71 32.18 -0.71 32.21 -0.71 32.23 -0.71 32.25 -0.71 32.28 -0.71 32.3 -0.71 32.32 -0.71 32.34 -0.71 32.37 -0.71 32.39 -0.71 32.41 -0.71 32.44 -0.71 32.46 -0.71 32.48 -0.71 32.51 -0.71 32.53 -0.71 32.55 -0.71 32.58 -0.71 32.6 -0.71 32.62 -0.71 32.65 -0.71 32.67 -0.71 32.69 -0.71 32.71 -0.71 32.74 -0.71 32.76 -0.71 32.78 -0.71 32.81 -0.71 32.83 -0.71 32.85 -0.71 32.88 -0.71 32.9 -0.71 32.92 -0.71 32.95 -0.71 32.97 -0.71 32.99 -0.71 33.01 -0.71 33.04 -0.71 33.06 -0.71 33.08 -0.71 33.11 -0.71 33.13 -0.71 33.15 -0.71 33.18 -0.71 33.2 -0.71 33.22 -0.71 33.25 -0.71 33.27 -0.71 33.29 -0.71 33.32 -0.71 33.34 -0.71 33.36 -0.71 33.38 -0.71 33.41 -0.71 33.43 -0.71 33.45 -0.71 33.48 -0.71 33.5 -0.71 33.52 -0.71 33.55 -0.71 33.57 -0.71 33.59 -0.71 33.62 -0.71 33.64 -0.71 33.66 -0.71 33.68 -0.71 33.71 -0.71 33.73 -0.71 33.75 -0.71 33.78 -0.71 33.8 -0.71 33.82 -0.71 33.85 -0.71 33.87 -0.71 33.89 -0.71 33.92 -0.71 33.94 -0.71 33.96 -0.71 33.99 -0.71 34.01 -0.71 34.03 -0.71 34.05 -0.71 34.08 -0.71 34.1 -0.71 34.12 -0.71 34.15 -0.71 34.17 -0.71 34.19 -0.71 34.22 -0.71 34.24 -0.71 34.26 -0.71 34.29 -0.71 34.31 -0.71 34.33 -0.71 34.35 -0.71 34.38 -0.71 34.4 -0.71 34.42 -0.71 34.45 -0.71 34.47 -0.71 34.49 -0.71 34.52 -0.71 34.54 -0.71 34.56 -0.71 34.59 -0.71 34.61 -0.71 34.63 -0.71 34.66 -0.71 34.68 -0.71 34.7 -0.71 34.72 -0.71 34.75 -0.71 34.77 -0.71 34.79 -0.71 34.82 -0.71 34.84 -0.71 34.86 -0.71 34.89 -0.71 34.91 -0.71 34.93 -0.71 34.96 -0.71 34.98 -0.71 35 -0.71 35.02 -0.71 35.05 -0.71 35.07 -0.71 35.09 -0.71 35.12 -0.71 35.14 -0.71 35.16 -0.71 35.19 -0.71 35.21 -0.71 35.23 -0.71 35.26 -0.71 35.28 -0.71 35.3 -0.71 35.33 -0.71 35.35 -0.71 35.37 -0.71 35.39 -0.71 35.42 -0.71 35.44 -0.71 35.46 -0.71 35.49 -0.71 35.51 -0.71 35.53 -0.71 35.56 -0.71 35.58 -0.71 35.6 -0.71 35.63 -0.71 35.65 -0.71 35.67 -0.71 35.69 -0.71 35.72 -0.71 35.74 -0.71 35.76 -0.71 35.79 -0.71 35.81 -0.71 35.83 -0.71 35.86 -0.71 35.88 -0.71 35.9 -0.71 35.93 -0.71 35.95 -0.71 35.97 -0.71 36 -0.71 36.02 -0.71 36.04 -0.71 36.06 -0.71 36.09 -0.71 36.11 -0.71 36.13 -0.71 36.16 -0.71 36.18 -0.71 36.2 -0.71 36.23 -0.71 36.25 -0.71 36.27 -0.71 36.3 -0.71 36.32 -0.71 36.34 -0.71 36.36 -0.71 36.39 -0.71 36.41 -0.71 36.43 -0.71 36.46 -0.71 36.48 -0.71 36.5 -0.71 36.53 -0.71 36.55 -0.71 36.57 -0.71 36.6 -0.71 36.62 -0.71 36.64 -0.71 36.67 -0.71 36.69 -0.71 36.71 -0.71 36.73 -0.71 36.76 -0.71 36.78 -0.71 36.8 -0.71 36.83 -0.71 36.85 -0.71 36.87 -0.71 36.9 -0.71 36.92 -0.71 36.94 -0.71 36.97 -0.71 36.99 -0.71 37.01 -0.71 37.03 -0.71 37.06 -0.71 37.08 -0.71 37.1 -0.71 37.13 -0.71 37.15 -0.71 37.17 -0.71 37.2 -0.71 37.22 -0.71 37.24 -0.71 37.27 -0.71 37.29 -0.71 37.31 -0.71 37.34 -0.71 37.36 -0.71 37.38 -0.71 37.4 -0.71 37.43 -0.71 37.45 -0.71 37.47 -0.71 37.5 -0.71 37.52 -0.71 37.54 -0.71 37.57 -0.71 37.59 -0.71 37.61 -0.71 37.64 -0.71 37.66 -0.71 37.68 -0.71 37.7 -0.71 37.73 -0.71 37.75 -0.71 37.77 -0.71 37.8 -0.71 37.82 -0.71 37.84 -0.71 37.87 -0.71 37.89 -0.71 37.91 -0.71 37.94 -0.71 37.96 -0.71 37.98 -0.71 38.01 -0.71 38.03 -0.71 38.05 -0.71 38.07 -0.71 38.1 -0.71 38.12 -0.71 38.14 -0.71 38.17 -0.71 38.19 -0.71 38.21 -0.71 38.24 -0.71 38.26 -0.71 38.28 -0.71 38.31 -0.71 38.33 -0.71 38.35 -0.71 38.37 -0.71 38.4 -0.71 38.42 -0.71 38.44 -0.71 38.47 -0.71 38.49 -0.71 38.51 -0.71 38.54 -0.71 38.56 -0.71 38.58 -0.71 38.61 -0.71 38.63 -0.71 38.65 -0.71 38.68 -0.71 38.7 -0.71 38.72 -0.71 38.74 -0.71 38.77 -0.71 38.79 -0.71 38.81 -0.71 38.84 -0.71 38.86 -0.71 38.88 -0.71 38.91 -0.71 38.93 -0.71 38.95 -0.71 38.98 -0.71 39 -0.71 39.02 -0.71 39.04 -0.71 39.07 -0.71 39.09 -0.71 39.11 -0.71 39.14 -0.71 39.16 -0.71 39.18 -0.71 39.21 -0.71 39.23 -0.71 39.25 -0.71 39.28 -0.71 39.3 -0.71 39.32 -0.71 39.35 -0.71 39.37 -0.71 39.39 -0.71 39.41 -0.71 39.44 -0.71 39.46 -0.71 39.48 -0.71 39.51 -0.71 39.53 -0.71 39.55 -0.71 39.58 -0.71 39.6 -0.71 39.62 -0.71 39.65 -0.71 39.67 -0.71 39.69 -0.71 39.71 -0.71 39.74 -0.71 39.76 -0.71 39.78 -0.71 39.81 -0.71 39.83 -0.71 39.85 -0.71 39.88 -0.71 39.9 -0.71 39.92 -0.71 39.95 -0.71 39.97 -0.71 39.99 -0.71 40.02 -0.71 40.04 -0.71 40.06 -0.71 40.08 -0.71 40.11 -0.71 40.13 -0.71 40.15 -0.71 40.18 -0.71 40.2 -0.71 40.22 -0.71 40.25 -0.71 40.27 -0.71 40.29 -0.71 40.32 -0.71 40.34 -0.71 40.36 -0.71 40.38 -0.71 40.41 -0.71 40.43 -0.71 40.45 -0.71 40.48 -0.71 40.5 -0.71 40.52 -0.71 40.55 -0.71 40.57 -0.71 40.59 -0.71 40.62 -0.71 40.64 -0.71 40.66 -0.71 40.69 -0.71 40.71 -0.71 40.73 -0.71 40.75 -0.71 40.78 -0.71 40.8 -0.71 40.82 -0.71 40.85 -0.71 40.87 -0.71 40.89 -0.71 40.92 -0.71 40.94 -0.71 40.96 -0.71 40.99 -0.71 41.01 -0.71 41.03 -0.71 41.05 -0.71 41.08 -0.71 41.1 -0.71 41.12 -0.71 41.15 -0.71 41.17 -0.71 41.19 -0.71 41.22 -0.71 41.24 -0.71 41.26 -0.71 41.29 -0.71 41.31 -0.71 41.33 -0.71 41.36 -0.71 41.38 -0.71 41.4 -0.71 41.42 -0.71 41.45 -0.71 41.47 -0.71 41.49 -0.71 41.52 -0.71 41.54 -0.71 41.56 -0.71 41.59 -0.71 41.61 -0.71 41.63 -0.71 41.66 -0.71 41.68 -0.71 41.7 -0.71 41.72 -0.71 41.75 -0.71 41.77 -0.71 41.79 -0.71 41.82 -0.71 41.84 -0.71 41.86 -0.71 41.89 -0.71 41.91 -0.71 41.93 -0.71 41.96 -0.71 41.98 -0.71 42 -0.71 42.03 -0.71 42.05 -0.71 42.07 -0.71 42.09 -0.71 42.12 -0.71 42.14 -0.71 42.16 -0.71" class="primitive"/>
          </g>
          <g transform="translate(70.89,78.36)" id="img-8311e7c5-234" class="geometry color_E_V" stroke="#D4CA3A">
            <path fill="none" d="M-42.16,0.36 L -42.14 0.36 -42.12 0.36 -42.09 0.36 -42.07 0.36 -42.05 0.35 -42.03 0.35 -42 0.35 -41.98 0.35 -41.96 0.35 -41.93 0.35 -41.91 0.35 -41.89 0.35 -41.86 0.35 -41.84 0.35 -41.82 0.35 -41.79 0.35 -41.77 0.34 -41.75 0.34 -41.72 0.34 -41.7 0.34 -41.68 0.34 -41.66 0.34 -41.63 0.34 -41.61 0.34 -41.59 0.34 -41.56 0.34 -41.54 0.34 -41.52 0.34 -41.49 0.34 -41.47 0.34 -41.45 0.34 -41.42 0.34 -41.4 0.34 -41.38 0.34 -41.36 0.34 -41.33 0.33 -41.31 0.33 -41.29 0.33 -41.26 0.33 -41.24 0.33 -41.22 0.33 -41.19 0.33 -41.17 0.33 -41.15 0.33 -41.12 0.33 -41.1 0.33 -41.08 0.33 -41.05 0.33 -41.03 0.33 -41.01 0.33 -40.99 0.33 -40.96 0.33 -40.94 0.33 -40.92 0.33 -40.89 0.33 -40.87 0.33 -40.85 0.33 -40.82 0.33 -40.8 0.33 -40.78 0.33 -40.75 0.32 -40.73 0.32 -40.71 0.32 -40.69 0.32 -40.66 0.32 -40.64 0.32 -40.62 0.32 -40.59 0.32 -40.57 0.32 -40.55 0.32 -40.52 0.32 -40.5 0.32 -40.48 0.32 -40.45 0.32 -40.43 0.32 -40.41 0.32 -40.38 0.32 -40.36 0.32 -40.34 0.32 -40.32 0.32 -40.29 0.32 -40.27 0.31 -40.25 0.31 -40.22 0.31 -40.2 0.31 -40.18 0.31 -40.15 0.31 -40.13 0.31 -40.11 0.31 -40.08 0.31 -40.06 0.31 -40.04 0.31 -40.02 0.31 -39.99 0.31 -39.97 0.31 -39.95 0.3 -39.92 0.3 -39.9 0.3 -39.88 0.3 -39.85 0.3 -39.83 0.3 -39.81 0.3 -39.78 0.3 -39.76 0.3 -39.74 0.3 -39.71 0.29 -39.69 0.29 -39.67 0.29 -39.65 0.29 -39.62 0.29 -39.6 0.29 -39.58 0.29 -39.55 0.29 -39.53 0.28 -39.51 0.28 -39.48 0.28 -39.46 0.28 -39.44 0.28 -39.41 0.28 -39.39 0.28 -39.37 0.27 -39.35 0.27 -39.32 0.27 -39.3 0.27 -39.28 0.27 -39.25 0.27 -39.23 0.26 -39.21 0.26 -39.18 0.26 -39.16 0.26 -39.14 0.26 -39.11 0.26 -39.09 0.25 -39.07 0.25 -39.04 0.25 -39.02 0.25 -39 0.25 -38.98 0.24 -38.95 0.24 -38.93 0.24 -38.91 0.24 -38.88 0.23 -38.86 0.23 -38.84 0.23 -38.81 0.23 -38.79 0.22 -38.77 0.22 -38.74 0.22 -38.72 0.22 -38.7 0.21 -38.68 0.21 -38.65 0.21 -38.63 0.21 -38.61 0.2 -38.58 0.2 -38.56 0.2 -38.54 0.19 -38.51 0.19 -38.49 0.19 -38.47 0.18 -38.44 0.18 -38.42 0.18 -38.4 0.17 -38.37 0.17 -38.35 0.17 -38.33 0.16 -38.31 0.16 -38.28 0.16 -38.26 0.15 -38.24 0.15 -38.21 0.14 -38.19 0.14 -38.17 0.14 -38.14 0.13 -38.12 0.13 -38.1 0.12 -38.07 0.12 -38.05 0.11 -38.03 0.11 -38.01 0.11 -37.98 0.1 -37.96 0.1 -37.94 0.09 -37.91 0.09 -37.89 0.08 -37.87 0.08 -37.84 0.07 -37.82 0.07 -37.8 0.06 -37.77 0.06 -37.75 0.05 -37.73 0.05 -37.7 0.04 -37.68 0.03 -37.66 0.03 -37.64 0.02 -37.61 0.02 -37.59 0.01 -37.57 0 -37.54 -0 -37.52 -0.01 -37.5 -0.01 -37.47 -0.02 -37.45 -0.03 -37.43 -0.03 -37.4 -0.04 -37.38 -0.05 -37.36 -0.05 -37.34 -0.06 -37.31 -0.07 -37.29 -0.08 -37.27 -0.08 -37.24 -0.09 -37.22 -0.1 -37.2 -0.11 -37.17 -0.11 -37.15 -0.12 -37.13 -0.13 -37.1 -0.14 -37.08 -0.15 -37.06 -0.15 -37.03 -0.16 -37.01 -0.17 -36.99 -0.18 -36.97 -0.19 -36.94 -0.2 -36.92 -0.21 -36.9 -0.21 -36.87 -0.22 -36.85 -0.23 -36.83 -0.24 -36.8 -0.25 -36.78 -0.26 -36.76 -0.27 -36.73 -0.28 -36.71 -0.29 -36.69 -0.3 -36.67 -0.31 -36.64 -0.32 -36.62 -0.33 -36.6 -0.34 -36.57 -0.35 -36.55 -0.36 -36.53 -0.37 -36.5 -0.38 -36.48 -0.4 -36.46 -0.41 -36.43 -0.42 -36.41 -0.43 -36.39 -0.44 -36.36 -0.45 -36.34 -0.46 -36.32 -0.48 -36.3 -0.49 -36.27 -0.5 -36.25 -0.51 -36.23 -0.52 -36.2 -0.54 -36.18 -0.55 -36.16 -0.56 -36.13 -0.57 -36.11 -0.58 -36.09 -0.6 -36.06 -0.61 -36.04 -0.62 -36.02 -0.64 -36 -0.65 -35.97 -0.66 -35.95 -0.67 -35.93 -0.69 -35.9 -0.7 -35.88 -0.71 -35.86 -0.73 -35.83 -0.74 -35.81 -0.75 -35.79 -0.77 -35.76 -0.78 -35.74 -0.79 -35.72 -0.81 -35.69 -0.82 -35.67 -0.84 -35.65 -0.85 -35.63 -0.86 -35.6 -0.88 -35.58 -0.89 -35.56 -0.9 -35.53 -0.92 -35.51 -0.93 -35.49 -0.95 -35.46 -0.96 -35.44 -0.97 -35.42 -0.99 -35.39 -1 -35.37 -1.02 -35.35 -1.03 -35.33 -1.04 -35.3 -1.06 -35.28 -1.07 -35.26 -1.08 -35.23 -1.1 -35.21 -1.11 -35.19 -1.12 -35.16 -1.14 -35.14 -1.15 -35.12 -1.17 -35.09 -1.18 -35.07 -1.19 -35.05 -1.21 -35.02 -1.22 -35 -1.23 -34.98 -1.24 -34.96 -1.26 -34.93 -1.27 -34.91 -1.28 -34.89 -1.3 -34.86 -1.31 -34.84 -1.32 -34.82 -1.33 -34.79 -1.34 -34.77 -1.36 -34.75 -1.37 -34.72 -1.38 -34.7 -1.39 -34.68 -1.4 -34.66 -1.41 -34.63 -1.43 -34.61 -1.44 -34.59 -1.45 -34.56 -1.46 -34.54 -1.47 -34.52 -1.48 -34.49 -1.49 -34.47 -1.5 -34.45 -1.51 -34.42 -1.52 -34.4 -1.53 -34.38 -1.54 -34.35 -1.55 -34.33 -1.56 -34.31 -1.57 -34.29 -1.57 -34.26 -1.58 -34.24 -1.59 -34.22 -1.6 -34.19 -1.61 -34.17 -1.61 -34.15 -1.62 -34.12 -1.63 -34.1 -1.64 -34.08 -1.64 -34.05 -1.65 -34.03 -1.66 -34.01 -1.66 -33.99 -1.67 -33.96 -1.67 -33.94 -1.68 -33.92 -1.69 -33.89 -1.69 -33.87 -1.7 -33.85 -1.7 -33.82 -1.7 -33.8 -1.71 -33.78 -1.71 -33.75 -1.72 -33.73 -1.72 -33.71 -1.72 -33.68 -1.73 -33.66 -1.73 -33.64 -1.73 -33.62 -1.74 -33.59 -1.74 -33.57 -1.74 -33.55 -1.74 -33.52 -1.75 -33.5 -1.75 -33.48 -1.75 -33.45 -1.75 -33.43 -1.75 -33.41 -1.75 -33.38 -1.75 -33.36 -1.75 -33.34 -1.75 -33.32 -1.75 -33.29 -1.75 -33.27 -1.75 -33.25 -1.75 -33.22 -1.75 -33.2 -1.75 -33.18 -1.75 -33.15 -1.75 -33.13 -1.75 -33.11 -1.75 -33.08 -1.75 -33.06 -1.74 -33.04 -1.74 -33.01 -1.74 -32.99 -1.74 -32.97 -1.74 -32.95 -1.73 -32.92 -1.73 -32.9 -1.73 -32.88 -1.73 -32.85 -1.72 -32.83 -1.72 -32.81 -1.72 -32.78 -1.71 -32.76 -1.71 -32.74 -1.71 -32.71 -1.7 -32.69 -1.7 -32.67 -1.69 -32.65 -1.69 -32.62 -1.69 -32.6 -1.68 -32.58 -1.68 -32.55 -1.67 -32.53 -1.67 -32.51 -1.66 -32.48 -1.66 -32.46 -1.65 -32.44 -1.65 -32.41 -1.64 -32.39 -1.64 -32.37 -1.63 -32.34 -1.63 -32.32 -1.62 -32.3 -1.61 -32.28 -1.61 -32.25 -1.6 -32.23 -1.6 -32.21 -1.59 -32.18 -1.58 -32.16 -1.58 -32.14 -1.57 -32.11 -1.57 -32.09 -1.56 -32.07 -1.55 -32.04 -1.55 -32.02 -1.54 -32 -1.53 -31.98 -1.53 -31.95 -1.52 -31.93 -1.51 -31.91 -1.51 -31.88 -1.5 -31.86 -1.49 -31.84 -1.49 -31.81 -1.48 -31.79 -1.47 -31.77 -1.47 -31.74 -1.46 -31.72 -1.45 -31.7 -1.44 -31.67 -1.44 -31.65 -1.43 -31.63 -1.42 -31.61 -1.42 -31.58 -1.41 -31.56 -1.4 -31.54 -1.39 -31.51 -1.39 -31.49 -1.38 -31.47 -1.37 -31.44 -1.36 -31.42 -1.36 -31.4 -1.35 -31.37 -1.34 -31.35 -1.33 -31.33 -1.33 -31.31 -1.32 -31.28 -1.31 -31.26 -1.31 -31.24 -1.3 -31.21 -1.29 -31.19 -1.28 -31.17 -1.28 -31.14 -1.27 -31.12 -1.26 -31.1 -1.25 -31.07 -1.25 -31.05 -1.24 -31.03 -1.23 -31 -1.22 -30.98 -1.22 -30.96 -1.21 -30.94 -1.2 -30.91 -1.19 -30.89 -1.19 -30.87 -1.18 -30.84 -1.17 -30.82 -1.16 -30.8 -1.16 -30.77 -1.15 -30.75 -1.14 -30.73 -1.13 -30.7 -1.13 -30.68 -1.12 -30.66 -1.11 -30.64 -1.11 -30.61 -1.1 -30.59 -1.09 -30.57 -1.08 -30.54 -1.08 -30.52 -1.07 -30.5 -1.06 -30.47 -1.06 -30.45 -1.05 -30.43 -1.04 -30.4 -1.03 -30.38 -1.03 -30.36 -1.02 -30.33 -1.01 -30.31 -1.01 -30.29 -1 -30.27 -0.99 -30.24 -0.98 -30.22 -0.98 -30.2 -0.97 -30.17 -0.96 -30.15 -0.96 -30.13 -0.95 -30.1 -0.94 -30.08 -0.94 -30.06 -0.93 -30.03 -0.92 -30.01 -0.92 -29.99 -0.91 -29.97 -0.9 -29.94 -0.9 -29.92 -0.89 -29.9 -0.88 -29.87 -0.88 -29.85 -0.87 -29.83 -0.86 -29.8 -0.86 -29.78 -0.85 -29.76 -0.85 -29.73 -0.84 -29.71 -0.83 -29.69 -0.83 -29.66 -0.82 -29.64 -0.81 -29.62 -0.81 -29.6 -0.8 -29.57 -0.8 -29.55 -0.79 -29.53 -0.78 -29.5 -0.78 -29.48 -0.77 -29.46 -0.77 -29.43 -0.76 -29.41 -0.75 -29.39 -0.75 -29.36 -0.74 -29.34 -0.74 -29.32 -0.73 -29.3 -0.72 -29.27 -0.72 -29.25 -0.71 -29.23 -0.71 -29.2 -0.7 -29.18 -0.7 -29.16 -0.69 -29.13 -0.68 -29.11 -0.68 -29.09 -0.67 -29.06 -0.67 -29.04 -0.66 -29.02 -0.66 -28.99 -0.65 -28.97 -0.65 -28.95 -0.64 -28.93 -0.64 -28.9 -0.63 -28.88 -0.63 -28.86 -0.62 -28.83 -0.61 -28.81 -0.61 -28.79 -0.6 -28.76 -0.6 -28.74 -0.59 -28.72 -0.59 -28.69 -0.58 -28.67 -0.58 -28.65 -0.57 -28.63 -0.57 -28.6 -0.56 -28.58 -0.56 -28.56 -0.55 -28.53 -0.55 -28.51 -0.55 -28.49 -0.54 -28.46 -0.54 -28.44 -0.53 -28.42 -0.53 -28.39 -0.52 -28.37 -0.52 -28.35 -0.51 -28.32 -0.51 -28.3 -0.5 -28.28 -0.5 -28.26 -0.49 -28.23 -0.49 -28.21 -0.49 -28.19 -0.48 -28.16 -0.48 -28.14 -0.47 -28.12 -0.47 -28.09 -0.46 -28.07 -0.46 -28.05 -0.46 -28.02 -0.45 -28 -0.45 -27.98 -0.44 -27.96 -0.44 -27.93 -0.43 -27.91 -0.43 -27.89 -0.43 -27.86 -0.42 -27.84 -0.42 -27.82 -0.41 -27.79 -0.41 -27.77 -0.41 -27.75 -0.4 -27.72 -0.4 -27.7 -0.39 -27.68 -0.39 -27.65 -0.39 -27.63 -0.38 -27.61 -0.38 -27.59 -0.38 -27.56 -0.37 -27.54 -0.37 -27.52 -0.36 -27.49 -0.36 -27.47 -0.36 -27.45 -0.35 -27.42 -0.35 -27.4 -0.35 -27.38 -0.34 -27.35 -0.34 -27.33 -0.34 -27.31 -0.33 -27.29 -0.33 -27.26 -0.33 -27.24 -0.32 -27.22 -0.32 -27.19 -0.32 -27.17 -0.31 -27.15 -0.31 -27.12 -0.31 -27.1 -0.3 -27.08 -0.3 -27.05 -0.3 -27.03 -0.29 -27.01 -0.29 -26.98 -0.29 -26.96 -0.28 -26.94 -0.28 -26.92 -0.28 -26.89 -0.27 -26.87 -0.27 -26.85 -0.27 -26.82 -0.27 -26.8 -0.26 -26.78 -0.26 -26.75 -0.26 -26.73 -0.25 -26.71 -0.25 -26.68 -0.25 -26.66 -0.24 -26.64 -0.24 -26.62 -0.24 -26.59 -0.24 -26.57 -0.23 -26.55 -0.23 -26.52 -0.23 -26.5 -0.23 -26.48 -0.22 -26.45 -0.22 -26.43 -0.22 -26.41 -0.21 -26.38 -0.21 -26.36 -0.21 -26.34 -0.21 -26.31 -0.2 -26.29 -0.2 -26.27 -0.2 -26.25 -0.2 -26.22 -0.19 -26.2 -0.19 -26.18 -0.19 -26.15 -0.19 -26.13 -0.18 -26.11 -0.18 -26.08 -0.18 -26.06 -0.18 -26.04 -0.17 -26.01 -0.17 -25.99 -0.17 -25.97 -0.17 -25.95 -0.16 -25.92 -0.16 -25.9 -0.16 -25.88 -0.16 -25.85 -0.16 -25.83 -0.15 -25.81 -0.15 -25.78 -0.15 -25.76 -0.15 -25.74 -0.14 -25.71 -0.14 -25.69 -0.14 -25.67 -0.14 -25.64 -0.14 -25.62 -0.13 -25.6 -0.13 -25.58 -0.13 -25.55 -0.13 -25.53 -0.13 -25.51 -0.12 -25.48 -0.12 -25.46 -0.12 -25.44 -0.12 -25.41 -0.12 -25.39 -0.11 -25.37 -0.11 -25.34 -0.11 -25.32 -0.11 -25.3 -0.11 -25.28 -0.1 -25.25 -0.1 -25.23 -0.1 -25.21 -0.1 -25.18 -0.1 -25.16 -0.09 -25.14 -0.09 -25.11 -0.09 -25.09 -0.09 -25.07 -0.09 -25.04 -0.09 -25.02 -0.08 -25 -0.08 -24.97 -0.08 -24.95 -0.08 -24.93 -0.08 -24.91 -0.08 -24.88 -0.07 -24.86 -0.07 -24.84 -0.07 -24.81 -0.07 -24.79 -0.07 -24.77 -0.07 -24.74 -0.06 -24.72 -0.06 -24.7 -0.06 -24.67 -0.06 -24.65 -0.06 -24.63 -0.06 -24.61 -0.05 -24.58 -0.05 -24.56 -0.05 -24.54 -0.05 -24.51 -0.05 -24.49 -0.05 -24.47 -0.04 -24.44 -0.04 -24.42 -0.04 -24.4 -0.04 -24.37 -0.04 -24.35 -0.04 -24.33 -0.04 -24.3 -0.03 -24.28 -0.03 -24.26 -0.03 -24.24 -0.03 -24.21 -0.03 -24.19 -0.03 -24.17 -0.03 -24.14 -0.02 -24.12 -0.02 -24.1 -0.02 -24.07 -0.02 -24.05 -0.02 -24.03 -0.02 -24 -0.02 -23.98 -0.02 -23.96 -0.01 -23.94 -0.01 -23.91 -0.01 -23.89 -0.01 -23.87 -0.01 -23.84 -0.01 -23.82 -0.01 -23.8 -0.01 -23.77 -0 -23.75 -0 -23.73 -0 -23.7 -0 -23.68 0 -23.66 0 -23.63 0 -23.61 0 -23.59 0 -23.57 0.01 -23.54 0.01 -23.52 0.01 -23.5 0.01 -23.47 0.01 -23.45 0.01 -23.43 0.01 -23.4 0.01 -23.38 0.01 -23.36 0.02 -23.33 0.02 -23.31 0.02 -23.29 0.02 -23.27 0.02 -23.24 0.02 -23.22 0.02 -23.2 0.02 -23.17 0.02 -23.15 0.03 -23.13 0.03 -23.1 0.03 -23.08 0.03 -23.06 0.03 -23.03 0.03 -23.01 0.03 -22.99 0.03 -22.96 0.03 -22.94 0.03 -22.92 0.03 -22.9 0.04 -22.87 0.04 -22.85 0.04 -22.83 0.04 -22.8 0.04 -22.78 0.04 -22.76 0.04 -22.73 0.04 -22.71 0.04 -22.69 0.04 -22.66 0.04 -22.64 0.05 -22.62 0.05 -22.6 0.05 -22.57 0.05 -22.55 0.05 -22.53 0.05 -22.5 0.05 -22.48 0.05 -22.46 0.05 -22.43 0.05 -22.41 0.05 -22.39 0.05 -22.36 0.06 -22.34 0.06 -22.32 0.06 -22.29 0.06 -22.27 0.06 -22.25 0.06 -22.23 0.06 -22.2 0.06 -22.18 0.06 -22.16 0.06 -22.13 0.06 -22.11 0.06 -22.09 0.06 -22.06 0.06 -22.04 0.07 -22.02 0.07 -21.99 0.07 -21.97 0.07 -21.95 0.07 -21.93 0.07 -21.9 0.07 -21.88 0.07 -21.86 0.07 -21.83 0.07 -21.81 0.07 -21.79 0.07 -21.76 0.07 -21.74 0.07 -21.72 0.07 -21.69 0.08 -21.67 0.08 -21.65 0.08 -21.62 0.08 -21.6 0.08 -21.58 0.08 -21.56 0.08 -21.53 0.08 -21.51 0.08 -21.49 0.08 -21.46 0.08 -21.44 0.08 -21.42 0.08 -21.39 0.08 -21.37 0.08 -21.35 0.08 -21.32 0.09 -21.3 0.09 -21.28 0.09 -21.26 0.09 -21.23 0.09 -21.21 0.09 -21.19 0.09 -21.16 0.09 -21.14 0.09 -21.12 0.09 -21.09 0.09 -21.07 0.09 -21.05 0.09 -21.02 0.09 -21 0.09 -20.98 0.09 -20.95 0.09 -20.93 0.09 -20.91 0.09 -20.89 0.1 -20.86 0.1 -20.84 0.1 -20.82 0.1 -20.79 0.1 -20.77 0.1 -20.75 0.1 -20.72 0.1 -20.7 0.1 -20.68 0.1 -20.65 0.1 -20.63 0.1 -20.61 0.1 -20.59 0.1 -20.56 0.1 -20.54 0.1 -20.52 0.1 -20.49 0.1 -20.47 0.1 -20.45 0.1 -20.42 0.1 -20.4 0.1 -20.38 0.11 -20.35 0.11 -20.33 0.11 -20.31 0.11 -20.28 0.11 -20.26 0.11 -20.24 0.11 -20.22 0.11 -20.19 0.11 -20.17 0.11 -20.15 0.11 -20.12 0.11 -20.1 0.11 -20.08 0.11 -20.05 0.11 -20.03 0.11 -20.01 0.11 -19.98 0.11 -19.96 0.11 -19.94 0.11 -19.92 0.11 -19.89 0.11 -19.87 0.11 -19.85 0.11 -19.82 0.11 -19.8 0.11 -19.78 0.11 -19.75 0.12 -19.73 0.12 -19.71 0.12 -19.68 0.12 -19.66 0.12 -19.64 0.12 -19.61 0.12 -19.59 0.12 -19.57 0.12 -19.55 0.12 -19.52 0.12 -19.5 0.12 -19.48 0.12 -19.45 0.12 -19.43 0.12 -19.41 0.12 -19.38 0.12 -19.36 0.12 -19.34 0.12 -19.31 0.12 -19.29 0.12 -19.27 0.12 -19.25 0.12 -19.22 0.12 -19.2 0.12 -19.18 0.12 -19.15 0.12 -19.13 0.12 -19.11 0.12 -19.08 0.12 -19.06 0.12 -19.04 0.12 -19.01 0.12 -18.99 0.13 -18.97 0.13 -18.94 0.13 -18.92 0.13 -18.9 0.13 -18.88 0.13 -18.85 0.13 -18.83 0.13 -18.81 0.13 -18.78 0.13 -18.76 0.13 -18.74 0.13 -18.71 0.13 -18.69 0.13 -18.67 0.13 -18.64 0.13 -18.62 0.13 -18.6 0.13 -18.58 0.13 -18.55 0.13 -18.53 0.13 -18.51 0.13 -18.48 0.13 -18.46 0.13 -18.44 0.13 -18.41 0.13 -18.39 0.13 -18.37 0.13 -18.34 0.13 -18.32 0.13 -18.3 0.13 -18.27 0.13 -18.25 0.13 -18.23 0.13 -18.21 0.13 -18.18 0.13 -18.16 0.13 -18.14 0.13 -18.11 0.13 -18.09 0.13 -18.07 0.13 -18.04 0.13 -18.02 0.13 -18 0.13 -17.97 0.14 -17.95 0.14 -17.93 0.14 -17.91 0.14 -17.88 0.14 -17.86 0.14 -17.84 0.14 -17.81 0.14 -17.79 0.14 -17.77 0.14 -17.74 0.14 -17.72 0.14 -17.7 0.14 -17.67 0.14 -17.65 0.14 -17.63 0.14 -17.6 0.14 -17.58 0.14 -17.56 0.14 -17.54 0.14 -17.51 0.14 -17.49 0.14 -17.47 0.14 -17.44 0.14 -17.42 0.14 -17.4 0.14 -17.37 0.14 -17.35 0.14 -17.33 0.14 -17.3 0.14 -17.28 0.14 -17.26 0.14 -17.24 0.14 -17.21 0.14 -17.19 0.14 -17.17 0.14 -17.14 0.14 -17.12 0.14 -17.1 0.14 -17.07 0.14 -17.05 0.14 -17.03 0.14 -17 0.14 -16.98 0.14 -16.96 0.14 -16.93 0.14 -16.91 0.14 -16.89 0.14 -16.87 0.14 -16.84 0.14 -16.82 0.14 -16.8 0.14 -16.77 0.14 -16.75 0.14 -16.73 0.14 -16.7 0.14 -16.68 0.14 -16.66 0.14 -16.63 0.14 -16.61 0.14 -16.59 0.14 -16.57 0.14 -16.54 0.14 -16.52 0.14 -16.5 0.14 -16.47 0.14 -16.45 0.14 -16.43 0.15 -16.4 0.15 -16.38 0.15 -16.36 0.15 -16.33 0.15 -16.31 0.15 -16.29 0.15 -16.26 0.15 -16.24 0.15 -16.22 0.15 -16.2 0.15 -16.17 0.15 -16.15 0.15 -16.13 0.15 -16.1 0.15 -16.08 0.15 -16.06 0.15 -16.03 0.15 -16.01 0.15 -15.99 0.15 -15.96 0.15 -15.94 0.15 -15.92 0.15 -15.9 0.15 -15.87 0.15 -15.85 0.15 -15.83 0.15 -15.8 0.15 -15.78 0.15 -15.76 0.15 -15.73 0.15 -15.71 0.15 -15.69 0.15 -15.66 0.15 -15.64 0.15 -15.62 0.15 -15.59 0.15 -15.57 0.15 -15.55 0.15 -15.53 0.15 -15.5 0.15 -15.48 0.15 -15.46 0.15 -15.43 0.15 -15.41 0.15 -15.39 0.15 -15.36 0.15 -15.34 0.15 -15.32 0.15 -15.29 0.15 -15.27 0.15 -15.25 0.15 -15.23 0.15 -15.2 0.15 -15.18 0.15 -15.16 0.15 -15.13 0.15 -15.11 0.15 -15.09 0.15 -15.06 0.15 -15.04 0.15 -15.02 0.15 -14.99 0.15 -14.97 0.15 -14.95 0.15 -14.92 0.15 -14.9 0.15 -14.88 0.15 -14.86 0.15 -14.83 0.15 -14.81 0.15 -14.79 0.15 -14.76 0.15 -14.74 0.15 -14.72 0.15 -14.69 0.15 -14.67 0.15 -14.65 0.15 -14.62 0.15 -14.6 0.15 -14.58 0.15 -14.56 0.15 -14.53 0.15 -14.51 0.15 -14.49 0.15 -14.46 0.15 -14.44 0.15 -14.42 0.15 -14.39 0.15 -14.37 0.15 -14.35 0.15 -14.32 0.15 -14.3 0.15 -14.28 0.15 -14.25 0.15 -14.23 0.15 -14.21 0.15 -14.19 0.15 -14.16 0.15 -14.14 0.15 -14.12 0.15 -14.09 0.15 -14.07 0.15 -14.05 0.15 -14.02 0.15 -14 0.15 -13.98 0.15 -13.95 0.15 -13.93 0.15 -13.91 0.15 -13.89 0.15 -13.86 0.15 -13.84 0.15 -13.82 0.15 -13.79 0.15 -13.77 0.15 -13.75 0.15 -13.72 0.15 -13.7 0.15 -13.68 0.15 -13.65 0.15 -13.63 0.15 -13.61 0.15 -13.58 0.15 -13.56 0.15 -13.54 0.15 -13.52 0.15 -13.49 0.15 -13.47 0.15 -13.45 0.15 -13.42 0.15 -13.4 0.15 -13.38 0.15 -13.35 0.15 -13.33 0.15 -13.31 0.15 -13.28 0.15 -13.26 0.15 -13.24 0.15 -13.22 0.15 -13.19 0.15 -13.17 0.15 -13.15 0.15 -13.12 0.15 -13.1 0.15 -13.08 0.15 -13.05 0.15 -13.03 0.15 -13.01 0.15 -12.98 0.15 -12.96 0.15 -12.94 0.15 -12.91 0.15 -12.89 0.15 -12.87 0.15 -12.85 0.15 -12.82 0.15 -12.8 0.15 -12.78 0.15 -12.75 0.15 -12.73 0.15 -12.71 0.16 -12.68 0.16 -12.66 0.16 -12.64 0.16 -12.61 0.16 -12.59 0.16 -12.57 0.16 -12.55 0.16 -12.52 0.16 -12.5 0.16 -12.48 0.16 -12.45 0.16 -12.43 0.16 -12.41 0.16 -12.38 0.16 -12.36 0.16 -12.34 0.16 -12.31 0.16 -12.29 0.16 -12.27 0.16 -12.24 0.16 -12.22 0.16 -12.2 0.16 -12.18 0.16 -12.15 0.16 -12.13 0.16 -12.11 0.16 -12.08 0.16 -12.06 0.16 -12.04 0.16 -12.01 0.16 -11.99 0.16 -11.97 0.16 -11.94 0.16 -11.92 0.16 -11.9 0.16 -11.88 0.16 -11.85 0.16 -11.83 0.16 -11.81 0.16 -11.78 0.16 -11.76 0.16 -11.74 0.16 -11.71 0.16 -11.69 0.16 -11.67 0.16 -11.64 0.16 -11.62 0.16 -11.6 0.16 -11.57 0.16 -11.55 0.16 -11.53 0.16 -11.51 0.16 -11.48 0.16 -11.46 0.16 -11.44 0.16 -11.41 0.16 -11.39 0.16 -11.37 0.16 -11.34 0.16 -11.32 0.16 -11.3 0.16 -11.27 0.16 -11.25 0.16 -11.23 0.16 -11.21 0.16 -11.18 0.16 -11.16 0.16 -11.14 0.16 -11.11 0.16 -11.09 0.16 -11.07 0.16 -11.04 0.16 -11.02 0.16 -11 0.16 -10.97 0.16 -10.95 0.16 -10.93 0.16 -10.9 0.16 -10.88 0.16 -10.86 0.16 -10.84 0.16 -10.81 0.16 -10.79 0.16 -10.77 0.16 -10.74 0.16 -10.72 0.16 -10.7 0.16 -10.67 0.16 -10.65 0.16 -10.63 0.16 -10.6 0.16 -10.58 0.16 -10.56 0.16 -10.54 0.16 -10.51 0.16 -10.49 0.16 -10.47 0.16 -10.44 0.16 -10.42 0.16 -10.4 0.16 -10.37 0.16 -10.35 0.16 -10.33 0.16 -10.3 0.16 -10.28 0.16 -10.26 0.16 -10.23 0.16 -10.21 0.16 -10.19 0.16 -10.17 0.16 -10.14 0.16 -10.12 0.16 -10.1 0.16 -10.07 0.16 -10.05 0.16 -10.03 0.16 -10 0.16 -9.98 0.16 -9.96 0.16 -9.93 0.16 -9.91 0.16 -9.89 0.16 -9.87 0.16 -9.84 0.16 -9.82 0.16 -9.8 0.16 -9.77 0.16 -9.75 0.16 -9.73 0.16 -9.7 0.16 -9.68 0.16 -9.66 0.16 -9.63 0.16 -9.61 0.16 -9.59 0.16 -9.56 0.16 -9.54 0.16 -9.52 0.16 -9.5 0.16 -9.47 0.16 -9.45 0.16 -9.43 0.16 -9.4 0.16 -9.38 0.16 -9.36 0.16 -9.33 0.16 -9.31 0.16 -9.29 0.16 -9.26 0.16 -9.24 0.16 -9.22 0.16 -9.2 0.16 -9.17 0.16 -9.15 0.16 -9.13 0.16 -9.1 0.16 -9.08 0.16 -9.06 0.16 -9.03 0.16 -9.01 0.16 -8.99 0.16 -8.96 0.16 -8.94 0.16 -8.92 0.16 -8.89 0.16 -8.87 0.16 -8.85 0.16 -8.83 0.16 -8.8 0.16 -8.78 0.16 -8.76 0.16 -8.73 0.16 -8.71 0.16 -8.69 0.16 -8.66 0.16 -8.64 0.16 -8.62 0.16 -8.59 0.16 -8.57 0.16 -8.55 0.16 -8.53 0.16 -8.5 0.16 -8.48 0.16 -8.46 0.16 -8.43 0.16 -8.41 0.16 -8.39 0.16 -8.36 0.16 -8.34 0.16 -8.32 0.16 -8.29 0.16 -8.27 0.16 -8.25 0.16 -8.22 0.16 -8.2 0.16 -8.18 0.16 -8.16 0.16 -8.13 0.16 -8.11 0.16 -8.09 0.16 -8.06 0.16 -8.04 0.16 -8.02 0.16 -7.99 0.16 -7.97 0.16 -7.95 0.16 -7.92 0.16 -7.9 0.16 -7.88 0.16 -7.86 0.16 -7.83 0.16 -7.81 0.16 -7.79 0.16 -7.76 0.16 -7.74 0.16 -7.72 0.16 -7.69 0.16 -7.67 0.16 -7.65 0.16 -7.62 0.16 -7.6 0.16 -7.58 0.16 -7.55 0.16 -7.53 0.16 -7.51 0.16 -7.49 0.16 -7.46 0.16 -7.44 0.16 -7.42 0.16 -7.39 0.16 -7.37 0.16 -7.35 0.16 -7.32 0.16 -7.3 0.16 -7.28 0.16 -7.25 0.16 -7.23 0.16 -7.21 0.16 -7.19 0.16 -7.16 0.16 -7.14 0.16 -7.12 0.16 -7.09 0.16 -7.07 0.16 -7.05 0.16 -7.02 0.16 -7 0.16 -6.98 0.16 -6.95 0.16 -6.93 0.16 -6.91 0.16 -6.88 0.16 -6.86 0.16 -6.84 0.16 -6.82 0.16 -6.79 0.16 -6.77 0.16 -6.75 0.16 -6.72 0.16 -6.7 0.16 -6.68 0.16 -6.65 0.16 -6.63 0.16 -6.61 0.16 -6.58 0.16 -6.56 0.16 -6.54 0.16 -6.52 0.16 -6.49 0.16 -6.47 0.16 -6.45 0.16 -6.42 0.16 -6.4 0.16 -6.38 0.16 -6.35 0.16 -6.33 0.16 -6.31 0.16 -6.28 0.16 -6.26 0.16 -6.24 0.16 -6.21 0.16 -6.19 0.16 -6.17 0.16 -6.15 0.16 -6.12 0.16 -6.1 0.16 -6.08 0.16 -6.05 0.16 -6.03 0.16 -6.01 0.16 -5.98 0.16 -5.96 0.16 -5.94 0.16 -5.91 0.16 -5.89 0.16 -5.87 0.16 -5.85 0.16 -5.82 0.16 -5.8 0.16 -5.78 0.16 -5.75 0.16 -5.73 0.16 -5.71 0.16 -5.68 0.16 -5.66 0.16 -5.64 0.16 -5.61 0.16 -5.59 0.16 -5.57 0.16 -5.54 0.16 -5.52 0.16 -5.5 0.16 -5.48 0.16 -5.45 0.16 -5.43 0.16 -5.41 0.16 -5.38 0.16 -5.36 0.16 -5.34 0.16 -5.31 0.16 -5.29 0.16 -5.27 0.16 -5.24 0.16 -5.22 0.16 -5.2 0.16 -5.18 0.16 -5.15 0.16 -5.13 0.16 -5.11 0.16 -5.08 0.16 -5.06 0.16 -5.04 0.16 -5.01 0.16 -4.99 0.16 -4.97 0.16 -4.94 0.16 -4.92 0.16 -4.9 0.16 -4.87 0.16 -4.85 0.16 -4.83 0.16 -4.81 0.16 -4.78 0.16 -4.76 0.16 -4.74 0.16 -4.71 0.16 -4.69 0.16 -4.67 0.16 -4.64 0.16 -4.62 0.16 -4.6 0.16 -4.57 0.16 -4.55 0.16 -4.53 0.16 -4.51 0.16 -4.48 0.16 -4.46 0.16 -4.44 0.16 -4.41 0.16 -4.39 0.16 -4.37 0.16 -4.34 0.16 -4.32 0.16 -4.3 0.16 -4.27 0.16 -4.25 0.16 -4.23 0.16 -4.2 0.16 -4.18 0.16 -4.16 0.16 -4.14 0.16 -4.11 0.16 -4.09 0.16 -4.07 0.16 -4.04 0.16 -4.02 0.16 -4 0.16 -3.97 0.16 -3.95 0.16 -3.93 0.16 -3.9 0.16 -3.88 0.16 -3.86 0.16 -3.84 0.16 -3.81 0.16 -3.79 0.16 -3.77 0.16 -3.74 0.16 -3.72 0.16 -3.7 0.16 -3.67 0.16 -3.65 0.16 -3.63 0.16 -3.6 0.16 -3.58 0.16 -3.56 0.16 -3.53 0.16 -3.51 0.16 -3.49 0.16 -3.47 0.16 -3.44 0.16 -3.42 0.16 -3.4 0.16 -3.37 0.16 -3.35 0.16 -3.33 0.16 -3.3 0.16 -3.28 0.16 -3.26 0.16 -3.23 0.16 -3.21 0.16 -3.19 0.16 -3.17 0.16 -3.14 0.16 -3.12 0.16 -3.1 0.16 -3.07 0.16 -3.05 0.16 -3.03 0.16 -3 0.16 -2.98 0.16 -2.96 0.16 -2.93 0.16 -2.91 0.16 -2.89 0.16 -2.86 0.16 -2.84 0.16 -2.82 0.16 -2.8 0.16 -2.77 0.16 -2.75 0.16 -2.73 0.16 -2.7 0.16 -2.68 0.16 -2.66 0.16 -2.63 0.16 -2.61 0.16 -2.59 0.16 -2.56 0.16 -2.54 0.16 -2.52 0.16 -2.5 0.16 -2.47 0.16 -2.45 0.16 -2.43 0.16 -2.4 0.16 -2.38 0.16 -2.36 0.16 -2.33 0.16 -2.31 0.16 -2.29 0.16 -2.26 0.16 -2.24 0.16 -2.22 0.16 -2.19 0.16 -2.17 0.16 -2.15 0.16 -2.13 0.16 -2.1 0.16 -2.08 0.16 -2.06 0.16 -2.03 0.16 -2.01 0.16 -1.99 0.16 -1.96 0.16 -1.94 0.16 -1.92 0.16 -1.89 0.16 -1.87 0.16 -1.85 0.15 -1.83 0.15 -1.8 0.15 -1.78 0.15 -1.76 0.15 -1.73 0.15 -1.71 0.15 -1.69 0.15 -1.66 0.15 -1.64 0.15 -1.62 0.15 -1.59 0.15 -1.57 0.15 -1.55 0.15 -1.52 0.15 -1.5 0.15 -1.48 0.15 -1.46 0.15 -1.43 0.15 -1.41 0.15 -1.39 0.15 -1.36 0.15 -1.34 0.15 -1.32 0.15 -1.29 0.15 -1.27 0.15 -1.25 0.15 -1.22 0.15 -1.2 0.15 -1.18 0.15 -1.16 0.15 -1.13 0.15 -1.11 0.15 -1.09 0.15 -1.06 0.15 -1.04 0.15 -1.02 0.15 -0.99 0.15 -0.97 0.15 -0.95 0.15 -0.92 0.15 -0.9 0.15 -0.88 0.15 -0.85 0.15 -0.83 0.15 -0.81 0.15 -0.79 0.15 -0.76 0.15 -0.74 0.15 -0.72 0.15 -0.69 0.15 -0.67 0.15 -0.65 0.15 -0.62 0.15 -0.6 0.15 -0.58 0.15 -0.55 0.15 -0.53 0.15 -0.51 0.15 -0.49 0.15 -0.46 0.15 -0.44 0.15 -0.42 0.15 -0.39 0.15 -0.37 0.15 -0.35 0.15 -0.32 0.15 -0.3 0.15 -0.28 0.15 -0.25 0.15 -0.23 0.15 -0.21 0.15 -0.18 0.15 -0.16 0.15 -0.14 0.15 -0.12 0.15 -0.09 0.15 -0.07 0.15 -0.05 0.15 -0.02 0.15 0 0.15 0.02 0.15 0.05 0.15 0.07 0.15 0.09 0.15 0.12 0.15 0.14 0.15 0.16 0.15 0.18 0.15 0.21 0.15 0.23 0.15 0.25 0.15 0.28 0.15 0.3 0.15 0.32 0.15 0.35 0.15 0.37 0.15 0.39 0.15 0.42 0.15 0.44 0.15 0.46 0.15 0.49 0.15 0.51 0.15 0.53 0.15 0.55 0.15 0.58 0.15 0.6 0.15 0.62 0.15 0.65 0.15 0.67 0.15 0.69 0.15 0.72 0.15 0.74 0.15 0.76 0.15 0.79 0.15 0.81 0.15 0.83 0.15 0.85 0.15 0.88 0.15 0.9 0.15 0.92 0.15 0.95 0.15 0.97 0.15 0.99 0.15 1.02 0.15 1.04 0.15 1.06 0.15 1.09 0.15 1.11 0.15 1.13 0.15 1.16 0.15 1.18 0.15 1.2 0.15 1.22 0.15 1.25 0.15 1.27 0.15 1.29 0.15 1.32 0.15 1.34 0.15 1.36 0.15 1.39 0.15 1.41 0.15 1.43 0.15 1.46 0.15 1.48 0.15 1.5 0.15 1.52 0.15 1.55 0.15 1.57 0.15 1.59 0.15 1.62 0.15 1.64 0.15 1.66 0.15 1.69 0.15 1.71 0.15 1.73 0.15 1.76 0.15 1.78 0.15 1.8 0.15 1.83 0.15 1.85 0.15 1.87 0.15 1.89 0.15 1.92 0.15 1.94 0.15 1.96 0.15 1.99 0.15 2.01 0.15 2.03 0.15 2.06 0.15 2.08 0.15 2.1 0.15 2.13 0.15 2.15 0.15 2.17 0.15 2.19 0.15 2.22 0.15 2.24 0.15 2.26 0.15 2.29 0.15 2.31 0.15 2.33 0.15 2.36 0.15 2.38 0.15 2.4 0.15 2.43 0.15 2.45 0.15 2.47 0.15 2.5 0.15 2.52 0.15 2.54 0.15 2.56 0.15 2.59 0.15 2.61 0.15 2.63 0.15 2.66 0.15 2.68 0.15 2.7 0.15 2.73 0.15 2.75 0.15 2.77 0.15 2.8 0.15 2.82 0.15 2.84 0.15 2.86 0.15 2.89 0.15 2.91 0.15 2.93 0.15 2.96 0.15 2.98 0.15 3 0.15 3.03 0.15 3.05 0.15 3.07 0.15 3.1 0.15 3.12 0.15 3.14 0.15 3.17 0.15 3.19 0.15 3.21 0.15 3.23 0.15 3.26 0.15 3.28 0.15 3.3 0.15 3.33 0.15 3.35 0.15 3.37 0.15 3.4 0.15 3.42 0.15 3.44 0.15 3.47 0.15 3.49 0.15 3.51 0.15 3.53 0.15 3.56 0.15 3.58 0.15 3.6 0.15 3.63 0.15 3.65 0.15 3.67 0.15 3.7 0.15 3.72 0.15 3.74 0.15 3.77 0.15 3.79 0.15 3.81 0.15 3.84 0.15 3.86 0.15 3.88 0.15 3.9 0.15 3.93 0.15 3.95 0.15 3.97 0.15 4 0.15 4.02 0.15 4.04 0.15 4.07 0.15 4.09 0.15 4.11 0.15 4.14 0.15 4.16 0.15 4.18 0.15 4.2 0.15 4.23 0.15 4.25 0.15 4.27 0.15 4.3 0.15 4.32 0.15 4.34 0.15 4.37 0.15 4.39 0.15 4.41 0.15 4.44 0.15 4.46 0.15 4.48 0.15 4.51 0.15 4.53 0.15 4.55 0.15 4.57 0.15 4.6 0.15 4.62 0.15 4.64 0.15 4.67 0.15 4.69 0.15 4.71 0.15 4.74 0.15 4.76 0.15 4.78 0.15 4.81 0.15 4.83 0.15 4.85 0.15 4.87 0.15 4.9 0.15 4.92 0.15 4.94 0.15 4.97 0.15 4.99 0.15 5.01 0.15 5.04 0.15 5.06 0.15 5.08 0.15 5.11 0.15 5.13 0.15 5.15 0.15 5.18 0.15 5.2 0.15 5.22 0.15 5.24 0.15 5.27 0.15 5.29 0.15 5.31 0.15 5.34 0.15 5.36 0.15 5.38 0.15 5.41 0.15 5.43 0.15 5.45 0.15 5.48 0.15 5.5 0.15 5.52 0.15 5.54 0.15 5.57 0.15 5.59 0.15 5.61 0.15 5.64 0.15 5.66 0.15 5.68 0.15 5.71 0.15 5.73 0.15 5.75 0.15 5.78 0.15 5.8 0.15 5.82 0.15 5.85 0.15 5.87 0.15 5.89 0.15 5.91 0.15 5.94 0.15 5.96 0.15 5.98 0.15 6.01 0.15 6.03 0.15 6.05 0.15 6.08 0.15 6.1 0.15 6.12 0.15 6.15 0.15 6.17 0.15 6.19 0.15 6.21 0.15 6.24 0.15 6.26 0.15 6.28 0.15 6.31 0.15 6.33 0.15 6.35 0.15 6.38 0.15 6.4 0.15 6.42 0.15 6.45 0.15 6.47 0.15 6.49 0.15 6.52 0.15 6.54 0.15 6.56 0.15 6.58 0.15 6.61 0.15 6.63 0.15 6.65 0.15 6.68 0.15 6.7 0.15 6.72 0.15 6.75 0.15 6.77 0.15 6.79 0.15 6.82 0.15 6.84 0.15 6.86 0.15 6.88 0.15 6.91 0.15 6.93 0.15 6.95 0.15 6.98 0.15 7 0.15 7.02 0.15 7.05 0.15 7.07 0.15 7.09 0.15 7.12 0.15 7.14 0.15 7.16 0.15 7.19 0.15 7.21 0.15 7.23 0.15 7.25 0.15 7.28 0.15 7.3 0.15 7.32 0.15 7.35 0.15 7.37 0.15 7.39 0.15 7.42 0.15 7.44 0.15 7.46 0.15 7.49 0.15 7.51 0.15 7.53 0.15 7.55 0.15 7.58 0.15 7.6 0.15 7.62 0.15 7.65 0.15 7.67 0.15 7.69 0.15 7.72 0.15 7.74 0.15 7.76 0.15 7.79 0.15 7.81 0.15 7.83 0.15 7.86 0.15 7.88 0.15 7.9 0.15 7.92 0.15 7.95 0.15 7.97 0.15 7.99 0.15 8.02 0.15 8.04 0.15 8.06 0.15 8.09 0.15 8.11 0.15 8.13 0.15 8.16 0.15 8.18 0.15 8.2 0.15 8.22 0.15 8.25 0.15 8.27 0.15 8.29 0.15 8.32 0.15 8.34 0.15 8.36 0.15 8.39 0.15 8.41 0.15 8.43 0.15 8.46 0.15 8.48 0.15 8.5 0.15 8.53 0.15 8.55 0.15 8.57 0.15 8.59 0.15 8.62 0.15 8.64 0.15 8.66 0.15 8.69 0.15 8.71 0.15 8.73 0.15 8.76 0.15 8.78 0.15 8.8 0.15 8.83 0.15 8.85 0.15 8.87 0.15 8.89 0.15 8.92 0.15 8.94 0.15 8.96 0.15 8.99 0.15 9.01 0.15 9.03 0.15 9.06 0.15 9.08 0.15 9.1 0.15 9.13 0.15 9.15 0.15 9.17 0.15 9.2 0.15 9.22 0.15 9.24 0.15 9.26 0.15 9.29 0.15 9.31 0.15 9.33 0.15 9.36 0.15 9.38 0.15 9.4 0.15 9.43 0.15 9.45 0.15 9.47 0.15 9.5 0.15 9.52 0.15 9.54 0.15 9.56 0.15 9.59 0.15 9.61 0.15 9.63 0.15 9.66 0.15 9.68 0.15 9.7 0.15 9.73 0.15 9.75 0.15 9.77 0.15 9.8 0.15 9.82 0.15 9.84 0.15 9.87 0.15 9.89 0.15 9.91 0.15 9.93 0.15 9.96 0.15 9.98 0.15 10 0.15 10.03 0.15 10.05 0.15 10.07 0.15 10.1 0.15 10.12 0.15 10.14 0.15 10.17 0.15 10.19 0.15 10.21 0.15 10.23 0.15 10.26 0.15 10.28 0.15 10.3 0.15 10.33 0.15 10.35 0.15 10.37 0.15 10.4 0.15 10.42 0.15 10.44 0.15 10.47 0.15 10.49 0.15 10.51 0.15 10.54 0.15 10.56 0.15 10.58 0.15 10.6 0.15 10.63 0.15 10.65 0.15 10.67 0.15 10.7 0.15 10.72 0.15 10.74 0.15 10.77 0.15 10.79 0.15 10.81 0.15 10.84 0.15 10.86 0.15 10.88 0.15 10.9 0.15 10.93 0.15 10.95 0.15 10.97 0.15 11 0.15 11.02 0.15 11.04 0.15 11.07 0.15 11.09 0.15 11.11 0.15 11.14 0.15 11.16 0.15 11.18 0.15 11.21 0.15 11.23 0.15 11.25 0.15 11.27 0.15 11.3 0.15 11.32 0.15 11.34 0.15 11.37 0.15 11.39 0.15 11.41 0.15 11.44 0.15 11.46 0.15 11.48 0.15 11.51 0.15 11.53 0.15 11.55 0.15 11.57 0.15 11.6 0.15 11.62 0.15 11.64 0.15 11.67 0.15 11.69 0.15 11.71 0.15 11.74 0.15 11.76 0.15 11.78 0.15 11.81 0.15 11.83 0.15 11.85 0.15 11.88 0.15 11.9 0.15 11.92 0.15 11.94 0.15 11.97 0.15 11.99 0.15 12.01 0.15 12.04 0.15 12.06 0.15 12.08 0.15 12.11 0.15 12.13 0.15 12.15 0.15 12.18 0.15 12.2 0.15 12.22 0.15 12.24 0.15 12.27 0.15 12.29 0.15 12.31 0.15 12.34 0.15 12.36 0.15 12.38 0.15 12.41 0.15 12.43 0.15 12.45 0.15 12.48 0.15 12.5 0.15 12.52 0.15 12.55 0.15 12.57 0.15 12.59 0.15 12.61 0.15 12.64 0.15 12.66 0.15 12.68 0.15 12.71 0.15 12.73 0.15 12.75 0.15 12.78 0.15 12.8 0.15 12.82 0.15 12.85 0.15 12.87 0.15 12.89 0.15 12.91 0.15 12.94 0.15 12.96 0.15 12.98 0.15 13.01 0.15 13.03 0.15 13.05 0.15 13.08 0.15 13.1 0.15 13.12 0.15 13.15 0.15 13.17 0.15 13.19 0.15 13.22 0.15 13.24 0.15 13.26 0.15 13.28 0.15 13.31 0.15 13.33 0.15 13.35 0.15 13.38 0.15 13.4 0.15 13.42 0.15 13.45 0.15 13.47 0.15 13.49 0.15 13.52 0.15 13.54 0.15 13.56 0.15 13.58 0.15 13.61 0.15 13.63 0.15 13.65 0.15 13.68 0.15 13.7 0.15 13.72 0.15 13.75 0.15 13.77 0.15 13.79 0.15 13.82 0.15 13.84 0.15 13.86 0.15 13.89 0.15 13.91 0.15 13.93 0.15 13.95 0.15 13.98 0.15 14 0.15 14.02 0.15 14.05 0.15 14.07 0.15 14.09 0.15 14.12 0.15 14.14 0.15 14.16 0.15 14.19 0.15 14.21 0.15 14.23 0.15 14.25 0.15 14.28 0.15 14.3 0.15 14.32 0.15 14.35 0.15 14.37 0.15 14.39 0.15 14.42 0.15 14.44 0.15 14.46 0.15 14.49 0.15 14.51 0.15 14.53 0.15 14.56 0.15 14.58 0.15 14.6 0.15 14.62 0.15 14.65 0.15 14.67 0.15 14.69 0.15 14.72 0.15 14.74 0.15 14.76 0.15 14.79 0.15 14.81 0.15 14.83 0.15 14.86 0.15 14.88 0.15 14.9 0.15 14.92 0.15 14.95 0.15 14.97 0.15 14.99 0.15 15.02 0.15 15.04 0.15 15.06 0.15 15.09 0.15 15.11 0.15 15.13 0.15 15.16 0.15 15.18 0.15 15.2 0.15 15.23 0.15 15.25 0.15 15.27 0.15 15.29 0.15 15.32 0.15 15.34 0.15 15.36 0.15 15.39 0.15 15.41 0.15 15.43 0.15 15.46 0.15 15.48 0.15 15.5 0.15 15.53 0.15 15.55 0.15 15.57 0.15 15.59 0.15 15.62 0.15 15.64 0.15 15.66 0.15 15.69 0.15 15.71 0.15 15.73 0.15 15.76 0.15 15.78 0.15 15.8 0.15 15.83 0.15 15.85 0.15 15.87 0.15 15.9 0.15 15.92 0.15 15.94 0.15 15.96 0.15 15.99 0.15 16.01 0.15 16.03 0.15 16.06 0.15 16.08 0.15 16.1 0.15 16.13 0.15 16.15 0.15 16.17 0.15 16.2 0.15 16.22 0.15 16.24 0.15 16.26 0.15 16.29 0.15 16.31 0.15 16.33 0.15 16.36 0.15 16.38 0.15 16.4 0.15 16.43 0.15 16.45 0.15 16.47 0.15 16.5 0.15 16.52 0.15 16.54 0.15 16.57 0.15 16.59 0.15 16.61 0.15 16.63 0.15 16.66 0.15 16.68 0.15 16.7 0.15 16.73 0.15 16.75 0.15 16.77 0.15 16.8 0.15 16.82 0.15 16.84 0.15 16.87 0.15 16.89 0.15 16.91 0.15 16.93 0.15 16.96 0.15 16.98 0.15 17 0.15 17.03 0.15 17.05 0.15 17.07 0.15 17.1 0.15 17.12 0.15 17.14 0.15 17.17 0.15 17.19 0.15 17.21 0.15 17.24 0.15 17.26 0.15 17.28 0.15 17.3 0.15 17.33 0.15 17.35 0.15 17.37 0.15 17.4 0.15 17.42 0.15 17.44 0.15 17.47 0.15 17.49 0.15 17.51 0.15 17.54 0.15 17.56 0.15 17.58 0.15 17.6 0.15 17.63 0.15 17.65 0.15 17.67 0.15 17.7 0.15 17.72 0.15 17.74 0.15 17.77 0.15 17.79 0.15 17.81 0.15 17.84 0.15 17.86 0.15 17.88 0.15 17.91 0.15 17.93 0.15 17.95 0.15 17.97 0.15 18 0.15 18.02 0.15 18.04 0.15 18.07 0.15 18.09 0.15 18.11 0.15 18.14 0.15 18.16 0.15 18.18 0.15 18.21 0.15 18.23 0.15 18.25 0.15 18.27 0.15 18.3 0.15 18.32 0.15 18.34 0.15 18.37 0.15 18.39 0.15 18.41 0.15 18.44 0.15 18.46 0.15 18.48 0.15 18.51 0.15 18.53 0.15 18.55 0.15 18.58 0.15 18.6 0.15 18.62 0.15 18.64 0.15 18.67 0.15 18.69 0.15 18.71 0.15 18.74 0.15 18.76 0.15 18.78 0.15 18.81 0.15 18.83 0.15 18.85 0.15 18.88 0.15 18.9 0.15 18.92 0.15 18.94 0.15 18.97 0.15 18.99 0.15 19.01 0.15 19.04 0.15 19.06 0.15 19.08 0.15 19.11 0.15 19.13 0.15 19.15 0.15 19.18 0.15 19.2 0.15 19.22 0.15 19.25 0.15 19.27 0.15 19.29 0.15 19.31 0.15 19.34 0.15 19.36 0.15 19.38 0.15 19.41 0.15 19.43 0.15 19.45 0.15 19.48 0.15 19.5 0.15 19.52 0.15 19.55 0.15 19.57 0.15 19.59 0.15 19.61 0.15 19.64 0.15 19.66 0.15 19.68 0.15 19.71 0.15 19.73 0.15 19.75 0.15 19.78 0.15 19.8 0.15 19.82 0.15 19.85 0.15 19.87 0.15 19.89 0.15 19.92 0.15 19.94 0.15 19.96 0.15 19.98 0.15 20.01 0.15 20.03 0.15 20.05 0.15 20.08 0.15 20.1 0.15 20.12 0.15 20.15 0.15 20.17 0.15 20.19 0.15 20.22 0.15 20.24 0.15 20.26 0.15 20.28 0.15 20.31 0.15 20.33 0.15 20.35 0.15 20.38 0.15 20.4 0.15 20.42 0.15 20.45 0.15 20.47 0.15 20.49 0.15 20.52 0.15 20.54 0.15 20.56 0.15 20.59 0.15 20.61 0.15 20.63 0.15 20.65 0.15 20.68 0.15 20.7 0.15 20.72 0.15 20.75 0.15 20.77 0.15 20.79 0.15 20.82 0.15 20.84 0.15 20.86 0.15 20.89 0.15 20.91 0.15 20.93 0.15 20.95 0.15 20.98 0.15 21 0.15 21.02 0.15 21.05 0.15 21.07 0.15 21.09 0.15 21.12 0.15 21.14 0.15 21.16 0.15 21.19 0.15 21.21 0.15 21.23 0.15 21.26 0.15 21.28 0.15 21.3 0.15 21.32 0.15 21.35 0.15 21.37 0.15 21.39 0.15 21.42 0.15 21.44 0.15 21.46 0.15 21.49 0.15 21.51 0.15 21.53 0.15 21.56 0.15 21.58 0.15 21.6 0.15 21.62 0.15 21.65 0.15 21.67 0.15 21.69 0.15 21.72 0.15 21.74 0.15 21.76 0.15 21.79 0.15 21.81 0.15 21.83 0.15 21.86 0.15 21.88 0.15 21.9 0.15 21.93 0.15 21.95 0.15 21.97 0.15 21.99 0.15 22.02 0.15 22.04 0.15 22.06 0.15 22.09 0.15 22.11 0.15 22.13 0.15 22.16 0.15 22.18 0.15 22.2 0.15 22.23 0.15 22.25 0.15 22.27 0.15 22.29 0.15 22.32 0.15 22.34 0.15 22.36 0.15 22.39 0.15 22.41 0.15 22.43 0.15 22.46 0.15 22.48 0.15 22.5 0.15 22.53 0.15 22.55 0.15 22.57 0.15 22.6 0.15 22.62 0.15 22.64 0.15 22.66 0.15 22.69 0.15 22.71 0.15 22.73 0.15 22.76 0.15 22.78 0.15 22.8 0.15 22.83 0.15 22.85 0.15 22.87 0.15 22.9 0.15 22.92 0.15 22.94 0.15 22.96 0.15 22.99 0.15 23.01 0.15 23.03 0.15 23.06 0.15 23.08 0.15 23.1 0.15 23.13 0.15 23.15 0.15 23.17 0.15 23.2 0.15 23.22 0.15 23.24 0.15 23.27 0.15 23.29 0.15 23.31 0.15 23.33 0.15 23.36 0.15 23.38 0.15 23.4 0.15 23.43 0.15 23.45 0.15 23.47 0.15 23.5 0.15 23.52 0.15 23.54 0.15 23.57 0.15 23.59 0.15 23.61 0.15 23.63 0.15 23.66 0.15 23.68 0.15 23.7 0.15 23.73 0.15 23.75 0.15 23.77 0.15 23.8 0.15 23.82 0.15 23.84 0.15 23.87 0.15 23.89 0.15 23.91 0.15 23.94 0.15 23.96 0.15 23.98 0.15 24 0.15 24.03 0.15 24.05 0.15 24.07 0.15 24.1 0.15 24.12 0.15 24.14 0.15 24.17 0.15 24.19 0.15 24.21 0.15 24.24 0.15 24.26 0.15 24.28 0.15 24.3 0.15 24.33 0.15 24.35 0.15 24.37 0.15 24.4 0.15 24.42 0.15 24.44 0.15 24.47 0.15 24.49 0.15 24.51 0.15 24.54 0.15 24.56 0.15 24.58 0.15 24.61 0.15 24.63 0.15 24.65 0.15 24.67 0.15 24.7 0.15 24.72 0.15 24.74 0.15 24.77 0.15 24.79 0.15 24.81 0.15 24.84 0.15 24.86 0.15 24.88 0.15 24.91 0.15 24.93 0.15 24.95 0.15 24.97 0.15 25 0.15 25.02 0.15 25.04 0.15 25.07 0.15 25.09 0.15 25.11 0.15 25.14 0.15 25.16 0.15 25.18 0.15 25.21 0.15 25.23 0.15 25.25 0.15 25.28 0.15 25.3 0.15 25.32 0.15 25.34 0.15 25.37 0.15 25.39 0.15 25.41 0.15 25.44 0.15 25.46 0.15 25.48 0.15 25.51 0.15 25.53 0.15 25.55 0.15 25.58 0.15 25.6 0.15 25.62 0.15 25.64 0.15 25.67 0.15 25.69 0.15 25.71 0.15 25.74 0.15 25.76 0.15 25.78 0.15 25.81 0.15 25.83 0.15 25.85 0.15 25.88 0.15 25.9 0.15 25.92 0.15 25.95 0.15 25.97 0.15 25.99 0.15 26.01 0.15 26.04 0.15 26.06 0.15 26.08 0.15 26.11 0.15 26.13 0.15 26.15 0.15 26.18 0.15 26.2 0.15 26.22 0.15 26.25 0.15 26.27 0.15 26.29 0.15 26.31 0.15 26.34 0.15 26.36 0.15 26.38 0.15 26.41 0.15 26.43 0.15 26.45 0.15 26.48 0.15 26.5 0.15 26.52 0.15 26.55 0.15 26.57 0.15 26.59 0.15 26.62 0.15 26.64 0.15 26.66 0.15 26.68 0.15 26.71 0.15 26.73 0.15 26.75 0.15 26.78 0.15 26.8 0.15 26.82 0.15 26.85 0.15 26.87 0.15 26.89 0.15 26.92 0.15 26.94 0.15 26.96 0.15 26.98 0.15 27.01 0.15 27.03 0.15 27.05 0.15 27.08 0.15 27.1 0.15 27.12 0.15 27.15 0.15 27.17 0.15 27.19 0.15 27.22 0.15 27.24 0.15 27.26 0.15 27.29 0.15 27.31 0.15 27.33 0.15 27.35 0.15 27.38 0.15 27.4 0.15 27.42 0.15 27.45 0.15 27.47 0.15 27.49 0.15 27.52 0.15 27.54 0.15 27.56 0.15 27.59 0.15 27.61 0.15 27.63 0.15 27.65 0.15 27.68 0.15 27.7 0.15 27.72 0.15 27.75 0.15 27.77 0.15 27.79 0.15 27.82 0.15 27.84 0.15 27.86 0.15 27.89 0.15 27.91 0.15 27.93 0.15 27.96 0.15 27.98 0.15 28 0.15 28.02 0.15 28.05 0.15 28.07 0.15 28.09 0.15 28.12 0.15 28.14 0.15 28.16 0.15 28.19 0.15 28.21 0.15 28.23 0.15 28.26 0.15 28.28 0.15 28.3 0.15 28.32 0.15 28.35 0.15 28.37 0.15 28.39 0.15 28.42 0.15 28.44 0.15 28.46 0.15 28.49 0.15 28.51 0.15 28.53 0.15 28.56 0.15 28.58 0.15 28.6 0.15 28.63 0.15 28.65 0.15 28.67 0.15 28.69 0.15 28.72 0.15 28.74 0.15 28.76 0.15 28.79 0.15 28.81 0.15 28.83 0.15 28.86 0.15 28.88 0.15 28.9 0.15 28.93 0.15 28.95 0.15 28.97 0.15 28.99 0.15 29.02 0.15 29.04 0.15 29.06 0.15 29.09 0.15 29.11 0.15 29.13 0.15 29.16 0.15 29.18 0.15 29.2 0.15 29.23 0.15 29.25 0.15 29.27 0.15 29.3 0.15 29.32 0.15 29.34 0.15 29.36 0.15 29.39 0.15 29.41 0.15 29.43 0.15 29.46 0.15 29.48 0.15 29.5 0.15 29.53 0.15 29.55 0.15 29.57 0.15 29.6 0.15 29.62 0.15 29.64 0.15 29.66 0.15 29.69 0.15 29.71 0.15 29.73 0.15 29.76 0.15 29.78 0.15 29.8 0.15 29.83 0.15 29.85 0.15 29.87 0.15 29.9 0.15 29.92 0.15 29.94 0.15 29.97 0.15 29.99 0.15 30.01 0.15 30.03 0.15 30.06 0.15 30.08 0.15 30.1 0.15 30.13 0.15 30.15 0.15 30.17 0.15 30.2 0.15 30.22 0.15 30.24 0.15 30.27 0.15 30.29 0.15 30.31 0.15 30.33 0.15 30.36 0.15 30.38 0.15 30.4 0.15 30.43 0.15 30.45 0.15 30.47 0.15 30.5 0.15 30.52 0.15 30.54 0.15 30.57 0.15 30.59 0.15 30.61 0.15 30.64 0.15 30.66 0.15 30.68 0.15 30.7 0.15 30.73 0.15 30.75 0.15 30.77 0.15 30.8 0.15 30.82 0.15 30.84 0.15 30.87 0.15 30.89 0.15 30.91 0.15 30.94 0.15 30.96 0.15 30.98 0.15 31 0.15 31.03 0.15 31.05 0.15 31.07 0.15 31.1 0.15 31.12 0.15 31.14 0.15 31.17 0.15 31.19 0.15 31.21 0.15 31.24 0.15 31.26 0.15 31.28 0.15 31.31 0.15 31.33 0.15 31.35 0.15 31.37 0.15 31.4 0.15 31.42 0.15 31.44 0.15 31.47 0.15 31.49 0.15 31.51 0.15 31.54 0.15 31.56 0.15 31.58 0.15 31.61 0.15 31.63 0.15 31.65 0.15 31.67 0.15 31.7 0.15 31.72 0.15 31.74 0.15 31.77 0.15 31.79 0.15 31.81 0.15 31.84 0.15 31.86 0.15 31.88 0.15 31.91 0.15 31.93 0.15 31.95 0.15 31.98 0.15 32 0.15 32.02 0.15 32.04 0.15 32.07 0.15 32.09 0.15 32.11 0.15 32.14 0.15 32.16 0.15 32.18 0.15 32.21 0.15 32.23 0.15 32.25 0.15 32.28 0.15 32.3 0.15 32.32 0.15 32.34 0.15 32.37 0.15 32.39 0.15 32.41 0.15 32.44 0.15 32.46 0.15 32.48 0.15 32.51 0.15 32.53 0.15 32.55 0.15 32.58 0.15 32.6 0.15 32.62 0.15 32.65 0.15 32.67 0.15 32.69 0.15 32.71 0.15 32.74 0.15 32.76 0.15 32.78 0.15 32.81 0.15 32.83 0.15 32.85 0.15 32.88 0.15 32.9 0.15 32.92 0.15 32.95 0.15 32.97 0.15 32.99 0.15 33.01 0.15 33.04 0.15 33.06 0.15 33.08 0.15 33.11 0.15 33.13 0.15 33.15 0.15 33.18 0.15 33.2 0.15 33.22 0.15 33.25 0.15 33.27 0.15 33.29 0.15 33.32 0.15 33.34 0.15 33.36 0.15 33.38 0.15 33.41 0.15 33.43 0.15 33.45 0.15 33.48 0.15 33.5 0.15 33.52 0.15 33.55 0.15 33.57 0.15 33.59 0.15 33.62 0.15 33.64 0.15 33.66 0.15 33.68 0.15 33.71 0.15 33.73 0.15 33.75 0.15 33.78 0.15 33.8 0.15 33.82 0.15 33.85 0.15 33.87 0.15 33.89 0.15 33.92 0.15 33.94 0.15 33.96 0.15 33.99 0.15 34.01 0.15 34.03 0.15 34.05 0.15 34.08 0.15 34.1 0.15 34.12 0.15 34.15 0.15 34.17 0.15 34.19 0.15 34.22 0.15 34.24 0.15 34.26 0.15 34.29 0.15 34.31 0.15 34.33 0.15 34.35 0.15 34.38 0.15 34.4 0.15 34.42 0.15 34.45 0.15 34.47 0.15 34.49 0.15 34.52 0.15 34.54 0.15 34.56 0.15 34.59 0.15 34.61 0.15 34.63 0.15 34.66 0.15 34.68 0.15 34.7 0.15 34.72 0.15 34.75 0.15 34.77 0.15 34.79 0.15 34.82 0.15 34.84 0.15 34.86 0.15 34.89 0.15 34.91 0.15 34.93 0.15 34.96 0.15 34.98 0.15 35 0.15 35.02 0.15 35.05 0.15 35.07 0.15 35.09 0.15 35.12 0.15 35.14 0.15 35.16 0.15 35.19 0.15 35.21 0.15 35.23 0.15 35.26 0.15 35.28 0.15 35.3 0.15 35.33 0.15 35.35 0.15 35.37 0.15 35.39 0.15 35.42 0.15 35.44 0.15 35.46 0.15 35.49 0.15 35.51 0.15 35.53 0.15 35.56 0.15 35.58 0.15 35.6 0.15 35.63 0.15 35.65 0.15 35.67 0.15 35.69 0.15 35.72 0.15 35.74 0.15 35.76 0.15 35.79 0.15 35.81 0.15 35.83 0.15 35.86 0.15 35.88 0.15 35.9 0.15 35.93 0.15 35.95 0.15 35.97 0.15 36 0.15 36.02 0.15 36.04 0.15 36.06 0.15 36.09 0.15 36.11 0.15 36.13 0.15 36.16 0.15 36.18 0.15 36.2 0.15 36.23 0.15 36.25 0.15 36.27 0.15 36.3 0.15 36.32 0.15 36.34 0.15 36.36 0.15 36.39 0.15 36.41 0.15 36.43 0.15 36.46 0.15 36.48 0.15 36.5 0.15 36.53 0.15 36.55 0.15 36.57 0.15 36.6 0.15 36.62 0.15 36.64 0.15 36.67 0.15 36.69 0.15 36.71 0.15 36.73 0.15 36.76 0.15 36.78 0.15 36.8 0.15 36.83 0.15 36.85 0.15 36.87 0.15 36.9 0.15 36.92 0.15 36.94 0.15 36.97 0.15 36.99 0.15 37.01 0.15 37.03 0.15 37.06 0.15 37.08 0.15 37.1 0.15 37.13 0.15 37.15 0.15 37.17 0.15 37.2 0.15 37.22 0.15 37.24 0.15 37.27 0.15 37.29 0.15 37.31 0.15 37.34 0.15 37.36 0.15 37.38 0.15 37.4 0.15 37.43 0.15 37.45 0.15 37.47 0.15 37.5 0.15 37.52 0.15 37.54 0.15 37.57 0.15 37.59 0.15 37.61 0.15 37.64 0.15 37.66 0.15 37.68 0.15 37.7 0.15 37.73 0.15 37.75 0.15 37.77 0.15 37.8 0.15 37.82 0.15 37.84 0.15 37.87 0.15 37.89 0.15 37.91 0.15 37.94 0.15 37.96 0.15 37.98 0.15 38.01 0.15 38.03 0.15 38.05 0.15 38.07 0.15 38.1 0.15 38.12 0.15 38.14 0.15 38.17 0.15 38.19 0.15 38.21 0.15 38.24 0.15 38.26 0.15 38.28 0.15 38.31 0.15 38.33 0.15 38.35 0.15 38.37 0.15 38.4 0.15 38.42 0.15 38.44 0.15 38.47 0.15 38.49 0.15 38.51 0.15 38.54 0.15 38.56 0.15 38.58 0.15 38.61 0.15 38.63 0.15 38.65 0.15 38.68 0.15 38.7 0.15 38.72 0.15 38.74 0.15 38.77 0.15 38.79 0.15 38.81 0.15 38.84 0.15 38.86 0.15 38.88 0.15 38.91 0.15 38.93 0.15 38.95 0.15 38.98 0.15 39 0.15 39.02 0.15 39.04 0.15 39.07 0.15 39.09 0.15 39.11 0.15 39.14 0.15 39.16 0.15 39.18 0.15 39.21 0.15 39.23 0.15 39.25 0.15 39.28 0.15 39.3 0.15 39.32 0.15 39.35 0.15 39.37 0.15 39.39 0.15 39.41 0.15 39.44 0.15 39.46 0.15 39.48 0.15 39.51 0.15 39.53 0.15 39.55 0.15 39.58 0.15 39.6 0.15 39.62 0.15 39.65 0.15 39.67 0.15 39.69 0.15 39.71 0.15 39.74 0.15 39.76 0.15 39.78 0.15 39.81 0.15 39.83 0.15 39.85 0.15 39.88 0.15 39.9 0.15 39.92 0.15 39.95 0.15 39.97 0.15 39.99 0.15 40.02 0.15 40.04 0.15 40.06 0.15 40.08 0.15 40.11 0.15 40.13 0.15 40.15 0.15 40.18 0.15 40.2 0.15 40.22 0.15 40.25 0.15 40.27 0.15 40.29 0.15 40.32 0.15 40.34 0.15 40.36 0.15 40.38 0.15 40.41 0.15 40.43 0.15 40.45 0.15 40.48 0.15 40.5 0.15 40.52 0.15 40.55 0.15 40.57 0.15 40.59 0.15 40.62 0.15 40.64 0.15 40.66 0.15 40.69 0.15 40.71 0.15 40.73 0.15 40.75 0.15 40.78 0.15 40.8 0.15 40.82 0.15 40.85 0.15 40.87 0.15 40.89 0.15 40.92 0.15 40.94 0.15 40.96 0.15 40.99 0.15 41.01 0.15 41.03 0.15 41.05 0.15 41.08 0.15 41.1 0.15 41.12 0.15 41.15 0.15 41.17 0.15 41.19 0.15 41.22 0.15 41.24 0.15 41.26 0.15 41.29 0.15 41.31 0.15 41.33 0.15 41.36 0.15 41.38 0.15 41.4 0.15 41.42 0.15 41.45 0.15 41.47 0.15 41.49 0.15 41.52 0.15 41.54 0.15 41.56 0.15 41.59 0.15 41.61 0.15 41.63 0.15 41.66 0.15 41.68 0.15 41.7 0.15 41.72 0.15 41.75 0.15 41.77 0.15 41.79 0.15 41.82 0.15 41.84 0.15 41.86 0.15 41.89 0.15 41.91 0.15 41.93 0.15 41.96 0.15 41.98 0.15 42 0.15 42.03 0.15 42.05 0.15 42.07 0.15 42.09 0.15 42.12 0.15 42.14 0.15 42.16 0.15" class="primitive"/>
          </g>
          <g transform="translate(70.89,77.25)" id="img-8311e7c5-235" class="geometry color_I_V" stroke="#FF6DAE">
            <path fill="none" d="M-42.16,1.47 L -42.14 1.47 -42.12 1.47 -42.09 1.47 -42.07 1.47 -42.05 1.47 -42.03 1.47 -42 1.47 -41.98 1.47 -41.96 1.47 -41.93 1.47 -41.91 1.47 -41.89 1.47 -41.86 1.47 -41.84 1.47 -41.82 1.47 -41.79 1.47 -41.77 1.47 -41.75 1.47 -41.72 1.46 -41.7 1.46 -41.68 1.46 -41.66 1.46 -41.63 1.46 -41.61 1.46 -41.59 1.46 -41.56 1.46 -41.54 1.46 -41.52 1.46 -41.49 1.46 -41.47 1.46 -41.45 1.46 -41.42 1.46 -41.4 1.46 -41.38 1.46 -41.36 1.46 -41.33 1.46 -41.31 1.46 -41.29 1.46 -41.26 1.46 -41.24 1.46 -41.22 1.46 -41.19 1.46 -41.17 1.46 -41.15 1.46 -41.12 1.46 -41.1 1.46 -41.08 1.46 -41.05 1.46 -41.03 1.46 -41.01 1.46 -40.99 1.45 -40.96 1.45 -40.94 1.45 -40.92 1.45 -40.89 1.45 -40.87 1.45 -40.85 1.45 -40.82 1.45 -40.8 1.45 -40.78 1.45 -40.75 1.45 -40.73 1.45 -40.71 1.45 -40.69 1.45 -40.66 1.45 -40.64 1.45 -40.62 1.45 -40.59 1.45 -40.57 1.45 -40.55 1.45 -40.52 1.45 -40.5 1.45 -40.48 1.44 -40.45 1.44 -40.43 1.44 -40.41 1.44 -40.38 1.44 -40.36 1.44 -40.34 1.44 -40.32 1.44 -40.29 1.44 -40.27 1.44 -40.25 1.44 -40.22 1.44 -40.2 1.44 -40.18 1.44 -40.15 1.44 -40.13 1.44 -40.11 1.44 -40.08 1.44 -40.06 1.43 -40.04 1.43 -40.02 1.43 -39.99 1.43 -39.97 1.43 -39.95 1.43 -39.92 1.43 -39.9 1.43 -39.88 1.43 -39.85 1.43 -39.83 1.43 -39.81 1.43 -39.78 1.43 -39.76 1.43 -39.74 1.43 -39.71 1.42 -39.69 1.42 -39.67 1.42 -39.65 1.42 -39.62 1.42 -39.6 1.42 -39.58 1.42 -39.55 1.42 -39.53 1.42 -39.51 1.42 -39.48 1.42 -39.46 1.41 -39.44 1.41 -39.41 1.41 -39.39 1.41 -39.37 1.41 -39.35 1.41 -39.32 1.41 -39.3 1.41 -39.28 1.41 -39.25 1.41 -39.23 1.4 -39.21 1.4 -39.18 1.4 -39.16 1.4 -39.14 1.4 -39.11 1.4 -39.09 1.4 -39.07 1.4 -39.04 1.39 -39.02 1.39 -39 1.39 -38.98 1.39 -38.95 1.39 -38.93 1.39 -38.91 1.39 -38.88 1.38 -38.86 1.38 -38.84 1.38 -38.81 1.38 -38.79 1.38 -38.77 1.38 -38.74 1.37 -38.72 1.37 -38.7 1.37 -38.68 1.37 -38.65 1.37 -38.63 1.37 -38.61 1.36 -38.58 1.36 -38.56 1.36 -38.54 1.36 -38.51 1.36 -38.49 1.35 -38.47 1.35 -38.44 1.35 -38.42 1.35 -38.4 1.35 -38.37 1.34 -38.35 1.34 -38.33 1.34 -38.31 1.34 -38.28 1.33 -38.26 1.33 -38.24 1.33 -38.21 1.33 -38.19 1.32 -38.17 1.32 -38.14 1.32 -38.12 1.32 -38.1 1.31 -38.07 1.31 -38.05 1.31 -38.03 1.3 -38.01 1.3 -37.98 1.3 -37.96 1.29 -37.94 1.29 -37.91 1.29 -37.89 1.28 -37.87 1.28 -37.84 1.28 -37.82 1.27 -37.8 1.27 -37.77 1.27 -37.75 1.26 -37.73 1.26 -37.7 1.26 -37.68 1.25 -37.66 1.25 -37.64 1.24 -37.61 1.24 -37.59 1.24 -37.57 1.23 -37.54 1.23 -37.52 1.22 -37.5 1.22 -37.47 1.21 -37.45 1.21 -37.43 1.2 -37.4 1.2 -37.38 1.19 -37.36 1.19 -37.34 1.18 -37.31 1.18 -37.29 1.17 -37.27 1.17 -37.24 1.16 -37.22 1.16 -37.2 1.15 -37.17 1.15 -37.15 1.14 -37.13 1.14 -37.1 1.13 -37.08 1.12 -37.06 1.12 -37.03 1.11 -37.01 1.1 -36.99 1.1 -36.97 1.09 -36.94 1.09 -36.92 1.08 -36.9 1.07 -36.87 1.06 -36.85 1.06 -36.83 1.05 -36.8 1.04 -36.78 1.04 -36.76 1.03 -36.73 1.02 -36.71 1.01 -36.69 1.01 -36.67 1 -36.64 0.99 -36.62 0.98 -36.6 0.97 -36.57 0.96 -36.55 0.96 -36.53 0.95 -36.5 0.94 -36.48 0.93 -36.46 0.92 -36.43 0.91 -36.41 0.9 -36.39 0.89 -36.36 0.88 -36.34 0.87 -36.32 0.86 -36.3 0.85 -36.27 0.84 -36.25 0.83 -36.23 0.82 -36.2 0.81 -36.18 0.8 -36.16 0.79 -36.13 0.78 -36.11 0.77 -36.09 0.76 -36.06 0.75 -36.04 0.73 -36.02 0.72 -36 0.71 -35.97 0.7 -35.95 0.69 -35.93 0.67 -35.9 0.66 -35.88 0.65 -35.86 0.64 -35.83 0.62 -35.81 0.61 -35.79 0.6 -35.76 0.58 -35.74 0.57 -35.72 0.56 -35.69 0.54 -35.67 0.53 -35.65 0.52 -35.63 0.5 -35.6 0.49 -35.58 0.47 -35.56 0.46 -35.53 0.44 -35.51 0.43 -35.49 0.41 -35.46 0.4 -35.44 0.38 -35.42 0.37 -35.39 0.35 -35.37 0.34 -35.35 0.32 -35.33 0.3 -35.3 0.29 -35.28 0.27 -35.26 0.25 -35.23 0.24 -35.21 0.22 -35.19 0.2 -35.16 0.19 -35.14 0.17 -35.12 0.15 -35.09 0.14 -35.07 0.12 -35.05 0.1 -35.02 0.08 -35 0.06 -34.98 0.05 -34.96 0.03 -34.93 0.01 -34.91 -0.01 -34.89 -0.03 -34.86 -0.05 -34.84 -0.07 -34.82 -0.08 -34.79 -0.1 -34.77 -0.12 -34.75 -0.14 -34.72 -0.16 -34.7 -0.18 -34.68 -0.2 -34.66 -0.22 -34.63 -0.24 -34.61 -0.26 -34.59 -0.28 -34.56 -0.3 -34.54 -0.32 -34.52 -0.34 -34.49 -0.36 -34.47 -0.38 -34.45 -0.4 -34.42 -0.42 -34.4 -0.44 -34.38 -0.46 -34.35 -0.48 -34.33 -0.5 -34.31 -0.52 -34.29 -0.55 -34.26 -0.57 -34.24 -0.59 -34.22 -0.61 -34.19 -0.63 -34.17 -0.65 -34.15 -0.67 -34.12 -0.69 -34.1 -0.71 -34.08 -0.74 -34.05 -0.76 -34.03 -0.78 -34.01 -0.8 -33.99 -0.82 -33.96 -0.84 -33.94 -0.86 -33.92 -0.89 -33.89 -0.91 -33.87 -0.93 -33.85 -0.95 -33.82 -0.97 -33.8 -0.99 -33.78 -1.01 -33.75 -1.03 -33.73 -1.06 -33.71 -1.08 -33.68 -1.1 -33.66 -1.12 -33.64 -1.14 -33.62 -1.16 -33.59 -1.18 -33.57 -1.2 -33.55 -1.23 -33.52 -1.25 -33.5 -1.27 -33.48 -1.29 -33.45 -1.31 -33.43 -1.33 -33.41 -1.35 -33.38 -1.37 -33.36 -1.39 -33.34 -1.41 -33.32 -1.43 -33.29 -1.45 -33.27 -1.47 -33.25 -1.5 -33.22 -1.52 -33.2 -1.54 -33.18 -1.56 -33.15 -1.58 -33.13 -1.6 -33.11 -1.62 -33.08 -1.64 -33.06 -1.66 -33.04 -1.67 -33.01 -1.69 -32.99 -1.71 -32.97 -1.73 -32.95 -1.75 -32.92 -1.77 -32.9 -1.79 -32.88 -1.81 -32.85 -1.83 -32.83 -1.85 -32.81 -1.87 -32.78 -1.88 -32.76 -1.9 -32.74 -1.92 -32.71 -1.94 -32.69 -1.96 -32.67 -1.97 -32.65 -1.99 -32.62 -2.01 -32.6 -2.03 -32.58 -2.05 -32.55 -2.06 -32.53 -2.08 -32.51 -2.1 -32.48 -2.11 -32.46 -2.13 -32.44 -2.15 -32.41 -2.16 -32.39 -2.18 -32.37 -2.2 -32.34 -2.21 -32.32 -2.23 -32.3 -2.24 -32.28 -2.26 -32.25 -2.28 -32.23 -2.29 -32.21 -2.31 -32.18 -2.32 -32.16 -2.34 -32.14 -2.35 -32.11 -2.37 -32.09 -2.38 -32.07 -2.4 -32.04 -2.41 -32.02 -2.42 -32 -2.44 -31.98 -2.45 -31.95 -2.47 -31.93 -2.48 -31.91 -2.49 -31.88 -2.51 -31.86 -2.52 -31.84 -2.53 -31.81 -2.55 -31.79 -2.56 -31.77 -2.57 -31.74 -2.58 -31.72 -2.6 -31.7 -2.61 -31.67 -2.62 -31.65 -2.63 -31.63 -2.64 -31.61 -2.66 -31.58 -2.67 -31.56 -2.68 -31.54 -2.69 -31.51 -2.7 -31.49 -2.71 -31.47 -2.72 -31.44 -2.73 -31.42 -2.74 -31.4 -2.76 -31.37 -2.77 -31.35 -2.78 -31.33 -2.79 -31.31 -2.8 -31.28 -2.81 -31.26 -2.81 -31.24 -2.82 -31.21 -2.83 -31.19 -2.84 -31.17 -2.85 -31.14 -2.86 -31.12 -2.87 -31.1 -2.88 -31.07 -2.89 -31.05 -2.89 -31.03 -2.9 -31 -2.91 -30.98 -2.92 -30.96 -2.93 -30.94 -2.93 -30.91 -2.94 -30.89 -2.95 -30.87 -2.96 -30.84 -2.96 -30.82 -2.97 -30.8 -2.98 -30.77 -2.98 -30.75 -2.99 -30.73 -3 -30.7 -3 -30.68 -3.01 -30.66 -3.02 -30.64 -3.02 -30.61 -3.03 -30.59 -3.03 -30.57 -3.04 -30.54 -3.05 -30.52 -3.05 -30.5 -3.06 -30.47 -3.06 -30.45 -3.07 -30.43 -3.07 -30.4 -3.08 -30.38 -3.08 -30.36 -3.09 -30.33 -3.09 -30.31 -3.09 -30.29 -3.1 -30.27 -3.1 -30.24 -3.11 -30.22 -3.11 -30.2 -3.11 -30.17 -3.12 -30.15 -3.12 -30.13 -3.13 -30.1 -3.13 -30.08 -3.13 -30.06 -3.14 -30.03 -3.14 -30.01 -3.14 -29.99 -3.14 -29.97 -3.15 -29.94 -3.15 -29.92 -3.15 -29.9 -3.15 -29.87 -3.16 -29.85 -3.16 -29.83 -3.16 -29.8 -3.16 -29.78 -3.16 -29.76 -3.17 -29.73 -3.17 -29.71 -3.17 -29.69 -3.17 -29.66 -3.17 -29.64 -3.17 -29.62 -3.18 -29.6 -3.18 -29.57 -3.18 -29.55 -3.18 -29.53 -3.18 -29.5 -3.18 -29.48 -3.18 -29.46 -3.18 -29.43 -3.18 -29.41 -3.18 -29.39 -3.18 -29.36 -3.18 -29.34 -3.18 -29.32 -3.18 -29.3 -3.18 -29.27 -3.18 -29.25 -3.18 -29.23 -3.18 -29.2 -3.18 -29.18 -3.18 -29.16 -3.18 -29.13 -3.18 -29.11 -3.18 -29.09 -3.18 -29.06 -3.18 -29.04 -3.18 -29.02 -3.18 -28.99 -3.18 -28.97 -3.18 -28.95 -3.17 -28.93 -3.17 -28.9 -3.17 -28.88 -3.17 -28.86 -3.17 -28.83 -3.17 -28.81 -3.17 -28.79 -3.16 -28.76 -3.16 -28.74 -3.16 -28.72 -3.16 -28.69 -3.16 -28.67 -3.16 -28.65 -3.15 -28.63 -3.15 -28.6 -3.15 -28.58 -3.15 -28.56 -3.14 -28.53 -3.14 -28.51 -3.14 -28.49 -3.14 -28.46 -3.13 -28.44 -3.13 -28.42 -3.13 -28.39 -3.13 -28.37 -3.12 -28.35 -3.12 -28.32 -3.12 -28.3 -3.12 -28.28 -3.11 -28.26 -3.11 -28.23 -3.11 -28.21 -3.1 -28.19 -3.1 -28.16 -3.1 -28.14 -3.09 -28.12 -3.09 -28.09 -3.09 -28.07 -3.08 -28.05 -3.08 -28.02 -3.08 -28 -3.07 -27.98 -3.07 -27.96 -3.07 -27.93 -3.06 -27.91 -3.06 -27.89 -3.05 -27.86 -3.05 -27.84 -3.05 -27.82 -3.04 -27.79 -3.04 -27.77 -3.03 -27.75 -3.03 -27.72 -3.03 -27.7 -3.02 -27.68 -3.02 -27.65 -3.01 -27.63 -3.01 -27.61 -3.01 -27.59 -3 -27.56 -3 -27.54 -2.99 -27.52 -2.99 -27.49 -2.98 -27.47 -2.98 -27.45 -2.97 -27.42 -2.97 -27.4 -2.96 -27.38 -2.96 -27.35 -2.96 -27.33 -2.95 -27.31 -2.95 -27.29 -2.94 -27.26 -2.94 -27.24 -2.93 -27.22 -2.93 -27.19 -2.92 -27.17 -2.92 -27.15 -2.91 -27.12 -2.91 -27.1 -2.9 -27.08 -2.9 -27.05 -2.89 -27.03 -2.89 -27.01 -2.88 -26.98 -2.88 -26.96 -2.87 -26.94 -2.86 -26.92 -2.86 -26.89 -2.85 -26.87 -2.85 -26.85 -2.84 -26.82 -2.84 -26.8 -2.83 -26.78 -2.83 -26.75 -2.82 -26.73 -2.82 -26.71 -2.81 -26.68 -2.8 -26.66 -2.8 -26.64 -2.79 -26.62 -2.79 -26.59 -2.78 -26.57 -2.78 -26.55 -2.77 -26.52 -2.77 -26.5 -2.76 -26.48 -2.75 -26.45 -2.75 -26.43 -2.74 -26.41 -2.74 -26.38 -2.73 -26.36 -2.72 -26.34 -2.72 -26.31 -2.71 -26.29 -2.71 -26.27 -2.7 -26.25 -2.7 -26.22 -2.69 -26.2 -2.68 -26.18 -2.68 -26.15 -2.67 -26.13 -2.67 -26.11 -2.66 -26.08 -2.65 -26.06 -2.65 -26.04 -2.64 -26.01 -2.64 -25.99 -2.63 -25.97 -2.62 -25.95 -2.62 -25.92 -2.61 -25.9 -2.6 -25.88 -2.6 -25.85 -2.59 -25.83 -2.59 -25.81 -2.58 -25.78 -2.57 -25.76 -2.57 -25.74 -2.56 -25.71 -2.56 -25.69 -2.55 -25.67 -2.54 -25.64 -2.54 -25.62 -2.53 -25.6 -2.52 -25.58 -2.52 -25.55 -2.51 -25.53 -2.51 -25.51 -2.5 -25.48 -2.49 -25.46 -2.49 -25.44 -2.48 -25.41 -2.47 -25.39 -2.47 -25.37 -2.46 -25.34 -2.45 -25.32 -2.45 -25.3 -2.44 -25.28 -2.44 -25.25 -2.43 -25.23 -2.42 -25.21 -2.42 -25.18 -2.41 -25.16 -2.4 -25.14 -2.4 -25.11 -2.39 -25.09 -2.38 -25.07 -2.38 -25.04 -2.37 -25.02 -2.36 -25 -2.36 -24.97 -2.35 -24.95 -2.35 -24.93 -2.34 -24.91 -2.33 -24.88 -2.33 -24.86 -2.32 -24.84 -2.31 -24.81 -2.31 -24.79 -2.3 -24.77 -2.29 -24.74 -2.29 -24.72 -2.28 -24.7 -2.27 -24.67 -2.27 -24.65 -2.26 -24.63 -2.25 -24.61 -2.25 -24.58 -2.24 -24.56 -2.24 -24.54 -2.23 -24.51 -2.22 -24.49 -2.22 -24.47 -2.21 -24.44 -2.2 -24.42 -2.2 -24.4 -2.19 -24.37 -2.18 -24.35 -2.18 -24.33 -2.17 -24.3 -2.16 -24.28 -2.16 -24.26 -2.15 -24.24 -2.14 -24.21 -2.14 -24.19 -2.13 -24.17 -2.13 -24.14 -2.12 -24.12 -2.11 -24.1 -2.11 -24.07 -2.1 -24.05 -2.09 -24.03 -2.09 -24 -2.08 -23.98 -2.07 -23.96 -2.07 -23.94 -2.06 -23.91 -2.05 -23.89 -2.05 -23.87 -2.04 -23.84 -2.04 -23.82 -2.03 -23.8 -2.02 -23.77 -2.02 -23.75 -2.01 -23.73 -2 -23.7 -2 -23.68 -1.99 -23.66 -1.98 -23.63 -1.98 -23.61 -1.97 -23.59 -1.96 -23.57 -1.96 -23.54 -1.95 -23.52 -1.95 -23.5 -1.94 -23.47 -1.93 -23.45 -1.93 -23.43 -1.92 -23.4 -1.91 -23.38 -1.91 -23.36 -1.9 -23.33 -1.9 -23.31 -1.89 -23.29 -1.88 -23.27 -1.88 -23.24 -1.87 -23.22 -1.86 -23.2 -1.86 -23.17 -1.85 -23.15 -1.84 -23.13 -1.84 -23.1 -1.83 -23.08 -1.83 -23.06 -1.82 -23.03 -1.81 -23.01 -1.81 -22.99 -1.8 -22.96 -1.79 -22.94 -1.79 -22.92 -1.78 -22.9 -1.78 -22.87 -1.77 -22.85 -1.76 -22.83 -1.76 -22.8 -1.75 -22.78 -1.75 -22.76 -1.74 -22.73 -1.73 -22.71 -1.73 -22.69 -1.72 -22.66 -1.71 -22.64 -1.71 -22.62 -1.7 -22.6 -1.7 -22.57 -1.69 -22.55 -1.68 -22.53 -1.68 -22.5 -1.67 -22.48 -1.67 -22.46 -1.66 -22.43 -1.65 -22.41 -1.65 -22.39 -1.64 -22.36 -1.64 -22.34 -1.63 -22.32 -1.62 -22.29 -1.62 -22.27 -1.61 -22.25 -1.61 -22.23 -1.6 -22.2 -1.59 -22.18 -1.59 -22.16 -1.58 -22.13 -1.58 -22.11 -1.57 -22.09 -1.57 -22.06 -1.56 -22.04 -1.55 -22.02 -1.55 -21.99 -1.54 -21.97 -1.54 -21.95 -1.53 -21.93 -1.52 -21.9 -1.52 -21.88 -1.51 -21.86 -1.51 -21.83 -1.5 -21.81 -1.5 -21.79 -1.49 -21.76 -1.48 -21.74 -1.48 -21.72 -1.47 -21.69 -1.47 -21.67 -1.46 -21.65 -1.45 -21.62 -1.45 -21.6 -1.44 -21.58 -1.44 -21.56 -1.43 -21.53 -1.43 -21.51 -1.42 -21.49 -1.42 -21.46 -1.41 -21.44 -1.4 -21.42 -1.4 -21.39 -1.39 -21.37 -1.39 -21.35 -1.38 -21.32 -1.38 -21.3 -1.37 -21.28 -1.37 -21.26 -1.36 -21.23 -1.35 -21.21 -1.35 -21.19 -1.34 -21.16 -1.34 -21.14 -1.33 -21.12 -1.33 -21.09 -1.32 -21.07 -1.32 -21.05 -1.31 -21.02 -1.3 -21 -1.3 -20.98 -1.29 -20.95 -1.29 -20.93 -1.28 -20.91 -1.28 -20.89 -1.27 -20.86 -1.27 -20.84 -1.26 -20.82 -1.26 -20.79 -1.25 -20.77 -1.25 -20.75 -1.24 -20.72 -1.24 -20.7 -1.23 -20.68 -1.22 -20.65 -1.22 -20.63 -1.21 -20.61 -1.21 -20.59 -1.2 -20.56 -1.2 -20.54 -1.19 -20.52 -1.19 -20.49 -1.18 -20.47 -1.18 -20.45 -1.17 -20.42 -1.17 -20.4 -1.16 -20.38 -1.16 -20.35 -1.15 -20.33 -1.15 -20.31 -1.14 -20.28 -1.14 -20.26 -1.13 -20.24 -1.13 -20.22 -1.12 -20.19 -1.12 -20.17 -1.11 -20.15 -1.11 -20.12 -1.1 -20.1 -1.1 -20.08 -1.09 -20.05 -1.09 -20.03 -1.08 -20.01 -1.08 -19.98 -1.07 -19.96 -1.07 -19.94 -1.06 -19.92 -1.06 -19.89 -1.05 -19.87 -1.05 -19.85 -1.04 -19.82 -1.04 -19.8 -1.03 -19.78 -1.03 -19.75 -1.02 -19.73 -1.02 -19.71 -1.01 -19.68 -1.01 -19.66 -1 -19.64 -1 -19.61 -0.99 -19.59 -0.99 -19.57 -0.99 -19.55 -0.98 -19.52 -0.98 -19.5 -0.97 -19.48 -0.97 -19.45 -0.96 -19.43 -0.96 -19.41 -0.95 -19.38 -0.95 -19.36 -0.94 -19.34 -0.94 -19.31 -0.93 -19.29 -0.93 -19.27 -0.92 -19.25 -0.92 -19.22 -0.92 -19.2 -0.91 -19.18 -0.91 -19.15 -0.9 -19.13 -0.9 -19.11 -0.89 -19.08 -0.89 -19.06 -0.88 -19.04 -0.88 -19.01 -0.88 -18.99 -0.87 -18.97 -0.87 -18.94 -0.86 -18.92 -0.86 -18.9 -0.85 -18.88 -0.85 -18.85 -0.84 -18.83 -0.84 -18.81 -0.84 -18.78 -0.83 -18.76 -0.83 -18.74 -0.82 -18.71 -0.82 -18.69 -0.81 -18.67 -0.81 -18.64 -0.81 -18.62 -0.8 -18.6 -0.8 -18.58 -0.79 -18.55 -0.79 -18.53 -0.78 -18.51 -0.78 -18.48 -0.78 -18.46 -0.77 -18.44 -0.77 -18.41 -0.76 -18.39 -0.76 -18.37 -0.75 -18.34 -0.75 -18.32 -0.75 -18.3 -0.74 -18.27 -0.74 -18.25 -0.73 -18.23 -0.73 -18.21 -0.73 -18.18 -0.72 -18.16 -0.72 -18.14 -0.71 -18.11 -0.71 -18.09 -0.71 -18.07 -0.7 -18.04 -0.7 -18.02 -0.69 -18 -0.69 -17.97 -0.69 -17.95 -0.68 -17.93 -0.68 -17.91 -0.67 -17.88 -0.67 -17.86 -0.67 -17.84 -0.66 -17.81 -0.66 -17.79 -0.65 -17.77 -0.65 -17.74 -0.65 -17.72 -0.64 -17.7 -0.64 -17.67 -0.64 -17.65 -0.63 -17.63 -0.63 -17.6 -0.62 -17.58 -0.62 -17.56 -0.62 -17.54 -0.61 -17.51 -0.61 -17.49 -0.6 -17.47 -0.6 -17.44 -0.6 -17.42 -0.59 -17.4 -0.59 -17.37 -0.59 -17.35 -0.58 -17.33 -0.58 -17.3 -0.58 -17.28 -0.57 -17.26 -0.57 -17.24 -0.56 -17.21 -0.56 -17.19 -0.56 -17.17 -0.55 -17.14 -0.55 -17.12 -0.55 -17.1 -0.54 -17.07 -0.54 -17.05 -0.54 -17.03 -0.53 -17 -0.53 -16.98 -0.52 -16.96 -0.52 -16.93 -0.52 -16.91 -0.51 -16.89 -0.51 -16.87 -0.51 -16.84 -0.5 -16.82 -0.5 -16.8 -0.5 -16.77 -0.49 -16.75 -0.49 -16.73 -0.49 -16.7 -0.48 -16.68 -0.48 -16.66 -0.48 -16.63 -0.47 -16.61 -0.47 -16.59 -0.47 -16.57 -0.46 -16.54 -0.46 -16.52 -0.46 -16.5 -0.45 -16.47 -0.45 -16.45 -0.45 -16.43 -0.44 -16.4 -0.44 -16.38 -0.44 -16.36 -0.43 -16.33 -0.43 -16.31 -0.43 -16.29 -0.42 -16.26 -0.42 -16.24 -0.42 -16.22 -0.41 -16.2 -0.41 -16.17 -0.41 -16.15 -0.4 -16.13 -0.4 -16.1 -0.4 -16.08 -0.39 -16.06 -0.39 -16.03 -0.39 -16.01 -0.38 -15.99 -0.38 -15.96 -0.38 -15.94 -0.38 -15.92 -0.37 -15.9 -0.37 -15.87 -0.37 -15.85 -0.36 -15.83 -0.36 -15.8 -0.36 -15.78 -0.35 -15.76 -0.35 -15.73 -0.35 -15.71 -0.34 -15.69 -0.34 -15.66 -0.34 -15.64 -0.34 -15.62 -0.33 -15.59 -0.33 -15.57 -0.33 -15.55 -0.32 -15.53 -0.32 -15.5 -0.32 -15.48 -0.32 -15.46 -0.31 -15.43 -0.31 -15.41 -0.31 -15.39 -0.3 -15.36 -0.3 -15.34 -0.3 -15.32 -0.29 -15.29 -0.29 -15.27 -0.29 -15.25 -0.29 -15.23 -0.28 -15.2 -0.28 -15.18 -0.28 -15.16 -0.28 -15.13 -0.27 -15.11 -0.27 -15.09 -0.27 -15.06 -0.26 -15.04 -0.26 -15.02 -0.26 -14.99 -0.26 -14.97 -0.25 -14.95 -0.25 -14.92 -0.25 -14.9 -0.24 -14.88 -0.24 -14.86 -0.24 -14.83 -0.24 -14.81 -0.23 -14.79 -0.23 -14.76 -0.23 -14.74 -0.23 -14.72 -0.22 -14.69 -0.22 -14.67 -0.22 -14.65 -0.22 -14.62 -0.21 -14.6 -0.21 -14.58 -0.21 -14.56 -0.2 -14.53 -0.2 -14.51 -0.2 -14.49 -0.2 -14.46 -0.19 -14.44 -0.19 -14.42 -0.19 -14.39 -0.19 -14.37 -0.18 -14.35 -0.18 -14.32 -0.18 -14.3 -0.18 -14.28 -0.17 -14.25 -0.17 -14.23 -0.17 -14.21 -0.17 -14.19 -0.16 -14.16 -0.16 -14.14 -0.16 -14.12 -0.16 -14.09 -0.15 -14.07 -0.15 -14.05 -0.15 -14.02 -0.15 -14 -0.14 -13.98 -0.14 -13.95 -0.14 -13.93 -0.14 -13.91 -0.13 -13.89 -0.13 -13.86 -0.13 -13.84 -0.13 -13.82 -0.13 -13.79 -0.12 -13.77 -0.12 -13.75 -0.12 -13.72 -0.12 -13.7 -0.11 -13.68 -0.11 -13.65 -0.11 -13.63 -0.11 -13.61 -0.1 -13.58 -0.1 -13.56 -0.1 -13.54 -0.1 -13.52 -0.1 -13.49 -0.09 -13.47 -0.09 -13.45 -0.09 -13.42 -0.09 -13.4 -0.08 -13.38 -0.08 -13.35 -0.08 -13.33 -0.08 -13.31 -0.08 -13.28 -0.07 -13.26 -0.07 -13.24 -0.07 -13.22 -0.07 -13.19 -0.06 -13.17 -0.06 -13.15 -0.06 -13.12 -0.06 -13.1 -0.06 -13.08 -0.05 -13.05 -0.05 -13.03 -0.05 -13.01 -0.05 -12.98 -0.04 -12.96 -0.04 -12.94 -0.04 -12.91 -0.04 -12.89 -0.04 -12.87 -0.03 -12.85 -0.03 -12.82 -0.03 -12.8 -0.03 -12.78 -0.03 -12.75 -0.02 -12.73 -0.02 -12.71 -0.02 -12.68 -0.02 -12.66 -0.02 -12.64 -0.01 -12.61 -0.01 -12.59 -0.01 -12.57 -0.01 -12.55 -0.01 -12.52 -0 -12.5 -0 -12.48 -0 -12.45 0 -12.43 0 -12.41 0.01 -12.38 0.01 -12.36 0.01 -12.34 0.01 -12.31 0.01 -12.29 0.02 -12.27 0.02 -12.24 0.02 -12.22 0.02 -12.2 0.02 -12.18 0.03 -12.15 0.03 -12.13 0.03 -12.11 0.03 -12.08 0.03 -12.06 0.03 -12.04 0.04 -12.01 0.04 -11.99 0.04 -11.97 0.04 -11.94 0.04 -11.92 0.05 -11.9 0.05 -11.88 0.05 -11.85 0.05 -11.83 0.05 -11.81 0.05 -11.78 0.06 -11.76 0.06 -11.74 0.06 -11.71 0.06 -11.69 0.06 -11.67 0.07 -11.64 0.07 -11.62 0.07 -11.6 0.07 -11.57 0.07 -11.55 0.07 -11.53 0.08 -11.51 0.08 -11.48 0.08 -11.46 0.08 -11.44 0.08 -11.41 0.08 -11.39 0.09 -11.37 0.09 -11.34 0.09 -11.32 0.09 -11.3 0.09 -11.27 0.09 -11.25 0.1 -11.23 0.1 -11.21 0.1 -11.18 0.1 -11.16 0.1 -11.14 0.1 -11.11 0.11 -11.09 0.11 -11.07 0.11 -11.04 0.11 -11.02 0.11 -11 0.11 -10.97 0.12 -10.95 0.12 -10.93 0.12 -10.9 0.12 -10.88 0.12 -10.86 0.12 -10.84 0.13 -10.81 0.13 -10.79 0.13 -10.77 0.13 -10.74 0.13 -10.72 0.13 -10.7 0.13 -10.67 0.14 -10.65 0.14 -10.63 0.14 -10.6 0.14 -10.58 0.14 -10.56 0.14 -10.54 0.15 -10.51 0.15 -10.49 0.15 -10.47 0.15 -10.44 0.15 -10.42 0.15 -10.4 0.15 -10.37 0.16 -10.35 0.16 -10.33 0.16 -10.3 0.16 -10.28 0.16 -10.26 0.16 -10.23 0.16 -10.21 0.17 -10.19 0.17 -10.17 0.17 -10.14 0.17 -10.12 0.17 -10.1 0.17 -10.07 0.17 -10.05 0.18 -10.03 0.18 -10 0.18 -9.98 0.18 -9.96 0.18 -9.93 0.18 -9.91 0.18 -9.89 0.19 -9.87 0.19 -9.84 0.19 -9.82 0.19 -9.8 0.19 -9.77 0.19 -9.75 0.19 -9.73 0.2 -9.7 0.2 -9.68 0.2 -9.66 0.2 -9.63 0.2 -9.61 0.2 -9.59 0.2 -9.56 0.2 -9.54 0.21 -9.52 0.21 -9.5 0.21 -9.47 0.21 -9.45 0.21 -9.43 0.21 -9.4 0.21 -9.38 0.21 -9.36 0.22 -9.33 0.22 -9.31 0.22 -9.29 0.22 -9.26 0.22 -9.24 0.22 -9.22 0.22 -9.2 0.22 -9.17 0.23 -9.15 0.23 -9.13 0.23 -9.1 0.23 -9.08 0.23 -9.06 0.23 -9.03 0.23 -9.01 0.23 -8.99 0.24 -8.96 0.24 -8.94 0.24 -8.92 0.24 -8.89 0.24 -8.87 0.24 -8.85 0.24 -8.83 0.24 -8.8 0.25 -8.78 0.25 -8.76 0.25 -8.73 0.25 -8.71 0.25 -8.69 0.25 -8.66 0.25 -8.64 0.25 -8.62 0.25 -8.59 0.26 -8.57 0.26 -8.55 0.26 -8.53 0.26 -8.5 0.26 -8.48 0.26 -8.46 0.26 -8.43 0.26 -8.41 0.26 -8.39 0.27 -8.36 0.27 -8.34 0.27 -8.32 0.27 -8.29 0.27 -8.27 0.27 -8.25 0.27 -8.22 0.27 -8.2 0.27 -8.18 0.28 -8.16 0.28 -8.13 0.28 -8.11 0.28 -8.09 0.28 -8.06 0.28 -8.04 0.28 -8.02 0.28 -7.99 0.28 -7.97 0.29 -7.95 0.29 -7.92 0.29 -7.9 0.29 -7.88 0.29 -7.86 0.29 -7.83 0.29 -7.81 0.29 -7.79 0.29 -7.76 0.29 -7.74 0.3 -7.72 0.3 -7.69 0.3 -7.67 0.3 -7.65 0.3 -7.62 0.3 -7.6 0.3 -7.58 0.3 -7.55 0.3 -7.53 0.3 -7.51 0.31 -7.49 0.31 -7.46 0.31 -7.44 0.31 -7.42 0.31 -7.39 0.31 -7.37 0.31 -7.35 0.31 -7.32 0.31 -7.3 0.31 -7.28 0.32 -7.25 0.32 -7.23 0.32 -7.21 0.32 -7.19 0.32 -7.16 0.32 -7.14 0.32 -7.12 0.32 -7.09 0.32 -7.07 0.32 -7.05 0.32 -7.02 0.33 -7 0.33 -6.98 0.33 -6.95 0.33 -6.93 0.33 -6.91 0.33 -6.88 0.33 -6.86 0.33 -6.84 0.33 -6.82 0.33 -6.79 0.33 -6.77 0.34 -6.75 0.34 -6.72 0.34 -6.7 0.34 -6.68 0.34 -6.65 0.34 -6.63 0.34 -6.61 0.34 -6.58 0.34 -6.56 0.34 -6.54 0.34 -6.52 0.34 -6.49 0.35 -6.47 0.35 -6.45 0.35 -6.42 0.35 -6.4 0.35 -6.38 0.35 -6.35 0.35 -6.33 0.35 -6.31 0.35 -6.28 0.35 -6.26 0.35 -6.24 0.35 -6.21 0.36 -6.19 0.36 -6.17 0.36 -6.15 0.36 -6.12 0.36 -6.1 0.36 -6.08 0.36 -6.05 0.36 -6.03 0.36 -6.01 0.36 -5.98 0.36 -5.96 0.36 -5.94 0.36 -5.91 0.37 -5.89 0.37 -5.87 0.37 -5.85 0.37 -5.82 0.37 -5.8 0.37 -5.78 0.37 -5.75 0.37 -5.73 0.37 -5.71 0.37 -5.68 0.37 -5.66 0.37 -5.64 0.37 -5.61 0.38 -5.59 0.38 -5.57 0.38 -5.54 0.38 -5.52 0.38 -5.5 0.38 -5.48 0.38 -5.45 0.38 -5.43 0.38 -5.41 0.38 -5.38 0.38 -5.36 0.38 -5.34 0.38 -5.31 0.38 -5.29 0.39 -5.27 0.39 -5.24 0.39 -5.22 0.39 -5.2 0.39 -5.18 0.39 -5.15 0.39 -5.13 0.39 -5.11 0.39 -5.08 0.39 -5.06 0.39 -5.04 0.39 -5.01 0.39 -4.99 0.39 -4.97 0.39 -4.94 0.4 -4.92 0.4 -4.9 0.4 -4.87 0.4 -4.85 0.4 -4.83 0.4 -4.81 0.4 -4.78 0.4 -4.76 0.4 -4.74 0.4 -4.71 0.4 -4.69 0.4 -4.67 0.4 -4.64 0.4 -4.62 0.4 -4.6 0.41 -4.57 0.41 -4.55 0.41 -4.53 0.41 -4.51 0.41 -4.48 0.41 -4.46 0.41 -4.44 0.41 -4.41 0.41 -4.39 0.41 -4.37 0.41 -4.34 0.41 -4.32 0.41 -4.3 0.41 -4.27 0.41 -4.25 0.41 -4.23 0.42 -4.2 0.42 -4.18 0.42 -4.16 0.42 -4.14 0.42 -4.11 0.42 -4.09 0.42 -4.07 0.42 -4.04 0.42 -4.02 0.42 -4 0.42 -3.97 0.42 -3.95 0.42 -3.93 0.42 -3.9 0.42 -3.88 0.42 -3.86 0.42 -3.84 0.43 -3.81 0.43 -3.79 0.43 -3.77 0.43 -3.74 0.43 -3.72 0.43 -3.7 0.43 -3.67 0.43 -3.65 0.43 -3.63 0.43 -3.6 0.43 -3.58 0.43 -3.56 0.43 -3.53 0.43 -3.51 0.43 -3.49 0.43 -3.47 0.43 -3.44 0.43 -3.42 0.43 -3.4 0.44 -3.37 0.44 -3.35 0.44 -3.33 0.44 -3.3 0.44 -3.28 0.44 -3.26 0.44 -3.23 0.44 -3.21 0.44 -3.19 0.44 -3.17 0.44 -3.14 0.44 -3.12 0.44 -3.1 0.44 -3.07 0.44 -3.05 0.44 -3.03 0.44 -3 0.44 -2.98 0.44 -2.96 0.44 -2.93 0.45 -2.91 0.45 -2.89 0.45 -2.86 0.45 -2.84 0.45 -2.82 0.45 -2.8 0.45 -2.77 0.45 -2.75 0.45 -2.73 0.45 -2.7 0.45 -2.68 0.45 -2.66 0.45 -2.63 0.45 -2.61 0.45 -2.59 0.45 -2.56 0.45 -2.54 0.45 -2.52 0.45 -2.5 0.45 -2.47 0.45 -2.45 0.46 -2.43 0.46 -2.4 0.46 -2.38 0.46 -2.36 0.46 -2.33 0.46 -2.31 0.46 -2.29 0.46 -2.26 0.46 -2.24 0.46 -2.22 0.46 -2.19 0.46 -2.17 0.46 -2.15 0.46 -2.13 0.46 -2.1 0.46 -2.08 0.46 -2.06 0.46 -2.03 0.46 -2.01 0.46 -1.99 0.46 -1.96 0.46 -1.94 0.46 -1.92 0.47 -1.89 0.47 -1.87 0.47 -1.85 0.47 -1.83 0.47 -1.8 0.47 -1.78 0.47 -1.76 0.47 -1.73 0.47 -1.71 0.47 -1.69 0.47 -1.66 0.47 -1.64 0.47 -1.62 0.47 -1.59 0.47 -1.57 0.47 -1.55 0.47 -1.52 0.47 -1.5 0.47 -1.48 0.47 -1.46 0.47 -1.43 0.47 -1.41 0.47 -1.39 0.47 -1.36 0.47 -1.34 0.48 -1.32 0.48 -1.29 0.48 -1.27 0.48 -1.25 0.48 -1.22 0.48 -1.2 0.48 -1.18 0.48 -1.16 0.48 -1.13 0.48 -1.11 0.48 -1.09 0.48 -1.06 0.48 -1.04 0.48 -1.02 0.48 -0.99 0.48 -0.97 0.48 -0.95 0.48 -0.92 0.48 -0.9 0.48 -0.88 0.48 -0.85 0.48 -0.83 0.48 -0.81 0.48 -0.79 0.48 -0.76 0.48 -0.74 0.48 -0.72 0.48 -0.69 0.49 -0.67 0.49 -0.65 0.49 -0.62 0.49 -0.6 0.49 -0.58 0.49 -0.55 0.49 -0.53 0.49 -0.51 0.49 -0.49 0.49 -0.46 0.49 -0.44 0.49 -0.42 0.49 -0.39 0.49 -0.37 0.49 -0.35 0.49 -0.32 0.49 -0.3 0.49 -0.28 0.49 -0.25 0.49 -0.23 0.49 -0.21 0.49 -0.18 0.49 -0.16 0.49 -0.14 0.49 -0.12 0.49 -0.09 0.49 -0.07 0.49 -0.05 0.49 -0.02 0.49 0 0.49 0.02 0.5 0.05 0.5 0.07 0.5 0.09 0.5 0.12 0.5 0.14 0.5 0.16 0.5 0.18 0.5 0.21 0.5 0.23 0.5 0.25 0.5 0.28 0.5 0.3 0.5 0.32 0.5 0.35 0.5 0.37 0.5 0.39 0.5 0.42 0.5 0.44 0.5 0.46 0.5 0.49 0.5 0.51 0.5 0.53 0.5 0.55 0.5 0.58 0.5 0.6 0.5 0.62 0.5 0.65 0.5 0.67 0.5 0.69 0.5 0.72 0.5 0.74 0.5 0.76 0.5 0.79 0.5 0.81 0.5 0.83 0.51 0.85 0.51 0.88 0.51 0.9 0.51 0.92 0.51 0.95 0.51 0.97 0.51 0.99 0.51 1.02 0.51 1.04 0.51 1.06 0.51 1.09 0.51 1.11 0.51 1.13 0.51 1.16 0.51 1.18 0.51 1.2 0.51 1.22 0.51 1.25 0.51 1.27 0.51 1.29 0.51 1.32 0.51 1.34 0.51 1.36 0.51 1.39 0.51 1.41 0.51 1.43 0.51 1.46 0.51 1.48 0.51 1.5 0.51 1.52 0.51 1.55 0.51 1.57 0.51 1.59 0.51 1.62 0.51 1.64 0.51 1.66 0.51 1.69 0.51 1.71 0.51 1.73 0.51 1.76 0.51 1.78 0.52 1.8 0.52 1.83 0.52 1.85 0.52 1.87 0.52 1.89 0.52 1.92 0.52 1.94 0.52 1.96 0.52 1.99 0.52 2.01 0.52 2.03 0.52 2.06 0.52 2.08 0.52 2.1 0.52 2.13 0.52 2.15 0.52 2.17 0.52 2.19 0.52 2.22 0.52 2.24 0.52 2.26 0.52 2.29 0.52 2.31 0.52 2.33 0.52 2.36 0.52 2.38 0.52 2.4 0.52 2.43 0.52 2.45 0.52 2.47 0.52 2.5 0.52 2.52 0.52 2.54 0.52 2.56 0.52 2.59 0.52 2.61 0.52 2.63 0.52 2.66 0.52 2.68 0.52 2.7 0.52 2.73 0.52 2.75 0.52 2.77 0.52 2.8 0.52 2.82 0.52 2.84 0.52 2.86 0.52 2.89 0.52 2.91 0.53 2.93 0.53 2.96 0.53 2.98 0.53 3 0.53 3.03 0.53 3.05 0.53 3.07 0.53 3.1 0.53 3.12 0.53 3.14 0.53 3.17 0.53 3.19 0.53 3.21 0.53 3.23 0.53 3.26 0.53 3.28 0.53 3.3 0.53 3.33 0.53 3.35 0.53 3.37 0.53 3.4 0.53 3.42 0.53 3.44 0.53 3.47 0.53 3.49 0.53 3.51 0.53 3.53 0.53 3.56 0.53 3.58 0.53 3.6 0.53 3.63 0.53 3.65 0.53 3.67 0.53 3.7 0.53 3.72 0.53 3.74 0.53 3.77 0.53 3.79 0.53 3.81 0.53 3.84 0.53 3.86 0.53 3.88 0.53 3.9 0.53 3.93 0.53 3.95 0.53 3.97 0.53 4 0.53 4.02 0.53 4.04 0.53 4.07 0.53 4.09 0.53 4.11 0.53 4.14 0.53 4.16 0.53 4.18 0.53 4.2 0.53 4.23 0.53 4.25 0.53 4.27 0.53 4.3 0.53 4.32 0.53 4.34 0.54 4.37 0.54 4.39 0.54 4.41 0.54 4.44 0.54 4.46 0.54 4.48 0.54 4.51 0.54 4.53 0.54 4.55 0.54 4.57 0.54 4.6 0.54 4.62 0.54 4.64 0.54 4.67 0.54 4.69 0.54 4.71 0.54 4.74 0.54 4.76 0.54 4.78 0.54 4.81 0.54 4.83 0.54 4.85 0.54 4.87 0.54 4.9 0.54 4.92 0.54 4.94 0.54 4.97 0.54 4.99 0.54 5.01 0.54 5.04 0.54 5.06 0.54 5.08 0.54 5.11 0.54 5.13 0.54 5.15 0.54 5.18 0.54 5.2 0.54 5.22 0.54 5.24 0.54 5.27 0.54 5.29 0.54 5.31 0.54 5.34 0.54 5.36 0.54 5.38 0.54 5.41 0.54 5.43 0.54 5.45 0.54 5.48 0.54 5.5 0.54 5.52 0.54 5.54 0.54 5.57 0.54 5.59 0.54 5.61 0.54 5.64 0.54 5.66 0.54 5.68 0.54 5.71 0.54 5.73 0.54 5.75 0.54 5.78 0.54 5.8 0.54 5.82 0.54 5.85 0.54 5.87 0.54 5.89 0.54 5.91 0.54 5.94 0.54 5.96 0.54 5.98 0.54 6.01 0.54 6.03 0.54 6.05 0.54 6.08 0.54 6.1 0.54 6.12 0.54 6.15 0.54 6.17 0.54 6.19 0.54 6.21 0.54 6.24 0.54 6.26 0.54 6.28 0.55 6.31 0.55 6.33 0.55 6.35 0.55 6.38 0.55 6.4 0.55 6.42 0.55 6.45 0.55 6.47 0.55 6.49 0.55 6.52 0.55 6.54 0.55 6.56 0.55 6.58 0.55 6.61 0.55 6.63 0.55 6.65 0.55 6.68 0.55 6.7 0.55 6.72 0.55 6.75 0.55 6.77 0.55 6.79 0.55 6.82 0.55 6.84 0.55 6.86 0.55 6.88 0.55 6.91 0.55 6.93 0.55 6.95 0.55 6.98 0.55 7 0.55 7.02 0.55 7.05 0.55 7.07 0.55 7.09 0.55 7.12 0.55 7.14 0.55 7.16 0.55 7.19 0.55 7.21 0.55 7.23 0.55 7.25 0.55 7.28 0.55 7.3 0.55 7.32 0.55 7.35 0.55 7.37 0.55 7.39 0.55 7.42 0.55 7.44 0.55 7.46 0.55 7.49 0.55 7.51 0.55 7.53 0.55 7.55 0.55 7.58 0.55 7.6 0.55 7.62 0.55 7.65 0.55 7.67 0.55 7.69 0.55 7.72 0.55 7.74 0.55 7.76 0.55 7.79 0.55 7.81 0.55 7.83 0.55 7.86 0.55 7.88 0.55 7.9 0.55 7.92 0.55 7.95 0.55 7.97 0.55 7.99 0.55 8.02 0.55 8.04 0.55 8.06 0.55 8.09 0.55 8.11 0.55 8.13 0.55 8.16 0.55 8.18 0.55 8.2 0.55 8.22 0.55 8.25 0.55 8.27 0.55 8.29 0.55 8.32 0.55 8.34 0.55 8.36 0.55 8.39 0.55 8.41 0.55 8.43 0.55 8.46 0.55 8.48 0.55 8.5 0.55 8.53 0.55 8.55 0.55 8.57 0.55 8.59 0.55 8.62 0.55 8.64 0.55 8.66 0.55 8.69 0.55 8.71 0.55 8.73 0.55 8.76 0.55 8.78 0.55 8.8 0.55 8.83 0.55 8.85 0.55 8.87 0.55 8.89 0.55 8.92 0.55 8.94 0.55 8.96 0.55 8.99 0.55 9.01 0.55 9.03 0.55 9.06 0.55 9.08 0.55 9.1 0.55 9.13 0.55 9.15 0.55 9.17 0.55 9.2 0.55 9.22 0.55 9.24 0.55 9.26 0.55 9.29 0.55 9.31 0.55 9.33 0.55 9.36 0.55 9.38 0.55 9.4 0.55 9.43 0.55 9.45 0.55 9.47 0.55 9.5 0.55 9.52 0.55 9.54 0.55 9.56 0.56 9.59 0.56 9.61 0.56 9.63 0.56 9.66 0.56 9.68 0.56 9.7 0.56 9.73 0.56 9.75 0.56 9.77 0.56 9.8 0.56 9.82 0.56 9.84 0.56 9.87 0.56 9.89 0.56 9.91 0.56 9.93 0.56 9.96 0.56 9.98 0.56 10 0.56 10.03 0.56 10.05 0.56 10.07 0.56 10.1 0.56 10.12 0.56 10.14 0.56 10.17 0.56 10.19 0.56 10.21 0.56 10.23 0.56 10.26 0.56 10.28 0.56 10.3 0.56 10.33 0.56 10.35 0.56 10.37 0.56 10.4 0.56 10.42 0.56 10.44 0.56 10.47 0.56 10.49 0.56 10.51 0.56 10.54 0.56 10.56 0.56 10.58 0.56 10.6 0.56 10.63 0.56 10.65 0.56 10.67 0.56 10.7 0.56 10.72 0.56 10.74 0.56 10.77 0.56 10.79 0.56 10.81 0.56 10.84 0.56 10.86 0.56 10.88 0.56 10.9 0.56 10.93 0.56 10.95 0.56 10.97 0.56 11 0.56 11.02 0.56 11.04 0.56 11.07 0.56 11.09 0.56 11.11 0.56 11.14 0.56 11.16 0.56 11.18 0.56 11.21 0.56 11.23 0.56 11.25 0.56 11.27 0.56 11.3 0.56 11.32 0.56 11.34 0.56 11.37 0.56 11.39 0.56 11.41 0.56 11.44 0.56 11.46 0.56 11.48 0.56 11.51 0.56 11.53 0.56 11.55 0.56 11.57 0.56 11.6 0.56 11.62 0.56 11.64 0.56 11.67 0.56 11.69 0.56 11.71 0.56 11.74 0.56 11.76 0.56 11.78 0.56 11.81 0.56 11.83 0.56 11.85 0.56 11.88 0.56 11.9 0.56 11.92 0.56 11.94 0.56 11.97 0.56 11.99 0.56 12.01 0.56 12.04 0.56 12.06 0.56 12.08 0.56 12.11 0.56 12.13 0.56 12.15 0.56 12.18 0.56 12.2 0.56 12.22 0.56 12.24 0.56 12.27 0.56 12.29 0.56 12.31 0.56 12.34 0.56 12.36 0.56 12.38 0.56 12.41 0.56 12.43 0.56 12.45 0.56 12.48 0.56 12.5 0.56 12.52 0.56 12.55 0.56 12.57 0.56 12.59 0.56 12.61 0.56 12.64 0.56 12.66 0.56 12.68 0.56 12.71 0.56 12.73 0.56 12.75 0.56 12.78 0.56 12.8 0.56 12.82 0.56 12.85 0.56 12.87 0.56 12.89 0.56 12.91 0.56 12.94 0.56 12.96 0.56 12.98 0.56 13.01 0.56 13.03 0.56 13.05 0.56 13.08 0.56 13.1 0.56 13.12 0.56 13.15 0.56 13.17 0.56 13.19 0.56 13.22 0.56 13.24 0.56 13.26 0.56 13.28 0.56 13.31 0.56 13.33 0.56 13.35 0.56 13.38 0.56 13.4 0.56 13.42 0.56 13.45 0.56 13.47 0.56 13.49 0.56 13.52 0.56 13.54 0.56 13.56 0.56 13.58 0.56 13.61 0.56 13.63 0.56 13.65 0.56 13.68 0.56 13.7 0.56 13.72 0.56 13.75 0.56 13.77 0.56 13.79 0.56 13.82 0.56 13.84 0.56 13.86 0.56 13.89 0.56 13.91 0.56 13.93 0.56 13.95 0.56 13.98 0.56 14 0.56 14.02 0.56 14.05 0.56 14.07 0.56 14.09 0.56 14.12 0.56 14.14 0.56 14.16 0.56 14.19 0.56 14.21 0.56 14.23 0.56 14.25 0.56 14.28 0.56 14.3 0.56 14.32 0.56 14.35 0.56 14.37 0.56 14.39 0.56 14.42 0.56 14.44 0.56 14.46 0.56 14.49 0.56 14.51 0.56 14.53 0.56 14.56 0.56 14.58 0.56 14.6 0.56 14.62 0.56 14.65 0.56 14.67 0.56 14.69 0.56 14.72 0.56 14.74 0.56 14.76 0.56 14.79 0.56 14.81 0.56 14.83 0.56 14.86 0.56 14.88 0.56 14.9 0.56 14.92 0.56 14.95 0.56 14.97 0.56 14.99 0.56 15.02 0.56 15.04 0.56 15.06 0.56 15.09 0.56 15.11 0.56 15.13 0.56 15.16 0.56 15.18 0.56 15.2 0.56 15.23 0.56 15.25 0.56 15.27 0.56 15.29 0.56 15.32 0.56 15.34 0.56 15.36 0.56 15.39 0.56 15.41 0.56 15.43 0.56 15.46 0.56 15.48 0.56 15.5 0.56 15.53 0.56 15.55 0.56 15.57 0.56 15.59 0.56 15.62 0.56 15.64 0.56 15.66 0.56 15.69 0.56 15.71 0.56 15.73 0.56 15.76 0.56 15.78 0.56 15.8 0.56 15.83 0.56 15.85 0.56 15.87 0.56 15.9 0.56 15.92 0.56 15.94 0.56 15.96 0.56 15.99 0.56 16.01 0.56 16.03 0.56 16.06 0.56 16.08 0.56 16.1 0.56 16.13 0.56 16.15 0.56 16.17 0.56 16.2 0.56 16.22 0.56 16.24 0.56 16.26 0.56 16.29 0.56 16.31 0.56 16.33 0.56 16.36 0.56 16.38 0.56 16.4 0.56 16.43 0.56 16.45 0.56 16.47 0.56 16.5 0.56 16.52 0.56 16.54 0.56 16.57 0.56 16.59 0.56 16.61 0.56 16.63 0.56 16.66 0.56 16.68 0.56 16.7 0.56 16.73 0.56 16.75 0.56 16.77 0.56 16.8 0.56 16.82 0.56 16.84 0.56 16.87 0.56 16.89 0.56 16.91 0.56 16.93 0.56 16.96 0.56 16.98 0.56 17 0.56 17.03 0.56 17.05 0.56 17.07 0.56 17.1 0.56 17.12 0.56 17.14 0.56 17.17 0.56 17.19 0.56 17.21 0.56 17.24 0.56 17.26 0.56 17.28 0.56 17.3 0.56 17.33 0.56 17.35 0.56 17.37 0.56 17.4 0.56 17.42 0.56 17.44 0.56 17.47 0.56 17.49 0.56 17.51 0.56 17.54 0.56 17.56 0.56 17.58 0.56 17.6 0.56 17.63 0.56 17.65 0.56 17.67 0.56 17.7 0.56 17.72 0.56 17.74 0.56 17.77 0.56 17.79 0.56 17.81 0.56 17.84 0.56 17.86 0.56 17.88 0.56 17.91 0.56 17.93 0.56 17.95 0.56 17.97 0.56 18 0.56 18.02 0.56 18.04 0.56 18.07 0.56 18.09 0.56 18.11 0.56 18.14 0.56 18.16 0.56 18.18 0.56 18.21 0.56 18.23 0.56 18.25 0.56 18.27 0.56 18.3 0.56 18.32 0.56 18.34 0.56 18.37 0.56 18.39 0.56 18.41 0.56 18.44 0.56 18.46 0.56 18.48 0.56 18.51 0.56 18.53 0.56 18.55 0.56 18.58 0.56 18.6 0.56 18.62 0.56 18.64 0.56 18.67 0.56 18.69 0.56 18.71 0.56 18.74 0.56 18.76 0.56 18.78 0.56 18.81 0.56 18.83 0.56 18.85 0.56 18.88 0.56 18.9 0.56 18.92 0.56 18.94 0.56 18.97 0.56 18.99 0.56 19.01 0.56 19.04 0.56 19.06 0.56 19.08 0.56 19.11 0.56 19.13 0.56 19.15 0.56 19.18 0.56 19.2 0.56 19.22 0.56 19.25 0.56 19.27 0.56 19.29 0.56 19.31 0.56 19.34 0.56 19.36 0.56 19.38 0.56 19.41 0.56 19.43 0.56 19.45 0.56 19.48 0.56 19.5 0.56 19.52 0.56 19.55 0.56 19.57 0.56 19.59 0.56 19.61 0.56 19.64 0.56 19.66 0.56 19.68 0.56 19.71 0.56 19.73 0.56 19.75 0.56 19.78 0.56 19.8 0.56 19.82 0.56 19.85 0.56 19.87 0.56 19.89 0.56 19.92 0.56 19.94 0.56 19.96 0.56 19.98 0.56 20.01 0.56 20.03 0.56 20.05 0.56 20.08 0.56 20.1 0.56 20.12 0.56 20.15 0.56 20.17 0.56 20.19 0.56 20.22 0.56 20.24 0.56 20.26 0.56 20.28 0.56 20.31 0.56 20.33 0.56 20.35 0.56 20.38 0.56 20.4 0.56 20.42 0.56 20.45 0.56 20.47 0.56 20.49 0.56 20.52 0.56 20.54 0.56 20.56 0.56 20.59 0.56 20.61 0.56 20.63 0.56 20.65 0.56 20.68 0.56 20.7 0.56 20.72 0.56 20.75 0.56 20.77 0.56 20.79 0.56 20.82 0.56 20.84 0.56 20.86 0.56 20.89 0.56 20.91 0.56 20.93 0.56 20.95 0.56 20.98 0.56 21 0.56 21.02 0.56 21.05 0.56 21.07 0.56 21.09 0.56 21.12 0.56 21.14 0.56 21.16 0.56 21.19 0.56 21.21 0.56 21.23 0.56 21.26 0.56 21.28 0.56 21.3 0.56 21.32 0.56 21.35 0.56 21.37 0.56 21.39 0.56 21.42 0.56 21.44 0.56 21.46 0.56 21.49 0.56 21.51 0.56 21.53 0.56 21.56 0.56 21.58 0.56 21.6 0.56 21.62 0.56 21.65 0.56 21.67 0.56 21.69 0.56 21.72 0.56 21.74 0.56 21.76 0.56 21.79 0.56 21.81 0.56 21.83 0.56 21.86 0.56 21.88 0.56 21.9 0.56 21.93 0.56 21.95 0.56 21.97 0.56 21.99 0.56 22.02 0.56 22.04 0.56 22.06 0.56 22.09 0.56 22.11 0.56 22.13 0.56 22.16 0.56 22.18 0.56 22.2 0.56 22.23 0.56 22.25 0.56 22.27 0.56 22.29 0.56 22.32 0.56 22.34 0.56 22.36 0.56 22.39 0.56 22.41 0.56 22.43 0.56 22.46 0.56 22.48 0.56 22.5 0.56 22.53 0.56 22.55 0.56 22.57 0.56 22.6 0.56 22.62 0.56 22.64 0.56 22.66 0.56 22.69 0.56 22.71 0.56 22.73 0.56 22.76 0.56 22.78 0.56 22.8 0.56 22.83 0.56 22.85 0.56 22.87 0.56 22.9 0.56 22.92 0.56 22.94 0.56 22.96 0.56 22.99 0.56 23.01 0.56 23.03 0.56 23.06 0.56 23.08 0.56 23.1 0.56 23.13 0.56 23.15 0.56 23.17 0.56 23.2 0.56 23.22 0.56 23.24 0.56 23.27 0.56 23.29 0.56 23.31 0.56 23.33 0.56 23.36 0.56 23.38 0.56 23.4 0.56 23.43 0.56 23.45 0.56 23.47 0.56 23.5 0.56 23.52 0.56 23.54 0.56 23.57 0.56 23.59 0.56 23.61 0.56 23.63 0.56 23.66 0.56 23.68 0.56 23.7 0.56 23.73 0.56 23.75 0.56 23.77 0.56 23.8 0.56 23.82 0.56 23.84 0.56 23.87 0.56 23.89 0.56 23.91 0.56 23.94 0.56 23.96 0.56 23.98 0.56 24 0.56 24.03 0.56 24.05 0.56 24.07 0.56 24.1 0.56 24.12 0.56 24.14 0.56 24.17 0.56 24.19 0.56 24.21 0.56 24.24 0.56 24.26 0.56 24.28 0.56 24.3 0.56 24.33 0.56 24.35 0.56 24.37 0.56 24.4 0.56 24.42 0.56 24.44 0.56 24.47 0.56 24.49 0.56 24.51 0.56 24.54 0.56 24.56 0.56 24.58 0.56 24.61 0.56 24.63 0.56 24.65 0.56 24.67 0.56 24.7 0.56 24.72 0.56 24.74 0.56 24.77 0.56 24.79 0.56 24.81 0.56 24.84 0.56 24.86 0.56 24.88 0.56 24.91 0.56 24.93 0.56 24.95 0.56 24.97 0.56 25 0.56 25.02 0.56 25.04 0.56 25.07 0.56 25.09 0.56 25.11 0.56 25.14 0.56 25.16 0.56 25.18 0.56 25.21 0.56 25.23 0.56 25.25 0.56 25.28 0.56 25.3 0.56 25.32 0.56 25.34 0.56 25.37 0.56 25.39 0.56 25.41 0.56 25.44 0.56 25.46 0.56 25.48 0.56 25.51 0.56 25.53 0.56 25.55 0.56 25.58 0.56 25.6 0.56 25.62 0.56 25.64 0.56 25.67 0.56 25.69 0.56 25.71 0.56 25.74 0.56 25.76 0.56 25.78 0.56 25.81 0.56 25.83 0.56 25.85 0.56 25.88 0.56 25.9 0.56 25.92 0.56 25.95 0.56 25.97 0.56 25.99 0.56 26.01 0.56 26.04 0.56 26.06 0.56 26.08 0.56 26.11 0.56 26.13 0.56 26.15 0.56 26.18 0.56 26.2 0.56 26.22 0.56 26.25 0.56 26.27 0.56 26.29 0.56 26.31 0.56 26.34 0.56 26.36 0.56 26.38 0.56 26.41 0.56 26.43 0.56 26.45 0.56 26.48 0.56 26.5 0.56 26.52 0.56 26.55 0.56 26.57 0.56 26.59 0.56 26.62 0.56 26.64 0.56 26.66 0.56 26.68 0.56 26.71 0.56 26.73 0.56 26.75 0.56 26.78 0.56 26.8 0.56 26.82 0.56 26.85 0.56 26.87 0.56 26.89 0.56 26.92 0.56 26.94 0.56 26.96 0.56 26.98 0.56 27.01 0.56 27.03 0.56 27.05 0.56 27.08 0.56 27.1 0.56 27.12 0.56 27.15 0.56 27.17 0.56 27.19 0.56 27.22 0.56 27.24 0.56 27.26 0.56 27.29 0.56 27.31 0.56 27.33 0.56 27.35 0.56 27.38 0.56 27.4 0.56 27.42 0.56 27.45 0.56 27.47 0.56 27.49 0.56 27.52 0.56 27.54 0.56 27.56 0.56 27.59 0.56 27.61 0.56 27.63 0.56 27.65 0.56 27.68 0.56 27.7 0.56 27.72 0.56 27.75 0.56 27.77 0.56 27.79 0.56 27.82 0.56 27.84 0.56 27.86 0.56 27.89 0.56 27.91 0.56 27.93 0.56 27.96 0.56 27.98 0.56 28 0.56 28.02 0.56 28.05 0.56 28.07 0.56 28.09 0.56 28.12 0.56 28.14 0.56 28.16 0.56 28.19 0.56 28.21 0.56 28.23 0.56 28.26 0.56 28.28 0.56 28.3 0.56 28.32 0.56 28.35 0.56 28.37 0.56 28.39 0.56 28.42 0.56 28.44 0.56 28.46 0.56 28.49 0.56 28.51 0.56 28.53 0.56 28.56 0.56 28.58 0.56 28.6 0.56 28.63 0.56 28.65 0.56 28.67 0.56 28.69 0.56 28.72 0.56 28.74 0.56 28.76 0.56 28.79 0.56 28.81 0.56 28.83 0.56 28.86 0.56 28.88 0.56 28.9 0.56 28.93 0.56 28.95 0.56 28.97 0.56 28.99 0.56 29.02 0.56 29.04 0.56 29.06 0.56 29.09 0.56 29.11 0.56 29.13 0.56 29.16 0.56 29.18 0.56 29.2 0.56 29.23 0.56 29.25 0.56 29.27 0.56 29.3 0.56 29.32 0.56 29.34 0.56 29.36 0.56 29.39 0.56 29.41 0.56 29.43 0.56 29.46 0.56 29.48 0.56 29.5 0.56 29.53 0.56 29.55 0.56 29.57 0.56 29.6 0.56 29.62 0.56 29.64 0.56 29.66 0.56 29.69 0.56 29.71 0.56 29.73 0.56 29.76 0.56 29.78 0.56 29.8 0.56 29.83 0.56 29.85 0.56 29.87 0.56 29.9 0.56 29.92 0.56 29.94 0.56 29.97 0.56 29.99 0.56 30.01 0.56 30.03 0.56 30.06 0.56 30.08 0.56 30.1 0.56 30.13 0.56 30.15 0.56 30.17 0.56 30.2 0.56 30.22 0.56 30.24 0.56 30.27 0.56 30.29 0.56 30.31 0.56 30.33 0.56 30.36 0.56 30.38 0.56 30.4 0.56 30.43 0.56 30.45 0.56 30.47 0.56 30.5 0.56 30.52 0.56 30.54 0.56 30.57 0.56 30.59 0.56 30.61 0.56 30.64 0.56 30.66 0.56 30.68 0.56 30.7 0.56 30.73 0.56 30.75 0.56 30.77 0.56 30.8 0.56 30.82 0.56 30.84 0.56 30.87 0.56 30.89 0.56 30.91 0.56 30.94 0.56 30.96 0.56 30.98 0.56 31 0.56 31.03 0.56 31.05 0.56 31.07 0.56 31.1 0.56 31.12 0.56 31.14 0.56 31.17 0.56 31.19 0.56 31.21 0.56 31.24 0.56 31.26 0.56 31.28 0.56 31.31 0.56 31.33 0.56 31.35 0.56 31.37 0.56 31.4 0.56 31.42 0.56 31.44 0.56 31.47 0.56 31.49 0.56 31.51 0.56 31.54 0.56 31.56 0.56 31.58 0.56 31.61 0.56 31.63 0.56 31.65 0.56 31.67 0.56 31.7 0.56 31.72 0.56 31.74 0.56 31.77 0.56 31.79 0.56 31.81 0.56 31.84 0.56 31.86 0.56 31.88 0.56 31.91 0.56 31.93 0.56 31.95 0.56 31.98 0.56 32 0.56 32.02 0.56 32.04 0.56 32.07 0.56 32.09 0.56 32.11 0.56 32.14 0.56 32.16 0.56 32.18 0.56 32.21 0.56 32.23 0.56 32.25 0.56 32.28 0.56 32.3 0.56 32.32 0.56 32.34 0.56 32.37 0.56 32.39 0.56 32.41 0.56 32.44 0.56 32.46 0.56 32.48 0.56 32.51 0.56 32.53 0.56 32.55 0.56 32.58 0.56 32.6 0.56 32.62 0.56 32.65 0.56 32.67 0.56 32.69 0.56 32.71 0.56 32.74 0.56 32.76 0.56 32.78 0.56 32.81 0.56 32.83 0.56 32.85 0.56 32.88 0.56 32.9 0.56 32.92 0.56 32.95 0.56 32.97 0.56 32.99 0.56 33.01 0.56 33.04 0.56 33.06 0.56 33.08 0.56 33.11 0.56 33.13 0.56 33.15 0.56 33.18 0.56 33.2 0.56 33.22 0.56 33.25 0.56 33.27 0.56 33.29 0.56 33.32 0.56 33.34 0.56 33.36 0.56 33.38 0.56 33.41 0.56 33.43 0.56 33.45 0.56 33.48 0.56 33.5 0.56 33.52 0.56 33.55 0.56 33.57 0.56 33.59 0.56 33.62 0.56 33.64 0.56 33.66 0.56 33.68 0.56 33.71 0.56 33.73 0.56 33.75 0.56 33.78 0.56 33.8 0.56 33.82 0.56 33.85 0.56 33.87 0.56 33.89 0.56 33.92 0.56 33.94 0.56 33.96 0.56 33.99 0.56 34.01 0.56 34.03 0.56 34.05 0.56 34.08 0.56 34.1 0.56 34.12 0.56 34.15 0.56 34.17 0.56 34.19 0.56 34.22 0.56 34.24 0.56 34.26 0.56 34.29 0.56 34.31 0.56 34.33 0.56 34.35 0.56 34.38 0.56 34.4 0.56 34.42 0.56 34.45 0.56 34.47 0.56 34.49 0.56 34.52 0.56 34.54 0.56 34.56 0.56 34.59 0.56 34.61 0.56 34.63 0.56 34.66 0.56 34.68 0.56 34.7 0.56 34.72 0.56 34.75 0.56 34.77 0.56 34.79 0.56 34.82 0.56 34.84 0.56 34.86 0.56 34.89 0.56 34.91 0.56 34.93 0.56 34.96 0.56 34.98 0.56 35 0.56 35.02 0.56 35.05 0.56 35.07 0.56 35.09 0.56 35.12 0.56 35.14 0.56 35.16 0.56 35.19 0.56 35.21 0.56 35.23 0.56 35.26 0.56 35.28 0.56 35.3 0.56 35.33 0.56 35.35 0.56 35.37 0.56 35.39 0.56 35.42 0.56 35.44 0.56 35.46 0.56 35.49 0.56 35.51 0.56 35.53 0.56 35.56 0.56 35.58 0.56 35.6 0.56 35.63 0.56 35.65 0.56 35.67 0.56 35.69 0.56 35.72 0.56 35.74 0.56 35.76 0.56 35.79 0.56 35.81 0.56 35.83 0.56 35.86 0.56 35.88 0.56 35.9 0.56 35.93 0.56 35.95 0.56 35.97 0.56 36 0.56 36.02 0.56 36.04 0.56 36.06 0.56 36.09 0.56 36.11 0.56 36.13 0.56 36.16 0.56 36.18 0.56 36.2 0.56 36.23 0.56 36.25 0.56 36.27 0.56 36.3 0.56 36.32 0.56 36.34 0.56 36.36 0.56 36.39 0.56 36.41 0.56 36.43 0.56 36.46 0.56 36.48 0.56 36.5 0.56 36.53 0.56 36.55 0.56 36.57 0.56 36.6 0.56 36.62 0.56 36.64 0.56 36.67 0.56 36.69 0.56 36.71 0.56 36.73 0.56 36.76 0.56 36.78 0.56 36.8 0.56 36.83 0.56 36.85 0.56 36.87 0.56 36.9 0.56 36.92 0.56 36.94 0.56 36.97 0.56 36.99 0.56 37.01 0.56 37.03 0.56 37.06 0.56 37.08 0.56 37.1 0.56 37.13 0.56 37.15 0.56 37.17 0.56 37.2 0.56 37.22 0.56 37.24 0.56 37.27 0.56 37.29 0.56 37.31 0.56 37.34 0.56 37.36 0.56 37.38 0.56 37.4 0.56 37.43 0.56 37.45 0.56 37.47 0.56 37.5 0.56 37.52 0.56 37.54 0.56 37.57 0.56 37.59 0.56 37.61 0.56 37.64 0.56 37.66 0.56 37.68 0.56 37.7 0.56 37.73 0.56 37.75 0.56 37.77 0.56 37.8 0.56 37.82 0.56 37.84 0.56 37.87 0.56 37.89 0.56 37.91 0.56 37.94 0.56 37.96 0.56 37.98 0.56 38.01 0.56 38.03 0.56 38.05 0.56 38.07 0.56 38.1 0.56 38.12 0.56 38.14 0.56 38.17 0.56 38.19 0.56 38.21 0.56 38.24 0.56 38.26 0.56 38.28 0.56 38.31 0.56 38.33 0.56 38.35 0.56 38.37 0.56 38.4 0.56 38.42 0.56 38.44 0.56 38.47 0.56 38.49 0.56 38.51 0.56 38.54 0.56 38.56 0.56 38.58 0.56 38.61 0.56 38.63 0.56 38.65 0.56 38.68 0.56 38.7 0.56 38.72 0.56 38.74 0.56 38.77 0.56 38.79 0.56 38.81 0.56 38.84 0.56 38.86 0.56 38.88 0.56 38.91 0.56 38.93 0.56 38.95 0.56 38.98 0.56 39 0.56 39.02 0.56 39.04 0.56 39.07 0.56 39.09 0.56 39.11 0.56 39.14 0.56 39.16 0.56 39.18 0.56 39.21 0.56 39.23 0.56 39.25 0.56 39.28 0.56 39.3 0.56 39.32 0.56 39.35 0.56 39.37 0.56 39.39 0.56 39.41 0.56 39.44 0.56 39.46 0.56 39.48 0.56 39.51 0.56 39.53 0.56 39.55 0.56 39.58 0.56 39.6 0.56 39.62 0.56 39.65 0.56 39.67 0.56 39.69 0.56 39.71 0.56 39.74 0.56 39.76 0.56 39.78 0.56 39.81 0.56 39.83 0.56 39.85 0.56 39.88 0.56 39.9 0.56 39.92 0.56 39.95 0.56 39.97 0.56 39.99 0.56 40.02 0.56 40.04 0.56 40.06 0.56 40.08 0.56 40.11 0.56 40.13 0.56 40.15 0.56 40.18 0.56 40.2 0.56 40.22 0.56 40.25 0.56 40.27 0.56 40.29 0.56 40.32 0.56 40.34 0.56 40.36 0.56 40.38 0.56 40.41 0.56 40.43 0.56 40.45 0.56 40.48 0.56 40.5 0.56 40.52 0.56 40.55 0.56 40.57 0.56 40.59 0.56 40.62 0.56 40.64 0.56 40.66 0.56 40.69 0.56 40.71 0.56 40.73 0.56 40.75 0.56 40.78 0.56 40.8 0.56 40.82 0.56 40.85 0.56 40.87 0.56 40.89 0.56 40.92 0.56 40.94 0.56 40.96 0.56 40.99 0.56 41.01 0.56 41.03 0.56 41.05 0.56 41.08 0.56 41.1 0.56 41.12 0.56 41.15 0.56 41.17 0.56 41.19 0.56 41.22 0.56 41.24 0.56 41.26 0.56 41.29 0.56 41.31 0.56 41.33 0.56 41.36 0.56 41.38 0.56 41.4 0.56 41.42 0.56 41.45 0.56 41.47 0.56 41.49 0.56 41.52 0.56 41.54 0.56 41.56 0.56 41.59 0.56 41.61 0.56 41.63 0.56 41.66 0.56 41.68 0.56 41.7 0.56 41.72 0.56 41.75 0.56 41.77 0.56 41.79 0.56 41.82 0.56 41.84 0.56 41.86 0.56 41.89 0.56 41.91 0.56 41.93 0.56 41.96 0.56 41.98 0.56 42 0.56 42.03 0.56 42.05 0.56 42.07 0.56 42.09 0.56 42.12 0.56 42.14 0.56 42.16 0.56" class="primitive"/>
          </g>
        </g>
      </g>
      <g opacity="0" class="guide zoomslider" stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-236">
        <g fill="#EAEAEA" stroke-width="0.3" stroke-opacity="0" stroke="#6A6A6A" id="img-8311e7c5-237">
          <g transform="translate(118.14,60)" id="img-8311e7c5-238">
            <path d="M-2,-2 L 2 -2 2 2 -2 2 z" class="primitive"/>
          </g>
          <g class="button_logo" fill="#6A6A6A" id="img-8311e7c5-239">
            <g transform="translate(118.14,60)" id="img-8311e7c5-240">
              <path d="M-1.2,-0.4 L -0.4 -0.4 -0.4 -1.2 0.4 -1.2 0.4 -0.4 1.2 -0.4 1.2 0.4 0.4 0.4 0.4 1.2 -0.4 1.2 -0.4 0.4 -1.2 0.4 z" class="primitive"/>
            </g>
          </g>
        </g>
        <g fill="#EAEAEA" id="img-8311e7c5-241">
          <g transform="translate(106.14,60)" id="img-8311e7c5-242">
            <path d="M-9.5,-2 L 9.5 -2 9.5 2 -9.5 2 z" class="primitive"/>
          </g>
        </g>
        <g class="zoomslider_thumb" fill="#6A6A6A" id="img-8311e7c5-243">
          <g transform="translate(106.14,60)" id="img-8311e7c5-244">
            <path d="M-1,-2 L 1 -2 1 2 -1 2 z" class="primitive"/>
          </g>
        </g>
        <g fill="#EAEAEA" stroke-width="0.3" stroke-opacity="0" stroke="#6A6A6A" id="img-8311e7c5-245">
          <g transform="translate(94.14,60)" id="img-8311e7c5-246">
            <path d="M-2,-2 L 2 -2 2 2 -2 2 z" class="primitive"/>
          </g>
          <g class="button_logo" fill="#6A6A6A" id="img-8311e7c5-247">
            <g transform="translate(94.14,60)" id="img-8311e7c5-248">
              <path d="M-1.2,-0.4 L 1.2 -0.4 1.2 0.4 -1.2 0.4 z" class="primitive"/>
            </g>
          </g>
        </g>
      </g>
    </g>
  </g>
  <g class="guide ylabels" font-size="2.82" font-family="'PT Sans Caption','Helvetica Neue','Helvetica',sans-serif" fill="#6C606B" id="img-8311e7c5-249">
    <g transform="translate(25.72,111.29)" id="img-8311e7c5-250" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-1.5×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,100.43)" id="img-8311e7c5-251" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-1.0×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,89.57)" id="img-8311e7c5-252" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-5.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,78.72)" id="img-8311e7c5-253" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">0</text>
      </g>
    </g>
    <g transform="translate(25.72,67.86)" id="img-8311e7c5-254" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">5.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,57)" id="img-8311e7c5-255" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.0×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,46.14)" id="img-8311e7c5-256" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.5×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,35.28)" id="img-8311e7c5-257" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">2.0×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,24.43)" id="img-8311e7c5-258" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">2.5×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,100.43)" id="img-8311e7c5-259" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-1.00×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,99.34)" id="img-8311e7c5-260" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-9.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,98.26)" id="img-8311e7c5-261" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-9.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,97.17)" id="img-8311e7c5-262" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-8.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,96.09)" id="img-8311e7c5-263" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-8.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,95)" id="img-8311e7c5-264" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-7.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,93.92)" id="img-8311e7c5-265" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-7.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,92.83)" id="img-8311e7c5-266" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-6.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,91.74)" id="img-8311e7c5-267" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-6.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,90.66)" id="img-8311e7c5-268" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-5.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,89.57)" id="img-8311e7c5-269" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-5.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,88.49)" id="img-8311e7c5-270" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-4.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,87.4)" id="img-8311e7c5-271" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-4.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,86.32)" id="img-8311e7c5-272" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-3.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,85.23)" id="img-8311e7c5-273" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-3.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,84.14)" id="img-8311e7c5-274" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-2.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,83.06)" id="img-8311e7c5-275" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-2.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,81.97)" id="img-8311e7c5-276" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-1.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,80.89)" id="img-8311e7c5-277" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-1.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,79.8)" id="img-8311e7c5-278" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-5.00×10²</text>
      </g>
    </g>
    <g transform="translate(25.72,78.72)" id="img-8311e7c5-279" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">0</text>
      </g>
    </g>
    <g transform="translate(25.72,77.63)" id="img-8311e7c5-280" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">5.00×10²</text>
      </g>
    </g>
    <g transform="translate(25.72,76.54)" id="img-8311e7c5-281" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,75.46)" id="img-8311e7c5-282" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,74.37)" id="img-8311e7c5-283" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">2.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,73.29)" id="img-8311e7c5-284" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">2.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,72.2)" id="img-8311e7c5-285" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">3.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,71.11)" id="img-8311e7c5-286" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">3.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,70.03)" id="img-8311e7c5-287" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">4.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,68.94)" id="img-8311e7c5-288" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">4.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,67.86)" id="img-8311e7c5-289" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">5.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,66.77)" id="img-8311e7c5-290" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">5.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,65.69)" id="img-8311e7c5-291" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">6.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,64.6)" id="img-8311e7c5-292" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">6.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,63.51)" id="img-8311e7c5-293" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">7.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,62.43)" id="img-8311e7c5-294" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">7.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,61.34)" id="img-8311e7c5-295" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">8.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,60.26)" id="img-8311e7c5-296" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">8.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,59.17)" id="img-8311e7c5-297" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">9.00×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,58.09)" id="img-8311e7c5-298" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">9.50×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,57)" id="img-8311e7c5-299" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.00×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,55.91)" id="img-8311e7c5-300" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.05×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,54.83)" id="img-8311e7c5-301" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.10×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,53.74)" id="img-8311e7c5-302" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.15×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,52.66)" id="img-8311e7c5-303" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.20×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,51.57)" id="img-8311e7c5-304" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.25×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,50.49)" id="img-8311e7c5-305" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.30×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,49.4)" id="img-8311e7c5-306" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.35×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,48.31)" id="img-8311e7c5-307" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.40×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,47.23)" id="img-8311e7c5-308" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.45×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,46.14)" id="img-8311e7c5-309" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.50×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,45.06)" id="img-8311e7c5-310" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.55×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,43.97)" id="img-8311e7c5-311" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.60×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,42.89)" id="img-8311e7c5-312" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.65×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,41.8)" id="img-8311e7c5-313" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.70×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,40.71)" id="img-8311e7c5-314" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.75×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,39.63)" id="img-8311e7c5-315" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.80×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,38.54)" id="img-8311e7c5-316" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.85×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,37.46)" id="img-8311e7c5-317" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.90×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,36.37)" id="img-8311e7c5-318" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.95×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,35.28)" id="img-8311e7c5-319" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">2.00×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,100.43)" id="img-8311e7c5-320" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-1×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,78.72)" id="img-8311e7c5-321" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">0</text>
      </g>
    </g>
    <g transform="translate(25.72,57)" id="img-8311e7c5-322" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,35.28)" id="img-8311e7c5-323" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">2×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,100.43)" id="img-8311e7c5-324" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-1.0×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,98.26)" id="img-8311e7c5-325" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-9.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,96.09)" id="img-8311e7c5-326" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-8.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,93.92)" id="img-8311e7c5-327" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-7.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,91.74)" id="img-8311e7c5-328" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-6.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,89.57)" id="img-8311e7c5-329" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-5.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,87.4)" id="img-8311e7c5-330" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-4.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,85.23)" id="img-8311e7c5-331" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-3.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,83.06)" id="img-8311e7c5-332" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-2.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,80.89)" id="img-8311e7c5-333" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-1.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,78.72)" id="img-8311e7c5-334" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">0</text>
      </g>
    </g>
    <g transform="translate(25.72,76.54)" id="img-8311e7c5-335" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,74.37)" id="img-8311e7c5-336" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">2.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,72.2)" id="img-8311e7c5-337" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">3.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,70.03)" id="img-8311e7c5-338" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">4.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,67.86)" id="img-8311e7c5-339" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">5.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,65.69)" id="img-8311e7c5-340" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">6.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,63.51)" id="img-8311e7c5-341" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">7.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,61.34)" id="img-8311e7c5-342" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">8.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,59.17)" id="img-8311e7c5-343" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">9.0×10³</text>
      </g>
    </g>
    <g transform="translate(25.72,57)" id="img-8311e7c5-344" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.0×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,54.83)" id="img-8311e7c5-345" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.1×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,52.66)" id="img-8311e7c5-346" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.2×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,50.49)" id="img-8311e7c5-347" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.3×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,48.31)" id="img-8311e7c5-348" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.4×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,46.14)" id="img-8311e7c5-349" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.5×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,43.97)" id="img-8311e7c5-350" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.6×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,41.8)" id="img-8311e7c5-351" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.7×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,39.63)" id="img-8311e7c5-352" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.8×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,37.46)" id="img-8311e7c5-353" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">1.9×10⁴</text>
      </g>
    </g>
    <g transform="translate(25.72,35.28)" id="img-8311e7c5-354" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">2.0×10⁴</text>
      </g>
    </g>
  </g>
  <g font-size="3.88" font-family="'PT Sans','Helvetica Neue','Helvetica',sans-serif" fill="#564A55" stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-355">
    <g transform="translate(8.81,65.86)" id="img-8311e7c5-356">
      <g class="primitive">
        <text text-anchor="middle" dy="0.35em" transform="rotate(-90,0, 2)">value</text>
      </g>
    </g>
  </g>
</g>
<g class="plotroot xscalable yscalable" id="img-8311e7c5-357">
  <g font-size="3.88" font-family="'PT Sans','Helvetica Neue','Helvetica',sans-serif" fill="#564A55" stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-358">
    <g transform="translate(71.39,38.39)" id="img-8311e7c5-359">
      <g class="primitive">
        <text text-anchor="middle" dy="0.6em">t</text>
      </g>
    </g>
  </g>
  <g class="guide xlabels" font-size="2.82" font-family="'PT Sans Caption','Helvetica Neue','Helvetica',sans-serif" fill="#6C606B" id="img-8311e7c5-360">
    <g transform="translate(-102.75,34.39)" id="img-8311e7c5-361" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-500</text>
      </g>
    </g>
    <g transform="translate(-77.88,34.39)" id="img-8311e7c5-362" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-400</text>
      </g>
    </g>
    <g transform="translate(-53,34.39)" id="img-8311e7c5-363" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-300</text>
      </g>
    </g>
    <g transform="translate(-28.12,34.39)" id="img-8311e7c5-364" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-200</text>
      </g>
    </g>
    <g transform="translate(-3.25,34.39)" id="img-8311e7c5-365" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-100</text>
      </g>
    </g>
    <g transform="translate(21.63,34.39)" id="img-8311e7c5-366" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">0</text>
      </g>
    </g>
    <g transform="translate(46.51,34.39)" id="img-8311e7c5-367" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">100</text>
      </g>
    </g>
    <g transform="translate(71.39,34.39)" id="img-8311e7c5-368" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">200</text>
      </g>
    </g>
    <g transform="translate(96.26,34.39)" id="img-8311e7c5-369" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">300</text>
      </g>
    </g>
    <g transform="translate(121.14,34.39)" id="img-8311e7c5-370" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="middle">400</text>
      </g>
    </g>
    <g transform="translate(146.02,34.39)" id="img-8311e7c5-371" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">500</text>
      </g>
    </g>
    <g transform="translate(170.89,34.39)" id="img-8311e7c5-372" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">600</text>
      </g>
    </g>
    <g transform="translate(195.77,34.39)" id="img-8311e7c5-373" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">700</text>
      </g>
    </g>
    <g transform="translate(220.65,34.39)" id="img-8311e7c5-374" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">800</text>
      </g>
    </g>
    <g transform="translate(245.52,34.39)" id="img-8311e7c5-375" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">900</text>
      </g>
    </g>
    <g transform="translate(-77.88,34.39)" id="img-8311e7c5-376" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-400</text>
      </g>
    </g>
    <g transform="translate(-72.9,34.39)" id="img-8311e7c5-377" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-380</text>
      </g>
    </g>
    <g transform="translate(-67.92,34.39)" id="img-8311e7c5-378" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-360</text>
      </g>
    </g>
    <g transform="translate(-62.95,34.39)" id="img-8311e7c5-379" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-340</text>
      </g>
    </g>
    <g transform="translate(-57.97,34.39)" id="img-8311e7c5-380" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-320</text>
      </g>
    </g>
    <g transform="translate(-53,34.39)" id="img-8311e7c5-381" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-300</text>
      </g>
    </g>
    <g transform="translate(-48.02,34.39)" id="img-8311e7c5-382" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-280</text>
      </g>
    </g>
    <g transform="translate(-43.05,34.39)" id="img-8311e7c5-383" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-260</text>
      </g>
    </g>
    <g transform="translate(-38.07,34.39)" id="img-8311e7c5-384" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-240</text>
      </g>
    </g>
    <g transform="translate(-33.1,34.39)" id="img-8311e7c5-385" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-220</text>
      </g>
    </g>
    <g transform="translate(-28.12,34.39)" id="img-8311e7c5-386" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-200</text>
      </g>
    </g>
    <g transform="translate(-23.15,34.39)" id="img-8311e7c5-387" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-180</text>
      </g>
    </g>
    <g transform="translate(-18.17,34.39)" id="img-8311e7c5-388" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-160</text>
      </g>
    </g>
    <g transform="translate(-13.2,34.39)" id="img-8311e7c5-389" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-140</text>
      </g>
    </g>
    <g transform="translate(-8.22,34.39)" id="img-8311e7c5-390" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-120</text>
      </g>
    </g>
    <g transform="translate(-3.25,34.39)" id="img-8311e7c5-391" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-100</text>
      </g>
    </g>
    <g transform="translate(1.73,34.39)" id="img-8311e7c5-392" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-80</text>
      </g>
    </g>
    <g transform="translate(6.71,34.39)" id="img-8311e7c5-393" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-60</text>
      </g>
    </g>
    <g transform="translate(11.68,34.39)" id="img-8311e7c5-394" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-40</text>
      </g>
    </g>
    <g transform="translate(16.66,34.39)" id="img-8311e7c5-395" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-20</text>
      </g>
    </g>
    <g transform="translate(21.63,34.39)" id="img-8311e7c5-396" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">0</text>
      </g>
    </g>
    <g transform="translate(26.61,34.39)" id="img-8311e7c5-397" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">20</text>
      </g>
    </g>
    <g transform="translate(31.58,34.39)" id="img-8311e7c5-398" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">40</text>
      </g>
    </g>
    <g transform="translate(36.56,34.39)" id="img-8311e7c5-399" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">60</text>
      </g>
    </g>
    <g transform="translate(41.53,34.39)" id="img-8311e7c5-400" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">80</text>
      </g>
    </g>
    <g transform="translate(46.51,34.39)" id="img-8311e7c5-401" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">100</text>
      </g>
    </g>
    <g transform="translate(51.48,34.39)" id="img-8311e7c5-402" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">120</text>
      </g>
    </g>
    <g transform="translate(56.46,34.39)" id="img-8311e7c5-403" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">140</text>
      </g>
    </g>
    <g transform="translate(61.43,34.39)" id="img-8311e7c5-404" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">160</text>
      </g>
    </g>
    <g transform="translate(66.41,34.39)" id="img-8311e7c5-405" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">180</text>
      </g>
    </g>
    <g transform="translate(71.39,34.39)" id="img-8311e7c5-406" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">200</text>
      </g>
    </g>
    <g transform="translate(76.36,34.39)" id="img-8311e7c5-407" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">220</text>
      </g>
    </g>
    <g transform="translate(81.34,34.39)" id="img-8311e7c5-408" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">240</text>
      </g>
    </g>
    <g transform="translate(86.31,34.39)" id="img-8311e7c5-409" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">260</text>
      </g>
    </g>
    <g transform="translate(91.29,34.39)" id="img-8311e7c5-410" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">280</text>
      </g>
    </g>
    <g transform="translate(96.26,34.39)" id="img-8311e7c5-411" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">300</text>
      </g>
    </g>
    <g transform="translate(101.24,34.39)" id="img-8311e7c5-412" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">320</text>
      </g>
    </g>
    <g transform="translate(106.21,34.39)" id="img-8311e7c5-413" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">340</text>
      </g>
    </g>
    <g transform="translate(111.19,34.39)" id="img-8311e7c5-414" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">360</text>
      </g>
    </g>
    <g transform="translate(116.16,34.39)" id="img-8311e7c5-415" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">380</text>
      </g>
    </g>
    <g transform="translate(121.14,34.39)" id="img-8311e7c5-416" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">400</text>
      </g>
    </g>
    <g transform="translate(126.11,34.39)" id="img-8311e7c5-417" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">420</text>
      </g>
    </g>
    <g transform="translate(131.09,34.39)" id="img-8311e7c5-418" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">440</text>
      </g>
    </g>
    <g transform="translate(136.06,34.39)" id="img-8311e7c5-419" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">460</text>
      </g>
    </g>
    <g transform="translate(141.04,34.39)" id="img-8311e7c5-420" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">480</text>
      </g>
    </g>
    <g transform="translate(146.02,34.39)" id="img-8311e7c5-421" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">500</text>
      </g>
    </g>
    <g transform="translate(150.99,34.39)" id="img-8311e7c5-422" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">520</text>
      </g>
    </g>
    <g transform="translate(155.97,34.39)" id="img-8311e7c5-423" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">540</text>
      </g>
    </g>
    <g transform="translate(160.94,34.39)" id="img-8311e7c5-424" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">560</text>
      </g>
    </g>
    <g transform="translate(165.92,34.39)" id="img-8311e7c5-425" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">580</text>
      </g>
    </g>
    <g transform="translate(170.89,34.39)" id="img-8311e7c5-426" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">600</text>
      </g>
    </g>
    <g transform="translate(175.87,34.39)" id="img-8311e7c5-427" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">620</text>
      </g>
    </g>
    <g transform="translate(180.84,34.39)" id="img-8311e7c5-428" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">640</text>
      </g>
    </g>
    <g transform="translate(185.82,34.39)" id="img-8311e7c5-429" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">660</text>
      </g>
    </g>
    <g transform="translate(190.79,34.39)" id="img-8311e7c5-430" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">680</text>
      </g>
    </g>
    <g transform="translate(195.77,34.39)" id="img-8311e7c5-431" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">700</text>
      </g>
    </g>
    <g transform="translate(200.74,34.39)" id="img-8311e7c5-432" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">720</text>
      </g>
    </g>
    <g transform="translate(205.72,34.39)" id="img-8311e7c5-433" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">740</text>
      </g>
    </g>
    <g transform="translate(210.7,34.39)" id="img-8311e7c5-434" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">760</text>
      </g>
    </g>
    <g transform="translate(215.67,34.39)" id="img-8311e7c5-435" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">780</text>
      </g>
    </g>
    <g transform="translate(220.65,34.39)" id="img-8311e7c5-436" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">800</text>
      </g>
    </g>
    <g transform="translate(-102.75,34.39)" id="img-8311e7c5-437" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-500</text>
      </g>
    </g>
    <g transform="translate(21.63,34.39)" id="img-8311e7c5-438" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">0</text>
      </g>
    </g>
    <g transform="translate(146.02,34.39)" id="img-8311e7c5-439" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">500</text>
      </g>
    </g>
    <g transform="translate(270.4,34.39)" id="img-8311e7c5-440" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">1000</text>
      </g>
    </g>
    <g transform="translate(-77.88,34.39)" id="img-8311e7c5-441" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-400</text>
      </g>
    </g>
    <g transform="translate(-65.44,34.39)" id="img-8311e7c5-442" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-350</text>
      </g>
    </g>
    <g transform="translate(-53,34.39)" id="img-8311e7c5-443" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-300</text>
      </g>
    </g>
    <g transform="translate(-40.56,34.39)" id="img-8311e7c5-444" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-250</text>
      </g>
    </g>
    <g transform="translate(-28.12,34.39)" id="img-8311e7c5-445" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-200</text>
      </g>
    </g>
    <g transform="translate(-15.68,34.39)" id="img-8311e7c5-446" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-150</text>
      </g>
    </g>
    <g transform="translate(-3.25,34.39)" id="img-8311e7c5-447" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-100</text>
      </g>
    </g>
    <g transform="translate(9.19,34.39)" id="img-8311e7c5-448" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">-50</text>
      </g>
    </g>
    <g transform="translate(21.63,34.39)" id="img-8311e7c5-449" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">0</text>
      </g>
    </g>
    <g transform="translate(34.07,34.39)" id="img-8311e7c5-450" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">50</text>
      </g>
    </g>
    <g transform="translate(46.51,34.39)" id="img-8311e7c5-451" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">100</text>
      </g>
    </g>
    <g transform="translate(58.95,34.39)" id="img-8311e7c5-452" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">150</text>
      </g>
    </g>
    <g transform="translate(71.39,34.39)" id="img-8311e7c5-453" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">200</text>
      </g>
    </g>
    <g transform="translate(83.82,34.39)" id="img-8311e7c5-454" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">250</text>
      </g>
    </g>
    <g transform="translate(96.26,34.39)" id="img-8311e7c5-455" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">300</text>
      </g>
    </g>
    <g transform="translate(108.7,34.39)" id="img-8311e7c5-456" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">350</text>
      </g>
    </g>
    <g transform="translate(121.14,34.39)" id="img-8311e7c5-457" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">400</text>
      </g>
    </g>
    <g transform="translate(133.58,34.39)" id="img-8311e7c5-458" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">450</text>
      </g>
    </g>
    <g transform="translate(146.02,34.39)" id="img-8311e7c5-459" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">500</text>
      </g>
    </g>
    <g transform="translate(158.45,34.39)" id="img-8311e7c5-460" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">550</text>
      </g>
    </g>
    <g transform="translate(170.89,34.39)" id="img-8311e7c5-461" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">600</text>
      </g>
    </g>
    <g transform="translate(183.33,34.39)" id="img-8311e7c5-462" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">650</text>
      </g>
    </g>
    <g transform="translate(195.77,34.39)" id="img-8311e7c5-463" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">700</text>
      </g>
    </g>
    <g transform="translate(208.21,34.39)" id="img-8311e7c5-464" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">750</text>
      </g>
    </g>
    <g transform="translate(220.65,34.39)" id="img-8311e7c5-465" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="middle">800</text>
      </g>
    </g>
  </g>
  <g class="guide colorkey" id="img-8311e7c5-466">
    <g fill="#4C404B" font-size="2.82" font-family="'PT Sans','Helvetica Neue','Helvetica',sans-serif" id="img-8311e7c5-467">
      <g transform="translate(126.95,14.23)" id="img-8311e7c5-468" class="color_S_H">
        <g class="primitive">
          <text dy="0.35em">S_H</text>
        </g>
      </g>
      <g transform="translate(126.95,17.86)" id="img-8311e7c5-469" class="color_E_H">
        <g class="primitive">
          <text dy="0.35em">E_H</text>
        </g>
      </g>
      <g transform="translate(126.95,21.48)" id="img-8311e7c5-470" class="color_I_H">
        <g class="primitive">
          <text dy="0.35em">I_H</text>
        </g>
      </g>
      <g transform="translate(126.95,25.11)" id="img-8311e7c5-471" class="color_R_H">
        <g class="primitive">
          <text dy="0.35em">R_H</text>
        </g>
      </g>
    </g>
    <g stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-472">
      <g transform="translate(125.05,14.23)" id="img-8311e7c5-473" fill="#00BFFF" class="color_S_H">
        <path d="M-0.91,-0.91 L 0.91 -0.91 0.91 0.91 -0.91 0.91 z" class="primitive"/>
      </g>
      <g transform="translate(125.05,17.86)" id="img-8311e7c5-474" fill="#D4CA3A" class="color_E_H">
        <path d="M-0.91,-0.91 L 0.91 -0.91 0.91 0.91 -0.91 0.91 z" class="primitive"/>
      </g>
      <g transform="translate(125.05,21.48)" id="img-8311e7c5-475" fill="#FF6DAE" class="color_I_H">
        <path d="M-0.91,-0.91 L 0.91 -0.91 0.91 0.91 -0.91 0.91 z" class="primitive"/>
      </g>
      <g transform="translate(125.05,25.11)" id="img-8311e7c5-476" fill="#00B78D" class="color_R_H">
        <path d="M-0.91,-0.91 L 0.91 -0.91 0.91 0.91 -0.91 0.91 z" class="primitive"/>
      </g>
    </g>
    <g fill="#362A35" font-size="3.88" font-family="'PT Sans','Helvetica Neue','Helvetica',sans-serif" stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-477">
      <g transform="translate(124.14,10.41)" id="img-8311e7c5-478">
        <g class="primitive">
          <text>variable</text>
        </g>
      </g>
    </g>
  </g>
  <g clip-path="url(#img-8311e7c5-479)">
    <g id="img-8311e7c5-480">
      <g pointer-events="visible" opacity="1" fill="#000000" fill-opacity="0.000" stroke="#000000" stroke-opacity="0.000" class="guide background" id="img-8311e7c5-481">
        <g transform="translate(71.39,17.86)" id="img-8311e7c5-482">
          <path d="M-51.75,-12.86 L 51.75 -12.86 51.75 12.86 -51.75 12.86 z" class="primitive"/>
        </g>
      </g>
      <g class="guide ygridlines xfixed" stroke-dasharray="0.5,0.5" stroke-width="0.2" stroke="#D0D0E0" id="img-8311e7c5-483">
        <g transform="translate(71.39,57.67)" id="img-8311e7c5-484" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,50.43)" id="img-8311e7c5-485" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,43.19)" id="img-8311e7c5-486" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,35.95)" id="img-8311e7c5-487" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,28.72)" id="img-8311e7c5-488" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,21.48)" id="img-8311e7c5-489" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,14.24)" id="img-8311e7c5-490" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,7)" id="img-8311e7c5-491" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-0.24)" id="img-8311e7c5-492" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-7.48)" id="img-8311e7c5-493" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-14.71)" id="img-8311e7c5-494" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-21.95)" id="img-8311e7c5-495" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,50.43)" id="img-8311e7c5-496" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,49.71)" id="img-8311e7c5-497" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,48.98)" id="img-8311e7c5-498" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,48.26)" id="img-8311e7c5-499" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,47.53)" id="img-8311e7c5-500" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,46.81)" id="img-8311e7c5-501" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,46.09)" id="img-8311e7c5-502" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,45.36)" id="img-8311e7c5-503" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,44.64)" id="img-8311e7c5-504" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,43.92)" id="img-8311e7c5-505" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,43.19)" id="img-8311e7c5-506" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,42.47)" id="img-8311e7c5-507" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,41.74)" id="img-8311e7c5-508" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,41.02)" id="img-8311e7c5-509" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,40.3)" id="img-8311e7c5-510" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,39.57)" id="img-8311e7c5-511" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,38.85)" id="img-8311e7c5-512" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,38.12)" id="img-8311e7c5-513" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,37.4)" id="img-8311e7c5-514" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,36.68)" id="img-8311e7c5-515" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,35.95)" id="img-8311e7c5-516" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,35.23)" id="img-8311e7c5-517" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,34.51)" id="img-8311e7c5-518" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,33.78)" id="img-8311e7c5-519" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,33.06)" id="img-8311e7c5-520" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,32.33)" id="img-8311e7c5-521" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,31.61)" id="img-8311e7c5-522" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,30.89)" id="img-8311e7c5-523" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,30.16)" id="img-8311e7c5-524" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,29.44)" id="img-8311e7c5-525" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,28.72)" id="img-8311e7c5-526" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,27.99)" id="img-8311e7c5-527" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,27.27)" id="img-8311e7c5-528" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,26.54)" id="img-8311e7c5-529" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,25.82)" id="img-8311e7c5-530" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,25.1)" id="img-8311e7c5-531" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,24.37)" id="img-8311e7c5-532" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,23.65)" id="img-8311e7c5-533" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,22.92)" id="img-8311e7c5-534" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,22.2)" id="img-8311e7c5-535" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,21.48)" id="img-8311e7c5-536" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,20.75)" id="img-8311e7c5-537" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,20.03)" id="img-8311e7c5-538" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,19.31)" id="img-8311e7c5-539" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,18.58)" id="img-8311e7c5-540" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,17.86)" id="img-8311e7c5-541" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,17.13)" id="img-8311e7c5-542" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,16.41)" id="img-8311e7c5-543" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,15.69)" id="img-8311e7c5-544" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,14.96)" id="img-8311e7c5-545" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,14.24)" id="img-8311e7c5-546" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,13.51)" id="img-8311e7c5-547" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,12.79)" id="img-8311e7c5-548" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,12.07)" id="img-8311e7c5-549" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,11.34)" id="img-8311e7c5-550" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,10.62)" id="img-8311e7c5-551" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,9.9)" id="img-8311e7c5-552" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,9.17)" id="img-8311e7c5-553" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,8.45)" id="img-8311e7c5-554" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,7.72)" id="img-8311e7c5-555" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,7)" id="img-8311e7c5-556" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,6.28)" id="img-8311e7c5-557" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,5.55)" id="img-8311e7c5-558" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,4.83)" id="img-8311e7c5-559" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,4.1)" id="img-8311e7c5-560" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,3.38)" id="img-8311e7c5-561" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,2.66)" id="img-8311e7c5-562" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,1.93)" id="img-8311e7c5-563" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,1.21)" id="img-8311e7c5-564" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,0.49)" id="img-8311e7c5-565" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-0.24)" id="img-8311e7c5-566" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-0.96)" id="img-8311e7c5-567" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-1.69)" id="img-8311e7c5-568" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-2.41)" id="img-8311e7c5-569" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-3.13)" id="img-8311e7c5-570" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-3.86)" id="img-8311e7c5-571" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-4.58)" id="img-8311e7c5-572" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-5.31)" id="img-8311e7c5-573" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-6.03)" id="img-8311e7c5-574" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-6.75)" id="img-8311e7c5-575" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-7.48)" id="img-8311e7c5-576" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-8.2)" id="img-8311e7c5-577" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-8.92)" id="img-8311e7c5-578" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-9.65)" id="img-8311e7c5-579" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-10.37)" id="img-8311e7c5-580" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-11.1)" id="img-8311e7c5-581" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-11.82)" id="img-8311e7c5-582" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-12.54)" id="img-8311e7c5-583" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-13.27)" id="img-8311e7c5-584" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-13.99)" id="img-8311e7c5-585" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-14.71)" id="img-8311e7c5-586" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,57.67)" id="img-8311e7c5-587" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,28.72)" id="img-8311e7c5-588" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-0.24)" id="img-8311e7c5-589" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-29.19)" id="img-8311e7c5-590" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,50.43)" id="img-8311e7c5-591" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,48.98)" id="img-8311e7c5-592" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,47.53)" id="img-8311e7c5-593" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,46.09)" id="img-8311e7c5-594" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,44.64)" id="img-8311e7c5-595" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,43.19)" id="img-8311e7c5-596" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,41.74)" id="img-8311e7c5-597" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,40.3)" id="img-8311e7c5-598" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,38.85)" id="img-8311e7c5-599" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,37.4)" id="img-8311e7c5-600" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,35.95)" id="img-8311e7c5-601" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,34.51)" id="img-8311e7c5-602" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,33.06)" id="img-8311e7c5-603" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,31.61)" id="img-8311e7c5-604" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,30.16)" id="img-8311e7c5-605" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,28.72)" id="img-8311e7c5-606" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,27.27)" id="img-8311e7c5-607" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,25.82)" id="img-8311e7c5-608" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,24.37)" id="img-8311e7c5-609" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,22.92)" id="img-8311e7c5-610" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,21.48)" id="img-8311e7c5-611" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,20.03)" id="img-8311e7c5-612" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,18.58)" id="img-8311e7c5-613" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,17.13)" id="img-8311e7c5-614" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,15.69)" id="img-8311e7c5-615" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,14.24)" id="img-8311e7c5-616" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,12.79)" id="img-8311e7c5-617" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,11.34)" id="img-8311e7c5-618" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,9.9)" id="img-8311e7c5-619" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,8.45)" id="img-8311e7c5-620" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,7)" id="img-8311e7c5-621" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,5.55)" id="img-8311e7c5-622" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,4.1)" id="img-8311e7c5-623" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,2.66)" id="img-8311e7c5-624" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,1.21)" id="img-8311e7c5-625" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-0.24)" id="img-8311e7c5-626" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-1.69)" id="img-8311e7c5-627" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-3.13)" id="img-8311e7c5-628" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-4.58)" id="img-8311e7c5-629" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-6.03)" id="img-8311e7c5-630" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-7.48)" id="img-8311e7c5-631" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-8.92)" id="img-8311e7c5-632" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-10.37)" id="img-8311e7c5-633" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-11.82)" id="img-8311e7c5-634" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-13.27)" id="img-8311e7c5-635" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
        <g transform="translate(71.39,-14.71)" id="img-8311e7c5-636" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M-51.75,0 L 51.75 0" class="primitive"/>
        </g>
      </g>
      <g class="guide xgridlines yfixed" stroke-dasharray="0.5,0.5" stroke-width="0.2" stroke="#D0D0E0" id="img-8311e7c5-637">
        <g transform="translate(-102.75,17.86)" id="img-8311e7c5-638" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-77.88,17.86)" id="img-8311e7c5-639" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-53,17.86)" id="img-8311e7c5-640" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-28.12,17.86)" id="img-8311e7c5-641" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-3.25,17.86)" id="img-8311e7c5-642" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(21.63,17.86)" id="img-8311e7c5-643" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(46.51,17.86)" id="img-8311e7c5-644" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(71.39,17.86)" id="img-8311e7c5-645" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(96.26,17.86)" id="img-8311e7c5-646" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(121.14,17.86)" id="img-8311e7c5-647" gadfly:scale="1.0" visibility="visible">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(146.02,17.86)" id="img-8311e7c5-648" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(170.89,17.86)" id="img-8311e7c5-649" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(195.77,17.86)" id="img-8311e7c5-650" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(220.65,17.86)" id="img-8311e7c5-651" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(245.52,17.86)" id="img-8311e7c5-652" gadfly:scale="1.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-77.88,17.86)" id="img-8311e7c5-653" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-72.9,17.86)" id="img-8311e7c5-654" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-67.92,17.86)" id="img-8311e7c5-655" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-62.95,17.86)" id="img-8311e7c5-656" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-57.97,17.86)" id="img-8311e7c5-657" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-53,17.86)" id="img-8311e7c5-658" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-48.02,17.86)" id="img-8311e7c5-659" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-43.05,17.86)" id="img-8311e7c5-660" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-38.07,17.86)" id="img-8311e7c5-661" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-33.1,17.86)" id="img-8311e7c5-662" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-28.12,17.86)" id="img-8311e7c5-663" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-23.15,17.86)" id="img-8311e7c5-664" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-18.17,17.86)" id="img-8311e7c5-665" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-13.2,17.86)" id="img-8311e7c5-666" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-8.22,17.86)" id="img-8311e7c5-667" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-3.25,17.86)" id="img-8311e7c5-668" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(1.73,17.86)" id="img-8311e7c5-669" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(6.71,17.86)" id="img-8311e7c5-670" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(11.68,17.86)" id="img-8311e7c5-671" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(16.66,17.86)" id="img-8311e7c5-672" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(21.63,17.86)" id="img-8311e7c5-673" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(26.61,17.86)" id="img-8311e7c5-674" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(31.58,17.86)" id="img-8311e7c5-675" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(36.56,17.86)" id="img-8311e7c5-676" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(41.53,17.86)" id="img-8311e7c5-677" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(46.51,17.86)" id="img-8311e7c5-678" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(51.48,17.86)" id="img-8311e7c5-679" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(56.46,17.86)" id="img-8311e7c5-680" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(61.43,17.86)" id="img-8311e7c5-681" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(66.41,17.86)" id="img-8311e7c5-682" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(71.39,17.86)" id="img-8311e7c5-683" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(76.36,17.86)" id="img-8311e7c5-684" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(81.34,17.86)" id="img-8311e7c5-685" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(86.31,17.86)" id="img-8311e7c5-686" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(91.29,17.86)" id="img-8311e7c5-687" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(96.26,17.86)" id="img-8311e7c5-688" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(101.24,17.86)" id="img-8311e7c5-689" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(106.21,17.86)" id="img-8311e7c5-690" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(111.19,17.86)" id="img-8311e7c5-691" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(116.16,17.86)" id="img-8311e7c5-692" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(121.14,17.86)" id="img-8311e7c5-693" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(126.11,17.86)" id="img-8311e7c5-694" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(131.09,17.86)" id="img-8311e7c5-695" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(136.06,17.86)" id="img-8311e7c5-696" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(141.04,17.86)" id="img-8311e7c5-697" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(146.02,17.86)" id="img-8311e7c5-698" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(150.99,17.86)" id="img-8311e7c5-699" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(155.97,17.86)" id="img-8311e7c5-700" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(160.94,17.86)" id="img-8311e7c5-701" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(165.92,17.86)" id="img-8311e7c5-702" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(170.89,17.86)" id="img-8311e7c5-703" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(175.87,17.86)" id="img-8311e7c5-704" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(180.84,17.86)" id="img-8311e7c5-705" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(185.82,17.86)" id="img-8311e7c5-706" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(190.79,17.86)" id="img-8311e7c5-707" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(195.77,17.86)" id="img-8311e7c5-708" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(200.74,17.86)" id="img-8311e7c5-709" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(205.72,17.86)" id="img-8311e7c5-710" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(210.7,17.86)" id="img-8311e7c5-711" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(215.67,17.86)" id="img-8311e7c5-712" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(220.65,17.86)" id="img-8311e7c5-713" gadfly:scale="10.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-102.75,17.86)" id="img-8311e7c5-714" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(21.63,17.86)" id="img-8311e7c5-715" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(146.02,17.86)" id="img-8311e7c5-716" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(270.4,17.86)" id="img-8311e7c5-717" gadfly:scale="0.5" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-77.88,17.86)" id="img-8311e7c5-718" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-65.44,17.86)" id="img-8311e7c5-719" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-53,17.86)" id="img-8311e7c5-720" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-40.56,17.86)" id="img-8311e7c5-721" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-28.12,17.86)" id="img-8311e7c5-722" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-15.68,17.86)" id="img-8311e7c5-723" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(-3.25,17.86)" id="img-8311e7c5-724" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(9.19,17.86)" id="img-8311e7c5-725" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(21.63,17.86)" id="img-8311e7c5-726" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(34.07,17.86)" id="img-8311e7c5-727" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(46.51,17.86)" id="img-8311e7c5-728" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(58.95,17.86)" id="img-8311e7c5-729" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(71.39,17.86)" id="img-8311e7c5-730" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(83.82,17.86)" id="img-8311e7c5-731" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(96.26,17.86)" id="img-8311e7c5-732" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(108.7,17.86)" id="img-8311e7c5-733" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(121.14,17.86)" id="img-8311e7c5-734" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(133.58,17.86)" id="img-8311e7c5-735" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(146.02,17.86)" id="img-8311e7c5-736" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(158.45,17.86)" id="img-8311e7c5-737" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(170.89,17.86)" id="img-8311e7c5-738" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(183.33,17.86)" id="img-8311e7c5-739" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(195.77,17.86)" id="img-8311e7c5-740" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(208.21,17.86)" id="img-8311e7c5-741" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
        <g transform="translate(220.65,17.86)" id="img-8311e7c5-742" gadfly:scale="5.0" visibility="hidden">
          <path fill="none" d="M0,-12.86 L 0 12.86" class="primitive"/>
        </g>
      </g>
      <g class="plotpanel" id="img-8311e7c5-743">
        <metadata>
          <boundingbox value="19.63166666666666mm 5.0mm 103.50718957064284mm 25.715000000000003mm"/>
          <unitbox value="-8.039620086265812 163.8153350218743 416.07924017253157 -177.63067004374858"/>
        </metadata>
        <g stroke-width="0.3" fill="#000000" fill-opacity="0.000" stroke-dasharray="none" id="img-8311e7c5-744">
          <g transform="translate(67.03,16.41)" id="img-8311e7c5-745" class="geometry color_R_H" stroke="#00B78D">
            <path fill="none" d="M-45.4,12.3 L -45.38 12.3 -45.35 12.3 -45.33 12.3 -45.3 12.3 -45.28 12.3 -45.25 12.3 -45.23 12.3 -45.2 12.29 -45.18 12.29 -45.15 12.29 -45.13 12.29 -45.1 12.29 -45.08 12.29 -45.05 12.29 -45.03 12.29 -45 12.29 -44.98 12.29 -44.95 12.28 -44.93 12.28 -44.9 12.28 -44.88 12.28 -44.85 12.28 -44.83 12.28 -44.8 12.28 -44.78 12.28 -44.75 12.28 -44.73 12.28 -44.7 12.28 -44.68 12.28 -44.65 12.27 -44.63 12.27 -44.6 12.27 -44.58 12.27 -44.55 12.27 -44.53 12.27 -44.5 12.27 -44.48 12.27 -44.45 12.27 -44.43 12.27 -44.41 12.27 -44.38 12.27 -44.36 12.26 -44.33 12.26 -44.31 12.26 -44.28 12.26 -44.26 12.26 -44.23 12.26 -44.21 12.26 -44.18 12.26 -44.16 12.26 -44.13 12.26 -44.11 12.26 -44.08 12.26 -44.06 12.25 -44.03 12.25 -44.01 12.25 -43.98 12.25 -43.96 12.25 -43.93 12.25 -43.91 12.25 -43.88 12.25 -43.86 12.25 -43.83 12.25 -43.81 12.24 -43.78 12.24 -43.76 12.24 -43.73 12.24 -43.71 12.24 -43.68 12.24 -43.66 12.24 -43.63 12.24 -43.61 12.24 -43.58 12.24 -43.56 12.23 -43.53 12.23 -43.51 12.23 -43.48 12.23 -43.46 12.23 -43.43 12.23 -43.41 12.23 -43.39 12.23 -43.36 12.22 -43.34 12.22 -43.31 12.22 -43.29 12.22 -43.26 12.22 -43.24 12.22 -43.21 12.22 -43.19 12.21 -43.16 12.21 -43.14 12.21 -43.11 12.21 -43.09 12.21 -43.06 12.21 -43.04 12.2 -43.01 12.2 -42.99 12.2 -42.96 12.2 -42.94 12.2 -42.91 12.2 -42.89 12.19 -42.86 12.19 -42.84 12.19 -42.81 12.19 -42.79 12.19 -42.76 12.18 -42.74 12.18 -42.71 12.18 -42.69 12.18 -42.66 12.17 -42.64 12.17 -42.61 12.17 -42.59 12.17 -42.56 12.16 -42.54 12.16 -42.51 12.16 -42.49 12.16 -42.46 12.15 -42.44 12.15 -42.41 12.15 -42.39 12.15 -42.37 12.14 -42.34 12.14 -42.32 12.14 -42.29 12.13 -42.27 12.13 -42.24 12.13 -42.22 12.12 -42.19 12.12 -42.17 12.12 -42.14 12.11 -42.12 12.11 -42.09 12.11 -42.07 12.1 -42.04 12.1 -42.02 12.09 -41.99 12.09 -41.97 12.09 -41.94 12.08 -41.92 12.08 -41.89 12.07 -41.87 12.07 -41.84 12.07 -41.82 12.06 -41.79 12.06 -41.77 12.05 -41.74 12.05 -41.72 12.04 -41.69 12.04 -41.67 12.03 -41.64 12.03 -41.62 12.02 -41.59 12.02 -41.57 12.01 -41.54 12.01 -41.52 12 -41.49 11.99 -41.47 11.99 -41.44 11.98 -41.42 11.98 -41.39 11.97 -41.37 11.96 -41.35 11.96 -41.32 11.95 -41.3 11.94 -41.27 11.94 -41.25 11.93 -41.22 11.92 -41.2 11.92 -41.17 11.91 -41.15 11.9 -41.12 11.89 -41.1 11.89 -41.07 11.88 -41.05 11.87 -41.02 11.86 -41 11.85 -40.97 11.85 -40.95 11.84 -40.92 11.83 -40.9 11.82 -40.87 11.81 -40.85 11.8 -40.82 11.79 -40.8 11.78 -40.77 11.77 -40.75 11.76 -40.72 11.75 -40.7 11.74 -40.67 11.73 -40.65 11.72 -40.62 11.71 -40.6 11.7 -40.57 11.69 -40.55 11.68 -40.52 11.67 -40.5 11.66 -40.47 11.64 -40.45 11.63 -40.42 11.62 -40.4 11.61 -40.38 11.6 -40.35 11.58 -40.33 11.57 -40.3 11.56 -40.28 11.54 -40.25 11.53 -40.23 11.52 -40.2 11.5 -40.18 11.49 -40.15 11.47 -40.13 11.46 -40.1 11.44 -40.08 11.43 -40.05 11.41 -40.03 11.4 -40 11.38 -39.98 11.36 -39.95 11.35 -39.93 11.33 -39.9 11.31 -39.88 11.3 -39.85 11.28 -39.83 11.26 -39.8 11.24 -39.78 11.23 -39.75 11.21 -39.73 11.19 -39.7 11.17 -39.68 11.15 -39.65 11.13 -39.63 11.11 -39.6 11.09 -39.58 11.07 -39.55 11.05 -39.53 11.03 -39.5 11 -39.48 10.98 -39.45 10.96 -39.43 10.94 -39.4 10.92 -39.38 10.89 -39.36 10.87 -39.33 10.85 -39.31 10.82 -39.28 10.8 -39.26 10.77 -39.23 10.75 -39.21 10.72 -39.18 10.7 -39.16 10.67 -39.13 10.64 -39.11 10.62 -39.08 10.59 -39.06 10.56 -39.03 10.54 -39.01 10.51 -38.98 10.48 -38.96 10.45 -38.93 10.42 -38.91 10.39 -38.88 10.36 -38.86 10.33 -38.83 10.3 -38.81 10.27 -38.78 10.24 -38.76 10.21 -38.73 10.18 -38.71 10.15 -38.68 10.11 -38.66 10.08 -38.63 10.05 -38.61 10.01 -38.58 9.98 -38.56 9.95 -38.53 9.91 -38.51 9.88 -38.48 9.84 -38.46 9.81 -38.43 9.77 -38.41 9.73 -38.38 9.7 -38.36 9.66 -38.34 9.62 -38.31 9.59 -38.29 9.55 -38.26 9.51 -38.24 9.47 -38.21 9.43 -38.19 9.39 -38.16 9.35 -38.14 9.31 -38.11 9.27 -38.09 9.23 -38.06 9.19 -38.04 9.15 -38.01 9.11 -37.99 9.07 -37.96 9.03 -37.94 8.98 -37.91 8.94 -37.89 8.9 -37.86 8.86 -37.84 8.81 -37.81 8.77 -37.79 8.73 -37.76 8.68 -37.74 8.64 -37.71 8.59 -37.69 8.55 -37.66 8.5 -37.64 8.46 -37.61 8.41 -37.59 8.37 -37.56 8.32 -37.54 8.28 -37.51 8.23 -37.49 8.18 -37.46 8.14 -37.44 8.09 -37.41 8.04 -37.39 8 -37.36 7.95 -37.34 7.9 -37.32 7.86 -37.29 7.81 -37.27 7.76 -37.24 7.71 -37.22 7.67 -37.19 7.62 -37.17 7.57 -37.14 7.52 -37.12 7.47 -37.09 7.42 -37.07 7.38 -37.04 7.33 -37.02 7.28 -36.99 7.23 -36.97 7.18 -36.94 7.13 -36.92 7.08 -36.89 7.04 -36.87 6.99 -36.84 6.94 -36.82 6.89 -36.79 6.84 -36.77 6.79 -36.74 6.74 -36.72 6.7 -36.69 6.65 -36.67 6.6 -36.64 6.55 -36.62 6.5 -36.59 6.45 -36.57 6.4 -36.54 6.36 -36.52 6.31 -36.49 6.26 -36.47 6.21 -36.44 6.16 -36.42 6.12 -36.39 6.07 -36.37 6.02 -36.35 5.97 -36.32 5.92 -36.3 5.88 -36.27 5.83 -36.25 5.78 -36.22 5.74 -36.2 5.69 -36.17 5.64 -36.15 5.6 -36.12 5.55 -36.1 5.5 -36.07 5.46 -36.05 5.41 -36.02 5.36 -36 5.32 -35.97 5.27 -35.95 5.23 -35.92 5.18 -35.9 5.14 -35.87 5.09 -35.85 5.05 -35.82 5 -35.8 4.96 -35.77 4.92 -35.75 4.87 -35.72 4.83 -35.7 4.78 -35.67 4.74 -35.65 4.7 -35.62 4.65 -35.6 4.61 -35.57 4.57 -35.55 4.53 -35.52 4.49 -35.5 4.44 -35.47 4.4 -35.45 4.36 -35.42 4.32 -35.4 4.28 -35.37 4.24 -35.35 4.2 -35.33 4.16 -35.3 4.12 -35.28 4.08 -35.25 4.04 -35.23 4 -35.2 3.96 -35.18 3.92 -35.15 3.88 -35.13 3.84 -35.1 3.8 -35.08 3.76 -35.05 3.73 -35.03 3.69 -35 3.65 -34.98 3.61 -34.95 3.58 -34.93 3.54 -34.9 3.5 -34.88 3.47 -34.85 3.43 -34.83 3.39 -34.8 3.36 -34.78 3.32 -34.75 3.29 -34.73 3.25 -34.7 3.22 -34.68 3.18 -34.65 3.15 -34.63 3.12 -34.6 3.08 -34.58 3.05 -34.55 3.01 -34.53 2.98 -34.5 2.95 -34.48 2.91 -34.45 2.88 -34.43 2.85 -34.4 2.82 -34.38 2.79 -34.35 2.75 -34.33 2.72 -34.31 2.69 -34.28 2.66 -34.26 2.63 -34.23 2.6 -34.21 2.57 -34.18 2.54 -34.16 2.51 -34.13 2.48 -34.11 2.45 -34.08 2.42 -34.06 2.39 -34.03 2.36 -34.01 2.33 -33.98 2.3 -33.96 2.27 -33.93 2.24 -33.91 2.22 -33.88 2.19 -33.86 2.16 -33.83 2.13 -33.81 2.11 -33.78 2.08 -33.76 2.05 -33.73 2.02 -33.71 2 -33.68 1.97 -33.66 1.94 -33.63 1.92 -33.61 1.89 -33.58 1.87 -33.56 1.84 -33.53 1.82 -33.51 1.79 -33.48 1.77 -33.46 1.74 -33.43 1.72 -33.41 1.69 -33.38 1.67 -33.36 1.64 -33.33 1.62 -33.31 1.6 -33.29 1.57 -33.26 1.55 -33.24 1.53 -33.21 1.5 -33.19 1.48 -33.16 1.46 -33.14 1.43 -33.11 1.41 -33.09 1.39 -33.06 1.37 -33.04 1.35 -33.01 1.32 -32.99 1.3 -32.96 1.28 -32.94 1.26 -32.91 1.24 -32.89 1.22 -32.86 1.2 -32.84 1.18 -32.81 1.15 -32.79 1.13 -32.76 1.11 -32.74 1.09 -32.71 1.07 -32.69 1.05 -32.66 1.03 -32.64 1.01 -32.61 1 -32.59 0.98 -32.56 0.96 -32.54 0.94 -32.51 0.92 -32.49 0.9 -32.46 0.88 -32.44 0.86 -32.41 0.84 -32.39 0.83 -32.36 0.81 -32.34 0.79 -32.31 0.77 -32.29 0.75 -32.27 0.74 -32.24 0.72 -32.22 0.7 -32.19 0.69 -32.17 0.67 -32.14 0.65 -32.12 0.63 -32.09 0.62 -32.07 0.6 -32.04 0.59 -32.02 0.57 -31.99 0.55 -31.97 0.54 -31.94 0.52 -31.92 0.5 -31.89 0.49 -31.87 0.47 -31.84 0.46 -31.82 0.44 -31.79 0.43 -31.77 0.41 -31.74 0.4 -31.72 0.38 -31.69 0.37 -31.67 0.35 -31.64 0.34 -31.62 0.32 -31.59 0.31 -31.57 0.29 -31.54 0.28 -31.52 0.27 -31.49 0.25 -31.47 0.24 -31.44 0.22 -31.42 0.21 -31.39 0.2 -31.37 0.18 -31.34 0.17 -31.32 0.16 -31.3 0.14 -31.27 0.13 -31.25 0.12 -31.22 0.11 -31.2 0.09 -31.17 0.08 -31.15 0.07 -31.12 0.05 -31.1 0.04 -31.07 0.03 -31.05 0.02 -31.02 0.01 -31 -0.01 -30.97 -0.02 -30.95 -0.03 -30.92 -0.04 -30.9 -0.05 -30.87 -0.07 -30.85 -0.08 -30.82 -0.09 -30.8 -0.1 -30.77 -0.11 -30.75 -0.12 -30.72 -0.13 -30.7 -0.14 -30.67 -0.16 -30.65 -0.17 -30.62 -0.18 -30.6 -0.19 -30.57 -0.2 -30.55 -0.21 -30.52 -0.22 -30.5 -0.23 -30.47 -0.24 -30.45 -0.25 -30.42 -0.26 -30.4 -0.27 -30.37 -0.28 -30.35 -0.29 -30.32 -0.3 -30.3 -0.31 -30.28 -0.32 -30.25 -0.33 -30.23 -0.34 -30.2 -0.35 -30.18 -0.36 -30.15 -0.37 -30.13 -0.38 -30.1 -0.39 -30.08 -0.4 -30.05 -0.4 -30.03 -0.41 -30 -0.42 -29.98 -0.43 -29.95 -0.44 -29.93 -0.45 -29.9 -0.46 -29.88 -0.47 -29.85 -0.47 -29.83 -0.48 -29.8 -0.49 -29.78 -0.5 -29.75 -0.51 -29.73 -0.52 -29.7 -0.53 -29.68 -0.53 -29.65 -0.54 -29.63 -0.55 -29.6 -0.56 -29.58 -0.57 -29.55 -0.57 -29.53 -0.58 -29.5 -0.59 -29.48 -0.6 -29.45 -0.6 -29.43 -0.61 -29.4 -0.62 -29.38 -0.63 -29.35 -0.63 -29.33 -0.64 -29.3 -0.65 -29.28 -0.66 -29.26 -0.66 -29.23 -0.67 -29.21 -0.68 -29.18 -0.68 -29.16 -0.69 -29.13 -0.7 -29.11 -0.7 -29.08 -0.71 -29.06 -0.72 -29.03 -0.73 -29.01 -0.73 -28.98 -0.74 -28.96 -0.74 -28.93 -0.75 -28.91 -0.76 -28.88 -0.76 -28.86 -0.77 -28.83 -0.78 -28.81 -0.78 -28.78 -0.79 -28.76 -0.8 -28.73 -0.8 -28.71 -0.81 -28.68 -0.81 -28.66 -0.82 -28.63 -0.83 -28.61 -0.83 -28.58 -0.84 -28.56 -0.84 -28.53 -0.85 -28.51 -0.86 -28.48 -0.86 -28.46 -0.87 -28.43 -0.87 -28.41 -0.88 -28.38 -0.88 -28.36 -0.89 -28.33 -0.89 -28.31 -0.9 -28.28 -0.9 -28.26 -0.91 -28.24 -0.92 -28.21 -0.92 -28.19 -0.93 -28.16 -0.93 -28.14 -0.94 -28.11 -0.94 -28.09 -0.95 -28.06 -0.95 -28.04 -0.96 -28.01 -0.96 -27.99 -0.97 -27.96 -0.97 -27.94 -0.98 -27.91 -0.98 -27.89 -0.99 -27.86 -0.99 -27.84 -0.99 -27.81 -1 -27.79 -1 -27.76 -1.01 -27.74 -1.01 -27.71 -1.02 -27.69 -1.02 -27.66 -1.03 -27.64 -1.03 -27.61 -1.04 -27.59 -1.04 -27.56 -1.04 -27.54 -1.05 -27.51 -1.05 -27.49 -1.06 -27.46 -1.06 -27.44 -1.07 -27.41 -1.07 -27.39 -1.07 -27.36 -1.08 -27.34 -1.08 -27.31 -1.09 -27.29 -1.09 -27.26 -1.09 -27.24 -1.1 -27.22 -1.1 -27.19 -1.11 -27.17 -1.11 -27.14 -1.11 -27.12 -1.12 -27.09 -1.12 -27.07 -1.12 -27.04 -1.13 -27.02 -1.13 -26.99 -1.14 -26.97 -1.14 -26.94 -1.14 -26.92 -1.15 -26.89 -1.15 -26.87 -1.15 -26.84 -1.16 -26.82 -1.16 -26.79 -1.16 -26.77 -1.17 -26.74 -1.17 -26.72 -1.17 -26.69 -1.18 -26.67 -1.18 -26.64 -1.18 -26.62 -1.19 -26.59 -1.19 -26.57 -1.19 -26.54 -1.2 -26.52 -1.2 -26.49 -1.2 -26.47 -1.21 -26.44 -1.21 -26.42 -1.21 -26.39 -1.21 -26.37 -1.22 -26.34 -1.22 -26.32 -1.22 -26.29 -1.23 -26.27 -1.23 -26.25 -1.23 -26.22 -1.24 -26.2 -1.24 -26.17 -1.24 -26.15 -1.24 -26.12 -1.25 -26.1 -1.25 -26.07 -1.25 -26.05 -1.25 -26.02 -1.26 -26 -1.26 -25.97 -1.26 -25.95 -1.27 -25.92 -1.27 -25.9 -1.27 -25.87 -1.27 -25.85 -1.28 -25.82 -1.28 -25.8 -1.28 -25.77 -1.28 -25.75 -1.29 -25.72 -1.29 -25.7 -1.29 -25.67 -1.29 -25.65 -1.3 -25.62 -1.3 -25.6 -1.3 -25.57 -1.3 -25.55 -1.31 -25.52 -1.31 -25.5 -1.31 -25.47 -1.31 -25.45 -1.31 -25.42 -1.32 -25.4 -1.32 -25.37 -1.32 -25.35 -1.32 -25.32 -1.33 -25.3 -1.33 -25.27 -1.33 -25.25 -1.33 -25.23 -1.33 -25.2 -1.34 -25.18 -1.34 -25.15 -1.34 -25.13 -1.34 -25.1 -1.34 -25.08 -1.35 -25.05 -1.35 -25.03 -1.35 -25 -1.35 -24.98 -1.35 -24.95 -1.36 -24.93 -1.36 -24.9 -1.36 -24.88 -1.36 -24.85 -1.36 -24.83 -1.37 -24.8 -1.37 -24.78 -1.37 -24.75 -1.37 -24.73 -1.37 -24.7 -1.37 -24.68 -1.38 -24.65 -1.38 -24.63 -1.38 -24.6 -1.38 -24.58 -1.38 -24.55 -1.39 -24.53 -1.39 -24.5 -1.39 -24.48 -1.39 -24.45 -1.39 -24.43 -1.39 -24.4 -1.4 -24.38 -1.4 -24.35 -1.4 -24.33 -1.4 -24.3 -1.4 -24.28 -1.4 -24.25 -1.4 -24.23 -1.41 -24.21 -1.41 -24.18 -1.41 -24.16 -1.41 -24.13 -1.41 -24.11 -1.41 -24.08 -1.42 -24.06 -1.42 -24.03 -1.42 -24.01 -1.42 -23.98 -1.42 -23.96 -1.42 -23.93 -1.42 -23.91 -1.43 -23.88 -1.43 -23.86 -1.43 -23.83 -1.43 -23.81 -1.43 -23.78 -1.43 -23.76 -1.43 -23.73 -1.44 -23.71 -1.44 -23.68 -1.44 -23.66 -1.44 -23.63 -1.44 -23.61 -1.44 -23.58 -1.44 -23.56 -1.44 -23.53 -1.45 -23.51 -1.45 -23.48 -1.45 -23.46 -1.45 -23.43 -1.45 -23.41 -1.45 -23.38 -1.45 -23.36 -1.45 -23.33 -1.45 -23.31 -1.46 -23.28 -1.46 -23.26 -1.46 -23.23 -1.46 -23.21 -1.46 -23.19 -1.46 -23.16 -1.46 -23.14 -1.46 -23.11 -1.46 -23.09 -1.47 -23.06 -1.47 -23.04 -1.47 -23.01 -1.47 -22.99 -1.47 -22.96 -1.47 -22.94 -1.47 -22.91 -1.47 -22.89 -1.47 -22.86 -1.48 -22.84 -1.48 -22.81 -1.48 -22.79 -1.48 -22.76 -1.48 -22.74 -1.48 -22.71 -1.48 -22.69 -1.48 -22.66 -1.48 -22.64 -1.48 -22.61 -1.48 -22.59 -1.49 -22.56 -1.49 -22.54 -1.49 -22.51 -1.49 -22.49 -1.49 -22.46 -1.49 -22.44 -1.49 -22.41 -1.49 -22.39 -1.49 -22.36 -1.49 -22.34 -1.49 -22.31 -1.5 -22.29 -1.5 -22.26 -1.5 -22.24 -1.5 -22.21 -1.5 -22.19 -1.5 -22.17 -1.5 -22.14 -1.5 -22.12 -1.5 -22.09 -1.5 -22.07 -1.5 -22.04 -1.5 -22.02 -1.5 -21.99 -1.51 -21.97 -1.51 -21.94 -1.51 -21.92 -1.51 -21.89 -1.51 -21.87 -1.51 -21.84 -1.51 -21.82 -1.51 -21.79 -1.51 -21.77 -1.51 -21.74 -1.51 -21.72 -1.51 -21.69 -1.51 -21.67 -1.51 -21.64 -1.52 -21.62 -1.52 -21.59 -1.52 -21.57 -1.52 -21.54 -1.52 -21.52 -1.52 -21.49 -1.52 -21.47 -1.52 -21.44 -1.52 -21.42 -1.52 -21.39 -1.52 -21.37 -1.52 -21.34 -1.52 -21.32 -1.52 -21.29 -1.52 -21.27 -1.52 -21.24 -1.53 -21.22 -1.53 -21.2 -1.53 -21.17 -1.53 -21.15 -1.53 -21.12 -1.53 -21.1 -1.53 -21.07 -1.53 -21.05 -1.53 -21.02 -1.53 -21 -1.53 -20.97 -1.53 -20.95 -1.53 -20.92 -1.53 -20.9 -1.53 -20.87 -1.53 -20.85 -1.53 -20.82 -1.53 -20.8 -1.53 -20.77 -1.54 -20.75 -1.54 -20.72 -1.54 -20.7 -1.54 -20.67 -1.54 -20.65 -1.54 -20.62 -1.54 -20.6 -1.54 -20.57 -1.54 -20.55 -1.54 -20.52 -1.54 -20.5 -1.54 -20.47 -1.54 -20.45 -1.54 -20.42 -1.54 -20.4 -1.54 -20.37 -1.54 -20.35 -1.54 -20.32 -1.54 -20.3 -1.54 -20.27 -1.54 -20.25 -1.54 -20.22 -1.54 -20.2 -1.55 -20.18 -1.55 -20.15 -1.55 -20.13 -1.55 -20.1 -1.55 -20.08 -1.55 -20.05 -1.55 -20.03 -1.55 -20 -1.55 -19.98 -1.55 -19.95 -1.55 -19.93 -1.55 -19.9 -1.55 -19.88 -1.55 -19.85 -1.55 -19.83 -1.55 -19.8 -1.55 -19.78 -1.55 -19.75 -1.55 -19.73 -1.55 -19.7 -1.55 -19.68 -1.55 -19.65 -1.55 -19.63 -1.55 -19.6 -1.55 -19.58 -1.55 -19.55 -1.55 -19.53 -1.55 -19.5 -1.55 -19.48 -1.56 -19.45 -1.56 -19.43 -1.56 -19.4 -1.56 -19.38 -1.56 -19.35 -1.56 -19.33 -1.56 -19.3 -1.56 -19.28 -1.56 -19.25 -1.56 -19.23 -1.56 -19.2 -1.56 -19.18 -1.56 -19.16 -1.56 -19.13 -1.56 -19.11 -1.56 -19.08 -1.56 -19.06 -1.56 -19.03 -1.56 -19.01 -1.56 -18.98 -1.56 -18.96 -1.56 -18.93 -1.56 -18.91 -1.56 -18.88 -1.56 -18.86 -1.56 -18.83 -1.56 -18.81 -1.56 -18.78 -1.56 -18.76 -1.56 -18.73 -1.56 -18.71 -1.56 -18.68 -1.56 -18.66 -1.56 -18.63 -1.56 -18.61 -1.56 -18.58 -1.56 -18.56 -1.56 -18.53 -1.56 -18.51 -1.56 -18.48 -1.56 -18.46 -1.56 -18.43 -1.56 -18.41 -1.56 -18.38 -1.56 -18.36 -1.57 -18.33 -1.57 -18.31 -1.57 -18.28 -1.57 -18.26 -1.57 -18.23 -1.57 -18.21 -1.57 -18.18 -1.57 -18.16 -1.57 -18.14 -1.57 -18.11 -1.57 -18.09 -1.57 -18.06 -1.57 -18.04 -1.57 -18.01 -1.57 -17.99 -1.57 -17.96 -1.57 -17.94 -1.57 -17.91 -1.57 -17.89 -1.57 -17.86 -1.57 -17.84 -1.57 -17.81 -1.57 -17.79 -1.57 -17.76 -1.57 -17.74 -1.57 -17.71 -1.57 -17.69 -1.57 -17.66 -1.57 -17.64 -1.57 -17.61 -1.57 -17.59 -1.57 -17.56 -1.57 -17.54 -1.57 -17.51 -1.57 -17.49 -1.57 -17.46 -1.57 -17.44 -1.57 -17.41 -1.57 -17.39 -1.57 -17.36 -1.57 -17.34 -1.57 -17.31 -1.57 -17.29 -1.57 -17.26 -1.57 -17.24 -1.57 -17.21 -1.57 -17.19 -1.57 -17.16 -1.57 -17.14 -1.57 -17.12 -1.57 -17.09 -1.57 -17.07 -1.57 -17.04 -1.57 -17.02 -1.57 -16.99 -1.57 -16.97 -1.57 -16.94 -1.57 -16.92 -1.57 -16.89 -1.57 -16.87 -1.57 -16.84 -1.57 -16.82 -1.57 -16.79 -1.57 -16.77 -1.57 -16.74 -1.57 -16.72 -1.57 -16.69 -1.57 -16.67 -1.57 -16.64 -1.57 -16.62 -1.57 -16.59 -1.57 -16.57 -1.57 -16.54 -1.57 -16.52 -1.57 -16.49 -1.57 -16.47 -1.57 -16.44 -1.57 -16.42 -1.57 -16.39 -1.57 -16.37 -1.57 -16.34 -1.57 -16.32 -1.57 -16.29 -1.57 -16.27 -1.57 -16.24 -1.57 -16.22 -1.57 -16.19 -1.57 -16.17 -1.57 -16.15 -1.57 -16.12 -1.57 -16.1 -1.57 -16.07 -1.57 -16.05 -1.57 -16.02 -1.57 -16 -1.57 -15.97 -1.57 -15.95 -1.57 -15.92 -1.57 -15.9 -1.57 -15.87 -1.57 -15.85 -1.57 -15.82 -1.57 -15.8 -1.57 -15.77 -1.57 -15.75 -1.57 -15.72 -1.57 -15.7 -1.57 -15.67 -1.57 -15.65 -1.57 -15.62 -1.57 -15.6 -1.57 -15.57 -1.57 -15.55 -1.57 -15.52 -1.57 -15.5 -1.57 -15.47 -1.57 -15.45 -1.57 -15.42 -1.57 -15.4 -1.57 -15.37 -1.57 -15.35 -1.57 -15.32 -1.57 -15.3 -1.57 -15.27 -1.57 -15.25 -1.57 -15.22 -1.57 -15.2 -1.57 -15.17 -1.57 -15.15 -1.57 -15.13 -1.57 -15.1 -1.57 -15.08 -1.57 -15.05 -1.57 -15.03 -1.57 -15 -1.57 -14.98 -1.57 -14.95 -1.57 -14.93 -1.57 -14.9 -1.57 -14.88 -1.57 -14.85 -1.57 -14.83 -1.57 -14.8 -1.57 -14.78 -1.57 -14.75 -1.57 -14.73 -1.57 -14.7 -1.57 -14.68 -1.57 -14.65 -1.57 -14.63 -1.57 -14.6 -1.57 -14.58 -1.57 -14.55 -1.57 -14.53 -1.57 -14.5 -1.57 -14.48 -1.57 -14.45 -1.57 -14.43 -1.57 -14.4 -1.57 -14.38 -1.57 -14.35 -1.57 -14.33 -1.57 -14.3 -1.57 -14.28 -1.57 -14.25 -1.57 -14.23 -1.57 -14.2 -1.57 -14.18 -1.57 -14.15 -1.57 -14.13 -1.57 -14.11 -1.57 -14.08 -1.57 -14.06 -1.57 -14.03 -1.57 -14.01 -1.57 -13.98 -1.57 -13.96 -1.57 -13.93 -1.57 -13.91 -1.57 -13.88 -1.57 -13.86 -1.57 -13.83 -1.57 -13.81 -1.57 -13.78 -1.57 -13.76 -1.57 -13.73 -1.57 -13.71 -1.57 -13.68 -1.57 -13.66 -1.57 -13.63 -1.57 -13.61 -1.57 -13.58 -1.57 -13.56 -1.57 -13.53 -1.57 -13.51 -1.57 -13.48 -1.57 -13.46 -1.57 -13.43 -1.57 -13.41 -1.57 -13.38 -1.57 -13.36 -1.57 -13.33 -1.57 -13.31 -1.57 -13.28 -1.57 -13.26 -1.57 -13.23 -1.57 -13.21 -1.57 -13.18 -1.57 -13.16 -1.57 -13.13 -1.57 -13.11 -1.57 -13.09 -1.57 -13.06 -1.57 -13.04 -1.57 -13.01 -1.57 -12.99 -1.57 -12.96 -1.57 -12.94 -1.57 -12.91 -1.57 -12.89 -1.57 -12.86 -1.57 -12.84 -1.57 -12.81 -1.57 -12.79 -1.57 -12.76 -1.57 -12.74 -1.57 -12.71 -1.57 -12.69 -1.57 -12.66 -1.57 -12.64 -1.57 -12.61 -1.57 -12.59 -1.57 -12.56 -1.57 -12.54 -1.57 -12.51 -1.57 -12.49 -1.57 -12.46 -1.57 -12.44 -1.57 -12.41 -1.57 -12.39 -1.57 -12.36 -1.57 -12.34 -1.57 -12.31 -1.57 -12.29 -1.57 -12.26 -1.57 -12.24 -1.57 -12.21 -1.57 -12.19 -1.57 -12.16 -1.56 -12.14 -1.56 -12.12 -1.56 -12.09 -1.56 -12.07 -1.56 -12.04 -1.56 -12.02 -1.56 -11.99 -1.56 -11.97 -1.56 -11.94 -1.56 -11.92 -1.56 -11.89 -1.56 -11.87 -1.56 -11.84 -1.56 -11.82 -1.56 -11.79 -1.56 -11.77 -1.56 -11.74 -1.56 -11.72 -1.56 -11.69 -1.56 -11.67 -1.56 -11.64 -1.56 -11.62 -1.56 -11.59 -1.56 -11.57 -1.56 -11.54 -1.56 -11.52 -1.56 -11.49 -1.56 -11.47 -1.56 -11.44 -1.56 -11.42 -1.56 -11.39 -1.56 -11.37 -1.56 -11.34 -1.56 -11.32 -1.56 -11.29 -1.56 -11.27 -1.56 -11.24 -1.56 -11.22 -1.56 -11.19 -1.56 -11.17 -1.56 -11.14 -1.56 -11.12 -1.56 -11.1 -1.56 -11.07 -1.56 -11.05 -1.56 -11.02 -1.56 -11 -1.56 -10.97 -1.56 -10.95 -1.56 -10.92 -1.56 -10.9 -1.56 -10.87 -1.56 -10.85 -1.56 -10.82 -1.56 -10.8 -1.56 -10.77 -1.56 -10.75 -1.56 -10.72 -1.56 -10.7 -1.56 -10.67 -1.56 -10.65 -1.56 -10.62 -1.56 -10.6 -1.56 -10.57 -1.56 -10.55 -1.56 -10.52 -1.56 -10.5 -1.56 -10.47 -1.56 -10.45 -1.56 -10.42 -1.56 -10.4 -1.56 -10.37 -1.56 -10.35 -1.56 -10.32 -1.56 -10.3 -1.56 -10.27 -1.56 -10.25 -1.56 -10.22 -1.56 -10.2 -1.56 -10.17 -1.56 -10.15 -1.56 -10.12 -1.56 -10.1 -1.56 -10.08 -1.56 -10.05 -1.56 -10.03 -1.56 -10 -1.56 -9.98 -1.56 -9.95 -1.56 -9.93 -1.56 -9.9 -1.56 -9.88 -1.55 -9.85 -1.55 -9.83 -1.55 -9.8 -1.55 -9.78 -1.55 -9.75 -1.55 -9.73 -1.55 -9.7 -1.55 -9.68 -1.55 -9.65 -1.55 -9.63 -1.55 -9.6 -1.55 -9.58 -1.55 -9.55 -1.55 -9.53 -1.55 -9.5 -1.55 -9.48 -1.55 -9.45 -1.55 -9.43 -1.55 -9.4 -1.55 -9.38 -1.55 -9.35 -1.55 -9.33 -1.55 -9.3 -1.55 -9.28 -1.55 -9.25 -1.55 -9.23 -1.55 -9.2 -1.55 -9.18 -1.55 -9.15 -1.55 -9.13 -1.55 -9.1 -1.55 -9.08 -1.55 -9.06 -1.55 -9.03 -1.55 -9.01 -1.55 -8.98 -1.55 -8.96 -1.55 -8.93 -1.55 -8.91 -1.55 -8.88 -1.55 -8.86 -1.55 -8.83 -1.55 -8.81 -1.55 -8.78 -1.55 -8.76 -1.55 -8.73 -1.55 -8.71 -1.55 -8.68 -1.55 -8.66 -1.55 -8.63 -1.55 -8.61 -1.55 -8.58 -1.55 -8.56 -1.55 -8.53 -1.55 -8.51 -1.55 -8.48 -1.55 -8.46 -1.55 -8.43 -1.55 -8.41 -1.55 -8.38 -1.55 -8.36 -1.55 -8.33 -1.55 -8.31 -1.55 -8.28 -1.55 -8.26 -1.55 -8.23 -1.55 -8.21 -1.55 -8.18 -1.55 -8.16 -1.55 -8.13 -1.55 -8.11 -1.55 -8.08 -1.55 -8.06 -1.55 -8.04 -1.55 -8.01 -1.55 -7.99 -1.55 -7.96 -1.55 -7.94 -1.55 -7.91 -1.55 -7.89 -1.55 -7.86 -1.55 -7.84 -1.54 -7.81 -1.54 -7.79 -1.54 -7.76 -1.54 -7.74 -1.54 -7.71 -1.54 -7.69 -1.54 -7.66 -1.54 -7.64 -1.54 -7.61 -1.54 -7.59 -1.54 -7.56 -1.54 -7.54 -1.54 -7.51 -1.54 -7.49 -1.54 -7.46 -1.54 -7.44 -1.54 -7.41 -1.54 -7.39 -1.54 -7.36 -1.54 -7.34 -1.54 -7.31 -1.54 -7.29 -1.54 -7.26 -1.54 -7.24 -1.54 -7.21 -1.54 -7.19 -1.54 -7.16 -1.54 -7.14 -1.54 -7.11 -1.54 -7.09 -1.54 -7.07 -1.54 -7.04 -1.54 -7.02 -1.54 -6.99 -1.54 -6.97 -1.54 -6.94 -1.54 -6.92 -1.54 -6.89 -1.54 -6.87 -1.54 -6.84 -1.54 -6.82 -1.54 -6.79 -1.54 -6.77 -1.54 -6.74 -1.54 -6.72 -1.54 -6.69 -1.54 -6.67 -1.54 -6.64 -1.54 -6.62 -1.54 -6.59 -1.54 -6.57 -1.54 -6.54 -1.54 -6.52 -1.54 -6.49 -1.54 -6.47 -1.54 -6.44 -1.54 -6.42 -1.54 -6.39 -1.54 -6.37 -1.54 -6.34 -1.54 -6.32 -1.54 -6.29 -1.54 -6.27 -1.54 -6.24 -1.54 -6.22 -1.54 -6.19 -1.54 -6.17 -1.54 -6.14 -1.54 -6.12 -1.54 -6.09 -1.54 -6.07 -1.54 -6.05 -1.54 -6.02 -1.54 -6 -1.54 -5.97 -1.54 -5.95 -1.54 -5.92 -1.54 -5.9 -1.54 -5.87 -1.54 -5.85 -1.54 -5.82 -1.53 -5.8 -1.53 -5.77 -1.53 -5.75 -1.53 -5.72 -1.53 -5.7 -1.53 -5.67 -1.53 -5.65 -1.53 -5.62 -1.53 -5.6 -1.53 -5.57 -1.53 -5.55 -1.53 -5.52 -1.53 -5.5 -1.53 -5.47 -1.53 -5.45 -1.53 -5.42 -1.53 -5.4 -1.53 -5.37 -1.53 -5.35 -1.53 -5.32 -1.53 -5.3 -1.53 -5.27 -1.53 -5.25 -1.53 -5.22 -1.53 -5.2 -1.53 -5.17 -1.53 -5.15 -1.53 -5.12 -1.53 -5.1 -1.53 -5.07 -1.53 -5.05 -1.53 -5.03 -1.53 -5 -1.53 -4.98 -1.53 -4.95 -1.53 -4.93 -1.53 -4.9 -1.53 -4.88 -1.53 -4.85 -1.53 -4.83 -1.53 -4.8 -1.53 -4.78 -1.53 -4.75 -1.53 -4.73 -1.53 -4.7 -1.53 -4.68 -1.53 -4.65 -1.53 -4.63 -1.53 -4.6 -1.53 -4.58 -1.53 -4.55 -1.53 -4.53 -1.53 -4.5 -1.53 -4.48 -1.53 -4.45 -1.53 -4.43 -1.53 -4.4 -1.53 -4.38 -1.53 -4.35 -1.53 -4.33 -1.53 -4.3 -1.53 -4.28 -1.53 -4.25 -1.53 -4.23 -1.53 -4.2 -1.53 -4.18 -1.53 -4.15 -1.53 -4.13 -1.53 -4.1 -1.53 -4.08 -1.53 -4.05 -1.53 -4.03 -1.53 -4.01 -1.53 -3.98 -1.53 -3.96 -1.53 -3.93 -1.53 -3.91 -1.53 -3.88 -1.53 -3.86 -1.53 -3.83 -1.53 -3.81 -1.53 -3.78 -1.53 -3.76 -1.53 -3.73 -1.53 -3.71 -1.53 -3.68 -1.52 -3.66 -1.52 -3.63 -1.52 -3.61 -1.52 -3.58 -1.52 -3.56 -1.52 -3.53 -1.52 -3.51 -1.52 -3.48 -1.52 -3.46 -1.52 -3.43 -1.52 -3.41 -1.52 -3.38 -1.52 -3.36 -1.52 -3.33 -1.52 -3.31 -1.52 -3.28 -1.52 -3.26 -1.52 -3.23 -1.52 -3.21 -1.52 -3.18 -1.52 -3.16 -1.52 -3.13 -1.52 -3.11 -1.52 -3.08 -1.52 -3.06 -1.52 -3.03 -1.52 -3.01 -1.52 -2.99 -1.52 -2.96 -1.52 -2.94 -1.52 -2.91 -1.52 -2.89 -1.52 -2.86 -1.52 -2.84 -1.52 -2.81 -1.52 -2.79 -1.52 -2.76 -1.52 -2.74 -1.52 -2.71 -1.52 -2.69 -1.52 -2.66 -1.52 -2.64 -1.52 -2.61 -1.52 -2.59 -1.52 -2.56 -1.52 -2.54 -1.52 -2.51 -1.52 -2.49 -1.52 -2.46 -1.52 -2.44 -1.52 -2.41 -1.52 -2.39 -1.52 -2.36 -1.52 -2.34 -1.52 -2.31 -1.52 -2.29 -1.52 -2.26 -1.52 -2.24 -1.52 -2.21 -1.52 -2.19 -1.52 -2.16 -1.52 -2.14 -1.52 -2.11 -1.52 -2.09 -1.52 -2.06 -1.52 -2.04 -1.52 -2.02 -1.52 -1.99 -1.52 -1.97 -1.52 -1.94 -1.52 -1.92 -1.52 -1.89 -1.52 -1.87 -1.52 -1.84 -1.52 -1.82 -1.52 -1.79 -1.52 -1.77 -1.52 -1.74 -1.52 -1.72 -1.52 -1.69 -1.52 -1.67 -1.52 -1.64 -1.52 -1.62 -1.52 -1.59 -1.52 -1.57 -1.52 -1.54 -1.52 -1.52 -1.52 -1.49 -1.52 -1.47 -1.52 -1.44 -1.52 -1.42 -1.52 -1.39 -1.52 -1.37 -1.52 -1.34 -1.52 -1.32 -1.52 -1.29 -1.51 -1.27 -1.51 -1.24 -1.51 -1.22 -1.51 -1.19 -1.51 -1.17 -1.51 -1.14 -1.51 -1.12 -1.51 -1.09 -1.51 -1.07 -1.51 -1.04 -1.51 -1.02 -1.51 -1 -1.51 -0.97 -1.51 -0.95 -1.51 -0.92 -1.51 -0.9 -1.51 -0.87 -1.51 -0.85 -1.51 -0.82 -1.51 -0.8 -1.51 -0.77 -1.51 -0.75 -1.51 -0.72 -1.51 -0.7 -1.51 -0.67 -1.51 -0.65 -1.51 -0.62 -1.51 -0.6 -1.51 -0.57 -1.51 -0.55 -1.51 -0.52 -1.51 -0.5 -1.51 -0.47 -1.51 -0.45 -1.51 -0.42 -1.51 -0.4 -1.51 -0.37 -1.51 -0.35 -1.51 -0.32 -1.51 -0.3 -1.51 -0.27 -1.51 -0.25 -1.51 -0.22 -1.51 -0.2 -1.51 -0.17 -1.51 -0.15 -1.51 -0.12 -1.51 -0.1 -1.51 -0.07 -1.51 -0.05 -1.51 -0.02 -1.51 0 -1.51 0.02 -1.51 0.05 -1.51 0.07 -1.51 0.1 -1.51 0.12 -1.51 0.15 -1.51 0.17 -1.51 0.2 -1.51 0.22 -1.51 0.25 -1.51 0.27 -1.51 0.3 -1.51 0.32 -1.51 0.35 -1.51 0.37 -1.51 0.4 -1.51 0.42 -1.51 0.45 -1.51 0.47 -1.51 0.5 -1.51 0.52 -1.51 0.55 -1.51 0.57 -1.51 0.6 -1.51 0.62 -1.51 0.65 -1.51 0.67 -1.51 0.7 -1.51 0.72 -1.51 0.75 -1.51 0.77 -1.51 0.8 -1.51 0.82 -1.51 0.85 -1.51 0.87 -1.51 0.9 -1.51 0.92 -1.51 0.95 -1.51 0.97 -1.51 1 -1.51 1.02 -1.51 1.04 -1.51 1.07 -1.51 1.09 -1.51 1.12 -1.51 1.14 -1.51 1.17 -1.51 1.19 -1.51 1.22 -1.51 1.24 -1.51 1.27 -1.51 1.29 -1.51 1.32 -1.51 1.34 -1.51 1.37 -1.51 1.39 -1.51 1.42 -1.51 1.44 -1.51 1.47 -1.51 1.49 -1.51 1.52 -1.51 1.54 -1.5 1.57 -1.5 1.59 -1.5 1.62 -1.5 1.64 -1.5 1.67 -1.5 1.69 -1.5 1.72 -1.5 1.74 -1.5 1.77 -1.5 1.79 -1.5 1.82 -1.5 1.84 -1.5 1.87 -1.5 1.89 -1.5 1.92 -1.5 1.94 -1.5 1.97 -1.5 1.99 -1.5 2.02 -1.5 2.04 -1.5 2.06 -1.5 2.09 -1.5 2.11 -1.5 2.14 -1.5 2.16 -1.5 2.19 -1.5 2.21 -1.5 2.24 -1.5 2.26 -1.5 2.29 -1.5 2.31 -1.5 2.34 -1.5 2.36 -1.5 2.39 -1.5 2.41 -1.5 2.44 -1.5 2.46 -1.5 2.49 -1.5 2.51 -1.5 2.54 -1.5 2.56 -1.5 2.59 -1.5 2.61 -1.5 2.64 -1.5 2.66 -1.5 2.69 -1.5 2.71 -1.5 2.74 -1.5 2.76 -1.5 2.79 -1.5 2.81 -1.5 2.84 -1.5 2.86 -1.5 2.89 -1.5 2.91 -1.5 2.94 -1.5 2.96 -1.5 2.99 -1.5 3.01 -1.5 3.03 -1.5 3.06 -1.5 3.08 -1.5 3.11 -1.5 3.13 -1.5 3.16 -1.5 3.18 -1.5 3.21 -1.5 3.23 -1.5 3.26 -1.5 3.28 -1.5 3.31 -1.5 3.33 -1.5 3.36 -1.5 3.38 -1.5 3.41 -1.5 3.43 -1.5 3.46 -1.5 3.48 -1.5 3.51 -1.5 3.53 -1.5 3.56 -1.5 3.58 -1.5 3.61 -1.5 3.63 -1.5 3.66 -1.5 3.68 -1.5 3.71 -1.5 3.73 -1.5 3.76 -1.5 3.78 -1.5 3.81 -1.5 3.83 -1.5 3.86 -1.5 3.88 -1.5 3.91 -1.5 3.93 -1.5 3.96 -1.5 3.98 -1.5 4.01 -1.5 4.03 -1.5 4.05 -1.5 4.08 -1.5 4.1 -1.5 4.13 -1.5 4.15 -1.5 4.18 -1.5 4.2 -1.5 4.23 -1.5 4.25 -1.5 4.28 -1.5 4.3 -1.5 4.33 -1.5 4.35 -1.5 4.38 -1.5 4.4 -1.5 4.43 -1.5 4.45 -1.5 4.48 -1.5 4.5 -1.5 4.53 -1.5 4.55 -1.5 4.58 -1.5 4.6 -1.5 4.63 -1.5 4.65 -1.5 4.68 -1.5 4.7 -1.5 4.73 -1.5 4.75 -1.5 4.78 -1.5 4.8 -1.5 4.83 -1.5 4.85 -1.5 4.88 -1.5 4.9 -1.5 4.93 -1.5 4.95 -1.5 4.98 -1.5 5 -1.5 5.03 -1.5 5.05 -1.5 5.07 -1.5 5.1 -1.5 5.12 -1.5 5.15 -1.5 5.17 -1.5 5.2 -1.5 5.22 -1.5 5.25 -1.5 5.27 -1.5 5.3 -1.5 5.32 -1.5 5.35 -1.5 5.37 -1.5 5.4 -1.49 5.42 -1.49 5.45 -1.49 5.47 -1.49 5.5 -1.49 5.52 -1.49 5.55 -1.49 5.57 -1.49 5.6 -1.49 5.62 -1.49 5.65 -1.49 5.67 -1.49 5.7 -1.49 5.72 -1.49 5.75 -1.49 5.77 -1.49 5.8 -1.49 5.82 -1.49 5.85 -1.49 5.87 -1.49 5.9 -1.49 5.92 -1.49 5.95 -1.49 5.97 -1.49 6 -1.49 6.02 -1.49 6.05 -1.49 6.07 -1.49 6.09 -1.49 6.12 -1.49 6.14 -1.49 6.17 -1.49 6.19 -1.49 6.22 -1.49 6.24 -1.49 6.27 -1.49 6.29 -1.49 6.32 -1.49 6.34 -1.49 6.37 -1.49 6.39 -1.49 6.42 -1.49 6.44 -1.49 6.47 -1.49 6.49 -1.49 6.52 -1.49 6.54 -1.49 6.57 -1.49 6.59 -1.49 6.62 -1.49 6.64 -1.49 6.67 -1.49 6.69 -1.49 6.72 -1.49 6.74 -1.49 6.77 -1.49 6.79 -1.49 6.82 -1.49 6.84 -1.49 6.87 -1.49 6.89 -1.49 6.92 -1.49 6.94 -1.49 6.97 -1.49 6.99 -1.49 7.02 -1.49 7.04 -1.49 7.07 -1.49 7.09 -1.49 7.11 -1.49 7.14 -1.49 7.16 -1.49 7.19 -1.49 7.21 -1.49 7.24 -1.49 7.26 -1.49 7.29 -1.49 7.31 -1.49 7.34 -1.49 7.36 -1.49 7.39 -1.49 7.41 -1.49 7.44 -1.49 7.46 -1.49 7.49 -1.49 7.51 -1.49 7.54 -1.49 7.56 -1.49 7.59 -1.49 7.61 -1.49 7.64 -1.49 7.66 -1.49 7.69 -1.49 7.71 -1.49 7.74 -1.49 7.76 -1.49 7.79 -1.49 7.81 -1.49 7.84 -1.49 7.86 -1.49 7.89 -1.49 7.91 -1.49 7.94 -1.49 7.96 -1.49 7.99 -1.49 8.01 -1.49 8.04 -1.49 8.06 -1.49 8.08 -1.49 8.11 -1.49 8.13 -1.49 8.16 -1.49 8.18 -1.49 8.21 -1.49 8.23 -1.49 8.26 -1.49 8.28 -1.49 8.31 -1.49 8.33 -1.49 8.36 -1.49 8.38 -1.49 8.41 -1.49 8.43 -1.49 8.46 -1.49 8.48 -1.49 8.51 -1.49 8.53 -1.49 8.56 -1.49 8.58 -1.49 8.61 -1.49 8.63 -1.49 8.66 -1.49 8.68 -1.49 8.71 -1.49 8.73 -1.49 8.76 -1.49 8.78 -1.49 8.81 -1.49 8.83 -1.49 8.86 -1.49 8.88 -1.49 8.91 -1.49 8.93 -1.49 8.96 -1.49 8.98 -1.49 9.01 -1.49 9.03 -1.49 9.06 -1.49 9.08 -1.49 9.1 -1.49 9.13 -1.49 9.15 -1.49 9.18 -1.49 9.2 -1.49 9.23 -1.49 9.25 -1.49 9.28 -1.49 9.3 -1.49 9.33 -1.49 9.35 -1.49 9.38 -1.49 9.4 -1.49 9.43 -1.49 9.45 -1.49 9.48 -1.49 9.5 -1.49 9.53 -1.49 9.55 -1.49 9.58 -1.49 9.6 -1.49 9.63 -1.49 9.65 -1.49 9.68 -1.49 9.7 -1.49 9.73 -1.49 9.75 -1.49 9.78 -1.49 9.8 -1.49 9.83 -1.49 9.85 -1.49 9.88 -1.49 9.9 -1.49 9.93 -1.49 9.95 -1.49 9.98 -1.49 10 -1.49 10.03 -1.49 10.05 -1.49 10.08 -1.49 10.1 -1.49 10.12 -1.49 10.15 -1.49 10.17 -1.49 10.2 -1.49 10.22 -1.49 10.25 -1.49 10.27 -1.49 10.3 -1.49 10.32 -1.49 10.35 -1.49 10.37 -1.49 10.4 -1.49 10.42 -1.49 10.45 -1.49 10.47 -1.49 10.5 -1.49 10.52 -1.49 10.55 -1.49 10.57 -1.49 10.6 -1.49 10.62 -1.49 10.65 -1.49 10.67 -1.49 10.7 -1.49 10.72 -1.49 10.75 -1.49 10.77 -1.49 10.8 -1.49 10.82 -1.49 10.85 -1.49 10.87 -1.49 10.9 -1.49 10.92 -1.49 10.95 -1.49 10.97 -1.49 11 -1.49 11.02 -1.49 11.05 -1.49 11.07 -1.49 11.1 -1.49 11.12 -1.49 11.14 -1.49 11.17 -1.49 11.19 -1.49 11.22 -1.49 11.24 -1.49 11.27 -1.49 11.29 -1.49 11.32 -1.49 11.34 -1.49 11.37 -1.49 11.39 -1.49 11.42 -1.49 11.44 -1.49 11.47 -1.49 11.49 -1.49 11.52 -1.49 11.54 -1.49 11.57 -1.49 11.59 -1.49 11.62 -1.49 11.64 -1.49 11.67 -1.49 11.69 -1.49 11.72 -1.49 11.74 -1.49 11.77 -1.49 11.79 -1.49 11.82 -1.49 11.84 -1.49 11.87 -1.49 11.89 -1.49 11.92 -1.49 11.94 -1.49 11.97 -1.49 11.99 -1.49 12.02 -1.49 12.04 -1.49 12.07 -1.49 12.09 -1.49 12.12 -1.49 12.14 -1.49 12.16 -1.49 12.19 -1.49 12.21 -1.49 12.24 -1.49 12.26 -1.49 12.29 -1.49 12.31 -1.49 12.34 -1.49 12.36 -1.49 12.39 -1.49 12.41 -1.49 12.44 -1.49 12.46 -1.49 12.49 -1.49 12.51 -1.49 12.54 -1.49 12.56 -1.49 12.59 -1.49 12.61 -1.49 12.64 -1.49 12.66 -1.49 12.69 -1.49 12.71 -1.49 12.74 -1.49 12.76 -1.49 12.79 -1.49 12.81 -1.48 12.84 -1.48 12.86 -1.48 12.89 -1.48 12.91 -1.48 12.94 -1.48 12.96 -1.48 12.99 -1.48 13.01 -1.48 13.04 -1.48 13.06 -1.48 13.09 -1.48 13.11 -1.48 13.13 -1.48 13.16 -1.48 13.18 -1.48 13.21 -1.48 13.23 -1.48 13.26 -1.48 13.28 -1.48 13.31 -1.48 13.33 -1.48 13.36 -1.48 13.38 -1.48 13.41 -1.48 13.43 -1.48 13.46 -1.48 13.48 -1.48 13.51 -1.48 13.53 -1.48 13.56 -1.48 13.58 -1.48 13.61 -1.48 13.63 -1.48 13.66 -1.48 13.68 -1.48 13.71 -1.48 13.73 -1.48 13.76 -1.48 13.78 -1.48 13.81 -1.48 13.83 -1.48 13.86 -1.48 13.88 -1.48 13.91 -1.48 13.93 -1.48 13.96 -1.48 13.98 -1.48 14.01 -1.48 14.03 -1.48 14.06 -1.48 14.08 -1.48 14.11 -1.48 14.13 -1.48 14.15 -1.48 14.18 -1.48 14.2 -1.48 14.23 -1.48 14.25 -1.48 14.28 -1.48 14.3 -1.48 14.33 -1.48 14.35 -1.48 14.38 -1.48 14.4 -1.48 14.43 -1.48 14.45 -1.48 14.48 -1.48 14.5 -1.48 14.53 -1.48 14.55 -1.48 14.58 -1.48 14.6 -1.48 14.63 -1.48 14.65 -1.48 14.68 -1.48 14.7 -1.48 14.73 -1.48 14.75 -1.48 14.78 -1.48 14.8 -1.48 14.83 -1.48 14.85 -1.48 14.88 -1.48 14.9 -1.48 14.93 -1.48 14.95 -1.48 14.98 -1.48 15 -1.48 15.03 -1.48 15.05 -1.48 15.08 -1.48 15.1 -1.48 15.13 -1.48 15.15 -1.48 15.17 -1.48 15.2 -1.48 15.22 -1.48 15.25 -1.48 15.27 -1.48 15.3 -1.48 15.32 -1.48 15.35 -1.48 15.37 -1.48 15.4 -1.48 15.42 -1.48 15.45 -1.48 15.47 -1.48 15.5 -1.48 15.52 -1.48 15.55 -1.48 15.57 -1.48 15.6 -1.48 15.62 -1.48 15.65 -1.48 15.67 -1.48 15.7 -1.48 15.72 -1.48 15.75 -1.48 15.77 -1.48 15.8 -1.48 15.82 -1.48 15.85 -1.48 15.87 -1.48 15.9 -1.48 15.92 -1.48 15.95 -1.48 15.97 -1.48 16 -1.48 16.02 -1.48 16.05 -1.48 16.07 -1.48 16.1 -1.48 16.12 -1.48 16.15 -1.48 16.17 -1.48 16.19 -1.48 16.22 -1.48 16.24 -1.48 16.27 -1.48 16.29 -1.48 16.32 -1.48 16.34 -1.48 16.37 -1.48 16.39 -1.48 16.42 -1.48 16.44 -1.48 16.47 -1.48 16.49 -1.48 16.52 -1.48 16.54 -1.48 16.57 -1.48 16.59 -1.48 16.62 -1.48 16.64 -1.48 16.67 -1.48 16.69 -1.48 16.72 -1.48 16.74 -1.48 16.77 -1.48 16.79 -1.48 16.82 -1.48 16.84 -1.48 16.87 -1.48 16.89 -1.48 16.92 -1.48 16.94 -1.48 16.97 -1.48 16.99 -1.48 17.02 -1.48 17.04 -1.48 17.07 -1.48 17.09 -1.48 17.12 -1.48 17.14 -1.48 17.16 -1.48 17.19 -1.48 17.21 -1.48 17.24 -1.48 17.26 -1.48 17.29 -1.48 17.31 -1.48 17.34 -1.48 17.36 -1.48 17.39 -1.48 17.41 -1.48 17.44 -1.48 17.46 -1.48 17.49 -1.48 17.51 -1.48 17.54 -1.48 17.56 -1.48 17.59 -1.48 17.61 -1.48 17.64 -1.48 17.66 -1.48 17.69 -1.48 17.71 -1.48 17.74 -1.48 17.76 -1.48 17.79 -1.48 17.81 -1.48 17.84 -1.48 17.86 -1.48 17.89 -1.48 17.91 -1.48 17.94 -1.48 17.96 -1.48 17.99 -1.48 18.01 -1.48 18.04 -1.48 18.06 -1.48 18.09 -1.48 18.11 -1.48 18.14 -1.48 18.16 -1.48 18.18 -1.48 18.21 -1.48 18.23 -1.48 18.26 -1.48 18.28 -1.48 18.31 -1.48 18.33 -1.48 18.36 -1.48 18.38 -1.48 18.41 -1.48 18.43 -1.48 18.46 -1.48 18.48 -1.48 18.51 -1.48 18.53 -1.48 18.56 -1.48 18.58 -1.48 18.61 -1.48 18.63 -1.48 18.66 -1.48 18.68 -1.48 18.71 -1.48 18.73 -1.48 18.76 -1.48 18.78 -1.48 18.81 -1.48 18.83 -1.48 18.86 -1.48 18.88 -1.48 18.91 -1.48 18.93 -1.48 18.96 -1.48 18.98 -1.48 19.01 -1.48 19.03 -1.48 19.06 -1.48 19.08 -1.48 19.11 -1.48 19.13 -1.48 19.16 -1.48 19.18 -1.48 19.2 -1.48 19.23 -1.48 19.25 -1.48 19.28 -1.48 19.3 -1.48 19.33 -1.48 19.35 -1.48 19.38 -1.48 19.4 -1.48 19.43 -1.48 19.45 -1.48 19.48 -1.48 19.5 -1.48 19.53 -1.48 19.55 -1.48 19.58 -1.48 19.6 -1.48 19.63 -1.48 19.65 -1.48 19.68 -1.48 19.7 -1.48 19.73 -1.48 19.75 -1.48 19.78 -1.48 19.8 -1.48 19.83 -1.48 19.85 -1.48 19.88 -1.48 19.9 -1.48 19.93 -1.48 19.95 -1.48 19.98 -1.48 20 -1.48 20.03 -1.48 20.05 -1.48 20.08 -1.48 20.1 -1.48 20.13 -1.48 20.15 -1.48 20.18 -1.48 20.2 -1.48 20.22 -1.48 20.25 -1.48 20.27 -1.48 20.3 -1.48 20.32 -1.48 20.35 -1.48 20.37 -1.48 20.4 -1.48 20.42 -1.48 20.45 -1.48 20.47 -1.48 20.5 -1.48 20.52 -1.48 20.55 -1.48 20.57 -1.48 20.6 -1.48 20.62 -1.48 20.65 -1.48 20.67 -1.48 20.7 -1.48 20.72 -1.48 20.75 -1.48 20.77 -1.48 20.8 -1.48 20.82 -1.48 20.85 -1.48 20.87 -1.48 20.9 -1.48 20.92 -1.48 20.95 -1.48 20.97 -1.48 21 -1.48 21.02 -1.48 21.05 -1.48 21.07 -1.48 21.1 -1.48 21.12 -1.48 21.15 -1.48 21.17 -1.48 21.2 -1.48 21.22 -1.48 21.24 -1.48 21.27 -1.48 21.29 -1.48 21.32 -1.48 21.34 -1.48 21.37 -1.48 21.39 -1.48 21.42 -1.48 21.44 -1.48 21.47 -1.48 21.49 -1.48 21.52 -1.48 21.54 -1.48 21.57 -1.48 21.59 -1.48 21.62 -1.48 21.64 -1.48 21.67 -1.48 21.69 -1.48 21.72 -1.48 21.74 -1.48 21.77 -1.48 21.79 -1.48 21.82 -1.48 21.84 -1.48 21.87 -1.48 21.89 -1.48 21.92 -1.48 21.94 -1.48 21.97 -1.48 21.99 -1.48 22.02 -1.48 22.04 -1.48 22.07 -1.48 22.09 -1.48 22.12 -1.48 22.14 -1.48 22.17 -1.48 22.19 -1.48 22.21 -1.48 22.24 -1.48 22.26 -1.48 22.29 -1.48 22.31 -1.48 22.34 -1.48 22.36 -1.48 22.39 -1.48 22.41 -1.48 22.44 -1.48 22.46 -1.48 22.49 -1.48 22.51 -1.48 22.54 -1.48 22.56 -1.48 22.59 -1.48 22.61 -1.48 22.64 -1.48 22.66 -1.48 22.69 -1.48 22.71 -1.48 22.74 -1.48 22.76 -1.48 22.79 -1.48 22.81 -1.48 22.84 -1.48 22.86 -1.48 22.89 -1.48 22.91 -1.48 22.94 -1.48 22.96 -1.48 22.99 -1.48 23.01 -1.48 23.04 -1.48 23.06 -1.48 23.09 -1.48 23.11 -1.48 23.14 -1.48 23.16 -1.48 23.19 -1.48 23.21 -1.48 23.23 -1.48 23.26 -1.48 23.28 -1.48 23.31 -1.48 23.33 -1.48 23.36 -1.48 23.38 -1.48 23.41 -1.48 23.43 -1.48 23.46 -1.48 23.48 -1.48 23.51 -1.48 23.53 -1.48 23.56 -1.48 23.58 -1.48 23.61 -1.48 23.63 -1.48 23.66 -1.48 23.68 -1.48 23.71 -1.48 23.73 -1.48 23.76 -1.48 23.78 -1.48 23.81 -1.48 23.83 -1.48 23.86 -1.48 23.88 -1.48 23.91 -1.48 23.93 -1.48 23.96 -1.48 23.98 -1.48 24.01 -1.48 24.03 -1.48 24.06 -1.48 24.08 -1.48 24.11 -1.48 24.13 -1.48 24.16 -1.48 24.18 -1.48 24.21 -1.48 24.23 -1.48 24.25 -1.48 24.28 -1.48 24.3 -1.48 24.33 -1.48 24.35 -1.48 24.38 -1.48 24.4 -1.48 24.43 -1.48 24.45 -1.48 24.48 -1.48 24.5 -1.48 24.53 -1.48 24.55 -1.48 24.58 -1.48 24.6 -1.48 24.63 -1.48 24.65 -1.48 24.68 -1.48 24.7 -1.48 24.73 -1.48 24.75 -1.48 24.78 -1.48 24.8 -1.48 24.83 -1.48 24.85 -1.48 24.88 -1.48 24.9 -1.48 24.93 -1.48 24.95 -1.48 24.98 -1.48 25 -1.48 25.03 -1.48 25.05 -1.48 25.08 -1.48 25.1 -1.48 25.13 -1.48 25.15 -1.48 25.18 -1.48 25.2 -1.48 25.23 -1.48 25.25 -1.48 25.27 -1.48 25.3 -1.48 25.32 -1.48 25.35 -1.48 25.37 -1.48 25.4 -1.48 25.42 -1.48 25.45 -1.48 25.47 -1.48 25.5 -1.48 25.52 -1.48 25.55 -1.48 25.57 -1.48 25.6 -1.48 25.62 -1.48 25.65 -1.48 25.67 -1.48 25.7 -1.48 25.72 -1.48 25.75 -1.48 25.77 -1.48 25.8 -1.48 25.82 -1.48 25.85 -1.48 25.87 -1.48 25.9 -1.48 25.92 -1.48 25.95 -1.48 25.97 -1.48 26 -1.48 26.02 -1.48 26.05 -1.48 26.07 -1.48 26.1 -1.48 26.12 -1.48 26.15 -1.48 26.17 -1.48 26.2 -1.48 26.22 -1.48 26.25 -1.48 26.27 -1.48 26.29 -1.48 26.32 -1.48 26.34 -1.48 26.37 -1.48 26.39 -1.48 26.42 -1.48 26.44 -1.48 26.47 -1.48 26.49 -1.48 26.52 -1.48 26.54 -1.48 26.57 -1.48 26.59 -1.48 26.62 -1.48 26.64 -1.48 26.67 -1.48 26.69 -1.48 26.72 -1.48 26.74 -1.48 26.77 -1.48 26.79 -1.48 26.82 -1.48 26.84 -1.48 26.87 -1.48 26.89 -1.48 26.92 -1.48 26.94 -1.48 26.97 -1.48 26.99 -1.48 27.02 -1.48 27.04 -1.48 27.07 -1.48 27.09 -1.48 27.12 -1.48 27.14 -1.48 27.17 -1.48 27.19 -1.48 27.22 -1.48 27.24 -1.48 27.26 -1.48 27.29 -1.48 27.31 -1.48 27.34 -1.48 27.36 -1.48 27.39 -1.48 27.41 -1.48 27.44 -1.48 27.46 -1.48 27.49 -1.48 27.51 -1.48 27.54 -1.48 27.56 -1.48 27.59 -1.48 27.61 -1.48 27.64 -1.48 27.66 -1.48 27.69 -1.48 27.71 -1.48 27.74 -1.48 27.76 -1.48 27.79 -1.48 27.81 -1.48 27.84 -1.48 27.86 -1.48 27.89 -1.48 27.91 -1.48 27.94 -1.48 27.96 -1.48 27.99 -1.48 28.01 -1.48 28.04 -1.48 28.06 -1.48 28.09 -1.48 28.11 -1.48 28.14 -1.48 28.16 -1.48 28.19 -1.48 28.21 -1.48 28.24 -1.48 28.26 -1.48 28.28 -1.48 28.31 -1.48 28.33 -1.48 28.36 -1.48 28.38 -1.48 28.41 -1.48 28.43 -1.48 28.46 -1.48 28.48 -1.48 28.51 -1.48 28.53 -1.48 28.56 -1.48 28.58 -1.48 28.61 -1.48 28.63 -1.48 28.66 -1.48 28.68 -1.48 28.71 -1.48 28.73 -1.48 28.76 -1.48 28.78 -1.48 28.81 -1.48 28.83 -1.48 28.86 -1.48 28.88 -1.48 28.91 -1.48 28.93 -1.48 28.96 -1.48 28.98 -1.48 29.01 -1.48 29.03 -1.48 29.06 -1.48 29.08 -1.48 29.11 -1.48 29.13 -1.48 29.16 -1.48 29.18 -1.48 29.21 -1.48 29.23 -1.48 29.26 -1.48 29.28 -1.48 29.3 -1.48 29.33 -1.48 29.35 -1.48 29.38 -1.48 29.4 -1.48 29.43 -1.48 29.45 -1.48 29.48 -1.48 29.5 -1.48 29.53 -1.48 29.55 -1.48 29.58 -1.48 29.6 -1.48 29.63 -1.48 29.65 -1.48 29.68 -1.48 29.7 -1.48 29.73 -1.48 29.75 -1.48 29.78 -1.48 29.8 -1.48 29.83 -1.48 29.85 -1.48 29.88 -1.48 29.9 -1.48 29.93 -1.48 29.95 -1.48 29.98 -1.48 30 -1.48 30.03 -1.48 30.05 -1.48 30.08 -1.48 30.1 -1.48 30.13 -1.48 30.15 -1.48 30.18 -1.48 30.2 -1.48 30.23 -1.48 30.25 -1.48 30.28 -1.48 30.3 -1.48 30.32 -1.48 30.35 -1.48 30.37 -1.48 30.4 -1.48 30.42 -1.48 30.45 -1.48 30.47 -1.48 30.5 -1.48 30.52 -1.48 30.55 -1.48 30.57 -1.48 30.6 -1.48 30.62 -1.48 30.65 -1.48 30.67 -1.48 30.7 -1.48 30.72 -1.48 30.75 -1.48 30.77 -1.48 30.8 -1.48 30.82 -1.48 30.85 -1.48 30.87 -1.48 30.9 -1.48 30.92 -1.48 30.95 -1.48 30.97 -1.48 31 -1.48 31.02 -1.48 31.05 -1.48 31.07 -1.48 31.1 -1.48 31.12 -1.48 31.15 -1.48 31.17 -1.48 31.2 -1.48 31.22 -1.48 31.25 -1.48 31.27 -1.48 31.3 -1.48 31.32 -1.48 31.34 -1.48 31.37 -1.48 31.39 -1.48 31.42 -1.48 31.44 -1.48 31.47 -1.48 31.49 -1.48 31.52 -1.48 31.54 -1.48 31.57 -1.48 31.59 -1.48 31.62 -1.48 31.64 -1.48 31.67 -1.48 31.69 -1.48 31.72 -1.48 31.74 -1.48 31.77 -1.48 31.79 -1.48 31.82 -1.48 31.84 -1.48 31.87 -1.48 31.89 -1.48 31.92 -1.48 31.94 -1.48 31.97 -1.48 31.99 -1.48 32.02 -1.48 32.04 -1.48 32.07 -1.48 32.09 -1.48 32.12 -1.48 32.14 -1.48 32.17 -1.48 32.19 -1.48 32.22 -1.48 32.24 -1.48 32.27 -1.48 32.29 -1.48 32.31 -1.48 32.34 -1.48 32.36 -1.48 32.39 -1.48 32.41 -1.48 32.44 -1.48 32.46 -1.48 32.49 -1.48 32.51 -1.48 32.54 -1.48 32.56 -1.48 32.59 -1.48 32.61 -1.48 32.64 -1.48 32.66 -1.48 32.69 -1.48 32.71 -1.48 32.74 -1.48 32.76 -1.48 32.79 -1.48 32.81 -1.48 32.84 -1.48 32.86 -1.48 32.89 -1.48 32.91 -1.48 32.94 -1.48 32.96 -1.48 32.99 -1.48 33.01 -1.48 33.04 -1.48 33.06 -1.48 33.09 -1.48 33.11 -1.48 33.14 -1.48 33.16 -1.48 33.19 -1.48 33.21 -1.48 33.24 -1.48 33.26 -1.48 33.29 -1.48 33.31 -1.48 33.33 -1.48 33.36 -1.48 33.38 -1.48 33.41 -1.48 33.43 -1.48 33.46 -1.48 33.48 -1.48 33.51 -1.48 33.53 -1.48 33.56 -1.48 33.58 -1.48 33.61 -1.48 33.63 -1.48 33.66 -1.48 33.68 -1.48 33.71 -1.48 33.73 -1.48 33.76 -1.48 33.78 -1.48 33.81 -1.48 33.83 -1.48 33.86 -1.48 33.88 -1.48 33.91 -1.48 33.93 -1.48 33.96 -1.48 33.98 -1.48 34.01 -1.48 34.03 -1.48 34.06 -1.48 34.08 -1.48 34.11 -1.48 34.13 -1.48 34.16 -1.48 34.18 -1.48 34.21 -1.48 34.23 -1.48 34.26 -1.48 34.28 -1.48 34.31 -1.48 34.33 -1.48 34.35 -1.48 34.38 -1.48 34.4 -1.48 34.43 -1.48 34.45 -1.48 34.48 -1.48 34.5 -1.48 34.53 -1.48 34.55 -1.48 34.58 -1.48 34.6 -1.48 34.63 -1.48 34.65 -1.48 34.68 -1.48 34.7 -1.48 34.73 -1.48 34.75 -1.48 34.78 -1.48 34.8 -1.48 34.83 -1.48 34.85 -1.48 34.88 -1.48 34.9 -1.48 34.93 -1.48 34.95 -1.48 34.98 -1.48 35 -1.48 35.03 -1.48 35.05 -1.48 35.08 -1.48 35.1 -1.48 35.13 -1.48 35.15 -1.48 35.18 -1.48 35.2 -1.48 35.23 -1.48 35.25 -1.48 35.28 -1.48 35.3 -1.48 35.33 -1.48 35.35 -1.48 35.37 -1.48 35.4 -1.48 35.42 -1.48 35.45 -1.48 35.47 -1.48 35.5 -1.48 35.52 -1.48 35.55 -1.48 35.57 -1.48 35.6 -1.48 35.62 -1.48 35.65 -1.48 35.67 -1.48 35.7 -1.48 35.72 -1.48 35.75 -1.48 35.77 -1.48 35.8 -1.48 35.82 -1.48 35.85 -1.48 35.87 -1.48 35.9 -1.48 35.92 -1.48 35.95 -1.48 35.97 -1.48 36 -1.48 36.02 -1.48 36.05 -1.48 36.07 -1.48 36.1 -1.48 36.12 -1.48 36.15 -1.48 36.17 -1.48 36.2 -1.48 36.22 -1.48 36.25 -1.48 36.27 -1.48 36.3 -1.48 36.32 -1.48 36.35 -1.48 36.37 -1.48 36.39 -1.48 36.42 -1.48 36.44 -1.48 36.47 -1.48 36.49 -1.48 36.52 -1.48 36.54 -1.48 36.57 -1.48 36.59 -1.48 36.62 -1.48 36.64 -1.48 36.67 -1.48 36.69 -1.48 36.72 -1.48 36.74 -1.48 36.77 -1.48 36.79 -1.48 36.82 -1.48 36.84 -1.48 36.87 -1.48 36.89 -1.48 36.92 -1.48 36.94 -1.48 36.97 -1.48 36.99 -1.48 37.02 -1.48 37.04 -1.48 37.07 -1.48 37.09 -1.48 37.12 -1.48 37.14 -1.48 37.17 -1.48 37.19 -1.48 37.22 -1.48 37.24 -1.48 37.27 -1.48 37.29 -1.48 37.32 -1.48 37.34 -1.48 37.36 -1.48 37.39 -1.48 37.41 -1.48 37.44 -1.48 37.46 -1.48 37.49 -1.48 37.51 -1.48 37.54 -1.48 37.56 -1.48 37.59 -1.48 37.61 -1.48 37.64 -1.48 37.66 -1.48 37.69 -1.48 37.71 -1.48 37.74 -1.48 37.76 -1.48 37.79 -1.48 37.81 -1.48 37.84 -1.48 37.86 -1.48 37.89 -1.48 37.91 -1.48 37.94 -1.48 37.96 -1.48 37.99 -1.48 38.01 -1.48 38.04 -1.48 38.06 -1.48 38.09 -1.48 38.11 -1.48 38.14 -1.48 38.16 -1.48 38.19 -1.48 38.21 -1.48 38.24 -1.48 38.26 -1.48 38.29 -1.48 38.31 -1.48 38.34 -1.48 38.36 -1.48 38.38 -1.48 38.41 -1.48 38.43 -1.48 38.46 -1.48 38.48 -1.48 38.51 -1.48 38.53 -1.48 38.56 -1.48 38.58 -1.48 38.61 -1.48 38.63 -1.48 38.66 -1.48 38.68 -1.48 38.71 -1.48 38.73 -1.48 38.76 -1.48 38.78 -1.48 38.81 -1.48 38.83 -1.48 38.86 -1.48 38.88 -1.48 38.91 -1.48 38.93 -1.48 38.96 -1.48 38.98 -1.48 39.01 -1.48 39.03 -1.48 39.06 -1.48 39.08 -1.48 39.11 -1.48 39.13 -1.48 39.16 -1.48 39.18 -1.48 39.21 -1.48 39.23 -1.48 39.26 -1.48 39.28 -1.48 39.31 -1.48 39.33 -1.48 39.36 -1.48 39.38 -1.48 39.4 -1.48 39.43 -1.48 39.45 -1.48 39.48 -1.48 39.5 -1.48 39.53 -1.48 39.55 -1.48 39.58 -1.48 39.6 -1.48 39.63 -1.48 39.65 -1.48 39.68 -1.48 39.7 -1.48 39.73 -1.48 39.75 -1.48 39.78 -1.48 39.8 -1.48 39.83 -1.48 39.85 -1.48 39.88 -1.48 39.9 -1.48 39.93 -1.48 39.95 -1.48 39.98 -1.48 40 -1.48 40.03 -1.48 40.05 -1.48 40.08 -1.48 40.1 -1.48 40.13 -1.48 40.15 -1.48 40.18 -1.48 40.2 -1.48 40.23 -1.48 40.25 -1.48 40.28 -1.48 40.3 -1.48 40.33 -1.48 40.35 -1.48 40.38 -1.48 40.4 -1.48 40.42 -1.48 40.45 -1.48 40.47 -1.48 40.5 -1.48 40.52 -1.48 40.55 -1.48 40.57 -1.48 40.6 -1.48 40.62 -1.48 40.65 -1.48 40.67 -1.48 40.7 -1.48 40.72 -1.48 40.75 -1.48 40.77 -1.48 40.8 -1.48 40.82 -1.48 40.85 -1.48 40.87 -1.48 40.9 -1.48 40.92 -1.48 40.95 -1.48 40.97 -1.48 41 -1.48 41.02 -1.48 41.05 -1.48 41.07 -1.48 41.1 -1.48 41.12 -1.48 41.15 -1.48 41.17 -1.48 41.2 -1.48 41.22 -1.48 41.25 -1.48 41.27 -1.48 41.3 -1.48 41.32 -1.48 41.35 -1.48 41.37 -1.48 41.39 -1.48 41.42 -1.48 41.44 -1.48 41.47 -1.48 41.49 -1.48 41.52 -1.48 41.54 -1.48 41.57 -1.48 41.59 -1.48 41.62 -1.48 41.64 -1.48 41.67 -1.48 41.69 -1.48 41.72 -1.48 41.74 -1.48 41.77 -1.48 41.79 -1.48 41.82 -1.48 41.84 -1.48 41.87 -1.48 41.89 -1.48 41.92 -1.48 41.94 -1.48 41.97 -1.48 41.99 -1.48 42.02 -1.48 42.04 -1.48 42.07 -1.48 42.09 -1.48 42.12 -1.48 42.14 -1.48 42.17 -1.48 42.19 -1.48 42.22 -1.48 42.24 -1.48 42.27 -1.48 42.29 -1.48 42.32 -1.48 42.34 -1.48 42.37 -1.48 42.39 -1.48 42.41 -1.48 42.44 -1.48 42.46 -1.48 42.49 -1.48 42.51 -1.48 42.54 -1.48 42.56 -1.48 42.59 -1.48 42.61 -1.48 42.64 -1.48 42.66 -1.48 42.69 -1.48 42.71 -1.48 42.74 -1.48 42.76 -1.48 42.79 -1.48 42.81 -1.48 42.84 -1.48 42.86 -1.48 42.89 -1.48 42.91 -1.48 42.94 -1.48 42.96 -1.48 42.99 -1.48 43.01 -1.48 43.04 -1.48 43.06 -1.48 43.09 -1.48 43.11 -1.48 43.14 -1.48 43.16 -1.48 43.19 -1.48 43.21 -1.48 43.24 -1.48 43.26 -1.48 43.29 -1.48 43.31 -1.48 43.34 -1.48 43.36 -1.48 43.39 -1.48 43.41 -1.48 43.43 -1.48 43.46 -1.48 43.48 -1.48 43.51 -1.48 43.53 -1.48 43.56 -1.48 43.58 -1.48 43.61 -1.48 43.63 -1.48 43.66 -1.48 43.68 -1.48 43.71 -1.48 43.73 -1.48 43.76 -1.48 43.78 -1.48 43.81 -1.48 43.83 -1.48 43.86 -1.48 43.88 -1.48 43.91 -1.48 43.93 -1.48 43.96 -1.48 43.98 -1.48 44.01 -1.48 44.03 -1.48 44.06 -1.48 44.08 -1.48 44.11 -1.48 44.13 -1.48 44.16 -1.48 44.18 -1.48 44.21 -1.48 44.23 -1.48 44.26 -1.48 44.28 -1.48 44.31 -1.48 44.33 -1.48 44.36 -1.48 44.38 -1.48 44.41 -1.48 44.43 -1.48 44.45 -1.48 44.48 -1.48 44.5 -1.48 44.53 -1.48 44.55 -1.48 44.58 -1.48 44.6 -1.48 44.63 -1.48 44.65 -1.48 44.68 -1.48 44.7 -1.48 44.73 -1.48 44.75 -1.48 44.78 -1.48 44.8 -1.48 44.83 -1.48 44.85 -1.48 44.88 -1.48 44.9 -1.48 44.93 -1.48 44.95 -1.48 44.98 -1.48 45 -1.48 45.03 -1.48 45.05 -1.48 45.08 -1.48 45.1 -1.48 45.13 -1.48 45.15 -1.48 45.18 -1.48 45.2 -1.48 45.23 -1.48 45.25 -1.48 45.28 -1.48 45.3 -1.48 45.33 -1.48 45.35 -1.48 45.38 -1.48 45.4 -1.48" class="primitive"/>
          </g>
          <g transform="translate(67.03,27.62)" id="img-8311e7c5-746" class="geometry color_S_H" stroke="#00BFFF">
            <path fill="none" d="M-45.4,-13.38 L -45.38 -13.38 -45.35 -13.38 -45.33 -13.38 -45.3 -13.38 -45.28 -13.38 -45.25 -13.38 -45.23 -13.38 -45.2 -13.38 -45.18 -13.38 -45.15 -13.38 -45.13 -13.38 -45.1 -13.38 -45.08 -13.38 -45.05 -13.38 -45.03 -13.38 -45 -13.38 -44.98 -13.38 -44.95 -13.38 -44.93 -13.38 -44.9 -13.38 -44.88 -13.38 -44.85 -13.38 -44.83 -13.38 -44.8 -13.37 -44.78 -13.37 -44.75 -13.37 -44.73 -13.37 -44.7 -13.37 -44.68 -13.37 -44.65 -13.37 -44.63 -13.36 -44.6 -13.36 -44.58 -13.36 -44.55 -13.36 -44.53 -13.36 -44.5 -13.35 -44.48 -13.35 -44.45 -13.35 -44.43 -13.35 -44.41 -13.34 -44.38 -13.34 -44.36 -13.34 -44.33 -13.34 -44.31 -13.33 -44.28 -13.33 -44.26 -13.33 -44.23 -13.32 -44.21 -13.32 -44.18 -13.31 -44.16 -13.31 -44.13 -13.31 -44.11 -13.3 -44.08 -13.3 -44.06 -13.29 -44.03 -13.29 -44.01 -13.28 -43.98 -13.28 -43.96 -13.27 -43.93 -13.27 -43.91 -13.26 -43.88 -13.26 -43.86 -13.25 -43.83 -13.25 -43.81 -13.24 -43.78 -13.24 -43.76 -13.23 -43.73 -13.22 -43.71 -13.22 -43.68 -13.21 -43.66 -13.2 -43.63 -13.2 -43.61 -13.19 -43.58 -13.18 -43.56 -13.18 -43.53 -13.17 -43.51 -13.16 -43.48 -13.15 -43.46 -13.15 -43.43 -13.14 -43.41 -13.13 -43.39 -13.12 -43.36 -13.11 -43.34 -13.1 -43.31 -13.09 -43.29 -13.08 -43.26 -13.08 -43.24 -13.07 -43.21 -13.06 -43.19 -13.05 -43.16 -13.04 -43.14 -13.03 -43.11 -13.02 -43.09 -13.01 -43.06 -12.99 -43.04 -12.98 -43.01 -12.97 -42.99 -12.96 -42.96 -12.95 -42.94 -12.94 -42.91 -12.93 -42.89 -12.91 -42.86 -12.9 -42.84 -12.89 -42.81 -12.88 -42.79 -12.86 -42.76 -12.85 -42.74 -12.84 -42.71 -12.82 -42.69 -12.81 -42.66 -12.79 -42.64 -12.78 -42.61 -12.76 -42.59 -12.75 -42.56 -12.73 -42.54 -12.72 -42.51 -12.7 -42.49 -12.69 -42.46 -12.67 -42.44 -12.65 -42.41 -12.64 -42.39 -12.62 -42.37 -12.6 -42.34 -12.58 -42.32 -12.56 -42.29 -12.55 -42.27 -12.53 -42.24 -12.51 -42.22 -12.49 -42.19 -12.47 -42.17 -12.45 -42.14 -12.43 -42.12 -12.41 -42.09 -12.39 -42.07 -12.37 -42.04 -12.34 -42.02 -12.32 -41.99 -12.3 -41.97 -12.28 -41.94 -12.25 -41.92 -12.23 -41.89 -12.21 -41.87 -12.18 -41.84 -12.16 -41.82 -12.13 -41.79 -12.1 -41.77 -12.08 -41.74 -12.05 -41.72 -12.02 -41.69 -12 -41.67 -11.97 -41.64 -11.94 -41.62 -11.91 -41.59 -11.88 -41.57 -11.85 -41.54 -11.82 -41.52 -11.79 -41.49 -11.76 -41.47 -11.73 -41.44 -11.7 -41.42 -11.66 -41.39 -11.63 -41.37 -11.6 -41.35 -11.56 -41.32 -11.53 -41.3 -11.49 -41.27 -11.46 -41.25 -11.42 -41.22 -11.38 -41.2 -11.34 -41.17 -11.31 -41.15 -11.27 -41.12 -11.23 -41.1 -11.19 -41.07 -11.15 -41.05 -11.11 -41.02 -11.06 -41 -11.02 -40.97 -10.98 -40.95 -10.93 -40.92 -10.89 -40.9 -10.84 -40.87 -10.8 -40.85 -10.75 -40.82 -10.7 -40.8 -10.66 -40.77 -10.61 -40.75 -10.56 -40.72 -10.51 -40.7 -10.46 -40.67 -10.41 -40.65 -10.35 -40.62 -10.3 -40.6 -10.25 -40.57 -10.19 -40.55 -10.14 -40.52 -10.08 -40.5 -10.03 -40.47 -9.97 -40.45 -9.91 -40.42 -9.85 -40.4 -9.79 -40.38 -9.73 -40.35 -9.67 -40.33 -9.61 -40.3 -9.55 -40.28 -9.49 -40.25 -9.42 -40.23 -9.36 -40.2 -9.29 -40.18 -9.23 -40.15 -9.16 -40.13 -9.09 -40.1 -9.02 -40.08 -8.95 -40.05 -8.88 -40.03 -8.81 -40 -8.74 -39.98 -8.67 -39.95 -8.6 -39.93 -8.52 -39.9 -8.45 -39.88 -8.37 -39.85 -8.3 -39.83 -8.22 -39.8 -8.14 -39.78 -8.07 -39.75 -7.99 -39.73 -7.91 -39.7 -7.83 -39.68 -7.75 -39.65 -7.67 -39.63 -7.58 -39.6 -7.5 -39.58 -7.42 -39.55 -7.33 -39.53 -7.25 -39.5 -7.16 -39.48 -7.08 -39.45 -6.99 -39.43 -6.91 -39.4 -6.82 -39.38 -6.73 -39.36 -6.64 -39.33 -6.55 -39.31 -6.47 -39.28 -6.38 -39.26 -6.29 -39.23 -6.2 -39.21 -6.11 -39.18 -6.01 -39.16 -5.92 -39.13 -5.83 -39.11 -5.74 -39.08 -5.65 -39.06 -5.56 -39.03 -5.46 -39.01 -5.37 -38.98 -5.28 -38.96 -5.18 -38.93 -5.09 -38.91 -5 -38.88 -4.9 -38.86 -4.81 -38.83 -4.72 -38.81 -4.62 -38.78 -4.53 -38.76 -4.44 -38.73 -4.35 -38.71 -4.25 -38.68 -4.16 -38.66 -4.07 -38.63 -3.98 -38.61 -3.88 -38.58 -3.79 -38.56 -3.7 -38.53 -3.61 -38.51 -3.52 -38.48 -3.43 -38.46 -3.34 -38.43 -3.25 -38.41 -3.16 -38.38 -3.08 -38.36 -2.99 -38.34 -2.9 -38.31 -2.82 -38.29 -2.73 -38.26 -2.65 -38.24 -2.56 -38.21 -2.48 -38.19 -2.4 -38.16 -2.32 -38.14 -2.23 -38.11 -2.15 -38.09 -2.08 -38.06 -2 -38.04 -1.92 -38.01 -1.84 -37.99 -1.77 -37.96 -1.69 -37.94 -1.62 -37.91 -1.55 -37.89 -1.48 -37.86 -1.41 -37.84 -1.34 -37.81 -1.27 -37.79 -1.2 -37.76 -1.14 -37.74 -1.07 -37.71 -1.01 -37.69 -0.95 -37.66 -0.89 -37.64 -0.83 -37.61 -0.77 -37.59 -0.71 -37.56 -0.65 -37.54 -0.6 -37.51 -0.54 -37.49 -0.49 -37.46 -0.44 -37.44 -0.39 -37.41 -0.34 -37.39 -0.29 -37.36 -0.24 -37.34 -0.2 -37.32 -0.15 -37.29 -0.11 -37.27 -0.07 -37.24 -0.02 -37.22 0.02 -37.19 0.05 -37.17 0.09 -37.14 0.13 -37.12 0.17 -37.09 0.2 -37.07 0.23 -37.04 0.27 -37.02 0.3 -36.99 0.33 -36.97 0.36 -36.94 0.39 -36.92 0.42 -36.89 0.44 -36.87 0.47 -36.84 0.49 -36.82 0.52 -36.79 0.54 -36.77 0.56 -36.74 0.58 -36.72 0.61 -36.69 0.63 -36.67 0.64 -36.64 0.66 -36.62 0.68 -36.59 0.7 -36.57 0.71 -36.54 0.73 -36.52 0.74 -36.49 0.76 -36.47 0.77 -36.44 0.79 -36.42 0.8 -36.39 0.81 -36.37 0.82 -36.35 0.83 -36.32 0.84 -36.3 0.85 -36.27 0.86 -36.25 0.87 -36.22 0.88 -36.2 0.89 -36.17 0.9 -36.15 0.91 -36.12 0.91 -36.1 0.92 -36.07 0.93 -36.05 0.93 -36.02 0.94 -36 0.94 -35.97 0.95 -35.95 0.96 -35.92 0.96 -35.9 0.96 -35.87 0.97 -35.85 0.97 -35.82 0.98 -35.8 0.98 -35.77 0.98 -35.75 0.99 -35.72 0.99 -35.7 0.99 -35.67 1 -35.65 1 -35.62 1 -35.6 1 -35.57 1.01 -35.55 1.01 -35.52 1.01 -35.5 1.01 -35.47 1.01 -35.45 1.02 -35.42 1.02 -35.4 1.02 -35.37 1.02 -35.35 1.02 -35.33 1.02 -35.3 1.02 -35.28 1.03 -35.25 1.03 -35.23 1.03 -35.2 1.03 -35.18 1.03 -35.15 1.03 -35.13 1.03 -35.1 1.03 -35.08 1.03 -35.05 1.03 -35.03 1.03 -35 1.04 -34.98 1.04 -34.95 1.04 -34.93 1.04 -34.9 1.04 -34.88 1.04 -34.85 1.04 -34.83 1.04 -34.8 1.04 -34.78 1.04 -34.75 1.04 -34.73 1.04 -34.7 1.04 -34.68 1.04 -34.65 1.04 -34.63 1.04 -34.6 1.04 -34.58 1.04 -34.55 1.04 -34.53 1.04 -34.5 1.04 -34.48 1.04 -34.45 1.04 -34.43 1.05 -34.4 1.05 -34.38 1.05 -34.35 1.05 -34.33 1.05 -34.31 1.05 -34.28 1.05 -34.26 1.05 -34.23 1.05 -34.21 1.05 -34.18 1.05 -34.16 1.05 -34.13 1.05 -34.11 1.05 -34.08 1.05 -34.06 1.05 -34.03 1.05 -34.01 1.05 -33.98 1.05 -33.96 1.05 -33.93 1.05 -33.91 1.05 -33.88 1.05 -33.86 1.05 -33.83 1.05 -33.81 1.05 -33.78 1.05 -33.76 1.05 -33.73 1.05 -33.71 1.05 -33.68 1.05 -33.66 1.05 -33.63 1.05 -33.61 1.05 -33.58 1.05 -33.56 1.05 -33.53 1.05 -33.51 1.05 -33.48 1.05 -33.46 1.05 -33.43 1.05 -33.41 1.05 -33.38 1.05 -33.36 1.05 -33.33 1.05 -33.31 1.05 -33.29 1.05 -33.26 1.05 -33.24 1.05 -33.21 1.05 -33.19 1.05 -33.16 1.05 -33.14 1.05 -33.11 1.05 -33.09 1.05 -33.06 1.05 -33.04 1.05 -33.01 1.05 -32.99 1.05 -32.96 1.05 -32.94 1.05 -32.91 1.05 -32.89 1.05 -32.86 1.05 -32.84 1.05 -32.81 1.05 -32.79 1.05 -32.76 1.05 -32.74 1.05 -32.71 1.05 -32.69 1.05 -32.66 1.05 -32.64 1.05 -32.61 1.05 -32.59 1.05 -32.56 1.05 -32.54 1.05 -32.51 1.05 -32.49 1.05 -32.46 1.05 -32.44 1.05 -32.41 1.05 -32.39 1.05 -32.36 1.05 -32.34 1.05 -32.31 1.05 -32.29 1.05 -32.27 1.05 -32.24 1.05 -32.22 1.05 -32.19 1.05 -32.17 1.05 -32.14 1.05 -32.12 1.05 -32.09 1.05 -32.07 1.05 -32.04 1.05 -32.02 1.05 -31.99 1.05 -31.97 1.05 -31.94 1.05 -31.92 1.05 -31.89 1.05 -31.87 1.05 -31.84 1.05 -31.82 1.06 -31.79 1.06 -31.77 1.06 -31.74 1.06 -31.72 1.06 -31.69 1.06 -31.67 1.06 -31.64 1.06 -31.62 1.06 -31.59 1.06 -31.57 1.06 -31.54 1.06 -31.52 1.06 -31.49 1.06 -31.47 1.06 -31.44 1.06 -31.42 1.06 -31.39 1.06 -31.37 1.06 -31.34 1.06 -31.32 1.06 -31.3 1.06 -31.27 1.06 -31.25 1.06 -31.22 1.06 -31.2 1.06 -31.17 1.06 -31.15 1.06 -31.12 1.06 -31.1 1.06 -31.07 1.06 -31.05 1.06 -31.02 1.06 -31 1.06 -30.97 1.06 -30.95 1.06 -30.92 1.06 -30.9 1.06 -30.87 1.06 -30.85 1.06 -30.82 1.06 -30.8 1.06 -30.77 1.06 -30.75 1.06 -30.72 1.05 -30.7 1.05 -30.67 1.05 -30.65 1.05 -30.62 1.05 -30.6 1.05 -30.57 1.05 -30.55 1.05 -30.52 1.05 -30.5 1.05 -30.47 1.05 -30.45 1.05 -30.42 1.05 -30.4 1.05 -30.37 1.05 -30.35 1.05 -30.32 1.05 -30.3 1.05 -30.28 1.05 -30.25 1.05 -30.23 1.05 -30.2 1.05 -30.18 1.05 -30.15 1.05 -30.13 1.05 -30.1 1.05 -30.08 1.05 -30.05 1.05 -30.03 1.05 -30 1.05 -29.98 1.05 -29.95 1.05 -29.93 1.05 -29.9 1.05 -29.88 1.05 -29.85 1.05 -29.83 1.05 -29.8 1.05 -29.78 1.05 -29.75 1.05 -29.73 1.05 -29.7 1.05 -29.68 1.05 -29.65 1.05 -29.63 1.05 -29.6 1.05 -29.58 1.05 -29.55 1.05 -29.53 1.05 -29.5 1.05 -29.48 1.05 -29.45 1.05 -29.43 1.05 -29.4 1.05 -29.38 1.05 -29.35 1.05 -29.33 1.05 -29.3 1.05 -29.28 1.05 -29.26 1.05 -29.23 1.05 -29.21 1.05 -29.18 1.05 -29.16 1.05 -29.13 1.05 -29.11 1.05 -29.08 1.05 -29.06 1.05 -29.03 1.05 -29.01 1.05 -28.98 1.05 -28.96 1.05 -28.93 1.05 -28.91 1.05 -28.88 1.05 -28.86 1.05 -28.83 1.05 -28.81 1.05 -28.78 1.05 -28.76 1.05 -28.73 1.05 -28.71 1.05 -28.68 1.05 -28.66 1.05 -28.63 1.05 -28.61 1.05 -28.58 1.05 -28.56 1.05 -28.53 1.05 -28.51 1.05 -28.48 1.05 -28.46 1.05 -28.43 1.05 -28.41 1.05 -28.38 1.05 -28.36 1.05 -28.33 1.05 -28.31 1.05 -28.28 1.05 -28.26 1.05 -28.24 1.05 -28.21 1.05 -28.19 1.05 -28.16 1.05 -28.14 1.05 -28.11 1.05 -28.09 1.05 -28.06 1.05 -28.04 1.05 -28.01 1.05 -27.99 1.05 -27.96 1.05 -27.94 1.05 -27.91 1.05 -27.89 1.05 -27.86 1.05 -27.84 1.05 -27.81 1.05 -27.79 1.05 -27.76 1.05 -27.74 1.05 -27.71 1.05 -27.69 1.05 -27.66 1.05 -27.64 1.05 -27.61 1.05 -27.59 1.05 -27.56 1.05 -27.54 1.05 -27.51 1.05 -27.49 1.05 -27.46 1.05 -27.44 1.05 -27.41 1.05 -27.39 1.05 -27.36 1.05 -27.34 1.05 -27.31 1.05 -27.29 1.05 -27.26 1.05 -27.24 1.05 -27.22 1.05 -27.19 1.05 -27.17 1.05 -27.14 1.05 -27.12 1.05 -27.09 1.05 -27.07 1.05 -27.04 1.05 -27.02 1.05 -26.99 1.05 -26.97 1.05 -26.94 1.05 -26.92 1.05 -26.89 1.05 -26.87 1.05 -26.84 1.05 -26.82 1.05 -26.79 1.05 -26.77 1.05 -26.74 1.05 -26.72 1.05 -26.69 1.05 -26.67 1.05 -26.64 1.05 -26.62 1.05 -26.59 1.05 -26.57 1.05 -26.54 1.05 -26.52 1.05 -26.49 1.05 -26.47 1.05 -26.44 1.05 -26.42 1.05 -26.39 1.05 -26.37 1.05 -26.34 1.05 -26.32 1.05 -26.29 1.05 -26.27 1.05 -26.25 1.05 -26.22 1.05 -26.2 1.05 -26.17 1.05 -26.15 1.05 -26.12 1.05 -26.1 1.05 -26.07 1.05 -26.05 1.05 -26.02 1.05 -26 1.05 -25.97 1.04 -25.95 1.04 -25.92 1.04 -25.9 1.04 -25.87 1.04 -25.85 1.04 -25.82 1.04 -25.8 1.04 -25.77 1.04 -25.75 1.04 -25.72 1.04 -25.7 1.04 -25.67 1.04 -25.65 1.04 -25.62 1.04 -25.6 1.04 -25.57 1.04 -25.55 1.04 -25.52 1.04 -25.5 1.04 -25.47 1.04 -25.45 1.04 -25.42 1.04 -25.4 1.04 -25.37 1.04 -25.35 1.04 -25.32 1.04 -25.3 1.04 -25.27 1.04 -25.25 1.04 -25.23 1.04 -25.2 1.04 -25.18 1.04 -25.15 1.04 -25.13 1.04 -25.1 1.04 -25.08 1.04 -25.05 1.04 -25.03 1.04 -25 1.04 -24.98 1.04 -24.95 1.04 -24.93 1.04 -24.9 1.04 -24.88 1.04 -24.85 1.04 -24.83 1.04 -24.8 1.04 -24.78 1.04 -24.75 1.04 -24.73 1.04 -24.7 1.04 -24.68 1.04 -24.65 1.04 -24.63 1.04 -24.6 1.04 -24.58 1.04 -24.55 1.04 -24.53 1.04 -24.5 1.04 -24.48 1.04 -24.45 1.04 -24.43 1.04 -24.4 1.04 -24.38 1.04 -24.35 1.04 -24.33 1.04 -24.3 1.04 -24.28 1.04 -24.25 1.04 -24.23 1.04 -24.21 1.04 -24.18 1.04 -24.16 1.04 -24.13 1.04 -24.11 1.04 -24.08 1.04 -24.06 1.04 -24.03 1.04 -24.01 1.04 -23.98 1.04 -23.96 1.04 -23.93 1.04 -23.91 1.04 -23.88 1.04 -23.86 1.04 -23.83 1.04 -23.81 1.04 -23.78 1.04 -23.76 1.04 -23.73 1.04 -23.71 1.04 -23.68 1.04 -23.66 1.04 -23.63 1.04 -23.61 1.04 -23.58 1.04 -23.56 1.04 -23.53 1.04 -23.51 1.04 -23.48 1.04 -23.46 1.04 -23.43 1.04 -23.41 1.04 -23.38 1.03 -23.36 1.03 -23.33 1.03 -23.31 1.03 -23.28 1.03 -23.26 1.03 -23.23 1.03 -23.21 1.03 -23.19 1.03 -23.16 1.03 -23.14 1.03 -23.11 1.03 -23.09 1.03 -23.06 1.03 -23.04 1.03 -23.01 1.03 -22.99 1.03 -22.96 1.03 -22.94 1.03 -22.91 1.03 -22.89 1.03 -22.86 1.03 -22.84 1.03 -22.81 1.03 -22.79 1.03 -22.76 1.03 -22.74 1.03 -22.71 1.03 -22.69 1.03 -22.66 1.03 -22.64 1.03 -22.61 1.03 -22.59 1.03 -22.56 1.03 -22.54 1.03 -22.51 1.03 -22.49 1.03 -22.46 1.03 -22.44 1.03 -22.41 1.03 -22.39 1.03 -22.36 1.03 -22.34 1.03 -22.31 1.03 -22.29 1.03 -22.26 1.03 -22.24 1.03 -22.21 1.03 -22.19 1.03 -22.17 1.03 -22.14 1.03 -22.12 1.03 -22.09 1.03 -22.07 1.03 -22.04 1.03 -22.02 1.03 -21.99 1.03 -21.97 1.03 -21.94 1.03 -21.92 1.03 -21.89 1.03 -21.87 1.03 -21.84 1.03 -21.82 1.03 -21.79 1.03 -21.77 1.03 -21.74 1.03 -21.72 1.03 -21.69 1.03 -21.67 1.03 -21.64 1.03 -21.62 1.03 -21.59 1.03 -21.57 1.03 -21.54 1.03 -21.52 1.03 -21.49 1.03 -21.47 1.03 -21.44 1.03 -21.42 1.03 -21.39 1.03 -21.37 1.03 -21.34 1.03 -21.32 1.03 -21.29 1.03 -21.27 1.02 -21.24 1.02 -21.22 1.02 -21.2 1.02 -21.17 1.02 -21.15 1.02 -21.12 1.02 -21.1 1.02 -21.07 1.02 -21.05 1.02 -21.02 1.02 -21 1.02 -20.97 1.02 -20.95 1.02 -20.92 1.02 -20.9 1.02 -20.87 1.02 -20.85 1.02 -20.82 1.02 -20.8 1.02 -20.77 1.02 -20.75 1.02 -20.72 1.02 -20.7 1.02 -20.67 1.02 -20.65 1.02 -20.62 1.02 -20.6 1.02 -20.57 1.02 -20.55 1.02 -20.52 1.02 -20.5 1.02 -20.47 1.02 -20.45 1.02 -20.42 1.02 -20.4 1.02 -20.37 1.02 -20.35 1.02 -20.32 1.02 -20.3 1.02 -20.27 1.02 -20.25 1.02 -20.22 1.02 -20.2 1.02 -20.18 1.02 -20.15 1.02 -20.13 1.02 -20.1 1.02 -20.08 1.02 -20.05 1.02 -20.03 1.02 -20 1.02 -19.98 1.02 -19.95 1.02 -19.93 1.02 -19.9 1.02 -19.88 1.02 -19.85 1.02 -19.83 1.02 -19.8 1.02 -19.78 1.02 -19.75 1.02 -19.73 1.02 -19.7 1.02 -19.68 1.02 -19.65 1.02 -19.63 1.02 -19.6 1.02 -19.58 1.02 -19.55 1.02 -19.53 1.02 -19.5 1.02 -19.48 1.02 -19.45 1.02 -19.43 1.01 -19.4 1.01 -19.38 1.01 -19.35 1.01 -19.33 1.01 -19.3 1.01 -19.28 1.01 -19.25 1.01 -19.23 1.01 -19.2 1.01 -19.18 1.01 -19.16 1.01 -19.13 1.01 -19.11 1.01 -19.08 1.01 -19.06 1.01 -19.03 1.01 -19.01 1.01 -18.98 1.01 -18.96 1.01 -18.93 1.01 -18.91 1.01 -18.88 1.01 -18.86 1.01 -18.83 1.01 -18.81 1.01 -18.78 1.01 -18.76 1.01 -18.73 1.01 -18.71 1.01 -18.68 1.01 -18.66 1.01 -18.63 1.01 -18.61 1.01 -18.58 1.01 -18.56 1.01 -18.53 1.01 -18.51 1.01 -18.48 1.01 -18.46 1.01 -18.43 1.01 -18.41 1.01 -18.38 1.01 -18.36 1.01 -18.33 1.01 -18.31 1.01 -18.28 1.01 -18.26 1.01 -18.23 1.01 -18.21 1.01 -18.18 1.01 -18.16 1.01 -18.14 1.01 -18.11 1.01 -18.09 1.01 -18.06 1.01 -18.04 1.01 -18.01 1.01 -17.99 1.01 -17.96 1.01 -17.94 1.01 -17.91 1.01 -17.89 1.01 -17.86 1.01 -17.84 1.01 -17.81 1.01 -17.79 1.01 -17.76 1.01 -17.74 1.01 -17.71 1.01 -17.69 1 -17.66 1 -17.64 1 -17.61 1 -17.59 1 -17.56 1 -17.54 1 -17.51 1 -17.49 1 -17.46 1 -17.44 1 -17.41 1 -17.39 1 -17.36 1 -17.34 1 -17.31 1 -17.29 1 -17.26 1 -17.24 1 -17.21 1 -17.19 1 -17.16 1 -17.14 1 -17.12 1 -17.09 1 -17.07 1 -17.04 1 -17.02 1 -16.99 1 -16.97 1 -16.94 1 -16.92 1 -16.89 1 -16.87 1 -16.84 1 -16.82 1 -16.79 1 -16.77 1 -16.74 1 -16.72 1 -16.69 1 -16.67 1 -16.64 1 -16.62 1 -16.59 1 -16.57 1 -16.54 1 -16.52 1 -16.49 1 -16.47 1 -16.44 1 -16.42 1 -16.39 1 -16.37 1 -16.34 1 -16.32 1 -16.29 1 -16.27 1 -16.24 1 -16.22 1 -16.19 1 -16.17 1 -16.15 1 -16.12 1 -16.1 1 -16.07 1 -16.05 0.99 -16.02 0.99 -16 0.99 -15.97 0.99 -15.95 0.99 -15.92 0.99 -15.9 0.99 -15.87 0.99 -15.85 0.99 -15.82 0.99 -15.8 0.99 -15.77 0.99 -15.75 0.99 -15.72 0.99 -15.7 0.99 -15.67 0.99 -15.65 0.99 -15.62 0.99 -15.6 0.99 -15.57 0.99 -15.55 0.99 -15.52 0.99 -15.5 0.99 -15.47 0.99 -15.45 0.99 -15.42 0.99 -15.4 0.99 -15.37 0.99 -15.35 0.99 -15.32 0.99 -15.3 0.99 -15.27 0.99 -15.25 0.99 -15.22 0.99 -15.2 0.99 -15.17 0.99 -15.15 0.99 -15.13 0.99 -15.1 0.99 -15.08 0.99 -15.05 0.99 -15.03 0.99 -15 0.99 -14.98 0.99 -14.95 0.99 -14.93 0.99 -14.9 0.99 -14.88 0.99 -14.85 0.99 -14.83 0.99 -14.8 0.99 -14.78 0.99 -14.75 0.99 -14.73 0.99 -14.7 0.99 -14.68 0.99 -14.65 0.99 -14.63 0.99 -14.6 0.99 -14.58 0.99 -14.55 0.99 -14.53 0.99 -14.5 0.99 -14.48 0.99 -14.45 0.98 -14.43 0.98 -14.4 0.98 -14.38 0.98 -14.35 0.98 -14.33 0.98 -14.3 0.98 -14.28 0.98 -14.25 0.98 -14.23 0.98 -14.2 0.98 -14.18 0.98 -14.15 0.98 -14.13 0.98 -14.11 0.98 -14.08 0.98 -14.06 0.98 -14.03 0.98 -14.01 0.98 -13.98 0.98 -13.96 0.98 -13.93 0.98 -13.91 0.98 -13.88 0.98 -13.86 0.98 -13.83 0.98 -13.81 0.98 -13.78 0.98 -13.76 0.98 -13.73 0.98 -13.71 0.98 -13.68 0.98 -13.66 0.98 -13.63 0.98 -13.61 0.98 -13.58 0.98 -13.56 0.98 -13.53 0.98 -13.51 0.98 -13.48 0.98 -13.46 0.98 -13.43 0.98 -13.41 0.98 -13.38 0.98 -13.36 0.98 -13.33 0.98 -13.31 0.98 -13.28 0.98 -13.26 0.98 -13.23 0.98 -13.21 0.98 -13.18 0.98 -13.16 0.98 -13.13 0.98 -13.11 0.98 -13.09 0.98 -13.06 0.98 -13.04 0.98 -13.01 0.98 -12.99 0.98 -12.96 0.98 -12.94 0.98 -12.91 0.98 -12.89 0.98 -12.86 0.98 -12.84 0.97 -12.81 0.97 -12.79 0.97 -12.76 0.97 -12.74 0.97 -12.71 0.97 -12.69 0.97 -12.66 0.97 -12.64 0.97 -12.61 0.97 -12.59 0.97 -12.56 0.97 -12.54 0.97 -12.51 0.97 -12.49 0.97 -12.46 0.97 -12.44 0.97 -12.41 0.97 -12.39 0.97 -12.36 0.97 -12.34 0.97 -12.31 0.97 -12.29 0.97 -12.26 0.97 -12.24 0.97 -12.21 0.97 -12.19 0.97 -12.16 0.97 -12.14 0.97 -12.12 0.97 -12.09 0.97 -12.07 0.97 -12.04 0.97 -12.02 0.97 -11.99 0.97 -11.97 0.97 -11.94 0.97 -11.92 0.97 -11.89 0.97 -11.87 0.97 -11.84 0.97 -11.82 0.97 -11.79 0.97 -11.77 0.97 -11.74 0.97 -11.72 0.97 -11.69 0.97 -11.67 0.97 -11.64 0.97 -11.62 0.97 -11.59 0.97 -11.57 0.97 -11.54 0.97 -11.52 0.97 -11.49 0.97 -11.47 0.97 -11.44 0.97 -11.42 0.97 -11.39 0.97 -11.37 0.97 -11.34 0.97 -11.32 0.97 -11.29 0.97 -11.27 0.97 -11.24 0.97 -11.22 0.97 -11.19 0.97 -11.17 0.96 -11.14 0.96 -11.12 0.96 -11.1 0.96 -11.07 0.96 -11.05 0.96 -11.02 0.96 -11 0.96 -10.97 0.96 -10.95 0.96 -10.92 0.96 -10.9 0.96 -10.87 0.96 -10.85 0.96 -10.82 0.96 -10.8 0.96 -10.77 0.96 -10.75 0.96 -10.72 0.96 -10.7 0.96 -10.67 0.96 -10.65 0.96 -10.62 0.96 -10.6 0.96 -10.57 0.96 -10.55 0.96 -10.52 0.96 -10.5 0.96 -10.47 0.96 -10.45 0.96 -10.42 0.96 -10.4 0.96 -10.37 0.96 -10.35 0.96 -10.32 0.96 -10.3 0.96 -10.27 0.96 -10.25 0.96 -10.22 0.96 -10.2 0.96 -10.17 0.96 -10.15 0.96 -10.12 0.96 -10.1 0.96 -10.08 0.96 -10.05 0.96 -10.03 0.96 -10 0.96 -9.98 0.96 -9.95 0.96 -9.93 0.96 -9.9 0.96 -9.88 0.96 -9.85 0.96 -9.83 0.96 -9.8 0.96 -9.78 0.96 -9.75 0.96 -9.73 0.96 -9.7 0.96 -9.68 0.96 -9.65 0.96 -9.63 0.96 -9.6 0.96 -9.58 0.96 -9.55 0.96 -9.53 0.96 -9.5 0.96 -9.48 0.96 -9.45 0.96 -9.43 0.96 -9.4 0.95 -9.38 0.95 -9.35 0.95 -9.33 0.95 -9.3 0.95 -9.28 0.95 -9.25 0.95 -9.23 0.95 -9.2 0.95 -9.18 0.95 -9.15 0.95 -9.13 0.95 -9.1 0.95 -9.08 0.95 -9.06 0.95 -9.03 0.95 -9.01 0.95 -8.98 0.95 -8.96 0.95 -8.93 0.95 -8.91 0.95 -8.88 0.95 -8.86 0.95 -8.83 0.95 -8.81 0.95 -8.78 0.95 -8.76 0.95 -8.73 0.95 -8.71 0.95 -8.68 0.95 -8.66 0.95 -8.63 0.95 -8.61 0.95 -8.58 0.95 -8.56 0.95 -8.53 0.95 -8.51 0.95 -8.48 0.95 -8.46 0.95 -8.43 0.95 -8.41 0.95 -8.38 0.95 -8.36 0.95 -8.33 0.95 -8.31 0.95 -8.28 0.95 -8.26 0.95 -8.23 0.95 -8.21 0.95 -8.18 0.95 -8.16 0.95 -8.13 0.95 -8.11 0.95 -8.08 0.95 -8.06 0.95 -8.04 0.95 -8.01 0.95 -7.99 0.95 -7.96 0.95 -7.94 0.95 -7.91 0.95 -7.89 0.95 -7.86 0.95 -7.84 0.95 -7.81 0.95 -7.79 0.95 -7.76 0.95 -7.74 0.95 -7.71 0.95 -7.69 0.95 -7.66 0.95 -7.64 0.95 -7.61 0.95 -7.59 0.95 -7.56 0.95 -7.54 0.95 -7.51 0.95 -7.49 0.94 -7.46 0.94 -7.44 0.94 -7.41 0.94 -7.39 0.94 -7.36 0.94 -7.34 0.94 -7.31 0.94 -7.29 0.94 -7.26 0.94 -7.24 0.94 -7.21 0.94 -7.19 0.94 -7.16 0.94 -7.14 0.94 -7.11 0.94 -7.09 0.94 -7.07 0.94 -7.04 0.94 -7.02 0.94 -6.99 0.94 -6.97 0.94 -6.94 0.94 -6.92 0.94 -6.89 0.94 -6.87 0.94 -6.84 0.94 -6.82 0.94 -6.79 0.94 -6.77 0.94 -6.74 0.94 -6.72 0.94 -6.69 0.94 -6.67 0.94 -6.64 0.94 -6.62 0.94 -6.59 0.94 -6.57 0.94 -6.54 0.94 -6.52 0.94 -6.49 0.94 -6.47 0.94 -6.44 0.94 -6.42 0.94 -6.39 0.94 -6.37 0.94 -6.34 0.94 -6.32 0.94 -6.29 0.94 -6.27 0.94 -6.24 0.94 -6.22 0.94 -6.19 0.94 -6.17 0.94 -6.14 0.94 -6.12 0.94 -6.09 0.94 -6.07 0.94 -6.05 0.94 -6.02 0.94 -6 0.94 -5.97 0.94 -5.95 0.94 -5.92 0.94 -5.9 0.94 -5.87 0.94 -5.85 0.94 -5.82 0.94 -5.8 0.94 -5.77 0.94 -5.75 0.94 -5.72 0.94 -5.7 0.94 -5.67 0.94 -5.65 0.94 -5.62 0.94 -5.6 0.94 -5.57 0.94 -5.55 0.94 -5.52 0.94 -5.5 0.94 -5.47 0.94 -5.45 0.94 -5.42 0.94 -5.4 0.94 -5.37 0.94 -5.35 0.94 -5.32 0.94 -5.3 0.93 -5.27 0.93 -5.25 0.93 -5.22 0.93 -5.2 0.93 -5.17 0.93 -5.15 0.93 -5.12 0.93 -5.1 0.93 -5.07 0.93 -5.05 0.93 -5.03 0.93 -5 0.93 -4.98 0.93 -4.95 0.93 -4.93 0.93 -4.9 0.93 -4.88 0.93 -4.85 0.93 -4.83 0.93 -4.8 0.93 -4.78 0.93 -4.75 0.93 -4.73 0.93 -4.7 0.93 -4.68 0.93 -4.65 0.93 -4.63 0.93 -4.6 0.93 -4.58 0.93 -4.55 0.93 -4.53 0.93 -4.5 0.93 -4.48 0.93 -4.45 0.93 -4.43 0.93 -4.4 0.93 -4.38 0.93 -4.35 0.93 -4.33 0.93 -4.3 0.93 -4.28 0.93 -4.25 0.93 -4.23 0.93 -4.2 0.93 -4.18 0.93 -4.15 0.93 -4.13 0.93 -4.1 0.93 -4.08 0.93 -4.05 0.93 -4.03 0.93 -4.01 0.93 -3.98 0.93 -3.96 0.93 -3.93 0.93 -3.91 0.93 -3.88 0.93 -3.86 0.93 -3.83 0.93 -3.81 0.93 -3.78 0.93 -3.76 0.93 -3.73 0.93 -3.71 0.93 -3.68 0.93 -3.66 0.93 -3.63 0.93 -3.61 0.93 -3.58 0.93 -3.56 0.93 -3.53 0.93 -3.51 0.93 -3.48 0.93 -3.46 0.93 -3.43 0.93 -3.41 0.93 -3.38 0.93 -3.36 0.93 -3.33 0.93 -3.31 0.93 -3.28 0.93 -3.26 0.93 -3.23 0.93 -3.21 0.93 -3.18 0.93 -3.16 0.93 -3.13 0.93 -3.11 0.93 -3.08 0.93 -3.06 0.93 -3.03 0.93 -3.01 0.93 -2.99 0.93 -2.96 0.93 -2.94 0.93 -2.91 0.93 -2.89 0.93 -2.86 0.93 -2.84 0.93 -2.81 0.93 -2.79 0.93 -2.76 0.93 -2.74 0.93 -2.71 0.93 -2.69 0.93 -2.66 0.92 -2.64 0.92 -2.61 0.92 -2.59 0.92 -2.56 0.92 -2.54 0.92 -2.51 0.92 -2.49 0.92 -2.46 0.92 -2.44 0.92 -2.41 0.92 -2.39 0.92 -2.36 0.92 -2.34 0.92 -2.31 0.92 -2.29 0.92 -2.26 0.92 -2.24 0.92 -2.21 0.92 -2.19 0.92 -2.16 0.92 -2.14 0.92 -2.11 0.92 -2.09 0.92 -2.06 0.92 -2.04 0.92 -2.02 0.92 -1.99 0.92 -1.97 0.92 -1.94 0.92 -1.92 0.92 -1.89 0.92 -1.87 0.92 -1.84 0.92 -1.82 0.92 -1.79 0.92 -1.77 0.92 -1.74 0.92 -1.72 0.92 -1.69 0.92 -1.67 0.92 -1.64 0.92 -1.62 0.92 -1.59 0.92 -1.57 0.92 -1.54 0.92 -1.52 0.92 -1.49 0.92 -1.47 0.92 -1.44 0.92 -1.42 0.92 -1.39 0.92 -1.37 0.92 -1.34 0.92 -1.32 0.92 -1.29 0.92 -1.27 0.92 -1.24 0.92 -1.22 0.92 -1.19 0.92 -1.17 0.92 -1.14 0.92 -1.12 0.92 -1.09 0.92 -1.07 0.92 -1.04 0.92 -1.02 0.92 -1 0.92 -0.97 0.92 -0.95 0.92 -0.92 0.92 -0.9 0.92 -0.87 0.92 -0.85 0.92 -0.82 0.92 -0.8 0.92 -0.77 0.92 -0.75 0.92 -0.72 0.92 -0.7 0.92 -0.67 0.92 -0.65 0.92 -0.62 0.92 -0.6 0.92 -0.57 0.92 -0.55 0.92 -0.52 0.92 -0.5 0.92 -0.47 0.92 -0.45 0.92 -0.42 0.92 -0.4 0.92 -0.37 0.92 -0.35 0.92 -0.32 0.92 -0.3 0.92 -0.27 0.92 -0.25 0.92 -0.22 0.92 -0.2 0.92 -0.17 0.92 -0.15 0.92 -0.12 0.92 -0.1 0.92 -0.07 0.92 -0.05 0.92 -0.02 0.92 0 0.92 0.02 0.92 0.05 0.92 0.07 0.92 0.1 0.92 0.12 0.92 0.15 0.92 0.17 0.92 0.2 0.92 0.22 0.92 0.25 0.92 0.27 0.92 0.3 0.92 0.32 0.92 0.35 0.92 0.37 0.92 0.4 0.92 0.42 0.92 0.45 0.92 0.47 0.92 0.5 0.92 0.52 0.92 0.55 0.92 0.57 0.92 0.6 0.92 0.62 0.92 0.65 0.92 0.67 0.92 0.7 0.92 0.72 0.92 0.75 0.92 0.77 0.92 0.8 0.92 0.82 0.92 0.85 0.92 0.87 0.92 0.9 0.92 0.92 0.91 0.95 0.91 0.97 0.91 1 0.91 1.02 0.91 1.04 0.91 1.07 0.91 1.09 0.91 1.12 0.91 1.14 0.91 1.17 0.91 1.19 0.91 1.22 0.91 1.24 0.91 1.27 0.91 1.29 0.91 1.32 0.91 1.34 0.91 1.37 0.91 1.39 0.91 1.42 0.91 1.44 0.91 1.47 0.91 1.49 0.91 1.52 0.91 1.54 0.91 1.57 0.91 1.59 0.91 1.62 0.91 1.64 0.91 1.67 0.91 1.69 0.91 1.72 0.91 1.74 0.91 1.77 0.91 1.79 0.91 1.82 0.91 1.84 0.91 1.87 0.91 1.89 0.91 1.92 0.91 1.94 0.91 1.97 0.91 1.99 0.91 2.02 0.91 2.04 0.91 2.06 0.91 2.09 0.91 2.11 0.91 2.14 0.91 2.16 0.91 2.19 0.91 2.21 0.91 2.24 0.91 2.26 0.91 2.29 0.91 2.31 0.91 2.34 0.91 2.36 0.91 2.39 0.91 2.41 0.91 2.44 0.91 2.46 0.91 2.49 0.91 2.51 0.91 2.54 0.91 2.56 0.91 2.59 0.91 2.61 0.91 2.64 0.91 2.66 0.91 2.69 0.91 2.71 0.91 2.74 0.91 2.76 0.91 2.79 0.91 2.81 0.91 2.84 0.91 2.86 0.91 2.89 0.91 2.91 0.91 2.94 0.91 2.96 0.91 2.99 0.91 3.01 0.91 3.03 0.91 3.06 0.91 3.08 0.91 3.11 0.91 3.13 0.91 3.16 0.91 3.18 0.91 3.21 0.91 3.23 0.91 3.26 0.91 3.28 0.91 3.31 0.91 3.33 0.91 3.36 0.91 3.38 0.91 3.41 0.91 3.43 0.91 3.46 0.91 3.48 0.91 3.51 0.91 3.53 0.91 3.56 0.91 3.58 0.91 3.61 0.91 3.63 0.91 3.66 0.91 3.68 0.91 3.71 0.91 3.73 0.91 3.76 0.91 3.78 0.91 3.81 0.91 3.83 0.91 3.86 0.91 3.88 0.91 3.91 0.91 3.93 0.91 3.96 0.91 3.98 0.91 4.01 0.91 4.03 0.91 4.05 0.91 4.08 0.91 4.1 0.91 4.13 0.91 4.15 0.91 4.18 0.91 4.2 0.91 4.23 0.91 4.25 0.91 4.28 0.91 4.3 0.91 4.33 0.91 4.35 0.91 4.38 0.91 4.4 0.91 4.43 0.91 4.45 0.91 4.48 0.91 4.5 0.91 4.53 0.91 4.55 0.91 4.58 0.91 4.6 0.91 4.63 0.91 4.65 0.91 4.68 0.91 4.7 0.91 4.73 0.91 4.75 0.91 4.78 0.91 4.8 0.91 4.83 0.91 4.85 0.91 4.88 0.91 4.9 0.91 4.93 0.91 4.95 0.91 4.98 0.91 5 0.91 5.03 0.91 5.05 0.91 5.07 0.91 5.1 0.91 5.12 0.91 5.15 0.91 5.17 0.91 5.2 0.91 5.22 0.91 5.25 0.91 5.27 0.91 5.3 0.91 5.32 0.91 5.35 0.91 5.37 0.91 5.4 0.91 5.42 0.91 5.45 0.91 5.47 0.91 5.5 0.91 5.52 0.91 5.55 0.91 5.57 0.91 5.6 0.91 5.62 0.91 5.65 0.91 5.67 0.91 5.7 0.91 5.72 0.91 5.75 0.91 5.77 0.91 5.8 0.91 5.82 0.91 5.85 0.91 5.87 0.91 5.9 0.91 5.92 0.91 5.95 0.91 5.97 0.91 6 0.91 6.02 0.91 6.05 0.91 6.07 0.91 6.09 0.91 6.12 0.91 6.14 0.91 6.17 0.91 6.19 0.91 6.22 0.91 6.24 0.91 6.27 0.91 6.29 0.91 6.32 0.91 6.34 0.91 6.37 0.91 6.39 0.91 6.42 0.91 6.44 0.91 6.47 0.91 6.49 0.91 6.52 0.91 6.54 0.91 6.57 0.91 6.59 0.91 6.62 0.91 6.64 0.91 6.67 0.91 6.69 0.91 6.72 0.91 6.74 0.91 6.77 0.91 6.79 0.91 6.82 0.91 6.84 0.91 6.87 0.91 6.89 0.91 6.92 0.91 6.94 0.91 6.97 0.91 6.99 0.91 7.02 0.91 7.04 0.91 7.07 0.91 7.09 0.91 7.11 0.91 7.14 0.91 7.16 0.91 7.19 0.91 7.21 0.91 7.24 0.91 7.26 0.91 7.29 0.91 7.31 0.91 7.34 0.91 7.36 0.91 7.39 0.91 7.41 0.91 7.44 0.91 7.46 0.91 7.49 0.91 7.51 0.91 7.54 0.91 7.56 0.91 7.59 0.91 7.61 0.91 7.64 0.91 7.66 0.91 7.69 0.91 7.71 0.91 7.74 0.91 7.76 0.91 7.79 0.91 7.81 0.91 7.84 0.9 7.86 0.9 7.89 0.9 7.91 0.9 7.94 0.9 7.96 0.9 7.99 0.9 8.01 0.9 8.04 0.9 8.06 0.9 8.08 0.9 8.11 0.9 8.13 0.9 8.16 0.9 8.18 0.9 8.21 0.9 8.23 0.9 8.26 0.9 8.28 0.9 8.31 0.9 8.33 0.9 8.36 0.9 8.38 0.9 8.41 0.9 8.43 0.9 8.46 0.9 8.48 0.9 8.51 0.9 8.53 0.9 8.56 0.9 8.58 0.9 8.61 0.9 8.63 0.9 8.66 0.9 8.68 0.9 8.71 0.9 8.73 0.9 8.76 0.9 8.78 0.9 8.81 0.9 8.83 0.9 8.86 0.9 8.88 0.9 8.91 0.9 8.93 0.9 8.96 0.9 8.98 0.9 9.01 0.9 9.03 0.9 9.06 0.9 9.08 0.9 9.1 0.9 9.13 0.9 9.15 0.9 9.18 0.9 9.2 0.9 9.23 0.9 9.25 0.9 9.28 0.9 9.3 0.9 9.33 0.9 9.35 0.9 9.38 0.9 9.4 0.9 9.43 0.9 9.45 0.9 9.48 0.9 9.5 0.9 9.53 0.9 9.55 0.9 9.58 0.9 9.6 0.9 9.63 0.9 9.65 0.9 9.68 0.9 9.7 0.9 9.73 0.9 9.75 0.9 9.78 0.9 9.8 0.9 9.83 0.9 9.85 0.9 9.88 0.9 9.9 0.9 9.93 0.9 9.95 0.9 9.98 0.9 10 0.9 10.03 0.9 10.05 0.9 10.08 0.9 10.1 0.9 10.12 0.9 10.15 0.9 10.17 0.9 10.2 0.9 10.22 0.9 10.25 0.9 10.27 0.9 10.3 0.9 10.32 0.9 10.35 0.9 10.37 0.9 10.4 0.9 10.42 0.9 10.45 0.9 10.47 0.9 10.5 0.9 10.52 0.9 10.55 0.9 10.57 0.9 10.6 0.9 10.62 0.9 10.65 0.9 10.67 0.9 10.7 0.9 10.72 0.9 10.75 0.9 10.77 0.9 10.8 0.9 10.82 0.9 10.85 0.9 10.87 0.9 10.9 0.9 10.92 0.9 10.95 0.9 10.97 0.9 11 0.9 11.02 0.9 11.05 0.9 11.07 0.9 11.1 0.9 11.12 0.9 11.14 0.9 11.17 0.9 11.19 0.9 11.22 0.9 11.24 0.9 11.27 0.9 11.29 0.9 11.32 0.9 11.34 0.9 11.37 0.9 11.39 0.9 11.42 0.9 11.44 0.9 11.47 0.9 11.49 0.9 11.52 0.9 11.54 0.9 11.57 0.9 11.59 0.9 11.62 0.9 11.64 0.9 11.67 0.9 11.69 0.9 11.72 0.9 11.74 0.9 11.77 0.9 11.79 0.9 11.82 0.9 11.84 0.9 11.87 0.9 11.89 0.9 11.92 0.9 11.94 0.9 11.97 0.9 11.99 0.9 12.02 0.9 12.04 0.9 12.07 0.9 12.09 0.9 12.12 0.9 12.14 0.9 12.16 0.9 12.19 0.9 12.21 0.9 12.24 0.9 12.26 0.9 12.29 0.9 12.31 0.9 12.34 0.9 12.36 0.9 12.39 0.9 12.41 0.9 12.44 0.9 12.46 0.9 12.49 0.9 12.51 0.9 12.54 0.9 12.56 0.9 12.59 0.9 12.61 0.9 12.64 0.9 12.66 0.9 12.69 0.9 12.71 0.9 12.74 0.9 12.76 0.9 12.79 0.9 12.81 0.9 12.84 0.9 12.86 0.9 12.89 0.9 12.91 0.9 12.94 0.9 12.96 0.9 12.99 0.9 13.01 0.9 13.04 0.9 13.06 0.9 13.09 0.9 13.11 0.9 13.13 0.9 13.16 0.9 13.18 0.9 13.21 0.9 13.23 0.9 13.26 0.9 13.28 0.9 13.31 0.9 13.33 0.9 13.36 0.9 13.38 0.9 13.41 0.9 13.43 0.9 13.46 0.9 13.48 0.9 13.51 0.9 13.53 0.9 13.56 0.9 13.58 0.9 13.61 0.9 13.63 0.9 13.66 0.9 13.68 0.9 13.71 0.9 13.73 0.9 13.76 0.9 13.78 0.9 13.81 0.9 13.83 0.9 13.86 0.9 13.88 0.9 13.91 0.9 13.93 0.9 13.96 0.9 13.98 0.9 14.01 0.9 14.03 0.9 14.06 0.9 14.08 0.9 14.11 0.9 14.13 0.9 14.15 0.9 14.18 0.9 14.2 0.9 14.23 0.9 14.25 0.9 14.28 0.9 14.3 0.9 14.33 0.9 14.35 0.9 14.38 0.9 14.4 0.9 14.43 0.9 14.45 0.9 14.48 0.9 14.5 0.9 14.53 0.9 14.55 0.9 14.58 0.9 14.6 0.9 14.63 0.9 14.65 0.9 14.68 0.9 14.7 0.9 14.73 0.9 14.75 0.9 14.78 0.9 14.8 0.9 14.83 0.9 14.85 0.9 14.88 0.9 14.9 0.9 14.93 0.9 14.95 0.9 14.98 0.9 15 0.9 15.03 0.9 15.05 0.9 15.08 0.9 15.1 0.9 15.13 0.9 15.15 0.9 15.17 0.9 15.2 0.9 15.22 0.9 15.25 0.9 15.27 0.9 15.3 0.9 15.32 0.9 15.35 0.9 15.37 0.9 15.4 0.9 15.42 0.9 15.45 0.9 15.47 0.9 15.5 0.9 15.52 0.9 15.55 0.9 15.57 0.9 15.6 0.9 15.62 0.9 15.65 0.9 15.67 0.9 15.7 0.9 15.72 0.9 15.75 0.9 15.77 0.9 15.8 0.9 15.82 0.9 15.85 0.9 15.87 0.9 15.9 0.9 15.92 0.9 15.95 0.9 15.97 0.9 16 0.9 16.02 0.9 16.05 0.9 16.07 0.9 16.1 0.9 16.12 0.9 16.15 0.9 16.17 0.9 16.19 0.9 16.22 0.9 16.24 0.9 16.27 0.9 16.29 0.9 16.32 0.9 16.34 0.9 16.37 0.9 16.39 0.9 16.42 0.9 16.44 0.9 16.47 0.9 16.49 0.9 16.52 0.9 16.54 0.9 16.57 0.9 16.59 0.9 16.62 0.9 16.64 0.9 16.67 0.9 16.69 0.9 16.72 0.9 16.74 0.9 16.77 0.9 16.79 0.9 16.82 0.9 16.84 0.9 16.87 0.9 16.89 0.9 16.92 0.9 16.94 0.9 16.97 0.9 16.99 0.9 17.02 0.9 17.04 0.9 17.07 0.9 17.09 0.9 17.12 0.9 17.14 0.9 17.16 0.9 17.19 0.9 17.21 0.9 17.24 0.9 17.26 0.9 17.29 0.9 17.31 0.9 17.34 0.9 17.36 0.9 17.39 0.9 17.41 0.9 17.44 0.9 17.46 0.9 17.49 0.9 17.51 0.9 17.54 0.9 17.56 0.9 17.59 0.9 17.61 0.9 17.64 0.9 17.66 0.9 17.69 0.9 17.71 0.9 17.74 0.9 17.76 0.9 17.79 0.9 17.81 0.9 17.84 0.9 17.86 0.9 17.89 0.9 17.91 0.9 17.94 0.9 17.96 0.9 17.99 0.9 18.01 0.9 18.04 0.9 18.06 0.9 18.09 0.9 18.11 0.9 18.14 0.9 18.16 0.9 18.18 0.9 18.21 0.9 18.23 0.9 18.26 0.9 18.28 0.9 18.31 0.9 18.33 0.9 18.36 0.9 18.38 0.9 18.41 0.9 18.43 0.9 18.46 0.9 18.48 0.9 18.51 0.9 18.53 0.9 18.56 0.9 18.58 0.9 18.61 0.9 18.63 0.9 18.66 0.9 18.68 0.9 18.71 0.9 18.73 0.9 18.76 0.9 18.78 0.9 18.81 0.9 18.83 0.9 18.86 0.9 18.88 0.9 18.91 0.9 18.93 0.9 18.96 0.9 18.98 0.9 19.01 0.9 19.03 0.9 19.06 0.9 19.08 0.9 19.11 0.9 19.13 0.9 19.16 0.9 19.18 0.9 19.2 0.9 19.23 0.9 19.25 0.9 19.28 0.9 19.3 0.9 19.33 0.9 19.35 0.9 19.38 0.9 19.4 0.9 19.43 0.9 19.45 0.9 19.48 0.9 19.5 0.9 19.53 0.9 19.55 0.9 19.58 0.9 19.6 0.9 19.63 0.9 19.65 0.9 19.68 0.9 19.7 0.9 19.73 0.9 19.75 0.9 19.78 0.9 19.8 0.9 19.83 0.9 19.85 0.9 19.88 0.9 19.9 0.9 19.93 0.9 19.95 0.9 19.98 0.9 20 0.9 20.03 0.9 20.05 0.9 20.08 0.9 20.1 0.9 20.13 0.9 20.15 0.9 20.18 0.9 20.2 0.9 20.22 0.9 20.25 0.9 20.27 0.9 20.3 0.9 20.32 0.9 20.35 0.9 20.37 0.9 20.4 0.9 20.42 0.9 20.45 0.9 20.47 0.9 20.5 0.9 20.52 0.9 20.55 0.9 20.57 0.9 20.6 0.9 20.62 0.9 20.65 0.9 20.67 0.9 20.7 0.9 20.72 0.9 20.75 0.9 20.77 0.9 20.8 0.9 20.82 0.9 20.85 0.9 20.87 0.9 20.9 0.9 20.92 0.9 20.95 0.9 20.97 0.9 21 0.9 21.02 0.9 21.05 0.9 21.07 0.9 21.1 0.9 21.12 0.9 21.15 0.9 21.17 0.9 21.2 0.9 21.22 0.9 21.24 0.9 21.27 0.9 21.29 0.9 21.32 0.9 21.34 0.9 21.37 0.9 21.39 0.9 21.42 0.9 21.44 0.9 21.47 0.9 21.49 0.9 21.52 0.9 21.54 0.9 21.57 0.9 21.59 0.9 21.62 0.9 21.64 0.9 21.67 0.9 21.69 0.9 21.72 0.9 21.74 0.9 21.77 0.9 21.79 0.9 21.82 0.9 21.84 0.9 21.87 0.9 21.89 0.9 21.92 0.9 21.94 0.9 21.97 0.9 21.99 0.9 22.02 0.9 22.04 0.9 22.07 0.9 22.09 0.9 22.12 0.9 22.14 0.9 22.17 0.9 22.19 0.9 22.21 0.9 22.24 0.9 22.26 0.9 22.29 0.9 22.31 0.9 22.34 0.9 22.36 0.9 22.39 0.9 22.41 0.9 22.44 0.9 22.46 0.9 22.49 0.9 22.51 0.9 22.54 0.9 22.56 0.9 22.59 0.9 22.61 0.9 22.64 0.9 22.66 0.9 22.69 0.9 22.71 0.9 22.74 0.9 22.76 0.9 22.79 0.9 22.81 0.9 22.84 0.9 22.86 0.9 22.89 0.9 22.91 0.9 22.94 0.9 22.96 0.9 22.99 0.9 23.01 0.9 23.04 0.9 23.06 0.9 23.09 0.9 23.11 0.9 23.14 0.9 23.16 0.9 23.19 0.9 23.21 0.9 23.23 0.9 23.26 0.9 23.28 0.9 23.31 0.9 23.33 0.9 23.36 0.9 23.38 0.9 23.41 0.9 23.43 0.9 23.46 0.9 23.48 0.9 23.51 0.9 23.53 0.9 23.56 0.9 23.58 0.9 23.61 0.9 23.63 0.9 23.66 0.9 23.68 0.9 23.71 0.9 23.73 0.9 23.76 0.9 23.78 0.9 23.81 0.9 23.83 0.9 23.86 0.9 23.88 0.9 23.91 0.9 23.93 0.9 23.96 0.9 23.98 0.9 24.01 0.9 24.03 0.9 24.06 0.9 24.08 0.9 24.11 0.9 24.13 0.9 24.16 0.9 24.18 0.9 24.21 0.9 24.23 0.9 24.25 0.9 24.28 0.9 24.3 0.9 24.33 0.9 24.35 0.9 24.38 0.9 24.4 0.9 24.43 0.9 24.45 0.9 24.48 0.9 24.5 0.9 24.53 0.9 24.55 0.9 24.58 0.9 24.6 0.9 24.63 0.9 24.65 0.9 24.68 0.9 24.7 0.9 24.73 0.9 24.75 0.9 24.78 0.9 24.8 0.9 24.83 0.9 24.85 0.9 24.88 0.9 24.9 0.9 24.93 0.9 24.95 0.9 24.98 0.9 25 0.9 25.03 0.9 25.05 0.9 25.08 0.9 25.1 0.9 25.13 0.9 25.15 0.9 25.18 0.9 25.2 0.9 25.23 0.9 25.25 0.9 25.27 0.9 25.3 0.9 25.32 0.9 25.35 0.9 25.37 0.9 25.4 0.9 25.42 0.9 25.45 0.9 25.47 0.9 25.5 0.9 25.52 0.9 25.55 0.9 25.57 0.9 25.6 0.9 25.62 0.9 25.65 0.9 25.67 0.9 25.7 0.9 25.72 0.9 25.75 0.9 25.77 0.9 25.8 0.9 25.82 0.9 25.85 0.9 25.87 0.9 25.9 0.9 25.92 0.9 25.95 0.9 25.97 0.9 26 0.9 26.02 0.9 26.05 0.9 26.07 0.9 26.1 0.9 26.12 0.9 26.15 0.9 26.17 0.9 26.2 0.9 26.22 0.9 26.25 0.9 26.27 0.9 26.29 0.9 26.32 0.9 26.34 0.9 26.37 0.9 26.39 0.9 26.42 0.9 26.44 0.9 26.47 0.9 26.49 0.9 26.52 0.9 26.54 0.9 26.57 0.9 26.59 0.9 26.62 0.9 26.64 0.9 26.67 0.9 26.69 0.9 26.72 0.9 26.74 0.9 26.77 0.9 26.79 0.9 26.82 0.9 26.84 0.9 26.87 0.9 26.89 0.9 26.92 0.9 26.94 0.9 26.97 0.9 26.99 0.9 27.02 0.9 27.04 0.9 27.07 0.9 27.09 0.9 27.12 0.9 27.14 0.9 27.17 0.9 27.19 0.9 27.22 0.9 27.24 0.9 27.26 0.9 27.29 0.9 27.31 0.9 27.34 0.9 27.36 0.9 27.39 0.9 27.41 0.9 27.44 0.9 27.46 0.9 27.49 0.9 27.51 0.9 27.54 0.9 27.56 0.9 27.59 0.9 27.61 0.9 27.64 0.9 27.66 0.9 27.69 0.9 27.71 0.9 27.74 0.9 27.76 0.9 27.79 0.9 27.81 0.9 27.84 0.9 27.86 0.9 27.89 0.9 27.91 0.9 27.94 0.9 27.96 0.9 27.99 0.9 28.01 0.9 28.04 0.9 28.06 0.9 28.09 0.9 28.11 0.9 28.14 0.9 28.16 0.9 28.19 0.9 28.21 0.9 28.24 0.9 28.26 0.9 28.28 0.9 28.31 0.9 28.33 0.9 28.36 0.9 28.38 0.9 28.41 0.9 28.43 0.9 28.46 0.9 28.48 0.9 28.51 0.9 28.53 0.9 28.56 0.9 28.58 0.9 28.61 0.9 28.63 0.9 28.66 0.9 28.68 0.9 28.71 0.9 28.73 0.9 28.76 0.9 28.78 0.9 28.81 0.9 28.83 0.9 28.86 0.9 28.88 0.9 28.91 0.9 28.93 0.9 28.96 0.9 28.98 0.9 29.01 0.9 29.03 0.9 29.06 0.9 29.08 0.9 29.11 0.9 29.13 0.9 29.16 0.9 29.18 0.9 29.21 0.9 29.23 0.9 29.26 0.9 29.28 0.9 29.3 0.9 29.33 0.9 29.35 0.9 29.38 0.9 29.4 0.9 29.43 0.9 29.45 0.9 29.48 0.9 29.5 0.9 29.53 0.9 29.55 0.9 29.58 0.9 29.6 0.9 29.63 0.9 29.65 0.9 29.68 0.9 29.7 0.9 29.73 0.9 29.75 0.9 29.78 0.9 29.8 0.9 29.83 0.9 29.85 0.9 29.88 0.9 29.9 0.9 29.93 0.9 29.95 0.9 29.98 0.9 30 0.9 30.03 0.9 30.05 0.9 30.08 0.9 30.1 0.9 30.13 0.9 30.15 0.9 30.18 0.9 30.2 0.9 30.23 0.9 30.25 0.9 30.28 0.9 30.3 0.9 30.32 0.9 30.35 0.9 30.37 0.9 30.4 0.9 30.42 0.9 30.45 0.9 30.47 0.9 30.5 0.9 30.52 0.9 30.55 0.9 30.57 0.9 30.6 0.9 30.62 0.9 30.65 0.9 30.67 0.9 30.7 0.9 30.72 0.9 30.75 0.9 30.77 0.9 30.8 0.9 30.82 0.9 30.85 0.9 30.87 0.9 30.9 0.9 30.92 0.9 30.95 0.9 30.97 0.9 31 0.9 31.02 0.9 31.05 0.9 31.07 0.9 31.1 0.9 31.12 0.9 31.15 0.9 31.17 0.9 31.2 0.9 31.22 0.9 31.25 0.9 31.27 0.9 31.3 0.9 31.32 0.9 31.34 0.9 31.37 0.9 31.39 0.9 31.42 0.9 31.44 0.9 31.47 0.9 31.49 0.9 31.52 0.9 31.54 0.9 31.57 0.9 31.59 0.9 31.62 0.9 31.64 0.9 31.67 0.9 31.69 0.9 31.72 0.9 31.74 0.9 31.77 0.9 31.79 0.9 31.82 0.9 31.84 0.9 31.87 0.9 31.89 0.9 31.92 0.9 31.94 0.9 31.97 0.9 31.99 0.9 32.02 0.9 32.04 0.9 32.07 0.9 32.09 0.9 32.12 0.9 32.14 0.9 32.17 0.9 32.19 0.9 32.22 0.9 32.24 0.9 32.27 0.9 32.29 0.9 32.31 0.9 32.34 0.9 32.36 0.9 32.39 0.9 32.41 0.9 32.44 0.9 32.46 0.9 32.49 0.9 32.51 0.9 32.54 0.9 32.56 0.9 32.59 0.9 32.61 0.9 32.64 0.9 32.66 0.9 32.69 0.9 32.71 0.9 32.74 0.9 32.76 0.9 32.79 0.9 32.81 0.9 32.84 0.9 32.86 0.9 32.89 0.9 32.91 0.9 32.94 0.9 32.96 0.9 32.99 0.9 33.01 0.9 33.04 0.9 33.06 0.9 33.09 0.9 33.11 0.9 33.14 0.9 33.16 0.9 33.19 0.9 33.21 0.9 33.24 0.9 33.26 0.9 33.29 0.9 33.31 0.9 33.33 0.9 33.36 0.9 33.38 0.9 33.41 0.9 33.43 0.9 33.46 0.9 33.48 0.9 33.51 0.9 33.53 0.9 33.56 0.9 33.58 0.9 33.61 0.9 33.63 0.9 33.66 0.9 33.68 0.9 33.71 0.9 33.73 0.9 33.76 0.9 33.78 0.9 33.81 0.9 33.83 0.9 33.86 0.9 33.88 0.9 33.91 0.9 33.93 0.9 33.96 0.9 33.98 0.9 34.01 0.9 34.03 0.9 34.06 0.9 34.08 0.9 34.11 0.9 34.13 0.9 34.16 0.9 34.18 0.9 34.21 0.9 34.23 0.9 34.26 0.9 34.28 0.9 34.31 0.9 34.33 0.9 34.35 0.9 34.38 0.9 34.4 0.9 34.43 0.9 34.45 0.9 34.48 0.9 34.5 0.9 34.53 0.9 34.55 0.9 34.58 0.9 34.6 0.9 34.63 0.9 34.65 0.9 34.68 0.9 34.7 0.9 34.73 0.9 34.75 0.9 34.78 0.9 34.8 0.9 34.83 0.9 34.85 0.9 34.88 0.9 34.9 0.9 34.93 0.9 34.95 0.9 34.98 0.9 35 0.9 35.03 0.9 35.05 0.9 35.08 0.9 35.1 0.9 35.13 0.9 35.15 0.9 35.18 0.9 35.2 0.9 35.23 0.9 35.25 0.9 35.28 0.9 35.3 0.9 35.33 0.9 35.35 0.9 35.37 0.9 35.4 0.9 35.42 0.9 35.45 0.9 35.47 0.9 35.5 0.9 35.52 0.9 35.55 0.9 35.57 0.9 35.6 0.9 35.62 0.9 35.65 0.9 35.67 0.9 35.7 0.9 35.72 0.9 35.75 0.9 35.77 0.9 35.8 0.9 35.82 0.9 35.85 0.9 35.87 0.9 35.9 0.9 35.92 0.9 35.95 0.9 35.97 0.9 36 0.9 36.02 0.9 36.05 0.9 36.07 0.9 36.1 0.9 36.12 0.9 36.15 0.9 36.17 0.9 36.2 0.9 36.22 0.9 36.25 0.9 36.27 0.9 36.3 0.9 36.32 0.9 36.35 0.9 36.37 0.9 36.39 0.9 36.42 0.9 36.44 0.9 36.47 0.9 36.49 0.9 36.52 0.9 36.54 0.9 36.57 0.9 36.59 0.9 36.62 0.9 36.64 0.9 36.67 0.9 36.69 0.9 36.72 0.9 36.74 0.9 36.77 0.9 36.79 0.9 36.82 0.9 36.84 0.9 36.87 0.9 36.89 0.9 36.92 0.9 36.94 0.9 36.97 0.9 36.99 0.9 37.02 0.9 37.04 0.9 37.07 0.9 37.09 0.9 37.12 0.9 37.14 0.9 37.17 0.9 37.19 0.9 37.22 0.9 37.24 0.9 37.27 0.9 37.29 0.9 37.32 0.9 37.34 0.9 37.36 0.9 37.39 0.9 37.41 0.9 37.44 0.9 37.46 0.9 37.49 0.9 37.51 0.9 37.54 0.9 37.56 0.9 37.59 0.9 37.61 0.9 37.64 0.9 37.66 0.9 37.69 0.9 37.71 0.9 37.74 0.9 37.76 0.9 37.79 0.9 37.81 0.9 37.84 0.9 37.86 0.9 37.89 0.9 37.91 0.9 37.94 0.9 37.96 0.9 37.99 0.9 38.01 0.9 38.04 0.9 38.06 0.9 38.09 0.9 38.11 0.9 38.14 0.9 38.16 0.9 38.19 0.9 38.21 0.9 38.24 0.9 38.26 0.9 38.29 0.9 38.31 0.9 38.34 0.9 38.36 0.9 38.38 0.9 38.41 0.9 38.43 0.9 38.46 0.9 38.48 0.9 38.51 0.9 38.53 0.9 38.56 0.9 38.58 0.9 38.61 0.9 38.63 0.9 38.66 0.9 38.68 0.9 38.71 0.9 38.73 0.9 38.76 0.9 38.78 0.9 38.81 0.9 38.83 0.9 38.86 0.9 38.88 0.9 38.91 0.9 38.93 0.9 38.96 0.9 38.98 0.9 39.01 0.9 39.03 0.9 39.06 0.9 39.08 0.9 39.11 0.9 39.13 0.9 39.16 0.9 39.18 0.9 39.21 0.9 39.23 0.9 39.26 0.9 39.28 0.9 39.31 0.9 39.33 0.9 39.36 0.9 39.38 0.9 39.4 0.9 39.43 0.9 39.45 0.9 39.48 0.9 39.5 0.9 39.53 0.9 39.55 0.9 39.58 0.9 39.6 0.9 39.63 0.9 39.65 0.9 39.68 0.9 39.7 0.9 39.73 0.9 39.75 0.9 39.78 0.9 39.8 0.9 39.83 0.9 39.85 0.9 39.88 0.9 39.9 0.9 39.93 0.9 39.95 0.9 39.98 0.9 40 0.9 40.03 0.9 40.05 0.9 40.08 0.9 40.1 0.9 40.13 0.9 40.15 0.9 40.18 0.9 40.2 0.9 40.23 0.9 40.25 0.9 40.28 0.9 40.3 0.9 40.33 0.9 40.35 0.9 40.38 0.9 40.4 0.9 40.42 0.9 40.45 0.9 40.47 0.9 40.5 0.9 40.52 0.9 40.55 0.9 40.57 0.9 40.6 0.9 40.62 0.9 40.65 0.9 40.67 0.9 40.7 0.9 40.72 0.9 40.75 0.9 40.77 0.9 40.8 0.9 40.82 0.9 40.85 0.9 40.87 0.9 40.9 0.9 40.92 0.9 40.95 0.9 40.97 0.9 41 0.9 41.02 0.9 41.05 0.9 41.07 0.9 41.1 0.9 41.12 0.9 41.15 0.9 41.17 0.9 41.2 0.9 41.22 0.9 41.25 0.9 41.27 0.9 41.3 0.9 41.32 0.9 41.35 0.9 41.37 0.9 41.39 0.9 41.42 0.9 41.44 0.9 41.47 0.9 41.49 0.9 41.52 0.9 41.54 0.9 41.57 0.9 41.59 0.9 41.62 0.9 41.64 0.9 41.67 0.9 41.69 0.9 41.72 0.9 41.74 0.9 41.77 0.9 41.79 0.9 41.82 0.9 41.84 0.9 41.87 0.9 41.89 0.9 41.92 0.9 41.94 0.9 41.97 0.9 41.99 0.9 42.02 0.9 42.04 0.9 42.07 0.9 42.09 0.9 42.12 0.9 42.14 0.9 42.17 0.9 42.19 0.9 42.22 0.9 42.24 0.9 42.27 0.9 42.29 0.9 42.32 0.9 42.34 0.9 42.37 0.9 42.39 0.9 42.41 0.9 42.44 0.9 42.46 0.9 42.49 0.9 42.51 0.9 42.54 0.9 42.56 0.9 42.59 0.9 42.61 0.9 42.64 0.9 42.66 0.9 42.69 0.9 42.71 0.9 42.74 0.9 42.76 0.9 42.79 0.9 42.81 0.9 42.84 0.9 42.86 0.9 42.89 0.9 42.91 0.9 42.94 0.9 42.96 0.9 42.99 0.9 43.01 0.9 43.04 0.9 43.06 0.9 43.09 0.9 43.11 0.9 43.14 0.9 43.16 0.9 43.19 0.9 43.21 0.9 43.24 0.9 43.26 0.9 43.29 0.9 43.31 0.9 43.34 0.9 43.36 0.9 43.39 0.9 43.41 0.9 43.43 0.9 43.46 0.9 43.48 0.9 43.51 0.9 43.53 0.9 43.56 0.9 43.58 0.9 43.61 0.9 43.63 0.9 43.66 0.9 43.68 0.9 43.71 0.9 43.73 0.9 43.76 0.9 43.78 0.9 43.81 0.9 43.83 0.9 43.86 0.9 43.88 0.9 43.91 0.9 43.93 0.9 43.96 0.9 43.98 0.9 44.01 0.9 44.03 0.9 44.06 0.9 44.08 0.9 44.11 0.9 44.13 0.9 44.16 0.9 44.18 0.9 44.21 0.9 44.23 0.9 44.26 0.9 44.28 0.9 44.31 0.9 44.33 0.9 44.36 0.9 44.38 0.9 44.41 0.9 44.43 0.9 44.45 0.9 44.48 0.9 44.5 0.9 44.53 0.9 44.55 0.9 44.58 0.9 44.6 0.9 44.63 0.9 44.65 0.9 44.68 0.9 44.7 0.9 44.73 0.9 44.75 0.9 44.78 0.9 44.8 0.9 44.83 0.9 44.85 0.9 44.88 0.9 44.9 0.9 44.93 0.9 44.95 0.9 44.98 0.9 45 0.9 45.03 0.9 45.05 0.9 45.08 0.9 45.1 0.9 45.13 0.9 45.15 0.9 45.18 0.9 45.2 0.9 45.23 0.9 45.25 0.9 45.28 0.9 45.3 0.9 45.33 0.9 45.35 0.9 45.38 0.9 45.4 0.9" class="primitive"/>
          </g>
          <g transform="translate(67.03,28.49)" id="img-8311e7c5-747" class="geometry color_E_H" stroke="#D4CA3A">
            <path fill="none" d="M-45.4,0.23 L -45.38 0.23 -45.35 0.23 -45.33 0.23 -45.3 0.23 -45.28 0.23 -45.25 0.23 -45.23 0.23 -45.2 0.23 -45.18 0.23 -45.15 0.23 -45.13 0.22 -45.1 0.22 -45.08 0.22 -45.05 0.22 -45.03 0.22 -45 0.22 -44.98 0.22 -44.95 0.22 -44.93 0.22 -44.9 0.22 -44.88 0.22 -44.85 0.22 -44.83 0.22 -44.8 0.22 -44.78 0.22 -44.75 0.22 -44.73 0.21 -44.7 0.21 -44.68 0.21 -44.65 0.21 -44.63 0.21 -44.6 0.21 -44.58 0.21 -44.55 0.21 -44.53 0.2 -44.5 0.2 -44.48 0.2 -44.45 0.2 -44.43 0.2 -44.41 0.2 -44.38 0.19 -44.36 0.19 -44.33 0.19 -44.31 0.19 -44.28 0.19 -44.26 0.18 -44.23 0.18 -44.21 0.18 -44.18 0.18 -44.16 0.18 -44.13 0.17 -44.11 0.17 -44.08 0.17 -44.06 0.17 -44.03 0.16 -44.01 0.16 -43.98 0.16 -43.96 0.15 -43.93 0.15 -43.91 0.15 -43.88 0.15 -43.86 0.14 -43.83 0.14 -43.81 0.14 -43.78 0.13 -43.76 0.13 -43.73 0.13 -43.71 0.13 -43.68 0.12 -43.66 0.12 -43.63 0.12 -43.61 0.11 -43.58 0.11 -43.56 0.11 -43.53 0.1 -43.51 0.1 -43.48 0.1 -43.46 0.09 -43.43 0.09 -43.41 0.08 -43.39 0.08 -43.36 0.08 -43.34 0.07 -43.31 0.07 -43.29 0.06 -43.26 0.06 -43.24 0.06 -43.21 0.05 -43.19 0.05 -43.16 0.04 -43.14 0.04 -43.11 0.04 -43.09 0.03 -43.06 0.03 -43.04 0.02 -43.01 0.02 -42.99 0.01 -42.96 0.01 -42.94 0.01 -42.91 0 -42.89 -0 -42.86 -0.01 -42.84 -0.01 -42.81 -0.02 -42.79 -0.02 -42.76 -0.03 -42.74 -0.03 -42.71 -0.04 -42.69 -0.04 -42.66 -0.05 -42.64 -0.06 -42.61 -0.06 -42.59 -0.07 -42.56 -0.07 -42.54 -0.08 -42.51 -0.08 -42.49 -0.09 -42.46 -0.1 -42.44 -0.1 -42.41 -0.11 -42.39 -0.11 -42.37 -0.12 -42.34 -0.13 -42.32 -0.13 -42.29 -0.14 -42.27 -0.15 -42.24 -0.15 -42.22 -0.16 -42.19 -0.17 -42.17 -0.17 -42.14 -0.18 -42.12 -0.19 -42.09 -0.2 -42.07 -0.2 -42.04 -0.21 -42.02 -0.22 -41.99 -0.23 -41.97 -0.23 -41.94 -0.24 -41.92 -0.25 -41.89 -0.26 -41.87 -0.27 -41.84 -0.28 -41.82 -0.28 -41.79 -0.29 -41.77 -0.3 -41.74 -0.31 -41.72 -0.32 -41.69 -0.33 -41.67 -0.34 -41.64 -0.35 -41.62 -0.36 -41.59 -0.37 -41.57 -0.38 -41.54 -0.39 -41.52 -0.4 -41.49 -0.41 -41.47 -0.42 -41.44 -0.43 -41.42 -0.44 -41.39 -0.45 -41.37 -0.47 -41.35 -0.48 -41.32 -0.49 -41.3 -0.5 -41.27 -0.51 -41.25 -0.52 -41.22 -0.54 -41.2 -0.55 -41.17 -0.56 -41.15 -0.57 -41.12 -0.59 -41.1 -0.6 -41.07 -0.61 -41.05 -0.63 -41.02 -0.64 -41 -0.66 -40.97 -0.67 -40.95 -0.68 -40.92 -0.7 -40.9 -0.71 -40.87 -0.73 -40.85 -0.74 -40.82 -0.76 -40.8 -0.77 -40.77 -0.79 -40.75 -0.8 -40.72 -0.82 -40.7 -0.84 -40.67 -0.85 -40.65 -0.87 -40.62 -0.89 -40.6 -0.9 -40.57 -0.92 -40.55 -0.94 -40.52 -0.95 -40.5 -0.97 -40.47 -0.99 -40.45 -1.01 -40.42 -1.02 -40.4 -1.04 -40.38 -1.06 -40.35 -1.08 -40.33 -1.1 -40.3 -1.12 -40.28 -1.14 -40.25 -1.15 -40.23 -1.17 -40.2 -1.19 -40.18 -1.21 -40.15 -1.23 -40.13 -1.25 -40.1 -1.27 -40.08 -1.29 -40.05 -1.31 -40.03 -1.33 -40 -1.35 -39.98 -1.37 -39.95 -1.39 -39.93 -1.41 -39.9 -1.43 -39.88 -1.45 -39.85 -1.47 -39.83 -1.49 -39.8 -1.51 -39.78 -1.53 -39.75 -1.56 -39.73 -1.58 -39.7 -1.6 -39.68 -1.62 -39.65 -1.64 -39.63 -1.66 -39.6 -1.68 -39.58 -1.7 -39.55 -1.72 -39.53 -1.74 -39.5 -1.76 -39.48 -1.78 -39.45 -1.8 -39.43 -1.82 -39.4 -1.84 -39.38 -1.86 -39.36 -1.88 -39.33 -1.9 -39.31 -1.92 -39.28 -1.94 -39.26 -1.96 -39.23 -1.98 -39.21 -1.99 -39.18 -2.01 -39.16 -2.03 -39.13 -2.05 -39.11 -2.06 -39.08 -2.08 -39.06 -2.1 -39.03 -2.11 -39.01 -2.13 -38.98 -2.15 -38.96 -2.16 -38.93 -2.18 -38.91 -2.19 -38.88 -2.2 -38.86 -2.22 -38.83 -2.23 -38.81 -2.24 -38.78 -2.26 -38.76 -2.27 -38.73 -2.28 -38.71 -2.29 -38.68 -2.3 -38.66 -2.31 -38.63 -2.32 -38.61 -2.33 -38.58 -2.34 -38.56 -2.34 -38.53 -2.35 -38.51 -2.36 -38.48 -2.36 -38.46 -2.37 -38.43 -2.37 -38.41 -2.38 -38.38 -2.38 -38.36 -2.38 -38.34 -2.38 -38.31 -2.39 -38.29 -2.39 -38.26 -2.39 -38.24 -2.39 -38.21 -2.38 -38.19 -2.38 -38.16 -2.38 -38.14 -2.37 -38.11 -2.37 -38.09 -2.37 -38.06 -2.36 -38.04 -2.35 -38.01 -2.35 -37.99 -2.34 -37.96 -2.33 -37.94 -2.32 -37.91 -2.31 -37.89 -2.3 -37.86 -2.29 -37.84 -2.28 -37.81 -2.26 -37.79 -2.25 -37.76 -2.24 -37.74 -2.22 -37.71 -2.21 -37.69 -2.19 -37.66 -2.18 -37.64 -2.16 -37.61 -2.14 -37.59 -2.12 -37.56 -2.11 -37.54 -2.09 -37.51 -2.07 -37.49 -2.05 -37.46 -2.03 -37.44 -2 -37.41 -1.98 -37.39 -1.96 -37.36 -1.94 -37.34 -1.92 -37.32 -1.89 -37.29 -1.87 -37.27 -1.85 -37.24 -1.82 -37.22 -1.8 -37.19 -1.77 -37.17 -1.75 -37.14 -1.72 -37.12 -1.7 -37.09 -1.67 -37.07 -1.65 -37.04 -1.62 -37.02 -1.59 -36.99 -1.57 -36.97 -1.54 -36.94 -1.51 -36.92 -1.49 -36.89 -1.46 -36.87 -1.43 -36.84 -1.41 -36.82 -1.38 -36.79 -1.36 -36.77 -1.33 -36.74 -1.3 -36.72 -1.28 -36.69 -1.25 -36.67 -1.22 -36.64 -1.2 -36.62 -1.17 -36.59 -1.15 -36.57 -1.12 -36.54 -1.1 -36.52 -1.07 -36.49 -1.05 -36.47 -1.02 -36.44 -1 -36.42 -0.97 -36.39 -0.95 -36.37 -0.93 -36.35 -0.9 -36.32 -0.88 -36.3 -0.86 -36.27 -0.83 -36.25 -0.81 -36.22 -0.79 -36.2 -0.77 -36.17 -0.75 -36.15 -0.73 -36.12 -0.71 -36.1 -0.69 -36.07 -0.67 -36.05 -0.65 -36.02 -0.63 -36 -0.61 -35.97 -0.59 -35.95 -0.57 -35.92 -0.55 -35.9 -0.54 -35.87 -0.52 -35.85 -0.5 -35.82 -0.49 -35.8 -0.47 -35.77 -0.46 -35.75 -0.44 -35.72 -0.43 -35.7 -0.41 -35.67 -0.4 -35.65 -0.38 -35.62 -0.37 -35.6 -0.36 -35.57 -0.34 -35.55 -0.33 -35.52 -0.32 -35.5 -0.31 -35.47 -0.29 -35.45 -0.28 -35.42 -0.27 -35.4 -0.26 -35.37 -0.25 -35.35 -0.24 -35.33 -0.23 -35.3 -0.22 -35.28 -0.21 -35.25 -0.2 -35.23 -0.19 -35.2 -0.18 -35.18 -0.17 -35.15 -0.16 -35.13 -0.16 -35.1 -0.15 -35.08 -0.14 -35.05 -0.13 -35.03 -0.13 -35 -0.12 -34.98 -0.11 -34.95 -0.11 -34.93 -0.1 -34.9 -0.09 -34.88 -0.09 -34.85 -0.08 -34.83 -0.08 -34.8 -0.07 -34.78 -0.06 -34.75 -0.06 -34.73 -0.05 -34.7 -0.05 -34.68 -0.04 -34.65 -0.04 -34.63 -0.04 -34.6 -0.03 -34.58 -0.03 -34.55 -0.02 -34.53 -0.02 -34.5 -0.01 -34.48 -0.01 -34.45 -0.01 -34.43 -0 -34.4 -0 -34.38 0 -34.35 0.01 -34.33 0.01 -34.31 0.01 -34.28 0.02 -34.26 0.02 -34.23 0.02 -34.21 0.02 -34.18 0.03 -34.16 0.03 -34.13 0.03 -34.11 0.03 -34.08 0.04 -34.06 0.04 -34.03 0.04 -34.01 0.04 -33.98 0.04 -33.96 0.05 -33.93 0.05 -33.91 0.05 -33.88 0.05 -33.86 0.05 -33.83 0.05 -33.81 0.06 -33.78 0.06 -33.76 0.06 -33.73 0.06 -33.71 0.06 -33.68 0.06 -33.66 0.06 -33.63 0.07 -33.61 0.07 -33.58 0.07 -33.56 0.07 -33.53 0.07 -33.51 0.07 -33.48 0.07 -33.46 0.07 -33.43 0.07 -33.41 0.08 -33.38 0.08 -33.36 0.08 -33.33 0.08 -33.31 0.08 -33.29 0.08 -33.26 0.08 -33.24 0.08 -33.21 0.08 -33.19 0.08 -33.16 0.08 -33.14 0.08 -33.11 0.09 -33.09 0.09 -33.06 0.09 -33.04 0.09 -33.01 0.09 -32.99 0.09 -32.96 0.09 -32.94 0.09 -32.91 0.09 -32.89 0.09 -32.86 0.09 -32.84 0.09 -32.81 0.09 -32.79 0.09 -32.76 0.09 -32.74 0.09 -32.71 0.09 -32.69 0.09 -32.66 0.09 -32.64 0.09 -32.61 0.09 -32.59 0.1 -32.56 0.1 -32.54 0.1 -32.51 0.1 -32.49 0.1 -32.46 0.1 -32.44 0.1 -32.41 0.1 -32.39 0.1 -32.36 0.1 -32.34 0.1 -32.31 0.1 -32.29 0.1 -32.27 0.1 -32.24 0.1 -32.22 0.1 -32.19 0.1 -32.17 0.1 -32.14 0.1 -32.12 0.1 -32.09 0.1 -32.07 0.1 -32.04 0.1 -32.02 0.1 -31.99 0.1 -31.97 0.1 -31.94 0.1 -31.92 0.1 -31.89 0.1 -31.87 0.1 -31.84 0.1 -31.82 0.1 -31.79 0.1 -31.77 0.1 -31.74 0.1 -31.72 0.1 -31.69 0.1 -31.67 0.1 -31.64 0.1 -31.62 0.1 -31.59 0.1 -31.57 0.1 -31.54 0.1 -31.52 0.1 -31.49 0.1 -31.47 0.1 -31.44 0.1 -31.42 0.1 -31.39 0.1 -31.37 0.1 -31.34 0.1 -31.32 0.1 -31.3 0.1 -31.27 0.1 -31.25 0.1 -31.22 0.1 -31.2 0.1 -31.17 0.1 -31.15 0.1 -31.12 0.11 -31.1 0.11 -31.07 0.11 -31.05 0.11 -31.02 0.11 -31 0.11 -30.97 0.11 -30.95 0.11 -30.92 0.11 -30.9 0.11 -30.87 0.11 -30.85 0.11 -30.82 0.11 -30.8 0.11 -30.77 0.11 -30.75 0.11 -30.72 0.11 -30.7 0.11 -30.67 0.11 -30.65 0.11 -30.62 0.11 -30.6 0.11 -30.57 0.11 -30.55 0.11 -30.52 0.11 -30.5 0.11 -30.47 0.11 -30.45 0.11 -30.42 0.11 -30.4 0.11 -30.37 0.11 -30.35 0.11 -30.32 0.11 -30.3 0.11 -30.28 0.11 -30.25 0.11 -30.23 0.11 -30.2 0.11 -30.18 0.11 -30.15 0.11 -30.13 0.11 -30.1 0.11 -30.08 0.11 -30.05 0.11 -30.03 0.11 -30 0.11 -29.98 0.11 -29.95 0.11 -29.93 0.11 -29.9 0.11 -29.88 0.11 -29.85 0.11 -29.83 0.11 -29.8 0.11 -29.78 0.11 -29.75 0.11 -29.73 0.11 -29.7 0.11 -29.68 0.11 -29.65 0.11 -29.63 0.11 -29.6 0.11 -29.58 0.11 -29.55 0.11 -29.53 0.11 -29.5 0.11 -29.48 0.11 -29.45 0.11 -29.43 0.11 -29.4 0.11 -29.38 0.11 -29.35 0.11 -29.33 0.11 -29.3 0.11 -29.28 0.11 -29.26 0.11 -29.23 0.11 -29.21 0.11 -29.18 0.11 -29.16 0.11 -29.13 0.11 -29.11 0.11 -29.08 0.11 -29.06 0.11 -29.03 0.11 -29.01 0.11 -28.98 0.11 -28.96 0.11 -28.93 0.11 -28.91 0.11 -28.88 0.11 -28.86 0.11 -28.83 0.11 -28.81 0.11 -28.78 0.11 -28.76 0.11 -28.73 0.11 -28.71 0.11 -28.68 0.11 -28.66 0.11 -28.63 0.11 -28.61 0.11 -28.58 0.11 -28.56 0.11 -28.53 0.11 -28.51 0.11 -28.48 0.11 -28.46 0.11 -28.43 0.11 -28.41 0.11 -28.38 0.11 -28.36 0.11 -28.33 0.11 -28.31 0.11 -28.28 0.11 -28.26 0.11 -28.24 0.11 -28.21 0.11 -28.19 0.11 -28.16 0.11 -28.14 0.11 -28.11 0.11 -28.09 0.11 -28.06 0.11 -28.04 0.11 -28.01 0.11 -27.99 0.11 -27.96 0.11 -27.94 0.11 -27.91 0.11 -27.89 0.11 -27.86 0.11 -27.84 0.11 -27.81 0.11 -27.79 0.11 -27.76 0.11 -27.74 0.11 -27.71 0.11 -27.69 0.11 -27.66 0.11 -27.64 0.11 -27.61 0.11 -27.59 0.11 -27.56 0.11 -27.54 0.11 -27.51 0.11 -27.49 0.11 -27.46 0.11 -27.44 0.11 -27.41 0.11 -27.39 0.11 -27.36 0.11 -27.34 0.11 -27.31 0.11 -27.29 0.11 -27.26 0.11 -27.24 0.11 -27.22 0.11 -27.19 0.11 -27.17 0.11 -27.14 0.11 -27.12 0.11 -27.09 0.11 -27.07 0.11 -27.04 0.11 -27.02 0.11 -26.99 0.11 -26.97 0.11 -26.94 0.11 -26.92 0.11 -26.89 0.11 -26.87 0.11 -26.84 0.11 -26.82 0.11 -26.79 0.11 -26.77 0.11 -26.74 0.11 -26.72 0.11 -26.69 0.11 -26.67 0.11 -26.64 0.11 -26.62 0.11 -26.59 0.11 -26.57 0.11 -26.54 0.11 -26.52 0.11 -26.49 0.11 -26.47 0.11 -26.44 0.11 -26.42 0.11 -26.39 0.11 -26.37 0.11 -26.34 0.11 -26.32 0.11 -26.29 0.11 -26.27 0.11 -26.25 0.11 -26.22 0.11 -26.2 0.11 -26.17 0.11 -26.15 0.11 -26.12 0.11 -26.1 0.11 -26.07 0.11 -26.05 0.11 -26.02 0.11 -26 0.11 -25.97 0.11 -25.95 0.11 -25.92 0.11 -25.9 0.11 -25.87 0.11 -25.85 0.11 -25.82 0.11 -25.8 0.11 -25.77 0.11 -25.75 0.11 -25.72 0.11 -25.7 0.11 -25.67 0.11 -25.65 0.11 -25.62 0.11 -25.6 0.11 -25.57 0.11 -25.55 0.11 -25.52 0.11 -25.5 0.11 -25.47 0.11 -25.45 0.11 -25.42 0.11 -25.4 0.11 -25.37 0.11 -25.35 0.11 -25.32 0.11 -25.3 0.11 -25.27 0.11 -25.25 0.11 -25.23 0.11 -25.2 0.11 -25.18 0.11 -25.15 0.11 -25.13 0.11 -25.1 0.11 -25.08 0.11 -25.05 0.11 -25.03 0.11 -25 0.11 -24.98 0.11 -24.95 0.11 -24.93 0.11 -24.9 0.11 -24.88 0.11 -24.85 0.11 -24.83 0.11 -24.8 0.11 -24.78 0.11 -24.75 0.11 -24.73 0.11 -24.7 0.11 -24.68 0.11 -24.65 0.11 -24.63 0.11 -24.6 0.11 -24.58 0.11 -24.55 0.11 -24.53 0.11 -24.5 0.11 -24.48 0.11 -24.45 0.11 -24.43 0.11 -24.4 0.11 -24.38 0.11 -24.35 0.11 -24.33 0.11 -24.3 0.11 -24.28 0.11 -24.25 0.11 -24.23 0.11 -24.21 0.11 -24.18 0.11 -24.16 0.11 -24.13 0.11 -24.11 0.11 -24.08 0.11 -24.06 0.11 -24.03 0.11 -24.01 0.11 -23.98 0.11 -23.96 0.11 -23.93 0.11 -23.91 0.11 -23.88 0.11 -23.86 0.11 -23.83 0.11 -23.81 0.11 -23.78 0.11 -23.76 0.11 -23.73 0.11 -23.71 0.11 -23.68 0.11 -23.66 0.11 -23.63 0.11 -23.61 0.11 -23.58 0.11 -23.56 0.11 -23.53 0.11 -23.51 0.11 -23.48 0.11 -23.46 0.11 -23.43 0.11 -23.41 0.11 -23.38 0.11 -23.36 0.11 -23.33 0.11 -23.31 0.11 -23.28 0.11 -23.26 0.11 -23.23 0.11 -23.21 0.11 -23.19 0.11 -23.16 0.11 -23.14 0.11 -23.11 0.11 -23.09 0.11 -23.06 0.11 -23.04 0.11 -23.01 0.11 -22.99 0.11 -22.96 0.11 -22.94 0.11 -22.91 0.11 -22.89 0.11 -22.86 0.11 -22.84 0.11 -22.81 0.11 -22.79 0.11 -22.76 0.11 -22.74 0.11 -22.71 0.11 -22.69 0.11 -22.66 0.11 -22.64 0.11 -22.61 0.11 -22.59 0.11 -22.56 0.11 -22.54 0.11 -22.51 0.11 -22.49 0.11 -22.46 0.11 -22.44 0.11 -22.41 0.11 -22.39 0.11 -22.36 0.11 -22.34 0.11 -22.31 0.11 -22.29 0.11 -22.26 0.11 -22.24 0.11 -22.21 0.11 -22.19 0.11 -22.17 0.11 -22.14 0.11 -22.12 0.11 -22.09 0.11 -22.07 0.11 -22.04 0.11 -22.02 0.11 -21.99 0.11 -21.97 0.11 -21.94 0.11 -21.92 0.11 -21.89 0.11 -21.87 0.11 -21.84 0.11 -21.82 0.11 -21.79 0.11 -21.77 0.11 -21.74 0.11 -21.72 0.11 -21.69 0.11 -21.67 0.11 -21.64 0.11 -21.62 0.11 -21.59 0.11 -21.57 0.11 -21.54 0.11 -21.52 0.11 -21.49 0.11 -21.47 0.11 -21.44 0.11 -21.42 0.11 -21.39 0.11 -21.37 0.11 -21.34 0.11 -21.32 0.11 -21.29 0.11 -21.27 0.11 -21.24 0.11 -21.22 0.11 -21.2 0.11 -21.17 0.11 -21.15 0.11 -21.12 0.11 -21.1 0.11 -21.07 0.11 -21.05 0.11 -21.02 0.11 -21 0.11 -20.97 0.11 -20.95 0.11 -20.92 0.11 -20.9 0.11 -20.87 0.11 -20.85 0.11 -20.82 0.11 -20.8 0.11 -20.77 0.11 -20.75 0.11 -20.72 0.11 -20.7 0.11 -20.67 0.11 -20.65 0.11 -20.62 0.11 -20.6 0.11 -20.57 0.11 -20.55 0.11 -20.52 0.11 -20.5 0.11 -20.47 0.11 -20.45 0.11 -20.42 0.11 -20.4 0.11 -20.37 0.11 -20.35 0.11 -20.32 0.11 -20.3 0.11 -20.27 0.11 -20.25 0.11 -20.22 0.11 -20.2 0.11 -20.18 0.11 -20.15 0.11 -20.13 0.11 -20.1 0.11 -20.08 0.11 -20.05 0.11 -20.03 0.11 -20 0.11 -19.98 0.11 -19.95 0.11 -19.93 0.11 -19.9 0.11 -19.88 0.11 -19.85 0.11 -19.83 0.11 -19.8 0.11 -19.78 0.11 -19.75 0.11 -19.73 0.11 -19.7 0.11 -19.68 0.11 -19.65 0.11 -19.63 0.11 -19.6 0.11 -19.58 0.11 -19.55 0.11 -19.53 0.11 -19.5 0.11 -19.48 0.11 -19.45 0.11 -19.43 0.11 -19.4 0.11 -19.38 0.11 -19.35 0.11 -19.33 0.11 -19.3 0.11 -19.28 0.11 -19.25 0.11 -19.23 0.11 -19.2 0.11 -19.18 0.11 -19.16 0.11 -19.13 0.11 -19.11 0.11 -19.08 0.11 -19.06 0.11 -19.03 0.11 -19.01 0.11 -18.98 0.11 -18.96 0.11 -18.93 0.11 -18.91 0.11 -18.88 0.11 -18.86 0.11 -18.83 0.11 -18.81 0.11 -18.78 0.11 -18.76 0.11 -18.73 0.11 -18.71 0.11 -18.68 0.11 -18.66 0.11 -18.63 0.11 -18.61 0.11 -18.58 0.11 -18.56 0.11 -18.53 0.11 -18.51 0.11 -18.48 0.11 -18.46 0.11 -18.43 0.11 -18.41 0.11 -18.38 0.11 -18.36 0.11 -18.33 0.11 -18.31 0.11 -18.28 0.11 -18.26 0.11 -18.23 0.11 -18.21 0.11 -18.18 0.11 -18.16 0.11 -18.14 0.11 -18.11 0.11 -18.09 0.11 -18.06 0.11 -18.04 0.11 -18.01 0.11 -17.99 0.11 -17.96 0.11 -17.94 0.11 -17.91 0.11 -17.89 0.11 -17.86 0.11 -17.84 0.11 -17.81 0.11 -17.79 0.11 -17.76 0.11 -17.74 0.11 -17.71 0.11 -17.69 0.11 -17.66 0.11 -17.64 0.11 -17.61 0.11 -17.59 0.11 -17.56 0.11 -17.54 0.11 -17.51 0.11 -17.49 0.11 -17.46 0.11 -17.44 0.11 -17.41 0.11 -17.39 0.11 -17.36 0.11 -17.34 0.11 -17.31 0.11 -17.29 0.11 -17.26 0.11 -17.24 0.11 -17.21 0.11 -17.19 0.11 -17.16 0.11 -17.14 0.11 -17.12 0.11 -17.09 0.11 -17.07 0.11 -17.04 0.11 -17.02 0.11 -16.99 0.11 -16.97 0.11 -16.94 0.11 -16.92 0.11 -16.89 0.11 -16.87 0.11 -16.84 0.11 -16.82 0.11 -16.79 0.11 -16.77 0.11 -16.74 0.11 -16.72 0.11 -16.69 0.11 -16.67 0.11 -16.64 0.11 -16.62 0.11 -16.59 0.11 -16.57 0.11 -16.54 0.11 -16.52 0.11 -16.49 0.11 -16.47 0.11 -16.44 0.11 -16.42 0.11 -16.39 0.11 -16.37 0.11 -16.34 0.11 -16.32 0.11 -16.29 0.11 -16.27 0.11 -16.24 0.11 -16.22 0.11 -16.19 0.11 -16.17 0.11 -16.15 0.11 -16.12 0.11 -16.1 0.11 -16.07 0.11 -16.05 0.11 -16.02 0.11 -16 0.11 -15.97 0.11 -15.95 0.11 -15.92 0.11 -15.9 0.11 -15.87 0.11 -15.85 0.11 -15.82 0.11 -15.8 0.11 -15.77 0.11 -15.75 0.11 -15.72 0.11 -15.7 0.11 -15.67 0.11 -15.65 0.11 -15.62 0.11 -15.6 0.11 -15.57 0.11 -15.55 0.11 -15.52 0.11 -15.5 0.11 -15.47 0.11 -15.45 0.11 -15.42 0.11 -15.4 0.11 -15.37 0.11 -15.35 0.11 -15.32 0.11 -15.3 0.11 -15.27 0.11 -15.25 0.11 -15.22 0.11 -15.2 0.11 -15.17 0.11 -15.15 0.11 -15.13 0.11 -15.1 0.11 -15.08 0.11 -15.05 0.11 -15.03 0.11 -15 0.11 -14.98 0.11 -14.95 0.11 -14.93 0.11 -14.9 0.11 -14.88 0.11 -14.85 0.11 -14.83 0.11 -14.8 0.11 -14.78 0.11 -14.75 0.11 -14.73 0.11 -14.7 0.11 -14.68 0.11 -14.65 0.11 -14.63 0.11 -14.6 0.11 -14.58 0.11 -14.55 0.11 -14.53 0.11 -14.5 0.11 -14.48 0.11 -14.45 0.11 -14.43 0.11 -14.4 0.11 -14.38 0.11 -14.35 0.11 -14.33 0.11 -14.3 0.11 -14.28 0.11 -14.25 0.11 -14.23 0.11 -14.2 0.11 -14.18 0.11 -14.15 0.11 -14.13 0.11 -14.11 0.11 -14.08 0.11 -14.06 0.11 -14.03 0.11 -14.01 0.11 -13.98 0.11 -13.96 0.11 -13.93 0.11 -13.91 0.11 -13.88 0.11 -13.86 0.11 -13.83 0.11 -13.81 0.11 -13.78 0.11 -13.76 0.11 -13.73 0.11 -13.71 0.11 -13.68 0.11 -13.66 0.11 -13.63 0.11 -13.61 0.11 -13.58 0.11 -13.56 0.11 -13.53 0.11 -13.51 0.11 -13.48 0.11 -13.46 0.11 -13.43 0.11 -13.41 0.11 -13.38 0.11 -13.36 0.11 -13.33 0.11 -13.31 0.11 -13.28 0.11 -13.26 0.11 -13.23 0.11 -13.21 0.11 -13.18 0.11 -13.16 0.11 -13.13 0.11 -13.11 0.11 -13.09 0.11 -13.06 0.11 -13.04 0.11 -13.01 0.11 -12.99 0.11 -12.96 0.11 -12.94 0.11 -12.91 0.11 -12.89 0.11 -12.86 0.11 -12.84 0.11 -12.81 0.11 -12.79 0.11 -12.76 0.11 -12.74 0.11 -12.71 0.11 -12.69 0.11 -12.66 0.11 -12.64 0.11 -12.61 0.11 -12.59 0.11 -12.56 0.11 -12.54 0.11 -12.51 0.11 -12.49 0.11 -12.46 0.11 -12.44 0.11 -12.41 0.11 -12.39 0.11 -12.36 0.11 -12.34 0.11 -12.31 0.11 -12.29 0.11 -12.26 0.11 -12.24 0.11 -12.21 0.11 -12.19 0.11 -12.16 0.11 -12.14 0.11 -12.12 0.11 -12.09 0.11 -12.07 0.11 -12.04 0.11 -12.02 0.11 -11.99 0.11 -11.97 0.11 -11.94 0.11 -11.92 0.11 -11.89 0.11 -11.87 0.11 -11.84 0.11 -11.82 0.11 -11.79 0.11 -11.77 0.11 -11.74 0.11 -11.72 0.11 -11.69 0.11 -11.67 0.11 -11.64 0.11 -11.62 0.11 -11.59 0.11 -11.57 0.11 -11.54 0.11 -11.52 0.11 -11.49 0.11 -11.47 0.11 -11.44 0.11 -11.42 0.11 -11.39 0.11 -11.37 0.11 -11.34 0.11 -11.32 0.11 -11.29 0.11 -11.27 0.11 -11.24 0.11 -11.22 0.11 -11.19 0.11 -11.17 0.11 -11.14 0.11 -11.12 0.11 -11.1 0.11 -11.07 0.11 -11.05 0.11 -11.02 0.11 -11 0.11 -10.97 0.11 -10.95 0.11 -10.92 0.11 -10.9 0.11 -10.87 0.11 -10.85 0.11 -10.82 0.11 -10.8 0.11 -10.77 0.11 -10.75 0.11 -10.72 0.11 -10.7 0.11 -10.67 0.11 -10.65 0.11 -10.62 0.11 -10.6 0.11 -10.57 0.11 -10.55 0.11 -10.52 0.11 -10.5 0.11 -10.47 0.11 -10.45 0.11 -10.42 0.11 -10.4 0.11 -10.37 0.11 -10.35 0.11 -10.32 0.11 -10.3 0.11 -10.27 0.11 -10.25 0.11 -10.22 0.11 -10.2 0.11 -10.17 0.11 -10.15 0.11 -10.12 0.11 -10.1 0.11 -10.08 0.11 -10.05 0.11 -10.03 0.11 -10 0.11 -9.98 0.11 -9.95 0.11 -9.93 0.11 -9.9 0.11 -9.88 0.11 -9.85 0.11 -9.83 0.11 -9.8 0.11 -9.78 0.11 -9.75 0.11 -9.73 0.11 -9.7 0.11 -9.68 0.11 -9.65 0.11 -9.63 0.11 -9.6 0.11 -9.58 0.11 -9.55 0.11 -9.53 0.11 -9.5 0.11 -9.48 0.11 -9.45 0.11 -9.43 0.11 -9.4 0.11 -9.38 0.11 -9.35 0.11 -9.33 0.11 -9.3 0.11 -9.28 0.11 -9.25 0.11 -9.23 0.11 -9.2 0.11 -9.18 0.11 -9.15 0.11 -9.13 0.11 -9.1 0.11 -9.08 0.11 -9.06 0.11 -9.03 0.11 -9.01 0.11 -8.98 0.11 -8.96 0.11 -8.93 0.11 -8.91 0.11 -8.88 0.11 -8.86 0.11 -8.83 0.11 -8.81 0.11 -8.78 0.11 -8.76 0.11 -8.73 0.11 -8.71 0.11 -8.68 0.11 -8.66 0.11 -8.63 0.11 -8.61 0.11 -8.58 0.11 -8.56 0.11 -8.53 0.11 -8.51 0.11 -8.48 0.11 -8.46 0.11 -8.43 0.11 -8.41 0.11 -8.38 0.11 -8.36 0.11 -8.33 0.11 -8.31 0.11 -8.28 0.11 -8.26 0.11 -8.23 0.11 -8.21 0.11 -8.18 0.11 -8.16 0.11 -8.13 0.11 -8.11 0.11 -8.08 0.11 -8.06 0.11 -8.04 0.11 -8.01 0.11 -7.99 0.11 -7.96 0.11 -7.94 0.11 -7.91 0.11 -7.89 0.11 -7.86 0.11 -7.84 0.11 -7.81 0.11 -7.79 0.11 -7.76 0.11 -7.74 0.11 -7.71 0.11 -7.69 0.11 -7.66 0.11 -7.64 0.11 -7.61 0.11 -7.59 0.11 -7.56 0.11 -7.54 0.11 -7.51 0.11 -7.49 0.11 -7.46 0.11 -7.44 0.11 -7.41 0.11 -7.39 0.11 -7.36 0.11 -7.34 0.11 -7.31 0.11 -7.29 0.11 -7.26 0.11 -7.24 0.11 -7.21 0.11 -7.19 0.11 -7.16 0.11 -7.14 0.11 -7.11 0.11 -7.09 0.11 -7.07 0.11 -7.04 0.11 -7.02 0.11 -6.99 0.11 -6.97 0.11 -6.94 0.11 -6.92 0.11 -6.89 0.11 -6.87 0.11 -6.84 0.11 -6.82 0.11 -6.79 0.11 -6.77 0.11 -6.74 0.11 -6.72 0.11 -6.69 0.11 -6.67 0.11 -6.64 0.11 -6.62 0.11 -6.59 0.11 -6.57 0.11 -6.54 0.11 -6.52 0.11 -6.49 0.11 -6.47 0.11 -6.44 0.11 -6.42 0.11 -6.39 0.11 -6.37 0.11 -6.34 0.11 -6.32 0.11 -6.29 0.11 -6.27 0.11 -6.24 0.11 -6.22 0.11 -6.19 0.11 -6.17 0.11 -6.14 0.11 -6.12 0.11 -6.09 0.11 -6.07 0.11 -6.05 0.11 -6.02 0.11 -6 0.11 -5.97 0.11 -5.95 0.11 -5.92 0.11 -5.9 0.11 -5.87 0.11 -5.85 0.11 -5.82 0.11 -5.8 0.11 -5.77 0.11 -5.75 0.11 -5.72 0.11 -5.7 0.11 -5.67 0.11 -5.65 0.11 -5.62 0.11 -5.6 0.11 -5.57 0.11 -5.55 0.11 -5.52 0.11 -5.5 0.11 -5.47 0.11 -5.45 0.11 -5.42 0.11 -5.4 0.11 -5.37 0.11 -5.35 0.11 -5.32 0.11 -5.3 0.11 -5.27 0.11 -5.25 0.11 -5.22 0.11 -5.2 0.11 -5.17 0.11 -5.15 0.11 -5.12 0.11 -5.1 0.11 -5.07 0.11 -5.05 0.11 -5.03 0.11 -5 0.11 -4.98 0.11 -4.95 0.11 -4.93 0.11 -4.9 0.11 -4.88 0.11 -4.85 0.11 -4.83 0.11 -4.8 0.11 -4.78 0.11 -4.75 0.11 -4.73 0.11 -4.7 0.11 -4.68 0.11 -4.65 0.11 -4.63 0.11 -4.6 0.11 -4.58 0.11 -4.55 0.11 -4.53 0.11 -4.5 0.11 -4.48 0.11 -4.45 0.11 -4.43 0.11 -4.4 0.11 -4.38 0.11 -4.35 0.11 -4.33 0.11 -4.3 0.11 -4.28 0.11 -4.25 0.11 -4.23 0.11 -4.2 0.11 -4.18 0.11 -4.15 0.11 -4.13 0.11 -4.1 0.11 -4.08 0.11 -4.05 0.11 -4.03 0.11 -4.01 0.11 -3.98 0.11 -3.96 0.11 -3.93 0.11 -3.91 0.11 -3.88 0.11 -3.86 0.11 -3.83 0.11 -3.81 0.11 -3.78 0.11 -3.76 0.11 -3.73 0.11 -3.71 0.11 -3.68 0.11 -3.66 0.11 -3.63 0.11 -3.61 0.11 -3.58 0.11 -3.56 0.11 -3.53 0.11 -3.51 0.11 -3.48 0.11 -3.46 0.11 -3.43 0.11 -3.41 0.11 -3.38 0.11 -3.36 0.11 -3.33 0.11 -3.31 0.11 -3.28 0.11 -3.26 0.11 -3.23 0.11 -3.21 0.11 -3.18 0.11 -3.16 0.11 -3.13 0.11 -3.11 0.11 -3.08 0.11 -3.06 0.11 -3.03 0.11 -3.01 0.11 -2.99 0.11 -2.96 0.11 -2.94 0.11 -2.91 0.11 -2.89 0.11 -2.86 0.11 -2.84 0.11 -2.81 0.11 -2.79 0.11 -2.76 0.11 -2.74 0.11 -2.71 0.11 -2.69 0.11 -2.66 0.11 -2.64 0.11 -2.61 0.11 -2.59 0.11 -2.56 0.11 -2.54 0.11 -2.51 0.11 -2.49 0.11 -2.46 0.11 -2.44 0.11 -2.41 0.11 -2.39 0.11 -2.36 0.11 -2.34 0.11 -2.31 0.11 -2.29 0.11 -2.26 0.11 -2.24 0.11 -2.21 0.11 -2.19 0.11 -2.16 0.11 -2.14 0.11 -2.11 0.11 -2.09 0.11 -2.06 0.11 -2.04 0.11 -2.02 0.11 -1.99 0.11 -1.97 0.11 -1.94 0.11 -1.92 0.11 -1.89 0.11 -1.87 0.11 -1.84 0.11 -1.82 0.11 -1.79 0.11 -1.77 0.11 -1.74 0.11 -1.72 0.11 -1.69 0.11 -1.67 0.11 -1.64 0.11 -1.62 0.11 -1.59 0.11 -1.57 0.11 -1.54 0.11 -1.52 0.11 -1.49 0.11 -1.47 0.11 -1.44 0.11 -1.42 0.11 -1.39 0.11 -1.37 0.11 -1.34 0.11 -1.32 0.11 -1.29 0.11 -1.27 0.11 -1.24 0.11 -1.22 0.11 -1.19 0.11 -1.17 0.11 -1.14 0.11 -1.12 0.11 -1.09 0.11 -1.07 0.11 -1.04 0.11 -1.02 0.11 -1 0.11 -0.97 0.11 -0.95 0.11 -0.92 0.11 -0.9 0.11 -0.87 0.11 -0.85 0.11 -0.82 0.11 -0.8 0.11 -0.77 0.11 -0.75 0.11 -0.72 0.11 -0.7 0.11 -0.67 0.11 -0.65 0.11 -0.62 0.11 -0.6 0.11 -0.57 0.11 -0.55 0.11 -0.52 0.11 -0.5 0.11 -0.47 0.11 -0.45 0.11 -0.42 0.11 -0.4 0.11 -0.37 0.11 -0.35 0.11 -0.32 0.11 -0.3 0.11 -0.27 0.11 -0.25 0.11 -0.22 0.11 -0.2 0.11 -0.17 0.11 -0.15 0.11 -0.12 0.11 -0.1 0.11 -0.07 0.11 -0.05 0.11 -0.02 0.11 0 0.11 0.02 0.11 0.05 0.11 0.07 0.11 0.1 0.11 0.12 0.11 0.15 0.11 0.17 0.11 0.2 0.11 0.22 0.11 0.25 0.11 0.27 0.11 0.3 0.11 0.32 0.11 0.35 0.11 0.37 0.11 0.4 0.11 0.42 0.11 0.45 0.11 0.47 0.11 0.5 0.11 0.52 0.11 0.55 0.11 0.57 0.11 0.6 0.11 0.62 0.11 0.65 0.11 0.67 0.11 0.7 0.11 0.72 0.11 0.75 0.11 0.77 0.11 0.8 0.11 0.82 0.11 0.85 0.11 0.87 0.11 0.9 0.11 0.92 0.11 0.95 0.11 0.97 0.11 1 0.11 1.02 0.11 1.04 0.11 1.07 0.11 1.09 0.11 1.12 0.11 1.14 0.11 1.17 0.11 1.19 0.11 1.22 0.11 1.24 0.11 1.27 0.11 1.29 0.11 1.32 0.11 1.34 0.11 1.37 0.11 1.39 0.11 1.42 0.11 1.44 0.11 1.47 0.11 1.49 0.11 1.52 0.11 1.54 0.11 1.57 0.11 1.59 0.11 1.62 0.11 1.64 0.11 1.67 0.11 1.69 0.11 1.72 0.11 1.74 0.11 1.77 0.11 1.79 0.11 1.82 0.11 1.84 0.11 1.87 0.11 1.89 0.11 1.92 0.11 1.94 0.11 1.97 0.11 1.99 0.11 2.02 0.11 2.04 0.11 2.06 0.11 2.09 0.11 2.11 0.11 2.14 0.11 2.16 0.11 2.19 0.11 2.21 0.11 2.24 0.11 2.26 0.11 2.29 0.11 2.31 0.11 2.34 0.11 2.36 0.11 2.39 0.11 2.41 0.11 2.44 0.11 2.46 0.11 2.49 0.11 2.51 0.11 2.54 0.11 2.56 0.11 2.59 0.11 2.61 0.11 2.64 0.11 2.66 0.11 2.69 0.11 2.71 0.11 2.74 0.11 2.76 0.11 2.79 0.11 2.81 0.11 2.84 0.11 2.86 0.11 2.89 0.11 2.91 0.11 2.94 0.11 2.96 0.11 2.99 0.11 3.01 0.11 3.03 0.11 3.06 0.11 3.08 0.11 3.11 0.11 3.13 0.11 3.16 0.11 3.18 0.11 3.21 0.11 3.23 0.11 3.26 0.11 3.28 0.11 3.31 0.11 3.33 0.11 3.36 0.11 3.38 0.11 3.41 0.11 3.43 0.11 3.46 0.11 3.48 0.11 3.51 0.11 3.53 0.11 3.56 0.11 3.58 0.11 3.61 0.11 3.63 0.11 3.66 0.11 3.68 0.11 3.71 0.11 3.73 0.11 3.76 0.11 3.78 0.11 3.81 0.11 3.83 0.11 3.86 0.11 3.88 0.11 3.91 0.11 3.93 0.11 3.96 0.11 3.98 0.11 4.01 0.11 4.03 0.11 4.05 0.11 4.08 0.11 4.1 0.11 4.13 0.11 4.15 0.11 4.18 0.11 4.2 0.11 4.23 0.11 4.25 0.11 4.28 0.11 4.3 0.11 4.33 0.11 4.35 0.11 4.38 0.11 4.4 0.11 4.43 0.11 4.45 0.11 4.48 0.11 4.5 0.11 4.53 0.11 4.55 0.11 4.58 0.11 4.6 0.11 4.63 0.11 4.65 0.11 4.68 0.11 4.7 0.11 4.73 0.11 4.75 0.11 4.78 0.11 4.8 0.11 4.83 0.11 4.85 0.11 4.88 0.11 4.9 0.11 4.93 0.11 4.95 0.11 4.98 0.11 5 0.11 5.03 0.11 5.05 0.11 5.07 0.11 5.1 0.11 5.12 0.11 5.15 0.11 5.17 0.11 5.2 0.11 5.22 0.11 5.25 0.11 5.27 0.11 5.3 0.11 5.32 0.11 5.35 0.11 5.37 0.11 5.4 0.11 5.42 0.11 5.45 0.11 5.47 0.11 5.5 0.11 5.52 0.11 5.55 0.11 5.57 0.11 5.6 0.11 5.62 0.11 5.65 0.11 5.67 0.11 5.7 0.11 5.72 0.11 5.75 0.11 5.77 0.11 5.8 0.11 5.82 0.11 5.85 0.11 5.87 0.11 5.9 0.11 5.92 0.11 5.95 0.11 5.97 0.11 6 0.11 6.02 0.11 6.05 0.11 6.07 0.11 6.09 0.11 6.12 0.11 6.14 0.11 6.17 0.11 6.19 0.11 6.22 0.11 6.24 0.11 6.27 0.11 6.29 0.11 6.32 0.11 6.34 0.11 6.37 0.11 6.39 0.11 6.42 0.11 6.44 0.11 6.47 0.11 6.49 0.11 6.52 0.11 6.54 0.11 6.57 0.11 6.59 0.11 6.62 0.11 6.64 0.11 6.67 0.11 6.69 0.11 6.72 0.11 6.74 0.11 6.77 0.11 6.79 0.11 6.82 0.11 6.84 0.11 6.87 0.11 6.89 0.11 6.92 0.11 6.94 0.11 6.97 0.11 6.99 0.11 7.02 0.11 7.04 0.11 7.07 0.11 7.09 0.11 7.11 0.11 7.14 0.11 7.16 0.11 7.19 0.11 7.21 0.11 7.24 0.11 7.26 0.11 7.29 0.11 7.31 0.11 7.34 0.11 7.36 0.11 7.39 0.11 7.41 0.11 7.44 0.11 7.46 0.11 7.49 0.11 7.51 0.11 7.54 0.11 7.56 0.11 7.59 0.11 7.61 0.11 7.64 0.11 7.66 0.11 7.69 0.11 7.71 0.11 7.74 0.11 7.76 0.11 7.79 0.11 7.81 0.11 7.84 0.11 7.86 0.11 7.89 0.11 7.91 0.11 7.94 0.11 7.96 0.11 7.99 0.11 8.01 0.11 8.04 0.11 8.06 0.11 8.08 0.11 8.11 0.11 8.13 0.11 8.16 0.11 8.18 0.11 8.21 0.11 8.23 0.11 8.26 0.11 8.28 0.11 8.31 0.11 8.33 0.11 8.36 0.11 8.38 0.11 8.41 0.11 8.43 0.11 8.46 0.11 8.48 0.11 8.51 0.11 8.53 0.11 8.56 0.11 8.58 0.11 8.61 0.11 8.63 0.11 8.66 0.11 8.68 0.11 8.71 0.11 8.73 0.11 8.76 0.11 8.78 0.11 8.81 0.11 8.83 0.11 8.86 0.11 8.88 0.11 8.91 0.11 8.93 0.11 8.96 0.11 8.98 0.11 9.01 0.11 9.03 0.11 9.06 0.11 9.08 0.11 9.1 0.11 9.13 0.11 9.15 0.11 9.18 0.11 9.2 0.11 9.23 0.11 9.25 0.11 9.28 0.11 9.3 0.11 9.33 0.11 9.35 0.11 9.38 0.11 9.4 0.11 9.43 0.11 9.45 0.11 9.48 0.11 9.5 0.11 9.53 0.11 9.55 0.11 9.58 0.11 9.6 0.11 9.63 0.11 9.65 0.11 9.68 0.11 9.7 0.11 9.73 0.11 9.75 0.11 9.78 0.11 9.8 0.11 9.83 0.11 9.85 0.11 9.88 0.11 9.9 0.11 9.93 0.11 9.95 0.11 9.98 0.11 10 0.11 10.03 0.11 10.05 0.11 10.08 0.11 10.1 0.11 10.12 0.11 10.15 0.11 10.17 0.11 10.2 0.11 10.22 0.11 10.25 0.11 10.27 0.11 10.3 0.11 10.32 0.11 10.35 0.11 10.37 0.11 10.4 0.11 10.42 0.11 10.45 0.11 10.47 0.11 10.5 0.11 10.52 0.11 10.55 0.11 10.57 0.11 10.6 0.11 10.62 0.11 10.65 0.11 10.67 0.11 10.7 0.11 10.72 0.11 10.75 0.11 10.77 0.11 10.8 0.11 10.82 0.11 10.85 0.11 10.87 0.11 10.9 0.11 10.92 0.11 10.95 0.11 10.97 0.11 11 0.11 11.02 0.11 11.05 0.11 11.07 0.11 11.1 0.11 11.12 0.11 11.14 0.11 11.17 0.11 11.19 0.11 11.22 0.11 11.24 0.11 11.27 0.11 11.29 0.11 11.32 0.11 11.34 0.11 11.37 0.11 11.39 0.11 11.42 0.11 11.44 0.11 11.47 0.11 11.49 0.11 11.52 0.11 11.54 0.11 11.57 0.11 11.59 0.11 11.62 0.11 11.64 0.11 11.67 0.11 11.69 0.11 11.72 0.11 11.74 0.11 11.77 0.11 11.79 0.11 11.82 0.11 11.84 0.11 11.87 0.11 11.89 0.11 11.92 0.11 11.94 0.11 11.97 0.11 11.99 0.11 12.02 0.11 12.04 0.11 12.07 0.11 12.09 0.11 12.12 0.11 12.14 0.11 12.16 0.11 12.19 0.11 12.21 0.11 12.24 0.11 12.26 0.11 12.29 0.11 12.31 0.11 12.34 0.11 12.36 0.11 12.39 0.11 12.41 0.11 12.44 0.11 12.46 0.11 12.49 0.11 12.51 0.11 12.54 0.11 12.56 0.11 12.59 0.11 12.61 0.11 12.64 0.11 12.66 0.11 12.69 0.11 12.71 0.11 12.74 0.11 12.76 0.11 12.79 0.11 12.81 0.11 12.84 0.11 12.86 0.11 12.89 0.11 12.91 0.11 12.94 0.11 12.96 0.11 12.99 0.11 13.01 0.11 13.04 0.11 13.06 0.11 13.09 0.11 13.11 0.11 13.13 0.11 13.16 0.11 13.18 0.11 13.21 0.11 13.23 0.11 13.26 0.11 13.28 0.11 13.31 0.11 13.33 0.11 13.36 0.11 13.38 0.11 13.41 0.11 13.43 0.11 13.46 0.11 13.48 0.11 13.51 0.11 13.53 0.11 13.56 0.11 13.58 0.11 13.61 0.11 13.63 0.11 13.66 0.11 13.68 0.11 13.71 0.11 13.73 0.11 13.76 0.11 13.78 0.11 13.81 0.11 13.83 0.11 13.86 0.11 13.88 0.11 13.91 0.11 13.93 0.11 13.96 0.11 13.98 0.11 14.01 0.11 14.03 0.11 14.06 0.11 14.08 0.11 14.11 0.11 14.13 0.11 14.15 0.11 14.18 0.11 14.2 0.11 14.23 0.11 14.25 0.11 14.28 0.11 14.3 0.11 14.33 0.11 14.35 0.11 14.38 0.11 14.4 0.11 14.43 0.11 14.45 0.11 14.48 0.11 14.5 0.11 14.53 0.11 14.55 0.11 14.58 0.11 14.6 0.11 14.63 0.11 14.65 0.11 14.68 0.11 14.7 0.11 14.73 0.11 14.75 0.11 14.78 0.11 14.8 0.11 14.83 0.11 14.85 0.11 14.88 0.11 14.9 0.11 14.93 0.11 14.95 0.11 14.98 0.11 15 0.11 15.03 0.11 15.05 0.11 15.08 0.11 15.1 0.11 15.13 0.11 15.15 0.11 15.17 0.11 15.2 0.11 15.22 0.11 15.25 0.11 15.27 0.11 15.3 0.11 15.32 0.11 15.35 0.11 15.37 0.11 15.4 0.11 15.42 0.11 15.45 0.11 15.47 0.11 15.5 0.11 15.52 0.11 15.55 0.11 15.57 0.11 15.6 0.11 15.62 0.11 15.65 0.11 15.67 0.11 15.7 0.11 15.72 0.11 15.75 0.11 15.77 0.11 15.8 0.11 15.82 0.11 15.85 0.11 15.87 0.11 15.9 0.11 15.92 0.11 15.95 0.11 15.97 0.11 16 0.11 16.02 0.11 16.05 0.11 16.07 0.11 16.1 0.11 16.12 0.11 16.15 0.11 16.17 0.11 16.19 0.11 16.22 0.11 16.24 0.11 16.27 0.11 16.29 0.11 16.32 0.11 16.34 0.11 16.37 0.11 16.39 0.11 16.42 0.11 16.44 0.11 16.47 0.11 16.49 0.11 16.52 0.11 16.54 0.11 16.57 0.11 16.59 0.11 16.62 0.11 16.64 0.11 16.67 0.11 16.69 0.11 16.72 0.11 16.74 0.11 16.77 0.11 16.79 0.11 16.82 0.11 16.84 0.11 16.87 0.11 16.89 0.11 16.92 0.11 16.94 0.11 16.97 0.11 16.99 0.11 17.02 0.11 17.04 0.11 17.07 0.11 17.09 0.11 17.12 0.11 17.14 0.11 17.16 0.11 17.19 0.11 17.21 0.11 17.24 0.11 17.26 0.11 17.29 0.11 17.31 0.11 17.34 0.11 17.36 0.11 17.39 0.11 17.41 0.11 17.44 0.11 17.46 0.11 17.49 0.11 17.51 0.11 17.54 0.11 17.56 0.11 17.59 0.11 17.61 0.11 17.64 0.11 17.66 0.11 17.69 0.11 17.71 0.11 17.74 0.11 17.76 0.11 17.79 0.11 17.81 0.11 17.84 0.11 17.86 0.11 17.89 0.11 17.91 0.11 17.94 0.11 17.96 0.11 17.99 0.11 18.01 0.11 18.04 0.11 18.06 0.11 18.09 0.11 18.11 0.11 18.14 0.11 18.16 0.11 18.18 0.11 18.21 0.11 18.23 0.11 18.26 0.11 18.28 0.11 18.31 0.11 18.33 0.11 18.36 0.11 18.38 0.11 18.41 0.11 18.43 0.11 18.46 0.11 18.48 0.11 18.51 0.11 18.53 0.11 18.56 0.11 18.58 0.11 18.61 0.11 18.63 0.11 18.66 0.11 18.68 0.11 18.71 0.11 18.73 0.11 18.76 0.11 18.78 0.11 18.81 0.11 18.83 0.11 18.86 0.11 18.88 0.11 18.91 0.11 18.93 0.11 18.96 0.11 18.98 0.11 19.01 0.11 19.03 0.11 19.06 0.11 19.08 0.11 19.11 0.11 19.13 0.11 19.16 0.11 19.18 0.11 19.2 0.11 19.23 0.11 19.25 0.11 19.28 0.11 19.3 0.11 19.33 0.11 19.35 0.11 19.38 0.11 19.4 0.11 19.43 0.11 19.45 0.11 19.48 0.11 19.5 0.11 19.53 0.11 19.55 0.11 19.58 0.11 19.6 0.11 19.63 0.11 19.65 0.11 19.68 0.11 19.7 0.11 19.73 0.11 19.75 0.11 19.78 0.11 19.8 0.11 19.83 0.11 19.85 0.11 19.88 0.11 19.9 0.11 19.93 0.11 19.95 0.11 19.98 0.11 20 0.11 20.03 0.11 20.05 0.11 20.08 0.11 20.1 0.11 20.13 0.11 20.15 0.11 20.18 0.11 20.2 0.11 20.22 0.11 20.25 0.11 20.27 0.11 20.3 0.11 20.32 0.11 20.35 0.11 20.37 0.11 20.4 0.11 20.42 0.11 20.45 0.11 20.47 0.11 20.5 0.11 20.52 0.11 20.55 0.11 20.57 0.11 20.6 0.11 20.62 0.11 20.65 0.11 20.67 0.11 20.7 0.11 20.72 0.11 20.75 0.11 20.77 0.11 20.8 0.11 20.82 0.11 20.85 0.11 20.87 0.11 20.9 0.11 20.92 0.11 20.95 0.11 20.97 0.11 21 0.11 21.02 0.11 21.05 0.11 21.07 0.11 21.1 0.11 21.12 0.11 21.15 0.11 21.17 0.11 21.2 0.11 21.22 0.11 21.24 0.11 21.27 0.11 21.29 0.11 21.32 0.11 21.34 0.11 21.37 0.11 21.39 0.11 21.42 0.11 21.44 0.11 21.47 0.11 21.49 0.11 21.52 0.11 21.54 0.11 21.57 0.11 21.59 0.11 21.62 0.11 21.64 0.11 21.67 0.11 21.69 0.11 21.72 0.11 21.74 0.11 21.77 0.11 21.79 0.11 21.82 0.11 21.84 0.11 21.87 0.11 21.89 0.11 21.92 0.11 21.94 0.11 21.97 0.11 21.99 0.11 22.02 0.11 22.04 0.11 22.07 0.11 22.09 0.11 22.12 0.11 22.14 0.11 22.17 0.11 22.19 0.11 22.21 0.11 22.24 0.11 22.26 0.11 22.29 0.11 22.31 0.11 22.34 0.11 22.36 0.11 22.39 0.11 22.41 0.11 22.44 0.11 22.46 0.11 22.49 0.11 22.51 0.11 22.54 0.11 22.56 0.11 22.59 0.11 22.61 0.11 22.64 0.11 22.66 0.11 22.69 0.11 22.71 0.11 22.74 0.11 22.76 0.11 22.79 0.11 22.81 0.11 22.84 0.11 22.86 0.11 22.89 0.11 22.91 0.11 22.94 0.11 22.96 0.11 22.99 0.11 23.01 0.11 23.04 0.11 23.06 0.11 23.09 0.11 23.11 0.11 23.14 0.11 23.16 0.11 23.19 0.11 23.21 0.11 23.23 0.11 23.26 0.11 23.28 0.11 23.31 0.11 23.33 0.11 23.36 0.11 23.38 0.11 23.41 0.11 23.43 0.11 23.46 0.11 23.48 0.11 23.51 0.11 23.53 0.11 23.56 0.11 23.58 0.11 23.61 0.11 23.63 0.11 23.66 0.11 23.68 0.11 23.71 0.11 23.73 0.11 23.76 0.11 23.78 0.11 23.81 0.11 23.83 0.11 23.86 0.11 23.88 0.11 23.91 0.11 23.93 0.11 23.96 0.11 23.98 0.11 24.01 0.11 24.03 0.11 24.06 0.11 24.08 0.11 24.11 0.11 24.13 0.11 24.16 0.11 24.18 0.11 24.21 0.11 24.23 0.11 24.25 0.11 24.28 0.11 24.3 0.11 24.33 0.11 24.35 0.11 24.38 0.11 24.4 0.11 24.43 0.11 24.45 0.11 24.48 0.11 24.5 0.11 24.53 0.11 24.55 0.11 24.58 0.11 24.6 0.11 24.63 0.11 24.65 0.11 24.68 0.11 24.7 0.11 24.73 0.11 24.75 0.11 24.78 0.11 24.8 0.11 24.83 0.11 24.85 0.11 24.88 0.11 24.9 0.11 24.93 0.11 24.95 0.11 24.98 0.11 25 0.11 25.03 0.11 25.05 0.11 25.08 0.11 25.1 0.11 25.13 0.11 25.15 0.11 25.18 0.11 25.2 0.11 25.23 0.11 25.25 0.11 25.27 0.11 25.3 0.11 25.32 0.11 25.35 0.11 25.37 0.11 25.4 0.11 25.42 0.11 25.45 0.11 25.47 0.11 25.5 0.11 25.52 0.11 25.55 0.11 25.57 0.11 25.6 0.11 25.62 0.11 25.65 0.11 25.67 0.11 25.7 0.11 25.72 0.11 25.75 0.11 25.77 0.11 25.8 0.11 25.82 0.11 25.85 0.11 25.87 0.11 25.9 0.11 25.92 0.11 25.95 0.11 25.97 0.11 26 0.11 26.02 0.11 26.05 0.11 26.07 0.11 26.1 0.11 26.12 0.11 26.15 0.11 26.17 0.11 26.2 0.11 26.22 0.11 26.25 0.11 26.27 0.11 26.29 0.11 26.32 0.11 26.34 0.11 26.37 0.11 26.39 0.11 26.42 0.11 26.44 0.11 26.47 0.11 26.49 0.11 26.52 0.11 26.54 0.11 26.57 0.11 26.59 0.11 26.62 0.11 26.64 0.11 26.67 0.11 26.69 0.11 26.72 0.11 26.74 0.11 26.77 0.11 26.79 0.11 26.82 0.11 26.84 0.11 26.87 0.11 26.89 0.11 26.92 0.11 26.94 0.11 26.97 0.11 26.99 0.11 27.02 0.11 27.04 0.11 27.07 0.11 27.09 0.11 27.12 0.11 27.14 0.11 27.17 0.11 27.19 0.11 27.22 0.11 27.24 0.11 27.26 0.11 27.29 0.11 27.31 0.11 27.34 0.11 27.36 0.11 27.39 0.11 27.41 0.11 27.44 0.11 27.46 0.11 27.49 0.11 27.51 0.11 27.54 0.11 27.56 0.11 27.59 0.11 27.61 0.11 27.64 0.11 27.66 0.11 27.69 0.11 27.71 0.11 27.74 0.11 27.76 0.11 27.79 0.11 27.81 0.11 27.84 0.11 27.86 0.11 27.89 0.11 27.91 0.11 27.94 0.11 27.96 0.11 27.99 0.11 28.01 0.11 28.04 0.11 28.06 0.11 28.09 0.11 28.11 0.11 28.14 0.11 28.16 0.11 28.19 0.11 28.21 0.11 28.24 0.11 28.26 0.11 28.28 0.11 28.31 0.11 28.33 0.11 28.36 0.11 28.38 0.11 28.41 0.11 28.43 0.11 28.46 0.11 28.48 0.11 28.51 0.11 28.53 0.11 28.56 0.11 28.58 0.11 28.61 0.11 28.63 0.11 28.66 0.11 28.68 0.11 28.71 0.11 28.73 0.11 28.76 0.11 28.78 0.11 28.81 0.11 28.83 0.11 28.86 0.11 28.88 0.11 28.91 0.11 28.93 0.11 28.96 0.11 28.98 0.11 29.01 0.11 29.03 0.11 29.06 0.11 29.08 0.11 29.11 0.11 29.13 0.11 29.16 0.11 29.18 0.11 29.21 0.11 29.23 0.11 29.26 0.11 29.28 0.11 29.3 0.11 29.33 0.11 29.35 0.11 29.38 0.11 29.4 0.11 29.43 0.11 29.45 0.11 29.48 0.11 29.5 0.11 29.53 0.11 29.55 0.11 29.58 0.11 29.6 0.11 29.63 0.11 29.65 0.11 29.68 0.11 29.7 0.11 29.73 0.11 29.75 0.11 29.78 0.11 29.8 0.11 29.83 0.11 29.85 0.11 29.88 0.11 29.9 0.11 29.93 0.11 29.95 0.11 29.98 0.11 30 0.11 30.03 0.11 30.05 0.11 30.08 0.11 30.1 0.11 30.13 0.11 30.15 0.11 30.18 0.11 30.2 0.11 30.23 0.11 30.25 0.11 30.28 0.11 30.3 0.11 30.32 0.11 30.35 0.11 30.37 0.11 30.4 0.11 30.42 0.11 30.45 0.11 30.47 0.11 30.5 0.11 30.52 0.11 30.55 0.11 30.57 0.11 30.6 0.11 30.62 0.11 30.65 0.11 30.67 0.11 30.7 0.11 30.72 0.11 30.75 0.11 30.77 0.11 30.8 0.11 30.82 0.11 30.85 0.11 30.87 0.11 30.9 0.11 30.92 0.11 30.95 0.11 30.97 0.11 31 0.11 31.02 0.11 31.05 0.11 31.07 0.11 31.1 0.11 31.12 0.11 31.15 0.11 31.17 0.11 31.2 0.11 31.22 0.11 31.25 0.11 31.27 0.11 31.3 0.11 31.32 0.11 31.34 0.11 31.37 0.11 31.39 0.11 31.42 0.11 31.44 0.11 31.47 0.11 31.49 0.11 31.52 0.11 31.54 0.11 31.57 0.11 31.59 0.11 31.62 0.11 31.64 0.11 31.67 0.11 31.69 0.11 31.72 0.11 31.74 0.11 31.77 0.11 31.79 0.11 31.82 0.11 31.84 0.11 31.87 0.11 31.89 0.11 31.92 0.11 31.94 0.11 31.97 0.11 31.99 0.11 32.02 0.11 32.04 0.11 32.07 0.11 32.09 0.11 32.12 0.11 32.14 0.11 32.17 0.11 32.19 0.11 32.22 0.11 32.24 0.11 32.27 0.11 32.29 0.11 32.31 0.11 32.34 0.11 32.36 0.11 32.39 0.11 32.41 0.11 32.44 0.11 32.46 0.11 32.49 0.11 32.51 0.11 32.54 0.11 32.56 0.11 32.59 0.11 32.61 0.11 32.64 0.11 32.66 0.11 32.69 0.11 32.71 0.11 32.74 0.11 32.76 0.11 32.79 0.11 32.81 0.11 32.84 0.11 32.86 0.11 32.89 0.11 32.91 0.11 32.94 0.11 32.96 0.11 32.99 0.11 33.01 0.11 33.04 0.11 33.06 0.11 33.09 0.11 33.11 0.11 33.14 0.11 33.16 0.11 33.19 0.11 33.21 0.11 33.24 0.11 33.26 0.11 33.29 0.11 33.31 0.11 33.33 0.11 33.36 0.11 33.38 0.11 33.41 0.11 33.43 0.11 33.46 0.11 33.48 0.11 33.51 0.11 33.53 0.11 33.56 0.11 33.58 0.11 33.61 0.11 33.63 0.11 33.66 0.11 33.68 0.11 33.71 0.11 33.73 0.11 33.76 0.11 33.78 0.11 33.81 0.11 33.83 0.11 33.86 0.11 33.88 0.11 33.91 0.11 33.93 0.11 33.96 0.11 33.98 0.11 34.01 0.11 34.03 0.11 34.06 0.11 34.08 0.11 34.11 0.11 34.13 0.11 34.16 0.11 34.18 0.11 34.21 0.11 34.23 0.11 34.26 0.11 34.28 0.11 34.31 0.11 34.33 0.11 34.35 0.11 34.38 0.11 34.4 0.11 34.43 0.11 34.45 0.11 34.48 0.11 34.5 0.11 34.53 0.11 34.55 0.11 34.58 0.11 34.6 0.11 34.63 0.11 34.65 0.11 34.68 0.11 34.7 0.11 34.73 0.11 34.75 0.11 34.78 0.11 34.8 0.11 34.83 0.11 34.85 0.11 34.88 0.11 34.9 0.11 34.93 0.11 34.95 0.11 34.98 0.11 35 0.11 35.03 0.11 35.05 0.11 35.08 0.11 35.1 0.11 35.13 0.11 35.15 0.11 35.18 0.11 35.2 0.11 35.23 0.11 35.25 0.11 35.28 0.11 35.3 0.11 35.33 0.11 35.35 0.11 35.37 0.11 35.4 0.11 35.42 0.11 35.45 0.11 35.47 0.11 35.5 0.11 35.52 0.11 35.55 0.11 35.57 0.11 35.6 0.11 35.62 0.11 35.65 0.11 35.67 0.11 35.7 0.11 35.72 0.11 35.75 0.11 35.77 0.11 35.8 0.11 35.82 0.11 35.85 0.11 35.87 0.11 35.9 0.11 35.92 0.11 35.95 0.11 35.97 0.11 36 0.11 36.02 0.11 36.05 0.11 36.07 0.11 36.1 0.11 36.12 0.11 36.15 0.11 36.17 0.11 36.2 0.11 36.22 0.11 36.25 0.11 36.27 0.11 36.3 0.11 36.32 0.11 36.35 0.11 36.37 0.11 36.39 0.11 36.42 0.11 36.44 0.11 36.47 0.11 36.49 0.11 36.52 0.11 36.54 0.11 36.57 0.11 36.59 0.11 36.62 0.11 36.64 0.11 36.67 0.11 36.69 0.11 36.72 0.11 36.74 0.11 36.77 0.11 36.79 0.11 36.82 0.11 36.84 0.11 36.87 0.11 36.89 0.11 36.92 0.11 36.94 0.11 36.97 0.11 36.99 0.11 37.02 0.11 37.04 0.11 37.07 0.11 37.09 0.11 37.12 0.11 37.14 0.11 37.17 0.11 37.19 0.11 37.22 0.11 37.24 0.11 37.27 0.11 37.29 0.11 37.32 0.11 37.34 0.11 37.36 0.11 37.39 0.11 37.41 0.11 37.44 0.11 37.46 0.11 37.49 0.11 37.51 0.11 37.54 0.11 37.56 0.11 37.59 0.11 37.61 0.11 37.64 0.11 37.66 0.11 37.69 0.11 37.71 0.11 37.74 0.11 37.76 0.11 37.79 0.11 37.81 0.11 37.84 0.11 37.86 0.11 37.89 0.11 37.91 0.11 37.94 0.11 37.96 0.11 37.99 0.11 38.01 0.11 38.04 0.11 38.06 0.11 38.09 0.11 38.11 0.11 38.14 0.11 38.16 0.11 38.19 0.11 38.21 0.11 38.24 0.11 38.26 0.11 38.29 0.11 38.31 0.11 38.34 0.11 38.36 0.11 38.38 0.11 38.41 0.11 38.43 0.11 38.46 0.11 38.48 0.11 38.51 0.11 38.53 0.11 38.56 0.11 38.58 0.11 38.61 0.11 38.63 0.11 38.66 0.11 38.68 0.11 38.71 0.11 38.73 0.11 38.76 0.11 38.78 0.11 38.81 0.11 38.83 0.11 38.86 0.11 38.88 0.11 38.91 0.11 38.93 0.11 38.96 0.11 38.98 0.11 39.01 0.11 39.03 0.11 39.06 0.11 39.08 0.11 39.11 0.11 39.13 0.11 39.16 0.11 39.18 0.11 39.21 0.11 39.23 0.11 39.26 0.11 39.28 0.11 39.31 0.11 39.33 0.11 39.36 0.11 39.38 0.11 39.4 0.11 39.43 0.11 39.45 0.11 39.48 0.11 39.5 0.11 39.53 0.11 39.55 0.11 39.58 0.11 39.6 0.11 39.63 0.11 39.65 0.11 39.68 0.11 39.7 0.11 39.73 0.11 39.75 0.11 39.78 0.11 39.8 0.11 39.83 0.11 39.85 0.11 39.88 0.11 39.9 0.11 39.93 0.11 39.95 0.11 39.98 0.11 40 0.11 40.03 0.11 40.05 0.11 40.08 0.11 40.1 0.11 40.13 0.11 40.15 0.11 40.18 0.11 40.2 0.11 40.23 0.11 40.25 0.11 40.28 0.11 40.3 0.11 40.33 0.11 40.35 0.11 40.38 0.11 40.4 0.11 40.42 0.11 40.45 0.11 40.47 0.11 40.5 0.11 40.52 0.11 40.55 0.11 40.57 0.11 40.6 0.11 40.62 0.11 40.65 0.11 40.67 0.11 40.7 0.11 40.72 0.11 40.75 0.11 40.77 0.11 40.8 0.11 40.82 0.11 40.85 0.11 40.87 0.11 40.9 0.11 40.92 0.11 40.95 0.11 40.97 0.11 41 0.11 41.02 0.11 41.05 0.11 41.07 0.11 41.1 0.11 41.12 0.11 41.15 0.11 41.17 0.11 41.2 0.11 41.22 0.11 41.25 0.11 41.27 0.11 41.3 0.11 41.32 0.11 41.35 0.11 41.37 0.11 41.39 0.11 41.42 0.11 41.44 0.11 41.47 0.11 41.49 0.11 41.52 0.11 41.54 0.11 41.57 0.11 41.59 0.11 41.62 0.11 41.64 0.11 41.67 0.11 41.69 0.11 41.72 0.11 41.74 0.11 41.77 0.11 41.79 0.11 41.82 0.11 41.84 0.11 41.87 0.11 41.89 0.11 41.92 0.11 41.94 0.11 41.97 0.11 41.99 0.11 42.02 0.11 42.04 0.11 42.07 0.11 42.09 0.11 42.12 0.11 42.14 0.11 42.17 0.11 42.19 0.11 42.22 0.11 42.24 0.11 42.27 0.11 42.29 0.11 42.32 0.11 42.34 0.11 42.37 0.11 42.39 0.11 42.41 0.11 42.44 0.11 42.46 0.11 42.49 0.11 42.51 0.11 42.54 0.11 42.56 0.11 42.59 0.11 42.61 0.11 42.64 0.11 42.66 0.11 42.69 0.11 42.71 0.11 42.74 0.11 42.76 0.11 42.79 0.11 42.81 0.11 42.84 0.11 42.86 0.11 42.89 0.11 42.91 0.11 42.94 0.11 42.96 0.11 42.99 0.11 43.01 0.11 43.04 0.11 43.06 0.11 43.09 0.11 43.11 0.11 43.14 0.11 43.16 0.11 43.19 0.11 43.21 0.11 43.24 0.11 43.26 0.11 43.29 0.11 43.31 0.11 43.34 0.11 43.36 0.11 43.39 0.11 43.41 0.11 43.43 0.11 43.46 0.11 43.48 0.11 43.51 0.11 43.53 0.11 43.56 0.11 43.58 0.11 43.61 0.11 43.63 0.11 43.66 0.11 43.68 0.11 43.71 0.11 43.73 0.11 43.76 0.11 43.78 0.11 43.81 0.11 43.83 0.11 43.86 0.11 43.88 0.11 43.91 0.11 43.93 0.11 43.96 0.11 43.98 0.11 44.01 0.11 44.03 0.11 44.06 0.11 44.08 0.11 44.11 0.11 44.13 0.11 44.16 0.11 44.18 0.11 44.21 0.11 44.23 0.11 44.26 0.11 44.28 0.11 44.31 0.11 44.33 0.11 44.36 0.11 44.38 0.11 44.41 0.11 44.43 0.11 44.45 0.11 44.48 0.11 44.5 0.11 44.53 0.11 44.55 0.11 44.58 0.11 44.6 0.11 44.63 0.11 44.65 0.11 44.68 0.11 44.7 0.11 44.73 0.11 44.75 0.11 44.78 0.11 44.8 0.11 44.83 0.11 44.85 0.11 44.88 0.11 44.9 0.11 44.93 0.11 44.95 0.11 44.98 0.11 45 0.11 45.03 0.11 45.05 0.11 45.08 0.11 45.1 0.11 45.13 0.11 45.15 0.11 45.18 0.11 45.2 0.11 45.23 0.11 45.25 0.11 45.28 0.11 45.3 0.11 45.33 0.11 45.35 0.11 45.38 0.11 45.4 0.11" class="primitive"/>
          </g>
          <g transform="translate(67.03,27.71)" id="img-8311e7c5-748" class="geometry color_I_H" stroke="#FF6DAE">
            <path fill="none" d="M-45.4,0.86 L -45.38 0.86 -45.35 0.86 -45.33 0.86 -45.3 0.86 -45.28 0.86 -45.25 0.86 -45.23 0.86 -45.2 0.86 -45.18 0.86 -45.15 0.87 -45.13 0.87 -45.1 0.87 -45.08 0.87 -45.05 0.87 -45.03 0.87 -45 0.87 -44.98 0.87 -44.95 0.87 -44.93 0.87 -44.9 0.87 -44.88 0.88 -44.85 0.88 -44.83 0.88 -44.8 0.88 -44.78 0.88 -44.75 0.88 -44.73 0.88 -44.7 0.88 -44.68 0.88 -44.65 0.88 -44.63 0.88 -44.6 0.88 -44.58 0.88 -44.55 0.88 -44.53 0.88 -44.5 0.88 -44.48 0.88 -44.45 0.88 -44.43 0.88 -44.41 0.88 -44.38 0.88 -44.36 0.88 -44.33 0.88 -44.31 0.88 -44.28 0.88 -44.26 0.88 -44.23 0.88 -44.21 0.88 -44.18 0.88 -44.16 0.88 -44.13 0.88 -44.11 0.88 -44.08 0.88 -44.06 0.87 -44.03 0.87 -44.01 0.87 -43.98 0.87 -43.96 0.87 -43.93 0.87 -43.91 0.87 -43.88 0.86 -43.86 0.86 -43.83 0.86 -43.81 0.86 -43.78 0.86 -43.76 0.86 -43.73 0.85 -43.71 0.85 -43.68 0.85 -43.66 0.85 -43.63 0.84 -43.61 0.84 -43.58 0.84 -43.56 0.84 -43.53 0.83 -43.51 0.83 -43.48 0.83 -43.46 0.82 -43.43 0.82 -43.41 0.82 -43.39 0.81 -43.36 0.81 -43.34 0.81 -43.31 0.8 -43.29 0.8 -43.26 0.8 -43.24 0.79 -43.21 0.79 -43.19 0.78 -43.16 0.78 -43.14 0.77 -43.11 0.77 -43.09 0.77 -43.06 0.76 -43.04 0.76 -43.01 0.75 -42.99 0.75 -42.96 0.74 -42.94 0.74 -42.91 0.73 -42.89 0.72 -42.86 0.72 -42.84 0.71 -42.81 0.71 -42.79 0.7 -42.76 0.69 -42.74 0.69 -42.71 0.68 -42.69 0.68 -42.66 0.67 -42.64 0.66 -42.61 0.65 -42.59 0.65 -42.56 0.64 -42.54 0.63 -42.51 0.63 -42.49 0.62 -42.46 0.61 -42.44 0.6 -42.41 0.59 -42.39 0.59 -42.37 0.58 -42.34 0.57 -42.32 0.56 -42.29 0.55 -42.27 0.54 -42.24 0.53 -42.22 0.53 -42.19 0.52 -42.17 0.51 -42.14 0.5 -42.12 0.49 -42.09 0.48 -42.07 0.47 -42.04 0.46 -42.02 0.45 -41.99 0.43 -41.97 0.42 -41.94 0.41 -41.92 0.4 -41.89 0.39 -41.87 0.38 -41.84 0.37 -41.82 0.35 -41.79 0.34 -41.77 0.33 -41.74 0.32 -41.72 0.3 -41.69 0.29 -41.67 0.28 -41.64 0.26 -41.62 0.25 -41.59 0.24 -41.57 0.22 -41.54 0.21 -41.52 0.19 -41.49 0.18 -41.47 0.16 -41.44 0.15 -41.42 0.13 -41.39 0.11 -41.37 0.1 -41.35 0.08 -41.32 0.07 -41.3 0.05 -41.27 0.03 -41.25 0.01 -41.22 -0 -41.2 -0.02 -41.17 -0.04 -41.15 -0.06 -41.12 -0.08 -41.1 -0.1 -41.07 -0.12 -41.05 -0.14 -41.02 -0.16 -41 -0.18 -40.97 -0.2 -40.95 -0.22 -40.92 -0.24 -40.9 -0.26 -40.87 -0.28 -40.85 -0.31 -40.82 -0.33 -40.8 -0.35 -40.77 -0.38 -40.75 -0.4 -40.72 -0.42 -40.7 -0.45 -40.67 -0.47 -40.65 -0.5 -40.62 -0.52 -40.6 -0.55 -40.57 -0.58 -40.55 -0.6 -40.52 -0.63 -40.5 -0.66 -40.47 -0.69 -40.45 -0.71 -40.42 -0.74 -40.4 -0.77 -40.38 -0.8 -40.35 -0.83 -40.33 -0.86 -40.3 -0.89 -40.28 -0.92 -40.25 -0.95 -40.23 -0.99 -40.2 -1.02 -40.18 -1.05 -40.15 -1.08 -40.13 -1.12 -40.1 -1.15 -40.08 -1.18 -40.05 -1.22 -40.03 -1.25 -40 -1.29 -39.98 -1.32 -39.95 -1.36 -39.93 -1.4 -39.9 -1.43 -39.88 -1.47 -39.85 -1.51 -39.83 -1.55 -39.8 -1.59 -39.78 -1.63 -39.75 -1.66 -39.73 -1.7 -39.7 -1.74 -39.68 -1.79 -39.65 -1.83 -39.63 -1.87 -39.6 -1.91 -39.58 -1.95 -39.55 -1.99 -39.53 -2.04 -39.5 -2.08 -39.48 -2.12 -39.45 -2.17 -39.43 -2.21 -39.4 -2.26 -39.38 -2.3 -39.36 -2.35 -39.33 -2.39 -39.31 -2.44 -39.28 -2.48 -39.26 -2.53 -39.23 -2.58 -39.21 -2.62 -39.18 -2.67 -39.16 -2.72 -39.13 -2.77 -39.11 -2.81 -39.08 -2.86 -39.06 -2.91 -39.03 -2.96 -39.01 -3.01 -38.98 -3.06 -38.96 -3.11 -38.93 -3.16 -38.91 -3.2 -38.88 -3.25 -38.86 -3.3 -38.83 -3.35 -38.81 -3.4 -38.78 -3.45 -38.76 -3.5 -38.73 -3.55 -38.71 -3.6 -38.68 -3.65 -38.66 -3.7 -38.63 -3.75 -38.61 -3.8 -38.58 -3.85 -38.56 -3.9 -38.53 -3.95 -38.51 -4 -38.48 -4.05 -38.46 -4.1 -38.43 -4.14 -38.41 -4.19 -38.38 -4.24 -38.36 -4.29 -38.34 -4.34 -38.31 -4.38 -38.29 -4.43 -38.26 -4.48 -38.24 -4.52 -38.21 -4.57 -38.19 -4.61 -38.16 -4.66 -38.14 -4.7 -38.11 -4.75 -38.09 -4.79 -38.06 -4.83 -38.04 -4.88 -38.01 -4.92 -37.99 -4.96 -37.96 -5 -37.94 -5.04 -37.91 -5.08 -37.89 -5.12 -37.86 -5.16 -37.84 -5.2 -37.81 -5.23 -37.79 -5.27 -37.76 -5.31 -37.74 -5.34 -37.71 -5.38 -37.69 -5.41 -37.66 -5.44 -37.64 -5.47 -37.61 -5.51 -37.59 -5.54 -37.56 -5.56 -37.54 -5.59 -37.51 -5.62 -37.49 -5.65 -37.46 -5.67 -37.44 -5.7 -37.41 -5.72 -37.39 -5.75 -37.36 -5.77 -37.34 -5.79 -37.32 -5.81 -37.29 -5.83 -37.27 -5.85 -37.24 -5.87 -37.22 -5.88 -37.19 -5.9 -37.17 -5.91 -37.14 -5.93 -37.12 -5.94 -37.09 -5.95 -37.07 -5.97 -37.04 -5.98 -37.02 -5.98 -36.99 -5.99 -36.97 -6 -36.94 -6.01 -36.92 -6.01 -36.89 -6.02 -36.87 -6.02 -36.84 -6.02 -36.82 -6.03 -36.79 -6.03 -36.77 -6.03 -36.74 -6.03 -36.72 -6.02 -36.69 -6.02 -36.67 -6.02 -36.64 -6.01 -36.62 -6.01 -36.59 -6 -36.57 -6 -36.54 -5.99 -36.52 -5.98 -36.49 -5.97 -36.47 -5.96 -36.44 -5.95 -36.42 -5.94 -36.39 -5.93 -36.37 -5.92 -36.35 -5.9 -36.32 -5.89 -36.3 -5.88 -36.27 -5.86 -36.25 -5.84 -36.22 -5.83 -36.2 -5.81 -36.17 -5.79 -36.15 -5.77 -36.12 -5.76 -36.1 -5.74 -36.07 -5.72 -36.05 -5.7 -36.02 -5.68 -36 -5.65 -35.97 -5.63 -35.95 -5.61 -35.92 -5.59 -35.9 -5.56 -35.87 -5.54 -35.85 -5.52 -35.82 -5.49 -35.8 -5.47 -35.77 -5.44 -35.75 -5.42 -35.72 -5.39 -35.7 -5.37 -35.67 -5.34 -35.65 -5.31 -35.62 -5.29 -35.6 -5.26 -35.57 -5.23 -35.55 -5.21 -35.52 -5.18 -35.5 -5.15 -35.47 -5.12 -35.45 -5.09 -35.42 -5.07 -35.4 -5.04 -35.37 -5.01 -35.35 -4.98 -35.33 -4.95 -35.3 -4.92 -35.28 -4.89 -35.25 -4.86 -35.23 -4.83 -35.2 -4.8 -35.18 -4.78 -35.15 -4.75 -35.13 -4.72 -35.1 -4.69 -35.08 -4.66 -35.05 -4.63 -35.03 -4.6 -35 -4.57 -34.98 -4.54 -34.95 -4.51 -34.93 -4.48 -34.9 -4.45 -34.88 -4.42 -34.85 -4.39 -34.83 -4.36 -34.8 -4.33 -34.78 -4.3 -34.75 -4.27 -34.73 -4.24 -34.7 -4.21 -34.68 -4.18 -34.65 -4.15 -34.63 -4.12 -34.6 -4.09 -34.58 -4.06 -34.55 -4.04 -34.53 -4.01 -34.5 -3.98 -34.48 -3.95 -34.45 -3.92 -34.43 -3.89 -34.4 -3.86 -34.38 -3.83 -34.35 -3.81 -34.33 -3.78 -34.31 -3.75 -34.28 -3.72 -34.26 -3.69 -34.23 -3.67 -34.21 -3.64 -34.18 -3.61 -34.16 -3.58 -34.13 -3.56 -34.11 -3.53 -34.08 -3.5 -34.06 -3.47 -34.03 -3.45 -34.01 -3.42 -33.98 -3.39 -33.96 -3.37 -33.93 -3.34 -33.91 -3.31 -33.88 -3.29 -33.86 -3.26 -33.83 -3.24 -33.81 -3.21 -33.78 -3.19 -33.76 -3.16 -33.73 -3.14 -33.71 -3.11 -33.68 -3.09 -33.66 -3.06 -33.63 -3.04 -33.61 -3.01 -33.58 -2.99 -33.56 -2.96 -33.53 -2.94 -33.51 -2.91 -33.48 -2.89 -33.46 -2.87 -33.43 -2.84 -33.41 -2.82 -33.38 -2.8 -33.36 -2.77 -33.33 -2.75 -33.31 -2.73 -33.29 -2.7 -33.26 -2.68 -33.24 -2.66 -33.21 -2.64 -33.19 -2.62 -33.16 -2.59 -33.14 -2.57 -33.11 -2.55 -33.09 -2.53 -33.06 -2.51 -33.04 -2.49 -33.01 -2.46 -32.99 -2.44 -32.96 -2.42 -32.94 -2.4 -32.91 -2.38 -32.89 -2.36 -32.86 -2.34 -32.84 -2.32 -32.81 -2.3 -32.79 -2.28 -32.76 -2.26 -32.74 -2.24 -32.71 -2.22 -32.69 -2.2 -32.66 -2.18 -32.64 -2.16 -32.61 -2.14 -32.59 -2.13 -32.56 -2.11 -32.54 -2.09 -32.51 -2.07 -32.49 -2.05 -32.46 -2.03 -32.44 -2.01 -32.41 -2 -32.39 -1.98 -32.36 -1.96 -32.34 -1.94 -32.31 -1.93 -32.29 -1.91 -32.27 -1.89 -32.24 -1.87 -32.22 -1.86 -32.19 -1.84 -32.17 -1.82 -32.14 -1.81 -32.12 -1.79 -32.09 -1.77 -32.07 -1.76 -32.04 -1.74 -32.02 -1.72 -31.99 -1.71 -31.97 -1.69 -31.94 -1.68 -31.92 -1.66 -31.89 -1.65 -31.87 -1.63 -31.84 -1.62 -31.82 -1.6 -31.79 -1.58 -31.77 -1.57 -31.74 -1.55 -31.72 -1.54 -31.69 -1.53 -31.67 -1.51 -31.64 -1.5 -31.62 -1.48 -31.59 -1.47 -31.57 -1.45 -31.54 -1.44 -31.52 -1.43 -31.49 -1.41 -31.47 -1.4 -31.44 -1.38 -31.42 -1.37 -31.39 -1.36 -31.37 -1.34 -31.34 -1.33 -31.32 -1.32 -31.3 -1.3 -31.27 -1.29 -31.25 -1.28 -31.22 -1.27 -31.2 -1.25 -31.17 -1.24 -31.15 -1.23 -31.12 -1.21 -31.1 -1.2 -31.07 -1.19 -31.05 -1.18 -31.02 -1.17 -31 -1.15 -30.97 -1.14 -30.95 -1.13 -30.92 -1.12 -30.9 -1.11 -30.87 -1.1 -30.85 -1.08 -30.82 -1.07 -30.8 -1.06 -30.77 -1.05 -30.75 -1.04 -30.72 -1.03 -30.7 -1.02 -30.67 -1.01 -30.65 -0.99 -30.62 -0.98 -30.6 -0.97 -30.57 -0.96 -30.55 -0.95 -30.52 -0.94 -30.5 -0.93 -30.47 -0.92 -30.45 -0.91 -30.42 -0.9 -30.4 -0.89 -30.37 -0.88 -30.35 -0.87 -30.32 -0.86 -30.3 -0.85 -30.28 -0.84 -30.25 -0.83 -30.23 -0.82 -30.2 -0.81 -30.18 -0.8 -30.15 -0.79 -30.13 -0.78 -30.1 -0.78 -30.08 -0.77 -30.05 -0.76 -30.03 -0.75 -30 -0.74 -29.98 -0.73 -29.95 -0.72 -29.93 -0.71 -29.9 -0.7 -29.88 -0.69 -29.85 -0.69 -29.83 -0.68 -29.8 -0.67 -29.78 -0.66 -29.75 -0.65 -29.73 -0.64 -29.7 -0.64 -29.68 -0.63 -29.65 -0.62 -29.63 -0.61 -29.6 -0.6 -29.58 -0.6 -29.55 -0.59 -29.53 -0.58 -29.5 -0.57 -29.48 -0.56 -29.45 -0.56 -29.43 -0.55 -29.4 -0.54 -29.38 -0.53 -29.35 -0.53 -29.33 -0.52 -29.3 -0.51 -29.28 -0.51 -29.26 -0.5 -29.23 -0.49 -29.21 -0.48 -29.18 -0.48 -29.16 -0.47 -29.13 -0.46 -29.11 -0.46 -29.08 -0.45 -29.06 -0.44 -29.03 -0.44 -29.01 -0.43 -28.98 -0.42 -28.96 -0.42 -28.93 -0.41 -28.91 -0.4 -28.88 -0.4 -28.86 -0.39 -28.83 -0.38 -28.81 -0.38 -28.78 -0.37 -28.76 -0.36 -28.73 -0.36 -28.71 -0.35 -28.68 -0.35 -28.66 -0.34 -28.63 -0.33 -28.61 -0.33 -28.58 -0.32 -28.56 -0.32 -28.53 -0.31 -28.51 -0.3 -28.48 -0.3 -28.46 -0.29 -28.43 -0.29 -28.41 -0.28 -28.38 -0.28 -28.36 -0.27 -28.33 -0.27 -28.31 -0.26 -28.28 -0.25 -28.26 -0.25 -28.24 -0.24 -28.21 -0.24 -28.19 -0.23 -28.16 -0.23 -28.14 -0.22 -28.11 -0.22 -28.09 -0.21 -28.06 -0.21 -28.04 -0.2 -28.01 -0.2 -27.99 -0.19 -27.96 -0.19 -27.94 -0.18 -27.91 -0.18 -27.89 -0.17 -27.86 -0.17 -27.84 -0.16 -27.81 -0.16 -27.79 -0.15 -27.76 -0.15 -27.74 -0.15 -27.71 -0.14 -27.69 -0.14 -27.66 -0.13 -27.64 -0.13 -27.61 -0.12 -27.59 -0.12 -27.56 -0.11 -27.54 -0.11 -27.51 -0.11 -27.49 -0.1 -27.46 -0.1 -27.44 -0.09 -27.41 -0.09 -27.39 -0.08 -27.36 -0.08 -27.34 -0.08 -27.31 -0.07 -27.29 -0.07 -27.26 -0.06 -27.24 -0.06 -27.22 -0.06 -27.19 -0.05 -27.17 -0.05 -27.14 -0.04 -27.12 -0.04 -27.09 -0.04 -27.07 -0.03 -27.04 -0.03 -27.02 -0.03 -26.99 -0.02 -26.97 -0.02 -26.94 -0.01 -26.92 -0.01 -26.89 -0.01 -26.87 -0 -26.84 -0 -26.82 0 -26.79 0.01 -26.77 0.01 -26.74 0.01 -26.72 0.02 -26.69 0.02 -26.67 0.02 -26.64 0.03 -26.62 0.03 -26.59 0.03 -26.57 0.04 -26.54 0.04 -26.52 0.04 -26.49 0.05 -26.47 0.05 -26.44 0.05 -26.42 0.06 -26.39 0.06 -26.37 0.06 -26.34 0.07 -26.32 0.07 -26.29 0.07 -26.27 0.07 -26.25 0.08 -26.22 0.08 -26.2 0.08 -26.17 0.09 -26.15 0.09 -26.12 0.09 -26.1 0.09 -26.07 0.1 -26.05 0.1 -26.02 0.1 -26 0.11 -25.97 0.11 -25.95 0.11 -25.92 0.11 -25.9 0.12 -25.87 0.12 -25.85 0.12 -25.82 0.12 -25.8 0.13 -25.77 0.13 -25.75 0.13 -25.72 0.13 -25.7 0.14 -25.67 0.14 -25.65 0.14 -25.62 0.14 -25.6 0.15 -25.57 0.15 -25.55 0.15 -25.52 0.15 -25.5 0.16 -25.47 0.16 -25.45 0.16 -25.42 0.16 -25.4 0.17 -25.37 0.17 -25.35 0.17 -25.32 0.17 -25.3 0.18 -25.27 0.18 -25.25 0.18 -25.23 0.18 -25.2 0.18 -25.18 0.19 -25.15 0.19 -25.13 0.19 -25.1 0.19 -25.08 0.19 -25.05 0.2 -25.03 0.2 -25 0.2 -24.98 0.2 -24.95 0.21 -24.93 0.21 -24.9 0.21 -24.88 0.21 -24.85 0.21 -24.83 0.22 -24.8 0.22 -24.78 0.22 -24.75 0.22 -24.73 0.22 -24.7 0.22 -24.68 0.23 -24.65 0.23 -24.63 0.23 -24.6 0.23 -24.58 0.23 -24.55 0.24 -24.53 0.24 -24.5 0.24 -24.48 0.24 -24.45 0.24 -24.43 0.24 -24.4 0.25 -24.38 0.25 -24.35 0.25 -24.33 0.25 -24.3 0.25 -24.28 0.25 -24.25 0.26 -24.23 0.26 -24.21 0.26 -24.18 0.26 -24.16 0.26 -24.13 0.26 -24.11 0.27 -24.08 0.27 -24.06 0.27 -24.03 0.27 -24.01 0.27 -23.98 0.27 -23.96 0.28 -23.93 0.28 -23.91 0.28 -23.88 0.28 -23.86 0.28 -23.83 0.28 -23.81 0.28 -23.78 0.29 -23.76 0.29 -23.73 0.29 -23.71 0.29 -23.68 0.29 -23.66 0.29 -23.63 0.29 -23.61 0.3 -23.58 0.3 -23.56 0.3 -23.53 0.3 -23.51 0.3 -23.48 0.3 -23.46 0.3 -23.43 0.3 -23.41 0.31 -23.38 0.31 -23.36 0.31 -23.33 0.31 -23.31 0.31 -23.28 0.31 -23.26 0.31 -23.23 0.32 -23.21 0.32 -23.19 0.32 -23.16 0.32 -23.14 0.32 -23.11 0.32 -23.09 0.32 -23.06 0.32 -23.04 0.32 -23.01 0.33 -22.99 0.33 -22.96 0.33 -22.94 0.33 -22.91 0.33 -22.89 0.33 -22.86 0.33 -22.84 0.33 -22.81 0.33 -22.79 0.34 -22.76 0.34 -22.74 0.34 -22.71 0.34 -22.69 0.34 -22.66 0.34 -22.64 0.34 -22.61 0.34 -22.59 0.34 -22.56 0.35 -22.54 0.35 -22.51 0.35 -22.49 0.35 -22.46 0.35 -22.44 0.35 -22.41 0.35 -22.39 0.35 -22.36 0.35 -22.34 0.35 -22.31 0.35 -22.29 0.36 -22.26 0.36 -22.24 0.36 -22.21 0.36 -22.19 0.36 -22.17 0.36 -22.14 0.36 -22.12 0.36 -22.09 0.36 -22.07 0.36 -22.04 0.36 -22.02 0.37 -21.99 0.37 -21.97 0.37 -21.94 0.37 -21.92 0.37 -21.89 0.37 -21.87 0.37 -21.84 0.37 -21.82 0.37 -21.79 0.37 -21.77 0.37 -21.74 0.37 -21.72 0.38 -21.69 0.38 -21.67 0.38 -21.64 0.38 -21.62 0.38 -21.59 0.38 -21.57 0.38 -21.54 0.38 -21.52 0.38 -21.49 0.38 -21.47 0.38 -21.44 0.38 -21.42 0.38 -21.39 0.39 -21.37 0.39 -21.34 0.39 -21.32 0.39 -21.29 0.39 -21.27 0.39 -21.24 0.39 -21.22 0.39 -21.2 0.39 -21.17 0.39 -21.15 0.39 -21.12 0.39 -21.1 0.39 -21.07 0.39 -21.05 0.4 -21.02 0.4 -21 0.4 -20.97 0.4 -20.95 0.4 -20.92 0.4 -20.9 0.4 -20.87 0.4 -20.85 0.4 -20.82 0.4 -20.8 0.4 -20.77 0.4 -20.75 0.4 -20.72 0.4 -20.7 0.4 -20.67 0.4 -20.65 0.41 -20.62 0.41 -20.6 0.41 -20.57 0.41 -20.55 0.41 -20.52 0.41 -20.5 0.41 -20.47 0.41 -20.45 0.41 -20.42 0.41 -20.4 0.41 -20.37 0.41 -20.35 0.41 -20.32 0.41 -20.3 0.41 -20.27 0.41 -20.25 0.41 -20.22 0.41 -20.2 0.42 -20.18 0.42 -20.15 0.42 -20.13 0.42 -20.1 0.42 -20.08 0.42 -20.05 0.42 -20.03 0.42 -20 0.42 -19.98 0.42 -19.95 0.42 -19.93 0.42 -19.9 0.42 -19.88 0.42 -19.85 0.42 -19.83 0.42 -19.8 0.42 -19.78 0.42 -19.75 0.42 -19.73 0.42 -19.7 0.42 -19.68 0.43 -19.65 0.43 -19.63 0.43 -19.6 0.43 -19.58 0.43 -19.55 0.43 -19.53 0.43 -19.5 0.43 -19.48 0.43 -19.45 0.43 -19.43 0.43 -19.4 0.43 -19.38 0.43 -19.35 0.43 -19.33 0.43 -19.3 0.43 -19.28 0.43 -19.25 0.43 -19.23 0.43 -19.2 0.43 -19.18 0.43 -19.16 0.43 -19.13 0.43 -19.11 0.43 -19.08 0.43 -19.06 0.44 -19.03 0.44 -19.01 0.44 -18.98 0.44 -18.96 0.44 -18.93 0.44 -18.91 0.44 -18.88 0.44 -18.86 0.44 -18.83 0.44 -18.81 0.44 -18.78 0.44 -18.76 0.44 -18.73 0.44 -18.71 0.44 -18.68 0.44 -18.66 0.44 -18.63 0.44 -18.61 0.44 -18.58 0.44 -18.56 0.44 -18.53 0.44 -18.51 0.44 -18.48 0.44 -18.46 0.44 -18.43 0.44 -18.41 0.44 -18.38 0.44 -18.36 0.44 -18.33 0.45 -18.31 0.45 -18.28 0.45 -18.26 0.45 -18.23 0.45 -18.21 0.45 -18.18 0.45 -18.16 0.45 -18.14 0.45 -18.11 0.45 -18.09 0.45 -18.06 0.45 -18.04 0.45 -18.01 0.45 -17.99 0.45 -17.96 0.45 -17.94 0.45 -17.91 0.45 -17.89 0.45 -17.86 0.45 -17.84 0.45 -17.81 0.45 -17.79 0.45 -17.76 0.45 -17.74 0.45 -17.71 0.45 -17.69 0.45 -17.66 0.45 -17.64 0.45 -17.61 0.45 -17.59 0.45 -17.56 0.45 -17.54 0.45 -17.51 0.45 -17.49 0.45 -17.46 0.45 -17.44 0.46 -17.41 0.46 -17.39 0.46 -17.36 0.46 -17.34 0.46 -17.31 0.46 -17.29 0.46 -17.26 0.46 -17.24 0.46 -17.21 0.46 -17.19 0.46 -17.16 0.46 -17.14 0.46 -17.12 0.46 -17.09 0.46 -17.07 0.46 -17.04 0.46 -17.02 0.46 -16.99 0.46 -16.97 0.46 -16.94 0.46 -16.92 0.46 -16.89 0.46 -16.87 0.46 -16.84 0.46 -16.82 0.46 -16.79 0.46 -16.77 0.46 -16.74 0.46 -16.72 0.46 -16.69 0.46 -16.67 0.46 -16.64 0.46 -16.62 0.46 -16.59 0.46 -16.57 0.46 -16.54 0.46 -16.52 0.46 -16.49 0.46 -16.47 0.46 -16.44 0.46 -16.42 0.46 -16.39 0.46 -16.37 0.46 -16.34 0.46 -16.32 0.46 -16.29 0.46 -16.27 0.46 -16.24 0.46 -16.22 0.47 -16.19 0.47 -16.17 0.47 -16.15 0.47 -16.12 0.47 -16.1 0.47 -16.07 0.47 -16.05 0.47 -16.02 0.47 -16 0.47 -15.97 0.47 -15.95 0.47 -15.92 0.47 -15.9 0.47 -15.87 0.47 -15.85 0.47 -15.82 0.47 -15.8 0.47 -15.77 0.47 -15.75 0.47 -15.72 0.47 -15.7 0.47 -15.67 0.47 -15.65 0.47 -15.62 0.47 -15.6 0.47 -15.57 0.47 -15.55 0.47 -15.52 0.47 -15.5 0.47 -15.47 0.47 -15.45 0.47 -15.42 0.47 -15.4 0.47 -15.37 0.47 -15.35 0.47 -15.32 0.47 -15.3 0.47 -15.27 0.47 -15.25 0.47 -15.22 0.47 -15.2 0.47 -15.17 0.47 -15.15 0.47 -15.13 0.47 -15.1 0.47 -15.08 0.47 -15.05 0.47 -15.03 0.47 -15 0.47 -14.98 0.47 -14.95 0.47 -14.93 0.47 -14.9 0.47 -14.88 0.47 -14.85 0.47 -14.83 0.47 -14.8 0.47 -14.78 0.47 -14.75 0.47 -14.73 0.47 -14.7 0.47 -14.68 0.47 -14.65 0.47 -14.63 0.47 -14.6 0.47 -14.58 0.47 -14.55 0.47 -14.53 0.47 -14.5 0.47 -14.48 0.47 -14.45 0.47 -14.43 0.47 -14.4 0.47 -14.38 0.48 -14.35 0.48 -14.33 0.48 -14.3 0.48 -14.28 0.48 -14.25 0.48 -14.23 0.48 -14.2 0.48 -14.18 0.48 -14.15 0.48 -14.13 0.48 -14.11 0.48 -14.08 0.48 -14.06 0.48 -14.03 0.48 -14.01 0.48 -13.98 0.48 -13.96 0.48 -13.93 0.48 -13.91 0.48 -13.88 0.48 -13.86 0.48 -13.83 0.48 -13.81 0.48 -13.78 0.48 -13.76 0.48 -13.73 0.48 -13.71 0.48 -13.68 0.48 -13.66 0.48 -13.63 0.48 -13.61 0.48 -13.58 0.48 -13.56 0.48 -13.53 0.48 -13.51 0.48 -13.48 0.48 -13.46 0.48 -13.43 0.48 -13.41 0.48 -13.38 0.48 -13.36 0.48 -13.33 0.48 -13.31 0.48 -13.28 0.48 -13.26 0.48 -13.23 0.48 -13.21 0.48 -13.18 0.48 -13.16 0.48 -13.13 0.48 -13.11 0.48 -13.09 0.48 -13.06 0.48 -13.04 0.48 -13.01 0.48 -12.99 0.48 -12.96 0.48 -12.94 0.48 -12.91 0.48 -12.89 0.48 -12.86 0.48 -12.84 0.48 -12.81 0.48 -12.79 0.48 -12.76 0.48 -12.74 0.48 -12.71 0.48 -12.69 0.48 -12.66 0.48 -12.64 0.48 -12.61 0.48 -12.59 0.48 -12.56 0.48 -12.54 0.48 -12.51 0.48 -12.49 0.48 -12.46 0.48 -12.44 0.48 -12.41 0.48 -12.39 0.48 -12.36 0.48 -12.34 0.48 -12.31 0.48 -12.29 0.48 -12.26 0.48 -12.24 0.48 -12.21 0.48 -12.19 0.48 -12.16 0.48 -12.14 0.48 -12.12 0.48 -12.09 0.48 -12.07 0.48 -12.04 0.48 -12.02 0.48 -11.99 0.48 -11.97 0.48 -11.94 0.48 -11.92 0.48 -11.89 0.48 -11.87 0.48 -11.84 0.48 -11.82 0.48 -11.79 0.48 -11.77 0.48 -11.74 0.48 -11.72 0.48 -11.69 0.48 -11.67 0.48 -11.64 0.48 -11.62 0.48 -11.59 0.48 -11.57 0.48 -11.54 0.48 -11.52 0.48 -11.49 0.48 -11.47 0.48 -11.44 0.48 -11.42 0.48 -11.39 0.48 -11.37 0.48 -11.34 0.48 -11.32 0.48 -11.29 0.48 -11.27 0.48 -11.24 0.48 -11.22 0.48 -11.19 0.48 -11.17 0.48 -11.14 0.48 -11.12 0.48 -11.1 0.48 -11.07 0.48 -11.05 0.48 -11.02 0.48 -11 0.48 -10.97 0.48 -10.95 0.48 -10.92 0.48 -10.9 0.48 -10.87 0.48 -10.85 0.48 -10.82 0.48 -10.8 0.48 -10.77 0.48 -10.75 0.48 -10.72 0.48 -10.7 0.48 -10.67 0.48 -10.65 0.48 -10.62 0.48 -10.6 0.48 -10.57 0.48 -10.55 0.48 -10.52 0.48 -10.5 0.48 -10.47 0.48 -10.45 0.48 -10.42 0.48 -10.4 0.48 -10.37 0.48 -10.35 0.48 -10.32 0.49 -10.3 0.49 -10.27 0.49 -10.25 0.49 -10.22 0.49 -10.2 0.49 -10.17 0.49 -10.15 0.49 -10.12 0.49 -10.1 0.49 -10.08 0.49 -10.05 0.49 -10.03 0.49 -10 0.49 -9.98 0.49 -9.95 0.49 -9.93 0.49 -9.9 0.49 -9.88 0.49 -9.85 0.49 -9.83 0.49 -9.8 0.49 -9.78 0.49 -9.75 0.49 -9.73 0.49 -9.7 0.49 -9.68 0.49 -9.65 0.49 -9.63 0.49 -9.6 0.49 -9.58 0.49 -9.55 0.49 -9.53 0.49 -9.5 0.49 -9.48 0.49 -9.45 0.49 -9.43 0.49 -9.4 0.49 -9.38 0.49 -9.35 0.49 -9.33 0.49 -9.3 0.49 -9.28 0.49 -9.25 0.49 -9.23 0.49 -9.2 0.49 -9.18 0.49 -9.15 0.49 -9.13 0.49 -9.1 0.49 -9.08 0.49 -9.06 0.49 -9.03 0.49 -9.01 0.49 -8.98 0.49 -8.96 0.49 -8.93 0.49 -8.91 0.49 -8.88 0.49 -8.86 0.49 -8.83 0.49 -8.81 0.49 -8.78 0.49 -8.76 0.49 -8.73 0.49 -8.71 0.49 -8.68 0.49 -8.66 0.49 -8.63 0.49 -8.61 0.49 -8.58 0.49 -8.56 0.49 -8.53 0.49 -8.51 0.49 -8.48 0.49 -8.46 0.49 -8.43 0.49 -8.41 0.49 -8.38 0.49 -8.36 0.49 -8.33 0.49 -8.31 0.49 -8.28 0.49 -8.26 0.49 -8.23 0.49 -8.21 0.49 -8.18 0.49 -8.16 0.49 -8.13 0.49 -8.11 0.49 -8.08 0.49 -8.06 0.49 -8.04 0.49 -8.01 0.49 -7.99 0.49 -7.96 0.49 -7.94 0.49 -7.91 0.49 -7.89 0.49 -7.86 0.49 -7.84 0.49 -7.81 0.49 -7.79 0.49 -7.76 0.49 -7.74 0.49 -7.71 0.49 -7.69 0.49 -7.66 0.49 -7.64 0.49 -7.61 0.49 -7.59 0.49 -7.56 0.49 -7.54 0.49 -7.51 0.49 -7.49 0.49 -7.46 0.49 -7.44 0.49 -7.41 0.49 -7.39 0.49 -7.36 0.49 -7.34 0.49 -7.31 0.49 -7.29 0.49 -7.26 0.49 -7.24 0.49 -7.21 0.49 -7.19 0.49 -7.16 0.49 -7.14 0.49 -7.11 0.49 -7.09 0.49 -7.07 0.49 -7.04 0.49 -7.02 0.49 -6.99 0.49 -6.97 0.49 -6.94 0.49 -6.92 0.49 -6.89 0.49 -6.87 0.49 -6.84 0.49 -6.82 0.49 -6.79 0.49 -6.77 0.49 -6.74 0.49 -6.72 0.49 -6.69 0.49 -6.67 0.49 -6.64 0.49 -6.62 0.49 -6.59 0.49 -6.57 0.49 -6.54 0.49 -6.52 0.49 -6.49 0.49 -6.47 0.49 -6.44 0.49 -6.42 0.49 -6.39 0.49 -6.37 0.49 -6.34 0.49 -6.32 0.49 -6.29 0.49 -6.27 0.49 -6.24 0.49 -6.22 0.49 -6.19 0.49 -6.17 0.49 -6.14 0.49 -6.12 0.49 -6.09 0.49 -6.07 0.49 -6.05 0.49 -6.02 0.49 -6 0.49 -5.97 0.49 -5.95 0.49 -5.92 0.49 -5.9 0.49 -5.87 0.49 -5.85 0.49 -5.82 0.49 -5.8 0.49 -5.77 0.49 -5.75 0.49 -5.72 0.49 -5.7 0.49 -5.67 0.49 -5.65 0.49 -5.62 0.49 -5.6 0.49 -5.57 0.49 -5.55 0.49 -5.52 0.49 -5.5 0.49 -5.47 0.49 -5.45 0.49 -5.42 0.49 -5.4 0.49 -5.37 0.49 -5.35 0.49 -5.32 0.49 -5.3 0.49 -5.27 0.49 -5.25 0.49 -5.22 0.49 -5.2 0.49 -5.17 0.49 -5.15 0.49 -5.12 0.49 -5.1 0.49 -5.07 0.49 -5.05 0.49 -5.03 0.49 -5 0.49 -4.98 0.49 -4.95 0.49 -4.93 0.49 -4.9 0.49 -4.88 0.49 -4.85 0.49 -4.83 0.49 -4.8 0.49 -4.78 0.49 -4.75 0.49 -4.73 0.49 -4.7 0.49 -4.68 0.49 -4.65 0.49 -4.63 0.49 -4.6 0.49 -4.58 0.49 -4.55 0.49 -4.53 0.49 -4.5 0.49 -4.48 0.49 -4.45 0.49 -4.43 0.49 -4.4 0.49 -4.38 0.49 -4.35 0.49 -4.33 0.49 -4.3 0.49 -4.28 0.49 -4.25 0.49 -4.23 0.49 -4.2 0.49 -4.18 0.49 -4.15 0.49 -4.13 0.49 -4.1 0.49 -4.08 0.49 -4.05 0.49 -4.03 0.49 -4.01 0.49 -3.98 0.49 -3.96 0.49 -3.93 0.49 -3.91 0.49 -3.88 0.49 -3.86 0.49 -3.83 0.49 -3.81 0.49 -3.78 0.49 -3.76 0.49 -3.73 0.49 -3.71 0.49 -3.68 0.49 -3.66 0.49 -3.63 0.49 -3.61 0.49 -3.58 0.49 -3.56 0.49 -3.53 0.49 -3.51 0.49 -3.48 0.49 -3.46 0.49 -3.43 0.49 -3.41 0.49 -3.38 0.49 -3.36 0.49 -3.33 0.49 -3.31 0.49 -3.28 0.49 -3.26 0.49 -3.23 0.49 -3.21 0.49 -3.18 0.49 -3.16 0.49 -3.13 0.49 -3.11 0.49 -3.08 0.49 -3.06 0.49 -3.03 0.49 -3.01 0.49 -2.99 0.49 -2.96 0.49 -2.94 0.49 -2.91 0.49 -2.89 0.49 -2.86 0.49 -2.84 0.49 -2.81 0.49 -2.79 0.49 -2.76 0.49 -2.74 0.49 -2.71 0.49 -2.69 0.49 -2.66 0.49 -2.64 0.49 -2.61 0.48 -2.59 0.48 -2.56 0.48 -2.54 0.48 -2.51 0.48 -2.49 0.48 -2.46 0.48 -2.44 0.48 -2.41 0.48 -2.39 0.48 -2.36 0.48 -2.34 0.48 -2.31 0.48 -2.29 0.48 -2.26 0.48 -2.24 0.48 -2.21 0.48 -2.19 0.48 -2.16 0.48 -2.14 0.48 -2.11 0.48 -2.09 0.48 -2.06 0.48 -2.04 0.48 -2.02 0.48 -1.99 0.48 -1.97 0.48 -1.94 0.48 -1.92 0.48 -1.89 0.48 -1.87 0.48 -1.84 0.48 -1.82 0.48 -1.79 0.48 -1.77 0.48 -1.74 0.48 -1.72 0.48 -1.69 0.48 -1.67 0.48 -1.64 0.48 -1.62 0.48 -1.59 0.48 -1.57 0.48 -1.54 0.48 -1.52 0.48 -1.49 0.48 -1.47 0.48 -1.44 0.48 -1.42 0.48 -1.39 0.48 -1.37 0.48 -1.34 0.48 -1.32 0.48 -1.29 0.48 -1.27 0.48 -1.24 0.48 -1.22 0.48 -1.19 0.48 -1.17 0.48 -1.14 0.48 -1.12 0.48 -1.09 0.48 -1.07 0.48 -1.04 0.48 -1.02 0.48 -1 0.48 -0.97 0.48 -0.95 0.48 -0.92 0.48 -0.9 0.48 -0.87 0.48 -0.85 0.48 -0.82 0.48 -0.8 0.48 -0.77 0.48 -0.75 0.48 -0.72 0.48 -0.7 0.48 -0.67 0.48 -0.65 0.48 -0.62 0.48 -0.6 0.48 -0.57 0.48 -0.55 0.48 -0.52 0.48 -0.5 0.48 -0.47 0.48 -0.45 0.48 -0.42 0.48 -0.4 0.48 -0.37 0.48 -0.35 0.48 -0.32 0.48 -0.3 0.48 -0.27 0.48 -0.25 0.48 -0.22 0.48 -0.2 0.48 -0.17 0.48 -0.15 0.48 -0.12 0.48 -0.1 0.48 -0.07 0.48 -0.05 0.48 -0.02 0.48 0 0.48 0.02 0.48 0.05 0.48 0.07 0.48 0.1 0.48 0.12 0.48 0.15 0.48 0.17 0.48 0.2 0.48 0.22 0.48 0.25 0.48 0.27 0.48 0.3 0.48 0.32 0.48 0.35 0.48 0.37 0.48 0.4 0.48 0.42 0.48 0.45 0.48 0.47 0.48 0.5 0.48 0.52 0.48 0.55 0.48 0.57 0.48 0.6 0.48 0.62 0.48 0.65 0.48 0.67 0.48 0.7 0.48 0.72 0.48 0.75 0.48 0.77 0.48 0.8 0.48 0.82 0.48 0.85 0.48 0.87 0.48 0.9 0.48 0.92 0.48 0.95 0.48 0.97 0.48 1 0.48 1.02 0.48 1.04 0.48 1.07 0.48 1.09 0.48 1.12 0.48 1.14 0.48 1.17 0.48 1.19 0.48 1.22 0.48 1.24 0.48 1.27 0.48 1.29 0.48 1.32 0.48 1.34 0.48 1.37 0.48 1.39 0.48 1.42 0.48 1.44 0.48 1.47 0.48 1.49 0.48 1.52 0.48 1.54 0.48 1.57 0.48 1.59 0.48 1.62 0.48 1.64 0.48 1.67 0.48 1.69 0.48 1.72 0.48 1.74 0.48 1.77 0.48 1.79 0.48 1.82 0.48 1.84 0.48 1.87 0.48 1.89 0.48 1.92 0.48 1.94 0.48 1.97 0.48 1.99 0.48 2.02 0.48 2.04 0.48 2.06 0.48 2.09 0.48 2.11 0.48 2.14 0.48 2.16 0.48 2.19 0.48 2.21 0.48 2.24 0.48 2.26 0.48 2.29 0.48 2.31 0.48 2.34 0.48 2.36 0.48 2.39 0.48 2.41 0.48 2.44 0.48 2.46 0.48 2.49 0.48 2.51 0.48 2.54 0.48 2.56 0.48 2.59 0.48 2.61 0.48 2.64 0.48 2.66 0.48 2.69 0.48 2.71 0.48 2.74 0.48 2.76 0.48 2.79 0.48 2.81 0.48 2.84 0.48 2.86 0.48 2.89 0.48 2.91 0.48 2.94 0.48 2.96 0.48 2.99 0.48 3.01 0.48 3.03 0.48 3.06 0.48 3.08 0.48 3.11 0.48 3.13 0.48 3.16 0.48 3.18 0.48 3.21 0.48 3.23 0.48 3.26 0.48 3.28 0.48 3.31 0.48 3.33 0.48 3.36 0.48 3.38 0.48 3.41 0.48 3.43 0.48 3.46 0.48 3.48 0.48 3.51 0.48 3.53 0.48 3.56 0.48 3.58 0.48 3.61 0.48 3.63 0.48 3.66 0.48 3.68 0.48 3.71 0.48 3.73 0.48 3.76 0.48 3.78 0.48 3.81 0.48 3.83 0.48 3.86 0.48 3.88 0.48 3.91 0.48 3.93 0.48 3.96 0.48 3.98 0.48 4.01 0.48 4.03 0.48 4.05 0.48 4.08 0.48 4.1 0.48 4.13 0.48 4.15 0.48 4.18 0.48 4.2 0.48 4.23 0.48 4.25 0.48 4.28 0.48 4.3 0.48 4.33 0.48 4.35 0.48 4.38 0.48 4.4 0.48 4.43 0.48 4.45 0.48 4.48 0.48 4.5 0.48 4.53 0.48 4.55 0.48 4.58 0.48 4.6 0.48 4.63 0.48 4.65 0.48 4.68 0.48 4.7 0.48 4.73 0.48 4.75 0.48 4.78 0.48 4.8 0.48 4.83 0.48 4.85 0.48 4.88 0.48 4.9 0.48 4.93 0.48 4.95 0.48 4.98 0.48 5 0.48 5.03 0.48 5.05 0.48 5.07 0.48 5.1 0.48 5.12 0.48 5.15 0.48 5.17 0.48 5.2 0.48 5.22 0.48 5.25 0.48 5.27 0.48 5.3 0.48 5.32 0.48 5.35 0.48 5.37 0.48 5.4 0.48 5.42 0.48 5.45 0.48 5.47 0.48 5.5 0.48 5.52 0.48 5.55 0.48 5.57 0.48 5.6 0.48 5.62 0.48 5.65 0.48 5.67 0.48 5.7 0.48 5.72 0.48 5.75 0.48 5.77 0.48 5.8 0.48 5.82 0.48 5.85 0.48 5.87 0.48 5.9 0.48 5.92 0.48 5.95 0.48 5.97 0.48 6 0.48 6.02 0.48 6.05 0.48 6.07 0.48 6.09 0.48 6.12 0.48 6.14 0.48 6.17 0.48 6.19 0.48 6.22 0.48 6.24 0.48 6.27 0.48 6.29 0.48 6.32 0.48 6.34 0.48 6.37 0.48 6.39 0.48 6.42 0.48 6.44 0.48 6.47 0.48 6.49 0.48 6.52 0.48 6.54 0.48 6.57 0.48 6.59 0.48 6.62 0.48 6.64 0.48 6.67 0.48 6.69 0.48 6.72 0.48 6.74 0.48 6.77 0.48 6.79 0.48 6.82 0.48 6.84 0.48 6.87 0.48 6.89 0.48 6.92 0.48 6.94 0.48 6.97 0.48 6.99 0.48 7.02 0.48 7.04 0.48 7.07 0.48 7.09 0.48 7.11 0.48 7.14 0.48 7.16 0.48 7.19 0.48 7.21 0.48 7.24 0.48 7.26 0.48 7.29 0.48 7.31 0.48 7.34 0.48 7.36 0.48 7.39 0.48 7.41 0.48 7.44 0.48 7.46 0.48 7.49 0.48 7.51 0.48 7.54 0.48 7.56 0.48 7.59 0.48 7.61 0.48 7.64 0.48 7.66 0.48 7.69 0.48 7.71 0.48 7.74 0.48 7.76 0.48 7.79 0.48 7.81 0.48 7.84 0.48 7.86 0.48 7.89 0.48 7.91 0.48 7.94 0.48 7.96 0.48 7.99 0.48 8.01 0.48 8.04 0.48 8.06 0.48 8.08 0.48 8.11 0.48 8.13 0.48 8.16 0.48 8.18 0.48 8.21 0.48 8.23 0.48 8.26 0.48 8.28 0.48 8.31 0.48 8.33 0.48 8.36 0.48 8.38 0.48 8.41 0.48 8.43 0.48 8.46 0.48 8.48 0.48 8.51 0.48 8.53 0.48 8.56 0.48 8.58 0.48 8.61 0.48 8.63 0.48 8.66 0.48 8.68 0.48 8.71 0.48 8.73 0.48 8.76 0.48 8.78 0.48 8.81 0.48 8.83 0.48 8.86 0.48 8.88 0.48 8.91 0.48 8.93 0.48 8.96 0.48 8.98 0.48 9.01 0.48 9.03 0.48 9.06 0.48 9.08 0.48 9.1 0.48 9.13 0.48 9.15 0.48 9.18 0.48 9.2 0.48 9.23 0.48 9.25 0.48 9.28 0.48 9.3 0.48 9.33 0.48 9.35 0.48 9.38 0.48 9.4 0.48 9.43 0.48 9.45 0.48 9.48 0.48 9.5 0.48 9.53 0.48 9.55 0.48 9.58 0.48 9.6 0.48 9.63 0.48 9.65 0.48 9.68 0.48 9.7 0.48 9.73 0.48 9.75 0.48 9.78 0.48 9.8 0.48 9.83 0.48 9.85 0.48 9.88 0.48 9.9 0.48 9.93 0.48 9.95 0.48 9.98 0.48 10 0.48 10.03 0.48 10.05 0.48 10.08 0.48 10.1 0.48 10.12 0.48 10.15 0.48 10.17 0.48 10.2 0.48 10.22 0.48 10.25 0.48 10.27 0.48 10.3 0.48 10.32 0.48 10.35 0.48 10.37 0.48 10.4 0.48 10.42 0.48 10.45 0.48 10.47 0.48 10.5 0.48 10.52 0.48 10.55 0.48 10.57 0.48 10.6 0.48 10.62 0.48 10.65 0.48 10.67 0.48 10.7 0.48 10.72 0.48 10.75 0.48 10.77 0.48 10.8 0.48 10.82 0.48 10.85 0.48 10.87 0.48 10.9 0.48 10.92 0.48 10.95 0.48 10.97 0.48 11 0.48 11.02 0.48 11.05 0.48 11.07 0.47 11.1 0.47 11.12 0.47 11.14 0.47 11.17 0.47 11.19 0.47 11.22 0.47 11.24 0.47 11.27 0.47 11.29 0.47 11.32 0.47 11.34 0.47 11.37 0.47 11.39 0.47 11.42 0.47 11.44 0.47 11.47 0.47 11.49 0.47 11.52 0.47 11.54 0.47 11.57 0.47 11.59 0.47 11.62 0.47 11.64 0.47 11.67 0.47 11.69 0.47 11.72 0.47 11.74 0.47 11.77 0.47 11.79 0.47 11.82 0.47 11.84 0.47 11.87 0.47 11.89 0.47 11.92 0.47 11.94 0.47 11.97 0.47 11.99 0.47 12.02 0.47 12.04 0.47 12.07 0.47 12.09 0.47 12.12 0.47 12.14 0.47 12.16 0.47 12.19 0.47 12.21 0.47 12.24 0.47 12.26 0.47 12.29 0.47 12.31 0.47 12.34 0.47 12.36 0.47 12.39 0.47 12.41 0.47 12.44 0.47 12.46 0.47 12.49 0.47 12.51 0.47 12.54 0.47 12.56 0.47 12.59 0.47 12.61 0.47 12.64 0.47 12.66 0.47 12.69 0.47 12.71 0.47 12.74 0.47 12.76 0.47 12.79 0.47 12.81 0.47 12.84 0.47 12.86 0.47 12.89 0.47 12.91 0.47 12.94 0.47 12.96 0.47 12.99 0.47 13.01 0.47 13.04 0.47 13.06 0.47 13.09 0.47 13.11 0.47 13.13 0.47 13.16 0.47 13.18 0.47 13.21 0.47 13.23 0.47 13.26 0.47 13.28 0.47 13.31 0.47 13.33 0.47 13.36 0.47 13.38 0.47 13.41 0.47 13.43 0.47 13.46 0.47 13.48 0.47 13.51 0.47 13.53 0.47 13.56 0.47 13.58 0.47 13.61 0.47 13.63 0.47 13.66 0.47 13.68 0.47 13.71 0.47 13.73 0.47 13.76 0.47 13.78 0.47 13.81 0.47 13.83 0.47 13.86 0.47 13.88 0.47 13.91 0.47 13.93 0.47 13.96 0.47 13.98 0.47 14.01 0.47 14.03 0.47 14.06 0.47 14.08 0.47 14.11 0.47 14.13 0.47 14.15 0.47 14.18 0.47 14.2 0.47 14.23 0.47 14.25 0.47 14.28 0.47 14.3 0.47 14.33 0.47 14.35 0.47 14.38 0.47 14.4 0.47 14.43 0.47 14.45 0.47 14.48 0.47 14.5 0.47 14.53 0.47 14.55 0.47 14.58 0.47 14.6 0.47 14.63 0.47 14.65 0.47 14.68 0.47 14.7 0.47 14.73 0.47 14.75 0.47 14.78 0.47 14.8 0.47 14.83 0.47 14.85 0.47 14.88 0.47 14.9 0.47 14.93 0.47 14.95 0.47 14.98 0.47 15 0.47 15.03 0.47 15.05 0.47 15.08 0.47 15.1 0.47 15.13 0.47 15.15 0.47 15.17 0.47 15.2 0.47 15.22 0.47 15.25 0.47 15.27 0.47 15.3 0.47 15.32 0.47 15.35 0.47 15.37 0.47 15.4 0.47 15.42 0.47 15.45 0.47 15.47 0.47 15.5 0.47 15.52 0.47 15.55 0.47 15.57 0.47 15.6 0.47 15.62 0.47 15.65 0.47 15.67 0.47 15.7 0.47 15.72 0.47 15.75 0.47 15.77 0.47 15.8 0.47 15.82 0.47 15.85 0.47 15.87 0.47 15.9 0.47 15.92 0.47 15.95 0.47 15.97 0.47 16 0.47 16.02 0.47 16.05 0.47 16.07 0.47 16.1 0.47 16.12 0.47 16.15 0.47 16.17 0.47 16.19 0.47 16.22 0.47 16.24 0.47 16.27 0.47 16.29 0.47 16.32 0.47 16.34 0.47 16.37 0.47 16.39 0.47 16.42 0.47 16.44 0.47 16.47 0.47 16.49 0.47 16.52 0.47 16.54 0.47 16.57 0.47 16.59 0.47 16.62 0.47 16.64 0.47 16.67 0.47 16.69 0.47 16.72 0.47 16.74 0.47 16.77 0.47 16.79 0.47 16.82 0.47 16.84 0.47 16.87 0.47 16.89 0.47 16.92 0.47 16.94 0.47 16.97 0.47 16.99 0.47 17.02 0.47 17.04 0.47 17.07 0.47 17.09 0.47 17.12 0.47 17.14 0.47 17.16 0.47 17.19 0.47 17.21 0.47 17.24 0.47 17.26 0.47 17.29 0.47 17.31 0.47 17.34 0.47 17.36 0.47 17.39 0.47 17.41 0.47 17.44 0.47 17.46 0.47 17.49 0.47 17.51 0.47 17.54 0.47 17.56 0.47 17.59 0.47 17.61 0.47 17.64 0.47 17.66 0.47 17.69 0.47 17.71 0.47 17.74 0.47 17.76 0.47 17.79 0.47 17.81 0.47 17.84 0.47 17.86 0.47 17.89 0.47 17.91 0.47 17.94 0.47 17.96 0.47 17.99 0.47 18.01 0.47 18.04 0.47 18.06 0.47 18.09 0.47 18.11 0.47 18.14 0.47 18.16 0.47 18.18 0.47 18.21 0.47 18.23 0.47 18.26 0.47 18.28 0.47 18.31 0.47 18.33 0.47 18.36 0.47 18.38 0.47 18.41 0.47 18.43 0.47 18.46 0.47 18.48 0.47 18.51 0.47 18.53 0.47 18.56 0.47 18.58 0.47 18.61 0.47 18.63 0.47 18.66 0.47 18.68 0.47 18.71 0.47 18.73 0.47 18.76 0.47 18.78 0.47 18.81 0.47 18.83 0.47 18.86 0.47 18.88 0.47 18.91 0.47 18.93 0.47 18.96 0.47 18.98 0.47 19.01 0.47 19.03 0.47 19.06 0.47 19.08 0.47 19.11 0.47 19.13 0.47 19.16 0.47 19.18 0.47 19.2 0.47 19.23 0.47 19.25 0.47 19.28 0.47 19.3 0.47 19.33 0.47 19.35 0.47 19.38 0.47 19.4 0.47 19.43 0.47 19.45 0.47 19.48 0.47 19.5 0.47 19.53 0.47 19.55 0.47 19.58 0.47 19.6 0.47 19.63 0.47 19.65 0.47 19.68 0.47 19.7 0.47 19.73 0.47 19.75 0.47 19.78 0.47 19.8 0.47 19.83 0.47 19.85 0.47 19.88 0.47 19.9 0.47 19.93 0.47 19.95 0.47 19.98 0.47 20 0.47 20.03 0.47 20.05 0.47 20.08 0.47 20.1 0.47 20.13 0.47 20.15 0.47 20.18 0.47 20.2 0.47 20.22 0.47 20.25 0.47 20.27 0.47 20.3 0.47 20.32 0.47 20.35 0.47 20.37 0.47 20.4 0.47 20.42 0.47 20.45 0.47 20.47 0.47 20.5 0.47 20.52 0.47 20.55 0.47 20.57 0.47 20.6 0.47 20.62 0.47 20.65 0.47 20.67 0.47 20.7 0.47 20.72 0.47 20.75 0.47 20.77 0.47 20.8 0.47 20.82 0.47 20.85 0.47 20.87 0.47 20.9 0.47 20.92 0.47 20.95 0.47 20.97 0.47 21 0.47 21.02 0.47 21.05 0.47 21.07 0.47 21.1 0.47 21.12 0.47 21.15 0.47 21.17 0.47 21.2 0.47 21.22 0.47 21.24 0.47 21.27 0.47 21.29 0.47 21.32 0.47 21.34 0.47 21.37 0.47 21.39 0.47 21.42 0.47 21.44 0.47 21.47 0.47 21.49 0.47 21.52 0.47 21.54 0.47 21.57 0.47 21.59 0.47 21.62 0.47 21.64 0.47 21.67 0.47 21.69 0.47 21.72 0.47 21.74 0.47 21.77 0.47 21.79 0.47 21.82 0.47 21.84 0.47 21.87 0.47 21.89 0.47 21.92 0.47 21.94 0.47 21.97 0.47 21.99 0.47 22.02 0.47 22.04 0.47 22.07 0.47 22.09 0.47 22.12 0.47 22.14 0.47 22.17 0.47 22.19 0.47 22.21 0.47 22.24 0.47 22.26 0.47 22.29 0.47 22.31 0.47 22.34 0.47 22.36 0.47 22.39 0.47 22.41 0.47 22.44 0.47 22.46 0.47 22.49 0.47 22.51 0.47 22.54 0.47 22.56 0.47 22.59 0.47 22.61 0.47 22.64 0.47 22.66 0.47 22.69 0.47 22.71 0.47 22.74 0.47 22.76 0.47 22.79 0.47 22.81 0.47 22.84 0.47 22.86 0.47 22.89 0.47 22.91 0.47 22.94 0.47 22.96 0.47 22.99 0.47 23.01 0.47 23.04 0.47 23.06 0.47 23.09 0.47 23.11 0.47 23.14 0.47 23.16 0.47 23.19 0.47 23.21 0.47 23.23 0.47 23.26 0.47 23.28 0.47 23.31 0.47 23.33 0.47 23.36 0.47 23.38 0.47 23.41 0.47 23.43 0.47 23.46 0.47 23.48 0.47 23.51 0.47 23.53 0.47 23.56 0.47 23.58 0.47 23.61 0.47 23.63 0.47 23.66 0.47 23.68 0.47 23.71 0.47 23.73 0.47 23.76 0.47 23.78 0.47 23.81 0.47 23.83 0.47 23.86 0.47 23.88 0.47 23.91 0.47 23.93 0.47 23.96 0.47 23.98 0.47 24.01 0.47 24.03 0.47 24.06 0.47 24.08 0.47 24.11 0.47 24.13 0.47 24.16 0.47 24.18 0.47 24.21 0.47 24.23 0.47 24.25 0.47 24.28 0.47 24.3 0.47 24.33 0.47 24.35 0.47 24.38 0.47 24.4 0.47 24.43 0.47 24.45 0.47 24.48 0.47 24.5 0.47 24.53 0.47 24.55 0.47 24.58 0.47 24.6 0.47 24.63 0.47 24.65 0.47 24.68 0.47 24.7 0.47 24.73 0.47 24.75 0.47 24.78 0.47 24.8 0.47 24.83 0.47 24.85 0.47 24.88 0.47 24.9 0.47 24.93 0.47 24.95 0.47 24.98 0.47 25 0.47 25.03 0.47 25.05 0.47 25.08 0.47 25.1 0.47 25.13 0.47 25.15 0.47 25.18 0.47 25.2 0.47 25.23 0.47 25.25 0.47 25.27 0.47 25.3 0.47 25.32 0.47 25.35 0.47 25.37 0.47 25.4 0.47 25.42 0.47 25.45 0.47 25.47 0.47 25.5 0.47 25.52 0.47 25.55 0.47 25.57 0.47 25.6 0.47 25.62 0.47 25.65 0.47 25.67 0.47 25.7 0.47 25.72 0.47 25.75 0.47 25.77 0.47 25.8 0.47 25.82 0.47 25.85 0.47 25.87 0.47 25.9 0.47 25.92 0.47 25.95 0.47 25.97 0.47 26 0.47 26.02 0.47 26.05 0.47 26.07 0.47 26.1 0.47 26.12 0.47 26.15 0.47 26.17 0.47 26.2 0.47 26.22 0.47 26.25 0.47 26.27 0.47 26.29 0.47 26.32 0.47 26.34 0.47 26.37 0.47 26.39 0.47 26.42 0.47 26.44 0.47 26.47 0.47 26.49 0.47 26.52 0.47 26.54 0.47 26.57 0.47 26.59 0.47 26.62 0.47 26.64 0.47 26.67 0.47 26.69 0.47 26.72 0.47 26.74 0.47 26.77 0.47 26.79 0.47 26.82 0.47 26.84 0.47 26.87 0.47 26.89 0.47 26.92 0.47 26.94 0.47 26.97 0.47 26.99 0.47 27.02 0.47 27.04 0.47 27.07 0.47 27.09 0.47 27.12 0.47 27.14 0.47 27.17 0.47 27.19 0.47 27.22 0.47 27.24 0.47 27.26 0.47 27.29 0.47 27.31 0.47 27.34 0.47 27.36 0.47 27.39 0.47 27.41 0.47 27.44 0.47 27.46 0.47 27.49 0.47 27.51 0.47 27.54 0.47 27.56 0.47 27.59 0.47 27.61 0.47 27.64 0.47 27.66 0.47 27.69 0.47 27.71 0.47 27.74 0.47 27.76 0.47 27.79 0.47 27.81 0.47 27.84 0.47 27.86 0.47 27.89 0.47 27.91 0.47 27.94 0.47 27.96 0.47 27.99 0.47 28.01 0.47 28.04 0.47 28.06 0.47 28.09 0.47 28.11 0.47 28.14 0.47 28.16 0.47 28.19 0.47 28.21 0.47 28.24 0.47 28.26 0.47 28.28 0.47 28.31 0.47 28.33 0.47 28.36 0.47 28.38 0.47 28.41 0.47 28.43 0.47 28.46 0.47 28.48 0.47 28.51 0.47 28.53 0.47 28.56 0.47 28.58 0.47 28.61 0.47 28.63 0.47 28.66 0.47 28.68 0.47 28.71 0.47 28.73 0.47 28.76 0.47 28.78 0.47 28.81 0.47 28.83 0.47 28.86 0.47 28.88 0.47 28.91 0.47 28.93 0.47 28.96 0.47 28.98 0.47 29.01 0.47 29.03 0.47 29.06 0.47 29.08 0.47 29.11 0.47 29.13 0.47 29.16 0.47 29.18 0.47 29.21 0.47 29.23 0.47 29.26 0.47 29.28 0.47 29.3 0.47 29.33 0.47 29.35 0.47 29.38 0.47 29.4 0.47 29.43 0.47 29.45 0.47 29.48 0.47 29.5 0.47 29.53 0.47 29.55 0.47 29.58 0.47 29.6 0.47 29.63 0.47 29.65 0.47 29.68 0.47 29.7 0.47 29.73 0.47 29.75 0.47 29.78 0.47 29.8 0.47 29.83 0.47 29.85 0.47 29.88 0.47 29.9 0.47 29.93 0.47 29.95 0.47 29.98 0.47 30 0.47 30.03 0.47 30.05 0.47 30.08 0.47 30.1 0.47 30.13 0.47 30.15 0.47 30.18 0.47 30.2 0.47 30.23 0.47 30.25 0.47 30.28 0.47 30.3 0.47 30.32 0.47 30.35 0.47 30.37 0.47 30.4 0.47 30.42 0.47 30.45 0.47 30.47 0.47 30.5 0.47 30.52 0.47 30.55 0.47 30.57 0.47 30.6 0.47 30.62 0.47 30.65 0.47 30.67 0.47 30.7 0.47 30.72 0.47 30.75 0.47 30.77 0.47 30.8 0.47 30.82 0.47 30.85 0.47 30.87 0.47 30.9 0.47 30.92 0.47 30.95 0.47 30.97 0.47 31 0.47 31.02 0.47 31.05 0.47 31.07 0.47 31.1 0.47 31.12 0.47 31.15 0.47 31.17 0.47 31.2 0.47 31.22 0.47 31.25 0.47 31.27 0.47 31.3 0.47 31.32 0.47 31.34 0.47 31.37 0.47 31.39 0.47 31.42 0.47 31.44 0.47 31.47 0.47 31.49 0.47 31.52 0.47 31.54 0.47 31.57 0.47 31.59 0.47 31.62 0.47 31.64 0.47 31.67 0.47 31.69 0.47 31.72 0.47 31.74 0.47 31.77 0.47 31.79 0.47 31.82 0.47 31.84 0.47 31.87 0.47 31.89 0.47 31.92 0.47 31.94 0.47 31.97 0.47 31.99 0.47 32.02 0.47 32.04 0.47 32.07 0.47 32.09 0.47 32.12 0.47 32.14 0.47 32.17 0.47 32.19 0.47 32.22 0.47 32.24 0.47 32.27 0.47 32.29 0.47 32.31 0.47 32.34 0.47 32.36 0.47 32.39 0.47 32.41 0.47 32.44 0.47 32.46 0.47 32.49 0.47 32.51 0.47 32.54 0.47 32.56 0.47 32.59 0.47 32.61 0.47 32.64 0.47 32.66 0.47 32.69 0.47 32.71 0.47 32.74 0.47 32.76 0.47 32.79 0.47 32.81 0.47 32.84 0.47 32.86 0.47 32.89 0.47 32.91 0.47 32.94 0.47 32.96 0.47 32.99 0.47 33.01 0.47 33.04 0.47 33.06 0.47 33.09 0.47 33.11 0.47 33.14 0.47 33.16 0.47 33.19 0.47 33.21 0.47 33.24 0.47 33.26 0.47 33.29 0.47 33.31 0.47 33.33 0.47 33.36 0.47 33.38 0.47 33.41 0.47 33.43 0.47 33.46 0.47 33.48 0.47 33.51 0.47 33.53 0.47 33.56 0.47 33.58 0.47 33.61 0.47 33.63 0.47 33.66 0.47 33.68 0.47 33.71 0.47 33.73 0.47 33.76 0.47 33.78 0.47 33.81 0.47 33.83 0.47 33.86 0.47 33.88 0.47 33.91 0.47 33.93 0.47 33.96 0.47 33.98 0.47 34.01 0.47 34.03 0.47 34.06 0.47 34.08 0.47 34.11 0.47 34.13 0.47 34.16 0.47 34.18 0.47 34.21 0.47 34.23 0.47 34.26 0.47 34.28 0.47 34.31 0.47 34.33 0.47 34.35 0.47 34.38 0.47 34.4 0.47 34.43 0.47 34.45 0.47 34.48 0.47 34.5 0.47 34.53 0.47 34.55 0.47 34.58 0.47 34.6 0.47 34.63 0.47 34.65 0.47 34.68 0.47 34.7 0.47 34.73 0.47 34.75 0.47 34.78 0.47 34.8 0.47 34.83 0.47 34.85 0.47 34.88 0.47 34.9 0.47 34.93 0.47 34.95 0.47 34.98 0.47 35 0.47 35.03 0.47 35.05 0.47 35.08 0.47 35.1 0.47 35.13 0.47 35.15 0.47 35.18 0.47 35.2 0.47 35.23 0.47 35.25 0.47 35.28 0.47 35.3 0.47 35.33 0.47 35.35 0.47 35.37 0.47 35.4 0.47 35.42 0.47 35.45 0.47 35.47 0.47 35.5 0.47 35.52 0.47 35.55 0.47 35.57 0.47 35.6 0.47 35.62 0.47 35.65 0.47 35.67 0.47 35.7 0.47 35.72 0.47 35.75 0.47 35.77 0.47 35.8 0.47 35.82 0.47 35.85 0.47 35.87 0.47 35.9 0.47 35.92 0.47 35.95 0.47 35.97 0.47 36 0.47 36.02 0.47 36.05 0.47 36.07 0.47 36.1 0.47 36.12 0.47 36.15 0.47 36.17 0.47 36.2 0.47 36.22 0.47 36.25 0.47 36.27 0.47 36.3 0.47 36.32 0.47 36.35 0.47 36.37 0.47 36.39 0.47 36.42 0.47 36.44 0.47 36.47 0.47 36.49 0.47 36.52 0.47 36.54 0.47 36.57 0.47 36.59 0.47 36.62 0.47 36.64 0.47 36.67 0.47 36.69 0.47 36.72 0.47 36.74 0.47 36.77 0.47 36.79 0.47 36.82 0.47 36.84 0.47 36.87 0.47 36.89 0.47 36.92 0.47 36.94 0.47 36.97 0.47 36.99 0.47 37.02 0.47 37.04 0.47 37.07 0.47 37.09 0.47 37.12 0.47 37.14 0.47 37.17 0.47 37.19 0.47 37.22 0.47 37.24 0.47 37.27 0.47 37.29 0.47 37.32 0.47 37.34 0.47 37.36 0.47 37.39 0.47 37.41 0.47 37.44 0.47 37.46 0.47 37.49 0.47 37.51 0.47 37.54 0.47 37.56 0.47 37.59 0.47 37.61 0.47 37.64 0.47 37.66 0.47 37.69 0.47 37.71 0.47 37.74 0.47 37.76 0.47 37.79 0.47 37.81 0.47 37.84 0.47 37.86 0.47 37.89 0.47 37.91 0.47 37.94 0.47 37.96 0.47 37.99 0.47 38.01 0.47 38.04 0.47 38.06 0.47 38.09 0.47 38.11 0.47 38.14 0.47 38.16 0.47 38.19 0.47 38.21 0.47 38.24 0.47 38.26 0.47 38.29 0.47 38.31 0.47 38.34 0.47 38.36 0.47 38.38 0.47 38.41 0.47 38.43 0.47 38.46 0.47 38.48 0.47 38.51 0.47 38.53 0.47 38.56 0.47 38.58 0.47 38.61 0.47 38.63 0.47 38.66 0.47 38.68 0.47 38.71 0.47 38.73 0.47 38.76 0.47 38.78 0.47 38.81 0.47 38.83 0.47 38.86 0.47 38.88 0.47 38.91 0.47 38.93 0.47 38.96 0.47 38.98 0.47 39.01 0.47 39.03 0.47 39.06 0.47 39.08 0.47 39.11 0.47 39.13 0.47 39.16 0.47 39.18 0.47 39.21 0.47 39.23 0.47 39.26 0.47 39.28 0.47 39.31 0.47 39.33 0.47 39.36 0.47 39.38 0.47 39.4 0.47 39.43 0.47 39.45 0.47 39.48 0.47 39.5 0.47 39.53 0.47 39.55 0.47 39.58 0.47 39.6 0.47 39.63 0.47 39.65 0.47 39.68 0.47 39.7 0.47 39.73 0.47 39.75 0.47 39.78 0.47 39.8 0.47 39.83 0.47 39.85 0.47 39.88 0.47 39.9 0.47 39.93 0.47 39.95 0.47 39.98 0.47 40 0.47 40.03 0.47 40.05 0.47 40.08 0.47 40.1 0.47 40.13 0.47 40.15 0.47 40.18 0.47 40.2 0.47 40.23 0.47 40.25 0.47 40.28 0.47 40.3 0.47 40.33 0.47 40.35 0.47 40.38 0.47 40.4 0.47 40.42 0.47 40.45 0.47 40.47 0.47 40.5 0.47 40.52 0.47 40.55 0.47 40.57 0.47 40.6 0.47 40.62 0.47 40.65 0.47 40.67 0.47 40.7 0.47 40.72 0.47 40.75 0.47 40.77 0.47 40.8 0.47 40.82 0.47 40.85 0.47 40.87 0.47 40.9 0.47 40.92 0.47 40.95 0.47 40.97 0.47 41 0.47 41.02 0.47 41.05 0.47 41.07 0.47 41.1 0.47 41.12 0.47 41.15 0.47 41.17 0.47 41.2 0.47 41.22 0.47 41.25 0.47 41.27 0.47 41.3 0.47 41.32 0.47 41.35 0.47 41.37 0.47 41.39 0.47 41.42 0.47 41.44 0.47 41.47 0.47 41.49 0.47 41.52 0.47 41.54 0.47 41.57 0.47 41.59 0.47 41.62 0.47 41.64 0.47 41.67 0.47 41.69 0.47 41.72 0.47 41.74 0.47 41.77 0.47 41.79 0.47 41.82 0.47 41.84 0.47 41.87 0.47 41.89 0.47 41.92 0.47 41.94 0.47 41.97 0.47 41.99 0.47 42.02 0.47 42.04 0.47 42.07 0.47 42.09 0.47 42.12 0.47 42.14 0.47 42.17 0.47 42.19 0.47 42.22 0.47 42.24 0.47 42.27 0.47 42.29 0.47 42.32 0.47 42.34 0.47 42.37 0.47 42.39 0.47 42.41 0.47 42.44 0.47 42.46 0.47 42.49 0.47 42.51 0.47 42.54 0.47 42.56 0.47 42.59 0.47 42.61 0.47 42.64 0.47 42.66 0.47 42.69 0.47 42.71 0.47 42.74 0.47 42.76 0.47 42.79 0.47 42.81 0.47 42.84 0.47 42.86 0.47 42.89 0.47 42.91 0.47 42.94 0.47 42.96 0.47 42.99 0.47 43.01 0.47 43.04 0.47 43.06 0.47 43.09 0.47 43.11 0.47 43.14 0.47 43.16 0.47 43.19 0.47 43.21 0.47 43.24 0.47 43.26 0.47 43.29 0.47 43.31 0.47 43.34 0.47 43.36 0.47 43.39 0.47 43.41 0.47 43.43 0.47 43.46 0.47 43.48 0.47 43.51 0.47 43.53 0.47 43.56 0.47 43.58 0.47 43.61 0.47 43.63 0.47 43.66 0.47 43.68 0.47 43.71 0.47 43.73 0.47 43.76 0.47 43.78 0.47 43.81 0.47 43.83 0.47 43.86 0.47 43.88 0.47 43.91 0.47 43.93 0.47 43.96 0.47 43.98 0.47 44.01 0.47 44.03 0.47 44.06 0.47 44.08 0.47 44.11 0.47 44.13 0.47 44.16 0.47 44.18 0.47 44.21 0.47 44.23 0.47 44.26 0.47 44.28 0.47 44.31 0.47 44.33 0.47 44.36 0.47 44.38 0.47 44.41 0.47 44.43 0.47 44.45 0.47 44.48 0.47 44.5 0.47 44.53 0.47 44.55 0.47 44.58 0.47 44.6 0.47 44.63 0.47 44.65 0.47 44.68 0.47 44.7 0.47 44.73 0.47 44.75 0.47 44.78 0.47 44.8 0.47 44.83 0.47 44.85 0.47 44.88 0.47 44.9 0.47 44.93 0.47 44.95 0.47 44.98 0.47 45 0.47 45.03 0.47 45.05 0.47 45.08 0.47 45.1 0.47 45.13 0.47 45.15 0.47 45.18 0.47 45.2 0.47 45.23 0.47 45.25 0.47 45.28 0.47 45.3 0.47 45.33 0.47 45.35 0.47 45.38 0.47 45.4 0.47" class="primitive"/>
          </g>
        </g>
      </g>
      <g opacity="0" class="guide zoomslider" stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-749">
        <g fill="#EAEAEA" stroke-width="0.3" stroke-opacity="0" stroke="#6A6A6A" id="img-8311e7c5-750">
          <g transform="translate(118.14,10)" id="img-8311e7c5-751">
            <path d="M-2,-2 L 2 -2 2 2 -2 2 z" class="primitive"/>
          </g>
          <g class="button_logo" fill="#6A6A6A" id="img-8311e7c5-752">
            <g transform="translate(118.14,10)" id="img-8311e7c5-753">
              <path d="M-1.2,-0.4 L -0.4 -0.4 -0.4 -1.2 0.4 -1.2 0.4 -0.4 1.2 -0.4 1.2 0.4 0.4 0.4 0.4 1.2 -0.4 1.2 -0.4 0.4 -1.2 0.4 z" class="primitive"/>
            </g>
          </g>
        </g>
        <g fill="#EAEAEA" id="img-8311e7c5-754">
          <g transform="translate(106.14,10)" id="img-8311e7c5-755">
            <path d="M-9.5,-2 L 9.5 -2 9.5 2 -9.5 2 z" class="primitive"/>
          </g>
        </g>
        <g class="zoomslider_thumb" fill="#6A6A6A" id="img-8311e7c5-756">
          <g transform="translate(106.14,10)" id="img-8311e7c5-757">
            <path d="M-1,-2 L 1 -2 1 2 -1 2 z" class="primitive"/>
          </g>
        </g>
        <g fill="#EAEAEA" stroke-width="0.3" stroke-opacity="0" stroke="#6A6A6A" id="img-8311e7c5-758">
          <g transform="translate(94.14,10)" id="img-8311e7c5-759">
            <path d="M-2,-2 L 2 -2 2 2 -2 2 z" class="primitive"/>
          </g>
          <g class="button_logo" fill="#6A6A6A" id="img-8311e7c5-760">
            <g transform="translate(94.14,10)" id="img-8311e7c5-761">
              <path d="M-1.2,-0.4 L 1.2 -0.4 1.2 0.4 -1.2 0.4 z" class="primitive"/>
            </g>
          </g>
        </g>
      </g>
    </g>
  </g>
  <g class="guide ylabels" font-size="2.82" font-family="'PT Sans Caption','Helvetica Neue','Helvetica',sans-serif" fill="#6C606B" id="img-8311e7c5-762">
    <g transform="translate(18.63,57.67)" id="img-8311e7c5-763" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-200</text>
      </g>
    </g>
    <g transform="translate(18.63,50.43)" id="img-8311e7c5-764" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-150</text>
      </g>
    </g>
    <g transform="translate(18.63,43.19)" id="img-8311e7c5-765" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-100</text>
      </g>
    </g>
    <g transform="translate(18.63,35.95)" id="img-8311e7c5-766" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-50</text>
      </g>
    </g>
    <g transform="translate(18.63,28.72)" id="img-8311e7c5-767" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">0</text>
      </g>
    </g>
    <g transform="translate(18.63,21.48)" id="img-8311e7c5-768" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">50</text>
      </g>
    </g>
    <g transform="translate(18.63,14.24)" id="img-8311e7c5-769" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">100</text>
      </g>
    </g>
    <g transform="translate(18.63,7)" id="img-8311e7c5-770" gadfly:scale="1.0" visibility="visible">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">150</text>
      </g>
    </g>
    <g transform="translate(18.63,-0.24)" id="img-8311e7c5-771" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">200</text>
      </g>
    </g>
    <g transform="translate(18.63,-7.48)" id="img-8311e7c5-772" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">250</text>
      </g>
    </g>
    <g transform="translate(18.63,-14.71)" id="img-8311e7c5-773" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">300</text>
      </g>
    </g>
    <g transform="translate(18.63,-21.95)" id="img-8311e7c5-774" gadfly:scale="1.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">350</text>
      </g>
    </g>
    <g transform="translate(18.63,50.43)" id="img-8311e7c5-775" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-150</text>
      </g>
    </g>
    <g transform="translate(18.63,49.71)" id="img-8311e7c5-776" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-145</text>
      </g>
    </g>
    <g transform="translate(18.63,48.98)" id="img-8311e7c5-777" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-140</text>
      </g>
    </g>
    <g transform="translate(18.63,48.26)" id="img-8311e7c5-778" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-135</text>
      </g>
    </g>
    <g transform="translate(18.63,47.53)" id="img-8311e7c5-779" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-130</text>
      </g>
    </g>
    <g transform="translate(18.63,46.81)" id="img-8311e7c5-780" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-125</text>
      </g>
    </g>
    <g transform="translate(18.63,46.09)" id="img-8311e7c5-781" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-120</text>
      </g>
    </g>
    <g transform="translate(18.63,45.36)" id="img-8311e7c5-782" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-115</text>
      </g>
    </g>
    <g transform="translate(18.63,44.64)" id="img-8311e7c5-783" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-110</text>
      </g>
    </g>
    <g transform="translate(18.63,43.92)" id="img-8311e7c5-784" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-105</text>
      </g>
    </g>
    <g transform="translate(18.63,43.19)" id="img-8311e7c5-785" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-100</text>
      </g>
    </g>
    <g transform="translate(18.63,42.47)" id="img-8311e7c5-786" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-95</text>
      </g>
    </g>
    <g transform="translate(18.63,41.74)" id="img-8311e7c5-787" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-90</text>
      </g>
    </g>
    <g transform="translate(18.63,41.02)" id="img-8311e7c5-788" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-85</text>
      </g>
    </g>
    <g transform="translate(18.63,40.3)" id="img-8311e7c5-789" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-80</text>
      </g>
    </g>
    <g transform="translate(18.63,39.57)" id="img-8311e7c5-790" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-75</text>
      </g>
    </g>
    <g transform="translate(18.63,38.85)" id="img-8311e7c5-791" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-70</text>
      </g>
    </g>
    <g transform="translate(18.63,38.12)" id="img-8311e7c5-792" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-65</text>
      </g>
    </g>
    <g transform="translate(18.63,37.4)" id="img-8311e7c5-793" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-60</text>
      </g>
    </g>
    <g transform="translate(18.63,36.68)" id="img-8311e7c5-794" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-55</text>
      </g>
    </g>
    <g transform="translate(18.63,35.95)" id="img-8311e7c5-795" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-50</text>
      </g>
    </g>
    <g transform="translate(18.63,35.23)" id="img-8311e7c5-796" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-45</text>
      </g>
    </g>
    <g transform="translate(18.63,34.51)" id="img-8311e7c5-797" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-40</text>
      </g>
    </g>
    <g transform="translate(18.63,33.78)" id="img-8311e7c5-798" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-35</text>
      </g>
    </g>
    <g transform="translate(18.63,33.06)" id="img-8311e7c5-799" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-30</text>
      </g>
    </g>
    <g transform="translate(18.63,32.33)" id="img-8311e7c5-800" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-25</text>
      </g>
    </g>
    <g transform="translate(18.63,31.61)" id="img-8311e7c5-801" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-20</text>
      </g>
    </g>
    <g transform="translate(18.63,30.89)" id="img-8311e7c5-802" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-15</text>
      </g>
    </g>
    <g transform="translate(18.63,30.16)" id="img-8311e7c5-803" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-10</text>
      </g>
    </g>
    <g transform="translate(18.63,29.44)" id="img-8311e7c5-804" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-5</text>
      </g>
    </g>
    <g transform="translate(18.63,28.72)" id="img-8311e7c5-805" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">0</text>
      </g>
    </g>
    <g transform="translate(18.63,27.99)" id="img-8311e7c5-806" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">5</text>
      </g>
    </g>
    <g transform="translate(18.63,27.27)" id="img-8311e7c5-807" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">10</text>
      </g>
    </g>
    <g transform="translate(18.63,26.54)" id="img-8311e7c5-808" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">15</text>
      </g>
    </g>
    <g transform="translate(18.63,25.82)" id="img-8311e7c5-809" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">20</text>
      </g>
    </g>
    <g transform="translate(18.63,25.1)" id="img-8311e7c5-810" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">25</text>
      </g>
    </g>
    <g transform="translate(18.63,24.37)" id="img-8311e7c5-811" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">30</text>
      </g>
    </g>
    <g transform="translate(18.63,23.65)" id="img-8311e7c5-812" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">35</text>
      </g>
    </g>
    <g transform="translate(18.63,22.92)" id="img-8311e7c5-813" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">40</text>
      </g>
    </g>
    <g transform="translate(18.63,22.2)" id="img-8311e7c5-814" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">45</text>
      </g>
    </g>
    <g transform="translate(18.63,21.48)" id="img-8311e7c5-815" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">50</text>
      </g>
    </g>
    <g transform="translate(18.63,20.75)" id="img-8311e7c5-816" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">55</text>
      </g>
    </g>
    <g transform="translate(18.63,20.03)" id="img-8311e7c5-817" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">60</text>
      </g>
    </g>
    <g transform="translate(18.63,19.31)" id="img-8311e7c5-818" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">65</text>
      </g>
    </g>
    <g transform="translate(18.63,18.58)" id="img-8311e7c5-819" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">70</text>
      </g>
    </g>
    <g transform="translate(18.63,17.86)" id="img-8311e7c5-820" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">75</text>
      </g>
    </g>
    <g transform="translate(18.63,17.13)" id="img-8311e7c5-821" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">80</text>
      </g>
    </g>
    <g transform="translate(18.63,16.41)" id="img-8311e7c5-822" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">85</text>
      </g>
    </g>
    <g transform="translate(18.63,15.69)" id="img-8311e7c5-823" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">90</text>
      </g>
    </g>
    <g transform="translate(18.63,14.96)" id="img-8311e7c5-824" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">95</text>
      </g>
    </g>
    <g transform="translate(18.63,14.24)" id="img-8311e7c5-825" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">100</text>
      </g>
    </g>
    <g transform="translate(18.63,13.51)" id="img-8311e7c5-826" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">105</text>
      </g>
    </g>
    <g transform="translate(18.63,12.79)" id="img-8311e7c5-827" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">110</text>
      </g>
    </g>
    <g transform="translate(18.63,12.07)" id="img-8311e7c5-828" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">115</text>
      </g>
    </g>
    <g transform="translate(18.63,11.34)" id="img-8311e7c5-829" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">120</text>
      </g>
    </g>
    <g transform="translate(18.63,10.62)" id="img-8311e7c5-830" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">125</text>
      </g>
    </g>
    <g transform="translate(18.63,9.9)" id="img-8311e7c5-831" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">130</text>
      </g>
    </g>
    <g transform="translate(18.63,9.17)" id="img-8311e7c5-832" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">135</text>
      </g>
    </g>
    <g transform="translate(18.63,8.45)" id="img-8311e7c5-833" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">140</text>
      </g>
    </g>
    <g transform="translate(18.63,7.72)" id="img-8311e7c5-834" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">145</text>
      </g>
    </g>
    <g transform="translate(18.63,7)" id="img-8311e7c5-835" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">150</text>
      </g>
    </g>
    <g transform="translate(18.63,6.28)" id="img-8311e7c5-836" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">155</text>
      </g>
    </g>
    <g transform="translate(18.63,5.55)" id="img-8311e7c5-837" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">160</text>
      </g>
    </g>
    <g transform="translate(18.63,4.83)" id="img-8311e7c5-838" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">165</text>
      </g>
    </g>
    <g transform="translate(18.63,4.1)" id="img-8311e7c5-839" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">170</text>
      </g>
    </g>
    <g transform="translate(18.63,3.38)" id="img-8311e7c5-840" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">175</text>
      </g>
    </g>
    <g transform="translate(18.63,2.66)" id="img-8311e7c5-841" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">180</text>
      </g>
    </g>
    <g transform="translate(18.63,1.93)" id="img-8311e7c5-842" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">185</text>
      </g>
    </g>
    <g transform="translate(18.63,1.21)" id="img-8311e7c5-843" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">190</text>
      </g>
    </g>
    <g transform="translate(18.63,0.49)" id="img-8311e7c5-844" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">195</text>
      </g>
    </g>
    <g transform="translate(18.63,-0.24)" id="img-8311e7c5-845" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">200</text>
      </g>
    </g>
    <g transform="translate(18.63,-0.96)" id="img-8311e7c5-846" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">205</text>
      </g>
    </g>
    <g transform="translate(18.63,-1.69)" id="img-8311e7c5-847" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">210</text>
      </g>
    </g>
    <g transform="translate(18.63,-2.41)" id="img-8311e7c5-848" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">215</text>
      </g>
    </g>
    <g transform="translate(18.63,-3.13)" id="img-8311e7c5-849" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">220</text>
      </g>
    </g>
    <g transform="translate(18.63,-3.86)" id="img-8311e7c5-850" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">225</text>
      </g>
    </g>
    <g transform="translate(18.63,-4.58)" id="img-8311e7c5-851" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">230</text>
      </g>
    </g>
    <g transform="translate(18.63,-5.31)" id="img-8311e7c5-852" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">235</text>
      </g>
    </g>
    <g transform="translate(18.63,-6.03)" id="img-8311e7c5-853" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">240</text>
      </g>
    </g>
    <g transform="translate(18.63,-6.75)" id="img-8311e7c5-854" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">245</text>
      </g>
    </g>
    <g transform="translate(18.63,-7.48)" id="img-8311e7c5-855" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">250</text>
      </g>
    </g>
    <g transform="translate(18.63,-8.2)" id="img-8311e7c5-856" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">255</text>
      </g>
    </g>
    <g transform="translate(18.63,-8.92)" id="img-8311e7c5-857" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">260</text>
      </g>
    </g>
    <g transform="translate(18.63,-9.65)" id="img-8311e7c5-858" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">265</text>
      </g>
    </g>
    <g transform="translate(18.63,-10.37)" id="img-8311e7c5-859" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">270</text>
      </g>
    </g>
    <g transform="translate(18.63,-11.1)" id="img-8311e7c5-860" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">275</text>
      </g>
    </g>
    <g transform="translate(18.63,-11.82)" id="img-8311e7c5-861" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">280</text>
      </g>
    </g>
    <g transform="translate(18.63,-12.54)" id="img-8311e7c5-862" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">285</text>
      </g>
    </g>
    <g transform="translate(18.63,-13.27)" id="img-8311e7c5-863" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">290</text>
      </g>
    </g>
    <g transform="translate(18.63,-13.99)" id="img-8311e7c5-864" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">295</text>
      </g>
    </g>
    <g transform="translate(18.63,-14.71)" id="img-8311e7c5-865" gadfly:scale="10.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">300</text>
      </g>
    </g>
    <g transform="translate(18.63,57.67)" id="img-8311e7c5-866" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-200</text>
      </g>
    </g>
    <g transform="translate(18.63,28.72)" id="img-8311e7c5-867" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">0</text>
      </g>
    </g>
    <g transform="translate(18.63,-0.24)" id="img-8311e7c5-868" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">200</text>
      </g>
    </g>
    <g transform="translate(18.63,-29.19)" id="img-8311e7c5-869" gadfly:scale="0.5" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">400</text>
      </g>
    </g>
    <g transform="translate(18.63,50.43)" id="img-8311e7c5-870" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-150</text>
      </g>
    </g>
    <g transform="translate(18.63,48.98)" id="img-8311e7c5-871" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-140</text>
      </g>
    </g>
    <g transform="translate(18.63,47.53)" id="img-8311e7c5-872" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-130</text>
      </g>
    </g>
    <g transform="translate(18.63,46.09)" id="img-8311e7c5-873" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-120</text>
      </g>
    </g>
    <g transform="translate(18.63,44.64)" id="img-8311e7c5-874" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-110</text>
      </g>
    </g>
    <g transform="translate(18.63,43.19)" id="img-8311e7c5-875" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-100</text>
      </g>
    </g>
    <g transform="translate(18.63,41.74)" id="img-8311e7c5-876" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-90</text>
      </g>
    </g>
    <g transform="translate(18.63,40.3)" id="img-8311e7c5-877" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-80</text>
      </g>
    </g>
    <g transform="translate(18.63,38.85)" id="img-8311e7c5-878" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-70</text>
      </g>
    </g>
    <g transform="translate(18.63,37.4)" id="img-8311e7c5-879" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-60</text>
      </g>
    </g>
    <g transform="translate(18.63,35.95)" id="img-8311e7c5-880" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-50</text>
      </g>
    </g>
    <g transform="translate(18.63,34.51)" id="img-8311e7c5-881" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-40</text>
      </g>
    </g>
    <g transform="translate(18.63,33.06)" id="img-8311e7c5-882" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-30</text>
      </g>
    </g>
    <g transform="translate(18.63,31.61)" id="img-8311e7c5-883" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-20</text>
      </g>
    </g>
    <g transform="translate(18.63,30.16)" id="img-8311e7c5-884" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">-10</text>
      </g>
    </g>
    <g transform="translate(18.63,28.72)" id="img-8311e7c5-885" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">0</text>
      </g>
    </g>
    <g transform="translate(18.63,27.27)" id="img-8311e7c5-886" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">10</text>
      </g>
    </g>
    <g transform="translate(18.63,25.82)" id="img-8311e7c5-887" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">20</text>
      </g>
    </g>
    <g transform="translate(18.63,24.37)" id="img-8311e7c5-888" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">30</text>
      </g>
    </g>
    <g transform="translate(18.63,22.92)" id="img-8311e7c5-889" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">40</text>
      </g>
    </g>
    <g transform="translate(18.63,21.48)" id="img-8311e7c5-890" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">50</text>
      </g>
    </g>
    <g transform="translate(18.63,20.03)" id="img-8311e7c5-891" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">60</text>
      </g>
    </g>
    <g transform="translate(18.63,18.58)" id="img-8311e7c5-892" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">70</text>
      </g>
    </g>
    <g transform="translate(18.63,17.13)" id="img-8311e7c5-893" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">80</text>
      </g>
    </g>
    <g transform="translate(18.63,15.69)" id="img-8311e7c5-894" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">90</text>
      </g>
    </g>
    <g transform="translate(18.63,14.24)" id="img-8311e7c5-895" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">100</text>
      </g>
    </g>
    <g transform="translate(18.63,12.79)" id="img-8311e7c5-896" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">110</text>
      </g>
    </g>
    <g transform="translate(18.63,11.34)" id="img-8311e7c5-897" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">120</text>
      </g>
    </g>
    <g transform="translate(18.63,9.9)" id="img-8311e7c5-898" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">130</text>
      </g>
    </g>
    <g transform="translate(18.63,8.45)" id="img-8311e7c5-899" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">140</text>
      </g>
    </g>
    <g transform="translate(18.63,7)" id="img-8311e7c5-900" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">150</text>
      </g>
    </g>
    <g transform="translate(18.63,5.55)" id="img-8311e7c5-901" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">160</text>
      </g>
    </g>
    <g transform="translate(18.63,4.1)" id="img-8311e7c5-902" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">170</text>
      </g>
    </g>
    <g transform="translate(18.63,2.66)" id="img-8311e7c5-903" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">180</text>
      </g>
    </g>
    <g transform="translate(18.63,1.21)" id="img-8311e7c5-904" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">190</text>
      </g>
    </g>
    <g transform="translate(18.63,-0.24)" id="img-8311e7c5-905" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">200</text>
      </g>
    </g>
    <g transform="translate(18.63,-1.69)" id="img-8311e7c5-906" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">210</text>
      </g>
    </g>
    <g transform="translate(18.63,-3.13)" id="img-8311e7c5-907" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">220</text>
      </g>
    </g>
    <g transform="translate(18.63,-4.58)" id="img-8311e7c5-908" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">230</text>
      </g>
    </g>
    <g transform="translate(18.63,-6.03)" id="img-8311e7c5-909" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">240</text>
      </g>
    </g>
    <g transform="translate(18.63,-7.48)" id="img-8311e7c5-910" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">250</text>
      </g>
    </g>
    <g transform="translate(18.63,-8.92)" id="img-8311e7c5-911" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">260</text>
      </g>
    </g>
    <g transform="translate(18.63,-10.37)" id="img-8311e7c5-912" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">270</text>
      </g>
    </g>
    <g transform="translate(18.63,-11.82)" id="img-8311e7c5-913" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">280</text>
      </g>
    </g>
    <g transform="translate(18.63,-13.27)" id="img-8311e7c5-914" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">290</text>
      </g>
    </g>
    <g transform="translate(18.63,-14.71)" id="img-8311e7c5-915" gadfly:scale="5.0" visibility="hidden">
      <g class="primitive">
        <text text-anchor="end" dy="0.35em">300</text>
      </g>
    </g>
  </g>
  <g font-size="3.88" font-family="'PT Sans','Helvetica Neue','Helvetica',sans-serif" fill="#564A55" stroke="#000000" stroke-opacity="0.000" id="img-8311e7c5-916">
    <g transform="translate(8.81,15.86)" id="img-8311e7c5-917">
      <g class="primitive">
        <text text-anchor="middle" dy="0.35em" transform="rotate(-90,0, 2)">value</text>
      </g>
    </g>
  </g>
</g>
<defs>
  <clipPath id="img-8311e7c5-15">
    <path d="M26.72,55 L 123.14 55 123.14 80.72 26.72 80.72" />
  </clipPath>
  <clipPath id="img-8311e7c5-479">
    <path d="M19.63,5 L 123.14 5 123.14 30.72 19.63 30.72" />
  </clipPath>
</defs>
<script> <![CDATA[
(function(N){var k=/[\.\/]/,L=/\s*,\s*/,C=function(a,d){return a-d},a,v,y={n:{}},M=function(){for(var a=0,d=this.length;a<d;a++)if("undefined"!=typeof this[a])return this[a]},A=function(){for(var a=this.length;--a;)if("undefined"!=typeof this[a])return this[a]},w=function(k,d){k=String(k);var f=v,n=Array.prototype.slice.call(arguments,2),u=w.listeners(k),p=0,b,q=[],e={},l=[],r=a;l.firstDefined=M;l.lastDefined=A;a=k;for(var s=v=0,x=u.length;s<x;s++)"zIndex"in u[s]&&(q.push(u[s].zIndex),0>u[s].zIndex&&
(e[u[s].zIndex]=u[s]));for(q.sort(C);0>q[p];)if(b=e[q[p++] ],l.push(b.apply(d,n)),v)return v=f,l;for(s=0;s<x;s++)if(b=u[s],"zIndex"in b)if(b.zIndex==q[p]){l.push(b.apply(d,n));if(v)break;do if(p++,(b=e[q[p] ])&&l.push(b.apply(d,n)),v)break;while(b)}else e[b.zIndex]=b;else if(l.push(b.apply(d,n)),v)break;v=f;a=r;return l};w._events=y;w.listeners=function(a){a=a.split(k);var d=y,f,n,u,p,b,q,e,l=[d],r=[];u=0;for(p=a.length;u<p;u++){e=[];b=0;for(q=l.length;b<q;b++)for(d=l[b].n,f=[d[a[u] ],d["*"] ],n=2;n--;)if(d=
f[n])e.push(d),r=r.concat(d.f||[]);l=e}return r};w.on=function(a,d){a=String(a);if("function"!=typeof d)return function(){};for(var f=a.split(L),n=0,u=f.length;n<u;n++)(function(a){a=a.split(k);for(var b=y,f,e=0,l=a.length;e<l;e++)b=b.n,b=b.hasOwnProperty(a[e])&&b[a[e] ]||(b[a[e] ]={n:{}});b.f=b.f||[];e=0;for(l=b.f.length;e<l;e++)if(b.f[e]==d){f=!0;break}!f&&b.f.push(d)})(f[n]);return function(a){+a==+a&&(d.zIndex=+a)}};w.f=function(a){var d=[].slice.call(arguments,1);return function(){w.apply(null,
[a,null].concat(d).concat([].slice.call(arguments,0)))}};w.stop=function(){v=1};w.nt=function(k){return k?(new RegExp("(?:\\.|\\/|^)"+k+"(?:\\.|\\/|$)")).test(a):a};w.nts=function(){return a.split(k)};w.off=w.unbind=function(a,d){if(a){var f=a.split(L);if(1<f.length)for(var n=0,u=f.length;n<u;n++)w.off(f[n],d);else{for(var f=a.split(k),p,b,q,e,l=[y],n=0,u=f.length;n<u;n++)for(e=0;e<l.length;e+=q.length-2){q=[e,1];p=l[e].n;if("*"!=f[n])p[f[n] ]&&q.push(p[f[n] ]);else for(b in p)p.hasOwnProperty(b)&&
q.push(p[b]);l.splice.apply(l,q)}n=0;for(u=l.length;n<u;n++)for(p=l[n];p.n;){if(d){if(p.f){e=0;for(f=p.f.length;e<f;e++)if(p.f[e]==d){p.f.splice(e,1);break}!p.f.length&&delete p.f}for(b in p.n)if(p.n.hasOwnProperty(b)&&p.n[b].f){q=p.n[b].f;e=0;for(f=q.length;e<f;e++)if(q[e]==d){q.splice(e,1);break}!q.length&&delete p.n[b].f}}else for(b in delete p.f,p.n)p.n.hasOwnProperty(b)&&p.n[b].f&&delete p.n[b].f;p=p.n}}}else w._events=y={n:{}}};w.once=function(a,d){var f=function(){w.unbind(a,f);return d.apply(this,
arguments)};return w.on(a,f)};w.version="0.4.2";w.toString=function(){return"You are running Eve 0.4.2"};"undefined"!=typeof module&&module.exports?module.exports=w:"function"===typeof define&&define.amd?define("eve",[],function(){return w}):N.eve=w})(this);
(function(N,k){"function"===typeof define&&define.amd?define("Snap.svg",["eve"],function(L){return k(N,L)}):k(N,N.eve)})(this,function(N,k){var L=function(a){var k={},y=N.requestAnimationFrame||N.webkitRequestAnimationFrame||N.mozRequestAnimationFrame||N.oRequestAnimationFrame||N.msRequestAnimationFrame||function(a){setTimeout(a,16)},M=Array.isArray||function(a){return a instanceof Array||"[object Array]"==Object.prototype.toString.call(a)},A=0,w="M"+(+new Date).toString(36),z=function(a){if(null==
a)return this.s;var b=this.s-a;this.b+=this.dur*b;this.B+=this.dur*b;this.s=a},d=function(a){if(null==a)return this.spd;this.spd=a},f=function(a){if(null==a)return this.dur;this.s=this.s*a/this.dur;this.dur=a},n=function(){delete k[this.id];this.update();a("mina.stop."+this.id,this)},u=function(){this.pdif||(delete k[this.id],this.update(),this.pdif=this.get()-this.b)},p=function(){this.pdif&&(this.b=this.get()-this.pdif,delete this.pdif,k[this.id]=this)},b=function(){var a;if(M(this.start)){a=[];
for(var b=0,e=this.start.length;b<e;b++)a[b]=+this.start[b]+(this.end[b]-this.start[b])*this.easing(this.s)}else a=+this.start+(this.end-this.start)*this.easing(this.s);this.set(a)},q=function(){var l=0,b;for(b in k)if(k.hasOwnProperty(b)){var e=k[b],f=e.get();l++;e.s=(f-e.b)/(e.dur/e.spd);1<=e.s&&(delete k[b],e.s=1,l--,function(b){setTimeout(function(){a("mina.finish."+b.id,b)})}(e));e.update()}l&&y(q)},e=function(a,r,s,x,G,h,J){a={id:w+(A++).toString(36),start:a,end:r,b:s,s:0,dur:x-s,spd:1,get:G,
set:h,easing:J||e.linear,status:z,speed:d,duration:f,stop:n,pause:u,resume:p,update:b};k[a.id]=a;r=0;for(var K in k)if(k.hasOwnProperty(K)&&(r++,2==r))break;1==r&&y(q);return a};e.time=Date.now||function(){return+new Date};e.getById=function(a){return k[a]||null};e.linear=function(a){return a};e.easeout=function(a){return Math.pow(a,1.7)};e.easein=function(a){return Math.pow(a,0.48)};e.easeinout=function(a){if(1==a)return 1;if(0==a)return 0;var b=0.48-a/1.04,e=Math.sqrt(0.1734+b*b);a=e-b;a=Math.pow(Math.abs(a),
1/3)*(0>a?-1:1);b=-e-b;b=Math.pow(Math.abs(b),1/3)*(0>b?-1:1);a=a+b+0.5;return 3*(1-a)*a*a+a*a*a};e.backin=function(a){return 1==a?1:a*a*(2.70158*a-1.70158)};e.backout=function(a){if(0==a)return 0;a-=1;return a*a*(2.70158*a+1.70158)+1};e.elastic=function(a){return a==!!a?a:Math.pow(2,-10*a)*Math.sin(2*(a-0.075)*Math.PI/0.3)+1};e.bounce=function(a){a<1/2.75?a*=7.5625*a:a<2/2.75?(a-=1.5/2.75,a=7.5625*a*a+0.75):a<2.5/2.75?(a-=2.25/2.75,a=7.5625*a*a+0.9375):(a-=2.625/2.75,a=7.5625*a*a+0.984375);return a};
return N.mina=e}("undefined"==typeof k?function(){}:k),C=function(){function a(c,t){if(c){if(c.tagName)return x(c);if(y(c,"array")&&a.set)return a.set.apply(a,c);if(c instanceof e)return c;if(null==t)return c=G.doc.querySelector(c),x(c)}return new s(null==c?"100%":c,null==t?"100%":t)}function v(c,a){if(a){"#text"==c&&(c=G.doc.createTextNode(a.text||""));"string"==typeof c&&(c=v(c));if("string"==typeof a)return"xlink:"==a.substring(0,6)?c.getAttributeNS(m,a.substring(6)):"xml:"==a.substring(0,4)?c.getAttributeNS(la,
a.substring(4)):c.getAttribute(a);for(var da in a)if(a[h](da)){var b=J(a[da]);b?"xlink:"==da.substring(0,6)?c.setAttributeNS(m,da.substring(6),b):"xml:"==da.substring(0,4)?c.setAttributeNS(la,da.substring(4),b):c.setAttribute(da,b):c.removeAttribute(da)}}else c=G.doc.createElementNS(la,c);return c}function y(c,a){a=J.prototype.toLowerCase.call(a);return"finite"==a?isFinite(c):"array"==a&&(c instanceof Array||Array.isArray&&Array.isArray(c))?!0:"null"==a&&null===c||a==typeof c&&null!==c||"object"==
a&&c===Object(c)||$.call(c).slice(8,-1).toLowerCase()==a}function M(c){if("function"==typeof c||Object(c)!==c)return c;var a=new c.constructor,b;for(b in c)c[h](b)&&(a[b]=M(c[b]));return a}function A(c,a,b){function m(){var e=Array.prototype.slice.call(arguments,0),f=e.join("\u2400"),d=m.cache=m.cache||{},l=m.count=m.count||[];if(d[h](f)){a:for(var e=l,l=f,B=0,H=e.length;B<H;B++)if(e[B]===l){e.push(e.splice(B,1)[0]);break a}return b?b(d[f]):d[f]}1E3<=l.length&&delete d[l.shift()];l.push(f);d[f]=c.apply(a,
e);return b?b(d[f]):d[f]}return m}function w(c,a,b,m,e,f){return null==e?(c-=b,a-=m,c||a?(180*I.atan2(-a,-c)/C+540)%360:0):w(c,a,e,f)-w(b,m,e,f)}function z(c){return c%360*C/180}function d(c){var a=[];c=c.replace(/(?:^|\s)(\w+)\(([^)]+)\)/g,function(c,b,m){m=m.split(/\s*,\s*|\s+/);"rotate"==b&&1==m.length&&m.push(0,0);"scale"==b&&(2<m.length?m=m.slice(0,2):2==m.length&&m.push(0,0),1==m.length&&m.push(m[0],0,0));"skewX"==b?a.push(["m",1,0,I.tan(z(m[0])),1,0,0]):"skewY"==b?a.push(["m",1,I.tan(z(m[0])),
0,1,0,0]):a.push([b.charAt(0)].concat(m));return c});return a}function f(c,t){var b=O(c),m=new a.Matrix;if(b)for(var e=0,f=b.length;e<f;e++){var h=b[e],d=h.length,B=J(h[0]).toLowerCase(),H=h[0]!=B,l=H?m.invert():0,E;"t"==B&&2==d?m.translate(h[1],0):"t"==B&&3==d?H?(d=l.x(0,0),B=l.y(0,0),H=l.x(h[1],h[2]),l=l.y(h[1],h[2]),m.translate(H-d,l-B)):m.translate(h[1],h[2]):"r"==B?2==d?(E=E||t,m.rotate(h[1],E.x+E.width/2,E.y+E.height/2)):4==d&&(H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.rotate(h[1],H,l)):m.rotate(h[1],
h[2],h[3])):"s"==B?2==d||3==d?(E=E||t,m.scale(h[1],h[d-1],E.x+E.width/2,E.y+E.height/2)):4==d?H?(H=l.x(h[2],h[3]),l=l.y(h[2],h[3]),m.scale(h[1],h[1],H,l)):m.scale(h[1],h[1],h[2],h[3]):5==d&&(H?(H=l.x(h[3],h[4]),l=l.y(h[3],h[4]),m.scale(h[1],h[2],H,l)):m.scale(h[1],h[2],h[3],h[4])):"m"==B&&7==d&&m.add(h[1],h[2],h[3],h[4],h[5],h[6])}return m}function n(c,t){if(null==t){var m=!0;t="linearGradient"==c.type||"radialGradient"==c.type?c.node.getAttribute("gradientTransform"):"pattern"==c.type?c.node.getAttribute("patternTransform"):
c.node.getAttribute("transform");if(!t)return new a.Matrix;t=d(t)}else t=a._.rgTransform.test(t)?J(t).replace(/\.{3}|\u2026/g,c._.transform||aa):d(t),y(t,"array")&&(t=a.path?a.path.toString.call(t):J(t)),c._.transform=t;var b=f(t,c.getBBox(1));if(m)return b;c.matrix=b}function u(c){c=c.node.ownerSVGElement&&x(c.node.ownerSVGElement)||c.node.parentNode&&x(c.node.parentNode)||a.select("svg")||a(0,0);var t=c.select("defs"),t=null==t?!1:t.node;t||(t=r("defs",c.node).node);return t}function p(c){return c.node.ownerSVGElement&&
x(c.node.ownerSVGElement)||a.select("svg")}function b(c,a,m){function b(c){if(null==c)return aa;if(c==+c)return c;v(B,{width:c});try{return B.getBBox().width}catch(a){return 0}}function h(c){if(null==c)return aa;if(c==+c)return c;v(B,{height:c});try{return B.getBBox().height}catch(a){return 0}}function e(b,B){null==a?d[b]=B(c.attr(b)||0):b==a&&(d=B(null==m?c.attr(b)||0:m))}var f=p(c).node,d={},B=f.querySelector(".svg---mgr");B||(B=v("rect"),v(B,{x:-9E9,y:-9E9,width:10,height:10,"class":"svg---mgr",
fill:"none"}),f.appendChild(B));switch(c.type){case "rect":e("rx",b),e("ry",h);case "image":e("width",b),e("height",h);case "text":e("x",b);e("y",h);break;case "circle":e("cx",b);e("cy",h);e("r",b);break;case "ellipse":e("cx",b);e("cy",h);e("rx",b);e("ry",h);break;case "line":e("x1",b);e("x2",b);e("y1",h);e("y2",h);break;case "marker":e("refX",b);e("markerWidth",b);e("refY",h);e("markerHeight",h);break;case "radialGradient":e("fx",b);e("fy",h);break;case "tspan":e("dx",b);e("dy",h);break;default:e(a,
b)}f.removeChild(B);return d}function q(c){y(c,"array")||(c=Array.prototype.slice.call(arguments,0));for(var a=0,b=0,m=this.node;this[a];)delete this[a++];for(a=0;a<c.length;a++)"set"==c[a].type?c[a].forEach(function(c){m.appendChild(c.node)}):m.appendChild(c[a].node);for(var h=m.childNodes,a=0;a<h.length;a++)this[b++]=x(h[a]);return this}function e(c){if(c.snap in E)return E[c.snap];var a=this.id=V(),b;try{b=c.ownerSVGElement}catch(m){}this.node=c;b&&(this.paper=new s(b));this.type=c.tagName;this.anims=
{};this._={transform:[]};c.snap=a;E[a]=this;"g"==this.type&&(this.add=q);if(this.type in{g:1,mask:1,pattern:1})for(var e in s.prototype)s.prototype[h](e)&&(this[e]=s.prototype[e])}function l(c){this.node=c}function r(c,a){var b=v(c);a.appendChild(b);return x(b)}function s(c,a){var b,m,f,d=s.prototype;if(c&&"svg"==c.tagName){if(c.snap in E)return E[c.snap];var l=c.ownerDocument;b=new e(c);m=c.getElementsByTagName("desc")[0];f=c.getElementsByTagName("defs")[0];m||(m=v("desc"),m.appendChild(l.createTextNode("Created with Snap")),
b.node.appendChild(m));f||(f=v("defs"),b.node.appendChild(f));b.defs=f;for(var ca in d)d[h](ca)&&(b[ca]=d[ca]);b.paper=b.root=b}else b=r("svg",G.doc.body),v(b.node,{height:a,version:1.1,width:c,xmlns:la});return b}function x(c){return!c||c instanceof e||c instanceof l?c:c.tagName&&"svg"==c.tagName.toLowerCase()?new s(c):c.tagName&&"object"==c.tagName.toLowerCase()&&"image/svg+xml"==c.type?new s(c.contentDocument.getElementsByTagName("svg")[0]):new e(c)}a.version="0.3.0";a.toString=function(){return"Snap v"+
this.version};a._={};var G={win:N,doc:N.document};a._.glob=G;var h="hasOwnProperty",J=String,K=parseFloat,U=parseInt,I=Math,P=I.max,Q=I.min,Y=I.abs,C=I.PI,aa="",$=Object.prototype.toString,F=/^\s*((#[a-f\d]{6})|(#[a-f\d]{3})|rgba?\(\s*([\d\.]+%?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+%?(?:\s*,\s*[\d\.]+%?)?)\s*\)|hsba?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\)|hsla?\(\s*([\d\.]+(?:deg|\xb0|%)?\s*,\s*[\d\.]+%?\s*,\s*[\d\.]+(?:%?\s*,\s*[\d\.]+)?%?)\s*\))\s*$/i;a._.separator=
RegExp("[,\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]+");var S=RegExp("[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*"),X={hs:1,rg:1},W=RegExp("([a-z])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)",
"ig"),ma=RegExp("([rstm])[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029,]*((-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*)+)","ig"),Z=RegExp("(-?\\d*\\.?\\d*(?:e[\\-+]?\\d+)?)[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*,?[\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*",
"ig"),na=0,ba="S"+(+new Date).toString(36),V=function(){return ba+(na++).toString(36)},m="http://www.w3.org/1999/xlink",la="http://www.w3.org/2000/svg",E={},ca=a.url=function(c){return"url('#"+c+"')"};a._.$=v;a._.id=V;a.format=function(){var c=/\{([^\}]+)\}/g,a=/(?:(?:^|\.)(.+?)(?=\[|\.|$|\()|\[('|")(.+?)\2\])(\(\))?/g,b=function(c,b,m){var h=m;b.replace(a,function(c,a,b,m,t){a=a||m;h&&(a in h&&(h=h[a]),"function"==typeof h&&t&&(h=h()))});return h=(null==h||h==m?c:h)+""};return function(a,m){return J(a).replace(c,
function(c,a){return b(c,a,m)})}}();a._.clone=M;a._.cacher=A;a.rad=z;a.deg=function(c){return 180*c/C%360};a.angle=w;a.is=y;a.snapTo=function(c,a,b){b=y(b,"finite")?b:10;if(y(c,"array"))for(var m=c.length;m--;){if(Y(c[m]-a)<=b)return c[m]}else{c=+c;m=a%c;if(m<b)return a-m;if(m>c-b)return a-m+c}return a};a.getRGB=A(function(c){if(!c||(c=J(c)).indexOf("-")+1)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};if("none"==c)return{r:-1,g:-1,b:-1,hex:"none",toString:ka};!X[h](c.toLowerCase().substring(0,
2))&&"#"!=c.charAt()&&(c=T(c));if(!c)return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka};var b,m,e,f,d;if(c=c.match(F)){c[2]&&(e=U(c[2].substring(5),16),m=U(c[2].substring(3,5),16),b=U(c[2].substring(1,3),16));c[3]&&(e=U((d=c[3].charAt(3))+d,16),m=U((d=c[3].charAt(2))+d,16),b=U((d=c[3].charAt(1))+d,16));c[4]&&(d=c[4].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b*=2.55),m=K(d[1]),"%"==d[1].slice(-1)&&(m*=2.55),e=K(d[2]),"%"==d[2].slice(-1)&&(e*=2.55),"rgba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),
d[3]&&"%"==d[3].slice(-1)&&(f/=100));if(c[5])return d=c[5].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsba"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsb2rgb(b,m,e,f);if(c[6])return d=c[6].split(S),b=K(d[0]),"%"==d[0].slice(-1)&&(b/=100),m=K(d[1]),"%"==d[1].slice(-1)&&(m/=100),e=K(d[2]),"%"==d[2].slice(-1)&&(e/=100),
"deg"!=d[0].slice(-3)&&"\u00b0"!=d[0].slice(-1)||(b/=360),"hsla"==c[1].toLowerCase().slice(0,4)&&(f=K(d[3])),d[3]&&"%"==d[3].slice(-1)&&(f/=100),a.hsl2rgb(b,m,e,f);b=Q(I.round(b),255);m=Q(I.round(m),255);e=Q(I.round(e),255);f=Q(P(f,0),1);c={r:b,g:m,b:e,toString:ka};c.hex="#"+(16777216|e|m<<8|b<<16).toString(16).slice(1);c.opacity=y(f,"finite")?f:1;return c}return{r:-1,g:-1,b:-1,hex:"none",error:1,toString:ka}},a);a.hsb=A(function(c,b,m){return a.hsb2rgb(c,b,m).hex});a.hsl=A(function(c,b,m){return a.hsl2rgb(c,
b,m).hex});a.rgb=A(function(c,a,b,m){if(y(m,"finite")){var e=I.round;return"rgba("+[e(c),e(a),e(b),+m.toFixed(2)]+")"}return"#"+(16777216|b|a<<8|c<<16).toString(16).slice(1)});var T=function(c){var a=G.doc.getElementsByTagName("head")[0]||G.doc.getElementsByTagName("svg")[0];T=A(function(c){if("red"==c.toLowerCase())return"rgb(255, 0, 0)";a.style.color="rgb(255, 0, 0)";a.style.color=c;c=G.doc.defaultView.getComputedStyle(a,aa).getPropertyValue("color");return"rgb(255, 0, 0)"==c?null:c});return T(c)},
qa=function(){return"hsb("+[this.h,this.s,this.b]+")"},ra=function(){return"hsl("+[this.h,this.s,this.l]+")"},ka=function(){return 1==this.opacity||null==this.opacity?this.hex:"rgba("+[this.r,this.g,this.b,this.opacity]+")"},D=function(c,b,m){null==b&&y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&(m=c.b,b=c.g,c=c.r);null==b&&y(c,string)&&(m=a.getRGB(c),c=m.r,b=m.g,m=m.b);if(1<c||1<b||1<m)c/=255,b/=255,m/=255;return[c,b,m]},oa=function(c,b,m,e){c=I.round(255*c);b=I.round(255*b);m=I.round(255*m);c={r:c,
g:b,b:m,opacity:y(e,"finite")?e:1,hex:a.rgb(c,b,m),toString:ka};y(e,"finite")&&(c.opacity=e);return c};a.color=function(c){var b;y(c,"object")&&"h"in c&&"s"in c&&"b"in c?(b=a.hsb2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):y(c,"object")&&"h"in c&&"s"in c&&"l"in c?(b=a.hsl2rgb(c),c.r=b.r,c.g=b.g,c.b=b.b,c.opacity=1,c.hex=b.hex):(y(c,"string")&&(c=a.getRGB(c)),y(c,"object")&&"r"in c&&"g"in c&&"b"in c&&!("error"in c)?(b=a.rgb2hsl(c),c.h=b.h,c.s=b.s,c.l=b.l,b=a.rgb2hsb(c),c.v=b.b):(c={hex:"none"},
c.r=c.g=c.b=c.h=c.s=c.v=c.l=-1,c.error=1));c.toString=ka;return c};a.hsb2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"b"in c&&(b=c.b,a=c.s,c=c.h,m=c.o);var e,h,d;c=360*c%360/60;d=b*a;a=d*(1-Y(c%2-1));b=e=h=b-d;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.hsl2rgb=function(c,a,b,m){y(c,"object")&&"h"in c&&"s"in c&&"l"in c&&(b=c.l,a=c.s,c=c.h);if(1<c||1<a||1<b)c/=360,a/=100,b/=100;var e,h,d;c=360*c%360/60;d=2*a*(0.5>b?b:1-b);a=d*(1-Y(c%2-1));b=e=
h=b-d/2;c=~~c;b+=[d,a,0,0,a,d][c];e+=[a,d,d,a,0,0][c];h+=[0,0,a,d,d,a][c];return oa(b,e,h,m)};a.rgb2hsb=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e;m=P(c,a,b);e=m-Q(c,a,b);c=((0==e?0:m==c?(a-b)/e:m==a?(b-c)/e+2:(c-a)/e+4)+360)%6*60/360;return{h:c,s:0==e?0:e/m,b:m,toString:qa}};a.rgb2hsl=function(c,a,b){b=D(c,a,b);c=b[0];a=b[1];b=b[2];var m,e,h;m=P(c,a,b);e=Q(c,a,b);h=m-e;c=((0==h?0:m==c?(a-b)/h:m==a?(b-c)/h+2:(c-a)/h+4)+360)%6*60/360;m=(m+e)/2;return{h:c,s:0==h?0:0.5>m?h/(2*m):h/(2-2*
m),l:m,toString:ra}};a.parsePathString=function(c){if(!c)return null;var b=a.path(c);if(b.arr)return a.path.clone(b.arr);var m={a:7,c:6,o:2,h:1,l:2,m:2,r:4,q:4,s:4,t:2,v:1,u:3,z:0},e=[];y(c,"array")&&y(c[0],"array")&&(e=a.path.clone(c));e.length||J(c).replace(W,function(c,a,b){var h=[];c=a.toLowerCase();b.replace(Z,function(c,a){a&&h.push(+a)});"m"==c&&2<h.length&&(e.push([a].concat(h.splice(0,2))),c="l",a="m"==a?"l":"L");"o"==c&&1==h.length&&e.push([a,h[0] ]);if("r"==c)e.push([a].concat(h));else for(;h.length>=
m[c]&&(e.push([a].concat(h.splice(0,m[c]))),m[c]););});e.toString=a.path.toString;b.arr=a.path.clone(e);return e};var O=a.parseTransformString=function(c){if(!c)return null;var b=[];y(c,"array")&&y(c[0],"array")&&(b=a.path.clone(c));b.length||J(c).replace(ma,function(c,a,m){var e=[];a.toLowerCase();m.replace(Z,function(c,a){a&&e.push(+a)});b.push([a].concat(e))});b.toString=a.path.toString;return b};a._.svgTransform2string=d;a._.rgTransform=RegExp("^[a-z][\t\n\x0B\f\r \u00a0\u1680\u180e\u2000\u2001\u2002\u2003\u2004\u2005\u2006\u2007\u2008\u2009\u200a\u202f\u205f\u3000\u2028\u2029]*-?\\.?\\d",
"i");a._.transform2matrix=f;a._unit2px=b;a._.getSomeDefs=u;a._.getSomeSVG=p;a.select=function(c){return x(G.doc.querySelector(c))};a.selectAll=function(c){c=G.doc.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};setInterval(function(){for(var c in E)if(E[h](c)){var a=E[c],b=a.node;("svg"!=a.type&&!b.ownerSVGElement||"svg"==a.type&&(!b.parentNode||"ownerSVGElement"in b.parentNode&&!b.ownerSVGElement))&&delete E[c]}},1E4);(function(c){function m(c){function a(c,
b){var m=v(c.node,b);(m=(m=m&&m.match(d))&&m[2])&&"#"==m.charAt()&&(m=m.substring(1))&&(f[m]=(f[m]||[]).concat(function(a){var m={};m[b]=ca(a);v(c.node,m)}))}function b(c){var a=v(c.node,"xlink:href");a&&"#"==a.charAt()&&(a=a.substring(1))&&(f[a]=(f[a]||[]).concat(function(a){c.attr("xlink:href","#"+a)}))}var e=c.selectAll("*"),h,d=/^\s*url\(("|'|)(.*)\1\)\s*$/;c=[];for(var f={},l=0,E=e.length;l<E;l++){h=e[l];a(h,"fill");a(h,"stroke");a(h,"filter");a(h,"mask");a(h,"clip-path");b(h);var t=v(h.node,
"id");t&&(v(h.node,{id:h.id}),c.push({old:t,id:h.id}))}l=0;for(E=c.length;l<E;l++)if(e=f[c[l].old])for(h=0,t=e.length;h<t;h++)e[h](c[l].id)}function e(c,a,b){return function(m){m=m.slice(c,a);1==m.length&&(m=m[0]);return b?b(m):m}}function d(c){return function(){var a=c?"<"+this.type:"",b=this.node.attributes,m=this.node.childNodes;if(c)for(var e=0,h=b.length;e<h;e++)a+=" "+b[e].name+'="'+b[e].value.replace(/"/g,'\\"')+'"';if(m.length){c&&(a+=">");e=0;for(h=m.length;e<h;e++)3==m[e].nodeType?a+=m[e].nodeValue:
1==m[e].nodeType&&(a+=x(m[e]).toString());c&&(a+="</"+this.type+">")}else c&&(a+="/>");return a}}c.attr=function(c,a){if(!c)return this;if(y(c,"string"))if(1<arguments.length){var b={};b[c]=a;c=b}else return k("snap.util.getattr."+c,this).firstDefined();for(var m in c)c[h](m)&&k("snap.util.attr."+m,this,c[m]);return this};c.getBBox=function(c){if(!a.Matrix||!a.path)return this.node.getBBox();var b=this,m=new a.Matrix;if(b.removed)return a._.box();for(;"use"==b.type;)if(c||(m=m.add(b.transform().localMatrix.translate(b.attr("x")||
0,b.attr("y")||0))),b.original)b=b.original;else var e=b.attr("xlink:href"),b=b.original=b.node.ownerDocument.getElementById(e.substring(e.indexOf("#")+1));var e=b._,h=a.path.get[b.type]||a.path.get.deflt;try{if(c)return e.bboxwt=h?a.path.getBBox(b.realPath=h(b)):a._.box(b.node.getBBox()),a._.box(e.bboxwt);b.realPath=h(b);b.matrix=b.transform().localMatrix;e.bbox=a.path.getBBox(a.path.map(b.realPath,m.add(b.matrix)));return a._.box(e.bbox)}catch(d){return a._.box()}};var f=function(){return this.string};
c.transform=function(c){var b=this._;if(null==c){var m=this;c=new a.Matrix(this.node.getCTM());for(var e=n(this),h=[e],d=new a.Matrix,l=e.toTransformString(),b=J(e)==J(this.matrix)?J(b.transform):l;"svg"!=m.type&&(m=m.parent());)h.push(n(m));for(m=h.length;m--;)d.add(h[m]);return{string:b,globalMatrix:c,totalMatrix:d,localMatrix:e,diffMatrix:c.clone().add(e.invert()),global:c.toTransformString(),total:d.toTransformString(),local:l,toString:f}}c instanceof a.Matrix?this.matrix=c:n(this,c);this.node&&
("linearGradient"==this.type||"radialGradient"==this.type?v(this.node,{gradientTransform:this.matrix}):"pattern"==this.type?v(this.node,{patternTransform:this.matrix}):v(this.node,{transform:this.matrix}));return this};c.parent=function(){return x(this.node.parentNode)};c.append=c.add=function(c){if(c){if("set"==c.type){var a=this;c.forEach(function(c){a.add(c)});return this}c=x(c);this.node.appendChild(c.node);c.paper=this.paper}return this};c.appendTo=function(c){c&&(c=x(c),c.append(this));return this};
c.prepend=function(c){if(c){if("set"==c.type){var a=this,b;c.forEach(function(c){b?b.after(c):a.prepend(c);b=c});return this}c=x(c);var m=c.parent();this.node.insertBefore(c.node,this.node.firstChild);this.add&&this.add();c.paper=this.paper;this.parent()&&this.parent().add();m&&m.add()}return this};c.prependTo=function(c){c=x(c);c.prepend(this);return this};c.before=function(c){if("set"==c.type){var a=this;c.forEach(function(c){var b=c.parent();a.node.parentNode.insertBefore(c.node,a.node);b&&b.add()});
this.parent().add();return this}c=x(c);var b=c.parent();this.node.parentNode.insertBefore(c.node,this.node);this.parent()&&this.parent().add();b&&b.add();c.paper=this.paper;return this};c.after=function(c){c=x(c);var a=c.parent();this.node.nextSibling?this.node.parentNode.insertBefore(c.node,this.node.nextSibling):this.node.parentNode.appendChild(c.node);this.parent()&&this.parent().add();a&&a.add();c.paper=this.paper;return this};c.insertBefore=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,
c.node);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.insertAfter=function(c){c=x(c);var a=this.parent();c.node.parentNode.insertBefore(this.node,c.node.nextSibling);this.paper=c.paper;a&&a.add();c.parent()&&c.parent().add();return this};c.remove=function(){var c=this.parent();this.node.parentNode&&this.node.parentNode.removeChild(this.node);delete this.paper;this.removed=!0;c&&c.add();return this};c.select=function(c){return x(this.node.querySelector(c))};c.selectAll=
function(c){c=this.node.querySelectorAll(c);for(var b=(a.set||Array)(),m=0;m<c.length;m++)b.push(x(c[m]));return b};c.asPX=function(c,a){null==a&&(a=this.attr(c));return+b(this,c,a)};c.use=function(){var c,a=this.node.id;a||(a=this.id,v(this.node,{id:a}));c="linearGradient"==this.type||"radialGradient"==this.type||"pattern"==this.type?r(this.type,this.node.parentNode):r("use",this.node.parentNode);v(c.node,{"xlink:href":"#"+a});c.original=this;return c};var l=/\S+/g;c.addClass=function(c){var a=(c||
"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h,d;if(a.length){for(e=0;d=a[e++];)h=m.indexOf(d),~h||m.push(d);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.removeClass=function(c){var a=(c||"").match(l)||[];c=this.node;var b=c.className.baseVal,m=b.match(l)||[],e,h;if(m.length){for(e=0;h=a[e++];)h=m.indexOf(h),~h&&m.splice(h,1);a=m.join(" ");b!=a&&(c.className.baseVal=a)}return this};c.hasClass=function(c){return!!~(this.node.className.baseVal.match(l)||[]).indexOf(c)};
c.toggleClass=function(c,a){if(null!=a)return a?this.addClass(c):this.removeClass(c);var b=(c||"").match(l)||[],m=this.node,e=m.className.baseVal,h=e.match(l)||[],d,f,E;for(d=0;E=b[d++];)f=h.indexOf(E),~f?h.splice(f,1):h.push(E);b=h.join(" ");e!=b&&(m.className.baseVal=b);return this};c.clone=function(){var c=x(this.node.cloneNode(!0));v(c.node,"id")&&v(c.node,{id:c.id});m(c);c.insertAfter(this);return c};c.toDefs=function(){u(this).appendChild(this.node);return this};c.pattern=c.toPattern=function(c,
a,b,m){var e=r("pattern",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,c=c.x);v(e.node,{x:c,y:a,width:b,height:m,patternUnits:"userSpaceOnUse",id:e.id,viewBox:[c,a,b,m].join(" ")});e.node.appendChild(this.node);return e};c.marker=function(c,a,b,m,e,h){var d=r("marker",u(this));null==c&&(c=this.getBBox());y(c,"object")&&"x"in c&&(a=c.y,b=c.width,m=c.height,e=c.refX||c.cx,h=c.refY||c.cy,c=c.x);v(d.node,{viewBox:[c,a,b,m].join(" "),markerWidth:b,markerHeight:m,
orient:"auto",refX:e||0,refY:h||0,id:d.id});d.node.appendChild(this.node);return d};var E=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);this.attr=c;this.dur=a;b&&(this.easing=b);m&&(this.callback=m)};a._.Animation=E;a.animation=function(c,a,b,m){return new E(c,a,b,m)};c.inAnim=function(){var c=[],a;for(a in this.anims)this.anims[h](a)&&function(a){c.push({anim:new E(a._attrs,a.dur,a.easing,a._callback),mina:a,curStatus:a.status(),status:function(c){return a.status(c)},stop:function(){a.stop()}})}(this.anims[a]);
return c};a.animate=function(c,a,b,m,e,h){"function"!=typeof e||e.length||(h=e,e=L.linear);var d=L.time();c=L(c,a,d,d+m,L.time,b,e);h&&k.once("mina.finish."+c.id,h);return c};c.stop=function(){for(var c=this.inAnim(),a=0,b=c.length;a<b;a++)c[a].stop();return this};c.animate=function(c,a,b,m){"function"!=typeof b||b.length||(m=b,b=L.linear);c instanceof E&&(m=c.callback,b=c.easing,a=b.dur,c=c.attr);var d=[],f=[],l={},t,ca,n,T=this,q;for(q in c)if(c[h](q)){T.equal?(n=T.equal(q,J(c[q])),t=n.from,ca=
n.to,n=n.f):(t=+T.attr(q),ca=+c[q]);var la=y(t,"array")?t.length:1;l[q]=e(d.length,d.length+la,n);d=d.concat(t);f=f.concat(ca)}t=L.time();var p=L(d,f,t,t+a,L.time,function(c){var a={},b;for(b in l)l[h](b)&&(a[b]=l[b](c));T.attr(a)},b);T.anims[p.id]=p;p._attrs=c;p._callback=m;k("snap.animcreated."+T.id,p);k.once("mina.finish."+p.id,function(){delete T.anims[p.id];m&&m.call(T)});k.once("mina.stop."+p.id,function(){delete T.anims[p.id]});return T};var T={};c.data=function(c,b){var m=T[this.id]=T[this.id]||
{};if(0==arguments.length)return k("snap.data.get."+this.id,this,m,null),m;if(1==arguments.length){if(a.is(c,"object")){for(var e in c)c[h](e)&&this.data(e,c[e]);return this}k("snap.data.get."+this.id,this,m[c],c);return m[c]}m[c]=b;k("snap.data.set."+this.id,this,b,c);return this};c.removeData=function(c){null==c?T[this.id]={}:T[this.id]&&delete T[this.id][c];return this};c.outerSVG=c.toString=d(1);c.innerSVG=d()})(e.prototype);a.parse=function(c){var a=G.doc.createDocumentFragment(),b=!0,m=G.doc.createElement("div");
c=J(c);c.match(/^\s*<\s*svg(?:\s|>)/)||(c="<svg>"+c+"</svg>",b=!1);m.innerHTML=c;if(c=m.getElementsByTagName("svg")[0])if(b)a=c;else for(;c.firstChild;)a.appendChild(c.firstChild);m.innerHTML=aa;return new l(a)};l.prototype.select=e.prototype.select;l.prototype.selectAll=e.prototype.selectAll;a.fragment=function(){for(var c=Array.prototype.slice.call(arguments,0),b=G.doc.createDocumentFragment(),m=0,e=c.length;m<e;m++){var h=c[m];h.node&&h.node.nodeType&&b.appendChild(h.node);h.nodeType&&b.appendChild(h);
"string"==typeof h&&b.appendChild(a.parse(h).node)}return new l(b)};a._.make=r;a._.wrap=x;s.prototype.el=function(c,a){var b=r(c,this.node);a&&b.attr(a);return b};k.on("snap.util.getattr",function(){var c=k.nt(),c=c.substring(c.lastIndexOf(".")+1),a=c.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});return pa[h](a)?this.node.ownerDocument.defaultView.getComputedStyle(this.node,null).getPropertyValue(a):v(this.node,c)});var pa={"alignment-baseline":0,"baseline-shift":0,clip:0,"clip-path":0,
"clip-rule":0,color:0,"color-interpolation":0,"color-interpolation-filters":0,"color-profile":0,"color-rendering":0,cursor:0,direction:0,display:0,"dominant-baseline":0,"enable-background":0,fill:0,"fill-opacity":0,"fill-rule":0,filter:0,"flood-color":0,"flood-opacity":0,font:0,"font-family":0,"font-size":0,"font-size-adjust":0,"font-stretch":0,"font-style":0,"font-variant":0,"font-weight":0,"glyph-orientation-horizontal":0,"glyph-orientation-vertical":0,"image-rendering":0,kerning:0,"letter-spacing":0,
"lighting-color":0,marker:0,"marker-end":0,"marker-mid":0,"marker-start":0,mask:0,opacity:0,overflow:0,"pointer-events":0,"shape-rendering":0,"stop-color":0,"stop-opacity":0,stroke:0,"stroke-dasharray":0,"stroke-dashoffset":0,"stroke-linecap":0,"stroke-linejoin":0,"stroke-miterlimit":0,"stroke-opacity":0,"stroke-width":0,"text-anchor":0,"text-decoration":0,"text-rendering":0,"unicode-bidi":0,visibility:0,"word-spacing":0,"writing-mode":0};k.on("snap.util.attr",function(c){var a=k.nt(),b={},a=a.substring(a.lastIndexOf(".")+
1);b[a]=c;var m=a.replace(/-(\w)/gi,function(c,a){return a.toUpperCase()}),a=a.replace(/[A-Z]/g,function(c){return"-"+c.toLowerCase()});pa[h](a)?this.node.style[m]=null==c?aa:c:v(this.node,b)});a.ajax=function(c,a,b,m){var e=new XMLHttpRequest,h=V();if(e){if(y(a,"function"))m=b,b=a,a=null;else if(y(a,"object")){var d=[],f;for(f in a)a.hasOwnProperty(f)&&d.push(encodeURIComponent(f)+"="+encodeURIComponent(a[f]));a=d.join("&")}e.open(a?"POST":"GET",c,!0);a&&(e.setRequestHeader("X-Requested-With","XMLHttpRequest"),
e.setRequestHeader("Content-type","application/x-www-form-urlencoded"));b&&(k.once("snap.ajax."+h+".0",b),k.once("snap.ajax."+h+".200",b),k.once("snap.ajax."+h+".304",b));e.onreadystatechange=function(){4==e.readyState&&k("snap.ajax."+h+"."+e.status,m,e)};if(4==e.readyState)return e;e.send(a);return e}};a.load=function(c,b,m){a.ajax(c,function(c){c=a.parse(c.responseText);m?b.call(m,c):b(c)})};a.getElementByPoint=function(c,a){var b,m,e=G.doc.elementFromPoint(c,a);if(G.win.opera&&"svg"==e.tagName){b=
e;m=b.getBoundingClientRect();b=b.ownerDocument;var h=b.body,d=b.documentElement;b=m.top+(g.win.pageYOffset||d.scrollTop||h.scrollTop)-(d.clientTop||h.clientTop||0);m=m.left+(g.win.pageXOffset||d.scrollLeft||h.scrollLeft)-(d.clientLeft||h.clientLeft||0);h=e.createSVGRect();h.x=c-m;h.y=a-b;h.width=h.height=1;b=e.getIntersectionList(h,null);b.length&&(e=b[b.length-1])}return e?x(e):null};a.plugin=function(c){c(a,e,s,G,l)};return G.win.Snap=a}();C.plugin(function(a,k,y,M,A){function w(a,d,f,b,q,e){null==
d&&"[object SVGMatrix]"==z.call(a)?(this.a=a.a,this.b=a.b,this.c=a.c,this.d=a.d,this.e=a.e,this.f=a.f):null!=a?(this.a=+a,this.b=+d,this.c=+f,this.d=+b,this.e=+q,this.f=+e):(this.a=1,this.c=this.b=0,this.d=1,this.f=this.e=0)}var z=Object.prototype.toString,d=String,f=Math;(function(n){function k(a){return a[0]*a[0]+a[1]*a[1]}function p(a){var d=f.sqrt(k(a));a[0]&&(a[0]/=d);a[1]&&(a[1]/=d)}n.add=function(a,d,e,f,n,p){var k=[[],[],[] ],u=[[this.a,this.c,this.e],[this.b,this.d,this.f],[0,0,1] ];d=[[a,
e,n],[d,f,p],[0,0,1] ];a&&a instanceof w&&(d=[[a.a,a.c,a.e],[a.b,a.d,a.f],[0,0,1] ]);for(a=0;3>a;a++)for(e=0;3>e;e++){for(f=n=0;3>f;f++)n+=u[a][f]*d[f][e];k[a][e]=n}this.a=k[0][0];this.b=k[1][0];this.c=k[0][1];this.d=k[1][1];this.e=k[0][2];this.f=k[1][2];return this};n.invert=function(){var a=this.a*this.d-this.b*this.c;return new w(this.d/a,-this.b/a,-this.c/a,this.a/a,(this.c*this.f-this.d*this.e)/a,(this.b*this.e-this.a*this.f)/a)};n.clone=function(){return new w(this.a,this.b,this.c,this.d,this.e,
this.f)};n.translate=function(a,d){return this.add(1,0,0,1,a,d)};n.scale=function(a,d,e,f){null==d&&(d=a);(e||f)&&this.add(1,0,0,1,e,f);this.add(a,0,0,d,0,0);(e||f)&&this.add(1,0,0,1,-e,-f);return this};n.rotate=function(b,d,e){b=a.rad(b);d=d||0;e=e||0;var l=+f.cos(b).toFixed(9);b=+f.sin(b).toFixed(9);this.add(l,b,-b,l,d,e);return this.add(1,0,0,1,-d,-e)};n.x=function(a,d){return a*this.a+d*this.c+this.e};n.y=function(a,d){return a*this.b+d*this.d+this.f};n.get=function(a){return+this[d.fromCharCode(97+
a)].toFixed(4)};n.toString=function(){return"matrix("+[this.get(0),this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)].join()+")"};n.offset=function(){return[this.e.toFixed(4),this.f.toFixed(4)]};n.determinant=function(){return this.a*this.d-this.b*this.c};n.split=function(){var b={};b.dx=this.e;b.dy=this.f;var d=[[this.a,this.c],[this.b,this.d] ];b.scalex=f.sqrt(k(d[0]));p(d[0]);b.shear=d[0][0]*d[1][0]+d[0][1]*d[1][1];d[1]=[d[1][0]-d[0][0]*b.shear,d[1][1]-d[0][1]*b.shear];b.scaley=f.sqrt(k(d[1]));
p(d[1]);b.shear/=b.scaley;0>this.determinant()&&(b.scalex=-b.scalex);var e=-d[0][1],d=d[1][1];0>d?(b.rotate=a.deg(f.acos(d)),0>e&&(b.rotate=360-b.rotate)):b.rotate=a.deg(f.asin(e));b.isSimple=!+b.shear.toFixed(9)&&(b.scalex.toFixed(9)==b.scaley.toFixed(9)||!b.rotate);b.isSuperSimple=!+b.shear.toFixed(9)&&b.scalex.toFixed(9)==b.scaley.toFixed(9)&&!b.rotate;b.noRotation=!+b.shear.toFixed(9)&&!b.rotate;return b};n.toTransformString=function(a){a=a||this.split();if(+a.shear.toFixed(9))return"m"+[this.get(0),
this.get(1),this.get(2),this.get(3),this.get(4),this.get(5)];a.scalex=+a.scalex.toFixed(4);a.scaley=+a.scaley.toFixed(4);a.rotate=+a.rotate.toFixed(4);return(a.dx||a.dy?"t"+[+a.dx.toFixed(4),+a.dy.toFixed(4)]:"")+(1!=a.scalex||1!=a.scaley?"s"+[a.scalex,a.scaley,0,0]:"")+(a.rotate?"r"+[+a.rotate.toFixed(4),0,0]:"")}})(w.prototype);a.Matrix=w;a.matrix=function(a,d,f,b,k,e){return new w(a,d,f,b,k,e)}});C.plugin(function(a,v,y,M,A){function w(h){return function(d){k.stop();d instanceof A&&1==d.node.childNodes.length&&
("radialGradient"==d.node.firstChild.tagName||"linearGradient"==d.node.firstChild.tagName||"pattern"==d.node.firstChild.tagName)&&(d=d.node.firstChild,b(this).appendChild(d),d=u(d));if(d instanceof v)if("radialGradient"==d.type||"linearGradient"==d.type||"pattern"==d.type){d.node.id||e(d.node,{id:d.id});var f=l(d.node.id)}else f=d.attr(h);else f=a.color(d),f.error?(f=a(b(this).ownerSVGElement).gradient(d))?(f.node.id||e(f.node,{id:f.id}),f=l(f.node.id)):f=d:f=r(f);d={};d[h]=f;e(this.node,d);this.node.style[h]=
x}}function z(a){k.stop();a==+a&&(a+="px");this.node.style.fontSize=a}function d(a){var b=[];a=a.childNodes;for(var e=0,f=a.length;e<f;e++){var l=a[e];3==l.nodeType&&b.push(l.nodeValue);"tspan"==l.tagName&&(1==l.childNodes.length&&3==l.firstChild.nodeType?b.push(l.firstChild.nodeValue):b.push(d(l)))}return b}function f(){k.stop();return this.node.style.fontSize}var n=a._.make,u=a._.wrap,p=a.is,b=a._.getSomeDefs,q=/^url\(#?([^)]+)\)$/,e=a._.$,l=a.url,r=String,s=a._.separator,x="";k.on("snap.util.attr.mask",
function(a){if(a instanceof v||a instanceof A){k.stop();a instanceof A&&1==a.node.childNodes.length&&(a=a.node.firstChild,b(this).appendChild(a),a=u(a));if("mask"==a.type)var d=a;else d=n("mask",b(this)),d.node.appendChild(a.node);!d.node.id&&e(d.node,{id:d.id});e(this.node,{mask:l(d.id)})}});(function(a){k.on("snap.util.attr.clip",a);k.on("snap.util.attr.clip-path",a);k.on("snap.util.attr.clipPath",a)})(function(a){if(a instanceof v||a instanceof A){k.stop();if("clipPath"==a.type)var d=a;else d=
n("clipPath",b(this)),d.node.appendChild(a.node),!d.node.id&&e(d.node,{id:d.id});e(this.node,{"clip-path":l(d.id)})}});k.on("snap.util.attr.fill",w("fill"));k.on("snap.util.attr.stroke",w("stroke"));var G=/^([lr])(?:\(([^)]*)\))?(.*)$/i;k.on("snap.util.grad.parse",function(a){a=r(a);var b=a.match(G);if(!b)return null;a=b[1];var e=b[2],b=b[3],e=e.split(/\s*,\s*/).map(function(a){return+a==a?+a:a});1==e.length&&0==e[0]&&(e=[]);b=b.split("-");b=b.map(function(a){a=a.split(":");var b={color:a[0]};a[1]&&
(b.offset=parseFloat(a[1]));return b});return{type:a,params:e,stops:b}});k.on("snap.util.attr.d",function(b){k.stop();p(b,"array")&&p(b[0],"array")&&(b=a.path.toString.call(b));b=r(b);b.match(/[ruo]/i)&&(b=a.path.toAbsolute(b));e(this.node,{d:b})})(-1);k.on("snap.util.attr.#text",function(a){k.stop();a=r(a);for(a=M.doc.createTextNode(a);this.node.firstChild;)this.node.removeChild(this.node.firstChild);this.node.appendChild(a)})(-1);k.on("snap.util.attr.path",function(a){k.stop();this.attr({d:a})})(-1);
k.on("snap.util.attr.class",function(a){k.stop();this.node.className.baseVal=a})(-1);k.on("snap.util.attr.viewBox",function(a){a=p(a,"object")&&"x"in a?[a.x,a.y,a.width,a.height].join(" "):p(a,"array")?a.join(" "):a;e(this.node,{viewBox:a});k.stop()})(-1);k.on("snap.util.attr.transform",function(a){this.transform(a);k.stop()})(-1);k.on("snap.util.attr.r",function(a){"rect"==this.type&&(k.stop(),e(this.node,{rx:a,ry:a}))})(-1);k.on("snap.util.attr.textpath",function(a){k.stop();if("text"==this.type){var d,
f;if(!a&&this.textPath){for(a=this.textPath;a.node.firstChild;)this.node.appendChild(a.node.firstChild);a.remove();delete this.textPath}else if(p(a,"string")?(d=b(this),a=u(d.parentNode).path(a),d.appendChild(a.node),d=a.id,a.attr({id:d})):(a=u(a),a instanceof v&&(d=a.attr("id"),d||(d=a.id,a.attr({id:d})))),d)if(a=this.textPath,f=this.node,a)a.attr({"xlink:href":"#"+d});else{for(a=e("textPath",{"xlink:href":"#"+d});f.firstChild;)a.appendChild(f.firstChild);f.appendChild(a);this.textPath=u(a)}}})(-1);
k.on("snap.util.attr.text",function(a){if("text"==this.type){for(var b=this.node,d=function(a){var b=e("tspan");if(p(a,"array"))for(var f=0;f<a.length;f++)b.appendChild(d(a[f]));else b.appendChild(M.doc.createTextNode(a));b.normalize&&b.normalize();return b};b.firstChild;)b.removeChild(b.firstChild);for(a=d(a);a.firstChild;)b.appendChild(a.firstChild)}k.stop()})(-1);k.on("snap.util.attr.fontSize",z)(-1);k.on("snap.util.attr.font-size",z)(-1);k.on("snap.util.getattr.transform",function(){k.stop();
return this.transform()})(-1);k.on("snap.util.getattr.textpath",function(){k.stop();return this.textPath})(-1);(function(){function b(d){return function(){k.stop();var b=M.doc.defaultView.getComputedStyle(this.node,null).getPropertyValue("marker-"+d);return"none"==b?b:a(M.doc.getElementById(b.match(q)[1]))}}function d(a){return function(b){k.stop();var d="marker"+a.charAt(0).toUpperCase()+a.substring(1);if(""==b||!b)this.node.style[d]="none";else if("marker"==b.type){var f=b.node.id;f||e(b.node,{id:b.id});
this.node.style[d]=l(f)}}}k.on("snap.util.getattr.marker-end",b("end"))(-1);k.on("snap.util.getattr.markerEnd",b("end"))(-1);k.on("snap.util.getattr.marker-start",b("start"))(-1);k.on("snap.util.getattr.markerStart",b("start"))(-1);k.on("snap.util.getattr.marker-mid",b("mid"))(-1);k.on("snap.util.getattr.markerMid",b("mid"))(-1);k.on("snap.util.attr.marker-end",d("end"))(-1);k.on("snap.util.attr.markerEnd",d("end"))(-1);k.on("snap.util.attr.marker-start",d("start"))(-1);k.on("snap.util.attr.markerStart",
d("start"))(-1);k.on("snap.util.attr.marker-mid",d("mid"))(-1);k.on("snap.util.attr.markerMid",d("mid"))(-1)})();k.on("snap.util.getattr.r",function(){if("rect"==this.type&&e(this.node,"rx")==e(this.node,"ry"))return k.stop(),e(this.node,"rx")})(-1);k.on("snap.util.getattr.text",function(){if("text"==this.type||"tspan"==this.type){k.stop();var a=d(this.node);return 1==a.length?a[0]:a}})(-1);k.on("snap.util.getattr.#text",function(){return this.node.textContent})(-1);k.on("snap.util.getattr.viewBox",
function(){k.stop();var b=e(this.node,"viewBox");if(b)return b=b.split(s),a._.box(+b[0],+b[1],+b[2],+b[3])})(-1);k.on("snap.util.getattr.points",function(){var a=e(this.node,"points");k.stop();if(a)return a.split(s)})(-1);k.on("snap.util.getattr.path",function(){var a=e(this.node,"d");k.stop();return a})(-1);k.on("snap.util.getattr.class",function(){return this.node.className.baseVal})(-1);k.on("snap.util.getattr.fontSize",f)(-1);k.on("snap.util.getattr.font-size",f)(-1)});C.plugin(function(a,v,y,
M,A){function w(a){return a}function z(a){return function(b){return+b.toFixed(3)+a}}var d={"+":function(a,b){return a+b},"-":function(a,b){return a-b},"/":function(a,b){return a/b},"*":function(a,b){return a*b}},f=String,n=/[a-z]+$/i,u=/^\s*([+\-\/*])\s*=\s*([\d.eE+\-]+)\s*([^\d\s]+)?\s*$/;k.on("snap.util.attr",function(a){if(a=f(a).match(u)){var b=k.nt(),b=b.substring(b.lastIndexOf(".")+1),q=this.attr(b),e={};k.stop();var l=a[3]||"",r=q.match(n),s=d[a[1] ];r&&r==l?a=s(parseFloat(q),+a[2]):(q=this.asPX(b),
a=s(this.asPX(b),this.asPX(b,a[2]+l)));isNaN(q)||isNaN(a)||(e[b]=a,this.attr(e))}})(-10);k.on("snap.util.equal",function(a,b){var q=f(this.attr(a)||""),e=f(b).match(u);if(e){k.stop();var l=e[3]||"",r=q.match(n),s=d[e[1] ];if(r&&r==l)return{from:parseFloat(q),to:s(parseFloat(q),+e[2]),f:z(r)};q=this.asPX(a);return{from:q,to:s(q,this.asPX(a,e[2]+l)),f:w}}})(-10)});C.plugin(function(a,v,y,M,A){var w=y.prototype,z=a.is;w.rect=function(a,d,k,p,b,q){var e;null==q&&(q=b);z(a,"object")&&"[object Object]"==
a?e=a:null!=a&&(e={x:a,y:d,width:k,height:p},null!=b&&(e.rx=b,e.ry=q));return this.el("rect",e)};w.circle=function(a,d,k){var p;z(a,"object")&&"[object Object]"==a?p=a:null!=a&&(p={cx:a,cy:d,r:k});return this.el("circle",p)};var d=function(){function a(){this.parentNode.removeChild(this)}return function(d,k){var p=M.doc.createElement("img"),b=M.doc.body;p.style.cssText="position:absolute;left:-9999em;top:-9999em";p.onload=function(){k.call(p);p.onload=p.onerror=null;b.removeChild(p)};p.onerror=a;
b.appendChild(p);p.src=d}}();w.image=function(f,n,k,p,b){var q=this.el("image");if(z(f,"object")&&"src"in f)q.attr(f);else if(null!=f){var e={"xlink:href":f,preserveAspectRatio:"none"};null!=n&&null!=k&&(e.x=n,e.y=k);null!=p&&null!=b?(e.width=p,e.height=b):d(f,function(){a._.$(q.node,{width:this.offsetWidth,height:this.offsetHeight})});a._.$(q.node,e)}return q};w.ellipse=function(a,d,k,p){var b;z(a,"object")&&"[object Object]"==a?b=a:null!=a&&(b={cx:a,cy:d,rx:k,ry:p});return this.el("ellipse",b)};
w.path=function(a){var d;z(a,"object")&&!z(a,"array")?d=a:a&&(d={d:a});return this.el("path",d)};w.group=w.g=function(a){var d=this.el("g");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.svg=function(a,d,k,p,b,q,e,l){var r={};z(a,"object")&&null==d?r=a:(null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l]));return this.el("svg",r)};w.mask=function(a){var d=
this.el("mask");1==arguments.length&&a&&!a.type?d.attr(a):arguments.length&&d.add(Array.prototype.slice.call(arguments,0));return d};w.ptrn=function(a,d,k,p,b,q,e,l){if(z(a,"object"))var r=a;else arguments.length?(r={},null!=a&&(r.x=a),null!=d&&(r.y=d),null!=k&&(r.width=k),null!=p&&(r.height=p),null!=b&&null!=q&&null!=e&&null!=l&&(r.viewBox=[b,q,e,l])):r={patternUnits:"userSpaceOnUse"};return this.el("pattern",r)};w.use=function(a){return null!=a?(make("use",this.node),a instanceof v&&(a.attr("id")||
a.attr({id:ID()}),a=a.attr("id")),this.el("use",{"xlink:href":a})):v.prototype.use.call(this)};w.text=function(a,d,k){var p={};z(a,"object")?p=a:null!=a&&(p={x:a,y:d,text:k||""});return this.el("text",p)};w.line=function(a,d,k,p){var b={};z(a,"object")?b=a:null!=a&&(b={x1:a,x2:k,y1:d,y2:p});return this.el("line",b)};w.polyline=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polyline",d)};
w.polygon=function(a){1<arguments.length&&(a=Array.prototype.slice.call(arguments,0));var d={};z(a,"object")&&!z(a,"array")?d=a:null!=a&&(d={points:a});return this.el("polygon",d)};(function(){function d(){return this.selectAll("stop")}function n(b,d){var f=e("stop"),k={offset:+d+"%"};b=a.color(b);k["stop-color"]=b.hex;1>b.opacity&&(k["stop-opacity"]=b.opacity);e(f,k);this.node.appendChild(f);return this}function u(){if("linearGradient"==this.type){var b=e(this.node,"x1")||0,d=e(this.node,"x2")||
1,f=e(this.node,"y1")||0,k=e(this.node,"y2")||0;return a._.box(b,f,math.abs(d-b),math.abs(k-f))}b=this.node.r||0;return a._.box((this.node.cx||0.5)-b,(this.node.cy||0.5)-b,2*b,2*b)}function p(a,d){function f(a,b){for(var d=(b-u)/(a-w),e=w;e<a;e++)h[e].offset=+(+u+d*(e-w)).toFixed(2);w=a;u=b}var n=k("snap.util.grad.parse",null,d).firstDefined(),p;if(!n)return null;n.params.unshift(a);p="l"==n.type.toLowerCase()?b.apply(0,n.params):q.apply(0,n.params);n.type!=n.type.toLowerCase()&&e(p.node,{gradientUnits:"userSpaceOnUse"});
var h=n.stops,n=h.length,u=0,w=0;n--;for(var v=0;v<n;v++)"offset"in h[v]&&f(v,h[v].offset);h[n].offset=h[n].offset||100;f(n,h[n].offset);for(v=0;v<=n;v++){var y=h[v];p.addStop(y.color,y.offset)}return p}function b(b,k,p,q,w){b=a._.make("linearGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{x1:k,y1:p,x2:q,y2:w});return b}function q(b,k,p,q,w,h){b=a._.make("radialGradient",b);b.stops=d;b.addStop=n;b.getBBox=u;null!=k&&e(b.node,{cx:k,cy:p,r:q});null!=w&&null!=h&&e(b.node,{fx:w,fy:h});
return b}var e=a._.$;w.gradient=function(a){return p(this.defs,a)};w.gradientLinear=function(a,d,e,f){return b(this.defs,a,d,e,f)};w.gradientRadial=function(a,b,d,e,f){return q(this.defs,a,b,d,e,f)};w.toString=function(){var b=this.node.ownerDocument,d=b.createDocumentFragment(),b=b.createElement("div"),e=this.node.cloneNode(!0);d.appendChild(b);b.appendChild(e);a._.$(e,{xmlns:"http://www.w3.org/2000/svg"});b=b.innerHTML;d.removeChild(d.firstChild);return b};w.clear=function(){for(var a=this.node.firstChild,
b;a;)b=a.nextSibling,"defs"!=a.tagName?a.parentNode.removeChild(a):w.clear.call({node:a}),a=b}})()});C.plugin(function(a,k,y,M){function A(a){var b=A.ps=A.ps||{};b[a]?b[a].sleep=100:b[a]={sleep:100};setTimeout(function(){for(var d in b)b[L](d)&&d!=a&&(b[d].sleep--,!b[d].sleep&&delete b[d])});return b[a]}function w(a,b,d,e){null==a&&(a=b=d=e=0);null==b&&(b=a.y,d=a.width,e=a.height,a=a.x);return{x:a,y:b,width:d,w:d,height:e,h:e,x2:a+d,y2:b+e,cx:a+d/2,cy:b+e/2,r1:F.min(d,e)/2,r2:F.max(d,e)/2,r0:F.sqrt(d*
d+e*e)/2,path:s(a,b,d,e),vb:[a,b,d,e].join(" ")}}function z(){return this.join(",").replace(N,"$1")}function d(a){a=C(a);a.toString=z;return a}function f(a,b,d,h,f,k,l,n,p){if(null==p)return e(a,b,d,h,f,k,l,n);if(0>p||e(a,b,d,h,f,k,l,n)<p)p=void 0;else{var q=0.5,O=1-q,s;for(s=e(a,b,d,h,f,k,l,n,O);0.01<Z(s-p);)q/=2,O+=(s<p?1:-1)*q,s=e(a,b,d,h,f,k,l,n,O);p=O}return u(a,b,d,h,f,k,l,n,p)}function n(b,d){function e(a){return+(+a).toFixed(3)}return a._.cacher(function(a,h,l){a instanceof k&&(a=a.attr("d"));
a=I(a);for(var n,p,D,q,O="",s={},c=0,t=0,r=a.length;t<r;t++){D=a[t];if("M"==D[0])n=+D[1],p=+D[2];else{q=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6]);if(c+q>h){if(d&&!s.start){n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c);O+=["C"+e(n.start.x),e(n.start.y),e(n.m.x),e(n.m.y),e(n.x),e(n.y)];if(l)return O;s.start=O;O=["M"+e(n.x),e(n.y)+"C"+e(n.n.x),e(n.n.y),e(n.end.x),e(n.end.y),e(D[5]),e(D[6])].join();c+=q;n=+D[5];p=+D[6];continue}if(!b&&!d)return n=f(n,p,D[1],D[2],D[3],D[4],D[5],D[6],h-c)}c+=q;n=+D[5];p=+D[6]}O+=
D.shift()+D}s.end=O;return n=b?c:d?s:u(n,p,D[0],D[1],D[2],D[3],D[4],D[5],1)},null,a._.clone)}function u(a,b,d,e,h,f,k,l,n){var p=1-n,q=ma(p,3),s=ma(p,2),c=n*n,t=c*n,r=q*a+3*s*n*d+3*p*n*n*h+t*k,q=q*b+3*s*n*e+3*p*n*n*f+t*l,s=a+2*n*(d-a)+c*(h-2*d+a),t=b+2*n*(e-b)+c*(f-2*e+b),x=d+2*n*(h-d)+c*(k-2*h+d),c=e+2*n*(f-e)+c*(l-2*f+e);a=p*a+n*d;b=p*b+n*e;h=p*h+n*k;f=p*f+n*l;l=90-180*F.atan2(s-x,t-c)/S;return{x:r,y:q,m:{x:s,y:t},n:{x:x,y:c},start:{x:a,y:b},end:{x:h,y:f},alpha:l}}function p(b,d,e,h,f,n,k,l){a.is(b,
"array")||(b=[b,d,e,h,f,n,k,l]);b=U.apply(null,b);return w(b.min.x,b.min.y,b.max.x-b.min.x,b.max.y-b.min.y)}function b(a,b,d){return b>=a.x&&b<=a.x+a.width&&d>=a.y&&d<=a.y+a.height}function q(a,d){a=w(a);d=w(d);return b(d,a.x,a.y)||b(d,a.x2,a.y)||b(d,a.x,a.y2)||b(d,a.x2,a.y2)||b(a,d.x,d.y)||b(a,d.x2,d.y)||b(a,d.x,d.y2)||b(a,d.x2,d.y2)||(a.x<d.x2&&a.x>d.x||d.x<a.x2&&d.x>a.x)&&(a.y<d.y2&&a.y>d.y||d.y<a.y2&&d.y>a.y)}function e(a,b,d,e,h,f,n,k,l){null==l&&(l=1);l=(1<l?1:0>l?0:l)/2;for(var p=[-0.1252,
0.1252,-0.3678,0.3678,-0.5873,0.5873,-0.7699,0.7699,-0.9041,0.9041,-0.9816,0.9816],q=[0.2491,0.2491,0.2335,0.2335,0.2032,0.2032,0.1601,0.1601,0.1069,0.1069,0.0472,0.0472],s=0,c=0;12>c;c++)var t=l*p[c]+l,r=t*(t*(-3*a+9*d-9*h+3*n)+6*a-12*d+6*h)-3*a+3*d,t=t*(t*(-3*b+9*e-9*f+3*k)+6*b-12*e+6*f)-3*b+3*e,s=s+q[c]*F.sqrt(r*r+t*t);return l*s}function l(a,b,d){a=I(a);b=I(b);for(var h,f,l,n,k,s,r,O,x,c,t=d?0:[],w=0,v=a.length;w<v;w++)if(x=a[w],"M"==x[0])h=k=x[1],f=s=x[2];else{"C"==x[0]?(x=[h,f].concat(x.slice(1)),
h=x[6],f=x[7]):(x=[h,f,h,f,k,s,k,s],h=k,f=s);for(var G=0,y=b.length;G<y;G++)if(c=b[G],"M"==c[0])l=r=c[1],n=O=c[2];else{"C"==c[0]?(c=[l,n].concat(c.slice(1)),l=c[6],n=c[7]):(c=[l,n,l,n,r,O,r,O],l=r,n=O);var z;var K=x,B=c;z=d;var H=p(K),J=p(B);if(q(H,J)){for(var H=e.apply(0,K),J=e.apply(0,B),H=~~(H/8),J=~~(J/8),U=[],A=[],F={},M=z?0:[],P=0;P<H+1;P++){var C=u.apply(0,K.concat(P/H));U.push({x:C.x,y:C.y,t:P/H})}for(P=0;P<J+1;P++)C=u.apply(0,B.concat(P/J)),A.push({x:C.x,y:C.y,t:P/J});for(P=0;P<H;P++)for(K=
0;K<J;K++){var Q=U[P],L=U[P+1],B=A[K],C=A[K+1],N=0.001>Z(L.x-Q.x)?"y":"x",S=0.001>Z(C.x-B.x)?"y":"x",R;R=Q.x;var Y=Q.y,V=L.x,ea=L.y,fa=B.x,ga=B.y,ha=C.x,ia=C.y;if(W(R,V)<X(fa,ha)||X(R,V)>W(fa,ha)||W(Y,ea)<X(ga,ia)||X(Y,ea)>W(ga,ia))R=void 0;else{var $=(R*ea-Y*V)*(fa-ha)-(R-V)*(fa*ia-ga*ha),aa=(R*ea-Y*V)*(ga-ia)-(Y-ea)*(fa*ia-ga*ha),ja=(R-V)*(ga-ia)-(Y-ea)*(fa-ha);if(ja){var $=$/ja,aa=aa/ja,ja=+$.toFixed(2),ba=+aa.toFixed(2);R=ja<+X(R,V).toFixed(2)||ja>+W(R,V).toFixed(2)||ja<+X(fa,ha).toFixed(2)||
ja>+W(fa,ha).toFixed(2)||ba<+X(Y,ea).toFixed(2)||ba>+W(Y,ea).toFixed(2)||ba<+X(ga,ia).toFixed(2)||ba>+W(ga,ia).toFixed(2)?void 0:{x:$,y:aa}}else R=void 0}R&&F[R.x.toFixed(4)]!=R.y.toFixed(4)&&(F[R.x.toFixed(4)]=R.y.toFixed(4),Q=Q.t+Z((R[N]-Q[N])/(L[N]-Q[N]))*(L.t-Q.t),B=B.t+Z((R[S]-B[S])/(C[S]-B[S]))*(C.t-B.t),0<=Q&&1>=Q&&0<=B&&1>=B&&(z?M++:M.push({x:R.x,y:R.y,t1:Q,t2:B})))}z=M}else z=z?0:[];if(d)t+=z;else{H=0;for(J=z.length;H<J;H++)z[H].segment1=w,z[H].segment2=G,z[H].bez1=x,z[H].bez2=c;t=t.concat(z)}}}return t}
function r(a){var b=A(a);if(b.bbox)return C(b.bbox);if(!a)return w();a=I(a);for(var d=0,e=0,h=[],f=[],l,n=0,k=a.length;n<k;n++)l=a[n],"M"==l[0]?(d=l[1],e=l[2],h.push(d),f.push(e)):(d=U(d,e,l[1],l[2],l[3],l[4],l[5],l[6]),h=h.concat(d.min.x,d.max.x),f=f.concat(d.min.y,d.max.y),d=l[5],e=l[6]);a=X.apply(0,h);l=X.apply(0,f);h=W.apply(0,h);f=W.apply(0,f);f=w(a,l,h-a,f-l);b.bbox=C(f);return f}function s(a,b,d,e,h){if(h)return[["M",+a+ +h,b],["l",d-2*h,0],["a",h,h,0,0,1,h,h],["l",0,e-2*h],["a",h,h,0,0,1,
-h,h],["l",2*h-d,0],["a",h,h,0,0,1,-h,-h],["l",0,2*h-e],["a",h,h,0,0,1,h,-h],["z"] ];a=[["M",a,b],["l",d,0],["l",0,e],["l",-d,0],["z"] ];a.toString=z;return a}function x(a,b,d,e,h){null==h&&null==e&&(e=d);a=+a;b=+b;d=+d;e=+e;if(null!=h){var f=Math.PI/180,l=a+d*Math.cos(-e*f);a+=d*Math.cos(-h*f);var n=b+d*Math.sin(-e*f);b+=d*Math.sin(-h*f);d=[["M",l,n],["A",d,d,0,+(180<h-e),0,a,b] ]}else d=[["M",a,b],["m",0,-e],["a",d,e,0,1,1,0,2*e],["a",d,e,0,1,1,0,-2*e],["z"] ];d.toString=z;return d}function G(b){var e=
A(b);if(e.abs)return d(e.abs);Q(b,"array")&&Q(b&&b[0],"array")||(b=a.parsePathString(b));if(!b||!b.length)return[["M",0,0] ];var h=[],f=0,l=0,n=0,k=0,p=0;"M"==b[0][0]&&(f=+b[0][1],l=+b[0][2],n=f,k=l,p++,h[0]=["M",f,l]);for(var q=3==b.length&&"M"==b[0][0]&&"R"==b[1][0].toUpperCase()&&"Z"==b[2][0].toUpperCase(),s,r,w=p,c=b.length;w<c;w++){h.push(s=[]);r=b[w];p=r[0];if(p!=p.toUpperCase())switch(s[0]=p.toUpperCase(),s[0]){case "A":s[1]=r[1];s[2]=r[2];s[3]=r[3];s[4]=r[4];s[5]=r[5];s[6]=+r[6]+f;s[7]=+r[7]+
l;break;case "V":s[1]=+r[1]+l;break;case "H":s[1]=+r[1]+f;break;case "R":for(var t=[f,l].concat(r.slice(1)),u=2,v=t.length;u<v;u++)t[u]=+t[u]+f,t[++u]=+t[u]+l;h.pop();h=h.concat(P(t,q));break;case "O":h.pop();t=x(f,l,r[1],r[2]);t.push(t[0]);h=h.concat(t);break;case "U":h.pop();h=h.concat(x(f,l,r[1],r[2],r[3]));s=["U"].concat(h[h.length-1].slice(-2));break;case "M":n=+r[1]+f,k=+r[2]+l;default:for(u=1,v=r.length;u<v;u++)s[u]=+r[u]+(u%2?f:l)}else if("R"==p)t=[f,l].concat(r.slice(1)),h.pop(),h=h.concat(P(t,
q)),s=["R"].concat(r.slice(-2));else if("O"==p)h.pop(),t=x(f,l,r[1],r[2]),t.push(t[0]),h=h.concat(t);else if("U"==p)h.pop(),h=h.concat(x(f,l,r[1],r[2],r[3])),s=["U"].concat(h[h.length-1].slice(-2));else for(t=0,u=r.length;t<u;t++)s[t]=r[t];p=p.toUpperCase();if("O"!=p)switch(s[0]){case "Z":f=+n;l=+k;break;case "H":f=s[1];break;case "V":l=s[1];break;case "M":n=s[s.length-2],k=s[s.length-1];default:f=s[s.length-2],l=s[s.length-1]}}h.toString=z;e.abs=d(h);return h}function h(a,b,d,e){return[a,b,d,e,d,
e]}function J(a,b,d,e,h,f){var l=1/3,n=2/3;return[l*a+n*d,l*b+n*e,l*h+n*d,l*f+n*e,h,f]}function K(b,d,e,h,f,l,n,k,p,s){var r=120*S/180,q=S/180*(+f||0),c=[],t,x=a._.cacher(function(a,b,c){var d=a*F.cos(c)-b*F.sin(c);a=a*F.sin(c)+b*F.cos(c);return{x:d,y:a}});if(s)v=s[0],t=s[1],l=s[2],u=s[3];else{t=x(b,d,-q);b=t.x;d=t.y;t=x(k,p,-q);k=t.x;p=t.y;F.cos(S/180*f);F.sin(S/180*f);t=(b-k)/2;v=(d-p)/2;u=t*t/(e*e)+v*v/(h*h);1<u&&(u=F.sqrt(u),e*=u,h*=u);var u=e*e,w=h*h,u=(l==n?-1:1)*F.sqrt(Z((u*w-u*v*v-w*t*t)/
(u*v*v+w*t*t)));l=u*e*v/h+(b+k)/2;var u=u*-h*t/e+(d+p)/2,v=F.asin(((d-u)/h).toFixed(9));t=F.asin(((p-u)/h).toFixed(9));v=b<l?S-v:v;t=k<l?S-t:t;0>v&&(v=2*S+v);0>t&&(t=2*S+t);n&&v>t&&(v-=2*S);!n&&t>v&&(t-=2*S)}if(Z(t-v)>r){var c=t,w=k,G=p;t=v+r*(n&&t>v?1:-1);k=l+e*F.cos(t);p=u+h*F.sin(t);c=K(k,p,e,h,f,0,n,w,G,[t,c,l,u])}l=t-v;f=F.cos(v);r=F.sin(v);n=F.cos(t);t=F.sin(t);l=F.tan(l/4);e=4/3*e*l;l*=4/3*h;h=[b,d];b=[b+e*r,d-l*f];d=[k+e*t,p-l*n];k=[k,p];b[0]=2*h[0]-b[0];b[1]=2*h[1]-b[1];if(s)return[b,d,k].concat(c);
c=[b,d,k].concat(c).join().split(",");s=[];k=0;for(p=c.length;k<p;k++)s[k]=k%2?x(c[k-1],c[k],q).y:x(c[k],c[k+1],q).x;return s}function U(a,b,d,e,h,f,l,k){for(var n=[],p=[[],[] ],s,r,c,t,q=0;2>q;++q)0==q?(r=6*a-12*d+6*h,s=-3*a+9*d-9*h+3*l,c=3*d-3*a):(r=6*b-12*e+6*f,s=-3*b+9*e-9*f+3*k,c=3*e-3*b),1E-12>Z(s)?1E-12>Z(r)||(s=-c/r,0<s&&1>s&&n.push(s)):(t=r*r-4*c*s,c=F.sqrt(t),0>t||(t=(-r+c)/(2*s),0<t&&1>t&&n.push(t),s=(-r-c)/(2*s),0<s&&1>s&&n.push(s)));for(r=q=n.length;q--;)s=n[q],c=1-s,p[0][q]=c*c*c*a+3*
c*c*s*d+3*c*s*s*h+s*s*s*l,p[1][q]=c*c*c*b+3*c*c*s*e+3*c*s*s*f+s*s*s*k;p[0][r]=a;p[1][r]=b;p[0][r+1]=l;p[1][r+1]=k;p[0].length=p[1].length=r+2;return{min:{x:X.apply(0,p[0]),y:X.apply(0,p[1])},max:{x:W.apply(0,p[0]),y:W.apply(0,p[1])}}}function I(a,b){var e=!b&&A(a);if(!b&&e.curve)return d(e.curve);var f=G(a),l=b&&G(b),n={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},k={x:0,y:0,bx:0,by:0,X:0,Y:0,qx:null,qy:null},p=function(a,b,c){if(!a)return["C",b.x,b.y,b.x,b.y,b.x,b.y];a[0]in{T:1,Q:1}||(b.qx=b.qy=null);
switch(a[0]){case "M":b.X=a[1];b.Y=a[2];break;case "A":a=["C"].concat(K.apply(0,[b.x,b.y].concat(a.slice(1))));break;case "S":"C"==c||"S"==c?(c=2*b.x-b.bx,b=2*b.y-b.by):(c=b.x,b=b.y);a=["C",c,b].concat(a.slice(1));break;case "T":"Q"==c||"T"==c?(b.qx=2*b.x-b.qx,b.qy=2*b.y-b.qy):(b.qx=b.x,b.qy=b.y);a=["C"].concat(J(b.x,b.y,b.qx,b.qy,a[1],a[2]));break;case "Q":b.qx=a[1];b.qy=a[2];a=["C"].concat(J(b.x,b.y,a[1],a[2],a[3],a[4]));break;case "L":a=["C"].concat(h(b.x,b.y,a[1],a[2]));break;case "H":a=["C"].concat(h(b.x,
b.y,a[1],b.y));break;case "V":a=["C"].concat(h(b.x,b.y,b.x,a[1]));break;case "Z":a=["C"].concat(h(b.x,b.y,b.X,b.Y))}return a},s=function(a,b){if(7<a[b].length){a[b].shift();for(var c=a[b];c.length;)q[b]="A",l&&(u[b]="A"),a.splice(b++,0,["C"].concat(c.splice(0,6)));a.splice(b,1);v=W(f.length,l&&l.length||0)}},r=function(a,b,c,d,e){a&&b&&"M"==a[e][0]&&"M"!=b[e][0]&&(b.splice(e,0,["M",d.x,d.y]),c.bx=0,c.by=0,c.x=a[e][1],c.y=a[e][2],v=W(f.length,l&&l.length||0))},q=[],u=[],c="",t="",x=0,v=W(f.length,
l&&l.length||0);for(;x<v;x++){f[x]&&(c=f[x][0]);"C"!=c&&(q[x]=c,x&&(t=q[x-1]));f[x]=p(f[x],n,t);"A"!=q[x]&&"C"==c&&(q[x]="C");s(f,x);l&&(l[x]&&(c=l[x][0]),"C"!=c&&(u[x]=c,x&&(t=u[x-1])),l[x]=p(l[x],k,t),"A"!=u[x]&&"C"==c&&(u[x]="C"),s(l,x));r(f,l,n,k,x);r(l,f,k,n,x);var w=f[x],z=l&&l[x],y=w.length,U=l&&z.length;n.x=w[y-2];n.y=w[y-1];n.bx=$(w[y-4])||n.x;n.by=$(w[y-3])||n.y;k.bx=l&&($(z[U-4])||k.x);k.by=l&&($(z[U-3])||k.y);k.x=l&&z[U-2];k.y=l&&z[U-1]}l||(e.curve=d(f));return l?[f,l]:f}function P(a,
b){for(var d=[],e=0,h=a.length;h-2*!b>e;e+=2){var f=[{x:+a[e-2],y:+a[e-1]},{x:+a[e],y:+a[e+1]},{x:+a[e+2],y:+a[e+3]},{x:+a[e+4],y:+a[e+5]}];b?e?h-4==e?f[3]={x:+a[0],y:+a[1]}:h-2==e&&(f[2]={x:+a[0],y:+a[1]},f[3]={x:+a[2],y:+a[3]}):f[0]={x:+a[h-2],y:+a[h-1]}:h-4==e?f[3]=f[2]:e||(f[0]={x:+a[e],y:+a[e+1]});d.push(["C",(-f[0].x+6*f[1].x+f[2].x)/6,(-f[0].y+6*f[1].y+f[2].y)/6,(f[1].x+6*f[2].x-f[3].x)/6,(f[1].y+6*f[2].y-f[3].y)/6,f[2].x,f[2].y])}return d}y=k.prototype;var Q=a.is,C=a._.clone,L="hasOwnProperty",
N=/,?([a-z]),?/gi,$=parseFloat,F=Math,S=F.PI,X=F.min,W=F.max,ma=F.pow,Z=F.abs;M=n(1);var na=n(),ba=n(0,1),V=a._unit2px;a.path=A;a.path.getTotalLength=M;a.path.getPointAtLength=na;a.path.getSubpath=function(a,b,d){if(1E-6>this.getTotalLength(a)-d)return ba(a,b).end;a=ba(a,d,1);return b?ba(a,b).end:a};y.getTotalLength=function(){if(this.node.getTotalLength)return this.node.getTotalLength()};y.getPointAtLength=function(a){return na(this.attr("d"),a)};y.getSubpath=function(b,d){return a.path.getSubpath(this.attr("d"),
b,d)};a._.box=w;a.path.findDotsAtSegment=u;a.path.bezierBBox=p;a.path.isPointInsideBBox=b;a.path.isBBoxIntersect=q;a.path.intersection=function(a,b){return l(a,b)};a.path.intersectionNumber=function(a,b){return l(a,b,1)};a.path.isPointInside=function(a,d,e){var h=r(a);return b(h,d,e)&&1==l(a,[["M",d,e],["H",h.x2+10] ],1)%2};a.path.getBBox=r;a.path.get={path:function(a){return a.attr("path")},circle:function(a){a=V(a);return x(a.cx,a.cy,a.r)},ellipse:function(a){a=V(a);return x(a.cx||0,a.cy||0,a.rx,
a.ry)},rect:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height,a.rx,a.ry)},image:function(a){a=V(a);return s(a.x||0,a.y||0,a.width,a.height)},line:function(a){return"M"+[a.attr("x1")||0,a.attr("y1")||0,a.attr("x2"),a.attr("y2")]},polyline:function(a){return"M"+a.attr("points")},polygon:function(a){return"M"+a.attr("points")+"z"},deflt:function(a){a=a.node.getBBox();return s(a.x,a.y,a.width,a.height)}};a.path.toRelative=function(b){var e=A(b),h=String.prototype.toLowerCase;if(e.rel)return d(e.rel);
a.is(b,"array")&&a.is(b&&b[0],"array")||(b=a.parsePathString(b));var f=[],l=0,n=0,k=0,p=0,s=0;"M"==b[0][0]&&(l=b[0][1],n=b[0][2],k=l,p=n,s++,f.push(["M",l,n]));for(var r=b.length;s<r;s++){var q=f[s]=[],x=b[s];if(x[0]!=h.call(x[0]))switch(q[0]=h.call(x[0]),q[0]){case "a":q[1]=x[1];q[2]=x[2];q[3]=x[3];q[4]=x[4];q[5]=x[5];q[6]=+(x[6]-l).toFixed(3);q[7]=+(x[7]-n).toFixed(3);break;case "v":q[1]=+(x[1]-n).toFixed(3);break;case "m":k=x[1],p=x[2];default:for(var c=1,t=x.length;c<t;c++)q[c]=+(x[c]-(c%2?l:
n)).toFixed(3)}else for(f[s]=[],"m"==x[0]&&(k=x[1]+l,p=x[2]+n),q=0,c=x.length;q<c;q++)f[s][q]=x[q];x=f[s].length;switch(f[s][0]){case "z":l=k;n=p;break;case "h":l+=+f[s][x-1];break;case "v":n+=+f[s][x-1];break;default:l+=+f[s][x-2],n+=+f[s][x-1]}}f.toString=z;e.rel=d(f);return f};a.path.toAbsolute=G;a.path.toCubic=I;a.path.map=function(a,b){if(!b)return a;var d,e,h,f,l,n,k;a=I(a);h=0;for(l=a.length;h<l;h++)for(k=a[h],f=1,n=k.length;f<n;f+=2)d=b.x(k[f],k[f+1]),e=b.y(k[f],k[f+1]),k[f]=d,k[f+1]=e;return a};
a.path.toString=z;a.path.clone=d});C.plugin(function(a,v,y,C){var A=Math.max,w=Math.min,z=function(a){this.items=[];this.bindings={};this.length=0;this.type="set";if(a)for(var f=0,n=a.length;f<n;f++)a[f]&&(this[this.items.length]=this.items[this.items.length]=a[f],this.length++)};v=z.prototype;v.push=function(){for(var a,f,n=0,k=arguments.length;n<k;n++)if(a=arguments[n])f=this.items.length,this[f]=this.items[f]=a,this.length++;return this};v.pop=function(){this.length&&delete this[this.length--];
return this.items.pop()};v.forEach=function(a,f){for(var n=0,k=this.items.length;n<k&&!1!==a.call(f,this.items[n],n);n++);return this};v.animate=function(d,f,n,u){"function"!=typeof n||n.length||(u=n,n=L.linear);d instanceof a._.Animation&&(u=d.callback,n=d.easing,f=n.dur,d=d.attr);var p=arguments;if(a.is(d,"array")&&a.is(p[p.length-1],"array"))var b=!0;var q,e=function(){q?this.b=q:q=this.b},l=0,r=u&&function(){l++==this.length&&u.call(this)};return this.forEach(function(a,l){k.once("snap.animcreated."+
a.id,e);b?p[l]&&a.animate.apply(a,p[l]):a.animate(d,f,n,r)})};v.remove=function(){for(;this.length;)this.pop().remove();return this};v.bind=function(a,f,k){var u={};if("function"==typeof f)this.bindings[a]=f;else{var p=k||a;this.bindings[a]=function(a){u[p]=a;f.attr(u)}}return this};v.attr=function(a){var f={},k;for(k in a)if(this.bindings[k])this.bindings[k](a[k]);else f[k]=a[k];a=0;for(k=this.items.length;a<k;a++)this.items[a].attr(f);return this};v.clear=function(){for(;this.length;)this.pop()};
v.splice=function(a,f,k){a=0>a?A(this.length+a,0):a;f=A(0,w(this.length-a,f));var u=[],p=[],b=[],q;for(q=2;q<arguments.length;q++)b.push(arguments[q]);for(q=0;q<f;q++)p.push(this[a+q]);for(;q<this.length-a;q++)u.push(this[a+q]);var e=b.length;for(q=0;q<e+u.length;q++)this.items[a+q]=this[a+q]=q<e?b[q]:u[q-e];for(q=this.items.length=this.length-=f-e;this[q];)delete this[q++];return new z(p)};v.exclude=function(a){for(var f=0,k=this.length;f<k;f++)if(this[f]==a)return this.splice(f,1),!0;return!1};
v.insertAfter=function(a){for(var f=this.items.length;f--;)this.items[f].insertAfter(a);return this};v.getBBox=function(){for(var a=[],f=[],k=[],u=[],p=this.items.length;p--;)if(!this.items[p].removed){var b=this.items[p].getBBox();a.push(b.x);f.push(b.y);k.push(b.x+b.width);u.push(b.y+b.height)}a=w.apply(0,a);f=w.apply(0,f);k=A.apply(0,k);u=A.apply(0,u);return{x:a,y:f,x2:k,y2:u,width:k-a,height:u-f,cx:a+(k-a)/2,cy:f+(u-f)/2}};v.clone=function(a){a=new z;for(var f=0,k=this.items.length;f<k;f++)a.push(this.items[f].clone());
return a};v.toString=function(){return"Snap\u2018s set"};v.type="set";a.set=function(){var a=new z;arguments.length&&a.push.apply(a,Array.prototype.slice.call(arguments,0));return a}});C.plugin(function(a,v,y,C){function A(a){var b=a[0];switch(b.toLowerCase()){case "t":return[b,0,0];case "m":return[b,1,0,0,1,0,0];case "r":return 4==a.length?[b,0,a[2],a[3] ]:[b,0];case "s":return 5==a.length?[b,1,1,a[3],a[4] ]:3==a.length?[b,1,1]:[b,1]}}function w(b,d,f){d=q(d).replace(/\.{3}|\u2026/g,b);b=a.parseTransformString(b)||
[];d=a.parseTransformString(d)||[];for(var k=Math.max(b.length,d.length),p=[],v=[],h=0,w,z,y,I;h<k;h++){y=b[h]||A(d[h]);I=d[h]||A(y);if(y[0]!=I[0]||"r"==y[0].toLowerCase()&&(y[2]!=I[2]||y[3]!=I[3])||"s"==y[0].toLowerCase()&&(y[3]!=I[3]||y[4]!=I[4])){b=a._.transform2matrix(b,f());d=a._.transform2matrix(d,f());p=[["m",b.a,b.b,b.c,b.d,b.e,b.f] ];v=[["m",d.a,d.b,d.c,d.d,d.e,d.f] ];break}p[h]=[];v[h]=[];w=0;for(z=Math.max(y.length,I.length);w<z;w++)w in y&&(p[h][w]=y[w]),w in I&&(v[h][w]=I[w])}return{from:u(p),
to:u(v),f:n(p)}}function z(a){return a}function d(a){return function(b){return+b.toFixed(3)+a}}function f(b){return a.rgb(b[0],b[1],b[2])}function n(a){var b=0,d,f,k,n,h,p,q=[];d=0;for(f=a.length;d<f;d++){h="[";p=['"'+a[d][0]+'"'];k=1;for(n=a[d].length;k<n;k++)p[k]="val["+b++ +"]";h+=p+"]";q[d]=h}return Function("val","return Snap.path.toString.call(["+q+"])")}function u(a){for(var b=[],d=0,f=a.length;d<f;d++)for(var k=1,n=a[d].length;k<n;k++)b.push(a[d][k]);return b}var p={},b=/[a-z]+$/i,q=String;
p.stroke=p.fill="colour";v.prototype.equal=function(a,b){return k("snap.util.equal",this,a,b).firstDefined()};k.on("snap.util.equal",function(e,k){var r,s;r=q(this.attr(e)||"");var x=this;if(r==+r&&k==+k)return{from:+r,to:+k,f:z};if("colour"==p[e])return r=a.color(r),s=a.color(k),{from:[r.r,r.g,r.b,r.opacity],to:[s.r,s.g,s.b,s.opacity],f:f};if("transform"==e||"gradientTransform"==e||"patternTransform"==e)return k instanceof a.Matrix&&(k=k.toTransformString()),a._.rgTransform.test(k)||(k=a._.svgTransform2string(k)),
w(r,k,function(){return x.getBBox(1)});if("d"==e||"path"==e)return r=a.path.toCubic(r,k),{from:u(r[0]),to:u(r[1]),f:n(r[0])};if("points"==e)return r=q(r).split(a._.separator),s=q(k).split(a._.separator),{from:r,to:s,f:function(a){return a}};aUnit=r.match(b);s=q(k).match(b);return aUnit&&aUnit==s?{from:parseFloat(r),to:parseFloat(k),f:d(aUnit)}:{from:this.asPX(e),to:this.asPX(e,k),f:z}})});C.plugin(function(a,v,y,C){var A=v.prototype,w="createTouch"in C.doc;v="click dblclick mousedown mousemove mouseout mouseover mouseup touchstart touchmove touchend touchcancel".split(" ");
var z={mousedown:"touchstart",mousemove:"touchmove",mouseup:"touchend"},d=function(a,b){var d="y"==a?"scrollTop":"scrollLeft",e=b&&b.node?b.node.ownerDocument:C.doc;return e[d in e.documentElement?"documentElement":"body"][d]},f=function(){this.returnValue=!1},n=function(){return this.originalEvent.preventDefault()},u=function(){this.cancelBubble=!0},p=function(){return this.originalEvent.stopPropagation()},b=function(){if(C.doc.addEventListener)return function(a,b,e,f){var k=w&&z[b]?z[b]:b,l=function(k){var l=
d("y",f),q=d("x",f);if(w&&z.hasOwnProperty(b))for(var r=0,u=k.targetTouches&&k.targetTouches.length;r<u;r++)if(k.targetTouches[r].target==a||a.contains(k.targetTouches[r].target)){u=k;k=k.targetTouches[r];k.originalEvent=u;k.preventDefault=n;k.stopPropagation=p;break}return e.call(f,k,k.clientX+q,k.clientY+l)};b!==k&&a.addEventListener(b,l,!1);a.addEventListener(k,l,!1);return function(){b!==k&&a.removeEventListener(b,l,!1);a.removeEventListener(k,l,!1);return!0}};if(C.doc.attachEvent)return function(a,
b,e,h){var k=function(a){a=a||h.node.ownerDocument.window.event;var b=d("y",h),k=d("x",h),k=a.clientX+k,b=a.clientY+b;a.preventDefault=a.preventDefault||f;a.stopPropagation=a.stopPropagation||u;return e.call(h,a,k,b)};a.attachEvent("on"+b,k);return function(){a.detachEvent("on"+b,k);return!0}}}(),q=[],e=function(a){for(var b=a.clientX,e=a.clientY,f=d("y"),l=d("x"),n,p=q.length;p--;){n=q[p];if(w)for(var r=a.touches&&a.touches.length,u;r--;){if(u=a.touches[r],u.identifier==n.el._drag.id||n.el.node.contains(u.target)){b=
u.clientX;e=u.clientY;(a.originalEvent?a.originalEvent:a).preventDefault();break}}else a.preventDefault();b+=l;e+=f;k("snap.drag.move."+n.el.id,n.move_scope||n.el,b-n.el._drag.x,e-n.el._drag.y,b,e,a)}},l=function(b){a.unmousemove(e).unmouseup(l);for(var d=q.length,f;d--;)f=q[d],f.el._drag={},k("snap.drag.end."+f.el.id,f.end_scope||f.start_scope||f.move_scope||f.el,b);q=[]};for(y=v.length;y--;)(function(d){a[d]=A[d]=function(e,f){a.is(e,"function")&&(this.events=this.events||[],this.events.push({name:d,
f:e,unbind:b(this.node||document,d,e,f||this)}));return this};a["un"+d]=A["un"+d]=function(a){for(var b=this.events||[],e=b.length;e--;)if(b[e].name==d&&(b[e].f==a||!a)){b[e].unbind();b.splice(e,1);!b.length&&delete this.events;break}return this}})(v[y]);A.hover=function(a,b,d,e){return this.mouseover(a,d).mouseout(b,e||d)};A.unhover=function(a,b){return this.unmouseover(a).unmouseout(b)};var r=[];A.drag=function(b,d,f,h,n,p){function u(r,v,w){(r.originalEvent||r).preventDefault();this._drag.x=v;
this._drag.y=w;this._drag.id=r.identifier;!q.length&&a.mousemove(e).mouseup(l);q.push({el:this,move_scope:h,start_scope:n,end_scope:p});d&&k.on("snap.drag.start."+this.id,d);b&&k.on("snap.drag.move."+this.id,b);f&&k.on("snap.drag.end."+this.id,f);k("snap.drag.start."+this.id,n||h||this,v,w,r)}if(!arguments.length){var v;return this.drag(function(a,b){this.attr({transform:v+(v?"T":"t")+[a,b]})},function(){v=this.transform().local})}this._drag={};r.push({el:this,start:u});this.mousedown(u);return this};
A.undrag=function(){for(var b=r.length;b--;)r[b].el==this&&(this.unmousedown(r[b].start),r.splice(b,1),k.unbind("snap.drag.*."+this.id));!r.length&&a.unmousemove(e).unmouseup(l);return this}});C.plugin(function(a,v,y,C){y=y.prototype;var A=/^\s*url\((.+)\)/,w=String,z=a._.$;a.filter={};y.filter=function(d){var f=this;"svg"!=f.type&&(f=f.paper);d=a.parse(w(d));var k=a._.id(),u=z("filter");z(u,{id:k,filterUnits:"userSpaceOnUse"});u.appendChild(d.node);f.defs.appendChild(u);return new v(u)};k.on("snap.util.getattr.filter",
function(){k.stop();var d=z(this.node,"filter");if(d)return(d=w(d).match(A))&&a.select(d[1])});k.on("snap.util.attr.filter",function(d){if(d instanceof v&&"filter"==d.type){k.stop();var f=d.node.id;f||(z(d.node,{id:d.id}),f=d.id);z(this.node,{filter:a.url(f)})}d&&"none"!=d||(k.stop(),this.node.removeAttribute("filter"))});a.filter.blur=function(d,f){null==d&&(d=2);return a.format('<feGaussianBlur stdDeviation="{def}"/>',{def:null==f?d:[d,f]})};a.filter.blur.toString=function(){return this()};a.filter.shadow=
function(d,f,k,u,p){"string"==typeof k&&(p=u=k,k=4);"string"!=typeof u&&(p=u,u="#000");null==k&&(k=4);null==p&&(p=1);null==d&&(d=0,f=2);null==f&&(f=d);u=a.color(u||"#000");return a.format('<feGaussianBlur in="SourceAlpha" stdDeviation="{blur}"/><feOffset dx="{dx}" dy="{dy}" result="offsetblur"/><feFlood flood-color="{color}"/><feComposite in2="offsetblur" operator="in"/><feComponentTransfer><feFuncA type="linear" slope="{opacity}"/></feComponentTransfer><feMerge><feMergeNode/><feMergeNode in="SourceGraphic"/></feMerge>',
{color:u,dx:d,dy:f,blur:k,opacity:p})};a.filter.shadow.toString=function(){return this()};a.filter.grayscale=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {b} {h} 0 0 0 0 0 1 0"/>',{a:0.2126+0.7874*(1-d),b:0.7152-0.7152*(1-d),c:0.0722-0.0722*(1-d),d:0.2126-0.2126*(1-d),e:0.7152+0.2848*(1-d),f:0.0722-0.0722*(1-d),g:0.2126-0.2126*(1-d),h:0.0722+0.9278*(1-d)})};a.filter.grayscale.toString=function(){return this()};a.filter.sepia=
function(d){null==d&&(d=1);return a.format('<feColorMatrix type="matrix" values="{a} {b} {c} 0 0 {d} {e} {f} 0 0 {g} {h} {i} 0 0 0 0 0 1 0"/>',{a:0.393+0.607*(1-d),b:0.769-0.769*(1-d),c:0.189-0.189*(1-d),d:0.349-0.349*(1-d),e:0.686+0.314*(1-d),f:0.168-0.168*(1-d),g:0.272-0.272*(1-d),h:0.534-0.534*(1-d),i:0.131+0.869*(1-d)})};a.filter.sepia.toString=function(){return this()};a.filter.saturate=function(d){null==d&&(d=1);return a.format('<feColorMatrix type="saturate" values="{amount}"/>',{amount:1-
d})};a.filter.saturate.toString=function(){return this()};a.filter.hueRotate=function(d){return a.format('<feColorMatrix type="hueRotate" values="{angle}"/>',{angle:d||0})};a.filter.hueRotate.toString=function(){return this()};a.filter.invert=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="table" tableValues="{amount} {amount2}"/><feFuncG type="table" tableValues="{amount} {amount2}"/><feFuncB type="table" tableValues="{amount} {amount2}"/></feComponentTransfer>',{amount:d,
amount2:1-d})};a.filter.invert.toString=function(){return this()};a.filter.brightness=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}"/><feFuncG type="linear" slope="{amount}"/><feFuncB type="linear" slope="{amount}"/></feComponentTransfer>',{amount:d})};a.filter.brightness.toString=function(){return this()};a.filter.contrast=function(d){null==d&&(d=1);return a.format('<feComponentTransfer><feFuncR type="linear" slope="{amount}" intercept="{amount2}"/><feFuncG type="linear" slope="{amount}" intercept="{amount2}"/><feFuncB type="linear" slope="{amount}" intercept="{amount2}"/></feComponentTransfer>',
{amount:d,amount2:0.5-d/2})};a.filter.contrast.toString=function(){return this()}});return C});

]]> </script>
<script> <![CDATA[

(function (glob, factory) {
    // AMD support
    if (typeof define === "function" && define.amd) {
        // Define as an anonymous module
        define("Gadfly", ["Snap.svg"], function (Snap) {
            return factory(Snap);
        });
    } else {
        // Browser globals (glob is window)
        // Snap adds itself to window
        glob.Gadfly = factory(glob.Snap);
    }
}(this, function (Snap) {

var Gadfly = {};

// Get an x/y coordinate value in pixels
var xPX = function(fig, x) {
    var client_box = fig.node.getBoundingClientRect();
    return x * fig.node.viewBox.baseVal.width / client_box.width;
};

var yPX = function(fig, y) {
    var client_box = fig.node.getBoundingClientRect();
    return y * fig.node.viewBox.baseVal.height / client_box.height;
};


Snap.plugin(function (Snap, Element, Paper, global) {
    // Traverse upwards from a snap element to find and return the first
    // note with the "plotroot" class.
    Element.prototype.plotroot = function () {
        var element = this;
        while (!element.hasClass("plotroot") && element.parent() != null) {
            element = element.parent();
        }
        return element;
    };

    Element.prototype.svgroot = function () {
        var element = this;
        while (element.node.nodeName != "svg" && element.parent() != null) {
            element = element.parent();
        }
        return element;
    };

    Element.prototype.plotbounds = function () {
        var root = this.plotroot()
        var bbox = root.select(".guide.background").node.getBBox();
        return {
            x0: bbox.x,
            x1: bbox.x + bbox.width,
            y0: bbox.y,
            y1: bbox.y + bbox.height
        };
    };

    Element.prototype.plotcenter = function () {
        var root = this.plotroot()
        var bbox = root.select(".guide.background").node.getBBox();
        return {
            x: bbox.x + bbox.width / 2,
            y: bbox.y + bbox.height / 2
        };
    };

    // Emulate IE style mouseenter/mouseleave events, since Microsoft always
    // does everything right.
    // See: http://www.dynamic-tools.net/toolbox/isMouseLeaveOrEnter/
    var events = ["mouseenter", "mouseleave"];

    for (i in events) {
        (function (event_name) {
            var event_name = events[i];
            Element.prototype[event_name] = function (fn, scope) {
                if (Snap.is(fn, "function")) {
                    var fn2 = function (event) {
                        if (event.type != "mouseover" && event.type != "mouseout") {
                            return;
                        }

                        var reltg = event.relatedTarget ? event.relatedTarget :
                            event.type == "mouseout" ? event.toElement : event.fromElement;
                        while (reltg && reltg != this.node) reltg = reltg.parentNode;

                        if (reltg != this.node) {
                            return fn.apply(this, event);
                        }
                    };

                    if (event_name == "mouseenter") {
                        this.mouseover(fn2, scope);
                    } else {
                        this.mouseout(fn2, scope);
                    }
                }
                return this;
            };
        })(events[i]);
    }


    Element.prototype.mousewheel = function (fn, scope) {
        if (Snap.is(fn, "function")) {
            var el = this;
            var fn2 = function (event) {
                fn.apply(el, [event]);
            };
        }

        this.node.addEventListener(
            /Firefox/i.test(navigator.userAgent) ? "DOMMouseScroll" : "mousewheel",
            fn2);

        return this;
    };


    // Snap's attr function can be too slow for things like panning/zooming.
    // This is a function to directly update element attributes without going
    // through eve.
    Element.prototype.attribute = function(key, val) {
        if (val === undefined) {
            return this.node.getAttribute(key);
        } else {
            this.node.setAttribute(key, val);
            return this;
        }
    };

    Element.prototype.init_gadfly = function() {
        this.mouseenter(Gadfly.plot_mouseover)
            .mouseleave(Gadfly.plot_mouseout)
            .dblclick(Gadfly.plot_dblclick)
            .mousewheel(Gadfly.guide_background_scroll)
            .drag(Gadfly.guide_background_drag_onmove,
                  Gadfly.guide_background_drag_onstart,
                  Gadfly.guide_background_drag_onend);
        this.mouseenter(function (event) {
            init_pan_zoom(this.plotroot());
        });
        return this;
    };
});


// When the plot is moused over, emphasize the grid lines.
Gadfly.plot_mouseover = function(event) {
    var root = this.plotroot();

    var keyboard_zoom = function(event) {
        if (event.which == 187) { // plus
            increase_zoom_by_position(root, 0.1, true);
        } else if (event.which == 189) { // minus
            increase_zoom_by_position(root, -0.1, true);
        }
    };
    root.data("keyboard_zoom", keyboard_zoom);
    window.addEventListener("keyup", keyboard_zoom);

    var xgridlines = root.select(".xgridlines"),
        ygridlines = root.select(".ygridlines");

    xgridlines.data("unfocused_strokedash",
                    xgridlines.attribute("stroke-dasharray").replace(/(\d)(,|$)/g, "$1mm$2"));
    ygridlines.data("unfocused_strokedash",
                    ygridlines.attribute("stroke-dasharray").replace(/(\d)(,|$)/g, "$1mm$2"));

    // emphasize grid lines
    var destcolor = root.data("focused_xgrid_color");
    xgridlines.attribute("stroke-dasharray", "none")
              .selectAll("path")
              .animate({stroke: destcolor}, 250);

    destcolor = root.data("focused_ygrid_color");
    ygridlines.attribute("stroke-dasharray", "none")
              .selectAll("path")
              .animate({stroke: destcolor}, 250);

    // reveal zoom slider
    root.select(".zoomslider")
        .animate({opacity: 1.0}, 250);
};

// Reset pan and zoom on double click
Gadfly.plot_dblclick = function(event) {
  set_plot_pan_zoom(this.plotroot(), 0.0, 0.0, 1.0);
};

// Unemphasize grid lines on mouse out.
Gadfly.plot_mouseout = function(event) {
    var root = this.plotroot();

    window.removeEventListener("keyup", root.data("keyboard_zoom"));
    root.data("keyboard_zoom", undefined);

    var xgridlines = root.select(".xgridlines"),
        ygridlines = root.select(".ygridlines");

    var destcolor = root.data("unfocused_xgrid_color");

    xgridlines.attribute("stroke-dasharray", xgridlines.data("unfocused_strokedash"))
              .selectAll("path")
              .animate({stroke: destcolor}, 250);

    destcolor = root.data("unfocused_ygrid_color");
    ygridlines.attribute("stroke-dasharray", ygridlines.data("unfocused_strokedash"))
              .selectAll("path")
              .animate({stroke: destcolor}, 250);

    // hide zoom slider
    root.select(".zoomslider")
        .animate({opacity: 0.0}, 250);
};


var set_geometry_transform = function(root, tx, ty, scale) {
    var xscalable = root.hasClass("xscalable"),
        yscalable = root.hasClass("yscalable");

    var old_scale = root.data("scale");

    var xscale = xscalable ? scale : 1.0,
        yscale = yscalable ? scale : 1.0;

    tx = xscalable ? tx : 0.0;
    ty = yscalable ? ty : 0.0;

    var t = new Snap.Matrix().translate(tx, ty).scale(xscale, yscale);

    root.selectAll(".geometry, image")
        .forEach(function (element, i) {
            element.transform(t);
        });

    bounds = root.plotbounds();

    if (yscalable) {
        var xfixed_t = new Snap.Matrix().translate(0, ty).scale(1.0, yscale);
        root.selectAll(".xfixed")
            .forEach(function (element, i) {
                element.transform(xfixed_t);
            });

        root.select(".ylabels")
            .transform(xfixed_t)
            .selectAll("text")
            .forEach(function (element, i) {
                if (element.attribute("gadfly:inscale") == "true") {
                    var cx = element.asPX("x"),
                        cy = element.asPX("y");
                    var st = element.data("static_transform");
                    unscale_t = new Snap.Matrix();
                    unscale_t.scale(1, 1/scale, cx, cy).add(st);
                    element.transform(unscale_t);

                    var y = cy * scale + ty;
                    element.attr("visibility",
                        bounds.y0 <= y && y <= bounds.y1 ? "visible" : "hidden");
                }
            });
    }

    if (xscalable) {
        var yfixed_t = new Snap.Matrix().translate(tx, 0).scale(xscale, 1.0);
        var xtrans = new Snap.Matrix().translate(tx, 0);
        root.selectAll(".yfixed")
            .forEach(function (element, i) {
                element.transform(yfixed_t);
            });

        root.select(".xlabels")
            .transform(yfixed_t)
            .selectAll("text")
            .forEach(function (element, i) {
                if (element.attribute("gadfly:inscale") == "true") {
                    var cx = element.asPX("x"),
                        cy = element.asPX("y");
                    var st = element.data("static_transform");
                    unscale_t = new Snap.Matrix();
                    unscale_t.scale(1/scale, 1, cx, cy).add(st);

                    element.transform(unscale_t);

                    var x = cx * scale + tx;
                    element.attr("visibility",
                        bounds.x0 <= x && x <= bounds.x1 ? "visible" : "hidden");
                    }
            });
    }

    // we must unscale anything that is scale invariance: widths, raiduses, etc.
    var size_attribs = ["font-size"];
    var unscaled_selection = ".geometry, .geometry *";
    if (xscalable) {
        size_attribs.push("rx");
        unscaled_selection += ", .xgridlines";
    }
    if (yscalable) {
        size_attribs.push("ry");
        unscaled_selection += ", .ygridlines";
    }

    root.selectAll(unscaled_selection)
        .forEach(function (element, i) {
            // circle need special help
            if (element.node.nodeName == "circle") {
                var cx = element.attribute("cx"),
                    cy = element.attribute("cy");
                unscale_t = new Snap.Matrix().scale(1/xscale, 1/yscale,
                                                        cx, cy);
                element.transform(unscale_t);
                return;
            }

            for (i in size_attribs) {
                var key = size_attribs[i];
                var val = parseFloat(element.attribute(key));
                if (val !== undefined && val != 0 && !isNaN(val)) {
                    element.attribute(key, val * old_scale / scale);
                }
            }
        });
};


// Find the most appropriate tick scale and update label visibility.
var update_tickscale = function(root, scale, axis) {
    if (!root.hasClass(axis + "scalable")) return;

    var tickscales = root.data(axis + "tickscales");
    var best_tickscale = 1.0;
    var best_tickscale_dist = Infinity;
    for (tickscale in tickscales) {
        var dist = Math.abs(Math.log(tickscale) - Math.log(scale));
        if (dist < best_tickscale_dist) {
            best_tickscale_dist = dist;
            best_tickscale = tickscale;
        }
    }

    if (best_tickscale != root.data(axis + "tickscale")) {
        root.data(axis + "tickscale", best_tickscale);
        var mark_inscale_gridlines = function (element, i) {
            var inscale = element.attr("gadfly:scale") == best_tickscale;
            element.attribute("gadfly:inscale", inscale);
            element.attr("visibility", inscale ? "visible" : "hidden");
        };

        var mark_inscale_labels = function (element, i) {
            var inscale = element.attr("gadfly:scale") == best_tickscale;
            element.attribute("gadfly:inscale", inscale);
            element.attr("visibility", inscale ? "visible" : "hidden");
        };

        root.select("." + axis + "gridlines").selectAll("path").forEach(mark_inscale_gridlines);
        root.select("." + axis + "labels").selectAll("text").forEach(mark_inscale_labels);
    }
};


var set_plot_pan_zoom = function(root, tx, ty, scale) {
    var old_scale = root.data("scale");
    var bounds = root.plotbounds();

    var width = bounds.x1 - bounds.x0,
        height = bounds.y1 - bounds.y0;

    // compute the viewport derived from tx, ty, and scale
    var x_min = -width * scale - (scale * width - width),
        x_max = width * scale,
        y_min = -height * scale - (scale * height - height),
        y_max = height * scale;

    var x0 = bounds.x0 - scale * bounds.x0,
        y0 = bounds.y0 - scale * bounds.y0;

    var tx = Math.max(Math.min(tx - x0, x_max), x_min),
        ty = Math.max(Math.min(ty - y0, y_max), y_min);

    tx += x0;
    ty += y0;

    // when the scale change, we may need to alter which set of
    // ticks is being displayed
    if (scale != old_scale) {
        update_tickscale(root, scale, "x");
        update_tickscale(root, scale, "y");
    }

    set_geometry_transform(root, tx, ty, scale);

    root.data("scale", scale);
    root.data("tx", tx);
    root.data("ty", ty);
};


var scale_centered_translation = function(root, scale) {
    var bounds = root.plotbounds();

    var width = bounds.x1 - bounds.x0,
        height = bounds.y1 - bounds.y0;

    var tx0 = root.data("tx"),
        ty0 = root.data("ty");

    var scale0 = root.data("scale");

    // how off from center the current view is
    var xoff = tx0 - (bounds.x0 * (1 - scale0) + (width * (1 - scale0)) / 2),
        yoff = ty0 - (bounds.y0 * (1 - scale0) + (height * (1 - scale0)) / 2);

    // rescale offsets
    xoff = xoff * scale / scale0;
    yoff = yoff * scale / scale0;

    // adjust for the panel position being scaled
    var x_edge_adjust = bounds.x0 * (1 - scale),
        y_edge_adjust = bounds.y0 * (1 - scale);

    return {
        x: xoff + x_edge_adjust + (width - width * scale) / 2,
        y: yoff + y_edge_adjust + (height - height * scale) / 2
    };
};


// Initialize data for panning zooming if it isn't already.
var init_pan_zoom = function(root) {
    if (root.data("zoompan-ready")) {
        return;
    }

    // The non-scaling-stroke trick. Rather than try to correct for the
    // stroke-width when zooming, we force it to a fixed value.
    var px_per_mm = root.node.getCTM().a;

    // Drag events report deltas in pixels, which we'd like to convert to
    // millimeters.
    root.data("px_per_mm", px_per_mm);

    root.selectAll("path")
        .forEach(function (element, i) {
        sw = element.asPX("stroke-width") * px_per_mm;
        if (sw > 0) {
            element.attribute("stroke-width", sw);
            element.attribute("vector-effect", "non-scaling-stroke");
        }
    });

    // Store ticks labels original tranformation
    root.selectAll(".xlabels > text, .ylabels > text")
        .forEach(function (element, i) {
            var lm = element.transform().localMatrix;
            element.data("static_transform",
                new Snap.Matrix(lm.a, lm.b, lm.c, lm.d, lm.e, lm.f));
        });

    var xgridlines = root.select(".xgridlines");
    var ygridlines = root.select(".ygridlines");
    var xlabels = root.select(".xlabels");
    var ylabels = root.select(".ylabels");

    if (root.data("tx") === undefined) root.data("tx", 0);
    if (root.data("ty") === undefined) root.data("ty", 0);
    if (root.data("scale") === undefined) root.data("scale", 1.0);
    if (root.data("xtickscales") === undefined) {

        // index all the tick scales that are listed
        var xtickscales = {};
        var ytickscales = {};
        var add_x_tick_scales = function (element, i) {
            xtickscales[element.attribute("gadfly:scale")] = true;
        };
        var add_y_tick_scales = function (element, i) {
            ytickscales[element.attribute("gadfly:scale")] = true;
        };

        if (xgridlines) xgridlines.selectAll("path").forEach(add_x_tick_scales);
        if (ygridlines) ygridlines.selectAll("path").forEach(add_y_tick_scales);
        if (xlabels) xlabels.selectAll("text").forEach(add_x_tick_scales);
        if (ylabels) ylabels.selectAll("text").forEach(add_y_tick_scales);

        root.data("xtickscales", xtickscales);
        root.data("ytickscales", ytickscales);
        root.data("xtickscale", 1.0);
    }

    var min_scale = 1.0, max_scale = 1.0;
    for (scale in xtickscales) {
        min_scale = Math.min(min_scale, scale);
        max_scale = Math.max(max_scale, scale);
    }
    for (scale in ytickscales) {
        min_scale = Math.min(min_scale, scale);
        max_scale = Math.max(max_scale, scale);
    }
    root.data("min_scale", min_scale);
    root.data("max_scale", max_scale);

    // store the original positions of labels
    if (xlabels) {
        xlabels.selectAll("text")
               .forEach(function (element, i) {
                   element.data("x", element.asPX("x"));
               });
    }

    if (ylabels) {
        ylabels.selectAll("text")
               .forEach(function (element, i) {
                   element.data("y", element.asPX("y"));
               });
    }

    // mark grid lines and ticks as in or out of scale.
    var mark_inscale = function (element, i) {
        element.attribute("gadfly:inscale", element.attribute("gadfly:scale") == 1.0);
    };

    if (xgridlines) xgridlines.selectAll("path").forEach(mark_inscale);
    if (ygridlines) ygridlines.selectAll("path").forEach(mark_inscale);
    if (xlabels) xlabels.selectAll("text").forEach(mark_inscale);
    if (ylabels) ylabels.selectAll("text").forEach(mark_inscale);

    // figure out the upper ond lower bounds on panning using the maximum
    // and minum grid lines
    var bounds = root.plotbounds();
    var pan_bounds = {
        x0: 0.0,
        y0: 0.0,
        x1: 0.0,
        y1: 0.0
    };

    if (xgridlines) {
        xgridlines
            .selectAll("path")
            .forEach(function (element, i) {
                if (element.attribute("gadfly:inscale") == "true") {
                    var bbox = element.node.getBBox();
                    if (bounds.x1 - bbox.x < pan_bounds.x0) {
                        pan_bounds.x0 = bounds.x1 - bbox.x;
                    }
                    if (bounds.x0 - bbox.x > pan_bounds.x1) {
                        pan_bounds.x1 = bounds.x0 - bbox.x;
                    }
                    element.attr("visibility", "visible");
                }
            });
    }

    if (ygridlines) {
        ygridlines
            .selectAll("path")
            .forEach(function (element, i) {
                if (element.attribute("gadfly:inscale") == "true") {
                    var bbox = element.node.getBBox();
                    if (bounds.y1 - bbox.y < pan_bounds.y0) {
                        pan_bounds.y0 = bounds.y1 - bbox.y;
                    }
                    if (bounds.y0 - bbox.y > pan_bounds.y1) {
                        pan_bounds.y1 = bounds.y0 - bbox.y;
                    }
                    element.attr("visibility", "visible");
                }
            });
    }

    // nudge these values a little
    pan_bounds.x0 -= 5;
    pan_bounds.x1 += 5;
    pan_bounds.y0 -= 5;
    pan_bounds.y1 += 5;
    root.data("pan_bounds", pan_bounds);

    root.data("zoompan-ready", true)
};


// drag actions, i.e. zooming and panning
var pan_action = {
    start: function(root, x, y, event) {
        root.data("dx", 0);
        root.data("dy", 0);
        root.data("tx0", root.data("tx"));
        root.data("ty0", root.data("ty"));
    },
    update: function(root, dx, dy, x, y, event) {
        var px_per_mm = root.data("px_per_mm");
        dx /= px_per_mm;
        dy /= px_per_mm;

        var tx0 = root.data("tx"),
            ty0 = root.data("ty");

        var dx0 = root.data("dx"),
            dy0 = root.data("dy");

        root.data("dx", dx);
        root.data("dy", dy);

        dx = dx - dx0;
        dy = dy - dy0;

        var tx = tx0 + dx,
            ty = ty0 + dy;

        set_plot_pan_zoom(root, tx, ty, root.data("scale"));
    },
    end: function(root, event) {

    },
    cancel: function(root) {
        set_plot_pan_zoom(root, root.data("tx0"), root.data("ty0"), root.data("scale"));
    }
};

var zoom_box;
var zoom_action = {
    start: function(root, x, y, event) {
        var bounds = root.plotbounds();
        var width = bounds.x1 - bounds.x0,
            height = bounds.y1 - bounds.y0;
        var ratio = width / height;
        var xscalable = root.hasClass("xscalable"),
            yscalable = root.hasClass("yscalable");
        var px_per_mm = root.data("px_per_mm");
        x = xscalable ? x / px_per_mm : bounds.x0;
        y = yscalable ? y / px_per_mm : bounds.y0;
        var w = xscalable ? 0 : width;
        var h = yscalable ? 0 : height;
        zoom_box = root.rect(x, y, w, h).attr({
            "fill": "#000",
            "opacity": 0.25
        });
        zoom_box.data("ratio", ratio);
    },
    update: function(root, dx, dy, x, y, event) {
        var xscalable = root.hasClass("xscalable"),
            yscalable = root.hasClass("yscalable");
        var px_per_mm = root.data("px_per_mm");
        var bounds = root.plotbounds();
        if (yscalable) {
            y /= px_per_mm;
            y = Math.max(bounds.y0, y);
            y = Math.min(bounds.y1, y);
        } else {
            y = bounds.y1;
        }
        if (xscalable) {
            x /= px_per_mm;
            x = Math.max(bounds.x0, x);
            x = Math.min(bounds.x1, x);
        } else {
            x = bounds.x1;
        }

        dx = x - zoom_box.attr("x");
        dy = y - zoom_box.attr("y");
        if (xscalable && yscalable) {
            var ratio = zoom_box.data("ratio");
            var width = Math.min(Math.abs(dx), ratio * Math.abs(dy));
            var height = Math.min(Math.abs(dy), Math.abs(dx) / ratio);
            dx = width * dx / Math.abs(dx);
            dy = height * dy / Math.abs(dy);
        }
        var xoffset = 0,
            yoffset = 0;
        if (dx < 0) {
            xoffset = dx;
            dx = -1 * dx;
        }
        if (dy < 0) {
            yoffset = dy;
            dy = -1 * dy;
        }
        if (isNaN(dy)) {
            dy = 0.0;
        }
        if (isNaN(dx)) {
            dx = 0.0;
        }
        zoom_box.transform("T" + xoffset + "," + yoffset);
        zoom_box.attr("width", dx);
        zoom_box.attr("height", dy);
    },
    end: function(root, event) {
        var xscalable = root.hasClass("xscalable"),
            yscalable = root.hasClass("yscalable");
        var zoom_bounds = zoom_box.getBBox();
        if (zoom_bounds.width * zoom_bounds.height <= 0) {
            return;
        }
        var plot_bounds = root.plotbounds();
        var zoom_factor = 1.0;
        if (yscalable) {
            zoom_factor = (plot_bounds.y1 - plot_bounds.y0) / zoom_bounds.height;
        } else {
            zoom_factor = (plot_bounds.x1 - plot_bounds.x0) / zoom_bounds.width;
        }
        var tx = (root.data("tx") - zoom_bounds.x) * zoom_factor + plot_bounds.x0,
            ty = (root.data("ty") - zoom_bounds.y) * zoom_factor + plot_bounds.y0;
        set_plot_pan_zoom(root, tx, ty, root.data("scale") * zoom_factor);
        zoom_box.remove();
    },
    cancel: function(root) {
        zoom_box.remove();
    }
};


Gadfly.guide_background_drag_onstart = function(x, y, event) {
    var root = this.plotroot();
    var scalable = root.hasClass("xscalable") || root.hasClass("yscalable");
    var zoomable = !event.altKey && !event.ctrlKey && event.shiftKey && scalable;
    var panable = !event.altKey && !event.ctrlKey && !event.shiftKey && scalable;
    var drag_action = zoomable ? zoom_action :
                      panable  ? pan_action :
                                 undefined;
    root.data("drag_action", drag_action);
    if (drag_action) {
        var cancel_drag_action = function(event) {
            if (event.which == 27) { // esc key
                drag_action.cancel(root);
                root.data("drag_action", undefined);
            }
        };
        window.addEventListener("keyup", cancel_drag_action);
        root.data("cancel_drag_action", cancel_drag_action);
        drag_action.start(root, x, y, event);
    }
};


Gadfly.guide_background_drag_onmove = function(dx, dy, x, y, event) {
    var root = this.plotroot();
    var drag_action = root.data("drag_action");
    if (drag_action) {
        drag_action.update(root, dx, dy, x, y, event);
    }
};


Gadfly.guide_background_drag_onend = function(event) {
    var root = this.plotroot();
    window.removeEventListener("keyup", root.data("cancel_drag_action"));
    root.data("cancel_drag_action", undefined);
    var drag_action = root.data("drag_action");
    if (drag_action) {
        drag_action.end(root, event);
    }
    root.data("drag_action", undefined);
};


Gadfly.guide_background_scroll = function(event) {
    if (event.shiftKey) {
        increase_zoom_by_position(this.plotroot(), 0.001 * event.wheelDelta);
        event.preventDefault();
    }
};


Gadfly.zoomslider_button_mouseover = function(event) {
    this.select(".button_logo")
         .animate({fill: this.data("mouseover_color")}, 100);
};


Gadfly.zoomslider_button_mouseout = function(event) {
     this.select(".button_logo")
         .animate({fill: this.data("mouseout_color")}, 100);
};


Gadfly.zoomslider_zoomout_click = function(event) {
    increase_zoom_by_position(this.plotroot(), -0.1, true);
};


Gadfly.zoomslider_zoomin_click = function(event) {
    increase_zoom_by_position(this.plotroot(), 0.1, true);
};


Gadfly.zoomslider_track_click = function(event) {
    // TODO
};


// Map slider position x to scale y using the function y = a*exp(b*x)+c.
// The constants a, b, and c are solved using the constraint that the function
// should go through the points (0; min_scale), (0.5; 1), and (1; max_scale).
var scale_from_slider_position = function(position, min_scale, max_scale) {
    var a = (1 - 2 * min_scale + min_scale * min_scale) / (min_scale + max_scale - 2),
        b = 2 * Math.log((max_scale - 1) / (1 - min_scale)),
        c = (min_scale * max_scale - 1) / (min_scale + max_scale - 2);
    return a * Math.exp(b * position) + c;
}

// inverse of scale_from_slider_position
var slider_position_from_scale = function(scale, min_scale, max_scale) {
    var a = (1 - 2 * min_scale + min_scale * min_scale) / (min_scale + max_scale - 2),
        b = 2 * Math.log((max_scale - 1) / (1 - min_scale)),
        c = (min_scale * max_scale - 1) / (min_scale + max_scale - 2);
    return 1 / b * Math.log((scale - c) / a);
}

var increase_zoom_by_position = function(root, delta_position, animate) {
    var scale = root.data("scale"),
        min_scale = root.data("min_scale"),
        max_scale = root.data("max_scale");
    var position = slider_position_from_scale(scale, min_scale, max_scale);
    position += delta_position;
    scale = scale_from_slider_position(position, min_scale, max_scale);
    set_zoom(root, scale, animate);
}

var set_zoom = function(root, scale, animate) {
    var min_scale = root.data("min_scale"),
        max_scale = root.data("max_scale"),
        old_scale = root.data("scale");
    var new_scale = Math.max(min_scale, Math.min(scale, max_scale));
    if (animate) {
        Snap.animate(
            old_scale,
            new_scale,
            function (new_scale) {
                update_plot_scale(root, new_scale);
            },
            200);
    } else {
        update_plot_scale(root, new_scale);
    }
}


var update_plot_scale = function(root, new_scale) {
    var trans = scale_centered_translation(root, new_scale);
    set_plot_pan_zoom(root, trans.x, trans.y, new_scale);

    root.selectAll(".zoomslider_thumb")
        .forEach(function (element, i) {
            var min_pos = element.data("min_pos"),
                max_pos = element.data("max_pos"),
                min_scale = root.data("min_scale"),
                max_scale = root.data("max_scale");
            var xmid = (min_pos + max_pos) / 2;
            var xpos = slider_position_from_scale(new_scale, min_scale, max_scale);
            element.transform(new Snap.Matrix().translate(
                Math.max(min_pos, Math.min(
                         max_pos, min_pos + (max_pos - min_pos) * xpos)) - xmid, 0));
    });
};


Gadfly.zoomslider_thumb_dragmove = function(dx, dy, x, y, event) {
    var root = this.plotroot();
    var min_pos = this.data("min_pos"),
        max_pos = this.data("max_pos"),
        min_scale = root.data("min_scale"),
        max_scale = root.data("max_scale"),
        old_scale = root.data("old_scale");

    var px_per_mm = root.data("px_per_mm");
    dx /= px_per_mm;
    dy /= px_per_mm;

    var xmid = (min_pos + max_pos) / 2;
    var xpos = slider_position_from_scale(old_scale, min_scale, max_scale) +
                   dx / (max_pos - min_pos);

    // compute the new scale
    var new_scale = scale_from_slider_position(xpos, min_scale, max_scale);
    new_scale = Math.min(max_scale, Math.max(min_scale, new_scale));

    update_plot_scale(root, new_scale);
    event.stopPropagation();
};


Gadfly.zoomslider_thumb_dragstart = function(x, y, event) {
    this.animate({fill: this.data("mouseover_color")}, 100);
    var root = this.plotroot();

    // keep track of what the scale was when we started dragging
    root.data("old_scale", root.data("scale"));
    event.stopPropagation();
};


Gadfly.zoomslider_thumb_dragend = function(event) {
    this.animate({fill: this.data("mouseout_color")}, 100);
    event.stopPropagation();
};


var toggle_color_class = function(root, color_class, ison) {
    var guides = root.selectAll(".guide." + color_class + ",.guide ." + color_class);
    var geoms = root.selectAll(".geometry." + color_class + ",.geometry ." + color_class);
    if (ison) {
        guides.animate({opacity: 0.5}, 250);
        geoms.animate({opacity: 0.0}, 250);
    } else {
        guides.animate({opacity: 1.0}, 250);
        geoms.animate({opacity: 1.0}, 250);
    }
};


Gadfly.colorkey_swatch_click = function(event) {
    var root = this.plotroot();
    var color_class = this.data("color_class");

    if (event.shiftKey) {
        root.selectAll(".colorkey text")
            .forEach(function (element) {
                var other_color_class = element.data("color_class");
                if (other_color_class != color_class) {
                    toggle_color_class(root, other_color_class,
                                       element.attr("opacity") == 1.0);
                }
            });
    } else {
        toggle_color_class(root, color_class, this.attr("opacity") == 1.0);
    }
};


return Gadfly;

}));


//@ sourceURL=gadfly.js

(function (glob, factory) {
    // AMD support
      if (typeof require === "function" && typeof define === "function" && define.amd) {
        require(["Snap.svg", "Gadfly"], function (Snap, Gadfly) {
            factory(Snap, Gadfly);
        });
      } else {
          factory(glob.Snap, glob.Gadfly);
      }
})(window, function (Snap, Gadfly) {
    var fig = Snap("#img-8311e7c5");
fig.select("#img-8311e7c5-4")
   .drag(function() {}, function() {}, function() {});
fig.select("#img-8311e7c5-6")
   .data("color_class", "color_S_V")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-7")
   .data("color_class", "color_E_V")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-8")
   .data("color_class", "color_I_V")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-10")
   .data("color_class", "color_S_V")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-11")
   .data("color_class", "color_E_V")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-12")
   .data("color_class", "color_I_V")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-16")
   .init_gadfly();
fig.select("#img-8311e7c5-19")
   .plotroot().data("unfocused_ygrid_color", "#D0D0E0")
;
fig.select("#img-8311e7c5-19")
   .plotroot().data("focused_ygrid_color", "#A0A0A0")
;
fig.select("#img-8311e7c5-125")
   .plotroot().data("unfocused_xgrid_color", "#D0D0E0")
;
fig.select("#img-8311e7c5-125")
   .plotroot().data("focused_xgrid_color", "#A0A0A0")
;
fig.select("#img-8311e7c5-237")
   .data("mouseover_color", "#CD5C5C")
;
fig.select("#img-8311e7c5-237")
   .data("mouseout_color", "#6A6A6A")
;
fig.select("#img-8311e7c5-237")
   .click(Gadfly.zoomslider_zoomin_click)
.mouseenter(Gadfly.zoomslider_button_mouseover)
.mouseleave(Gadfly.zoomslider_button_mouseout)
;
fig.select("#img-8311e7c5-241")
   .data("max_pos", 107.14)
;
fig.select("#img-8311e7c5-241")
   .data("min_pos", 90.14)
;
fig.select("#img-8311e7c5-241")
   .click(Gadfly.zoomslider_track_click);
fig.select("#img-8311e7c5-243")
   .data("max_pos", 107.14)
;
fig.select("#img-8311e7c5-243")
   .data("min_pos", 90.14)
;
fig.select("#img-8311e7c5-243")
   .data("mouseover_color", "#CD5C5C")
;
fig.select("#img-8311e7c5-243")
   .data("mouseout_color", "#6A6A6A")
;
fig.select("#img-8311e7c5-243")
   .drag(Gadfly.zoomslider_thumb_dragmove,
     Gadfly.zoomslider_thumb_dragstart,
     Gadfly.zoomslider_thumb_dragend)
;
fig.select("#img-8311e7c5-245")
   .data("mouseover_color", "#CD5C5C")
;
fig.select("#img-8311e7c5-245")
   .data("mouseout_color", "#6A6A6A")
;
fig.select("#img-8311e7c5-245")
   .click(Gadfly.zoomslider_zoomout_click)
.mouseenter(Gadfly.zoomslider_button_mouseover)
.mouseleave(Gadfly.zoomslider_button_mouseout)
;
fig.select("#img-8311e7c5-466")
   .drag(function() {}, function() {}, function() {});
fig.select("#img-8311e7c5-468")
   .data("color_class", "color_S_H")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-469")
   .data("color_class", "color_E_H")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-470")
   .data("color_class", "color_I_H")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-471")
   .data("color_class", "color_R_H")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-473")
   .data("color_class", "color_S_H")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-474")
   .data("color_class", "color_E_H")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-475")
   .data("color_class", "color_I_H")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-476")
   .data("color_class", "color_R_H")
.click(Gadfly.colorkey_swatch_click)
;
fig.select("#img-8311e7c5-480")
   .init_gadfly();
fig.select("#img-8311e7c5-483")
   .plotroot().data("unfocused_ygrid_color", "#D0D0E0")
;
fig.select("#img-8311e7c5-483")
   .plotroot().data("focused_ygrid_color", "#A0A0A0")
;
fig.select("#img-8311e7c5-637")
   .plotroot().data("unfocused_xgrid_color", "#D0D0E0")
;
fig.select("#img-8311e7c5-637")
   .plotroot().data("focused_xgrid_color", "#A0A0A0")
;
fig.select("#img-8311e7c5-750")
   .data("mouseover_color", "#CD5C5C")
;
fig.select("#img-8311e7c5-750")
   .data("mouseout_color", "#6A6A6A")
;
fig.select("#img-8311e7c5-750")
   .click(Gadfly.zoomslider_zoomin_click)
.mouseenter(Gadfly.zoomslider_button_mouseover)
.mouseleave(Gadfly.zoomslider_button_mouseout)
;
fig.select("#img-8311e7c5-754")
   .data("max_pos", 107.14)
;
fig.select("#img-8311e7c5-754")
   .data("min_pos", 90.14)
;
fig.select("#img-8311e7c5-754")
   .click(Gadfly.zoomslider_track_click);
fig.select("#img-8311e7c5-756")
   .data("max_pos", 107.14)
;
fig.select("#img-8311e7c5-756")
   .data("min_pos", 90.14)
;
fig.select("#img-8311e7c5-756")
   .data("mouseover_color", "#CD5C5C")
;
fig.select("#img-8311e7c5-756")
   .data("mouseout_color", "#6A6A6A")
;
fig.select("#img-8311e7c5-756")
   .drag(Gadfly.zoomslider_thumb_dragmove,
     Gadfly.zoomslider_thumb_dragstart,
     Gadfly.zoomslider_thumb_dragend)
;
fig.select("#img-8311e7c5-758")
   .data("mouseover_color", "#CD5C5C")
;
fig.select("#img-8311e7c5-758")
   .data("mouseout_color", "#6A6A6A")
;
fig.select("#img-8311e7c5-758")
   .click(Gadfly.zoomslider_zoomout_click)
.mouseenter(Gadfly.zoomslider_button_mouseover)
.mouseleave(Gadfly.zoomslider_button_mouseout)
;
    });
]]> </script>
</svg>

</div>


