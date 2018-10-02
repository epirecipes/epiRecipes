// URL: https://beta.observablehq.com/@epichef/forced-sir
// Title: Forced SIR
// Author: epichef (@epichef)
// Version: 132
// Runtime version: 1

const m0 = {
  id: "baf02793a26a64c3@132",
  variables: [
    {
      inputs: ["md"],
      value: (function(md){return(
md`# Forced SIR`
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`Author: Cameron Lack @c-lack

Date: 2018-10-02`
)})
    },
    {
      inputs: ["sir_sol","width","DOM","Plotly"],
      value: (function(sir_sol,width,DOM,Plotly)
{
  let t = sir_sol.t
  let I = sir_sol.y.map((x) => x[1])
  console.log(I)
  var itrace = {
    type: "scatter",
    mode: "step",
    name: 'I',
    x: t,
    y: I,
    line: {color: 'red'}
  }
  
  var data = [itrace];
  
  var layout = {
    width: width
  };
    
  const div = DOM.element('div');
  Plotly.newPlot(div, data, layout);
  return div;
}
)
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Parameters`
)})
    },
    {
      name: "viewof le",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 70,
  min: 50,
  max: 100,
  precision: 0,
  description: "Life expectancy"
})
)})
    },
    {
      name: "le",
      inputs: ["Generators","viewof le"],
      value: (G, _) => G.input(_)
    },
    {
      name: "mu",
      inputs: ["le"],
      value: (function(le){return(
1/(le*365)
)})
    },
    {
      name: "viewof beta0",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 1.4,
  min: 0,
  max: 2,
  precision: 1,
  description: "Average infectivity"
})
)})
    },
    {
      name: "beta0",
      inputs: ["Generators","viewof beta0"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof beta1",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 0.05,
  min: 0.01,
  max: 0.1,
  precision: 2,
  description: "Seasonal variation"
})
)})
    },
    {
      name: "beta1",
      inputs: ["Generators","viewof beta1"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof f",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 365,
  min: 1,
  max: 730,
  precision: 0,
  description: "Length of season (days)"
})
)})
    },
    {
      name: "f",
      inputs: ["Generators","viewof f"],
      value: (G, _) => G.input(_)
    },
    {
      name: "omega",
      inputs: ["f"],
      value: (function(f){return(
2*Math.PI/f
)})
    },
    {
      name: "viewof gamma",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 0.14,
  min: 0,
  max: 1,
  precision: 2,
  description: "Recovery rate"
})
)})
    },
    {
      name: "gamma",
      inputs: ["Generators","viewof gamma"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof tmax_years",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 2,
  min: 0,
  max: 5,
  precision: 1,
  description: "Length of simulation (years)"
})
)})
    },
    {
      name: "tmax_years",
      inputs: ["Generators","viewof tmax_years"],
      value: (G, _) => G.input(_)
    },
    {
      name: "tmax",
      inputs: ["tmax_years"],
      value: (function(tmax_years){return(
tmax_years*365
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Model code`
)})
    },
    {
      name: "sir_forced",
      inputs: ["beta0","beta1","omega","mu","gamma"],
      value: (function(beta0,beta1,omega,mu,gamma){return(
(dydt, y, t) => {
  let N = y[0] + y[1] + y[2]
  let force_inf = beta0*(1+beta1*Math.sin(omega*t))*y[0]*y[1]/N
  dydt[0] = mu*N - force_inf - mu*y[0]
  dydt[1] = force_inf - (gamma + mu)*y[1]
  dydt[2] = gamma*y[1] - mu*y[2]
}
)})
    },
    {
      name: "sir_sol",
      inputs: ["simulate","sir_forced","tmax"],
      value: (function(simulate,sir_forced,tmax){return(
simulate(sir_forced,0,[99999, 1, 0],tmax)
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Helper code`
)})
    },
    {
      name: "copy",
      value: (function(){return(
(x) => {
  return Object.assign({},x)
}
)})
    },
    {
      name: "simulate",
      inputs: ["rk4","copy"],
      value: (function(rk4,copy){return(
(f,t0,y0,tmax) => {
  let step = 1
  let integrator = rk4(y0, f, t0, step)
  let t = t0
  let y = copy(y0)
  let ts = []
  let ys = []
  // Simulate to equilibrium
  while (t < 1000*365) {
    t = t+step
    integrator = integrator.step()
  }
  t = 0
  ts.push(t)
  ys.push(copy(integrator.y))
  while (t < tmax) {
    t = t+step
    integrator = integrator.step()
    ys.push(copy(integrator.y))
    ts.push(t)
  }
  return {t:ts, y:ys}
}
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Imports`
)})
    },
    {
      name: "rk4",
      inputs: ["require"],
      value: (function(require){return(
require('https://bundle.run/ode-rk4')
)})
    },
    {
      name: "Plotly",
      inputs: ["require"],
      value: (function(require){return(
require("https://cdn.plot.ly/plotly-latest.min.js")
)})
    },
    {
      from: "@jashkenas/inputs",
      name: "slider",
      remote: "slider"
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Test`
)})
    },
    {
      value: (function()
{
  return 'Default parameters: le=70, beta0 = 1.4, beta0 = 0.05, f = 365, gamma = 0.14, tmax_years = 2. I_final = 14.646455760129378'
}
)
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Benchmark`
)})
    },
    {
      inputs: ["simulate","sir_forced","tmax"],
      value: (function(simulate,sir_forced,tmax)
{
  let start = new Date()
  let sir_sol = simulate(sir_forced,0,[99999, 1, 0],tmax)
  let time = new Date() - start
  return `Runs in ${time} milliseconds`
}
)
    }
  ]
};

const m1 = {
  id: "@jashkenas/inputs",
  variables: [
    {
      name: "slider",
      inputs: ["input"],
      value: (function(input){return(
function slider(config = {}) {
  let {value, min = 0, max = 1, step = "any", precision = 2, title, description, format, submit} = config;
  if (typeof config == "number") value = config;
  if (value == null) value = (max + min) / 2;
  precision = Math.pow(10, precision);
  return input({
    type: "range", title, description, submit, format,
    attributes: {min, max, step, value},
    getValue: input => Math.round(input.valueAsNumber * precision) / precision
  });
}
)})
    },
    {
      name: "input",
      inputs: ["html","d3format"],
      value: (function(html,d3format){return(
function input(config) {
  let {form, type = "text", attributes = {}, action, getValue, title, description, format, submit, options} = config;
  if (!form) form = html`<form>
	<input name=input type=${type} />
  </form>`;
  const input = form.input;
  Object.keys(attributes).forEach(key => {
    const val = attributes[key];
    if (val != null) input.setAttribute(key, val);
  });
  if (submit) form.append(html`<input name=submit type=submit style="margin: 0 0.75em" value="${typeof submit == 'string' ? submit : 'Submit'}" />`);
  form.append(html`<output name=output style="font: 14px Menlo, Consolas, monospace; margin-left: 0.5em;"></output>`);
  if (title) form.prepend(html`<div style="font: 700 0.9rem sans-serif;">${title}</div>`);
  if (description) form.append(html`<div style="font-size: 0.85rem; font-style: italic;">${description}</div>`);
  if (format) format = d3format.format(format);
  if (action) {
    action(form);
  } else {
    const verb = submit ? "onsubmit" : type == "button" ? "onclick" : type == "checkbox" || type == "radio" ? "onchange" : "oninput";
    form[verb] = (e) => {
      e && e.preventDefault();
      const value = getValue ? getValue(input) : input.value;
      if (form.output) form.output.value = format ? format(value) : value;
      form.value = value;
      if (verb !== "oninput") form.dispatchEvent(new CustomEvent("input"));
    };
    if (verb !== "oninput") input.oninput = e => e && e.stopPropagation() && e.preventDefault();
    if (verb !== "onsubmit") form.onsubmit = (e) => e && e.preventDefault();
    form[verb]();
  }
  return form;
}
)})
    },
    {
      name: "d3format",
      inputs: ["require"],
      value: (function(require){return(
require("d3-format")
)})
    }
  ]
};

const notebook = {
  id: "baf02793a26a64c3@132",
  modules: [m0,m1]
};

export default notebook;
