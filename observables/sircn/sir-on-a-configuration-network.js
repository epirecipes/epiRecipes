// URL: https://beta.observablehq.com/@epichef/sir-on-a-configuration-network
// Title: SIR on a configuration network
// Author: epichef (@epichef)
// Version: 94
// Runtime version: 1

const m0 = {
  id: "f7f315e722948697@94",
  variables: [
    {
      inputs: ["md"],
      value: (function(md){return(
md`# SIR on a configuration network`
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
      inputs: ["sir_sol","psi","k","width","DOM","Plotly"],
      value: (function(sir_sol,psi,k,width,DOM,Plotly)
{
  let t = sir_sol.t
  let theta = sir_sol.y.map((x) => x[0])
  let S = sir_sol.y.map((x) => psi(x[0],k))
  let I = sir_sol.y.map((x) => 1 - psi(x[0],k) - x[1])
  let R = sir_sol.y.map((x) => x[1])
  var strace = {
    type: "scatter",
    mode: "step",
    name: 'S',
    x: t,
    y: S,
    line: {color: 'green'}
  }
  var itrace = {
    type: "scatter",
    mode: "step",
    name: 'I',
    x: t,
    y: I,
    line: {color: 'red'}
  }
  var thetatrace = {
    type: "scatter",
    mode: "step",
    name: 'theta',
    x: t,
    y: theta,
    line: {color: 'yellow'}
  }
  var rtrace = {
    type: "scatter",
    mode: "step",
    name: 'R',
    x: t,
    y: R,
    line: {color: 'blue'}
  }
  
  var data = [strace,itrace,rtrace,thetatrace];
  
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
      name: "viewof beta",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 0.1,
  min: 0,
  max: 1,
  precision: 2,
  description: "Infectivity"
})
)})
    },
    {
      name: "beta",
      inputs: ["Generators","viewof beta"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof gamma",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 0.05,
  min: 0,
  max: 1,
  precision: 2,
  description: "Recovery"
})
)})
    },
    {
      name: "gamma",
      inputs: ["Generators","viewof gamma"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof k",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 5,
  min: 1,
  max: 10,
  precision: 0,
  description: "Degree"
})
)})
    },
    {
      name: "k",
      inputs: ["Generators","viewof k"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof tmax",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 125,
  min: 0,
  max: 150,
  precision: 0,
  description: "Simulationm time"
})
)})
    },
    {
      name: "tmax",
      inputs: ["Generators","viewof tmax"],
      value: (G, _) => G.input(_)
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Model code`
)})
    },
    {
      name: "sir_config_net",
      inputs: ["beta","dpsi","k","gamma","psi"],
      value: (function(beta,dpsi,k,gamma,psi){return(
(dydt, y, t) => {
  dydt[0] = -beta*y[0] + beta*dpsi(y[0],k)/dpsi(1,k) + gamma*(1-y[0])
  let S = psi(y[0],k)
  let I = 1-S-y[1]
  dydt[1] = gamma*I
}
)})
    },
    {
      name: "psi",
      value: (function(){return(
(theta,k) => {
  return Math.pow(theta,k)
}
)})
    },
    {
      name: "dpsi",
      value: (function(){return(
(theta,k) => {
  return k*Math.pow(theta,k-1)
}
)})
    },
    {
      name: "sir_sol",
      inputs: ["simulate","sir_config_net","tmax"],
      value: (function(simulate,sir_config_net,tmax){return(
simulate(sir_config_net,0,[0.999,0],tmax)
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Helper functions`
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
  let step = tmax/200
  let integrator = rk4(y0, f, t0, step)
  let ts = []
  let ys = []
  ts.push(t0)
  ys.push(copy(y0))
  let t = t0
  while (t < tmax) {
    integrator = integrator.step()
    t = t + step
    ts.push(t)
    ys.push(copy(integrator.y))
  }
  return {t:ts, y:ys}
}
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`## Dependencies`
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
  return 'Default parameters: beta = 0.1, gamma = 0.05, k = 5, tmax = 125. theta_final = 0.342509'
}
)
    },
    {
      inputs: ["simulate","sir_config_net","tmax"],
      value: (function(simulate,sir_config_net,tmax)
{
  let start = new Date()
  let sir_sol = simulate(sir_config_net,0,[0.999, 0],tmax)
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
  id: "f7f315e722948697@94",
  modules: [m0,m1]
};

export default notebook;
