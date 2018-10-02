// URL: https://beta.observablehq.com/@epichef/sis-model
// Title: SIS Model
// Author: epichef (@epichef)
// Version: 483
// Runtime version: 1

const m0 = {
  id: "b86dec6a4c1ccf7d@483",
  variables: [
    {
      inputs: ["md"],
      value: (function(md){return(
md`# SIS Model`
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`### Plot`
)})
    },
    {
      inputs: ["sis_sol","width","beta","mu","DOM","Plotly"],
      value: (function(sis_sol,width,beta,mu,DOM,Plotly)
{
  var t = sis_sol.time
  var S = sis_sol.SI.map((x)=>{return x[0]})
  var I = sis_sol.SI.map((x)=>{return x[1]})
  var strace = {
    type: "scatter",
    mode: "step",
    name: 'S',
    x: t,
    y: S,
    line: {color: 'blue'}
  }
  var itrace = {
    type: "scatter",
    mode: "step",
    name: 'I',
    x: t,
    y: I,
    line: {color: 'red'}
  }
  var data = [strace,itrace];
  
  var layout = {
    width: width
  };
  var r = Math.round((beta/mu)*100)/100;
  var s_star = Math.round(S.slice(-1)[0]*10)/10;
  var i_star = Math.round(I.slice(-1)[0]*10)/10;
  // var brn = tex.block`R_0`
  var layout = {
    title: `R0 = ` + r + `, Endemic state: (S*, I*) = ` +`(` + s_star + ` , ` + i_star + `)`,
    xaxis: {
    title: 'Time',
    titlefont: {
      size: 18,
    }
    },
    yaxis: {
    title: 'Number',
    titlefont: {
      size: 18,
    }
  }};
  const div = DOM.element('div');
  Plotly.newPlot(div, data, layout);
  return div;

}

)
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`### Parameters`
)})
    },
    {
      name: "viewof mu",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 0.1,
  min: 0,
  max: 1,
  precision: 2,
  description: "Recovery Rate"
})
)})
    },
    {
      name: "mu",
      inputs: ["Generators","viewof mu"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof beta",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 0.4,
  min: 0,
  max: 5,
  precision: 2,
  description: "Infectivity parameter"
})
)})
    },
    {
      name: "beta",
      inputs: ["Generators","viewof beta"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof N",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 1000,
  min: 1,
  max: 5000,
  step: 1,
  precision: 2,
  description: "Population Size"
})
)})
    },
    {
      name: "N",
      inputs: ["Generators","viewof N"],
      value: (G, _) => G.input(_)
    },
    {
      name: "viewof I0",
      inputs: ["slider"],
      value: (function(slider){return(
slider({
  value : 1,
  min: 0,
  max: 1000,
  step: 1,
  precision: 2,
  description: "Initial Number of Infecteds"
})
)})
    },
    {
      name: "I0",
      inputs: ["Generators","viewof I0"],
      value: (G, _) => G.input(_)
    },
    {
      name: "tmax",
      value: (function(){return(
100
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`### Model Code`
)})
    },
    {
      name: "SIS",
      inputs: ["beta","N","mu"],
      value: (function(beta,N,mu){return(
function SIS(dydt, y, t) {
  dydt[0] = -beta*y[0]*y[1]/N + mu*y[1];
  dydt[1] = beta*y[0]*y[1]/N - mu*y[1];
}
)})
    },
    {
      name: "sis_sol",
      inputs: ["simulate","SIS","N","I0"],
      value: (function(simulate,SIS,N,I0){return(
simulate(SIS, 0.01, [N-I0,I0], 0, 100)
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`### Utility Functions`
)})
    },
    {
      name: "simulate",
      inputs: ["rk4","copy"],
      value: (function(rk4,copy){return(
function simulate(func, step, initial, t0, maxtime){
  var integrator = rk4( initial, func, t0, step )
  var time = []
  var SI = []
  time.push(t0)
  SI.push(initial)
  var t = t0
  while(t<maxtime){
    t += step
    integrator = integrator.step()
    SI.push(copy(integrator.y))
    time.push(t)}

  return {time, SI};
}
)})
    },
    {
      name: "copy",
      value: (function(){return(
function copy(x) {
  return Object.assign({},x)
}
)})
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`### Dependencies`
)})
    },
    {
      name: "rk4",
      inputs: ["require"],
      value: (function(require){return(
require(`https://bundle.run/ode-rk4`)
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
md`### Test + Benchmark`
)})
    },
    {
      inputs: ["simulate","SIS","N","I0"],
      value: (function(simulate,SIS,N,I0)
{let times = 0;
 let i;
 for (i = 0; i < 1; i++) {
   const start = new Date();
   let sis_sol = simulate(SIS, 0.01, [N-I0,I0], 0, 100);
   const elapsed = new Date() - start;
   times += elapsed;
 }
 return 'Average time elapsed over 1 replicates: ' + times/1 + ' ms';
}
)
    },
    {
      inputs: ["md"],
      value: (function(md){return(
md`Fix parameters to beta = 1, mu = 0.1, I0 = 1, N = 1000, tmax = 100`
)})
    },
    {
      inputs: ["simulate","SIS","N","I0","tmax"],
      value: (function(simulate,SIS,N,I0,tmax)
{
 const beta = 1;
 const mu = 0.01
 const test_sol = simulate(SIS, 0.01,[N-I0, I0], 0, tmax)
 return `Final number of infected individuals: ` + Math.round(test_sol.SI.map((x)=>{return x[1]}).slice(-1))}
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
  id: "b86dec6a4c1ccf7d@483",
  modules: [m0,m1]
};

export default notebook;
