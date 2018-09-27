---
interact_link: notebooks/sir/cpp.ipynb
title: 'C++'
permalink: 'chapters/sir/cpp'
previouschapter:
  url: chapters/sir/scilab
  title: 'Scilab'
nextchapter:
  url: chapters/sir/xpp
  title: 'xppaut'
redirect_from:
  - 'chapters/sir/cpp'
---

## SIR model in C++ using Boost odeint


{:.input_area}
```python
%%writefile sir.cpp
#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

const double b = 0.1;
const double g = 0.05;

typedef boost::array< double , 3 > state_type;

void sir( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = -b * x[0] * x[1];
    dxdt[1] = b * x[0] * x[1] - g * x[1];
    dxdt[2] = g * x[1];
}

void write_sir(const state_type &x , const double t )
{
    cout << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << endl;
}

int main(int argc, char **argv)
{
    state_type x = { 0.99 , 0.01 , 0.0 }; // initial conditions
    integrate( sir , x , 0.0 , 200.0 , 0.1 , write_sir );
}
```

{:.output_stream}
```
Writing sir.cpp

```

We compile the code first.


{:.input_area}
```python
!g++ -O3 sir.cpp -o sir
```

We run the file and redirect the output to a file


{:.input_area}
```python
!./sir > sir_cpp.out
```

### Visualisation

We can graph the results in the output file e.g. using Python and Matplotlib.


{:.input_area}
```python
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("ggplot")
```


{:.input_area}
```python
sir_out = pd.read_csv("sir_cpp.out",sep=" ",header=None,names=["t","S","I","R"],index_col=False)
```


{:.input_area}
```python
sline = plt.plot("t","S","",data=sir_out,color="red",linewidth=2)
iline = plt.plot("t","I","",data=sir_out,color="green",linewidth=2)
rline = plt.plot("t","R","",data=sir_out,color="blue",linewidth=2)
plt.xlabel("Time",fontweight="bold")
plt.ylabel("Number",fontweight="bold")
legend = plt.legend(title="Population",loc=5,bbox_to_anchor=(1.25,0.5))
frame = legend.get_frame()
frame.set_facecolor("white")
frame.set_linewidth(0)
```


![png](../../images/chapters/sir/cpp_10_0.png)

