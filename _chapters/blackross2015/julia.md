---
interact_link: notebooks/blackross2015/julia.ipynb
title: 'Julia'
permalink: 'chapters/blackross2015/julia'
previouschapter:
  url: chapters/blackross2015/matlab
  title: 'Original Matlab code'
nextchapter:
  url: chapters/blackross2015/octave
  title: 'Octave'
redirect_from:
  - 'chapters/blackross2015/julia'
---

## Final size for SIR model


{:.input_area}
```julia
using Plots
```


{:.input_area}
```julia
@fastmath function SIR_fs(N,bet,gamm)
    final_size = zeros(N+1,1);
    final_size[2] = 1;
    for Z2 in 0:N
        @inbounds for Z1 in Z2+1:N-1
            p1 = 1 / ( 1 + gamm/(bet*(N-Z1)));
            final_size[Z1+2] = final_size[Z1+2] + final_size[Z1+1]*p1; 
            final_size[Z1+1] = final_size[Z1+1]*(1-p1);
        end
    end
    return final_size;
end
```




{:.output_data_text}
```
SIR_fs (generic function with 1 method)
```




{:.input_area}
```julia
N = 20;                       
bet = 2/(N-1);
gamm = 1.0;
```


{:.input_area}
```julia
final_size = SIR_fs(N,bet,gamm);
```


{:.input_area}
```julia
@time final_size = SIR_fs(N,bet,gamm);
```

{:.output_stream}
```
  0.000007 seconds (84 allocations: 6.342 KiB)

```


{:.input_area}
```julia
bar(0:N,final_size)
```




![svg](../../images/chapters/blackross2015/julia_6_0.svg)



### Final size for SI(4)R model


{:.input_area}
```julia
@fastmath function SI4R_fs(N::Int64,bet::Float64,gamm::Float64)
    Psi = (N+1)*(N+2)*(N+3)*(N+4)/24;
    final_size = zeros(N+1,1);
    p_vec = zeros(Psi,1);
    p_vec[2]=1;
    @inbounds for Z5 in 0:N
        w::Int64 = Psi - (N-Z5+1)*(N-Z5+2)*(N-Z5+3)*(N-Z5+4)/24 + 1; 
        @inbounds for Z4 in Z5:N
            a5 = 4*gamm*(Z4-Z5);
            @inbounds for Z3 in Z4:N
                a4 = 4*gamm*(Z3-Z4);
                @inbounds for Z2 in Z3:N;
                    a3 = 4*gamm*(Z2-Z3);            
                    @inbounds for Z1 in Z2:N
                    a1 = bet*(N-Z1)*(Z1-Z5);   
                    a2 = 4*gamm*(Z1-Z2);
                    tot = a1+a2+a3+a4+a5;
                    if Z1-Z5 == 0
                        final_size[Z5+1] = p_vec[w];
                    end
                    if a1 > 0
                       p_vec[w+1] = p_vec[w+1]+ p_vec[w]*a1/tot;
                    end
                    if a2 > 0
                        p_vec[w+N-Z2] = p_vec[w+N-Z2]+ p_vec[w]*a2/tot;
                    end
                    if a3 > 0
                        place3::Int64 = (N-Z3)*(N-Z3+1)/2;
                        p_vec[w+place3] = p_vec[w+place3]+ p_vec[w]*a3/tot;
                    end
                    if a4 > 0
                        place4::Int64 = (N-Z4)*(N-Z4+1)*(N-Z4+2)/6; 
                        p_vec[w+place4] = p_vec[w+place4] + p_vec[w]*a4/tot;
                    end
                    if a5 > 0
                        p_vec[w] = p_vec[w]*a5/tot;
                    end
                    w = w + 1;
                end
            end
        end
    end
end
    return final_size
end
```




{:.output_data_text}
```
SI4R_fs (generic function with 1 method)
```




{:.input_area}
```julia
final_size = SI4R_fs(N,bet,gamm);
```


{:.input_area}
```julia
@time final_size = SI4R_fs(N,bet,gamm);
```

{:.output_stream}
```
  0.000751 seconds (7 allocations: 83.547 KiB)

```


{:.input_area}
```julia
bar(0:N,final_size)
```




![svg](../../images/chapters/blackross2015/julia_11_0.svg)


