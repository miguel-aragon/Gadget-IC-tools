# Init

## Descrition

Codes to lower the resolution and fix periodic bondaries of gadget initial conditions files.

### lower_res

Adds adjacent particles in the initial conditions grid to lower the resolution by a factor of two in each dimension

## fix_periodic

Fixes particles outside the computational box of Gadget initial conditions. This can happen after applying Zel'dovich. In any case the initial conditions code should already account for this. 


