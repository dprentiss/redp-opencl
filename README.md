# Regional Evacuation Design Problem
This repository hosts the code associated with "Evacuation Route Modeling and Planning with General Purpose GPU Computing". [(here)](https://drum.lib.umd.edu/handle/1903/15954)

## Abstract
> This work introduces a bilevel, stochastic optimization problem aimed at robust, regional evacuation network design and shelter location under uncertain hazards. A regional planner, acting as a Stackelberg leader, chooses among evacuation-route contraflow operation and shelter location to minimize the expected risk exposure to evacuees. Evacuees then seek an equilibrium with respect to risk exposure in the lower level. An example network is solved exactly with a strategy that takes advantage of a fast, low-memory, equilibrium algorithm and general purpose computing on graphical processing units.

## Code
* /algorithmB

  An OpenCL implementation of Robert Dial's 'Algorithm B'

* /python

  A python script to generate input text files for the REDP solvers

* /redp

  A single cpu version of the REDP solver in C++

* /redpcl

  A GPU version of the REDP solver in OpenCL/C++
