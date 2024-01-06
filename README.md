## Project Title 
Atomistic Monte Carlo (MC) Simulation & Structual Analysis via Radial Distribution Function (RDF)

## Motivation
Course project in ChE696: Simulation of Condensed Matter System

## Content
Monte_Carlo-Inline.cpp simulates particles in an orthogonal box using Markov Chain Monte Carlo (MCMC) method and widom insertion to calculating system thermodynamic properties. 
RDF_Obj_Oriented.cpp calculates the RDF of a system once given the particle trajectory.


## Build & Execute
To build the MC code: g++ -O3 -o Monte_Carlo Monte_Carlo-Inline.cpp
To execute the MC code: ./Monte_Carlo <random_seed> <reduced_density> | tee Monte_Carlo.log

To build the RDF code: g++ -O3 -o RDF_Obj_Oriented RDF_Obj_Oriented.cpp
To execute the RDF code: ./RDF_Obj_Oriented <trajectory_file>

## Credit
The skeleton of the code was provided by the instructor Dr. Rebecca Lindsey.
