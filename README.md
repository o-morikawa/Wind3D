# Wind3D

We consider the winding number on numerical simulations.
We propose an effective approach to give an approximate integer for
```math
\int \tr(g^{-1}dg)^{3},
```
where g is a smooth map: X -> U(N).
To this end, we utilize a Julia repository, o-morikawa/Gaugefields.jl,
and formulate a gradient-flow method even on a course lattice.

- src: main function on Julia (o-morikawa/Gaugefields.jl)
- output: data and simple figures
- m_nb: Mathematica notebooks for small lattice size
