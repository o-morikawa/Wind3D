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

## Sample code
```julia
using Random
using Gaugefields
using LinearAlgebra
using Wilsonloop

using Plots

function UN_test_3D(NX,NY,NT,NC,η;λ=1,rε=0.1)

    Dim = 3
    L = NX

    n = 4
    println("Run Configurations: n=$n")

    eps = 0.01
    flow_number = 6000
    step = 10
    println("L=$L, eps=$eps,flow=$(eps*flow_number)")

    eta = η
    randscale = λ
    rand_eps = rε
    println("eta: $eta,  randscale: $randscale,  rand_eps: $(rand_eps)")

    for i = 1:n

        Random.seed!(123)

        if i == 1
            println("Configuration 1: Cold start")
            U = Initialize_3D_UN_Gaugefields(
                NC,NX,NY,NT,
                condition = "cold",
                randomnumber="Random"
            )
        elseif i == 2
            println("Configuration 2: Hot start")
            U = Initialize_3D_UN_Gaugefields(
                NC,NX,NY,NT,
                condition = "hot",
                randomnumber="Random",
                randscale=randscale,
            )
        elseif i == 3
            println("Configuration 3: Test mapping as T^3->SU(2)")
            U = Initialize_3D_UN_Gaugefields(
                NC,NX,NY,NT,
                condition = "test_map",
                m = -1, # -1, 1, 3, 5, 7
            )
        elseif i == 4
            println("Configuration 4: Test mapping as T^3->SU(2) with Random noise")
            U = Initialize_3D_UN_Gaugefields(
                NC,NX,NY,NT,
                condition = "test_map_rand",
                m = get_mass(i),
                randomnumber="Random",
                reps = rand_eps,
            )
        end
        println(typeof(U))

        temps = Temporalfields(U, num=9)
        println(typeof(temps))

        println(winding_UN_3D(U,temps))

        if i==1 || i==2
            g = Gradientflow_eta_3D(U, eta, eps=eps)
        else
            g = Gradientflow_TA_eta_3D(U, eta, eps=eps)
        end
        flownumber = flow_number

        W = winding_UN_3D(U,temps)
        Wh =winding_UN_3D(U,1,temps)
        S = calc_gdgaction_3D(U,eta,temps)
        D = det_unitary(U)

        println("flow time: 0")
        println("W=$(W), Wh=$(Wh), S=$(S), Det=$(D)")

        for iflow = 1:flownumber
            flow!(U, g)
            if iflow%step==0
                W = winding_UN_3D(U,temps)
                Wh =winding_UN_3D(U,1,temps)
                S = calc_gdgaction_3D(U,eta,temps)
                D = det_unitary(U)
                println("flow time: $(iflow*eps)")
                println("W=$(W), Wh=$(Wh), S=$(S), Det=$(D)")
            end
        end
    end

end


```
