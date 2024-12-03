#using LatticeQCD

using Gaugefields
using LinearAlgebra
using Wilsonloop

using Plots

function UN_msearch_3D(NX,NY,NT,NC)

    Dim = 3
    L = NX

    println("Test mapping configuration (improved action)")
    n = 1

    println("L=$L")

    m = -1.05:0.1:7.05
    ml = length(m)
    η = [1, 0, -10, -20]
    s = zeros(Float64, ml, 4)

    for i = 1:ml
        println("Run with $(m[i])...")

        U = Initialize_3D_UN_Gaugefields(
            NC,NX,NY,NT,
            condition = "test_map",
            m = m[i],
        )

        temps = Temporalfields(U, num=9)

        for j = 1:4
            s[i,j] = calc_gdgaction_3D(U,η[j],temps)
        end
    end

    println("Output...")
    open("./msearch_L$(L).csv", "w") do f
        write(f, "m, seta1, seta2, seta3, seta4\n")
        for i = 1:ml
            write(f, "$(m[i]), $(s[i,1]), $(s[i,2]), $(s[i,3]), $(s[i,4])\n")
        end
    end

    plt = plot(m, s[:,1], label="η=-1")
    plot!(plt, m, s[:,2], label="η=0")
    plot!(plt, m, s[:,3], label="η=-10")
    plot!(plt, m, s[:,4], label="η=-20")
    savefig("./msearch_L$(L).png")

end


function main()
    L = 30
    
    NX = L
    NY = L
    NT = L
    NC = 2
    
    @time UN_msearch_3D(NX,NY,NT,NC)
end
main()



