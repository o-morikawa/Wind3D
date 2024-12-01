#using LatticeQCD

using Random
using Dates
using Gaugefields
using LinearAlgebra
using Wilsonloop

using Plots

function UN_adm_3D(NX,NY,NT,NC,η,m)

    Dim = 3
    L = NX

    println("Test mapping configuration (improved action)")
    n = 1

    eps = 0.01
    if L==10
        flow_number = 4000
    elseif L==20
        flow_number = 8000
    elseif L==30
        flow_number = 16000
    elseif L==40
        flow_number = 27000
    elseif L==50
        flow_number = 42000
    end
    step = 100

    println("L=$L, eps=$eps,flow=$(eps*flow_number)")

    eta = η
    println("eta: $eta")

    s = zeros(Float64, NX, NY, Int(flow_number/step)+1)
    tslice = div(NT, 2)

    for i = 1:n

        #Random.seed!(123)
        t0 = Dates.DateTime(2024,1,1,16,10,7)
        t  = Dates.now()
        Random.seed!(Dates.value(t-t0))

        U = Initialize_3D_UN_Gaugefields(
            NC,NX,NY,NT,
            condition = "test_map",
            m = m,
        )
        println(typeof(U))

        temps = Temporalfields(U, num=9)
        println(typeof(temps))

        println(winding_UN_3D(U,temps))

        g = Gradientflow_TA_eta_3D(U, eta, eps=eps)
        flownumber = flow_number

        j = 1
        for iy = 1:NY
            for ix = 1:NX
                s[ix,iy,j] = calc_gdgactiondensity_3D(U,ix,iy,tslice,eta,temps)
            end
        end

        for iflow = 1:flownumber
            flow!(U, g)
            if iflow%step==0
                j += 1
                for iy = 1:NY
                    for ix = 1:NX
                        s[ix,iy,j] = calc_gdgactiondensity_3D(U,ix,iy,tslice,eta,temps)
                    end
                end
            end
        end
    end

    #println(w)

    flow = 0:(eps*step):(eps*flow_number)
    lt = length(flow)

    for it=1:lt
        open("./adm_map_L$(L)_eta$(eta)_flow$(flow[it]).csv", "w") do f
            write(f, "x, y, Adm\n")
            for iy = 1:NY
                for ix = 1:NX
                    write(f, "$ix, $iy, $(s[ix,iy,it])\n")
                end
                write(f, "\n")
            end
        end
    end

end


function main()
    L = 10
    
    NX = L
    NY = L
    NT = L
    NC = 2

    η = -20
    m = 1
    
    @time UN_adm_3D(NX,NY,NT,NC,η,m)
end
main()



