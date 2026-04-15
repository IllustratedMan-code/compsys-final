import DifferentialEquations as DE
using Plots
using LaTeXStrings
import BifurcationKit as BK
import BifurcationKit: @optic


# rate constants
k = (1.8, 1.8, 0.05, 0.23, 0.27, 0.27, 0.5, 1, 40, 0.1, 0.05, 0.02, 50, 1, 5, 0.12, 1.4, 50, 1.4)

# these are 0 indexed because the authors hate me :),
# so add 1 to all indicies compared to the paper. e.g. K = K[1] K_1 = K[2] etc
K = (1.25, 3, 1, 10, 3)


k0 = (one = 0.0, three = 0.0, four = 0.0)

const DUI = (
    frq_mRNA = 1, FRQ_c = 2, FRQ_n = 3, 
    wc1_mRNA = 4, WC1_c = 5, WC1_n = 6, 
    FRQ_n_WC1_n = 7, csp1_mRNA = 8, CSP1 = 9
)
# function model(du, u, p, t)
#     k = p.k
#     k0 = p.k0
#     K = p.K
    
#     dui = Dict("frq mRNA" => 1,
#               "FRQ_c" => 2,
#               "FRQ_n" => 3,
#               "wc-1 mRNA" => 4,
#               "WC-1_c" => 5,
#               "WC-1_n" => 6,
#               "FRQ_n:WC-1_n" => 7,
#               "csp-1 mRNA" => 8,
#               "CSP-1" => 9)

#     u = Dict{String, Any}("frq mRNA" => u[1],
#               "FRQ_c" => u[2],
#               "FRQ_n" => u[3],
#               "wc-1 mRNA" => u[4],
#               "WC-1_c" => u[5],
#               "WC-1_n" => u[6],
#               "FRQ_n:WC-1_n" => u[7],
#               "csp-1 mRNA" => u[8],
#              "CSP-1" => u[9])

#     du[dui["frq mRNA"]] = k[1] * ((u["WC-1_n"]^6)/(K[1] + u["WC-1_n"]^6)) - k[4] * u["frq mRNA"] + k0[1]

#     du[dui["FRQ_c"]] = k[2] * u["frq mRNA"] - (k[3] + k[5]) * u["FRQ_c"]

#     du[dui["FRQ_n"]] = k[3] * u["FRQ_c"] + k[14] * u["FRQ_n:WC-1_n"] - u["FRQ_n"] * (k[6] + k[13] * u["WC-1_n"])

#     du[dui["wc-1 mRNA"]] = k[7] * (K[2]/(K[2] + u["CSP-1"])) - k[10] * u["wc-1 mRNA"] + k0[3]

#     du[dui["WC-1_c"]] = k[8] * ((u["FRQ_c"]^2)/(K[3] + u["FRQ_c"]^2)) * (u["wc-1 mRNA"]/(K[4]+ u["wc-1 mRNA"])) - (k[9] + k[11]) * u["WC-1_c"]

#     du[dui["WC-1_n"]] = k[9] * u["WC-1_c"] - u["WC-1_n"] * (k[12] + k[13] * u["FRQ_n"]) + k[14] * u["FRQ_n:WC-1_n"]

#     du[dui["FRQ_n:WC-1_n"]] = k[13] * u["FRQ_n"] * u["WC-1_n"] - (k[14] + k[15]) * u["FRQ_n:WC-1_n"]

#     du[dui["csp-1 mRNA"]] = k[16] * u["WC-1_n"] * (K[5]/(K[5] + u["CSP-1"])) - k[17] * u["csp-1 mRNA"] + k0[4]

#     du[dui["CSP-1"]] = k[18] * u["csp-1 mRNA"] - k[19] * u["CSP-1"]
    

# end

function model(du, u, p, t)
    k, k0, K = p.k, p.k0, p.K
    
    # Map u to local variables (zero allocation, AD-safe)
    f_mRNA = u[1]; FRQ_c = u[2]; FRQ_n = u[3]
    w_mRNA = u[4]; WC1_c = u[5]; WC1_n = u[6]
    Complex = u[7]; c_mRNA = u[8]; CSP1 = u[9]

    # Math logic (Types will promote naturally here)
    du[DUI.frq_mRNA] = k[1] * ((WC1_n^6)/(K[1] + WC1_n^6)) - k[4] * f_mRNA + k0.one

    du[DUI.FRQ_c] = k[2] * f_mRNA - (k[3] + k[5]) * FRQ_c

    du[DUI.FRQ_n] = k[3] * FRQ_c + k[14] * Complex - FRQ_n * (k[6] + k[13] * WC1_n)

    du[DUI.wc1_mRNA] = k[7] * (K[2]/(K[2] + CSP1)) - k[10] * w_mRNA + k0.three

    du[DUI.WC1_c] = k[8] * ((FRQ_c^2)/(K[3] + FRQ_c^2)) * (w_mRNA/(K[4] + w_mRNA)) - (k[9] + k[11]) * WC1_c

    du[DUI.WC1_n] = k[9] * WC1_c - WC1_n * (k[12] + k[13] * FRQ_n) + k[14] * Complex

    du[DUI.FRQ_n_WC1_n] = k[13] * FRQ_n * WC1_n - (k[14] + k[15]) * Complex

    du[DUI.csp1_mRNA] = k[16] * WC1_n * (K[5]/(K[5] + CSP1)) - k[17] * c_mRNA + k0.four

    du[DUI.CSP1] = k[18] * c_mRNA - k[19] * CSP1
end


u0 = fill(10.0, 9) # need to give the initial condition at least a 1 to occilate
tspan = (0.0, 3000)

params = (k=k, K=K, k0=k0)

problem = DE.ODEProblem(model, u0, tspan, params)

sol = DE.solve(problem)

sol = sol(sol.t[sol.t .>= 500])


solution = Dict("frq mRNA" => sol[1, :],
                "FRQ_c" => sol[2, :],
                "FRQ_n" => sol[3, :],
                "wc-1 mRNA" => sol[4, :],
                "WC-1_c" => sol[5, :],
                "WC-1_n" => sol[6, :],
                "FRQ_n:WC-1_n" => sol[7, :],
                "csp-1 mRNA" => sol[8, :],
                "CSP-1" => sol[9, :],
                "WC-1_tot" => sol[5, :] + sol[6, :] + sol[7, :],
                "FRQ_tot" => sol[2, :] + sol[3, :] + sol[7, :],
                "time" => sol.t .- 500
                )

plot(solution["time"], solution["WC-1_tot"], label=L"WC-1_{tot}", xlims=(0, 48))
plot!(solution["time"], solution["FRQ_tot"], label=L"[FRQ_{tot}]", xlims=(0, 48))
xlabel!("Time (h)")
ylabel!(L"$FRQ_{tot}$, $WC-1_{tot}$, a.u.")
savefig("figure2A.svg")

plot(solution["time"], solution["WC-1_n"], label=L"WC-1_n", xlims=(0, 48))
plot!(solution["time"], solution["FRQ_n"], label=L"[FRQ_n]", xlims=(0, 48))
xlabel!("Time (h)")
ylabel!(L"$WC-1_n$, [$FRQ_n$], a.u.")
savefig("figure2B.svg")

plot(solution["time"], solution["frq mRNA"], label=L"$frq$ mRNA", xlims=(0, 48))
plot!(solution["time"], solution["wc-1 mRNA"], label=L"$wc-1$ mRNA", xlims=(0, 48))
plot!(solution["time"], solution["csp-1 mRNA"], label=L"$csp-1$ mRNA", xlims=(0, 48))
xlabel!("Time (h)")
ylabel!(L"[$frq$ mRNA], [$FRQ_n$], a.u.")
savefig("figure2C.svg")

# For Figures 3A-C we're using Bifurcationkit
# see: https://docs.sciml.ai/BifurcationKit/stable/tutorials/tutorials3/#brusauto
# Use a NamedTuple for indices instead of a Dict for better performance/type safety




function bifurcation_model(u, p)
    param_type = typeof(first(values(p.k0)))
    
    # 2. Promote between the state 'u' and the parameter type
    T = promote_type(eltype(u), param_type)
    
    # 3. Create du with the promoted type
    du = similar(u, T)

    model(du, u, p, 0.0)
    return du
end


prob = BK.BifurcationProblem(bifurcation_model, u0, params, (@optic _.k0.one))

k0_range = BK.ContinuationPar(p_min=0.0, p_max=1.0)

br = BK.continuation(prob, BK.PALC(), k0_range)

opt_po = BK.NewtonPar(tol = 1e-5, verbose = true, max_iterations = 1000)
opts_po_cont = BK.ContinuationPar(
    dsmin = 0.001,
    dsmax = 0.04, ds = 0.001,
    p_max = 2.2,
    max_steps = 3000,
    newton_options = opt_po,
    plot_every_step = 1,
    nev = 11,
    tol_stability = 1e-6,
		)

#period_prob = BK.PeriodicOrbitTrapProblem(M=3000)
period_prob = BK.PeriodicOrbitOCollProblem(250, 4)


period = BK.continuation(br, 1, opts_po_cont, period_prob, plot=true, 
                         ampfactor = 0.001,
record_from_solution = (x, p; k...) -> begin
        # Extract period from the internal state x
        return (period = BK.getperiod(p.prob, x, p.p),
                max_x = maximum(x)) # or whatever variable you want
    end
                         )

plot(period)
