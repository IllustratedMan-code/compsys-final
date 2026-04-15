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

original_sol = DE.solve(problem)

sol = original_sol(original_sol.t[original_sol.t .>= 500])


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

# plot(solution["time"], solution["WC-1_tot"], label=L"WC-1_{tot}", xlims=(0, 48))
# plot!(solution["time"], solution["FRQ_tot"], label=L"[FRQ_{tot}]", xlims=(0, 48))
plot(solution["time"], solution["WC-1_tot"], label=L"WC-1_{tot}")
plot!(solution["time"], solution["FRQ_tot"], label=L"[FRQ_{tot}]")
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

# For Figures 3A-C and 4A-C we're using Bifurcationkit
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


u_end = sol.u[end-1] # this doesn't give hopf points in figure4A for some reason
u_end = fill(1.0, 9) # have to just set everything to one instead



function calculate_period(optic, params=params)
    # define the bifurcation based on our parameter (i.e. k01)
    prob = BK.BifurcationProblem(bifurcation_model, u_end, params, optic)

    # Establish the sweep for the k0 param
    k0_range = BK.ContinuationPar(
        p_min = 0.0, 
        p_max = 2.0,
        ds = 0.01,
        dsmax = 0.01 # ensures that we get lots of points for plotting
    )

    # calculate hopf points
    br = BK.continuation(prob, BK.PALC(), k0_range)
    plot(br)
    savefig("br.svg")
    print(br)
    if length(br.specialpoint) < 2
        # ideally we only get two special points, the end point and the hopf point
        # if not, then return this object for debugging
        print("There are no hopf points")
        return br
    end

    # trapezoidal estimation of period orbits
    # this isn't exactly equivalent to xppaut, but close enough
    # and it's much faster to compute
    period_prob = BK.PeriodicOrbitTrapProblem(M=200)
    
    # compute the period vs param based on the hopf point
    period_branch = BK.continuation(
        br, 1, # Ensure '1' is the index of the hopf point
        k0_range, 
        period_prob;
    )
    
    return period_branch

end


period_branch = calculate_period((@optic _.k0.one))
plot(period_branch, vars = (:param, :period),
     xlabel = L"rate of $frq$ overexpression $k_{01}$",
     ylabel = "period, h",
     xlims = (0, 0.21),
     ylims = (16, 28)
     )

savefig("figure3A.svg")

#https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/tutorials/ode/tutorialsCodim2PO/#Periodic-predator-prey-model
#There is no hopf point for this one, need to use the above to take our established period
# from the ode solution

one_period_problem = DE.ODEProblem(model, original_sol.u[end], (0, 25), params)

one_period_sol = DE.solve(one_period_problem)

period_branch = calculate_period((@optic _.k0.three))
plot(period_branch)
savefig("figure2b_no_hopf_points.svg")

prob_bif = BK.ODEBifProblem(bifurcation_model, u_end, params, (@optic _.k0.three);)

k03_range = BK.ContinuationPar(
    p_min = 0.0, 
    p_max = 1.0,
    ds = 0.01,
    dsmax = 0.01
)

period_prob = BK.PeriodicOrbitTrapProblem(M=200)
probtrap, ci = BK.generate_ci_problem(period_prob,
	                              prob_bif, original_sol, (0, 22))

opts_po_cont = BK.ContinuationPar(k03_range, max_steps = 1800, tol_stability = 1e-5)
brpo_fold = BK.continuation(probtrap, ci, BK.PALC(), opts_po_cont;
	verbosity = 3, plot = true, argspo...)
plot(brpo_fold, vars = (:param, :period),
     xlabel = L"rate of $wc-1$ overexpression $k_{03}$",
     ylabel = "period, h",
     xlims = (0, 1),
     ylims = (19, 24)
     )

# savefig("figure3B.svg")


# need to set csp-1 expression constant k16=0
params = (k=k, K=K, k0=k0)
period_branch = calculate_period((@optic _.k0.four))
