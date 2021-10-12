using DrWatson
@quickactivate "OptimalEstimationProject"
using MATLAB
using IndirectTrajOpt
using StaticArrays
using LinearAlgebra
using DataFrames

function main()

    # File to save data 
    fid = open(datadir() * "\\MFSolutionTrajTOFContinuation.txt","w")

    # Data info
    tf      = 33.1

    # ME Initial CoStates
    λ0 = [5.963713007087376, 12.99948134789672, 1.6584130970338573, -0.04062357611978672, 0.018558730103857064, -0.0004993971529679529, 0.1220053310292062] 

    try
        # Init parameters
        ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
        μ = ps.crp.μ
        LU = ps.crp.LU
        TU = ps.crp.TU

        # Trajectory time span
        tspan = (0.0, tf)

        # Initial conditions
        ics = [-0.019488511458668, -0.016033479812051, 0.0,
                8.918881923678198, -4.081793688818725, 0.0,
                1.0] 

        # Target final conditions
        fcs = [(1.0 - μ) + 6060.483/LU, 19452.284/LU, -34982.968/LU, 0.082677*(TU/LU), 0.006820*(TU/LU), -0.368434*(TU/LU), 0.0]
        xf = @SVector [fcs[1], fcs[2], fcs[3], fcs[4], fcs[5], fcs[6]]

        # Solve for increasing tofs
        results = DataFrame(TOF = Vector{Float64}(undef, 0), 
                            Solution = Vector{Vector}(undef, 0),
                            Δm = Vector{Float64}(undef, 0))
        stop = false
        while !stop
            tspan = (tspan[1], tspan[2] + 0.1)
            prob = IndirectOptimizationProblem("Low Thrust 10 CR3BP", ics, fcs, tspan)
            fss = FSSSolver(λ0, ics, fcs, prob.BVPFunc, prob.BVPWithSTMFunc; 
                homotopy = true, homotopyParamVec = [0.0])
            solve!(fss; factor = 3.0, ftol = 1e-10)

            # Integrate target trajectory
            solTarg = integrate(xf, (0.0, 1.6), "Low Thrust 10 CR3BP", CR3BP(), FullSolutionHistoryNoControl())
            rxTarg = [u[1] for u in solTarg.u]
            ryTarg = [u[2] for u in solTarg.u]
            rzTarg = [u[3] for u in solTarg.u]

            # Set updates
            λ0 = GetHomotopySolutionVector(fss)[1]
            ps.ϵ = 0.0

            # Initialize final integration inputs
            y0 = @SVector [ics[1], ics[2], ics[3], ics[4], ics[5], ics[6], ics[7], 
                        λ0[1], λ0[2], λ0[3], λ0[4], λ0[5], λ0[6], λ0[7]]

            # Integrate 
            tspanInteg = (tspan[1], tspan[2] * 3600*24/TU)
            sol = integrate(y0, tspanInteg, ps, CR3BP(), FullSolutionHistory(), MEMF())

            # Compute thrusting and coasting arcs
            cSc = ps.sp.c * ps.crp.TU / (ps.crp.LU * 1000.0)
            rxt = Vector{Float64}(undef, length(sol)); fill!(rxt, NaN)
            ryt = Vector{Float64}(undef, length(sol)); fill!(ryt, NaN)
            rzt = Vector{Float64}(undef, length(sol)); fill!(rzt, NaN)
            rxc = Vector{Float64}(undef, length(sol)); fill!(rxc, NaN)
            ryc = Vector{Float64}(undef, length(sol)); fill!(ryc, NaN)
            rzc = Vector{Float64}(undef, length(sol)); fill!(rzc, NaN)
            Sold = 0.0
            for i in 1:length(sol)
                λv = norm(@view(sol.u[i][11:13]))
                S = IndirectTrajOpt.computeS(sol.u[i], λv, cSc)
                #println(S)
                i == 1 ? Sold = S : ()
                if S > ps.ϵ || Sold > ps.ϵ
                    rxc[i] = sol.u[i][1]
                    ryc[i] = sol.u[i][2]
                    rzc[i] = sol.u[i][3]
                end
                if S <= ps.ϵ || Sold <= ps.ϵ
                    rxt[i] = sol.u[i][1]
                    ryt[i] = sol.u[i][2]
                    rzt[i] = sol.u[i][3]
                end
                global Sold = S
            end

            mat"""
            clf
            plot3($rxc,$ryc,$rzc, 'b')
            hold on 
            plot3($rxt,$ryt,$rzt, 'r')
            plot3($rxTarg, $ryTarg, $rzTarg, '--k')
            axis equal
            grid on
            """

            # Push results to solution vector and print results 
            Δm = (1.0 - sol.u[end][7])*1500
            push!(results, Dict("TOF" => tspan[2],
                                "Solution" => λ0,
                                "Δm" => Δm))
            println("Fuel used: " * string(Δm)) 

            # Print results to file 
            println(fid, "TOF:\t\t" * string(tspan[2]))
            lineStr = "Solution:\t["
            for i in 1:7
                lineStr *= string(λ0[i]) * " "
            end
            println(fid, lineStr * "]")
            println(fid, "Fuel Used:\t" * string(Δm))
            println(fid, "Converged:\t" * (GetHomotopyConverged(fss) ? "Yes" : "No"))
            flush(fid)
        end
    finally
        close(fid)
    end
end
main()