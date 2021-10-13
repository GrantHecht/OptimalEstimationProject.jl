using DrWatson
@quickactivate "OptimalEstimationProject"
using MATLAB
using IndirectTrajOpt
using StaticArrays
using LinearAlgebra

# Data info
tf      = 33.1

# Init parameters
ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
μ = ps.crp.μ
LU = ps.crp.LU
TU = ps.crp.TU

# Trajectory time span
tspan = (0.0, tf)
tspanPlot = (0.0, tspan[2] * 86400 / TU)

# Initial conditions
ics = [-0.019488511458668, -0.016033479812051, 0.0,
        8.918881923678198, -4.081793688818725, 0.0,
        1.0] 

# Target final conditions
fcs = [(1.0 - μ) + 6060.483/LU, 19452.284/LU, -34982.968/LU, 0.082677*(TU/LU), 0.006820*(TU/LU), -0.368434*(TU/LU), 0.0]
xf = @SVector [fcs[1], fcs[2], fcs[3], fcs[4], fcs[5], fcs[6]]

# Integrate target trajectory
solTarg = integrate(xf, (0.0, 1.6), "Low Thrust 10 CR3BP", CR3BP(), FullSolutionHistoryNoControl())
rxTarg = [u[1] for u in solTarg.u]
ryTarg = [u[2] for u in solTarg.u]
rzTarg = [u[3] for u in solTarg.u]

#for i in 1:6
#    fcs[i] = solTarg(0.16)[i]
#end

# Grab requirements from DataFrame
#λ0 = [3.0872240069228414, 7.585751775356703, -0.8828499338288641, -0.022292151740247245, 0.01139245771460113, -0.00015851392516777712, 0.13268783730938283]
#ps.ϵ = 0.16662720666508826
#λ0 = [5.963713007087376, 12.99948134789672, 1.6584130970338573, -0.04062357611978672, 0.018558730103857064, -0.0004993971529679529, 0.1220053310292062] 
λ0 = [5.963713007227692, 12.999481348176722, 1.6584130970663207, -0.04062357612068792, 0.018558730104232632, -0.0004993971529784666, 0.12200533102974148]
ps.ϵ = 0.0
#λ0 = GetHomotopySolutionVector(fss)[1]

prob = IndirectOptimizationProblem("Low Thrust 10 CR3BP", ics, fcs, tspan)
fss = FSSSolver(λ0, ics, fcs, prob.BVPFunc, prob.BVPWithSTMFunc;
                homotopy = true, homotopyParamVec = [0.0])
#solve!(fss; factor = 3.0, ftol = 1e-10)
#λ0 = GetHomotopySolutionVector(fss)[1]

# Initialize final integration inputs
y0 = @SVector [ics[1], ics[2], ics[3], ics[4], ics[5], ics[6], ics[7], 
            λ0[1], λ0[2], λ0[3], λ0[4], λ0[5], λ0[6], λ0[7]]

# Integrate 
sol = integrate(y0, tspanPlot, ps, CR3BP(), FullSolutionHistory(), MEMF())

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
figure
plot3(-$μ*$LU, 0, 0, 'xk')
hold on
plot3((1-$μ)*$LU, 0, 0, 'xk')
plot3($rxc*$LU,$ryc*$LU,$rzc*$LU, 'b')
plot3($rxt*$LU,$ryt*$LU,$rzt*$LU, 'r')
plot3($rxTarg*$LU, $ryTarg*$LU, $rzTarg*$LU, '--k')
axis equal
grid on
xlabel("x, km")
ylabel("y, km")
zlabel("z, km")
"""