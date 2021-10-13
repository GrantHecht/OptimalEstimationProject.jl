
using DrWatson
@quickactivate "OptimalEstimationProject"
using MATLAB
using IndirectTrajOpt
using StaticArrays
using LinearAlgebra

function main()
# Data info
tf      = 13.0
folder  = "numParticles2000_tf"*replace(string(tf), '.' => 'p')
trial   = 10
hStep   = 2

# Load data
df = readTextData(datadir(folder))

# Init parameters
ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
μ = ps.crp.μ
LU = ps.crp.LU
TU = ps.crp.TU

# Trajectory time span
tspan = (0.0, tf * 86400 / TU)

# Initial conditions
ics = [-0.019488511458668, -0.016033479812051, 0.0,
        8.918881923678198, -4.081793688818725, 0.0,
        1.0] 

# Target final conditions
xf = @SVector [(1.0 - μ) + 6060.483/LU, 19452.284/LU, -34982.968/LU, 0.082677*(TU/LU), 0.006820*(TU/LU), -0.368434*(TU/LU)]

# Integrate target trajectory
solTarg = integrate(xf, (0.0, 1.6), "Low Thrust 10 CR3BP", CR3BP(), FullSolutionHistoryNoControl())
rxTarg = [u[1] for u in solTarg.u]
ryTarg = [u[2] for u in solTarg.u]
rzTarg = [u[3] for u in solTarg.u]

# Grab requirements from DataFrame
λ0 = df.HomotopySolutions[trial][hStep]
ps.ϵ = (hStep == 1 ? 1.0 : df.HomotopyParams[trial][hStep - 1])

# Initialize final integration inputs
y0 = @SVector [ics[1], ics[2], ics[3], ics[4], ics[5], ics[6], ics[7], 
            λ0[1], λ0[2], λ0[3], λ0[4], λ0[5], λ0[6], λ0[7]]

# Integrate 
sol = integrate(y0, tspan, ps, CR3BP(), FullSolutionHistory(), MEMF())

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
    Sold = S
end

mat"""
figure
plot3($rxc,$ryc,$rzc, 'b')
hold on 
plot3($rxt,$ryt,$rzt, 'r')
plot3($rxTarg, $ryTarg, $rzTarg, '--k')
axis equal
grid on
"""
end
main()