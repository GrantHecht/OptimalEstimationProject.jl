using DrWatson
@quickactivate "OptimalEstimationProject"
using MATLAB
using IndirectTrajOpt
using StaticArrays
using LinearAlgebra
using DataFrames

# Trajectory time span
tspan = (0.0, 26.0)

# ME Initial CoStates
λ0 = [-0.8778526555553886,-1.5568310993815173,-0.4120913077312581,0.005878876387420351,-0.002181056020767603,-0.00035144027756129617,0.0941110364795588]

# Initial conditions
ics = [-0.019488511458668, -0.016033479812051, 0.0,
        8.918881923678198, -4.081793688818725, 0.0,
        1.0] 

# Init parameters
ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
μ = ps.crp.μ
LU = ps.crp.LU
TU = ps.crp.TU

# Target final conditions
fcs = [(1.0 - μ) + 6060.483/LU, 19452.284/LU, -34982.968/LU, 0.082677*(TU/LU), 0.006820*(TU/LU), -0.368434*(TU/LU), 0.0]
xf = @SVector [fcs[1], fcs[2], fcs[3], fcs[4], fcs[5], fcs[6]]

# Solve with homotopy
N = 5000
prob = IndirectOptimizationProblem("Low Thrust 10 CR3BP", ics, fcs, tspan)
fss = FSSSolver(λ0, ics, fcs, prob.BVPFunc, prob.BVPWithSTMFunc;
    homotopy = true, homotopyParamVec = [(j^2 - 1.)/(N^2 - 1.) for j in N:-1:1])
solve!(fss; factor = 1.0, ftol = 1e-10)