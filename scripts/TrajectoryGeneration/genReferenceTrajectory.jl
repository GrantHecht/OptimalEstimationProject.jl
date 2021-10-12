
using DrWatson
@quickactivate "OptimalEstimationProject"
using IndirectTrajOpt

# This script attemplts to find an optimal minimum fuel trajectory with initial 
# and final conditions ics and fcs respectively. It does so with an indirect 
# approaching, first searching for a good selection for the initial co-states
# corresponding to the minimum energy problem. A homotopy continuation is 
# then performed to transition PSO initialized co-states from the minumum 
# energy problem to that of mimimum fuel.

# Ref: Hecht, G., Botta, E. 2021. "Co-State Initialization with Particle Swarm 
# Optimization for Low-Thrust Mimumum-Fuel Trajectory Optimization," AAS/AIAA
# Astrodynamic Specialist Conference.

# IndirectTrajOpt.jl Julia package and its dependancies can be found at:
#   IndirectTrajOpt.jl      - https://github.com/GrantHecht/IndirectTrajOpt.jl
#   IndirectCoStateInit.jl  - https://github.com/GrantHecht/IndirectCoStateInit.jl
#   IndirectShooting.jl     - https://github.com/GrantHecht/IndirectShooting.jl
#   Heuristics.jl           - https://github.com/GrantHecht/Heuristics.jl

function genRefTrajs(swarmSize, numTrials, tf)
    # Set tspan in days
    tspan = (0.0, tf) # [Days]

    # Grab parameters for scaling
    pTemp = IndirectTrajOpt.initCR3BPIndirectParams("Low Thrust 10 CR3BP")
    μ  = pTemp.crp.μ
    LU = pTemp.crp.LU 
    TU = pTemp.crp.TU

    # Set initial and final conditions 
    ics = [-0.019488511458668, -0.016033479812051, 0.0,
            8.918881923678198, -4.081793688818725, 0.0,
            1.0] 
    fcs = [(1.0 - μ) + 6060.483/LU, 19452.284/LU, -34982.968/LU, 0.082677*(TU/LU), 0.006820*(TU/LU), -0.368434*(TU/LU), 0.0]

    # Homotopy continuation 
    N = 20
    hVec = [(j^2 - 1.)/(N^2 - 1.) for j in N:-1:1]

    # Instantiate indirect optimization problem 
    prob = IndirectOptimizationProblem("Low Thrust 10 CR3BP", ics, fcs, tspan)

    # Create vector of co-state initializers
    itoVec = [IndirectTrajOptimizer(
        prob;
        initOptimizer    = :PSO,
        initCostFunc     = :WSS,
        numParticles     = swarmSize,
        numSwarms        = 8,
        UBs              = [20., 20., 5., 0.1, 0.1, 0.01, 1.0],
        LBs              = [-20., -20., -5., -0.1, -0.1, -0.01, 0.0],
        iUBs             = [1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1.0],
        iLBs             = [-1e-6, -1e-6, -1e-6, -1e-6, -1e-6, -1e-6, 1e-6],
        weights          = [10, 10, 10, 1, 1, 1, 1],
        MFD              = 1.0,
        displayInterval  = 5,
        funcTol          = 1e-5,
        minNeighborhoodFraction = 0.05,
        maxStallIters    = 50,
        maxIters         = 5000,
        maxTime          = 1800.0,
        homotopyParamVec = hVec,
        dataFolder       = datadir("numParticles"*string(swarmSize)*"_tf"*replace(string(tf), '.' => 'p')),
        writeData        = true) for i in 1:numTrials];
    
    # Initialize all co-states 
    for i in 1:numTrials
        initialize!(itoVec[i])
        solve!(itoVec[i], factor = 3.0, convergenceAttempts = 4, ftol = 1e-11)
    end
end
