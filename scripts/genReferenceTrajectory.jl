
using DrWatson
@quickactivate "OptimalEstimationProject"
using IndirectTrajOpt

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
   prob = IndirectOptimizationProblem("Low Thrust 10 CR3BP", ics, fcs, tspan, 
       initializationFlag = Initialization())

   # Create vector of co-state initializers
   itoVec = [IndirectTrajOptimizer(
       prob;
       initOptimizer    = :PSO,
       initCostFunc     = :WSS,
       numParticles     = swarmSize,
       numSwarms        = 8,
       UBs              = [150, 150, 20, 20, 20, 20, 20],
       LBs              = [-150, -150, -20, -20, -20, -20, 0],
       #iUBs             = [40, 40, 2, 2, 2, 2, 2],
       #iLBs             = [-40, -40, -2, -2, -2, -2, 0],
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
