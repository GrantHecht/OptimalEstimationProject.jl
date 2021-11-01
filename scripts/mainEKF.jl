using OptimalEstimationProject
using LinearAlgebra
using AstroTime

# Spacecraft simulator
t0        = TAIEpoch("2021-08-08T12:00:00")
scsim     = OptimalEstimationProject.SpacecraftSim(t0.second)

# Initial state error variance
σr      = 1e-3 # [km]
σv      = 1e-5  # [km / s]
σm      = 1e-6  # [kg]
sqrtP0  = Diagonal([σr, σr, σr, σv, σv, σv, σm])
P0      = 4*(sqrtP0^2)

# Initial state estimate
(y0true, us) = OptimalEstimationProject.GetStateAndControl(scsim, t0.second)
xhat0   = y0true[1:7] + sqrtP0*randn(7) 

# Measurement time stamps
Δt      = 60 # [sec]
ts      = range(t0.second + Δt; step = Δt, stop = t0.second + scsim.ts[2]*86400)

# Measurement statistics
σρ      = 1.0e-3 # [km] Pseudorange noise standard deviation
σr      = 5.0e-3 # [km] GPS broadcast ephemeris standard deviation

# Filter measurement and process noise covariance
R       = Diagonal((σρ^2 + 3*σr^2)*ones(32)) 
Q       = Diagonal([0.0, 0.0, 0.0, 1e-1, 1e-1, 1e-2, 1e-8])

# GPS Simulation Span
startWeek   = 2170
startDay    = 0
endWeek     = 2174
endDay      = 6

# Create GPS simulation object
gpssim = GPSSim(startWeek, startDay, endWeek, endDay; σρ = σρ, σr = σr)

# Create EKF
ekf = OptimalEstimationProject.EKF(xhat0, P0, R, Q, ts, gpssim, scsim);

OptimalEstimationProject.propagate!(ekf)
OptimalEstimationProject.update!(ekf)