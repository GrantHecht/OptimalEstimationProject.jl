using DrWatson
using OptimalEstimationProject
using LinearAlgebra
using AstroTime
using DataFrames
using CSV

# Activate project environment
@quickactivate "OptimalEstimationProject"

# ===== Compute EKF MC Statistics
(ts,em,ev) = OptimalEstimationProject.computeMCStats(datadir("UKF_MC"))

# ===== Run single UKF trial
# Spacecraft simulator
t0      = TAIEpoch("2021-08-08T12:00:00")
scsim   = SpacecraftSim(t0.second)

# Initial state error variance
σr      = 0.1   # [km]
σv      = 0.01  # [km / s]
σm      = 1e-6  # [kg]
sqrtP0  = Diagonal([σr, σr, σr, σv, σv, σv, σm])
P0      = sqrtP0.^2

# Initial state estimate
(y0true, us) = OptimalEstimationProject.GetStateAndControl(scsim, t0.second)
xhat0   = y0true[1:7] + sqrtP0*randn(7) 

# Measurement time stamps
Δt      = 60.0        # [sec]
gpsΔt   = 15*60     # [sec]
ts      = range(t0.second + Δt; step = Δt, stop = t0.second + scsim.ts[2]*86400)

# Measurement statistics
σρ      = 10.0e-3 # [km]     Pseudorange noise standard deviation
σr      = 5.0e-3  # [km]     GPS broadcast ephemeris standard deviation
σa      = 1.0e-6  # [km/s^2] Accelerometer noise standard deviation 

# Process noise covariance
Q       = Diagonal([1.0e-18, 1.0e-18, 1.0e-18, 1.0e-12, 1.0e-12, 1.0e-12, 5.0e-3])

# GPS Simulation Span
startWeek   = 2170
startDay    = 0
endWeek     = 2174
endDay      = 6

# Create GPS simulation object
gpssim = GPSSim(startWeek, startDay, endWeek, endDay; σρ = σρ, σr = σr)

# Create IMU simulation object
imusim = IMUSim(σa)

# Create EKF
ukf = UKF(xhat0, P0, Q, (σρ^2 + 3*σr^2), σa^2, ts, gpsΔt, gpssim, imusim, scsim; 
    lunaPerts = true, α = 1.0, β = 0.0, κ = 3.0+14.0);

runFilter!(ukf)

# ===== MC Plots
plotMC(ukf, em, ev, length(ukf.txp))