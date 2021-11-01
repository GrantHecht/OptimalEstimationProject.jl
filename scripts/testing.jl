
using OptimalEstimationProject
using AstroTime

# GPS Simulation Span
startWeek   = 2170
startDay    = 0
endWeek     = 2174
endDay      = 6

# Create GPS simulation object
gpssim = GPSSim(startWeek, startDay, endWeek, endDay)

# Generate measurement
measEpoch = TAIEpoch("2021-08-25T12:00:00")
σρ        = 0.1/1000;
r         = [6700.0, 0.0, 0.0]
mSecTIA   = Float64(measEpoch.second)
meas = OptimalEstimationProject.gpsMeasurement(gpssim, r, mSecTIA, σρ)

# Spacecraft simulator
t0        = TAIEpoch("2021-08-08T12:00:00")
scsim     = OptimalEstimationProject.SpacecraftSim(t0.second)

# Get state and control
(xs,us)   = OptimalEstimationProject.GetStateAndControl(scsim, t0.second .+ [t for t in 0:10:1000])

# Plotting
OptimalEstimationProject.PlotTrajectory(scsim; frame = :inertial)
#OptimalEstimationProject.PlotTrajectory(scsim, range(t0.second + 400; stop = t0.second + 800, length = 100))