using DrWatson
using OptimalEstimationProject
using LinearAlgebra
using AstroTime
using DataFrames
using CSV

# Number of trials
nTrials = 10

function main(nTrials)
    # Spacecraft simulator
    t0      = TAIEpoch("2021-08-08T12:00:00")
    scsim   = SpacecraftSim(t0.second)

    # Initial state error variance
    σr      = 0.1   # [km]
    σv      = 0.01  # [km / s]
    σm      = 1e-6  # [kg]
    sqrtP0  = Diagonal([σr, σr, σr, σv, σv, σv, σm])
    P0      = sqrtP0.^2

    # Initial true state
    (y0true, us) = OptimalEstimationProject.GetStateAndControl(scsim, t0.second)

    # Measurement time stamps
    Δt      = 30        # [sec]
    gpsΔt   = 15*60     # [sec]
    ts      = range(t0.second + Δt; step = Δt, stop = t0.second + scsim.ts[2]*86400)

    # Measurement statistics
    σρ      = 1.0e-3 # [km]     Pseudorange noise standard deviation
    σr      = 5.0e-3 # [km]     GPS broadcast ephemeris standard deviation
    σa      = 1.0e-5 # [km/s^2] Accelerometer noise standard deviation 

    # Process noise covariance
    Q       = Diagonal([0.0, 0.0, 0.0, 1e-10, 1e-10, 1e-10, 5e-4])

    # GPS Simulation Span
    startWeek   = 2170
    startDay    = 0
    endWeek     = 2174
    endDay      = 6

    # Create GPS simulation object
    gpssim = GPSSim(startWeek, startDay, endWeek, endDay; σρ = σρ, σr = σr)

    # Create IMU simulation object
    imusim = IMUSim(σa)

    # MC Loop
    Threads.@threads for i in 1:nTrials
        # Initial state estimate
        xhat0   = y0true[1:7] + sqrtP0*randn(7) 

        # Create EKF
        ekf = EKF(xhat0, P0, Q, (σρ^2 + 3*σr^2), σa^2, ts, gpsΔt, 
            deepcopy(gpssim), deepcopy(imusim), deepcopy(scsim); steps2save = 2, lunaPerts = true);

        # Run filter
        runFilter!(ekf)

        # Write results to file (Just looking at filter uncertainty quantification so only need time and error)
        df = DataFrame(
            Time    = ekf.txp,
            erx     = ekf.es[:,1],
            ery     = ekf.es[:,2],
            erz     = ekf.es[:,3],
            evx     = ekf.es[:,4],
            evy     = ekf.es[:,5],
            evz     = ekf.es[:,6],
            em      = ekf.es[:,7])

        CSV.write(datadir("EKF_MC", "trial_" * string(i) * ".csv"), df)
    end

end

main(nTrials)