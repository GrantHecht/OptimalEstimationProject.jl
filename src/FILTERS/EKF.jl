
mutable struct EKF
    # Measurement times
    ts::Vector{Float64}

    # Measurement residuals
    rs::Matrix{Float64}

    # State estimates and covariance matrix diags
    txp::Vector{Float64} # Epochs accosiated with entries in xhats and Ps
    ixp::Int # Index associated with latest entry
    xhats::Matrix{Float64}
    es::Matrix{Float64}
    Ps::Matrix{Float64}

    # Preallocate matricies
    P⁺::Matrix{Float64}
    P⁻::Matrix{Float64}
    H::Matrix{Float64}
    HPH::Matrix{Float64}
    PH::Matrix{Float64}
    KH::Matrix{Float64}
    F::Matrix{Float64}
    R::Diagonal{Float64, Vector{Float64}}
    Q::Diagonal{Float64, Vector{Float64}}
    y::Vector{Float64}

    # Simulators
    gpsSim::GPSSim
    scSim::SpacecraftSim

    # Steps to save per propagation
    stp::Int64

    # Step Index
    k::Int64

    # Update flag
    updated::Bool

    # Finished flag
    finished::Bool
end

function EKF(xhat0, P0, R, Q, ts, gpsSim, scSim; steps2save = 20)
    # Instantiate preallocated matriciesAdd
    ixp     = 1
    rs      = zeros(length(ts), 32)
    txp     = zeros((steps2save + 1)*length(ts) + 1); txp[1] = scSim.t0
    xhats   = zeros((steps2save + 1)*length(ts) + 1, 7); xhats[1, :] .= xhat0
    es      = zeros((steps2save + 1)*length(ts) + 1, 7); xhats[1, :]
    Ps      = zeros((steps2save + 1)*length(ts) + 1, 7); Ps[1, :] .= diag(P0)
    P⁻      = zeros(7,7)
    P⁺      = deepcopy(P0)
    H       = zeros(32, 7)
    HPH     = zeros(32, 32)
    PH      = zeros(7, 32)
    KH      = zeros(7,7)
    F       = zeros(7,7)
    y       = zeros(7 + 49)

    # Compute initial estimates error
    (ytrue, u) = GetStateAndControl(scSim, scSim.t0)
    @views es[1,:] .= xhat0 .- ytrue[1:7]

    # Instantiate EKF
    EKF(ts,rs,txp,ixp,xhats,es,Ps,P⁺,P⁻,H,HPH,PH,KH,F,R,Q,y,gpsSim,scSim,steps2save,0,true,false)
end

function propagate!(ekf::EKF)
    if ekf.updated == false
        throw(ErrorException("Cannot propagate filter before updating."))
    else
        ekf.k += 1
        ekf.updated = false
    end

    # Get time 
    if ekf.k == 1
        t0 = ekf.scSim.t0
        tf = ekf.ts[1]
    else
        t0 = ekf.ts[ekf.k - 1]
        tf = ekf.ts[ekf.k]
    end

    # Generate vector of times to save at 
    saveat = range(t0; length = ekf.stp, stop = tf)

    # Construct integration vector
    if ekf.k == 1
        @views ekf.y[1:7]   .= ekf.xhats[1, :]
    else
        @views ekf.y[1:7]   .= ekf.xhats[ekf.ixp, :]
    end
    pVec = reshape(ekf.P⁺, (49, 1))
    for i in 1:49
        ekf.y[7 + i] = pVec[i]
    end

    # Create ODE Problem 
    prob = ODEProblem(ekfNoLunarPertEOM!, ekf.y, (t0, tf), ekf)

    # Solve ODE 
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat=saveat)

    # Set P⁻
    ekf.P⁻ .= reshape(sol.u[end][8:end], (7, 7))

    # Update storage matricies
    start = ekf.ixp + 1
    stop  = ekf.ixp + ekf.stp
    step  = 0
    for i in start:stop
        step += 1

        # Time 
        ekf.txp[i] = sol.t[step]

        # State
        @views ekf.xhats[i,:] .= sol.u[step][1:7]

        # Error
        (ytrue, u) = GetStateAndControl(ekf.scSim, ekf.scSim.t0)
        @views ekf.es[i,:] .= ekf.xhats[i,:] .- ytrue[1:7]

        # Covariance
        @views ekf.Ps[i,:] .= diag(reshape(sol.u[step][8:end], (7, 7)))
        ekf.ixp = i
    end

    return nothing
end

function update!(ekf::EKF)
    if ekf.updated == true
        throw(ErrorException("Cannot update before propagating."))
    else
        ekf.updated = true
    end

    # Get current time
    t = ekf.ts[ekf.k]

    # Get state estimate
    @views xhat = ekf.xhats[ekf.ixp, :]

    # Get true state 
    (ytrue, u) = GetStateAndControl(ekf.scSim, t)

    # Generate measurements
    @views meas = gpsMeasurement(ekf.gpsSim, ytrue[1:3], t)

    # Get satelite ids
    @views sats = meas[2:2:end]
    numSats = length(sats)

    # Get expected measurements
    @views (expMeas, ps) = gpsMeasurement(ekf.gpsSim, xhat[1:3], sats, t; type = :expected)

    # Compute residual
    @views ekf.rs[ekf.k, 1:numSats] .= meas[3:2:end] .- expMeas[3:2:end]

    # Compute measurement Jacobian
    ekf.H .*= 0.0
    for i in 1:numSats
        @views ekf.H[i, 1:3] .= ps[i,1:3]
        @views ekf.H[i, 1:3] ./= norm(ps[i,1:3])
    end

    # Compute gain
    @views mul!(ekf.PH[:,1:numSats], ekf.P⁻, transpose(ekf.H[1:numSats, :]))
    @views mul!(ekf.HPH[1:numSats,1:numSats], ekf.H[1:numSats, :], ekf.PH[:,1:numSats])
    display(ekf.HPH)
    for i in 1:numSats; ekf.HPH[i,i] += ekf.R[i,i]; end
    @views HPHpRFact = factorize(ekf.HPH[1:numSats,1:numSats])
    @views rdiv!(ekf.PH[:,1:numSats], HPHpRFact) # Stores K in PH matrix
    display(ekf.PH)

    # Compute state update
    ekf.ixp += 1
    @views mul!(ekf.xhats[ekf.ixp, :], ekf.PH[:,1:numSats], ekf.rs[ekf.k, 1:numSats])
    display(ekf.xhats[ekf.ixp,:])
    display(xhat)
    ekf.xhats[ekf.ixp, :] .+= xhat

    # Set time in storage vector
    ekf.txp[ekf.ixp] = t

    # Compute error
    (ytrue, u) = GetStateAndControl(ekf.scSim, ekf.scSim.t0)
    @views ekf.es[ekf.ixp, :] .= ekf.xhats[ekf.ixp, :] .- ytrue[1:7]

    # Compute covariance update
    @views mul!(ekf.KH, ekf.PH[:,1:numSats], ekf.H[1:numSats,:], -1.0, 0.0)
    for i in 1:7; ekf.KH[i,i] += 1.0; end
    mul!(ekf.P⁺, ekf.KH, ekf.P⁻)
    ekf.Ps[ekf.ixp, :] .= diag(ekf.P⁺)

    return nothing
end

function runFilter(ekf::EKF)
    for i in 1:length(ekf.ts)
        propagate!(ekf)
        update!(ekf)
    end
end

function ekfNoLunarPertEOM!(dy, y, ekf::EKF, t)
    # Earth gravitational parameter
    μ = 3.986004418e5 # [km^3/s^2]

    # Thrust magnitued
    T = 10 # [N]

    # Get control
    u = GetControl(ekf.scSim, t)

    # States
    @views r        = norm(y[1:3])
    @views dy[1:3] .= y[4:6]
    @views dy[4:6] .= -(μ/r^3).*y[1:3] .+ (T*u[1] / y[7]).*u[2:4]
    dy[7]           = -u[1]*T / ekf.scSim.ps.sp.c

    # Dynamics Jacobian
    # dr/dv
    ekf.F[1,4] = 1.0
    ekf.F[2,5] = 1.0
    ekf.F[3,6] = 1.0

    # dv/dr
    @views mul!(ekf.F[4:6, 1:3], y[1:3], transpose(y[1:3]))
    @views ekf.F[4:6, 1:3] .*= (3.0 / r^2)
    for i in 1:3; ekf.F[3 + i, i] -= 1.0; end
    @views ekf.F[4:6, 1:3] .*= (μ/r^3)

    # dv/dm
    @views ekf.F[4:6, 7] .= -(u[1]*T / y[7]^2) .* u[2:4]

    # Covariance Dynamics
    @views P   = reshape(y[8:end], (7, 7))
    @views dP  = reshape(dy[8:end], (7, 7))
    mul!(dP, ekf.F, P)
    mul!(dP, P, transpose(ekf.F), 1.0, 1.0)
    dP .+= ekf.Q
end