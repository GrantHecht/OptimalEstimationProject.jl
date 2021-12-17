
mutable struct EKF
    # Measurement times
    ts::Vector{Float64}
    imuΔt::Float64
    gpsΔt::Float64

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
    Q::Diagonal{Float64, Vector{Float64}}
    y::Vector{Float64}
    dy::Vector{Float64}

    # Simulators
    gpsSim::GPSSim
    imuSim::IMUSim
    scSim::SpacecraftSim

    # Measurement covariance terms
    σ2GPS::Float64
    σ2Acc::Float64

    # Steps to save per propagation
    stp::Int64

    # Step Index
    k::Int64

    # Update flag
    updated::Bool

    # Finished flag
    finished::Bool

    # Luna perturbation flag
    lunaPerts::Bool
end

function EKF(xhat0, P0, Q, σ2GPS, σ2Acc, ts, gpsΔt, gpsSim, imuSim, scSim; lunaPerts = false, steps2save = 5)
    # Instantiate preallocated matriciesAdd
    ixp     = 1
    rs      = zeros(length(ts), 32 + 3)
    txp     = zeros((steps2save + 1)*length(ts) + 1); txp[1] = scSim.t0
    xhats   = zeros((steps2save + 1)*length(ts) + 1, 7); xhats[1, :] .= xhat0
    es      = zeros((steps2save + 1)*length(ts) + 1, 7);
    Ps      = zeros((steps2save + 1)*length(ts) + 1, 7); Ps[1, :] .= diag(P0)
    P⁻      = zeros(7,7)
    P⁺      = deepcopy(P0)
    H       = zeros(32 + 3, 7)
    HPH     = zeros(32 + 3, 32 + 3)
    PH      = zeros(7, 32 + 3)
    KH      = zeros(7,7)
    F       = zeros(7,7)
    y       = zeros(7 + 49)
    dy      = deepcopy(y)

    # Compute IMU Δt from ts vector
    imuΔt   = ts[2] - ts[1]

    # Compute initial estimates error
    (ytrue, u) = GetStateAndControl(scSim, scSim.t0)
    @views es[1,:] .= xhat0 .- ytrue[1:7]

    # Instantiate EKF
    EKF(ts,imuΔt,gpsΔt,rs,txp,ixp,xhats,es,Ps,P⁺,P⁻,H,HPH,PH,KH,F,Q,y,dy,gpsSim,imuSim,scSim,σ2GPS,σ2Acc,steps2save,0,true,false,lunaPerts)
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
    if ekf.lunaPerts == false
        prob = ODEProblem(ekfNoLunarPertEOM!, ekf.y, (t0, tf), ekf)
    else
        prob = ODEProblem(ekfWithLunarPertEOM!, ekf.y, (t0, tf), ekf)
    end

    # Solve ODE 
    sol = solve(prob, Vern7(), reltol=1e-10, abstol=1e-10, saveat=saveat)

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
        (ytrue, u) = GetStateAndControl(ekf.scSim, ekf.txp[i])
        @views ekf.es[i,:] .= ekf.xhats[i,:] .- ytrue[1:7]

        # Covariance
        @views ekf.Ps[i,:] .= diag(reshape(sol.u[step][8:end], (7, 7)))
        ekf.ixp = i
    end

    return nothing
end

function updateGPS!(ekf::EKF)
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

    # Generate gps measurements
    @views measGPS = gpsMeasurement(ekf.gpsSim, ytrue[1:3], t)

    # Get satelite ids
    @views sats = measGPS[2:2:end]
    numSats = length(sats)

    # Get expected gps measurements
    @views (expMeasGPS, ps) = gpsMeasurement(ekf.gpsSim, xhat[1:3], sats, t; type = :expected)

    # Compute residuals and measurement partials
    numRejected = 0
    rejectTol   = 10.0
    for i in 1:numSats
        r = measGPS[3 + 2*(i - 1)] - expMeasGPS[3 + 2*(i - 1)]
        if abs(r) < rejectTol
            # Compute residual while catching bad measurements
            @views ekf.rs[ekf.k, i - numRejected] = r

            # Compute measurement Jacobian
            ekf.H .*= 0.0
            @views ekf.H[i - numRejected, 1:3] .= -ps[i,1:3]
            @views ekf.H[i - numRejected, 1:3] ./= norm(ps[i,1:3])
        else
            numRejected += 1
        end
    end
    numSats -= numRejected
    
    # Compute gain
    numMeas = numSats
    @views mul!(ekf.PH[:,1:numMeas], ekf.P⁻, transpose(ekf.H[1:numMeas, :]))
    @views mul!(ekf.HPH[1:numMeas,1:numMeas], ekf.H[1:numMeas, :], ekf.PH[:,1:numMeas])
    for i in 1:numSats; ekf.HPH[i,i] += ekf.σ2GPS; end
    @views HPHpRFact = factorize(ekf.HPH[1:numMeas,1:numMeas])
    @views rdiv!(ekf.PH[:,1:numMeas], HPHpRFact) # Stores K in PH matrix

    # Compute state update
    ekf.ixp += 1
    @views mul!(ekf.xhats[ekf.ixp, :], ekf.PH[:,1:numMeas], ekf.rs[ekf.k, 1:numMeas])
    ekf.xhats[ekf.ixp, :] .+= ekf.xhats[ekf.ixp - 1, :]

    # Set time in storage vector
    ekf.txp[ekf.ixp] = t

    # Compute error
    (ytrue, u) = GetStateAndControl(ekf.scSim, ekf.txp[ekf.ixp])
    @views ekf.es[ekf.ixp, :] .= ekf.xhats[ekf.ixp, :] .- ytrue[1:7]

    # Compute covariance update
    @views mul!(ekf.KH, ekf.PH[:,1:numMeas], ekf.H[1:numMeas,:], -1.0, 0.0)
    for i in 1:7; ekf.KH[i,i] += 1.0; end
    mul!(ekf.P⁺, ekf.KH, ekf.P⁻)

    ekf.Ps[ekf.ixp, :] .= diag(ekf.P⁺)

    return nothing
end

function updateIMU!(ekf::EKF)
    if ekf.updated == true
        throw(ErrorException("Cannot update before propagating."))
    else
        ekf.updated = true
    end

    # Get current time
    t = ekf.ts[ekf.k]

    # Get state estimate
    @views xhat = ekf.xhats[ekf.ixp, :]

    # Generate accelerometer measurements
    (dy, _us) = GetStateAndControl(ekf.scSim, t; derivs = true)
    @views measAcc = ComputeMeasurement(ekf.imuSim, dy[4:6])

    # Compute expected accelerometer measurement
    @views ekf.y[1:7] .= xhat
    ekf.y[8:end] .= 0.0
    for i in 1:7; ekf.y[7 + 7*(i - 1) + i] = 1.0; end
    if ekf.lunaPerts == false
        ekfNoLunarPertEOM!(ekf.dy, ekf.y, ekf, t)
    else
        ekfWithLunarPertEOM!(ekf.dy, ekf.y, ekf, t)
    end   
    @views expMeasAcc = ComputeMeasurement(ekf.imuSim, ekf.dy[4:6]; type = :expected)

    # Compute accelerometer residuals and partials
    ekf.rs[ekf.k, 1:3]   .= measAcc .- expMeasAcc
    @views ekf.H[1:3, 1:3] .= ekf.F[4:6, 1:3]
    @views ekf.H[1:3, 7] .= ekf.F[4:6, 7]

    # Compute gain
    numMeas = 3
    @views mul!(ekf.PH[:,1:numMeas], ekf.P⁻, transpose(ekf.H[1:numMeas, :]))
    @views mul!(ekf.HPH[1:numMeas,1:numMeas], ekf.H[1:numMeas, :], ekf.PH[:,1:numMeas])
    for i in 1:3; ekf.HPH[i,i] += ekf.σ2Acc; end
    @views HPHpRFact = factorize(ekf.HPH[1:numMeas,1:numMeas])
    @views rdiv!(ekf.PH[:,1:numMeas], HPHpRFact) # Stores K in PH matrix

    # Compute state update
    ekf.ixp += 1
    @views mul!(ekf.xhats[ekf.ixp, :], ekf.PH[:,1:numMeas], ekf.rs[ekf.k, 1:numMeas])
    ekf.xhats[ekf.ixp, :] .+= ekf.xhats[ekf.ixp - 1, :]

    # Set time in storage vector
    ekf.txp[ekf.ixp] = t

    # Compute error
    (ytrue, u) = GetStateAndControl(ekf.scSim, ekf.txp[ekf.ixp])
    @views ekf.es[ekf.ixp, :] .= ekf.xhats[ekf.ixp, :] .- ytrue[1:7]

    # Compute covariance update
    @views mul!(ekf.KH, ekf.PH[:,1:numMeas], ekf.H[1:numMeas,:], -1.0, 0.0)
    for i in 1:7; ekf.KH[i,i] += 1.0; end
    mul!(ekf.P⁺, ekf.KH, ekf.P⁻)

    ekf.Ps[ekf.ixp, :] .= diag(ekf.P⁺)

    return nothing
end


function updateGPSIMU!(ekf::EKF)
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

    # Generate gps measurements
    @views measGPS = gpsMeasurement(ekf.gpsSim, ytrue[1:3], t)

    # Get satelite ids
    @views sats = measGPS[2:2:end]
    numSats = length(sats)

    # Get expected gps measurements
    @views (expMeasGPS, ps) = gpsMeasurement(ekf.gpsSim, xhat[1:3], sats, t; type = :expected)

    # Compute residuals and measurement partials
    numRejected = 0
    rejectTol   = 1.0
    for i in 1:numSats
        r = measGPS[3 + 2*(i - 1)] - expMeasGPS[3 + 2*(i - 1)]
        if abs(r) < rejectTol
            # Compute residual while catching bad measurements
            @views ekf.rs[ekf.k, i - numRejected] = r

            # Compute measurement Jacobian
            ekf.H .*= 0.0
            @views ekf.H[i - numRejected, 1:3] .= -ps[i,1:3]
            @views ekf.H[i - numRejected, 1:3] ./= norm(ps[i,1:3])
        else
            numRejected += 1
        end
    end
    numSats -= numRejected

    # Generate accelerometer measurements
    (dy, _us) = GetStateAndControl(ekf.scSim, t; derivs = true)
    @views measAcc = ComputeMeasurement(ekf.imuSim, dy[4:6])

    # Compute expected accelerometer measurement
    @views ekf.y[1:7] .= xhat
    ekf.y[8:end] .= 0.0
    for i in 1:7; ekf.y[7 + 7*(i - 1) + i] = 1.0; end
    if ekf.lunaPerts == false
        ekfNoLunarPertEOM!(ekf.dy, ekf.y, ekf, t)
    else
        ekfWithLunarPertEOM!(ekf.dy, ekf.y, ekf, t)
    end   
    @views expMeasAcc = ComputeMeasurement(ekf.imuSim, ekf.dy[4:6]; type = :expected)

    # Compute accelerometer residuals and partials
    ekf.rs[ekf.k, numSats + 1:numSats + 3]   .= measAcc .- expMeasAcc
    @views ekf.H[numSats + 1:numSats + 3, 1:3] .= ekf.F[4:6, 1:3]
    @views ekf.H[numSats + 1:numSats + 3, 7] .= ekf.F[4:6, 7]

    # Compute gain
    numMeas = numSats + 3
    @views mul!(ekf.PH[:,1:numMeas], ekf.P⁻, transpose(ekf.H[1:numMeas, :]))
    @views mul!(ekf.HPH[1:numMeas,1:numMeas], ekf.H[1:numMeas, :], ekf.PH[:,1:numMeas])
    for i in 1:numSats; ekf.HPH[i,i] += ekf.σ2GPS; end
    for i in numSats + 1:numMeas; ekf.HPH[i,i] += ekf.σ2Acc; end
    @views HPHpRFact = factorize(ekf.HPH[1:numMeas,1:numMeas])
    @views rdiv!(ekf.PH[:,1:numMeas], HPHpRFact) # Stores K in PH matrix

    # Compute state update
    ekf.ixp += 1
    @views mul!(ekf.xhats[ekf.ixp, :], ekf.PH[:,1:numMeas], ekf.rs[ekf.k, 1:numMeas])
    ekf.xhats[ekf.ixp, :] .+= ekf.xhats[ekf.ixp - 1, :]

    # Set time in storage vector
    ekf.txp[ekf.ixp] = t

    # Compute error
    (ytrue, u) = GetStateAndControl(ekf.scSim, ekf.txp[ekf.ixp])
    @views ekf.es[ekf.ixp, :] .= ekf.xhats[ekf.ixp, :] .- ytrue[1:7]

    # Compute covariance update
    @views mul!(ekf.KH, ekf.PH[:,1:numMeas], ekf.H[1:numMeas,:], -1.0, 0.0)
    for i in 1:7; ekf.KH[i,i] += 1.0; end
    mul!(ekf.P⁺, ekf.KH, ekf.P⁻)

    ekf.Ps[ekf.ixp, :] .= diag(ekf.P⁺)

    return nothing
end

function runFilter!(ekf::EKF)
    for i in 1:length(ekf.ts)
        propagate!(ekf)
        if i*ekf.imuΔt % ekf.gpsΔt == 0
            updateGPSIMU!(ekf)
        else
            updateIMU!(ekf)
        end
    end
end

function ekfNoLunarPertEOM!(dy, y, ekf::EKF, t)
    # Zero F 
    ekf.F .= 0.0

    # Earth gravitational parameter
    μ = 3.986004418e5 # [km^3/s^2]

    # Thrust magnitued
    T = 10e-3 # [kN]

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

function ekfWithLunarPertEOM!(dy, y, ekf::EKF, t)
    # Zero F 
    ekf.F .= 0.0

    # Earth gravitational parameter
    μ = 5.9742e24*6.673e-20

    # Luna gravitational parameter
    μl = 7.3483e22*6.673e-20

    # Thrust magnitued
    T = 10e-3 # [kN]

    # Get control
    (u, rLuna) = GetControlAndLunaPos(ekf.scSim, t)

    # Luna gravitational perturbation
    rls         = @SVector [rLuna[1] - y[1], rLuna[2] - y[2], rLuna[3] - y[3]]
    μldRls3     = μl / (norm(rls)^3)
    μldRluna3   = μl / (norm(rLuna)^3)
    al          = @SVector [rls[1]*μldRls3 - rLuna[1]*μldRluna3,
                            rls[2]*μldRls3 - rLuna[2]*μldRluna3,
                            rls[3]*μldRls3 - rLuna[3]*μldRluna3]
    # States
    @views r        = norm(y[1:3])
    @views dy[1:3] .= y[4:6]
    @views dy[4:6] .= -(μ/r^3).*y[1:3] .+ al .+ (T*u[1] / y[7]).*u[2:4]
    dy[7]           = -u[1]*T / ekf.scSim.ps.sp.c

    # Dynamics Jacobian
    # dr/dv
    ekf.F[1,4] = 1.0
    ekf.F[2,5] = 1.0
    ekf.F[3,6] = 1.0

    # dv/dr keplarian
    @views mul!(ekf.F[4:6, 1:3], y[1:3], transpose(y[1:3]))
    @views ekf.F[4:6, 1:3] .*= (3.0 / r^2)
    for i in 1:3; ekf.F[3 + i, i] -= 1.0; end
    @views ekf.F[4:6, 1:3] .*= (μ/r^3)

    # dv/dr luna perturbation (stealing memory from ekf.H)
    @views mul!(ekf.H[1:3, 1:3], rls, transpose(rls))
    @views ekf.H[1:3, 1:3] .*= (3.0 / norm(rls)^2)
    for i in 1:3; ekf.H[i,i] -= 1.0; end
    @views ekf.H[1:3, 1:3] .*= (μl / norm(rls)^3)

    # Combine dv/dr partials
    @views ekf.F[4:6, 1:3] .+= ekf.H[1:3, 1:3]

    # dv/dm
    @views ekf.F[4:6, 7] .= -(u[1]*T / y[7]^2) .* u[2:4]

    # Covariance Dynamics
    @views P   = reshape(y[8:end], (7, 7))
    @views dP  = reshape(dy[8:end], (7, 7))
    mul!(dP, ekf.F, P)
    mul!(dP, P, transpose(ekf.F), 1.0, 1.0)
    dP .+= ekf.Q
end

function plot(ekf::EKF, xtrue, n)
    ts   = ekf.txp[1:n]
    xhat = ekf.xhats[1:n, :]
    xt   = xtrue[1:n, :]
    es   = ekf.es[1:n, :]
    σ311 = 3*sqrt.(ekf.Ps[1:n,1])
    σ322 = 3*sqrt.(ekf.Ps[1:n,2])
    σ333 = 3*sqrt.(ekf.Ps[1:n,3])
    σ344 = 3*sqrt.(ekf.Ps[1:n,4])
    σ355 = 3*sqrt.(ekf.Ps[1:n,5])
    σ366 = 3*sqrt.(ekf.Ps[1:n,6])
    σ377 = 3*sqrt.(ekf.Ps[1:n,7])

    mat"""
    ts      = $ts;
    ts      = (ts - $(ekf.scSim.t0))/86400;
    xhat    = $xhat;
    xt      = $xt;
    es      = $es;
    s311    = $σ311;
    s322    = $σ322;
    s333    = $σ333;
    s344    = $σ344;
    s355    = $σ355;
    s366    = $σ366;
    s377    = $σ377;

    figure()
    subplot(3,1,1)
    plot(ts, es(:,1), 'k')
    hold on
    plot(ts, s311, 'r')
    plot(ts, -s311, 'r')
    legend('Estimation Error', 'EKF \$3\\sigma\$', 'Interpreter', 'latex')
    grid on

    subplot(3,1,2)
    plot(ts, es(:,1), 'k')
    hold on
    plot(ts, s322, 'r')
    plot(ts, -s322, 'r')
    grid on

    subplot(3,1,3)
    plot(ts, es(:,3), 'k')
    hold on
    plot(ts, s333, 'r')
    plot(ts, -s333, 'r')
    xlabel("Time, days")
    grid on

    figure()
    subplot(3,1,1)
    plot(ts, es(:,4), 'k')
    hold on
    plot(ts, s344, 'r')
    plot(ts, -s344, 'r')

    subplot(3,1,2)
    plot(ts, es(:,5), 'k')
    hold on
    plot(ts, s355, 'r')
    plot(ts, -s355, 'r')

    subplot(3,1,3)
    plot(ts, es(:,6), 'k')
    hold on
    plot(ts, s366, 'r')
    plot(ts, -s366, 'r')

    figure()
    plot(ts, es(:,7), 'k')
    hold on
    plot(ts, s377, 'r')
    plot(ts, -s377, 'r')

    figure()
    plot3(xhat(:,1),xhat(:,2), xhat(:,3), 'r')
    hold on
    plot3(xt(:,1),xt(:,2),xt(:,3), 'b')
    grid on
    axis equal
    """
end

function plotMC(ekf::EKF, em, ev, n)
    ts      = ekf.txp[1:n]
    xhat    = ekf.xhats[1:n, :]
    es      = ekf.es[1:n, :]
    σ311    = 3*sqrt.(ekf.Ps[1:n,1])
    σ311mc  = 3*sqrt.(ev[:,1])
    σ322    = 3*sqrt.(ekf.Ps[1:n,2])
    σ322mc  = 3*sqrt.(ev[:,2])
    σ333    = 3*sqrt.(ekf.Ps[1:n,3])
    σ333mc  = 3*sqrt.(ev[:,3])
    σ344    = 3*sqrt.(ekf.Ps[1:n,4])
    σ344mc  = 3*sqrt.(ev[:,4])
    σ355    = 3*sqrt.(ekf.Ps[1:n,5])
    σ355mc  = 3*sqrt.(ev[:,5])
    σ366    = 3*sqrt.(ekf.Ps[1:n,6])
    σ366mc  = 3*sqrt.(ev[:,6])
    σ377    = 3*sqrt.(ekf.Ps[1:n,7])
    σ377mc  = 3*sqrt.(ev[:,7])

    mat"""
    ts      = $ts;
    ts      = (ts - $(ekf.scSim.t0))/86400;
    xhat    = $xhat;
    es      = $es;
    s311    = $σ311;
    s311mc  = $σ311mc;
    s322    = $σ322;
    s322mc  = $σ322mc;
    s333    = $σ333;
    s333mc  = $σ333mc;
    s344    = $σ344;
    s344mc  = $σ344mc;
    s355    = $σ355;
    s355mc  = $σ355mc;
    s366    = $σ366;
    s366mc  = $σ366mc;
    s377    = $σ377;
    s377mc  = $σ377mc;
    em      = $em;

    figure()
    subplot(3,1,1)
    plot(ts, es(:,1).*1000, 'k')
    hold on
    plot(ts, s311.*1000, 'r')
    plot(ts, -s311.*1000, 'r')
    plot(ts, s311mc.*1000, '--r')
    plot(ts, -s311mc.*1000, '--r')
    legend('Estimation Error', 'EKF \$3\\sigma\$', 'MC \$3\\sigma\$','Interpreter', 'latex')
    ylabel('\$e_{r_x}\$, m', 'Interpreter', 'latex')

%    ylim([-1000, 1000])
    grid on

    subplot(3,1,2)
    plot(ts, es(:,1).*1000, 'k')
    hold on
    plot(ts, s322.*1000, 'r')
    plot(ts, -s322.*1000, 'r')
    plot(ts, s322mc.*1000, '--r')
    plot(ts, -s322mc.*1000, '--r')
    ylabel('\$e_{r_y}\$, m', 'Interpreter', 'latex')

%    ylim([-1000, 1000])
    grid on

    subplot(3,1,3)
    plot(ts, es(:,3).*1000, 'k')
    hold on
    plot(ts, s333.*1000, 'r')
    plot(ts, -s333.*1000, 'r')
    plot(ts, s333mc.*1000, '--r')
    plot(ts, -s333mc.*1000, '--r')
    xlabel('Time, days')
    ylabel('\$e_{r_z}\$, m', 'Interpreter', 'latex')

    grid on

    figure()
    subplot(3,1,1)
    plot(ts, es(:,4).*1000, 'k')
    hold on
    plot(ts, s344.*1000, 'r')
    plot(ts, -s344.*1000, 'r')
    plot(ts, s344mc.*1000, '--r')
    plot(ts, -s344mc.*1000, '--r')
    legend('Estimation Error', 'UKF \$3\\sigma\$', 'MC \$3\\sigma\$','Interpreter', 'latex')
    ylabel('\$e_{v_x}\$, m/s', 'Interpreter', 'latex')

    ylim([-1, 1])
    grid on

    subplot(3,1,2)
    plot(ts, es(:,5).*1000, 'k')
    hold on
    plot(ts, s355.*1000, 'r')
    plot(ts, -s355.*1000, 'r')
    plot(ts, s355mc.*1000, '--r')
    plot(ts, -s355mc.*1000, '--r')
    ylabel('\$e_{v_y}\$, m/s', 'Interpreter', 'latex')

    ylim([-1.5, 1.5])
    grid on

    subplot(3,1,3)
    plot(ts, es(:,6).*1000, 'k')
    hold on
    plot(ts, s366.*1000, 'r')
    plot(ts, -s366.*1000, 'r')
    plot(ts, s366mc.*1000, '--r')
    plot(ts, -s366mc.*1000, '--r')
    xlabel('Time, days')
    ylabel('\$e_{v_z}\$, m/s', 'Interpreter', 'latex')

    ylim([-1, 1])
    grid on

    figure()
    plot(ts, es(:,7), 'k')
    hold on
    plot(ts, s377, 'r')
    plot(ts, -s377, 'r')
    plot(ts, s377mc, '--r')
    plot(ts, -s377mc, '--r')
    xlabel('Time, days')
    ylabel('\$e_m\$, kg', 'Interpreter', 'latex')
    legend('Estimation Error', 'UKF \$3\\sigma\$', 'MC \$3\\sigma\$','Interpreter', 'latex')

    grid on

    figure()
    histogram(em(:,1).*1000, 20)
    hold on
    xline(mean(em(:,1))*1000, '--k', 'LineWidth', 3) 
    xlabel('\$e_{r_x}\$, m', 'Interpreter', 'latex')
    ylabel('# Occurances')
    grid on

    figure()
    histogram(em(:,2).*1000, 20)
    hold on
    xline(mean(em(:,2))*1000, '--k', 'LineWidth', 3)
    xlabel('\$e_{r_y}\$, m', 'Interpreter', 'latex')
    ylabel('# Occurances')
    grid on

    figure()
    histogram(em(:,3).*1000, 20)
    hold on
    xline(mean(em(:,3))*1000, '--k', 'LineWidth', 3)
    xlabel('\$e_{r_z}\$, m', 'Interpreter', 'latex')
    ylabel('# Occurances')
    grid on

    figure()
    histogram(em(:,4).*1000, 20)
    hold on
    xline(mean(em(:,4))*1000, '--k', 'LineWidth', 3) 
    xlabel('\$e_{v_x}\$, m/s', 'Interpreter', 'latex')
    ylabel('# Occurances')
    grid on

    figure()
    histogram(em(:,5).*1000, 20)
    hold on
    xline(mean(em(:,5))*1000, '--k', 'LineWidth', 3)
    xlabel('\$e_{v_y}\$, m/s', 'Interpreter', 'latex')
    ylabel('# Occurances')
    grid on

    figure()
    histogram(em(:,6).*1000, 30)
    hold on
    xline(mean(em(:,6))*1000, '--k', 'LineWidth', 3)
    xlabel('\$e_{v_z}\$, m/s', 'Interpreter', 'latex')
    ylabel('# Occurances')
    grid on

    figure()
    histogram(em(:,7), 20)
    hold on
    xline(mean(em(:,7)), '--k', 'LineWidth', 3)
    xlabel('\$e_m\$, m/s', 'Interpreter', 'latex')
    ylabel('# Occurances')
    grid on

    """
end

function saveMCmat(ekf::EKF, em, ev, n)
    ts      = ekf.txp[1:n]
    xhat    = ekf.xhats[1:n, :]
    es      = ekf.es[1:n, :]
    σ311    = 3*sqrt.(ekf.Ps[1:n,1])
    σ311mc  = 3*sqrt.(ev[:,1])
    σ322    = 3*sqrt.(ekf.Ps[1:n,2])
    σ322mc  = 3*sqrt.(ev[:,2])
    σ333    = 3*sqrt.(ekf.Ps[1:n,3])
    σ333mc  = 3*sqrt.(ev[:,3])
    σ344    = 3*sqrt.(ekf.Ps[1:n,4])
    σ344mc  = 3*sqrt.(ev[:,4])
    σ355    = 3*sqrt.(ekf.Ps[1:n,5])
    σ355mc  = 3*sqrt.(ev[:,5])
    σ366    = 3*sqrt.(ekf.Ps[1:n,6])
    σ366mc  = 3*sqrt.(ev[:,6])
    σ377    = 3*sqrt.(ekf.Ps[1:n,7])
    σ377mc  = 3*sqrt.(ev[:,7])

    fname   = datadir("MATLABdata","ekf.mat")

    mat"""
    ts      = $ts;
    ts      = (ts - $(ekf.scSim.t0))/86400;
    xhat    = $xhat;
    es      = $es;
    s311    = $σ311;
    s311mc  = $σ311mc;
    s322    = $σ322;
    s322mc  = $σ322mc;
    s333    = $σ333;
    s333mc  = $σ333mc;
    s344    = $σ344;
    s344mc  = $σ344mc;
    s355    = $σ355;
    s355mc  = $σ355mc;
    s366    = $σ366;
    s366mc  = $σ366mc;
    s377    = $σ377;
    s377mc  = $σ377mc;
    em      = $em;

    save($fname)
    """

    return nothing
end