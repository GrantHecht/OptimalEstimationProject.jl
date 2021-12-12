
mutable struct UKF 
    # Measurement times
    ts::Vector{Float64}
    imuΔt::Float64
    gpsΔt::Float64

    # Measurement residuals 
    rs::Matrix{Float64}

    # State estimates and covariance matrix diags
    txp::Vector{Float64}    # Epochs associated with entries in xhats and Ps
    ixp::Int                # Index associated with latest entry
    xhats::Matrix{Float64}  # Matrix of state estimates
    es::Matrix{Float64}     # Matrix of state estimation errors
    Ps::Matrix{Float64}     # Matrix of covariance matrix diagonals

    # Parameters 
    α::Float64
    β::Float64
    κ::Float64
    λ::Float64
    γ::Float64

    # Preallocate matricies
    χ::Matrix{Float64}
    χd::Matrix{Float64}
    P⁺::Matrix{Float64}
    P⁻::Matrix{Float64}
    Py::Matrix{Float64}
    Pyt::Matrix{Float64}
    Pxy::Matrix{Float64}
    Q::Diagonal{Float64, Vector{Float64}}
    expmeas::Matrix{Float64}
    expmeasd::Matrix{Float64}
    yhat::Vector{Float64}

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

    # Resample sigma points for update
    resamp::Bool
end

function UKF(xhat0, P0, Q, σ2GPS, σ2Acc, ts, gpsΔt, gpsSim, imuSim, scSim; 
    α = 1.0, β = 2.0, κ = 1.0, lunaPerts = false, steps2save = 5, resample = false)

    # Instantiate preallocated matricies
    ixp     = 1
    rs      = zeros(length(ts), 35)
    txp     = zeros((steps2save + 1)*length(ts) + 1); txp[1] = scSim.t0
    xhats   = zeros((steps2save + 1)*length(ts) + 1, 7); xhats[1, :] .= xhat0
    es      = zeros((steps2save + 1)*length(ts) + 1, 7);
    Ps      = zeros((steps2save + 1)*length(ts) + 1, 7); Ps[1, :] .= diag(P0)
    χ       = zeros(7,15)
    χd      = zeros(7,15)
    P⁻      = zeros(7,7)
    P⁺      = deepcopy(P0) 
    Py      = zeros(35, 35)
    Pyt     = zeros(35, 35)
    Pxy     = zeros(7, 35)
    expmeas = zeros(35, 15)
    expmeasd= zeros(35, 15)
    yhat    = zeros(35)
   
    # Compute IMU Δt from ts vector 
    imuΔt   = ts[2] - ts[1]

    # Compute initial estimates error 
    (ytrue, u) = GetStateAndControl(scSim, scSim.t0)
    @views es[1,:] .= xhat0 .- ytrue[1:7]

    # Compute gamma 
    λ       = α^2*(7.0 + κ) - 7.0
    γ       = sqrt(7.0 + λ)

    # Instantiate UKF
    UKF(ts,imuΔt,gpsΔt,rs,txp,ixp,xhats,es,Ps,α,β,κ,λ,γ,χ,χd,P⁺,P⁻,Py,Pyt,Pxy,Q,expmeas,expmeasd,yhat,
        gpsSim,imuSim,scSim,σ2GPS,σ2Acc,steps2save,0,true,false,lunaPerts,resample)
end

function propagate!(ukf::UKF)
    if ukf.updated == false 
        throw(ErrorException("Cannot propagate filter before updating."))
    else
        ukf.k += 1
        ukf.updated = false
    end

    # Get time 
    if ukf.k == 1
        t0  = ukf.scSim.t0 
        tf  = ukf.ts[1]
    else
        t0  = ukf.ts[ukf.k - 1]
        tf  = ukf.ts[ukf.k]
    end

    # Generate vector of times to save at 
    saveat  = range(t0; length = ukf.stp, stop = tf)

    # Compute sigma points
    C       = cholesky(ukf.P⁺)
    for i in 1:7
        @views ukf.χ[:,2*i - 1] .= ukf.xhats[ukf.ixp, :] .+ ukf.γ.*C.L[:,i]
        @views ukf.χ[:,2*i]     .= ukf.xhats[ukf.ixp, :] .- ukf.γ.*C.L[:,i]
    end
    ukf.χ[:,15] .= ukf.xhats[ukf.ixp, :]

    # Propagate sigma points 
    if ukf.lunaPerts == false
        @views sols = [solve(ODEProblem(ukfNoLunarPertEOM, SVector{7}(ukf.χ[:,i]), (t0, tf), ukf), 
            Vern7(); reltol=1e-10, abstol=1e-10, saveat=saveat) for i in 1:15]
    else
        @views sols = [solve(ODEProblem(ukfWithLunarPertEOM, SVector{7}(ukf.χ[:,i]), (t0, tf), ukf), 
            Vern7(); reltol=1e-10, abstol=1e-10, saveat=saveat) for i in 1:15]
    end

    # Update storage matricies
    start   = ukf.ixp + 1
    stop    = ukf.ixp + ukf.stp 
    step    = 0
    for i in start:stop
        step += 1

        # Time 
        ukf.txp[i] = sols[1].t[step]

        # Fill matrix of sigma point states at step "step"
        for j in 1:15
            for k in 1:7
                ukf.χ[k,j] = sols[j].u[step][k]
            end
        end

        # Compute weights
        w0m     = ukf.λ / (7.0 + ukf.λ)
        w0c     = w0m + (1 - ukf.α^2 + ukf.β)
        wi      = 1.0 / (2.0*(7 + ukf.λ))

        # State estimate
        ukf.xhats[i,:] .*= 0.0
        for j in 1:14
            @views ukf.xhats[i,:] .+= wi.*ukf.χ[:,j]
        end
        @views ukf.xhats[i,:] .+= w0m.*ukf.χ[:,15]

        # Error 
        (ytrue, u) = GetStateAndControl(ukf.scSim, ukf.txp[i])
        @views ukf.es[i,:] .= ukf.xhats[i,:] .- ytrue[1:7]

        # Covariance
        ukf.P⁻ .*= 0.0
        for j in 1:14
            @views ukf.χd[:,j] .= ukf.χ[:,j] .- ukf.xhats[i,:]
            @views mul!(ukf.P⁺, ukf.χd[:,j], transpose(ukf.χd[:,j])) # Using P⁺ for storage of data
            ukf.P⁻ .+= wi.*ukf.P⁺
        end
        @views ukf.χd[:,15] .= ukf.χ[:,15] .- ukf.xhats[i,:]
        @views mul!(ukf.P⁺, ukf.χd[:,15], transpose(ukf.χd[:,15]))
        ukf.P⁻ .+= w0c.*ukf.P⁺

        @views ukf.Ps[i,:] .= diag(ukf.P⁻)
        ukf.ixp = i
    end
end

function updateGPS!(ukf::UKF)
    if ukf.updated == true
        throw(ErrorException("Cannot update again before propagating."))
    else
        ukf.updated = true 
    end

    # Get current time 
    t   = ukf.ts[ukf.k]

    # Get true state 
    (ytrue, u) = GetStateAndControl(ukf.scSim, t)

    # Generate GPS measurements 
    @views measGPS = gpsMeasurement(ukf.gpsSim, ytrue[1:3], t)

    # Get satelite ids 
    @views sats = measGPS[2:2:end]
    numSats = length(sats)

    # Resample sigma points if desired 
    if ukf.resamp
        # Compute sigma points
        C       = cholesky(ukf.P⁻)
        for i in 1:7
            @views ukf.χ[:,2*i - 1] .= ukf.xhats[ukf.ixp, :] .+ ukf.γ.*C.L[:,i]
            @views ukf.χ[:,2*i]     .= ukf.xhats[ukf.ixp, :] .- ukf.γ.*C.L[:,i]
        end
        ukf.χ[:,15] .= ukf.xhats[ukf.ixp, :]

        # Compute sigma point diffs
        for i in 1:15
            @views ukf.χd[:,i] .= ukf.χ[:,i] .- ukf.xhats[ukf.ixp, :]
        end
    end

    # Compute expected measurements
    ukf.expmeas .*= 0.0
    for i in 1:15
        # Get expected gps measurement
        @views (expMeasGPS, ps) = gpsMeasurement(ukf.gpsSim, ukf.χ[1:3,i], sats, t; type = :expected)
        for j in 1:numSats 
            ukf.expmeas[j,i] = expMeasGPS[3 + 2*(j - 1)]
        end
    end

    # Compute weights
    w0m     = ukf.λ / (7.0 + ukf.λ)
    w0c     = w0m + (1 - ukf.α^2 + ukf.β)
    wi      = 1.0 / (2.0*(7 + ukf.λ))

    # Compute mean observation
    ukf.yhat .*= 0.0
    for i in 1:14
        @views ukf.yhat[1:numSats] .+= wi.*ukf.expmeas[1:numSats,i]
    end
    @views ukf.yhat[1:numSats] .+= w0m.*ukf.expmeas[1:numSats,15]

    # Compute residuals with measurement rejection
    accMeas     = Vector{Bool}(undef, 35); accMeas .= false
    rejectTol   = 1e9
    numRejected = 0
    for i in 1:numSats 
        r = measGPS[3 + 2*(i - 1)] - ukf.yhat[i]
        if r < rejectTol
            accMeas[i] = true
            ukf.rs[ukf.k, i - numRejected] = r
        else
            numRejected += 1
        end
    end

    # Compute output covariance
    ukf.Py .= 0.0
    for j in 1:14
        @views ukf.expmeasd[1:numSats,j] .= ukf.expmeas[1:numSats,j] .- ukf.yhat[1:numSats]
        @views mul!(ukf.Py[1:numSats-numRejected,1:numSats-numRejected], 
            ukf.expmeasd[accMeas,j], transpose(ukf.expmeasd[accMeas,j]), wi, 1.0)
        #@views ukf.Py[1:numSats-numRejected,1:numSats-numRejected] .+= 
        #    wi.*ukf.Pyt[1:numSats-numRejected,1:numSats-numRejected]
    end
    @views ukf.expmeasd[1:numSats,15] .= ukf.expmeas[1:numSats,15] .- ukf.yhat[1:numSats]
    @views mul!(ukf.Py[1:numSats-numRejected,1:numSats-numRejected], 
        ukf.expmeasd[accMeas,15], transpose(ukf.expmeasd[accMeas,15]), w0c, 1.0)
    #@views ukf.Py[1:numSats-numRejected,1:numSats-numRejected] .+= 
    #    w0c.*ukf.Pyt[1:numSats-numRejected,1:numSats-numRejected]
    for i in 1:numSats-numRejected; ukf.Py[i,i] += ukf.σ2GPS; end

    # Compute cross covariance
    ukf.Pxy .= 0.0
    for j in 1:14
        @views mul!(ukf.Pxy[1:7,1:numSats-numRejected], 
            ukf.χd[:,j], transpose(ukf.expmeasd[accMeas,j]), wi, 1.0)
        #@views ukf.Pxy[1:7,1:numSats-numRejected] .+= 
        #    wi.*ukf.Pyt[1:7,1:numSats-numRejected]
    end
    @views mul!(ukf.Pxy[1:7,1:numSats-numRejected],
        ukf.χd[:,15], transpose(ukf.expmeasd[accMeas,15]), w0c, 1.0)
    #@views ukf.Pxy[1:7,1:numSats-numRejected] .+= 
    #    w0c.*ukf.Pyt[1:7,1:numSats-numRejected]

    # Compute gain 
    @views Kfact = factorize(ukf.Py[1:numSats-numRejected,1:numSats-numRejected])
    @views rdiv!(ukf.Pxy[1:7,1:numSats-numRejected], Kfact)

    # Compute state update
    ukf.ixp += 1
    @views mul!(ukf.xhats[ukf.ixp, :], ukf.Pxy[1:7,1:numSats-numRejected], 
        ukf.rs[1:numSats-numRejected])
    ukf.xhats[ukf.ixp,:] .+= ukf.xhats[ukf.ixp - 1, :]

    # Compute covariance update
    ukf.P⁺ .= ukf.P⁻
    @views mul!(ukf.Pyt[1:7,1:numSats-numRejected], ukf.Pxy[1:7,1:numSats-numRejected], 
        ukf.Py[1:numSats-numRejected,1:numSats-numRejected])
    @views mul!(ukf.P⁺, ukf.Pyt[1:7,1:numSats-numRejected], 
        transpose(ukf.Pxy[1:7,1:numSats-numRejected]), -1.0, 1.0)

    # Set time in storage vector
    ukf.txp[ukf.ixp] = t

    # Compute error
    (ytrue,u) = GetStateAndControl(ukf.scSim, t)
    @views ukf.es[ukf.ixp,:] .= ukf.xhats[ukf.ixp,:] .- ytrue[1:7]

    # Set covariance diag 
    ukf.Ps[ukf.ixp,:] .= diag(ukf.P⁺)

    # Force covariance symmetry
    for r in 1:6
        for c in r+1:7
            ukf.P⁺[r,c] = ukf.P⁺[c,r]
        end
    end
end

function runFilter!(ukf::UKF)
    for i in 1:length(ukf.ts)
        propagate!(ukf)
        updateGPS!(ukf)
    end
end

function ukfNoLunarPertEOM(y, ukf::UKF, t)
    # Earth gravitational parameter
    μ = 3.986004418e5 # [km^3/s^2]

    # Thrust magnitued
    T = 10e-3 # [kN]

    # Get control
    u = GetControl(ukf.scSim, t)

    # States
    @views r    = norm(y[1:3])
    nμr3        = -μ/r^3
    at          = T*u[1] / y[7] 
    dy          = @SVector [y[4], y[5], y[6], 
                            nμr3*y[1] + at*u[2],
                            nμr3*y[2] + at*u[3],
                            nμr3*y[3] + at*u[4],
                            -u[1]*T / ukf.scSim.ps.sp.c]
                            
    return dy
end

function ukfWithLunarPertEOM(y, ukf::UKF, t)
    # Earth gravitational parameter
    μ = 5.9742e24*6.673e-20

    # Luna gravitational parameter
    μl = 7.3483e22*6.673e-20

    # Thrust magnitued
    T = 10e-3 # [kN]

    # Get control
    (u, rLuna) = GetControlAndLunaPos(ukf.scSim, t)

    # Luna gravitational perturbation
    rls         = @SVector [rLuna[1] - y[1], rLuna[2] - y[2], rLuna[3] - y[3]]
    μldRls3     = μl / (norm(rls)^3)
    μldRluna3   = μl / (norm(rLuna)^3)
    al          = @SVector [rls[1]*μldRls3 - rLuna[1]*μldRluna3,
                            rls[2]*μldRls3 - rLuna[2]*μldRluna3,
                            rls[3]*μldRls3 - rLuna[3]*μldRluna3]
    # States
    @views r    = norm(y[1:3])
    nμr3        = -μ/r^3
    at          = T*u[1] / y[7] 
    dy          = @SVector [y[4], y[5], y[6], 
                            nμr3*y[1] + al[1] + at*u[2],
                            nμr3*y[2] + al[2] + at*u[3],
                            nμr3*y[3] + al[3] + at*u[4],
                            -u[1]*T / ukf.scSim.ps.sp.c]
    
    return dy
end