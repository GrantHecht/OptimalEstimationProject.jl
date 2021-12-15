
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

    # Augmented state vector size 
    L::Int

    # Preallocate matricies
    χ::Matrix{Float64}
    χd::Matrix{Float64}
    Pa::SparseMatrixCSC{Float64, Int64}
    C::SuiteSparse.CHOLMOD.Factor{Float64}
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

    # Step Index
    k::Int64

    # Update flag
    updated::Bool

    # Finished flag
    finished::Bool

    # Luna perturbation flag
    lunaPerts::Bool

end

function UKF(xhat0, P0, Q, σ2GPS, σ2Acc, ts, gpsΔt, gpsSim, imuSim, scSim; 
    α = 1.0, β = 2.0, κ = 1.0, lunaPerts = false)

    # Instantiate preallocated matricies
    ixp     = 1
    rs      = zeros(length(ts), 35)
    txp     = zeros(2*length(ts) + 1); txp[1] = scSim.t0
    xhats   = zeros(2*length(ts) + 1, 7); xhats[1, :] .= xhat0
    es      = zeros(2*length(ts) + 1, 7);
    Ps      = zeros(2*length(ts) + 1, 7); Ps[1, :] .= diag(P0)
    χ       = zeros(14,29)
    χd      = zeros(7,29)
    P⁻      = zeros(7,7)
    P⁺      = deepcopy(P0) 
    Py      = zeros(35, 35)
    Pyt     = zeros(35, 35)
    Pxy     = zeros(7, 35)
    expmeas = zeros(35, 29)
    expmeasd= zeros(35, 29)
    yhat    = zeros(35)
    L       = 14

    # Construct sparse augmented covariance
    Pa      = sparse(zeros(14,14))
    for c in 1:7
        for r in 1:7
            Pa[r,c]     = P⁺[r,c]
            if c == r
                Pa[7 + r,7 + c]     = Q[r,c]
            end
        end
    end
    C       = cholesky(Pa; perm = 1:14)
   
    # Compute IMU Δt from ts vector 
    imuΔt   = ts[2] - ts[1]

    # Compute initial estimates error 
    (ytrue, u) = GetStateAndControl(scSim, scSim.t0)
    @views es[1,:] .= xhat0 .- ytrue[1:7]

    # Compute gamma 
    λ       = α^2*(L + κ) - L
    γ       = sqrt(L + λ)

    # Instantiate UKF
    UKF(ts,imuΔt,gpsΔt,rs,txp,ixp,xhats,es,Ps,α,β,κ,λ,γ,L,χ,χd,Pa,C,P⁺,P⁻,Py,Pyt,Pxy,Q,expmeas,expmeasd,yhat,
        gpsSim,imuSim,scSim,σ2GPS,σ2Acc,0,true,false,lunaPerts)
end

function propagate!(ukf::UKF)
    if ukf.updated == false 
        throw(ErrorException("Cannot propagate filter before updating."))
    else
        ukf.k += 1
        ukf.updated = false
    end

    @inbounds begin
    # Get time 
    if ukf.k == 1
        t0  = ukf.scSim.t0 
        tf  = ukf.ts[ukf.k]
    else
        t0  = ukf.ts[ukf.k - 1]
        tf  = ukf.ts[ukf.k]
    end

    # Generate vector of times to save at 
    saveat = range(t0; length = 2, stop = tf)

    # Compute sigma points
    cholesky!(ukf.C, ukf.Pa)
    L       = sparse(ukf.C.L)
    ukf.χ  .= 0.0
    for i in 1:ukf.L
        @views ukf.χ[1:7,2*i - 1]   .= ukf.xhats[ukf.ixp, :] .+ ukf.γ.*L[1:7,i]
        @views ukf.χ[8:14,2*i - 1]  .= ukf.γ.*L[8:14,i]
        @views ukf.χ[1:7,2*i]       .= ukf.xhats[ukf.ixp, :] .- ukf.γ.*L[1:7,i]
        @views ukf.χ[8:14,2*i]      .= -ukf.γ.*L[8:14,i]
    end
    ukf.χ[1:7,2*ukf.L + 1] .= ukf.xhats[ukf.ixp, :]

    # Propagate sigma points 
    if ukf.lunaPerts == false
        #@views sols = [rkProp(ukfNoLunarPertEOM, ukf.χ[:,i], saveat, ukf) for i in 1:2*ukf.L + 1]
        @views sols = [solve(
            ODEProblem{false}((u,p,t)->ukfNoLunarPertEOM(u,p,t,ukf.χ[8:14,i]), 
            SVector{7}(ukf.χ[1:7,i]), (t0, tf), ukf), Vern7();
            reltol=1e-10, abstol=1e-10, saveat=saveat) for i in 1:2*ukf.L + 1]
    else
        #@views sols = [rkProp(ukfWithLunarPertEOM, ukf.χ[:,i], saveat, ukf) for i in 1:2*ukf.L + 1]
        @views sols = [solve(
            ODEProblem{false}((u,p,t)->ukfWithLunarPertEOM(u,p,t,ukf.χ[8:14,i]), 
            SVector{7}(ukf.χ[1:7,i]), (t0, tf), ukf), Vern7();
            reltol=1e-10, abstol=1e-10, saveat=saveat) for i in 1:2*ukf.L + 1]
    end

    # Update storage matricies
    ukf.ixp += 1

    # Time 
    ukf.txp[ukf.ixp] = tf

    # Fill matrix of sigma point states 
    for j in 1:2*ukf.L + 1
        for k in 1:7
            ukf.χ[k,j] = sols[j].u[end][k]
        end
    end

    # Compute weights
    w0m     = ukf.λ / (ukf.L + ukf.λ)
    w0c     = w0m + (1 - ukf.α^2 + ukf.β)
    wi      = 1.0 / (2.0*(ukf.L + ukf.λ))

    # State estimate
    ukf.xhats[ukf.ixp,:] .*= 0.0
    for j in 1:2*ukf.L
        @views ukf.xhats[ukf.ixp,:] .+= wi.*ukf.χ[1:7,j]
    end
    @views ukf.xhats[ukf.ixp,:] .+= w0m.*ukf.χ[1:7,2*ukf.L + 1]

    # Error 
    (ytrue, u) = GetStateAndControl(ukf.scSim, tf)
    @views ukf.es[ukf.ixp,:] .= ukf.xhats[ukf.ixp,:] .- ytrue[1:7]

    # Covariance
    ukf.P⁻ .*= 0.0
    for j in 1:2*ukf.L
        @views ukf.χd[:,j] .= ukf.χ[1:7,j] .- ukf.xhats[ukf.ixp,:]
        @views mul!(ukf.P⁻, ukf.χd[:,j], transpose(ukf.χd[:,j]), wi, 1.0) # Using P⁺ for storage of data
    end
    @views ukf.χd[:,2*ukf.L + 1] .= ukf.χ[1:7,2*ukf.L + 1] .- ukf.xhats[ukf.ixp,:]
    @views mul!(ukf.P⁻, ukf.χd[:,2*ukf.L + 1], transpose(ukf.χd[:,2*ukf.L + 1]), w0c, 1.0)

    @views ukf.Ps[ukf.ixp,:] .= diag(ukf.P⁻)
    end
    nothing
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

    # Compute expected measurements
    ukf.expmeas .*= 0.0
    for i in 1:2*ukf.L + 1
        # Get expected gps measurement
        @views (expMeasGPS, ps) = gpsMeasurement(ukf.gpsSim, ukf.χ[1:3,i], sats, t; type = :expected)
        for j in 1:numSats 
            ukf.expmeas[j,i] = expMeasGPS[3 + 2*(j - 1)]
        end
    end

    # Compute weights
    w0m     = ukf.λ / (ukf.L + ukf.λ)
    w0c     = w0m + (1 - ukf.α^2 + ukf.β)
    wi      = 1.0 / (2.0*(ukf.L + ukf.λ))

    # Compute mean observation
    ukf.yhat .= 0.0
    for i in 1:2*ukf.L
        @views ukf.yhat[1:numSats] .+= wi.*ukf.expmeas[1:numSats,i]
    end
    @views ukf.yhat[1:numSats] .+= w0m.*ukf.expmeas[1:numSats,2*ukf.L + 1]

    # Compute residuals with measurement rejection
    accMeas     = Vector{Bool}(undef, 35); accMeas .= false
    rejectTol   = 10.0
    numRejected = 0
    for i in 1:numSats 
        r = measGPS[3 + 2*(i - 1)] - ukf.yhat[i]
        if abs(r) < rejectTol
            accMeas[i] = true
            ukf.rs[ukf.k, i - numRejected] = r
        else
            numRejected += 1
        end
    end

    # Compute output covariance
    ukf.Py      .= 0.0
    numMeas  = numSats - numRejected
    for j in 1:2*ukf.L
        @views ukf.expmeasd[1:numMeas,j] .= ukf.expmeas[accMeas,j] .- ukf.yhat[accMeas]
        @views mul!(ukf.Py[1:numMeas,1:numMeas], ukf.expmeasd[1:numMeas,j], 
            transpose(ukf.expmeasd[1:numMeas,j]), wi, 1.0)
    end
    @views ukf.expmeasd[1:numMeas,2*ukf.L + 1] .= ukf.expmeas[accMeas,2*ukf.L + 1] .- ukf.yhat[accMeas]
    @views mul!(ukf.Py[1:numMeas,1:numMeas], ukf.expmeasd[1:numMeas,2*ukf.L + 1], 
        transpose(ukf.expmeasd[1:numMeas,2*ukf.L + 1]), w0c, 1.0)

    # Update output covariance to innovations covariance (assuming linearly appearing measurement noise)
    for i in 1:numMeas; ukf.Py[i,i] += ukf.σ2GPS; end

    # Compute cross correlation matrix
    ukf.Pxy .= 0.0
    for j in 1:2*ukf.L 
        @views mul!(ukf.Pxy[1:7,1:numMeas], ukf.χd[:,j], transpose(ukf.expmeasd[1:numMeas,j]), wi, 1.0)
    end
    @views mul!(ukf.Pxy[1:7,1:numMeas], ukf.χd[:,2*ukf.L + 1], transpose(ukf.expmeasd[1:numMeas,2*ukf.L + 1]), w0c, 1.0)

    # Compute gain 
    @views Pyfact = factorize(ukf.Py[1:numMeas,1:numMeas])
    @views rdiv!(ukf.Pxy[1:7,1:numMeas], Pyfact)

    # Update storage matricies
    # Time
    ukf.ixp += 1
    ukf.txp[ukf.ixp] = t

    # Compute state update and save
    @views mul!(ukf.xhats[ukf.ixp, :], ukf.Pxy[1:7,1:numMeas], ukf.rs[ukf.k, 1:numMeas])
    ukf.xhats[ukf.ixp,:] .+= ukf.xhats[ukf.ixp - 1, :]

    # Compute covariance update
    ukf.P⁺ .= ukf.P⁻
    @views mul!(ukf.Pyt[1:7,1:numMeas], ukf.Pxy[1:7,1:numMeas], ukf.Py[1:numMeas,1:numMeas])
    @views mul!(ukf.P⁺, ukf.Pyt[1:7,1:numMeas], transpose(ukf.Pxy[1:7,1:numMeas]), -1.0, 1.0)

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

    # Update augmented covariance
    for r in 1:7
        for c in 1:7
            ukf.Pa[r,c] = ukf.P⁺[r,c]       
        end
    end
end

function updateIMU!(ukf::UKF)
    if ukf.updated == true
        throw(ErrorException("Cannot update again before propagating."))
    else
        ukf.updated = true 
    end

    @inbounds begin
    # Get current time 
    t   = ukf.ts[ukf.k]

    # Get true state 
    (ytrue, u)  = GetStateAndControl(ukf.scSim, t)

    # Generate accelerometer measurements 
    (dy, _us)   = GetStateAndControl(ukf.scSim, t; derivs = true)
    @views measAcc = ComputeMeasurement(ukf.imuSim, dy[4:6])

    # Compute expected measurements
    ukf.yhat    .= 0.0
    ukf.expmeas .= 0.0
    for i in 1:2*ukf.L + 1
        if ukf.lunaPerts == false
            dy = ukfNoLunarPertEOM(ukf.χ[:,i], ukf, t, ukf.yhat)
        else
            dy = ukfWithLunarPertEOM(ukf.χ[:,i], ukf, t, ukf.yhat)
        end
        @views expMeasAcc = ComputeMeasurement(ukf.imuSim, dy[4:6]; type = :expected)
        for j in 1:3
            ukf.expmeas[j,i] = expMeasAcc[j]
        end
    end

    # Number of measurements
    numMeas = 3

    # Compute weights
    w0m     = ukf.λ / (ukf.L + ukf.λ)
    w0c     = w0m + (1 - ukf.α^2 + ukf.β)
    wi      = 1.0 / (2.0*(ukf.L + ukf.λ))

    # Compute mean observation
    ukf.yhat .= 0.0
    for i in 1:2*ukf.L
        @views ukf.yhat[1:numMeas] .+= wi.*ukf.expmeas[1:numMeas,i]
    end
    @views ukf.yhat[1:numMeas] .+= w0m.*ukf.expmeas[1:numMeas,2*ukf.L + 1]

    # Compute residuals
    for i in 1:3
        ukf.rs[ukf.k, i] = measAcc[i] - ukf.yhat[i]
    end

    # Compute output covariance
    ukf.Py      .= 0.0
    for j in 1:2*ukf.L
        @views ukf.expmeasd[1:numMeas,j] .= ukf.expmeas[1:numMeas,j] .- ukf.yhat[1:numMeas]
        @views mul!(ukf.Py[1:numMeas,1:numMeas], ukf.expmeasd[1:numMeas,j], 
            transpose(ukf.expmeasd[1:numMeas,j]), wi, 1.0)
    end
    @views ukf.expmeasd[1:numMeas,2*ukf.L + 1] .= ukf.expmeas[1:numMeas,2*ukf.L + 1] .- ukf.yhat[1:numMeas]
    @views mul!(ukf.Py[1:numMeas,1:numMeas], ukf.expmeasd[1:numMeas,2*ukf.L + 1], 
        transpose(ukf.expmeasd[1:numMeas,2*ukf.L + 1]), w0c, 1.0)

    # Update output covariance to innovations covariance (assuming linearly appearing measurement noise)
    for i in 1:numMeas; ukf.Py[i,i] += ukf.σ2Acc; end

    # Compute cross correlation matrix
    ukf.Pxy .= 0.0
    for j in 1:2*ukf.L 
        @views mul!(ukf.Pxy[1:7,1:numMeas], ukf.χd[:,j], transpose(ukf.expmeasd[1:numMeas,j]), wi, 1.0)
    end
    @views mul!(ukf.Pxy[1:7,1:numMeas], ukf.χd[:,2*ukf.L + 1], transpose(ukf.expmeasd[1:numMeas,2*ukf.L + 1]), w0c, 1.0)

    # Compute gain 
    @views Pyfact = factorize(ukf.Py[1:numMeas,1:numMeas])
    @views rdiv!(ukf.Pxy[1:7,1:numMeas], Pyfact)

    # Update storage matricies
    # Time
    ukf.ixp += 1
    ukf.txp[ukf.ixp] = t

    # Compute state update and save
    @views mul!(ukf.xhats[ukf.ixp, :], ukf.Pxy[1:7,1:numMeas], ukf.rs[ukf.k, 1:numMeas])
    ukf.xhats[ukf.ixp,:] .+= ukf.xhats[ukf.ixp - 1, :]

    # Compute covariance update
    ukf.P⁺ .= ukf.P⁻
    @views mul!(ukf.Pyt[1:7,1:numMeas], ukf.Pxy[1:7,1:numMeas], ukf.Py[1:numMeas,1:numMeas])
    @views mul!(ukf.P⁺, ukf.Pyt[1:7,1:numMeas], transpose(ukf.Pxy[1:7,1:numMeas]), -1.0, 1.0)

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

    # Update augmented covariance
    for r in 1:7
        for c in 1:7
            ukf.Pa[r,c] = ukf.P⁺[r,c]       
        end
    end

    end
    nothing
end

function updateGPSIMU!(ukf::UKF)
    if ukf.updated == true
        throw(ErrorException("Cannot update again before propagating."))
    else
        ukf.updated = true 
    end

    @inbounds begin
    # Get current time 
    t   = ukf.ts[ukf.k]

    # Get true state 
    (ytrue, u) = GetStateAndControl(ukf.scSim, t)

    ## === GPS ===

    # Generate GPS measurements 
    @views measGPS = gpsMeasurement(ukf.gpsSim, ytrue[1:3], t)

    # Get satelite ids 
    @views sats = measGPS[2:2:end]
    numSats = length(sats)

    # Compute expected measurements
    ukf.expmeas .= 0.0
    for i in 1:2*ukf.L + 1
        # Get expected gps measurement
        @views (expMeasGPS, ps) = gpsMeasurement(ukf.gpsSim, ukf.χ[1:3,i], sats, t; type = :expected)
        for j in 1:numSats 
            ukf.expmeas[j,i] = expMeasGPS[3 + 2*(j - 1)]
        end
    end

    ## === IMU ===

    # Generate accelerometer measurement
    (dy, _us) = GetStateAndControl(ukf.scSim, t; derivs = true)
    @views measAcc = ComputeMeasurement(ukf.imuSim, dy[4:6])

    # Compute expected measurements
    ukf.yhat .= 0.0
    for i in 1:2*ukf.L + 1
        if ukf.lunaPerts == false
            dy = ukfNoLunarPertEOM(ukf.χ[:,i], ukf, t, ukf.yhat)
        else
            dy = ukfWithLunarPertEOM(ukf.χ[:,i], ukf, t, ukf.yhat)
        end
        @views expMeasAcc = ComputeMeasurement(ukf.imuSim, dy[4:6]; type = :expected)
        for j in 1:3
            ukf.expmeas[numSats + j,i] = expMeasAcc[j]
        end
    end   

    # Compute weights
    w0m     = ukf.λ / (ukf.L + ukf.λ)
    w0c     = w0m + (1 - ukf.α^2 + ukf.β)
    wi      = 1.0 / (2.0*(ukf.L + ukf.λ))

    # Compute mean observation
    ukf.yhat .= 0.0
    for i in 1:2*ukf.L
        @views ukf.yhat[1:numSats + 3] .+= wi.*ukf.expmeas[1:numSats + 3,i]
    end
    @views ukf.yhat[1:numSats + 3] .+= w0m.*ukf.expmeas[1:numSats + 3,2*ukf.L + 1]

    # Compute residuals with GPS pseudorange measurement rejection
    accMeas     = Vector{Bool}(undef, 35); accMeas .= false
    rejectTol   = 1.0
    numRejected = 0
    for i in 1:numSats 
        r = measGPS[3 + 2*(i - 1)] - ukf.yhat[i]
        if abs(r) < rejectTol
            accMeas[i] = true
            ukf.rs[ukf.k, i - numRejected] = r
        else
            numRejected += 1
        end
    end
    j = 0
    for i in numSats - numRejected + 1:numSats - numRejected + 3
        j += 1
        ukf.rs[ukf.k, i]    = measAcc[j] - ukf.yhat[i + numRejected]
        accMeas[i + numRejected]    = true
    end

    # Compute output covariance
    ukf.Py      .= 0.0
    numMeas  = numSats - numRejected + 3
    for j in 1:2*ukf.L
        @views ukf.expmeasd[1:numMeas,j] .= ukf.expmeas[accMeas,j] .- ukf.yhat[accMeas]
        @views mul!(ukf.Py[1:numMeas,1:numMeas], ukf.expmeasd[1:numMeas,j], 
            transpose(ukf.expmeasd[1:numMeas,j]), wi, 1.0)
    end
    @views ukf.expmeasd[1:numMeas,2*ukf.L + 1] .= ukf.expmeas[accMeas,2*ukf.L + 1] .- ukf.yhat[accMeas]
    @views mul!(ukf.Py[1:numMeas,1:numMeas], ukf.expmeasd[1:numMeas,2*ukf.L + 1], 
        transpose(ukf.expmeasd[1:numMeas,2*ukf.L + 1]), w0c, 1.0)

    # Update output covariance to innovations covariance (assuming linearly appearing measurement noise)
    for i in 1:numSats - numRejected; ukf.Py[i,i] += ukf.σ2GPS; end
    for i in numSats - numRejected + 1:numSats - numRejected + 3; ukf.Py[i,i] += ukf.σ2Acc; end

    # Compute cross correlation matrix
    ukf.Pxy .= 0.0
    for j in 1:2*ukf.L 
        @views mul!(ukf.Pxy[1:7,1:numMeas], ukf.χd[:,j], transpose(ukf.expmeasd[1:numMeas,j]), wi, 1.0)
    end
    @views mul!(ukf.Pxy[1:7,1:numMeas], ukf.χd[:,2*ukf.L + 1], transpose(ukf.expmeasd[1:numMeas,2*ukf.L + 1]), w0c, 1.0)

    # Compute gain 
    @views Pyfact = factorize(ukf.Py[1:numMeas,1:numMeas])
    @views rdiv!(ukf.Pxy[1:7,1:numMeas], Pyfact)

    # Update storage matricies
    # Time
    ukf.ixp += 1
    ukf.txp[ukf.ixp] = t

    # Compute state update and save
    @views mul!(ukf.xhats[ukf.ixp, :], ukf.Pxy[1:7,1:numMeas], ukf.rs[ukf.k, 1:numMeas])
    ukf.xhats[ukf.ixp,:] .+= ukf.xhats[ukf.ixp - 1, :]

    # Compute covariance update
    ukf.P⁺ .= ukf.P⁻
    @views mul!(ukf.Pyt[1:7,1:numMeas], ukf.Pxy[1:7,1:numMeas], ukf.Py[1:numMeas,1:numMeas])
    @views mul!(ukf.P⁺, ukf.Pyt[1:7,1:numMeas], transpose(ukf.Pxy[1:7,1:numMeas]), -1.0, 1.0)

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

    # Update augmented covariance
    for r in 1:7
        for c in 1:7
            ukf.Pa[r,c] = ukf.P⁺[r,c]       
        end
    end

    end
    nothing
end

function runFilter!(ukf::UKF)
    for i in 1:length(ukf.ts)
        propagate!(ukf)
        if i*ukf.imuΔt % ukf.gpsΔt == 0
            updateGPSIMU!(ukf)
        else
            updateIMU!(ukf)
        end
    end
end

function rkProp(f, x0, ts, ukf::UKF)
    h   = ts[2] - ts[1]
    χw  = view(x0, 8:14)
    xin = @SVector [x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7]]
    for i in 1:length(ts) - 1
        xi1 = @SVector [xin[1], xin[2], xin[3], xin[4], xin[5], xin[6], xin[7]]
        k1  = f(xi1, ukf, ts[i], χw)
        xi2 = @SVector [xi1[1] + 0.5*h*k1[1],
                        xi1[2] + 0.5*h*k1[2],
                        xi1[3] + 0.5*h*k1[3],
                        xi1[4] + 0.5*h*k1[4],
                        xi1[5] + 0.5*h*k1[5],
                        xi1[6] + 0.5*h*k1[6],
                        xi1[7] + 0.5*h*k1[7]]
        k2  = f(xi2, ukf, ts[i] + 0.5*h, χw)
        xi3 = @SVector [xi1[1] + 0.5*h*k2[1],
                        xi1[2] + 0.5*h*k2[2],
                        xi1[3] + 0.5*h*k2[3],
                        xi1[4] + 0.5*h*k2[4],
                        xi1[5] + 0.5*h*k2[5],
                        xi1[6] + 0.5*h*k2[6],
                        xi1[7] + 0.5*h*k2[7]]
        k3  = f(xi3, ukf, ts[i] + 0.5*h, χw)
        xi4 = @SVector [xi1[1] + h*k3[1],
                        xi1[2] + h*k3[2],
                        xi1[3] + h*k3[3],
                        xi1[4] + h*k3[4],
                        xi1[5] + h*k3[5],
                        xi1[6] + h*k3[6],
                        xi1[7] + h*k3[7]]
        k4  = f(xi4, ukf, ts[i] + h, χw)
        xin = @SVector [xi1[1] + h*(k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])/6.0,
                        xi1[2] + h*(k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2])/6.0,
                        xi1[3] + h*(k1[3] + 2.0*k2[3] + 2.0*k3[3] + k4[3])/6.0,
                        xi1[4] + h*(k1[4] + 2.0*k2[4] + 2.0*k3[4] + k4[4])/6.0,
                        xi1[5] + h*(k1[5] + 2.0*k2[5] + 2.0*k3[5] + k4[5])/6.0,
                        xi1[6] + h*(k1[6] + 2.0*k2[6] + 2.0*k3[6] + k4[6])/6.0,
                        xi1[7] + h*(k1[7] + 2.0*k2[7] + 2.0*k3[7] + k4[7])/6.0]
    end
    return xin
end

function ukfNoLunarPertEOM(y, ukf::UKF, t, χw)
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
    dy          = @SVector [y[4] + χw[1], y[5] + χw[2], y[6] + χw[3], 
                            nμr3*y[1] + at*u[2] + χw[4],
                            nμr3*y[2] + at*u[3] + χw[5],
                            nμr3*y[3] + at*u[4] + χw[6],
                            -u[1]*T / ukf.scSim.ps.sp.c + χw[7]]
                            
    return dy
end

function ukfWithLunarPertEOM(y, ukf::UKF, t, χw)
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
    dy          = @SVector [y[4] + χw[1], y[5] + χw[2], y[6] + χw[3], 
                            nμr3*y[1] + al[1] + at*u[2] + χw[4],
                            nμr3*y[2] + al[2] + at*u[3] + χw[5],
                            nμr3*y[3] + al[3] + at*u[4] + χw[6],
                            -u[1]*T / ukf.scSim.ps.sp.c + χw[7]]
    
    return dy
end

function plot(ukf::UKF, xtrue, n)
    ts   = ukf.txp[1:n]
    xhat = ukf.xhats[1:n, :]
    xt   = xtrue[1:n, :]
    es   = ukf.es[1:n, :]
    σ311 = 3*sqrt.(ukf.Ps[1:n,1])
    σ322 = 3*sqrt.(ukf.Ps[1:n,2])
    σ333 = 3*sqrt.(ukf.Ps[1:n,3])
    σ344 = 3*sqrt.(ukf.Ps[1:n,4])
    σ355 = 3*sqrt.(ukf.Ps[1:n,5])
    σ366 = 3*sqrt.(ukf.Ps[1:n,6])
    σ377 = 3*sqrt.(ukf.Ps[1:n,7])

    mat"""
    ts      = $ts;
    ts      = (ts - $(ukf.scSim.t0))/86400;
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
    legend('Estimation Error', 'UKF \$3\\sigma\$', 'Interpreter', 'latex')
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

function plotMC(ukf::UKF, em, ev, n)
    ts      = ukf.txp[1:n]
    xhat    = ukf.xhats[1:n, :]
    es      = ukf.es[1:n, :]
    σ311    = 3*sqrt.(ukf.Ps[1:n,1])
    σ311mc  = 3*sqrt.(ev[:,1])
    σ322    = 3*sqrt.(ukf.Ps[1:n,2])
    σ322mc  = 3*sqrt.(ev[:,2])
    σ333    = 3*sqrt.(ukf.Ps[1:n,3])
    σ333mc  = 3*sqrt.(ev[:,3])
    σ344    = 3*sqrt.(ukf.Ps[1:n,4])
    σ344mc  = 3*sqrt.(ev[:,4])
    σ355    = 3*sqrt.(ukf.Ps[1:n,5])
    σ355mc  = 3*sqrt.(ev[:,5])
    σ366    = 3*sqrt.(ukf.Ps[1:n,6])
    σ366mc  = 3*sqrt.(ev[:,6])
    σ377    = 3*sqrt.(ukf.Ps[1:n,7])
    σ377mc  = 3*sqrt.(ev[:,7])

    mat"""
    ts      = $ts;
    ts      = (ts - $(ukf.scSim.t0))/86400;
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

    figure()
    subplot(3,1,1)
    plot(ts, es(:,1), 'k')
    hold on
    plot(ts, s311, 'r')
    plot(ts, -s311, 'r')
    plot(ts, s311mc, '--r')
    plot(ts, -s311mc, '--r')
    legend('Estimation Error', 'UKF \$3\\sigma\$', 'Interpreter', 'latex')
    grid on

    subplot(3,1,2)
    plot(ts, es(:,1), 'k')
    hold on
    plot(ts, s322, 'r')
    plot(ts, -s322, 'r')
    plot(ts, s322mc, '--r')
    plot(ts, -s322mc, '--r')
    grid on

    subplot(3,1,3)
    plot(ts, es(:,3), 'k')
    hold on
    plot(ts, s333, 'r')
    plot(ts, -s333, 'r')
    plot(ts, s333mc, '--r')
    plot(ts, -s333mc, '--r')
    xlabel("Time, days")
    grid on

    figure()
    subplot(3,1,1)
    plot(ts, es(:,4), 'k')
    hold on
    plot(ts, s344, 'r')
    plot(ts, -s344, 'r')
    plot(ts, s344mc, '--r')
    plot(ts, -s344mc, '--r')

    subplot(3,1,2)
    plot(ts, es(:,5), 'k')
    hold on
    plot(ts, s355, 'r')
    plot(ts, -s355, 'r')
    plot(ts, s355mc, '--r')
    plot(ts, -s355mc, '--r')

    subplot(3,1,3)
    plot(ts, es(:,6), 'k')
    hold on
    plot(ts, s366, 'r')
    plot(ts, -s366, 'r')
    plot(ts, s366mc, '--r')
    plot(ts, -s366mc, '--r')

    figure()
    plot(ts, es(:,7), 'k')
    hold on
    plot(ts, s377, 'r')
    plot(ts, -s377, 'r')
    plot(ts, s377mc, '--r')
    plot(ts, -s377mc, '--r')

    """
end