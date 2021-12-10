
mutable struct UKF 
    # Measurement times
    ts::Vector{Float64}
    imuΔt::Float64
    gpsΔt::Float64

    # State estimates and covariance matrix diags
    txp::Vector{Float64} # Epochs associated with entries in xhats and Ps
    ixp::Int # Index associated with latest entry
    xhats::Matrix{Float64}
    es::Matrix{Float64}
    Ps::Matrix{Float64}

    # Parameters 
    α::Float64
    κ::Float64

    # Preallocate matricies
    Paacc::SparseMatrix{Float64, Int}
    Pagps::SparseMatrix{Float64, Int}
    χ::Matrix{Float64}
    C::Cholesky{Float64, Matrix{Float64}}
    P⁺::Matrix{Float64}
    P⁻::Matrix{Float64}

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

function UKF(xhat0, P0, Q, σ2GPS, σ2Acc, ts, gpsΔt, gpsSim, imuSim, scSim; 
    α = 1e-3, κ = 1.0, linaPerts = false, steps2save = 5)

    # Instantiate preallocated matricies
    Pa = sparse([P0 zeros(6,6); zeros(6,6) Q]) 

end