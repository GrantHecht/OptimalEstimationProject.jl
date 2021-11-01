mutable struct SpacecraftSim{DEST}
    # Initial Epoch (TAI Seconds)
    t0::Float64

    # Initial Co-States Defining Primer Vector Control 
    λ0::Vector{Float64}

    # Time span of trajectory
    ts::Tuple{Float64, Float64}

    # Indirect optimal control spacecraft parameters
    ps::IndirectTrajOpt.CR3BPIndirectParams

    # DifferentialEquations.jl solution object (providing interpolated traj.)
    interp::DEST

    # Constructor
    function SpacecraftSim(t0)
        # Initializing time span, initial co-states, and parameters
        λ0      = [ 5.963713007227692, 12.999481348176722, 1.6584130970663207, 
                    -0.04062357612068792, 0.018558730104232632, -0.0004993971529784666, 0.12200533102974148]
        ps      = initCR3BPIndirectParams("Low Thrust 10 CR3BP"); ps.ϵ = 0.0
        ts      = (0.0, 33.1) # [days]
        tsND    = (0.0, ts[2] * 86400 / ps.crp.TU) # [n.d. time units]

        # Boundary conditions and initial state/co-state vector
        μ   = ps.crp.μ; LU = ps.crp.LU; TU = ps.crp.TU
        ics = [-0.019488511458668, -0.016033479812051, 0.0, 8.918881923678198, -4.081793688818725, 0.0, 1.0] 
        fcs = [(1.0 - μ) + 6060.483/LU, 19452.284/LU, -34982.968/LU, 
                0.082677*(TU/LU), 0.006820*(TU/LU), -0.368434*(TU/LU), 0.0]
        y0  = @SVector [ics[1], ics[2], ics[3], ics[4], ics[5], ics[6], ics[7], 
                        λ0[1],  λ0[2],  λ0[3],  λ0[4],  λ0[5],  λ0[6],  λ0[7]]

        # Integrate state and co-state differential equations 
        sol = integrate(y0, tsND, ps, CR3BP(), FullSolutionHistory(), MEMF())

        # Create spacecraft sim object
        new{typeof(sol)}(t0, λ0, ts, ps, sol)
    end
end

function GetStateAndControl(ss::SpacecraftSim, t::Real; frame = :inertial)
    # Check that t ∈ t0 <= t <= t0 + ts[2]
    if t < ss.t0 || t > ss.t0 + ss.ts[2]*86400
        throw(ArgumentError("t does not fall within tranjectory time span."))
    end

    # Compute non-dimentioanl time
    tnd     = (t - ss.t0) / ss.ps.crp.TU

    # Get state in synodic reference frame 
    ysyn    = ss.interp(tnd)

    # Scale physical states to units of km, km/s, and kg
    xsyns   = @SVector [(ysyn[1] + ss.ps.crp.μ)*ss.ps.crp.LU, ysyn[2]*ss.ps.crp.LU, ysyn[3]*ss.ps.crp.LU,
                        ysyn[4]*ss.ps.crp.VU, ysyn[5]*ss.ps.crp.VU, ysyn[6]*ss.ps.crp.VU, ysyn[7]*ss.ps.crp.MU]

    # Rotate to inertial if desired
    if frame == :inertial
        # Create rotation matrix to inertial frame
        ω       = -1.0              # [rad / n.d. time]
        ωs      = ω / ss.ps.crp.TU  # [rad / second]
        θ       = ωs*(ss.t0 + ss.ts[2]*86400 - t)    # [rad]
        RsI     = @SMatrix [cos(θ) -sin(θ) 0;
                            sin(θ)  cos(θ) 0;
                            0       0      1]
        dRsI    = @SMatrix [-ωs*sin(θ) -ωs*cos(θ) 0;
                            ωs*cos(θ) -ωs*sin(θ) 0;
                            0         0        0]

        # Rotate to inertial frame
        xs      = zeros(7); xs[7] = xsyns[7]
        @views mul!(xs[1:3], RsI,  xsyns[1:3])
        @views mul!(xs[4:6], RsI,  xsyns[4:6])
        @views mul!(xs[4:6], dRsI, xsyns[1:3], 1.0, 1.0)
    else
        xs = xsyns
    end

    # Compute controls
    us      = zeros(4);
    
    # Scaled exaust velocity
    cs      = ss.ps.sp.c * ss.ps.crp.TU / (ss.ps.crp.LU * 1000.0)

    # Switching function
    λv      = norm(view(ysyn, 11:13))
    S       = IndirectTrajOpt.computeS(ysyn, λv, cs)

    # Throttling factor 
    us[1]   = S > 0.0 ? 0.0 : 1.0

    # Thrust direction
    ussyn   = @SVector [-ysyn[11] / λv, -ysyn[12] / λv, -ysyn[13] / λv]

    # Rotate to necisary if desired
    if frame == :inertial
        @views mul!(us[2:4], RsI, ussyn)
    else
        us[2:4] .= ussyn
    end
    
    # Return state and control
    return (xs, us)
end

function GetStateAndControl(ss::SpacecraftSim, ts::AbstractVector; frame = :inertial)
    xIs = zeros(length(ts), 7)
    us  = zeros(length(ts), 4)
    for i in 1:length(ts)
        # Get time 
        t = ts[i]

        # Get state and control
        sc = GetStateAndControl(ss, t; frame = frame)

        # Place in data matricies
        xIs[i,:] .= sc[1]
        us[i,:]  .= sc[2]
    end

    return (xIs, us)
end

function GetControl(ss::SpacecraftSim, t::Real; frame = :inertial)
    # Check that t ∈ t0 <= t <= t0 + ts[2]
    if t < ss.t0 || t > ss.t0 + ss.ts[2]*86400
        throw(ArgumentError("t does not fall within tranjectory time span."))
    end

    # Compute non-dimentioanl time
    tnd     = (t - ss.t0) / ss.ps.crp.TU

    # Get state in synodic reference frame 
    ysyn    = ss.interp(tnd)

    # Compute controls
    us      = zeros(4);
    
    # Scaled exaust velocity
    cs      = ss.ps.sp.c * ss.ps.crp.TU / (ss.ps.crp.LU * 1000.0)

    # Switching function
    λv      = norm(view(ysyn, 11:13))
    S       = IndirectTrajOpt.computeS(ysyn, λv, cs)

    # Throttling factor 
    us[1]   = S > 0.0 ? 0.0 : 1.0

    # Thrust direction
    ussyn   = @SVector [-ysyn[11] / λv, -ysyn[12] / λv, -ysyn[13] / λv]

    # Rotate to necisary if desired
    if frame == :inertial
        # Create rotation matrix to inertial frame
        ω       = -1.0              # [rad / n.d. time]
        ωs      = ω / ss.ps.crp.TU  # [rad / second]
        θ       = ωs*(ss.t0 + ss.ts[2]*86400 - t)    # [rad]
        RsI     = @SMatrix [cos(θ) -sin(θ) 0;
                            sin(θ)  cos(θ) 0;
                            0       0      1]

        @views mul!(us[2:4], RsI, ussyn)
    else
        us[2:4] .= ussyn
    end
    
    # Return state and control
    return us
end

function PlotTrajectory(ss::SpacecraftSim, ts; frame = :inertial)
    # Get states and controls
    (xs,us) = GetStateAndControl(ss, ts; frame = frame)

    # Get target trajectory 
    xt = generateTargetTrajectory()

    # Plot with MATLAB 
    mat"""
        % Bulk Data
        ts = $ts;
        xs = $xs;
        us = $us;
        xh = $xt;

        % Split into thrust and coast
        xc = nan(size(xs));
        xt = nan(size(xs));
        ut = nan(size(us));
        lu = us(1,1);
        for i = 1:length(ts)
            if us(i,1) == 1.0 || (us(i,1) == 0.0 && lu == 1.0)
                xt(i,:) = xs(i,:);
                ut(i,:) = us(i,:);
            end 
            if us(i,1) == 0.0 || (us(i,1) == 1.0 && lu == 0.0)
                xc(i,:) = xs(i,:);
            end
            lu = us(i,1);
        end

        figure()
        plot3(xt(:,1), xt(:,2), xt(:,3), 'r')
        hold on
        plot3(xc(:,1), xc(:,2), xs(:,3), 'b')
        %quiver3(xt(:,1), xt(:,2), xt(:,3), ut(:,2), ut(:,3), ut(:,4), 'r')
        plot3(xh(:,1), xh(:,2), xh(:, 3), '--k')

        axis equal 
        grid on
    """
end

function PlotTrajectory(ss::SpacecraftSim; frame = :inertial)
    # Plot full trajectory
    ts = range(ss.t0; stop = ss.t0 + ss.ts[2]*86400, length = 10000)
    PlotTrajectory(ss, ts; frame = frame)
end

function generateTargetTrajectory()
    # Parameters 
    ps = initCR3BPIndirectParams("Low Thrust 10 CR3BP")
    μ  = ps.crp.μ
    LU = ps.crp.LU
    TU = ps.crp.TU

    # Target final conditions
    xf = @SVector [(1.0 - μ) + 6060.483/LU, 19452.284/LU, -34982.968/LU, 0.082677*(TU/LU), 0.006820*(TU/LU), -0.368434*(TU/LU)]

    # Integrate target trajectory
    solTarg = integrate(xf, (0.0, 1.6), "Low Thrust 10 CR3BP", CR3BP(), FullSolutionHistoryNoControl())
    rt = zeros(length(solTarg.t), 3)
    for i in 1:length(solTarg.t)
        rt[i,1] = (solTarg.u[i][1] + μ)*LU
        rt[i,2] = solTarg.u[i][2]*LU
        rt[i,3] = solTarg.u[i][3]*LU
    end

    return rt
end