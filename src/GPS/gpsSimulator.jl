
mutable struct GPSSim
    # Ephemeris data initial and final epochs  
    startGPSWeek::Int
    startWeekDay::Int 
    endGPSWeek::Int
    endWeekDay::Int

    # Raw gps data with ephemeris in IGS reference frame (currently IGS14)
    rawSp3Data::DataFrame

    # GPS Processed Data
    data::Matrix{Float64}

    # Interpolators
    interps::DataFrame

    # ERPs
    erps::Matrix{Float64}

    # Measurement statistics
    σρ::Float64
    σr::Float64

    function GPSSim(startGPSWeek::Integer, startWeekDay::Integer, endGPSWeek::Integer, endWeekDay::Integer; 
                    σρ::Float64 = 1.0e-3, σr::Float64 = 5e-3)
        # Read in sp3 unprocessed data
        rawSp3Data = readSP3s(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay)

        # Read in Earth Orientation Parameters
        erps = readERPs(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay)

        # Rotate data from ITRF14 to GCRF and data required for gps simulation in fast to read matrix
        data = rotateData(rawSp3Data, erps)

        # Create interpolators
        interps = genInterpolants(data)

        new(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay, rawSp3Data, data, interps, erps, σρ, σr)
    end
end

function gpsMeasurement(gps::GPSSim, r::AbstractVector, t; type = :true)
    # Grab vectors of interpolators
    Xs = gps.interps.X 
    Ys = gps.interps.Y
    Zs = gps.interps.Z

    # Compute time of transmition from each satelite
    tts = zeros(32)
    for i in 1:32
        try
            func(tt) = transmissionTimeFunction(tt, t, r, Xs[i], Ys[i], Zs[i])
            tts[i] = find_zero(func, t)
            #tts[i] = transmissionTimeFunction(t, r, Xs[i], Ys[i], Zs[i])
        catch
            tts[i] = 0.0
        end
    end

    # Determine visable GPS satelites
    visable = Vector{Bool}(undef, 32); visable .= false
    for i in 1:32
        if tts[i] != 0
            # Get vector to GPS Sat if
            rsat = @SVector [Xs[i](tts[i]), Ys[i](tts[i]), Zs[i](tts[i])]

            # Compute vector from r to rsat 
            rdiff = @SVector [rsat[1] - r[1], rsat[2] - r[2], rsat[3] - r[3]]

            # Check for intersection with sphere by solving quadratic equation
            nrDiff = norm(rdiff)
            u = @SVector [rdiff[1] / nrDiff, rdiff[2] / nrDiff, rdiff[3] / nrDiff]
            b = 2.0*transpose(u)*r
            c = transpose(r)*r - 6378.137^2
            d = b^2 - 4*c
            # If d is negative, rdiff does not intersect the sphere of the Earth
            if d < 0.0
                visable[i] = true
            end
        end
    end

    # Generate measurements from visable satelites
    meas    = zeros(2*sum(visable) + 1)
    meas[1] = t
    count   = 0
    for i in 1:32
        if visable[i] == true
            count += 1

            # Set sat num
            meas[2*(count - 1) + 2] = i

            # Get vector to GPS Sat
            if type == :true # Use true position of satelite
                rsat = @SVector [Xs[i](tts[i]), Ys[i](tts[i]), Zs[i](tts[i])]
            else # If generating expected measurements, corrupt position to simulate ArgumentError
                 # in GPS broadcast ephemeris
                rsat = @SVector [Xs[i](tts[i]) + gps.σr*randn(), 
                                 Ys[i](tts[i]) + gps.σr*randn(), 
                                 Zs[i](tts[i]) + gps.σr*randn()]
            end

            # Compute vector from r to rsat 
            rdiff = @SVector [rsat[1] - r[1], rsat[2] - r[2], rsat[3] - r[3]]

            # Get GPS clock bias and corrupt with noise
            Δtsat = 0.0
            Δtm = 0.0
            n, m  = size(gps.data)
            for j in 32 + i:31:n
                if gps.data[j - 31,1] < t && gps.data[j, 1] > t 
                    #t1      = gps.data[j - 31,1]
                    #t2      = gps.data[j,1]
                    Δt1     = gps.data[j - 31,6] 
                    #Δt2     = gps.data[j,6]
                    Δtσ1    = gps.data[j - 31,10]
                    #Δtσ2    = gps.data[j,10]
                    #Δttrue  = (Δt1 + (t - t1)*(Δt2 - Δt1)/(t2 - t1)) * 1e-6
                    #Δtσ     = Δtσ1 + (t - t1)*(Δtσ2 - Δtσ1)/(t2 - t1) * 1e-12
                    Δtm     = Δt1*1e-6
                    Δtσ     = Δtσ1*1e-12
                    Δtsat   = Δtm + Δtσ*randn()
                    if isnan(Δtsat) # Currently data in GPSSim.data is getting set to NaN so this is a workaround for now.
                        Δtsat = 0.0
                    end
                end
            end

            # Compute pseudorange
            if type == :true
                meas[2*(count - 1) + 3] = norm(rdiff) - 299792.458*Δtsat + gps.σρ*randn()
            else
                meas[2*(count - 1) + 3] = norm(rdiff) - 299792.458*Δtm
            end
        end
    end

    return meas
end

function gpsMeasurement(gps::GPSSim, r::AbstractVector, gpsSats, t; type = :true)
    # Grab vectors of interpolators
    Xs = gps.interps.X 
    Ys = gps.interps.Y
    Zs = gps.interps.Z

    # Compute time of transmition from each satelite
    tts = zeros(32)
    for ifloat in gpsSats
        i = Int(ifloat)
        try
            #func(tt) = transmissionTimeFunction(tt, t, r, Xs[i], Ys[i], Zs[i])
            #tts[i] = find_zero(func, t)
            tts[i] = transmissionTimeFunction(t,r,Xs[i],Ys[i],Zs[i])
        catch
            tts[i] = 0.0
        end
    end

    # Generate measurements from visable satelites
    meas    = zeros(2*length(gpsSats) + 1)
    ps      = zeros(length(gpsSats), 3)
    meas[1] = t
    count   = 0
    for ifloat in gpsSats
        i = Int(ifloat)
        count += 1

        # Set sat num
        meas[2*(count - 1) + 2] = i

        # Get vector to GPS Sat
        if type == :true # Use true position of satelite
            rsat = @SVector [Xs[i](tts[i]), Ys[i](tts[i]), Zs[i](tts[i])]
        else # If generating expected measurements, corrupt position to simulate ArgumentError
                # in GPS broadcast ephemeris
            rsat = @SVector [Xs[i](tts[i]) + gps.σr*randn(), 
                             Ys[i](tts[i]) + gps.σr*randn(), 
                             Zs[i](tts[i]) + gps.σr*randn()]
        end

        # Compute vector from r to rsat 
        rdiff = @SVector [rsat[1] - r[1], rsat[2] - r[2], rsat[3] - r[3]]

        # Get GPS clock bias and corrupt with noise
        Δtsat = 0.0
        Δtm   = 0.0
        n, m  = size(gps.data)
        for j in 32 + i:31:n
            if gps.data[j - 31,1] < t && gps.data[j, 1] > t 
                #t1      = gps.data[j - 31,1]
                #t2      = gps.data[j,1]
                Δt1     = gps.data[j - 31,6] 
                #Δt2     = gps.data[j,6]
                Δtσ1    = gps.data[j - 31,10]
                #Δtσ2    = gps.data[j,10]
                #Δttrue  = (Δt1 + (t - t1)*(Δt2 - Δt1)/(t2 - t1)) * 1e-3
                #Δtσ     = Δtσ1 + (t - t1)*(Δtσ2 - Δtσ1)/(t2 - t1) * 1e-12
                Δtm     = Δt1*1e-6
                Δtσ     = Δtσ1*1e-12
                Δtsat   = Δtm + Δtσ*randn()
            end
        end

        # Compute pseudorange
        if type == :true
            meas[2*(count - 1) + 3] = norm(rdiff) - 299792.458*Δtsat + gps.σρ*randn()
        else
            meas[2*(count - 1) + 3] = norm(rdiff) - 299792.458*Δtm
        end

        # Set GPS Satelite position
        ps[count, :] .= rdiff
    end

    return (meas, ps)
end

function plotGPS(gps::GPSSim, sats::AbstractVector)
    mat"""
    data = $(gps.data);
    sats = $sats;
    satData = nan(floor(size(data,1)/30),3);

    figure()
    hold on
    for i = 1:length(sats)
        satData = satData .* NaN;
        sat     = sats(i);
        count   = 0;
        for j = 1:size(data,1)
            if data(j,2) == sat
                count = count + 1;
                satData(count,1:3)  = data(j,3:5);
            end
        end
        plot3(satData(:,1),satData(:,2),satData(:,3))
    end

    % Plot Earth
    [A,B,C] = ellipsoid(0.0, 0.0, 0.0, 6378, 6378, 6378);
    surf(A,B,-C, imread("./figures/WM.jpg"), "FaceColor", "texturemap", "EdgeColor", "none")

    grid on
    axis equal
    view(3)

    xlabel("X, km", "Interpreter", "latex")
    ylabel("Y, km", "Interpreter", "latex")
    zlabel("Z, km", "Interpreter", "latex")

    set(gca, "fontname", "Times New Roman", "fontsize", 10)
    set(gcf, "PaperUnits", "inches", "PaperPosition", [0.25, 0.25, 5.0, 4.0])
    set(gcf, "PaperPositionMode", "Manual")
    print("./figures/GPSConst.eps", "-depsc", "-r300")
    """
end

plotGPS(gps::GPSSim) = plotGPS(gps, 1:32)