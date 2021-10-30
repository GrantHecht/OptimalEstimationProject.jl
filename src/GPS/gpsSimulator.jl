
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

    function GPSSim(startGPSWeek::Integer, startWeekDay::Integer, endGPSWeek::Integer, endWeekDay::Integer)
        # Read in sp3 unprocessed data
        rawSp3Data = readSP3s(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay)

        # Read in Earth Orientation Parameters
        erps = readERPs(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay)

        # Rotate data from ITRF14 to GCRF and data required for gps simulation in fast to read matrix
        data = rotateData(rawSp3Data, erps)

        # Create interpolators
        interps = genInterpolants(data)

        new(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay, rawSp3Data, data, interps, erps)
    end
end

function gpsMeasurement(gps::GPSSim, r::AbstractVector, t, σρ)
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
        else
            visable[i] = false
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

            # Get vector to GPS Sat if
            rsat = @SVector [Xs[i](tts[i]), Ys[i](tts[i]), Zs[i](tts[i])]

            # Compute vector from r to rsat 
            rdiff = @SVector [rsat[1] - r[1], rsat[2] - r[2], rsat[3] - r[3]]

            # Get GPS clock bias and corrupt with noise
            Δtsat = 0.0
            n, m  = size(gps.data)
            for j in 32 + i:31:n
                if gps.data[j - 31,1] < t && gps.data[j, 1] > t 
                    t1      = gps.data[j - 31,1]
                    t2      = gps.data[j,1]
                    Δt1     = gps.data[j - 31,6] 
                    Δt2     = gps.data[j,6]
                    Δtσ1    = gps.data[j - 31,10]
                    Δtσ2    = gps.data[j,10]
                    Δttrue  = (Δt1 + (t - t1)*(Δt2 - Δt1)/(t2 - t1)) * 1e-3
                    Δtσ     = Δtσ1 + (t - t1)*(Δtσ2 - Δtσ1)/(t2 - t1) * 1e-12
                    Δtsat   = Δttrue + Δtσ*randn()
                end
            end

            # Compute pseudorange
            meas[2*(count - 1) + 3] = norm(rdiff) - 299792.458*Δtsat
        end
    end

    return meas
end

function plotGPS(gps::GPSSim, sats::AbstractVector)
    mat"""
    data = $(gps.data);
    sats = $sats;
    satData = nan(floor(size(data,1)/30),3);

    figure()
    hold on
    for i = 1:length(sats)
        sat   = sats(i);
        count = 0;
        for j = 1:size(data,1)
            if data(j,2) == sat
                count = count + 1;
                satData(count,1:3)  = data(j,3:5);
            end
        end
        plot3(satData(:,1),satData(:,2),satData(:,3))
    end
    grid on
    axis equal
    """
end