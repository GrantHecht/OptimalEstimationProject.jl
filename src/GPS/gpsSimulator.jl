
mutable struct gpsSim 
    # Ephemeris data initial and final epochs  
    startGPSWeek::Int
    startWeekDay::Int 
    endGPSWeek::Int
    endWeekDay::Int

    # Raw gps data with ephemeris in IGS reference frame (currently IGS14)
    rawSp3Data::DataFrame

    # GPS Processed Data
    data::Matrix{Float64}

    # ERPs
    erps::Matrix{Float64}

    function gpsSim(startGPSWeek::Integer, startWeekDay::Integer, endGPSWeek::Integer, endWeekDay::Integer)

        # Read in sp3 unprocessed data
        rawSp3Data = readSP3s(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay)

        # Read in Earth Orientation Parameters
        erps = readERPs(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay)

        # Rotate data from ITRF14 to GCRF and data required for gps simulation in fast to read matrix
        #data = rotateData(rawsp3Data, eops)
        data = zeros(1,1)

        new(startGPSWeek, startWeekDay, endGPSWeek, endWeekDay, rawSp3Data, data, erps)
    end
end