# This function reads in the necessary Earth orientation parameters for 
# rotating gps ephemeris from the ITRF14 reference frame to the GCRF
# (Geocentric Celestial Reference Frame) often just refered to as 
# the J2000 reference frame (i.e. in SPICE)
function readERPs(gpsWeekStart, weekDayStart, gpsWeekEnd, weekDayEnd)
    # If weekDayStart is 0, we need the final row of last weeks eop in addition to gpsWeekStart week
    needLastWeek = (weekDayStart == 0)

    # If weekDayEnd is 6, we need the first row in next weeks eop in addition to gpsWeekEnd week
    needNextWeek = (weekDayEnd == 6)

    # Compute number of weeks in requested data 
    numReqWeeks = gpsWeekEnd - gpsWeekStart + 1

    # Size EOP Matrix 
    eops = zeros(numReqWeeks*7 + Int(needLastWeek) + Int(needNextWeek), 18)

    # Read in data
    lineNum = 0

    # Read in last week data if necessary 
    if needLastWeek
        # Open file
        erpFile = datadir("igs") * (Sys.isunix() ? "/igs" : "\\igs") * string(gpsWeekStart - 1) * "7.erp"
        if !isfile(erpFile)
            throw(ErrorException("The file " * erpFile * " does not exist and must be downloaded from CDDIS.\n" * 
            "The file can be found at 'https://cddis.nasa.gov/archive/gnss/products/'"))
        end
        f = open(erpFile)

        # Get last line in file
        lastLine = ""
        for line in readlines(f)
            if length(line) < 2
                break
            end
            lastLine = line
        end

        # Put data in Matrix
        lineNum += 1
        splitLine = split(lastLine, " "; keepempty = false)
        for i in 1:18
            eops[lineNum, i] = parse(Float64, splitLine[i]) 
        end
        close(f)
    end

    # Read in span of weeks 
    weeks = gpsWeekStart:gpsWeekEnd
    for week in weeks
        # Open file
        erpFile = datadir("igs") * (Sys.isunix() ? "/igs" : "\\igs") * string(week) * "7.erp"
        if !isfile(erpFile)
            throw(ErrorException("The file 'igs"*string(week)*"7.erp' does not exist in the folder ./data/igs and must be downloaded from CDDIS.\n" * 
            "The file can be found at 'https://cddis.nasa.gov/archive/gnss/products/'"))
        end
        f = open(erpFile)

        # Read in data
        for line in readlines(f)
            if length(line) < 2
                break
            elseif isdigit(line[1])
                lineNum += 1
                splitLine = split(line, " "; keepempty = false)
                for i in 1:18
                   eops[lineNum, i] = parse(Float64, splitLine[i]) 
                end
            end
        end
        close(f)
    end

    # Read in last week data if necessary 
    if needNextWeek
        # Open file
        erpFile = datadir("igs") * (Sys.isunix() ? "/igs" : "\\igs") * string(gpsWeekEnd + 1) * "7.erp"
        if !isfile(erpFile)
            throw(ErrorException("The file 'igs"*string(gpsWeekEnd + 1)*"7.erp' does not exist in the folder ./data/igs and must be downloaded from CDDIS.\n" * 
            "The file can be found at 'https://cddis.nasa.gov/archive/gnss/products/'"))
        end
        f = open(erpFile)

        # Get first line of data in file
        firstDataLine = ""
        for line in readlines(f)
            if isdigit(line[1])
                firstDataLine = line
                break
            end
        end

        # Put data in Matrix
        lineNum += 1
        splitLine = split(firstDataLine, " "; keepempty = false)
        for i in 1:18
            eops[lineNum, i] = parse(Float64, splitLine[i]) 
        end
        close(f)
    end

    return eops
end