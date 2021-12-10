
function readSP3(gpsWeek, day)

    # Construct string
    #sp3File = datadir("igs") * "\\igs" * string(gpsWeek) * string(day) * ".sp3"
    sp3File = datadir("igs", "igs" * string(gpsWeek) * string(day) * ".sp3")

    if !isfile(sp3File)
        throw(ErrorException("The file " * sp3File * " does not exist and must be downloaded from CDDIS.\n" * 
        "The file can be found at 'https://cddis.nasa.gov/archive/gnss/products/'"))
    end

    # Instantiate Data Frame 
    df = DataFrame(
        MJD     = Float64[],
        Year    = Int[],
        Month   = Int[],
        Day     = Int[],
        Hour    = Int[],
        Min     = Int[],
        Sec     = Float64[],
        Sat     = Symbol[],
        X       = Float64[],
        Y       = Float64[],
        Z       = Float64[],
        Clk     = Float64[],
        Xstd    = Float64[],
        Ystd    = Float64[],
        Zstd    = Float64[],
        Clkstd  = Float64[]
        )

    # Open file 
    f = open(sp3File)

    # Base std errors
    basePos = 0.0
    baseClk = 0.0

    # Date and time 
    MJD     = 0.0
    Year    = 0
    Month   = 0
    Day     = 0
    Hour    = 0
    Min     = 0
    Sec     = 0.0

    # Begin looping through lines in file 
    lineNum = 0
    for line in readlines(f)
        lineNum += 1
        if lineNum > 22
            if line[1] == '*' # Line designates time
                Year    = parse(Int, line[4:7])
                Month   = parse(Int, line[9:10])
                Day     = parse(Int, line[12:13])
                Hour    = parse(Int, line[15:16])
                Min     = parse(Int, line[18:19])
                Sec     = parse(Float64, line[21:31])
            elseif line[1] == 'P'
                Sat         = Symbol(line[2:4])
                X           = parse(Float64, line[5:18])
                Y           = parse(Float64, line[19:32])
                Z           = parse(Float64, line[33:46])
                Clk         = parse(Float64, line[47:60])
                if length(line) > 60
                    XstdExp     = tryparse(Int, line[62:63])
                    XstdExp === nothing ? XstdExp = NaN : ()
                    YstdExp     = tryparse(Int, line[65:66])
                    YstdExp === nothing ? YstdExp = NaN : ()
                    ZstdExp     = tryparse(Int, line[68:69])
                    ZstdExp === nothing ? ZstdExp = NaN : ()
                    ClkstdExp   = tryparse(Int, line[71:73])
                    ClkstdExp === nothing ? ClkstdExp = NaN : ()
                else
                    XstdExp     = NaN
                    YstdExp     = NaN
                    ZstdExp     = NaN
                    ClkstdExp   = NaN
                end

                # Compute MJD 
                dayWithFrac = day + ((Sec/60.0 + Min)/60.0 + Hour)/24.0
                MJD = gps2MJD(gpsWeek, dayWithFrac)

                # Push to Data Frame
                push!(df, Dict(
                    "MJD"       => MJD,
                    "Year"      => Year,
                    "Month"     => Month,
                    "Day"       => Day,
                    "Hour"      => Hour,
                    "Min"       => Min,
                    "Sec"       => Sec,
                    "Sat"       => Sat,
                    "X"         => X,
                    "Y"         => Y,
                    "Z"         => Z,
                    "Clk"       => Clk,
                    "Xstd"      => basePos^XstdExp,
                    "Ystd"      => basePos^YstdExp,
                    "Zstd"      => basePos^ZstdExp,
                    "Clkstd"    => baseClk^ClkstdExp
                ))
            end
        else
            if lineNum == 15
                str = split(line, " "; keepempty = false)[2:3]
                basePos = parse(Float64, str[1])
                baseClk = parse(Float64, str[2])
            end
        end
    end

    close(f)
    return df
end

function readSP3s(gpsWeekStart, weekDayStart, gpsWeekEnd, weekDayEnd)
    # Handle gps weeks and week days
    gpsWeeks        = gpsWeekStart:gpsWeekEnd
    startWeekDays   = weekDayStart:6
    endWeekDays     = 0:weekDayEnd

    # Instantiate Data Frame 
    df = DataFrame(
        MJD     = Float64[],
        Year    = Int[],
        Month   = Int[],
        Day     = Int[],
        Hour    = Int[],
        Min     = Int[],
        Sec     = Float64[],
        Sat     = Symbol[],
        X       = Float64[],
        Y       = Float64[],
        Z       = Float64[],
        Clk     = Float64[],
        Xstd    = Float64[],
        Ystd    = Float64[],
        Zstd    = Float64[],
        Clkstd  = Float64[]
        )

    @showprogress "Reading in sp3 files..." for week in gpsWeeks 
        if week == gpsWeekStart 
            for day in startWeekDays
                df = vcat(df, readSP3(week, day))
            end
        elseif week == gpsWeekEnd 
            for day in endWeekDays
                df = vcat(df, readSP3(week, day))
            end
        else
            for day in 1:6
                df = vcat(df, readSP3(week, day))
            end
        end
    end

    return df
end