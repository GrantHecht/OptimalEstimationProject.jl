
# This function rotates igs GPS ephemeris from the ITRF14 reference frame to the GCRF reference frame
# using igs provided Earth rotation parameters

function rotateData(rawData::DataFrame, erps::Matrix)

    # Instantiate data matrix
    data = zeros(length(rawData.MJD), 10)

    # Read in nutation data
    #Ax = CSV.File(datadir("tab5")*"\\Ax.txt"; delim = " ", ignorerepeated = true, header = 0, types = Float64) |> Tables.matrix
    #Ay = CSV.File(datadir("tab5")*"\\Ay.txt"; delim = " ", ignorerepeated = true, header = 0, types = Float64) |> Tables.matrix
    #As = CSV.File(datadir("tab5")*"\\As.txt"; delim = " ", ignorerepeated = true, header = 0, types = Float64) |> Tables.matrix
    Ax = CSV.File(datadir("tab5", "Ax.txt"); delim = " ", ignorerepeated = true, header = 0, types = Float64) |> Tables.matrix
    Ay = CSV.File(datadir("tab5", "Ay.txt"); delim = " ", ignorerepeated = true, header = 0, types = Float64) |> Tables.matrix
    As = CSV.File(datadir("tab5", "As.txt"); delim = " ", ignorerepeated = true, header = 0, types = Float64) |> Tables.matrix

    # Loop through data and rotate
    @showprogress "Rotating ephemeris to J2000..." for i in 1:length(rawData.MJD)
        # ===== Get epoch in TAI 
        tai     = TAIEpoch(rawData.MJD[i]*days, origin=:modified_julian)

        # ===== First push data that does not require rotation to data matrix 
        data[i,1]   = tai.second
        data[i,2]   = parse(Float64, replace(string(rawData.Sat[i]), "G" => ""))
        data[i,6]   = rawData.Clk[i]
        data[i,10]  = rawData.Clkstd[i]

        # ===== Compute polar motion rotation matrix (W)
        ac = 0.26 # Chandler Wobble (arcseconds)
        aa = 0.12 # Annual Wobble (arcseconds)
        arcSec2Rad = π/(180.0*3600.0)

        # Get pole locations
        refMJD  = floor(data[i,1])
        dayFrac = data[i,1] % refMJD
        dayFrac < 0.5 ? (refMJD -= 0.5) : (refMJD += 0.5)
        refIdx  = findmin(abs.(view(erps,:,1) .- refMJD))[2]

        xp = (erps[refIdx,2] + dayFrac*erps[refIdx, 13]) * 10^-6 # (arcseconds)
        yp = (erps[refIdx,3] + dayFrac*erps[refIdx, 14]) * 10^-6 # (arcseconds)

        # Compute T_TT
        tt      = TTEpoch(tai)
        jdtt    = julian(tt)
        T_TT    = (jdtt.second / 86400.0 - 2451545.0) / 36525.0

        # Compute s'
        sp =  -0.0015*(ac^2/1.2 + aa^2)*T_TT

        # Compute W
        cxp = cos(xp*arcSec2Rad)
        sxp = sin(xp*arcSec2Rad)
        cyp = cos(yp*arcSec2Rad)
        syp = sin(yp*arcSec2Rad)
        csp = cos(sp*arcSec2Rad)
        ssp = sin(sp*arcSec2Rad)

        W11 = cxp*csp
        W12 = -cyp*ssp + syp*sxp*csp 
        W13 = -syp*ssp - cyp*sxp*csp
        W21 = cxp*ssp 
        W22 = cyp*csp + syp*sxp*ssp 
        W23 = syp*csp - cyp*sxp*ssp 
        W31 = sxp 
        W32 = -syp*cxp 
        W33 = cyp*cxp
        W   = @SMatrix [W11 W12 W13; W21 W22 W23; W31 W32 W33]

        # ===== Compute Earth rotation matrix (R)
        ut1 = UT1Epoch(tai)

        # Get Julian date 
        JDut1 = julian(ut1)
    
        θ = 2*π*(0.7790572732640 + 1.00273781191135448*(JDut1.second/86400.0 - 2451545.0))
        R = @SMatrix [ cos(-θ) sin(-θ) 0;
                      -sin(-θ) cos(-θ) 0;
                       0 0 1]

        # ===== Compute Precession-Nutation matrix (PN)
        PN = computePrecessionNutation(T_TT, Ax, Ay, As)

        # ===== Rotate Vectors
        T       = SMatrix{3,3}(PN*R*W) # Could improve this with BLAS
        r       = @SVector [rawData.X[i], rawData.Y[i], rawData.Z[i]]
        rstd    = @SVector [rawData.Xstd[i], rawData.Ystd[i], rawData.Zstd[i]]
        data[i, 3:5] .= T*r
        data[i, 7:9] .= T*rstd
    end

    return data
end