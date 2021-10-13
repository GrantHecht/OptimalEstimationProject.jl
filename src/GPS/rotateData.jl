
# This function rotates igs GPS ephemeris from the ITRF14 reference frame to the GCRF reference frame
# using igs provided Earth rotation parameters

function rotateData(rawData::DataFrame, erps::Matrix)

    # Instantiate data matrix
    data = zeros(length(rawData.MJD), 10)

    # Loop through data and rotate
    for i in 1:length(rawData.MJD)
        # ===== First push data that does not require rotation to data matrix 
        data[i,1]   = rawData.MDJ[i]
        data[i,2]   = parse(Float64,replace(rawData.Sat[i], "G" => ""))
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
        refIdx  = findmin(abs.(view(erps,:,1) - refMJD))[2]

        xp = (erps[refIdx,2] + dayFrac*erps[refIdx, 13]) * 10^-6 # (arcseconds)
        yp = (erps[refIDx,3] + dayFrac*erps[refIdx, 14]) * 10^-6 # (arcseconds)

        # Compute T_TT (See Valaddo p195 - This may not be exactly correct, but plenty correct for the project... Correct 
        # later if necessary)
        T_TT = (data[i,1] + (2400000.5 - 2451545.0)) / 36525.0

        # Compute s'
        sp =  -0.000047*T_TT # This is an approximation (See Valaddo p212 for exact)

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
        W   = @SVector [W11 W12 W13; W21 W22 W23; W31 W32 W33]

        # ===== Compute Earth rotation matrix (R)
        UTCmUTC = erps[refIdx, 4] + dayFrac*erps[refIdx, 5] # μs

        # Approximating here... need to more precicelly convert JDTT to JDUTC (if JDUTC contains leapseconds, need to include those)
        JDUTC1day = floor(data[i,1])
        JDUTC1frac = data[i,1] % JDUTC1day + UTCmUTC*1e-6/86400
        if JDUTC1frac > 1
            JDUTC1day += floor(JDUTC1frac)
            JDUTC1frac -= floor(JDUTC1frac)
        end
    
        θ = 2*π*(0.7790572732640 + 1.00273781191135448*(JDUT1day - 2451545) + 1.00273781191135448*JDUT1frac)
        R = @SVector [ cos(-θ) sin(-θ) 0;
                      -sin(-θ) cos(-θ) 0;
                      0 0 1]

        # ===== Compute Precession-Nutation matrix (PN)

        # DO THE ROTATION!!!!!!!
    end

    return data
end