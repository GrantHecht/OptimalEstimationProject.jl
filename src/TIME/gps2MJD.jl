# This function converts GPS week and day to T_TAI (MJD referenced to International Atomic Time)
function gps2MJD(week, day)
    # GPS Epoch in Mod. Julian Date
    gpsMJDEpoch = 2444244.5 - 2400000.5

    # Add GPS week and days
    MJD = gpsMJDEpoch + week*7.0 + day + 19.0 / 86400.0

    return MJD
end