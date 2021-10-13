# This function converts GPS week and day to MJD
function gps2MJD(week, day)
    # GPS Epoch in Mod. Julian Date
    gpsMJDEpoch = 2444244.5 - 2400000.5

    # Add GPS week and days
    MJD = gpsMJDEpoch + week*7.0 + day 

    return MJD
end