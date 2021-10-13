
function downloadGPSData(gpsWeek::Int, day::Int)
    baseUrl = "https://cddis.nasa.gov/archive/gnss/products/"
    sp3Url  = baseUrl * string(gpsWeek) * "\\igs" * string(gpsWeek) * string(day) * ".sp3"

    fileName = datadir("igs") * "\\igs" * string(gpsWeek) * string(day) * ".sp3"
    download(sp3Url, fileName*".Z")
    
    try
        stream = DeflateDecompressorStream(open(fileName*".Z"))
        for line in eachline(stream)

        end
    finally
        close(stream)
    end
end