function genInterpolants(data)
    # Get rows of data
    rows    = size(data, 1)

    # Allocate matrix for satelite data
    satData = zeros(Int(floor(rows / 30)), 4)

    # Generate interpolants for each GPS Satelite
    local df
    @showprogress "Generating spline interpolants..." for i in 1:32
        # Cound epochs of data for sat i and put in satData
        count = 0
        idxs  = [1,3,4,5]
        for j in 1:rows
            if data[j, 2] == i
                count += 1
                @views satData[count, :] .= data[j, idxs]
            end
        end

        if i == 1
            @views x  = CubicSpline(satData[1:count,2], satData[1:count,1])
            @views y  = CubicSpline(satData[1:count,3], satData[1:count,1])
            @views z  = CubicSpline(satData[1:count,4], satData[1:count,1])

            df = DataFrame(SAT = i, X = [x], Y = [y], Z = [z])
        else
            @views x  = CubicSpline(satData[1:count,2], satData[1:count,1])
            @views y  = CubicSpline(satData[1:count,3], satData[1:count,1])
            @views z  = CubicSpline(satData[1:count,4], satData[1:count,1])
            push!(df, Dict(
                "SAT" => i, "X" => x, "Y" => y, "Z" => z
            ))
        end
    end

    return df
end