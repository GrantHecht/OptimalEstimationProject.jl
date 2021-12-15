
function computeMCStats(folderName)
    # Get files in folder 
    files = readdir(folderName)

    # Read first file to get time
    f  = CSV.File(folderName * "/" * files[1])
    ts = Vector{Float64}(undef, 0) 
    for row in f
        push!(ts, row.Time)    
    end

    # Create matricies for data
    erx     = zeros(length(ts), length(files)) 
    ery     = zeros(length(ts), length(files)) 
    erz     = zeros(length(ts), length(files)) 
    evx     = zeros(length(ts), length(files)) 
    evy     = zeros(length(ts), length(files)) 
    evz     = zeros(length(ts), length(files)) 
    em      = zeros(length(ts), length(files)) 

    # Load data
    for i in 1:length(files)
        f = CSV.File(folderName * "/" * files[i])
        j = 0
        for row in f
            j += 1
            erx[j, i] = row.erx
            ery[j, i] = row.ery
            erz[j, i] = row.erz
            evx[j, i] = row.evx
            evy[j, i] = row.evy
            evz[j, i] = row.evz
            em[j, i]  = row.em
        end
    end

    # Compute mean and variance of MC errors
    errMean = zeros(length(ts), 7)
    errVar  = zeros(length(ts), 7)
    for i in 1:length(ts)
        @views errMean[i,1] = mean(erx[i,:])
        @views errMean[i,2] = mean(ery[i,:])
        @views errMean[i,3] = mean(erz[i,:])
        @views errMean[i,4] = mean(evx[i,:])
        @views errMean[i,5] = mean(evy[i,:])
        @views errMean[i,6] = mean(evz[i,:])
        @views errMean[i,7] = mean(em[i,:])
        @views errVar[i,1]  = var(erx[i,:]; mean = errMean[i,1])
        @views errVar[i,2]  = var(ery[i,:]; mean = errMean[i,2])
        @views errVar[i,3]  = var(erz[i,:]; mean = errMean[i,3])
        @views errVar[i,4]  = var(evx[i,:]; mean = errMean[i,4])
        @views errVar[i,5]  = var(evy[i,:]; mean = errMean[i,5])
        @views errVar[i,6]  = var(evz[i,:]; mean = errMean[i,6])
        @views errVar[i,7]  = var(em[i,:]; mean = errMean[i,7])
    end

    return (ts, errMean, errVar)
end