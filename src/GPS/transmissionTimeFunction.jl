# To determine time of signal transmission, this function must be solved

# Returns the difference between the distance light travels in signal reception time (t) - signal transmission time (tt)
# and the geometric distance between the two satelites
function transmissionTimeFunction(tt, t, r, X, Y, Z)
    # Speed of light
    c = 299792.458 # [km / s]

    rsat = @SVector [X(tt), Y(tt), Z(tt)]
    func = c*(t - tt) - norm(r .- rsat) 

    return func
end

# Determin transmission time through itteration
#  https://gssc.esa.int/navipedia/index.php/Emission_Time_Computation
function transmissionTimeFunction(t,r,X,Y,Z)
    # Speed of light
    c = 299792.458 # [km / s]

    # Convergence criteria
    iMax = 100
    tol  = 1e-12

    rsat    = @SVector [X(t), Y(t), Z(t)]
    stop    = false
    iter    = 0
    tt      = 0.0
    while stop == false && iter < iMax
        iter += 1

        Δt = sqrt((r[1] - rsat[1])^2 + (r[2] - rsat[2])^2 + (r[3] - rsat[3])^2) / c
        tt = t - Δt

        rdiff = sqrt((rsat[1] - X(tt))^2 + (rsat[2] - Y(tt))^2 + (rsat[3] - Z(tt))^2)
        if rdiff < tol
            stop = true
        else
            rsat = @SVector [X(tt), Y(tt), Z(tt)]
        end
    end

    if iter == iMax
        throw(ErrorException("Transmission time function did not converge."))
    end

    return tt
end
