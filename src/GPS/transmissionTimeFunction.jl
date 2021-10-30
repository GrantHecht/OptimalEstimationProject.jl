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