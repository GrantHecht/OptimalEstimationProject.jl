
# Very simple IMU accelerometer simulation
mutable struct IMUSim
    σa::Float64
end

function ComputeMeasurement(imu::IMUSim, acc; type = :true) 
    if type == :true
        meas = @SVector [acc[1] + imu.σa*randn(), 
                         acc[2] + imu.σa*randn(),
                         acc[3] + imu.σa*randn()]
    else
        meas = @SVector [acc[1], acc[2], acc[3]]
    end

    return meas
end