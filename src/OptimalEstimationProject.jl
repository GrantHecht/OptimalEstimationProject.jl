module OptimalEstimationProject

using DataFrames
using StaticArrays
using LinearAlgebra
using DrWatson
using ProgressMeter

# Time Utilities
include("TIME/gps2MJD.jl")

# GPS Simulation and Dependancies
include("GPS/readSP3.jl")
include("GPS/readERP.jl")
include("GPS/gpsSimulator.jl")

export readSP3
export readSP3s
export gpsSim

end
