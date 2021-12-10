module OptimalEstimationProject

using DataFrames
using IndirectTrajOpt
using StaticArrays
using LinearAlgebra
using SparseArrays
using DifferentialEquations
using Roots
using DataInterpolations
using CSV
using Tables
using DrWatson
using ProgressMeter
using AstroTime; AstroTime.update()
using EarthOrientation
using MATLAB

# Time Utilities
include("TIME/gps2MJD.jl")

# GPS Simulation and Dependancies
include("GPS/readSP3.jl")
include("GPS/readERP.jl")
include("GPS/computePrecessionNutation.jl")
include("GPS/rotateData.jl")
include("GPS/genInterpolants.jl")
include("GPS/gpsSimulator.jl")
include("GPS/transmissionTimeFunction.jl")

# IMU Simulation
include("IMU/IMUSim.jl")

# True Trajectory Dependancies
include("TRUTH/SpacecraftSim.jl")

# Filters
include("FILTERS/EKF.jl")

# Precompilation
precompile(GPSSim, (Int64, Int64, Int64, Int64))
precompile(gpsMeasurement, (GPSSim,Int64,Int64,Int64,Int64))

export GPSSim
export IMUSim
export SpacecraftSim
export EKF
export runFilter!
export plotGPS

end
