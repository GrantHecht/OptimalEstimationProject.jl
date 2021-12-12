module OptimalEstimationProject

using DataFrames
using IndirectTrajOpt
using StaticArrays
using LinearAlgebra
using SparseArrays
using SuiteSparse
using DifferentialEquations
using Statistics
using StatsBase
using Roots
using DataInterpolations
using CSV
using Tables
using DrWatson
using ProgressMeter
using AstroTime; AstroTime.update()
using EarthOrientation
using MATLAB

# Linear Algebra Utilities
#using LoopVectorization # Using BLAS routines rather than pure Julian linear algebra 
#include("UTILS/vecOuterProd.jl")

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
include("FILTERS/UKF.jl")

# Precompilation
precompile(GPSSim, (Int64, Int64, Int64, Int64))
precompile(gpsMeasurement, (GPSSim,Int64,Int64,Int64,Int64))

export GPSSim
export IMUSim
export SpacecraftSim
export EKF
export UKF
export runFilter!
export plotGPS
export plotEKF

end
