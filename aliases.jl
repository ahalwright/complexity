export Population 

using Distributions
using Random
using QuadGK
import Random.seed!

const DIST_TYPE = Dict
const IDIST_TYPE = Dict{Int64,Float64}
const SDIST_TYPE = Dict{String,Float64}
const Population = Vector
const IPopulation = Array{Int64,1}
const SPopulation = Array{String,1}
const PopList = Array{Population,1}

