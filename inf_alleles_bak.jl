include("types.jl")
#const Population = Array{Int64,1}
#const PopList = Array{Population,1}
using QuadGK
using Printf

@doc """ function simple_poplist( N::Int64, mu_per_locus::Float64, ngens::Int64, burn_in::Int64 )
No-frills implmentation of Wright-Fisher infinite alleles model.
N is the population size
N_mu is the per-population  mutation rate,  the per individual mutation rate is N_mu/N.
If combine==true, returns a single combined population of the populations from all ngens generations after burn in.
If combine==false, returns a list of ngens populations.
"""

function simple_poplist( N::Int64, N_mu::Float64, ngens::Int64; burn_in::Float64=0.0, combine::Bool=false)
  int_burn_in = burn_in > 0.0 ? Int(round(burn_in*N/N_mu)) : 0
  mu = N_mu/N
  println("int_burn_in: ",int_burn_in,"  mu: ",mu)
  base_pop = collect(1:N)
  poplist= Population[ collect(1:N) ]
  if combine
    pop_result = Population()
  end
  new_id = N+1
  g = 1
  print("g: ",g,"  pop: ",base_pop)
  @printf(" re: %4.2f  ent: %4.2f\n",relative_entropy(base_pop,base_pop),entropy(base_pop))
  for g = 2:(ngens+int_burn_in)
    #println("generation: ",g)
    pctr = pop_counter( poplist[g-1] )
    result = zeros(Int64,N);
    new_base_pop = Int64[]
    for i = 1:N
      if rand() < mu
        result[i] = new_id
        push!(new_base_pop,new_id)
        new_id += 1
      else  # Choose a random element of the previous population
        result[i] = poplist[g-1][rand(1:N)]
        push!(new_base_pop,result[i])
      end
    end
    #=
    base_set = Set{Int64}(base_pop)
    base_list = [s for s in base_set]
    for i = (length(new_base_pop)+1):N
      b = pop!(base_list) in new_base_pop && i < N end
      while b in new_base_pop
        b = pop!(base_list)
      end
      push!(new_base_pop,b)
    end
    =#
      
    print("g: ",g,"  pop: ",result)
    @printf(" re: %4.2f  ent: %4.2f\n",relative_entropy(result,base_pop),entropy(result))
    Base.push!(poplist,result)
    if combine && g >= int_burn_in+1
      pop_result = vcat( pop_result, result )
    end
  end
  if combine
    return [pop_result]
  else
    return poplist[int_burn_in+1:end]
  end
end 

@doc """ function pop_counts64( pop::Population )
Returns the sorted frequencies of the alleles of Population pop.
Example:  If pop = [5, 7, 9, 5, 4, 5, 7], then the returned list is [3, 2, 1, 1]
   because there are 3 5's, 2 7's, 1 9, and 1 4.  So the sum of the returned
   list is the length of the population.
"""
function pop_counts64( pop::Population )
  c = Dict{Int64,Int64}()
  for x in pop
    c[x] = get( c, x, 0 ) + 1
  end
  map( x->c[x], sort( unique(pop), by=x->c[x], rev=true ) )
end

function pop_counter( pop::Population )
  c = Dict{Int64,Int64}()
  for x in pop
    c[x] = get( c, x, 0 ) + 1
  end
  c
end

#=
integrand(x,theta) = theta/x*(1-x)^(theta-1)

function expected_richness( N::Int64=5, N_mu::Float64=1.0 )
  theta = 2.0*N_mu
  println("N: ",N,"  N_mu: ",N_mu,"  integrand(1/N,theta): ",integrand(1/2/N,theta))
  result, err = quadgk( x->integrand(x,theta), 1.0/N, 1.0 )
  result += theta
  println("add_expected richness result: ",result,"  err: ",err)
  result
end

# Return the pair of expected and observed richness
function exp_observed_richness(N::Int64=5, N_mu::Float64=1.0, ngens::Int64=5, burn_in::Float64=2.0 )
  exp_richness = expected_richness( N, N_mu ) 
  pop_result_list = simple_poplist(N,N_mu,ngens)
  #println("pop_result: ",pop_result_list)
  return exp_richness, mean(map(x->length(pop_counter(x)),pop_result_list))
end
=#
