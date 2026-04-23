using CSV
using DataFrames
include("../src/aliases.jl")


function dist_check( p::DIST_TYPE )
  sum = 0.0
  for x in keys(p)
    if p[x] < 0.0
      error(" value of distribution ", p, " is negative ")
    end
    sum += p[x]
  end
  if !isapprox(sum,1.0)
      error(" sum of distribution ", p, " is not 1.0")
  end
end

function pop_to_dist( pop::IPopulation )
  N = length(pop)
  result = DIST_TYPE()
  for i in pop
    #pval = get(result,i,0.0)
    result[i] = get(result,i,0.0) + 1.0/N
    #println("i: ",i,"  pval: ",pval,"  result: ",result)
  end
  dist_check(result)
  result
end

function pop_to_dist( pop::SPopulation )
  N = length(pop)
  result = DIST_TYPE()
  for i in pop
    #pval = get(result,i,0.0)
    result[i] = get(result,i,0.0) + 1.0/N
    #println("i: ",i,"  pval: ",pval,"  result: ",result)
  end
  dist_check(result)
  result
end

# There is one column per "allele" of the combined population
# The sum of the entries in the result is 1.0.
# The sum of the entries in row i is length(P[i])/(sum(length(p) for p in P)
function pops_to_tbl( P::Vector{Population} )
  #P = map(x->convert(Vector{Any},x),P)
  m = length(P)   # the number of rows in result
  combined_pop = reduce(vcat, P[i] for i = 1:m )
  N = length(combined_pop)
  allele_list = unique(combined_pop)
  n = length(allele_list)   # The number columns in result
  result = zeros(Float64,m,n)
  for i = 1:length(P)
    dist = Dict( allele_list .=> zeros(Float64,n))
    for p in P[i]
      dist[p] += 1.0/N
    end
    result[i,:] = [ dist[a] for a in allele_list ]
  end
  result
end

function pops_to_tbl( P::Vector{IPopulation} )
  pops_to_tbl( convert(Vector{Population},P) )
end

function pops_to_tbl( P::Vector{SPopulation} )
  pops_to_tbl( convert(Vector{Population},P) )
end

# D is a DataFrame with 2 columns per population.
# Each odd numbered column is a list of Strings which are "allele names".
# The next even numbered column is an integer count for the corresponding allele name.
# Each pair of columns produces a column of the table.
function pop_counts_to_tbl( D::DataFrame, columns::Vector{Int64}=collect(1:(div(size(D)[2],2))))
  #m = div(length(names(D)),2)   # the number of rows in result
  m = length(columns)
  #println("m: ",m)
  # Convert columns of data frame D from CSV columns to columns of Strings or columns of Int64.
  for i = 1:m
    D[!,2*i-1] = convert(Vector{typeof(D[!,1][1])} ,D[!,2*i-1])
    D[!,2*i] = convert(Vector{Int64},D[!,2*i])
  end
  combined_pop = reduce(vcat, D[!,2*i-1] for i in columns)   # combine allele names
  combined_counts = reduce(vcat, D[!,2*i] for i in columns)  # combine counts
  N = length(combined_pop)
  combined_sum = sum(combined_counts)
  allele_list = unique(combined_pop)
  n = length(allele_list)   # The number columns in result
  #println("m: ",m,"  N: ",N,"  n: ",n,"  combined_sum: ",combined_sum)
  result = zeros(Float64,m,n)
  for i = 1:m
    ii = columns[i]
    #println("i: ",i,"  ii: ",ii)
    dist = Dict( allele_list .=> zeros(Float64,n))
    for j = 1:length(D[!,2*ii-1])
      dist[D[!,2*ii-1][j]] += D[!,2*ii][j]/combined_sum
    end
    result[i,:] = [ dist[a] for a in allele_list ]
  end
  result
end

# Converts a row of a joint probability distribution table to a DIST
# The DIST is normalized so that the probabilities sum to 1.0
function table_row_to_dist( tbl::Array{Float64,2}, row_index::Int64 )
  result = DIST_TYPE()
  (m,n) = size(tbl)
  ssum = sum(tbl[row_index,:])
  for j = 1:n
    result[j] = get(result,j,0.0) + tbl[row_index,j]/ssum
  end
  result
end

function row_marginal( tbl::Array{Float64,2} )
  [ sum(tbl[i,:]) for i = 1:size(tbl)[1]]
end

function column_marginal( tbl::Array{Float64,2} )
  [ sum(tbl[:,j]) for j = 1:size(tbl)[2]]
end

function entropy( p::DIST_TYPE; base::Float64=2.0 )
  result = 0.0
  for x in keys(p)
    result += p[x] > 0.0 ? -p[x]*log(base, p[x] ) : 0.0
    #println("x: ",x, "result_inc: ", p[x] * log(base, p[x]),"  result: ",result)
  end
  result
end

function entropy( p::IPopulation; base::Float64=2.0 )
  entropy( pop_to_dist(p), base=base )
end

function entropy( p::SPopulation; base::Float64=2.0 )
  entropy( pop_to_dist(p), base=base )
end

function entropy( p::Vector{Float64}; base::Float64=2.0 )
  entf(x) = x > 0.0 ? -x*log(base,x) : 0.0
  reduce(+,map(entf,p))
end

function entropy( tbl::Array{Float64,2}, row_index::Int64; base::Float64=2.0 )
  entropy( table_row_to_dist( tbl, row_index ), base = base )
end

# relative_entropy()
# Also Kullback Leibler divergence D( q || p )  
# Note that if there are any keys x of p such that q[x] == 0, the result will be NaN
#   indicating that the result is not defined.
# Reversed p and q to agree with Cover & Thomas on 9/16/19
function relative_entroqy( p::DIST_TYPE, q::DIST_TYPE; base::Float64=2.0 )
  result = 0.0
  for x in keys(p)  # Assume that any keys of q that aren't in p have value 0 for p
    try qval = q[x]
    catch
      return NaN
    end
    #pval = get(q,x,0.0)
    pval = get(p,x,0.0)
    result += pval > 0.0 ? pval*(log(base,pval)-log(base,q[x])) : 0.0
    #println("x: ",x,"  pval: ",pval,"  q[x]: ",q[x],"  result: ",result)
  end
  result
end

function relative_entropy( q::IPopulation, p::IPopulation; base::Float64=2.0 )
  relative_entropy( pop_to_dist(q), pop_to_dist(p), base=base )
end

function relative_entropy( q::SPopulation, p::SPopulation; base::Float64=2.0 )
  relative_entropy( pop_to_dist(q), pop_to_dist(p), base=base )
end

function relative_entropy( p::Vector{Float64}, q::Vector{Float64}; base::Float64=2.0 )
  result = 0.0
  for i = 1:length(p)
    if q[i] == 0.0
      return NaN
    end
    result += p[i]*log(base,p[i]/q[i])
  end
  result
end

# Conditional entropy of a table.
# H(Y|X)  where  X  is the row-marginal of the table, and Y is the column marginal of the table
# The rows argument specifies the rows of the table that are used.  Order of rows is not important.
# Verification in the file /home/evotech/information_theory/test/conditional_entropy_test.jl
function conditional_entropy( tbl::Array{Float64,2}, rows::Vector{Int64}=collect(1:size(tbl)[1]); base::Float64=2.0 )
  #println("rows: ",rows)
  row_sums = map(x->sum(tbl[x,:]),rows)
  tbl_sum = sum(row_sums)   # sum of the part of the table that we are using.
  result = 0.0
  for x = 1:length(rows)
    #sum_px = sum(tbl[rows[x],:])
    for y = 1:size(tbl)[2]
      pyx = tbl[rows[x],y]/row_sums[x]
      result -= tbl[rows[x],y]/tbl_sum > 0.0 ? tbl[rows[x],y]/tbl_sum*log(base,pyx) : 0.0
      #println("(x,y): ",(x,y),"  pyx: ",pyx,"  inc: ",tbl[rows[x],y]>0.0 ? -tbl[rows[x],y]*log(2,pyx) : 0.0)
    end
  end
  result
end

function joint_entropy( tbl::Array{Float64,2}; base::Float64=2.0 )
  if ! reduce( &, tbl .>= 0.0 )
    error("all entries of tbl must be non-netative in function joint entropy")
  end
  if !(sum(tbl) ≈ 1.0)
    tbl = tbl/sum(tbl)   # Normalize to sum 1
  end
  (m,n) = size(tbl)
  entf(x) = x != 0 ? -x*log(base,x) : 0.0
  reduce(+, map(entf,tbl))
end


# Mutual information by the "standard" defintion of eq. 2.28 of Cover and Thomas.
#   This is the mutual information of a joint distribution.
function mutual_information( tbl::Array{Float64,2}; base::Float64=2.0 )
  if ! reduce( &, tbl .>= 0.0 )
    error("all entries of tbl must be non-netative in function mutual information")
  end
  if !(sum(tbl) ≈ 1.0)
    tbl = tbl/sum(tbl)   # Normalize to sum 1
  end
  (m,n) = size(tbl)
  mif(i,j) = tbl[i,j] != 0.0 ? tbl[i,j]*log(base,tbl[i,j]/mr1[i]/mr2[j]) : 0.0
  mr1 = [ sum( tbl[i,:]) for i = 1:m ] 
  mr2 = [ sum( tbl[:,j]) for j = 1:n ]
  reduce( +, [ mif(i,j) for i=1:m, j=1:n ] )
end

# tbl is the probability distribution for a total population which is made up of m subpopulations
#    all over the same allele set of size n.
# Row i of tbl corresponds to a probability distribution over subpopulation i.
function sherwin_mutual_information( tbl::Array{Float64,2}; base::Float64=2.0 )
  entf(x) = x != 0 ? -x*log(base,x) : 0.0
  (m,n) = size(tbl)
  # weight[i] is the fraction of the total population corresponding to subpopulation i.
  weights = [ sum(tbl[i,:]) for i = 1:m ]
  s = [ sum( map(entf,tbl[i,:]/weights[i])) for i = 1:m ]
  ent_total = sum(map(entf,[sum(tbl[:,j]) for j = 1:n]))
  ent = ent_total - sum(s[i]*weights[i] for i = 1:m)
  #println(" -- sh ent: ",ent)
  ent
end

function mutual_information( D::DataFrame, columns::Vector{Int64}=collect(1:(div(size(D)[2],2))); base::Float64=2.0 )
  mutual_information(pop_counts_to_tbl( D, columns ), base=base )
end

function mutual_information( P::Vector{IPopulation}; base::Float64=2.0 )
  sherwin_mutual_information( pops_to_tbl( P ), base=base )
end

function mutual_information( P::Vector{SPopulation}; base::Float64=2.0 )
  sherwin_mutual_information( pops_to_tbl( P ), base=base )
end

