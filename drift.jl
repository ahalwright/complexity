# 
include("types.jl")
const Population = Array{Int64,1}

# Initial population size
N = 5
p = DIST_TYPE( i=> 1.0/N for i = 1:N )
pop = [1,2,3,4]

# pop can include repeated elements.
# If i is repeated n times in pop, result[i] is Float64(n)/N.
function pop_to_dist( pop::Population )
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

function rel_entropy( q::Population, p::Population )
  relative_entropy( pop_to_dist(q), pop_to_dist(p) )
end

function drift( N::Int64, p::DIST_TYPE )
  q = deepcopy(p)
  for i = 1:N
    q = deepcopy(p)
    n = length(keys(q))
    k = rand(1:n)
    while q[k] == 0.0
      k = rand(1:n)
    end  
    for x in keys(q)
      if x != k
        q[x] += 1.0/n/(n-1)
        println("ne x: ",x,"  q[x]: ",q[x])
      else
        if 1.0/n == 0.
          delete!(q,x)
        else
          q[x] -= 1.0/n
        end
        println("eq x: ",x,"  q[x]: ",q[x])
      end
    end
    println("i: ",i,"  n: ",n,"  k:",k,"   q: ",q)
    dist_check(q)
    r = relative_entropy(q,p)
    println("i: ",i,"  rel entropy(q,p): ",r)
    p = q
  end
  q
end
    
