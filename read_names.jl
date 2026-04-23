const NAMES_DIST_TYPE = Dict{String,Float64}

# boys = CSV.read("../baby_names/boys_names.csv")

function names_to_dist( names, counts )
  N = length(names)
  num_counts = map( x->parse(Int64,replace(x,","=>"")),counts)
  sum_num_counts = sum(num_counts)
  d = NAMES_DIST_TYPE()
  for i = 1:N
    d[names[i]] = num_counts[i]/sum_num_counts
  end
  d
end

function dist_check( d::NAMES_DIST_TYPE )
  sum = 0.0
  for n in keys(d)
    if d[n] < 0.0
      error(" value of distribution ", d, " is negative ")
    end
    sum += d[n]
  end
  if !isapprox(sum,1.0)
      error(" sum of distribution ", d, " is not 1.0")
  end
end

function entropy( d::NAMES_DIST_TYPE )
  result = 0.0
  for n in keys(d)
    result += -d[n]*log2(d[n])
  end
  result
end

function combine_dists( dist1::NAMES_DIST_TYPE, dist2::NAMES_DIST_TYPE )
  d = NAMES_DIST_TYPE()
  for n in keys(dist1)
    d[n] = 0.5*dist1[n]
  end
  for n in keys(dist2)
    val = get(d,n,0.0)
    d[n] = val + 0.5*dist2[n]
  end
  println("A len(d): ",length(d),"  sum(d): ",sum_dist(d))
  d
end

function mutual_information( D::Vector{NAMES_DIST_TYPE} )
  combined_dist = reduce( combine_dists, D[i] for i = 1:length(D) )
  entropy( combined_dist ) - mean( map( entropy, D ) )
end

function sum_dist( dist::NAMES_DIST_TYPE )
  sum = 0.0
  for k in keys(dist)
    sum += dist[k]
  end
  sum
end
