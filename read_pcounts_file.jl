
function read_pcount_file( fname )
  open( fname ) do f
    map( x->parse(Int64,x),readlines(f))
  end
end
    
