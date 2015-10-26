julia> ex = :(f(x) = [sqrt(x[1]), 2*x[2]])
:(f(x) = [sqrt(x[1]),2 * x[2]])

julia> body = ex.args[2]
:([sqrt(x[1]),2 * x[2]])

julia> eval(body.args[1])
ERROR: x not defined

julia> body.args
2-element Array{Any,1}:
 :(sqrt(x[1]))
 :(2 * x[2])  

julia> typeof(ans
       )
Array{Any,1}

julia> f2(x) = eval(body.args[2])
f2 (generic function with 1 method)

julia> f2(3)
ERROR: x not defined
 in f2 at none:1

julia> f2([2,3])
ERROR: x not defined
 in f2 at none:1

julia> x = [2,3]
2-element Array{Int32,1}:
 2
 3

julia> f2(x)
6

