rect = [Interval(-0.5-0.0001, -0.5+0.0001), Interval(prec(pi)/6-0.0001, prec(pi)/6+0.0001)]

f(x) = path_general(x, c, [1,2,3,1])
g(x) = f(x) - x
krawczyk2d(g, rect, 64)

#=
Unique zero in [IntervalArithmetic.Interval(-0.5000000000001089,-0.49999999999999994),IntervalArithmetic.Interval(0.5235987755982217,0.523598775598299)]
Unique zero in [IntervalArithmetic.Interval(-0.5000000000000817,-0.49999999999999994),IntervalArithmetic.Interval(0.523598775598299,0.5235987755983784)]
Unique zero in [IntervalArithmetic.Interval(-0.49999999999999994,-0.4999999999999189),IntervalArithmetic.Interval(0.5235987755982198,0.523598775598299)]
Unique zero in [IntervalArithmetic.Interval(-0.49999999999999994,-0.4999999999999415),IntervalArithmetic.Interval(0.523598775598299,0.523598775598378)]
4-element Array{Array{IntervalArithmetic.Interval,1},1}:
 [IntervalArithmetic.Interval(-0.5000000000001089,-0.49999999999999994),IntervalArithmetic.Interval(0.5235987755982217,0.523598775598299)]
 [IntervalArithmetic.Interval(-0.5000000000000817,-0.49999999999999994),IntervalArithmetic.Interval(0.523598775598299,0.5235987755983784)]
 [IntervalArithmetic.Interval(-0.49999999999999994,-0.4999999999999189),IntervalArithmetic.Interval(0.5235987755982198,0.523598775598299)]
 [IntervalArithmetic.Interval(-0.49999999999999994,-0.4999999999999415),IntervalArithmetic.Interval(0.523598775598299,0.523598775598378)] 
=#

# The first solution is true, other 3 false! Check by substitution:

#=
julia> g(sol_kraw1[1])
2-element Array{IntervalArithmetic.Interval,1}:
 IntervalArithmetic.Interval(-8.872230727874353e-11,2.423063416578941e-10)  
 IntervalArithmetic.Interval(-1.4215018051544348e-10,2.9175617477505966e-10)

julia> g(sol_kraw1[2])
2-element Array{IntervalArithmetic.Interval,1}:
 IntervalArithmetic.Interval(-1.4020024030614309e-10,1.4034068351875817e-10)
 IntervalArithmetic.Interval(-1.8386203670672785e-10,1.839027818917316e-10) 

julia> g(sol_kraw1[3])
2-element Array{IntervalArithmetic.Interval,1}:
 IntervalArithmetic.Interval(-1.39411204802542e-10,1.3941403587125478e-10)
 IntervalArithmetic.Interval(-1.82645787383251e-10,1.8287005243422527e-10)

julia> g(sol_kraw1[4])
2-element Array{IntervalArithmetic.Interval,1}:
 IntervalArithmetic.Interval(-1.7383139372384449e-10,5.99419958113856e-11)
 IntervalArithmetic.Interval(-2.08681738556038e-10,9.780898313493935e-11)
=#
 
# krawczyk is not to blame - probably the sin function is! Look:

using KrawczykMethod
h(x) = sin(x) - x
krawczyk(h, Interval(-0.5, 0.5), 64)
# Solution: 0. Spews out lots of false solutions.

h(x) = sin(x) - 1
krawczyk(h, Interval(1, 2), 64)
# Solution: pi/2. Finds nothing.

