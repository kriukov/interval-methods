## Automatic differentiation 

#=
type Dual{T <: Real}
    fun::T
    der::T
end
=#

type Dual
    fun
    der
end


# Basic arithmetic
+(u::Dual, w::Dual) = Dual(u.fun + w.fun, u.der + w.der)
-(u::Dual, w::Dual) = Dual(u.fun - w.fun, u.der - w.der)
*(u::Dual, w::Dual) = Dual(u.fun*w.fun, u.der*w.fun + u.fun*w.der)
function /(u::Dual, w::Dual)
    ff = u.fun/w.fun
    dd = (u.der - ff*w.der)/w.fun
    Dual(ff, dd)
end

# Power; two functions for power because it is less costly
^(u::Dual, a::Integer)= Dual(u.fun^a, a*u.fun^(a-1)*u.der)

function ^(u::Dual, a::Real)
    ff = u.fun^(a-1)    
    Dual(ff*u.fun, a*ff*u.der)
end

# Enable dual-real operations as normal
for op in (:+, :-, :*, :/)
    @eval begin
        $op(u::Dual, c::Real) = $op(u, Dual(c))
        $op(c::Real, u::Dual) = $op(Dual(c), u)
    end
end
+(u::Dual) = u
-(u::Dual) = Dual(-u.fun, -u.der)


#= Trying to add intervals - unnecessary???
for op in (:+, :-, :*, :/)
    @eval begin
        $op(u::Dual, c::Interval) = $op(u, Dual(c))
        $op(c::Interval, u::Dual) = $op(Dual(c), u)
    end
end
=#

# Dual(constant)
Dual{T<:Real}(c::T) = Dual{T}(c, zero(T))

# Elementary functions and their derivatives
for (ff, dd) in ((:(Base.exp), :(Base.exp)), (:(Base.log), :(x->1/x)), (:(Base.sin), :(Base.cos)), (:(Base.cos), :(x->-sin(x))), (:(Base.tan), :(x->(sec(x))^2)))
    @eval begin
        function $ff(u::Dual)
            Dual($ff(u.fun), u.der*$dd(u.fun))
        end
    end
end

# Differentiating function f at point a

#=
function differentiate(f, a)
	deriv(a) = f(Dual(a, one(a))).der
	if typeof(a) != Interval
		return deriv(a)
	else return Interval(min(deriv(a.lo), deriv(a.hi)), max(deriv(a.lo), deriv(a.hi)))
	end
end
=#

differentiate(f, a) = f(Dual(a, one(a))).der


#diff(f, a) = 

#diff(f, a::Interval) = Interval(diff(f, a.lo), diff(a.hi))

	
	
	
	
