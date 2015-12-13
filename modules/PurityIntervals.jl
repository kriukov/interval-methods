module PurityIntervals
export PurityInterval, purity

using IntervalArithmetic

type PurityInterval
    interval#::Interval
    flag::Int
end

# 1 = full domain ("clean");
# 0 = partially defined ("unclean")
# -1 = undefined ("dirty")

function PurityInterval(x)
    if !isnan(x)
        return PurityInterval(x, 1)
    else
        return PurityInterval(x, -1)
    end
end


#PurityInterval(x::IntUnion) = IntUnion(map(PurityInterval, x.union))

import Base.+, Base.-, Base.*, Base./
for op in (:+, :-, :*, :/)
	@eval begin
    	$op(x::PurityInterval, y::PurityInterval) = PurityInterval($op(x.interval, y.interval), min(x.flag, y.flag))
	    $op(x::PurityInterval, c) = $op(x, PurityInterval(Interval(c)))
	    $op(c, x::PurityInterval) = $op(PurityInterval(Interval(c)), x)
	end
end

-(x::PurityInterval) = PurityInterval(-x.interval, x.flag)

# Embedding a PurityInterval into a PurityInterval should be idempotent
PurityInterval(x::PurityInterval, y::Int) = PurityInterval(x.interval, min(y, x.flag))

import Base.sin
sin(x::PurityInterval) = PurityInterval(sin(x.interval), x.flag)
import Base.cos
cos(x::PurityInterval) = PurityInterval(cos(x.interval), x.flag)
import Base.mod
mod(x::PurityInterval, y::Real) = PurityInterval(mod(x.interval, y), x.flag)


# On "empty-set" interval [Inf, Inf]: if you try to do [Inf, -Inf], then its subtraction may yield the whole real line!

import Base.sqrt
function sqrt(x::PurityInterval)
	domain = Interval(0, Inf)
	if inside(x.interval, domain)
		return PurityInterval(sqrt(x.interval), 1)
	else
		y = isect(x.interval, domain)
		if y != false
			return PurityInterval(sqrt(y), 0)
		else
			return PurityInterval(Interval(Inf, Inf), -1)
		end
	end
end

import Base.asin
function asin(x::PurityInterval)
	domain = Interval(-1, 1)
	if inside(x.interval, domain)
		return PurityInterval(asin(x.interval), 1)
	else
		y = isect(x.interval, domain)
		if y != false
			return PurityInterval(asin(y), 0)
		else
			return PurityInterval(Interval(Inf, Inf), -1)
		end
	end
end

purity(f, x::Interval) = f(PurityInterval(x)).flag

function purity(f, x::MultiDimInterval)
	p = f([PurityInterval(x[1]), PurityInterval(x[2])])
	min(p[1].flag, p[2].flag)
end

function purity(f, x::Array{MultiDimInterval, 1})
	p = f([PurityInterval(x[1]), PurityInterval(x[2])])
	min(p[1].flag, p[2].flag)
end

#=
function purity(f, x::IntUnion2D)
    if typeof(x[1]) == IntUnion || typeof(x[2]) == Interval
        for i = 1:length(x[1].union)
            return purity(f, [x[1].union[i], x[2]])
        end
    elseif typeof(x[1]) == Interval || typeof(x[2]) == IntUnion
        for i = 1:length(x[2].union)
            return purity(f, [x[1], x[2].union[i]])
        end  
    end
end
=#

# end of module
end
