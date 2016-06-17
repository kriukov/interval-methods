module PurityIntervalsVN
export PurityInterval, purity, purify, unpurify

using ValidatedNumerics

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
	if in(x.interval, domain)
		return PurityInterval(sqrt(x.interval), 1)
	else
		y = intersect(x.interval, domain)
		if y != ∅
			return PurityInterval(sqrt(y), 0)
		else
			return PurityInterval(∅, -1)
		end
	end
end

import Base.asin
function asin(x::PurityInterval)
	domain = Interval(-1, 1)
	if in(x.interval, domain)
		return PurityInterval(asin(x.interval), 1)
	else
		y = intersect(x.interval, domain)
		if y != ∅
			return PurityInterval(asin(y), 0)
		else
			return PurityInterval(∅, -1)
		end
	end
end

purity(f, x::Interval) = f(PurityInterval(x)).flag

function purity(f, x::IntervalBox)
	p = f([PurityInterval(x[1]), PurityInterval(x[2])])
	min(p[1].flag, p[2].flag)
end

#=
function purity(f, x::IntUnion)
    p = Int[]
    l = length(x.union)
    for i = 1:l
        push!(p, f(PurityInterval(x.union[i])).flag)
    end
	minimum(p)
end
=#



# Convert a 2D rectangle into one with purity 1 and extract rectangle from a PurityInterval
purify(x::IntervalBox) = IntervalBox(PurityInterval(x[1], 1), PurityInterval(x[2], 1))
unpurify(x::Array{PurityInterval, 1}) = IntervalBox(x[1].interval, x[2].interval)


# end of module
end
