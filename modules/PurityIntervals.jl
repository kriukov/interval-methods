module PurityIntervals
export PurityInterval

using IntervalArithmetic

type PurityInterval
    interval::Interval
    flag::Int
end

PurityInterval(x) = PurityInterval(x, 1)

import Base.+, Base.-, Base.*, Base./
for op in (:+, :-, :*, :/)
	@eval begin
    	$op(x::PurityInterval, y::PurityInterval) = PurityInterval($op(x.interval, y.interval), min(x.flag, y.flag))
	    $op(x::PurityInterval, c) = $op(x, PurityInterval(Interval(c)))
	    $op(c, x::PurityInterval) = $op(PurityInterval(Interval(c)), x)
	end
end

import Base.sqrt
sqrt(x::PurityInterval) = PurityInterval(sqrt_d(x.interval), domaincheck(sqrt1, x.interval))
import Base.asin
asin(x::PurityInterval) = PurityInterval(arcsin_d(x.interval), domaincheck(arcsin, x.interval))

# end of module
end
