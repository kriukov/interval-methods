## Automatic differentiation (interval version)

type ad
    u
    up
end

ad(c) = ad(c, Interval(0.0))

function +(x::ad, y::ad)
    ad(x.u + y.u, x.up + y.up)
end

function -(x::ad, y::ad)
    ad(x.u - y.u, x.up - y.up)
end

function *(x::ad, y::ad)
    ad(x.u*y.u, x.u*y.up + y.u*x.up)
end

function /(x::ad, y::ad)
    ad(x.u/y.u, (x.up - (x.u/y.u)*y.up)/y.u)
end

# Multiplying by numbers (intervals)

function *(x::ad, y::Interval)
    ad(x.u*y, x.up*y)
end

function *(x::Interval, y::ad)
    ad(x*y.u, x*y.up)
end

function /(x::ad, y::Interval)
    ad(x.u/y, x.up/y)
end

function /(x::Interval, y::ad)
    ad(x/y.u, x/y.up)
end


function +(x::Interval, y::ad)
    ad(x) + y
end

function +(x::ad, y::Interval)
    x + ad(y)
end

function -(x::Interval, y::ad)
    ad(x) - y
end

function -(x::ad, y::Interval)
    x - ad(y)
end

# Elementary functions

sin(x::ad) = ad(sin(x.u), x.up*cos(x.u))
cos(x::ad) = ad(cos(x.u), -x.up*sin(x.u))
e^(x::ad) = ad(e^x.u, x.up*e^x.u)
exp(x::ad) = ad(exp(x.u), x.up*exp(x.u))
log(x::ad) = ad(log(x.u), x.up/x.u)
(x::ad)^y::Interval = ad(x.u^y, x.up*y*x.u^(y-1))

differentiate(f, a) = f(ad(a, 1.)).up
