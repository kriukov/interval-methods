
function frac(alpha, epsilon)
	
	irr0 = alpha

	arr_a = Integer[]

	n = 0
	while n <= 38 # change!
		a = int(floor(alpha))
		print("$a, ")
		push!(arr_a, a)
		alpha = 1/(alpha - a)
		n += 1
	end


	h = Integer[]
	push!(h, 0, 1)
	for i = 3:40 # change!
		push!(h, arr_a[i-2]*h[i-1] + h[i-2])
	end
	shift!(h)
	shift!(h)


	k = Integer[]
	push!(k, 1, 0)
	for i = 3:40 # change!
		push!(k, arr_a[i-2]*k[i-1] + k[i-2])
	end
	shift!(k)
	shift!(k)


	function first_solution(alpha, epsilon)
		n = 1
		while abs(alpha*k[n] - h[n]) < epsilon
			n +=1	
		end
		(k[n], h[n])
	end

	first_solution(irr0, epsilon)
end

function eff(m, b, epsilon)
	bb = b
	kn = 0
	while bb > epsilon && 1 - bb > epsilon
		if bb < 0.5
			(q, p) = frac(m, 2bb)
		else
			(q, p) = frac(m, 2*(1 - bb))
		end
		bb = m*q + b - int(m*q + b)
		kn = kn + q		
	end
	q = kn
	p = int(q*m) + 1
	return (q, p)	
end




