
function first_collision(irr, delta)
	
	irr0 = irr

	arr_a = Integer[]

	n = 0
	while n <= 38
		a = floor(irr)
		#println(a)
		@show push!(arr_a, a)
		@show irr = 1/(irr - a)
		n += 1
	end


	h = Integer[]
	push!(h, 0, 1)
	for i = 3:40
		push!(h, arr_a[i-2]*h[i-1] + h[i-2])
	end
	shift!(h)
	shift!(h)


	k = Integer[]
	push!(k, 1, 0)
	for i = 3:40
		push!(k, arr_a[i-2]*k[i-1] + k[i-2])
	end
	shift!(k)
	shift!(k)

	#for i = 1:length(h)
	#	println("$(abs(irr0 - h[i]/k[i])), $(1/(sqrt(5)*(k[i])^2))")
	#	@show abs(irr0 - h[i]/k[i]) <= 1/(sqrt(5)*(k[i])^2)
	#end


	function first_solution(alpha, delta)
		n = 1
		while k[n] < 1/(sqrt(1 + alpha^2)*sqrt(5)*delta)
			n +=1	
		end
		[k[n], h[n]]
	end

	first_solution(irr0, delta)
end
