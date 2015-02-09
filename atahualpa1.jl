function fractions(alfa, ap)
	p=0;q=0
	#real*8 ap, alfa, mm(10000), prueba, pap
	#real rand
	#integer(8) i,j,h(10000), k(10000), a(10000), p,q
	
	mm = zeros(10000)
	a = zeros(Int, 10000)
	h = zeros(Int, 10000)
	k = zeros(Int, 10000)
	
	#write(*,*) alfa
	h[1]=0
	h[2]=1
	k[1]=1
	k[2]=0
	mm[1]=alfa
	a[1]=ifloor(mm[1])
	prueba = 10000.
	i=1
	while abs(k[i+1]*mm[1]-h[i+1])>ap

		#write(*,*) k(i+1), ap, mm(i), i
		mm[i+1]=1/(mm[i]-a[i])
		a[i+1]=ifloor(mm[i+1])
		h[i+2]=a[i]*h[i+1]+h[i]
		k[i+2]=a[i]*k[i+1]+k[i]
		p=h[i+2]
		q=k[i+2]

		i += 1
		
		if prueba<abs(h[i+2]/(k[i+2])-alfa)
			#write(*,*) ":-o", h(i+2)/(k(i+2)), h(i+2), (k(i+2)) 
			continue #goto 220
		end
		#prueba=abs(h[i+2]/(k[i+2])-m)
	end
	#220 continue
	p,q
end 




b = sqrt(rand())
m = sqrt(rand())


theta = rand()*pi/2
v = [cos(theta), sin(theta)] # Already normalized
m = v[1]/v[2] # Why? It should be v2/v1
delta = 100
i = 0

while delta>0.0000001*v[1]
	i += 1

	y = m*i + b
	hn = ifloor(y)
	kn = ifloor(y)+1

	delta1 = y - hn
	delta2 = kn - y

	if delta1 < delta
		delta = delta1
		p = hn
	end

	if delta2 < delta
		delta = delta2
		p = kn
	end
end


#write(*,*) sqrt(i**2.+p**2.), i, p
while b > 0.001

	if b<0.5
	ap1=b

	ap2=(1-b)
	alfa=m
	p3,i3 = fractions(alfa,ap1)

	pru=ifloor(i3*alfa)+1
	end
	#write(*,*) pru, i3, p3
	if pru != p3
		i3=-1 
	end
	p4,i4 = fractions(alfa,ap2)
	
	
	pru=ifloor(i4*alfa)

	#write(*,*) pru, i4, p4
	if pru != p4 
		i4=-1
	else

	ap1=(1-b)
	ap2=(2-2b)

	p3,i3 = fractions(alfa,ap1)

	p4,i4 = fractions(alfa,ap2)
	end
	
	
	if i3>0 && i4>0
	b=max(alfa*i3+b-p3, alfa*i4+b-p4)
	elseif i3<0 && i4>0
	b=alfa*i4+b-p4
	println("falla i3")
	elseif i3>0 && i4<0
	b=alfa*i3+b-p3
	println("falla i4")
	elseif i3<0 && i4<0
	println("falla i3 and i4")
	end
	
	println(b)
	
end
