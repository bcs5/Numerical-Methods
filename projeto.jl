using SymEngine

function euler(f, y, t, h)
	y = y + h*f(t, y)
	return y
end

function euler_inverso(f, y, t, h)
	y1 = euler(f, y, t, h)
	y = y + h*f(t+h, y1)
	return y
end

function euler_aprimorado(f, y, t, h)
	y1 = euler(f, y, t, h)
	y = y + ((f(t, y) + f(t+h, y1))*h)/2
end

function runge_kutta(f, y, t, h)
	k1 = f(t, y)
	k2 = f(t+0.5*h, y+0.5*h*k1)
	k3 = f(t+0.5*h, y+0.5*h*k2)
	k4 = f(t+h, y+h*k3)
	return y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
end

function singleStep(f, h, n, t, y, g, method)
	m = getfield(Main, Symbol(method))
	for i = 0:n
		push!(g, (t, y))
		y = m(f, y, t, h)
		t = t + h
	end
end

const ABcoef = [
	[1, 0, 0, 0, 0, 0, 0, 0],
	[3/2, -1/2, 0, 0, 0, 0, 0, 0],
	[23/12, -4/3, 5/12, 0, 0, 0, 0, 0],
	[55/24, -59/24, 37/24, -3/8, 0, 0, 0, 0],
	[1901/720, -1387/360, 109/30, -637/360, 251/720, 0, 0, 0],
	[4277/1440, -2641/480, 4991/720, -3649/720, 959/480, -95/288, 0, 0],
	[198721/60480, -18637/2520, 235183/20160, -10754/945, 135713/20160, -5603/2520, 19087/60480, 0],
	[16083/4480, -1152169/120960, 242653/13440, -296053/13440, 2102243/120960, -115747/13440, 32863/13440, -5257/17280]
]

const AMcoef = [
	[1, 0, 0, 0, 0, 0, 0, 0],
	[1/2, 1/2, 0, 0, 0, 0, 0, 0],
	[5/12, 2/3, 1/12, 0, 0, 0, 0, 0],
	[3/8, 19/24, -5/24, 1/24, 0, 0, 0, 0],
	[251/720, 323/360, -11/30, 53/360, -19/720, 0, 0, 0],
	[95/288, 1427/1440, -133/240, 241/720, -173/1440, 3/160, 0, 0],
	[19087/60480, 2713/2520, -15487/20160, 586/945, -6737/20160, 263/2520, -863/60480, 0],
	[5257/17280, 139849/120960, -4511/4480, 123133/120960, -88547/120960, 1537/4480, -11351/120960, 275/24192]
]

const DIcoef = [
	[1, 1, 0, 0, 0, 0, 0, 0],
	[4/3, -1/3, 2/3, 0, 0, 0, 0, 0],
	[18/11, -9/11, 2/11, 6/11, 0, 0, 0, 0], 
	[48/25, -36/25, 16/25, -3/25, 12/25, 0, 0, 0],
	[300/137, -300/137, 200/137, -75/137,  12/137, 60/137, 0, 0],
	[360/147, -450/147, 400/147, -225/147, 72/147, -10/147, 60/147, 0]
]

function adamsBF(f, h, g, idx, ordem) 		# ordem N, N pontos
	res = 0
	for i = 1:ordem
		res += ABcoef[ordem][i]*f(g[idx-i][1], g[idx-i][2])
	end
	return g[idx-1][2]+h*res
end


function adamsMU(f, h, g, idx, ordem) 		# ordem N, N-1 pontos, previsao bor adamsBf ordem N-1
	res = AMcoef[ordem][1]*f(g[end][1]+h, adamsBF(f, h, g, idx, ordem-1))
	for i = 1:(ordem-1)
		res += AMcoef[ordem][i+1]*f(g[idx-i][1], g[idx-i][2])
	end
	return g[idx-1][2]+h*res
end

function difInversa(f, h, g, idx, ordem) 	# ordem N, N-1 pontos, previsao bor adamsBf ordem N-1
	res = DIcoef[ordem-1][ordem]*h*f(g[end][1]+h, adamsBF(f, h, g, idx, ordem-1))
	for i = 1:(ordem-1)
		res += DIcoef[ordem-1][i]*g[idx-i][2]
	end
	return res
end

function adamBashExec(f, h, n, t, ordem, g)		# ordem N, N pontos
	for i = 0:n
		if (i+1 > ordem)
			y = adamsBF(f, h, g, i+1, ordem)
			push!(g, (t, y))
		end
		t = t+h
	end
end

function adamMultonExec(f, h, n, t, ordem, g)	# ordem N, N-1 pontos
	for i = 0:n
		if (i+1 >= ordem)
			y = adamsMU(f, h, g, i+1, ordem)
			push!(g, (t, y))
		end
		t = t+h
	end
end

function difInversaExec(f, h, n, t, ordem, g)	# ordem N, N-1 pontos
	for i = 0:n
		if (i+1 >= ordem)
			y = difInversa(f, h, g, i+1, ordem)
			push!(g, (t, y))
		end
		t = t+h
	end
end

function run(input)
	h = 0.1
	n = 0
	t1 = 0.0
	y1 = 0.0
	ordem = 0
	funcao = ""
	t, y = symbols("t y")
	f = ""

	arr = split(input)
	s = arr[1]
	g = []
	
	if (s[1] != 'a' && s[1] != 'f')
		y1 = parse(Float64, arr[2])
		t1 = parse(Float64, arr[3])
		h = parse(Float64, arr[4])
		n =  parse(Int, arr[5])
		funcao = Meta.parse(arr[6])
		f = lambdify(funcao, (t, y))
		singleStep(f, h, n, t1, y1, g, s)
	elseif (s[1] == 'a')
		if (s[6] == 'b')
			if (s == "adam_bashforth")
				ordem = parse(Int, arr[end])
				t1 = parse(Float64, arr[ordem+2])
				h = parse(Float64, arr[ordem+3])
				n =  parse(Int, arr[ordem+4])
				funcao = Meta.parse(arr[ordem+5])
				f = lambdify(funcao, (t, y))
				for i = 2:ordem+1
					push!(g, (t1+h*(i-2), parse(Float64, arr[i])))
				end
			else
				ordem = parse(Int, arr[end])
				y1 = parse(Float64, arr[2])
				t1 = parse(Float64, arr[3])
				h = parse(Float64, arr[4])
				n =  parse(Int, arr[5])
				funcao = Meta.parse(arr[6])
				f = lambdify(funcao, (t, y))
				singleStep(f, h, ordem-1, t1, y1, g, s[19:end])
			end
			adamBashExec(f, h, n, t1, ordem, g)
		else
			if (s == "adam_multon")
				ordem = parse(Int, arr[end])
				t1 = parse(Float64, arr[ordem+1])
				h = parse(Float64, arr[ordem+2])
				n =  parse(Int, arr[ordem+3])
				funcao = Meta.parse(arr[ordem+4])
				f = lambdify(funcao, (t, y))
				for i = 2:ordem
					push!(g, (t1+h*(i-2), parse(Float64, arr[i])))
				end
			else
				ordem = parse(Int, arr[end])
				y1 = parse(Float64, arr[2])
				t1 = parse(Float64, arr[3])
				h = parse(Float64, arr[4])
				n =  parse(Int, arr[5])
				funcao = Meta.parse(arr[6])
				f = lambdify(funcao, (t, y))
				singleStep(f, h, ordem-2, t1, y1, g, s[16:end])
			end
			adamMultonExec(f, h, n, t1, ordem, g)
		end
	elseif (s[1] == 'f')
		if (s == "formula_inversa")
			ordem = parse(Int, arr[end])
			t1 = parse(Float64, arr[ordem+1])
			h = parse(Float64, arr[ordem+2])
			n =  parse(Int, arr[ordem+3])
			funcao = Meta.parse(arr[ordem+4])
			f = lambdify(funcao, (t, y))
			for i = 2:ordem
				push!(g, (t1+h*(i-2), parse(Float64, arr[i])))
			end
		else
			ordem = parse(Int, arr[end])
			y1 = parse(Float64, arr[2])
			t1 = parse(Float64, arr[3])
			h = parse(Float64, arr[4])
			n =  parse(Int, arr[5])
			funcao = Meta.parse(arr[6])
			f = lambdify(funcao, (t, y))
			singleStep(f, h, ordem-2, t1, y1, g, s[20:end])
		end
		difInversaExec(f, h, n, t1, ordem, g)
	end
	println(s)
	println("y( ", t1, " ) = ", y1)
	println("h = ", h)
	for i = 0:n
		println(g[i+1][2])
	end
	
end

fin = open("entrada.txt")
for s in eachline(fin)
	run(s)
	println()
end









