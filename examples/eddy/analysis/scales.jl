function τ(B, d, ν=1e-6)
	return (1+2*B)*d^2 / (36*ν)
end

function w0(B, d, g=9.8, ν=1e-6)
	return (1-B)*d^2*g / (18*ν)
end

function t_shear(L, uz, B, d, g=9.8, ν=1e-6)
	return sqrt(2*L / (w0(B, d, g, ν)*uz))
end

function t_stoch(L, κ)
	return L^2 / κ
end

function t_attr(L, r, u, B, d, g=9.8, ν=1e-6)
	return L*g*r / (w0(B, d, g, ν)*u^2)
end

function t_rise(H, B, d, g=9.8, ν=1e-6)
	return H / w0(B, d, g, ν)
end
