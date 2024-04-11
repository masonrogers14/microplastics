#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Mon May 16 2022

@author Mason Rogers

crit_var.jl computes the (unstable) limiting variance of the problem for r≪a
=#

using Printf

#read parameters
io = open("../input/kv_param.py", "r")
param_str_list = readlines(io)
close(io)
for param_str in param_str_list
    param_expr = Meta.parse(param_str)
    if !(typeof(param_expr) == Nothing)
        if param_expr.head == :(=)
            eval(param_expr)
        end
    end
end

ϵ = ((1+2*B)*d^2*Us)/(36*ν*Ls)
c = Ls*ϵ/Us * 2*(1-B)/(1+2*B) * (Γ/(2*π))^2 * a^-4
Σc = 2*κ/c
σc = sqrt(Σc)
s = 4*κ
τ = 1/(2*c)

println("Relevant Parameters")
println("-------------------")
println("Units of L²:")
@printf "%-13s %.5f\n" "Σc (R² ≪ a²)" Σc
@printf "%-13s %.5f\n" "a²" a^2
@printf "%-13s %.5f\n" "Δx²" dx^2
println("Units of L²/T:")
@printf "%-13s %.5f\n" "s (R² ≫ a²)" s
@printf "%-13s %.5f\n" "s (R² = 0)" (4*κ/c - .001)/τ
println("Units of T:")
@printf "%-13s %.5f\n" "τ" τ

