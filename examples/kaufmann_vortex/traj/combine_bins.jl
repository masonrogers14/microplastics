#!/usr/bin/env julia
# -*- coding utf-8 -*-
#=
Created on Thu Apr 18 2024

@author Mason Rogers

combine_bins.jl combines single-processor binaries
=#

#imports
using Pkg
Pkg.activate(".")
using StatsBase

#combine trajectories by concatenation
function concat_bins(fnames::Array{T, 1} where T,
                     si::Array{T, 1} where T<:Integer)
    n = length(fnames)
    so = copy(si)
    so[end] *= n
    i = zeros(si...)
    o = zeros(so...)
    for (fname, j) in zip(fnames, collect(1:si[end]:so[end]))
        read!(fname, i) 
        o[:, j:j+si[end]-1] = i
    end
    return o
end

#combine histograms by adding
function add_bins(fnames::Array{T, 1} where T,
                  si::Array{T, 1} where T<:Integer)
    n = length(fnames)
    i = zeros(si...)
    o = zeros(si...)
    for fname in fnames
        read!(fname, i) 
        o .+= i
    end
    return o
end

#combine trajectories from a file to a histogram
function compute_histogram(arr)
    v = .~ isnan.(arr[1,:])
    l = size(arr)[1]%3==0 ? 3 : 2
    h = fit(Histogram, ([arr[i,v] for i in 1:l]...,), edges).weights
    if flip_last
        h = reverse(h, dims=length(size(h)))
    end
    return h
end

#take per-proc trajectories and compute combined histogram
function multi_traj_to_single_hist(nProc)
    for iter in 0 : Int(round(wFreq/dt)) : Int(round(tStop/dt))
        hist = zeros(nx, ny)
        fnames = [
            Printf.format(Printf.Format("%s.%010d_%04d.bin"), t_prefix, iter, i) 
            for i in 1:nProc
        ]
        traj = zeros(nDim, nTraj÷nProc)
        for fname in fnames
            read!(fname, traj)
            hist = hist + compute_histogram(traj)
        end
        hist = hist / nTraj / vol
        fout = Printf.format(Printf.Format("%s.%010d.data"), t_prefix, iter)
        io = open(fout, "w")
        write(io, hton.(convert(Array{Float32, 2}, hist)))
        close(io)
        println(sum(hist))
    end
end
