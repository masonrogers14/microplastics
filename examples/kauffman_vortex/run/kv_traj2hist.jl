#grid parameters
const nx = Int(round(2*R/dx)) + 4
const ny = Int(round(2*R/dy)) + 4
const nz = 1
const dt = 1e-3
const vol = dx*dy*dz
const XG = dx .* collect(-nx/2:nx/2-1)
const XC = XG .+ dx/2
const YG = dy .* collect(-ny/2:ny/2-1)
const YC = YG .+ dy/2
const RF = dz .* collect(nz:-1:0) #length nz+1
const RC = (RF[1:end-1] .+ RF[2:end]) ./ 2

#histogram bin edges and flips
const edges = (vcat(XG, maximum(XG)+dx), vcat(YG, maximum(YG)+dy))
const flip_last = false

#save trajectories
function save_trajectories()
    io = open(t_fname, "w")
    write(io, solu_arr)
    close(io)
end

#compute and save histogram data in MITgcm format
function save_histogram()
    for j in 1:nOuts
        v = .~ isnan.(solu_arr[:,j,1])
        h_j = fit(Histogram, ([solu_arr[v,j,i] for i in 1:size(solu_arr)[3]]...,), edges).weights
        if flip_last
            h_j = reverse(h_j, dims=length(size(h_j)))
	end
        h_suffix = @sprintf ".%010d.data" Int(wFreq * (j-1) / dt)
        io = open(h_prefix*h_suffix, "w")
        write(io, hton.(convert(Array{Float32,2}, h_j/(vol*nTraj))))
        close(io)
    end
end

#save grid if necessary
function pack_grid()
    #parameters----------------------------------------
    ρ₀ = 1000

    #positions-----------------------------------------
    io = open(dir*"RC.data", "w")
    write(io, hton.(convert(Array{Float32,1}, RC)))
    close(io)

    io = open(dir*"RF.data", "w")
    write(io, hton.(convert(Array{Float32,1}, RF)))
    close(io)

    io = open(dir*"XC.data", "w")
    write(io, hton.(convert(Array{Float32,2}, XC * ones(ny)')))
    close(io)

    io = open(dir*"XG.data", "w")
    write(io, hton.(convert(Array{Float32,2}, XG * ones(ny)')))
    close(io)

    io = open(dir*"YC.data", "w")
    write(io, hton.(convert(Array{Float32,2}, ones(nx) * YC')))
    close(io)

    io = open(dir*"YG.data", "w")
    write(io, hton.(convert(Array{Float32,2}, ones(nx) * YC')))
    close(io)

    io = open(dir*"RC.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n     1,    1,    1,\n     1,    1,    1,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nz nz)
    close(io)

    io = open(dir*"RF.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n     1,    1,    1,\n     1,    1,    1,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nz+1 nz+1)
    close(io)

    io = open(dir*"XC.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"XG.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"YC.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"YG.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    #areas---------------------------------------------
    io = open(dir*"RAC.data", "w")
    write(io, hton.(convert(Array{Float32,2}, dx*dy*ones(nx,ny))))
    close(io)

    io = open(dir*"RAS.data", "w")
    write(io, hton.(convert(Array{Float32,2}, dx*dy*ones(nx,ny))))
    close(io)

    io = open(dir*"RAW.data", "w")
    write(io, hton.(convert(Array{Float32,2}, dx*dy*ones(nx,ny))))
    close(io)

    io = open(dir*"RAZ.data", "w")
    write(io, hton.(convert(Array{Float32,2}, dx*dy*ones(nx,ny))))
    close(io)

    io = open(dir*"RAC.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"RAS.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"RAW.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"RAZ.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    #distances-----------------------------------------
    io = open(dir*"DRC.data", "w")
    write(io, hton.(convert(Array{Float32,1}, dz*ones(nz+1))))
    close(io)

    io = open(dir*"DRF.data", "w")
    write(io, hton.(convert(Array{Float32,1}, dz*ones(nz))))
    close(io)

    io = open(dir*"DXC.data", "w")
    write(io, hton.(convert(Array{Float32,2}, dx*ones(nx,ny))))
    close(io)

    io = open(dir*"DXG.data", "w")
    write(io, hton.(convert(Array{Float32,2}, dx*ones(nx,ny))))
    close(io)

    io = open(dir*"DYC.data", "w")
    write(io, hton.(convert(Array{Float32,2}, dy*ones(nx,ny))))
    close(io)

    io = open(dir*"DYG.data", "w")
    write(io, hton.(convert(Array{Float32,2}, dy*ones(nx,ny))))
    close(io)

    io = open(dir*"DRC.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n     1,    1,    1,\n     1,    1,    1,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nz+1 nz+1)
    close(io)

    io = open(dir*"DRF.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n     1,    1,    1,\n     1,    1,    1,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nz nz)
    close(io)

    io = open(dir*"DXC.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"DXG.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"DYC.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"DYG.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    #misc----------------------------------------------
    io = open(dir*"Depth.data", "w")
    write(io, hton.(convert(Array{Float32,2}, ones(nx,ny))))
    close(io)

    io = open(dir*"PHrefC.data", "w")
    write(io, hton.(convert(Array{Float32,1}, -g*ρ₀*RC)))
    close(io)

    io = open(dir*"PHrefF.data", "w")
    write(io, hton.(convert(Array{Float32,1}, -g*ρ₀*RF)))
    close(io)

    io = open(dir*"hFacC.data", "w")
    write(io, hton.(convert(Array{Float32,3}, ones(nx,ny,nz))))
    close(io)

    io = open(dir*"hFacS.data", "w")
    write(io, hton.(convert(Array{Float32,3}, ones(nx,ny,nz))))
    close(io)

    io = open(dir*"hFacW.data", "w")
    write(io, hton.(convert(Array{Float32,3}, ones(nx,ny,nz))))
    close(io)

    io = open(dir*"Depth.meta", "w")
    write(io, @sprintf "nDims = [2]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny)
    close(io)

    io = open(dir*"PHrefC.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n     1,    1,    1,\n     1,    1,    1,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nz nz)
    close(io)

    io = open(dir*"PHrefF.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n     1,    1,    1,\n     1,    1,    1,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nz+1 nz+1)
    close(io)

    io = open(dir*"hFacC.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n %5d,    1, %4d,\n %5d,    1, %4d,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny nz nz)
    close(io)

    io = open(dir*"hFacS.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n %5d,    1, %4d,\n %5d,    1, %4d,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny nz nz)
    close(io)

    io = open(dir*"hFacW.meta", "w")
    write(io, @sprintf "nDims = [3]; \ndimList = [\n %5d,    1, %4d,\n %5d,    1, %4d,\n %5d,    1, %4d\n]; \ndataprec = ['float32']; \nnrecords = [1];" nx nx ny ny nz nz)
    close(io)

    #non-grid stuff------------------------------------
    io = open(dir*"available_diagnostics.log", "w")
    write(io, @sprintf "%-4d|%-8s|%-4d|       |SMR     MR|1               |plastic concentration\n" 1 "TRAC01" nz)
    write(io, @sprintf "%-4d|%-8s|  1 |       |SM      M1|m               |Surface Height Anomaly\n" 2 "ETAN")
    write(io, @sprintf "%-4d|%-8s|  1 |       |SM      M1|m^2/s^2         |Bottom Pressure Pot.(p/rho) Anomaly\n" 3 "PHIBOT")
    write(io, @sprintf "%-4d|%-8s|%-4d|       |SMR     MR|psu             |Salinity\n" 4 "SALT" nz)
    write(io, @sprintf "%-4d|%-8s|%-4d|       |SMR     MR|1               |Potential Temperature\n" 5 "THETA" nz)
    close(io)

    for j in 1:nOuts
        nt = Int(wFreq*(j-1)/dt)
        h_suffix = @sprintf ".%010d.meta" nt
        io = open(h_prefix*h_suffix, "w")
        write(io, @sprintf "nDims = [3]; \ndimList = [\n  %d, 1, %d,\n  %d, 1, %d,\n  %d, 1, %d \n]; \ndataprec = ['float32']; \nnrecords = [1]; \ntimeStepNumber = [0]; \ntimeInterval = [%.8e]; \nnFlds = [1]; \nfldList = { \n'TRAC01  ' \n};" nx nx ny ny nz nz nt*dt)
        close(io)
    end
end
