# Source-free TOML representation. No QCL/Julia package/model side effects.
module DiagnosticSerialization
using TOML
export plain,small_toml
plain(x::NamedTuple)=Dict(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractDict)=Dict(string(k)=>plain(v) for (k,v) in x)
plain(x::Complex)=Dict("format"=>"complex-scalar-v1","real"=>real(x),"imag"=>imag(x))
plain(x::AbstractMatrix{<:Complex})=Dict("format"=>"complex-matrix-rows-v1","shape"=>collect(size(x)),
    "real_rows"=>[[real(x[i,j]) for j in axes(x,2)] for i in axes(x,1)],
    "imag_rows"=>[[imag(x[i,j]) for j in axes(x,2)] for i in axes(x,1)])
plain(x::AbstractMatrix)=Dict("format"=>"matrix-rows-v1","shape"=>collect(size(x)),
    "rows"=>[[plain(x[i,j]) for j in axes(x,2)] for i in axes(x,1)])
plain(x::AbstractVector)=map(plain,x)
plain(x::AbstractArray)=throw(ArgumentError("TOML tensor not supported; use declared H5 raw writer"))
plain(x)=x
function small_toml(path,x)
    ispath(path) && error("refusing overwrite $path")
    open(path,"w") do io;TOML.print(io,plain(x));flush(io);end
end
end
