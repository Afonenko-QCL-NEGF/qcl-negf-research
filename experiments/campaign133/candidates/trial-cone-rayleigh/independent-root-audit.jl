using HDF5, LinearAlgebra, TOML
base=joinpath(@__DIR__,"attempt01")
summary=TOML.parsefile(joinpath(base,"summary.toml"))
@assert summary["execution"]=="complete"
for (name,record) in summary["blocks"]
    X,v=h5open(joinpath(base,name*".h5"),"r") do f
        complex.(read(f["raw_real"]),read(f["raw_imag"])),
        complex.(read(f["float_eigenvector_real"]),read(f["float_eigenvector_imag"]))
    end
    # Base conversion independently replaces the helper's IEEE bit decoder.
    Q=Complex{Rational{BigInt}}.(X)
    w=Complex{Rational{BigInt}}.(v)
    result=dot(w,((Q+Q')/2)*w)/dot(w,w)
    @assert iszero(imag(result))
    q=real(result)
    expected=record["exact_Hermitian_part_raw_entries"]
    @assert q==parse(BigInt,expected["numerator"])//parse(BigInt,expected["denominator"])
    @assert q<0
    @assert q < -Rational{BigInt}(record["guard_allowance"])
    println(name,": independent Base-rational matrix quotient matches; negative beyond recorded allowance")
end
