# This fixture needs only stdlib TOML/Test; standalone runs no Core/arrays/H5 model.
using TOML,Test
isdefined(@__MODULE__,:DiagnosticSerialization) || include("serialization.jl")
function serialization_checks()
    D=DiagnosticSerialization
    fixture=(Hbar=[1.0 -0.0 1e-300;2e200 -3.5 7.25],
        complex_field=ComplexF64[1+2im -3+4im;5-6im -0.0-0.0im],
        records=[(absolute=1.25,relative=0.0),(absolute=3.5,relative=2.0)],
        Ritz_real=[.5,.5],Ritz_imag=[.2,-.2],shape_fixture=[2,3])
    function decoded(record)
        shape=Tuple(record["shape"])
        rowmatrix(rows)=permutedims(hcat(rows...))
        out=record["format"]=="complex-matrix-rows-v1" ?
            complex.(rowmatrix(record["real_rows"]),rowmatrix(record["imag_rows"])) : rowmatrix(record["rows"])
        @test size(out)==shape
        return out
    end
    @testset "source-free exact matrix/vector TOML roundtrip" begin
        encoded=D.plain(fixture);parsed=TOML.parse(sprint(io->TOML.print(io,encoded)))
        for key in ("Hbar","complex_field")
            original=getproperty(fixture,Symbol(key));restored=decoded(parsed[key])
            @test reinterpret(UInt64,vec(restored))==reinterpret(UInt64,vec(original))
        end
        @test parsed["records"][2]["absolute"]==3.5
        @test parsed["Ritz_imag"]==[.2,-.2]
        mktempdir() do dir
            path=joinpath(dir,"fixture.toml");D.small_toml(path,fixture)
            @test TOML.parsefile(path)==parsed
            @test_throws ErrorException D.small_toml(path,fixture)
        end
        @test_throws ArgumentError D.plain(zeros(1,1,1))
    end
end
if abspath(PROGRAM_FILE)==@__FILE__
    ARGS==["--tiny-serialization"] || error("usage: --tiny-serialization")
    serialization_checks()
end
