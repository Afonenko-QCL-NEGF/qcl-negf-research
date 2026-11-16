# Source-free metadata controls only. Does NOT manufacture a scientific proof.
using Test, TOML, SHA
function certificate_gate_controls(check)
    @test_throws ErrorException check(Dict("execution_authorized"=>false),".",Dict())
    mktempdir() do dir
        cert=joinpath(dir,"cert.toml"); proof=joinpath(dir,"proof.toml")
        open(proof,"w") do io; TOML.print(io,Dict("mode"=>"all")); end
        h(p)=bytes2hex(open(sha256,p))
        c=Dict{String,Any}("execution_authorized"=>true,"pd_certificate"=>cert,"pd_certificate_input"=>proof,
            "pd_certificate_input_sha256"=>h(proof),"expected_NE"=>30704,"expected_Nk"=>11,"expected_Nb"=>5)
        for record in (Dict("execution"=>"complete","mode"=>"benchmark5000","all_blocks_strict_PD_certified"=>false),
                       Dict("execution"=>"bounded_incomplete","mode"=>"all","all_blocks_strict_PD_certified"=>false),
                       Dict("execution"=>"failed","mode"=>"all","all_blocks_strict_PD_certified"=>true),
                       Dict("execution"=>"complete","mode"=>"all","all_blocks_strict_PD_certified"=>true,
                           "expected_all_blocks"=>337744,"requested_blocks"=>5000),
                       Dict("execution"=>"complete","mode"=>"all","all_blocks_strict_PD_certified"=>true,
                           "expected_all_blocks"=>337744,"requested_blocks"=>337744,
                           "results"=>Dict("checked_blocks"=>5000,"counts"=>Dict("pd"=>5000))))
            open(cert,"w") do io; TOML.print(io,record); end
            c["pd_certificate_sha256"]=h(cert)
            @test_throws ErrorException check(c,dir,Dict())
        end
    end
end
