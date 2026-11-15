using HDF5, TOML
length(ARGS)==2 || error("usage: read_markers.jl STATE.h5 NEW_OUTPUT_DIRECTORY")
ispath(ARGS[2]) && error("output exists")
mkpath(ARGS[2])
h5open(ARGS[1],"r") do f
    g=f["physical_markers"]
    fields=sort([k for k in keys(g) if !(k in ("channels","collisions"))])
    arrays=Dict(k=>read(g[k]) for k in fields)
    selected=findall(!iszero,arrays["available"])
    open(joinpath(ARGS[2],"physical_markers.csv"),"w") do io
        println(io,join(fields,","))
        for row in selected
            println(io,join([string(arrays[k][row]) for k in fields],","))
        end
    end
    lastrow=last(selected)
    final=Dict(k=>(arrays[k][lastrow] isa AbstractString ? String(arrays[k][lastrow]) : arrays[k][lastrow]) for k in fields)
    open(joinpath(ARGS[2],"final_markers.toml"),"w") do io
        TOML.print(io,final)
    end
    for name in ("channels","collisions")
        child=g[name]
        columns=sort(collect(keys(child)))
        vals=Dict(k=>read(child[k]) for k in columns)
        open(joinpath(ARGS[2],name*".csv"),"w") do io
            println(io,join(columns,","))
            for row in eachindex(vals[first(columns)])
                println(io,join([string(vals[k][row]) for k in columns],","))
            end
        end
    end
end
println("MarkerExtractionOK")
