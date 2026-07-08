convert_data(data::N, backend) where {names,N<:NamedTuple{names}} =
    NamedTuple{names}(convert_array(d, backend) for d in data)

function parse_ac_power_data(filename, T = Float64; from = nothing)
    path = isfile(filename) ? filename : joinpath(ExaPowerIO.get_path(:pglib), filename)
    @info "Loading power case file"
    return PowerIO.parse_ac_power_data(path; from = from, T = T)
end

# Resolve a case name to a parsed network. The same network feeds both the
# powerdata NamedTuple and the LoadSeries so their bus order aligns.
function parse_mp_network(filename; from = nothing)
    path = isfile(filename) ? filename : joinpath(ExaPowerIO.get_path(:pglib), filename)
    @info "Loading power case file"
    return PowerIO.parse_file(path; from = from)
end
