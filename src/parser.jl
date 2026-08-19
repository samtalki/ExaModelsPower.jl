# `map`, not `NamedTuple{names}(generator)`. The generator form builds the
# tuple from an iterator whose element types inference cannot pin one by one,
# so the result is a `NamedTuple{names}` with no value types — concrete names,
# abstract contents. `map` over a NamedTuple is inferable elementwise.
# Recurses into a nested NamedTuple. The multi-period arguments group their
# N-expanded bound matrices under one field, because there are enough of them
# to push a flat tuple past the 31-field limit above.
convert_data(data::NamedTuple, backend) = map(d -> _convert_field(d, backend), data)
_convert_field(d, backend) = convert_array(d, backend)
_convert_field(d::NamedTuple, backend) = convert_data(d, backend)

# A case name that is not a path on disk is looked up in the bundled PGLib-OPF
# directory, which is the one thing ExaPowerIO is still here for.
_case_path(filename) =
    isfile(filename) ? filename : joinpath(ExaPowerIO.get_path(:pglib), filename)

# `T` is a TYPE PARAMETER, not a value of type `Type`. As a plain positional
# argument it arrives typed `Type` — abstract — and then the parse cannot
# specialize, so its result is an unparameterized NamedTuple and every field
# read off it widens to `AbstractArray{BusRow{T}} where T`. Written `::Type{T}`
# it is a static parameter and the whole chain stays concrete.
parse_ac_power_data(filename) = parse_ac_power_data(filename, Float64)
function parse_ac_power_data(filename, ::Type{T}; from = nothing) where {T}
    @info "Loading power case file"
    return PowerIO.parse_ac_power_data(_case_path(filename); from = from, T = T)
end

# Resolve a case name to a parsed network. The same network feeds both the
# powerdata NamedTuple and the LoadSeries, so their bus order aligns by
# construction rather than by assumption.
function parse_mp_network(filename; from = nothing)
    @info "Loading power case file"
    return PowerIO.parse_file(_case_path(filename); from = from)
end
