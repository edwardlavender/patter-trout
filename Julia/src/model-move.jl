struct StateXYD <: Patter.State
    map_value::Float64
    x::Float64
    y::Float64
    angle::Float64  
end

struct ModelMoveXYD{T, U, V, W} <: ModelMove
    map::T
    mobility::U
    dbn_length::V
    dbn_angle_delta::W
end 

function Patter.simulate_step(state::StateXYD, model_move::ModelMoveXYD, t::Int64)
    length = rand(model_move.dbn_length)
    angle  = state.angle + rand(model_move.dbn_angle_delta)
    x      = state.x + length * cos(angle)
    y      = state.y + length * sin(angle)
    map_value = Patter.extract(model_move.map, x, y)
    StateXYD(map_value, x, y, angle)
end 

function Patter.logpdf_step(state_from::StateXYD, state_to::StateXYD, model_move::ModelMoveXYD, t::Int64, length::Float64, angle::Float64) 
    angle_delta = abs_angle_difference(angle, state_from.angle)
    logpdf(model_move.dbn_length, length) + logpdf(model_move.dbn_angle_delta, angle_delta)
end 

function Patter.states_init(state_type::Type{StateXYD}, coords)
    N = nrow(coords)
    coords.angle = rand(N) .* 2 .* π
    return coords
end 