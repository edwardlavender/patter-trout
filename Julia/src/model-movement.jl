# State and movement model structures for 2D-CRW 
# (e.g., as implemented by glatos::crw())

# Custom state (X, Y, D)
struct StateXYD <: Patter.State
    # Map value
    map_value::Float64   
    # Coordinates 
    x::Float64
    y::Float64
    # Direction
    angle::Float64  
end 

# Custom movement model structure for fixed step length
struct ModelMoveXYD{T, U, V} <: Patter.ModelMove
    # Environment
    map::T
    # Step length (fixed)
    length::U
    # Change in turning angle
    dbn_angle_delta::V
end 

# simulate_step() function for the particle filter
function Patter.simulate_step(state::StateXYD, model_move::ModelMoveXYD, t::Int64)
    # Simulate length (fixed)
    length = model_move.length
    # Simulate angle 
    # * Note that the change in angle is converted from degrees to radians
    angle  = state.angle + (rand(model_move.dbn_angle_delta) * (π / 180))
    x      = state.x + length * cos(angle)
    y      = state.y + length * sin(angle)
    map_value = extract(model_move.map, x, y)
    StateXYD(map_value, x, y, angle)
end 

# logpdf_step() function for the particle smoother
function Patter.logpdf_step(state_from::StateXYZD, state_to::StateXYZD, model_move::ModelMoveXYZD, t::Int64, length::Float64, angle::Float64) 
    # Compute change in angle 
    angle_delta = abs_angle_difference(angle, state_from.angle)
    # Sum up logpdfs
    # (fixed step length + change in angle)
    0.0 + logpdf(model_move.dbn_angle_delta, angle_delta)
end 