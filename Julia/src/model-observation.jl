# (TEMPORARY) ModelObsAcousticContainer structure & methods

struct ModelObsAcousticContainer <: Patter.ModelObs
  sensor_id::Int64
  receiver_x::Float64
  receiver_y::Float64
  max_dist::Float64
end

function Patter.logpdf_obs(state::State, model::ModelObsAcousticContainer, t::Int64, obs::Int64)
  # Calculate distance between particle (state) and receiver
  dist = Patter.distance(state.x, state.y, model.receiver_x, model.receiver_y)
  # Only particles within max_dist are permitted
  # * max_dist is a pre-calculated field in model
  # * (max_dist = receiver_gamma + (receiver_timestep - t) * mobility)
  ifelse(dist <= model.max_dist, 0.0, -Inf)
end
