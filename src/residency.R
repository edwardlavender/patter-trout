if (Sys.getenv("JULIA_SESSION") == "FALSE") {
  
  #' Quick overall residency estimation
  # * id:    individual_id
  # * map:   raster with regions defined
  # * paths: data.table of simulated paths
  
  qresidency <- function(id, map, paths) {
    
    path <- paths[individual_id == id, ]
    smo  <- here_output("particles", glue("smo-{id}.qs"))
    if (!file.exists(smo)) {
      return(NULL)
    }
    smo  <- qs::qread(smo)
    
    #### True residency in each region
    path_res <- 
      path |> 
      # Replace region column for consistency
      mutate(region = terra::extract(map, cbind(x, y))[, 1]) |>
      group_by(region) |> 
      summarise(n = n()) |> 
      ungroup() |>
      mutate(perc = (n / sum(n)) * 100, 
             estimate = "path", 
             individual_id = id) |> 
      select("individual_id", "region", "estimate", count = "n", "perc") |>
      as.data.table()
    
    #### Particle residency estimates
    part_res <- 
      smo$states |> 
      mutate(region = terra::extract(map, cbind(x, y))[, 1]) |>
      group_by(region) |> 
      summarise(n = n()) |> 
      ungroup() |>
      mutate(perc = (n / sum(n)) * 100, 
             estimate = "smoother", 
             individual_id = id) |> 
      select("individual_id", "region", "estimate", count = "n", "perc") |>
      as.data.table()
    
    #### Collect data
    res <- rbind(path_res, part_res)
    qs::qsave(res, here_output("residency", "qresidency", glue("{id}.qs")))
    
    NULL
    
  }
  
}