fnSimulate <- function(scenario, N, var, scale, nu, intercept, pnum, anum, unif){
  
  for(i in 1:N){
    print(i)
    
    # Simulate continuous surface
    r <- latt_generation(xlim, ylim, by, intercept, nu, scale, var)
    truesurface <- raster::extract(r, as.matrix(st_coordinates(dppoint)))
    
    # Sample observations in points and areas
    if(unif){p1 <- punifsample(pnum, r)}
    if(!unif){
      loct <- pcoxsample(pnum, r)
    p1 <- datagenerator(loct, r)
    p1 <- p1 %>% st_as_sf(coords = c("x", "y"), dim = "XY") %>% st_cast("POINT")
    }
    a1 <- areasample(anum, r)
    
    # Fit model with melding and downscaler approaches
    mesh <- fnCreateMesh(p1, boundaryregion)
    
    meld <- fnPredictMelding(depoint = p1, dearea = a1, dppoint = dppoint, dparea = NULL,
                             boundaryregion = boundaryregion, mesh = mesh)
    
    down <- fnPredictDown(depoint = p1, dearea = a1, dppoint = dppoint, dparea = NULL,
                          boundaryregion = boundaryregion, mesh = mesh)
    
    # Save truesurface
    write.csv(cbind(true = truesurface,meld = meld[[1]]$pred_mean,down = down[[1]]$pred_mean), paste0("results/", scenario, "results", i, ".csv"), row.names = FALSE)
    
    
  }}

fnSimulateCov <- function(scenario, N, var, scale, nu, intercept,fix_betax, pnum, anum, unif){
  for(i in 1:N){
    print(i)
    
    # Simulate continuous surface
    r_random <- latt_generation(xlim, ylim, by, intercept, nu, scale, var)
    r_fix <- raster(matrix(rep(x*fix_betax, length(y)), byrow = T, nrow = length(y)))
    r <- r_random + r_fix
    truesurface <- raster::extract(r, as.matrix(st_coordinates(dppoint)))
    
    # Sample observations in points and areas
    if(unif){p1 <- punifsample(pnum, r)}
    if(!unif){
      loct <- pcoxsample(pnum, r)
      p1 <- datagenerator(loct, r)
      p1 <- p1 %>% st_as_sf(coords = c("x", "y"), dim = "XY") %>% st_cast("POINT")
    }
    a1 <- areasample(anum, r)
    
    # Fit model with melding and downscaler approaches
    mesh <- fnCreateMesh(p1, boundaryregion)
    
    meld <- fnPredictMeldingCov(depoint = p1, dearea = a1, dppoint = dppoint, dparea = NULL,
                             boundaryregion = boundaryregion, mesh = mesh)
    
    down <- fnPredictDownCov(depoint = p1, dearea = a1, dppoint = dppoint, dparea = NULL,
                          boundaryregion = boundaryregion, mesh = mesh)
    
    # Save truesurface
    write.csv(cbind(true = truesurface,meld = meld[[1]]$pred_mean,down = down[[1]]$pred_mean), paste0("results/", scenario, "results", i, ".csv"), row.names = FALSE)
    
    
  }
}

 
fnCalculateMSE <- function(scenario, N){
  mse_meld <- NULL
  mse_down <- NULL
  for(i in 1:N){
    print(i)
    results <- read.csv(paste0("results/", scenario, "results", i, ".csv"))
    mse_meld <- c(mse_meld, Metrics::mse(results[ ,1], results[ ,2]))
    mse_down <- c(mse_down, Metrics::mse(results[ ,1], results[ ,3]))
  }
  data <- as.data.frame(cbind(MSE = c(mse_meld, mse_down), Model = rep(c('Meld', 'Downscaler'), each = N)))
  write.csv(data, paste0("results/mse/", scenario, "MSE", N, ".csv"), row.names = FALSE)
}
