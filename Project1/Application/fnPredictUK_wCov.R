fnPredictMeldingUK_1<- function( depoint = NULL, dearea = NULL, dppoint = NULL, dparea = NULL, boundaryregion,
                               mesh = NULL, priorspdesigma = NULL, priorspderange = NULL){
  formula <- y ~ 0 + Intercept + temp_mn + precptn + rod_scl + rod_prp + f(s,model = spde)
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Check inputs
  fnCheckInputsMelding(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  de1ToF <- !is.null(de1)
  de2ToF <- !is.null(de2)
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  # Mesh
  if(is.null(mesh)){
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  
 
    indexs <- inla.spde.make.index('s', spde$n.spde)
  
  
  
  # Projection matrices for points (estimation point and prediction point)
  if(de1ToF){Ae1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(de1)[ , c(1,2)]))}
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dp1)[ , c(1,2)]))}
  
  # Create projection matrix A for areas (estimation area and prediction area)
  if(de2ToF){Ae2 <- fnProjectionMatrixArea(de2, mesh)}
  if(dp2ToF){Ap2 <- fnProjectionMatrixArea(dp2, mesh)}
  
  # Specify formula melding
  
  
  # Create stk.full
  stk.e1 = NULL
  stk.e2 = NULL
  stk.p1 = NULL
  stk.p2 = NULL
  
  # estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
  if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$value), A = list(1,1,1,1,1, Ae1), effects = list(data.frame(Intercept = rep(1, nrow(de1))),temp_mn = de1$temp_mn ,
                                                                                                                   precptn = de1$precptn, rod_scl = de1$rod_scl, rod_prp = de1$rod_prp, s = indexs$s))}
  if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$value), A = list(1,1,1,1,1, Ae2), effects = list(data.frame(Intercept = rep(1, nrow(de2))),temp_mn = de2$temp_mn,
                                                                                                                                                                        precptn = de2$precptn, rod_scl = de2$rod_scl, rod_prp = de2$rod_prp,
                                                                                                                                                                        s = indexs$s))}
  
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA),A = list(1,1,1,1,1, Ap1), effects = list(data.frame(Intercept = rep(1, nrow(dp1))),temp_mn = dp1$temp_mn,
                                                                                                              precptn = dp1$precptn, rod_scl = dp1$rod_scl, rod_prp = dp1$rod_prp,
                                                                                                              s = indexs$s))}
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(rep(1,fix_num), Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2))), cov = dp2$cov, s = indexs$s))}
  
  # construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
  stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])
  # Call inla()
  
  res <- inla(formula, family = "gaussian", data = inla.stack.data(stk.full),
              control.compute = list(dic = TRUE),
              control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))
  dic <- res$dic$dic
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions(stk.full, res, "pred1", dp1)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions(stk.full, res, "pred2", dp2)}
  
  #return(list(dp1,dp2,dic))
  return(res)
}

fnPredictMeldingUK_2<- function( depoint = NULL, dearea = NULL, dppoint = NULL, dparea = NULL, boundaryregion,
                               mesh = NULL, priorspdesigma = NULL, priorspderange = NULL){
  formula <- y ~ 0 + Intercept + temp_mn + precptn  + rod_prp + f(s,model = spde)
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Check inputs
  fnCheckInputsMelding(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  de1ToF <- !is.null(de1)
  de2ToF <- !is.null(de2)
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  # Mesh
  if(is.null(mesh)){
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  
  
  indexs <- inla.spde.make.index('s', spde$n.spde)
  
  
  
  # Projection matrices for points (estimation point and prediction point)
  if(de1ToF){Ae1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(de1)[ , c(1,2)]))}
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dp1)[ , c(1,2)]))}
  
  # Create projection matrix A for areas (estimation area and prediction area)
  if(de2ToF){Ae2 <- fnProjectionMatrixArea(de2, mesh)}
  if(dp2ToF){Ap2 <- fnProjectionMatrixArea(dp2, mesh)}
  
  # Specify formula melding
  
  
  # Create stk.full
  stk.e1 = NULL
  stk.e2 = NULL
  stk.p1 = NULL
  stk.p2 = NULL
  
  # estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
  if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$value), A = list(1,1,1,1, Ae1), effects = list(data.frame(Intercept = rep(1, nrow(de1))),temp_mn = de1$temp_mn ,
                                                                                                                     precptn = de1$precptn, rod_prp = de1$rod_prp, s = indexs$s))}
  if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$value), A = list(1,1,1,1, Ae2), effects = list(data.frame(Intercept = rep(1, nrow(de2))),temp_mn = de2$temp_mn,
                                                                                                                     precptn = de2$precptn,  rod_prp = de2$rod_prp,
                                                                                                                     s = indexs$s))}
  
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA),A = list(1,1,1,1, Ap1), effects = list(data.frame(Intercept = rep(1, nrow(dp1))),temp_mn = dp1$temp_mn,
                                                                                                              precptn = dp1$precptn,  rod_prp = dp1$rod_prp,
                                                                                                              s = indexs$s))}
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(rep(1,fix_num), Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2))), cov = dp2$cov, s = indexs$s))}
  
  # construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
  stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])
  # Call inla()
  
  res <- inla(formula, family = "gaussian", data = inla.stack.data(stk.full),
              control.compute = list(dic = TRUE),
              control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))
  dic <- res$dic$dic
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions(stk.full, res, "pred1", dp1)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions(stk.full, res, "pred2", dp2)}
  
  #return(list(dp1,dp2,dic))
  return(res)
}

fnPredictMeldingUK_3<- function( depoint = NULL, dearea = NULL, dppoint = NULL, dparea = NULL, boundaryregion,
                                 mesh = NULL, priorspdesigma = NULL, priorspderange = NULL){
  formula <- y ~ 0 + Intercept + temp_mn + precptn + rod_scl + f(s,model = spde)
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Check inputs
  fnCheckInputsMelding(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  de1ToF <- !is.null(de1)
  de2ToF <- !is.null(de2)
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  # Mesh
  if(is.null(mesh)){
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  
  
  indexs <- inla.spde.make.index('s', spde$n.spde)
  
  
  
  # Projection matrices for points (estimation point and prediction point)
  if(de1ToF){Ae1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(de1)[ , c(1,2)]))}
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dp1)[ , c(1,2)]))}
  
  # Create projection matrix A for areas (estimation area and prediction area)
  if(de2ToF){Ae2 <- fnProjectionMatrixArea(de2, mesh)}
  if(dp2ToF){Ap2 <- fnProjectionMatrixArea(dp2, mesh)}
  
  # Specify formula melding
  
  
  # Create stk.full
  stk.e1 = NULL
  stk.e2 = NULL
  stk.p1 = NULL
  stk.p2 = NULL
  
  # estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
  if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$value), A = list(1,1,1,1, Ae1), effects = list(data.frame(Intercept = rep(1, nrow(de1))),temp_mn = de1$temp_mn ,
                                                                                                                     precptn = de1$precptn, rod_scl = de1$rod_scl, s = indexs$s))}
  if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$value), A = list(1,1,1,1, Ae2), effects = list(data.frame(Intercept = rep(1, nrow(de2))),temp_mn = de2$temp_mn,
                                                                                                                     precptn = de2$precptn, rod_scl = de2$rod_scl, 
                                                                                                                     s = indexs$s))}
  
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA),A = list(1,1,1,1, Ap1), effects = list(data.frame(Intercept = rep(1, nrow(dp1))),temp_mn = dp1$temp_mn,
                                                                                                              precptn = dp1$precptn, rod_scl = dp1$rod_scl,
                                                                                                              s = indexs$s))}
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(rep(1,fix_num), Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2))), cov = dp2$cov, s = indexs$s))}
  
  # construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
  stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])
  # Call inla()
  
  res <- inla(formula, family = "gaussian", data = inla.stack.data(stk.full),
              control.compute = list(dic = TRUE),
              control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))
  dic <- res$dic$dic
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions(stk.full, res, "pred1", dp1)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions(stk.full, res, "pred2", dp2)}
  
  return(list(dp1,dp2,dic))
}

fnPredictMeldingUK_4<- function( depoint = NULL, dearea = NULL, dppoint = NULL, dparea = NULL, boundaryregion,
                                 mesh = NULL, priorspdesigma = NULL, priorspderange = NULL){
  formula <- y ~ 0 + Intercept + temp_mn + precptn + f(s,model = spde)
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Check inputs
  fnCheckInputsMelding(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  de1ToF <- !is.null(de1)
  de2ToF <- !is.null(de2)
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  # Mesh
  if(is.null(mesh)){
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  
  
  indexs <- inla.spde.make.index('s', spde$n.spde)
  
  
  
  # Projection matrices for points (estimation point and prediction point)
  if(de1ToF){Ae1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(de1)[ , c(1,2)]))}
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dp1)[ , c(1,2)]))}
  
  # Create projection matrix A for areas (estimation area and prediction area)
  if(de2ToF){Ae2 <- fnProjectionMatrixArea(de2, mesh)}
  if(dp2ToF){Ap2 <- fnProjectionMatrixArea(dp2, mesh)}
  
  # Specify formula melding
  
  
  # Create stk.full
  stk.e1 = NULL
  stk.e2 = NULL
  stk.p1 = NULL
  stk.p2 = NULL
  
  # estimation point (de1), estimation area (de2), prediction point (dp1), prediction area (dp2)
  if(de1ToF){stk.e1 <- inla.stack(tag = "est1", data = list(y = de1$value), A = list(1,1,1, Ae1), effects = list(data.frame(Intercept = rep(1, nrow(de1))),temp_mn = de1$temp_mn ,
                                                                                                                     precptn = de1$precptn, s = indexs$s))}
  if(de2ToF){stk.e2 <- inla.stack(tag = "est2", data = list(y = de2$value), A = list(1,1,1, Ae2), effects = list(data.frame(Intercept = rep(1, nrow(de2))),temp_mn = de2$temp_mn,
                                                                                                                     precptn = de2$precptn, 
                                                                                                                     s = indexs$s))}
  
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA),A = list(1,1,1, Ap1), effects = list(data.frame(Intercept = rep(1, nrow(dp1))),temp_mn = dp1$temp_mn,
                                                                                                              precptn = dp1$precptn, 
                                                                                                              s = indexs$s))}
  if(dp2ToF){stk.p2 <- inla.stack(tag = "pred2", data = list(y = NA), A = list(rep(1,fix_num), Ap2), effects = list(data.frame(b0 = rep(1, nrow(dp2))), cov = dp2$cov, s = indexs$s))}
  
  # construct stack full with the data we have # stk.full <- inla.stack(stk.e1, stk.p1)
  stk.full <- do.call(inla.stack, list(stk.e1, stk.e2, stk.p1, stk.p2)[c(de1ToF, de2ToF, dp1ToF, dp2ToF)])
  # Call inla()
  
  res <- inla(formula, family = "gaussian", data = inla.stack.data(stk.full),
              control.compute = list(dic = TRUE),
              control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full)))
  dic <- res$dic$dic
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions(stk.full, res, "pred1", dp1)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions(stk.full, res, "pred2", dp2)}
  
  return(list(dp1,dp2,dic))
}

fnPredictDownUK_1 <- function(depoint, dearea, dppoint = NULL, dparea = NULL, boundaryregion,
                          mesh = NULL, priorspdesigma = NULL, priorspderange = NULL){
  
  
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Set avalue and pvalue instead of value
  de1$pvalue <- de1$value
  de2$avalue <- de2$value
  dp1$pvalue <- dp1$value
  dp2$avalue <- dp2$value  
  
  # Check inputs
  
  # check function
  # TO DO : check crs of input is crsproj
  fnCheckInputsDown(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  
  # Mesh
  if (is.null(mesh) == T) {
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  indexs <- inla.spde.make.index("s", spde$n.spde)
  indexs1 <- inla.spde.make.index("s1", spde$n.spde)
  
  
  
  # Match Areal data and point data, Match dppoint and dearea
  locin <- st_join(de1, de2, left = F)
  
  ### PAULA: Keep both points and areas are inside study region. Should join by de2 or boundaryregion?
  # If TRUE, predict in points
  if(dp1ToF){
    locin_pred <- st_join(dp1, de2, left = FALSE)
  }
  
  # If TRUE, predict in areas. Construct prediction surface
  if(dp2ToF){
    # Bounding box
    bb <- unname(attributes(st_geometry(boundaryregion))$bbox)
    # Grid
    x <- seq(bb[1] - 1, bb[3] + 1, length.out = 50)
    y <- seq(bb[2] - 1, bb[4] + 1, length.out = 50)
    coop <- expand.grid(x, y)
    coop_sf <- sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = st_crs(de1))
    # Keep points inside
    dpcontsurface <- coop_sf %>% st_join(boundaryregion, left = FALSE) %>% st_join(de2, left = FALSE)
  }
  
  
  
  # Create projection matrices A
  # Matrix A estimation
  Ae <- inla.spde.make.A(mesh, loc = as.matrix(st_coordinates(locin[, 1]))[, c(1, 2)])
  # Matrix A points prediction
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(locin_pred)[, c(1,2)]))}
  # Matrix A areas prediction (need to predict in a continuous surface first, dense points)
  if(dp2ToF){Ap2 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dpcontsurface[, 1])))}
  
  # Create stk.full
  # btilde0(x) = b0 + b0(x) (Ap1 + indexs1$s1)
  # btilde1(x) = b1 + b1(x) (Ap1 * locin_pred$avalue + indexs$s)
  # stack estimation
  
  stk.p1 = NULL
  stk.p2 = NULL
  stk.e <- inla.stack(tag = "est", data = list(y = locin$pvalue),
                      A = list(1,1,1,1,1, Ae * locin$avalue, Ae),
                      effects = list(data.frame(Intercept = rep(1, nrow(locin)), X = locin$avalue),temp_mn = de1$temp_mn ,
                                     precptn = de1$precptn, rod_scl = de1$rod_scl, rod_prp = de1$rod_prp, s = indexs$s, s1 = indexs1$s1))
  
 

  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA),
                                  A = list(1,1,1,1,1, Ap1 * locin_pred$avalue, Ap1),
                                  effects = list(data.frame(Intercept = rep(1, nrow(locin_pred)), X = locin_pred$avalue),temp_mn = dp1$temp_mn,
                                                 precptn = dp1$precptn, rod_scl = dp1$rod_scl, rod_prp = dp1$rod_prp, s = indexs$s, s1 = indexs1$s1))}
  # stack prediction areas
  
  # construct stack full with the data we have
  stk.full <-  do.call(inla.stack, list(stk.e, stk.p1, stk.p2)[c( T,dp1ToF,dp2ToF)])
  
  
  # Specify formula downscaling
  formula <- y ~ 0 + Intercept + X + temp_mn + precptn + rod_scl + rod_prp + f(s,model = spde)+ f(s1, model = spde)
  
  # Call inla()
  res <- inla(formula, data = inla.stack.data(stk.full), 
              control.compute=list(dic=TRUE),
              control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full)))
  
  dic <- res$dic$dic
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions_p(stk.full, res, "pred1", locin_pred)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions_a(stk.full, res, "pred2", dpcontsurface, dp2)}
  
  return(list(dp1, dp2,dic))
}

fnPredictDownUK_2 <- function(depoint, dearea, dppoint = NULL, dparea = NULL, boundaryregion,
                              mesh = NULL, priorspdesigma = NULL, priorspderange = NULL){
  
  
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Set avalue and pvalue instead of value
  de1$pvalue <- de1$value
  de2$avalue <- de2$value
  dp1$pvalue <- dp1$value
  dp2$avalue <- dp2$value  
  
  # Check inputs
  
  # check function
  # TO DO : check crs of input is crsproj
  fnCheckInputsDown(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  
  # Mesh
  if (is.null(mesh) == T) {
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  indexs <- inla.spde.make.index("s", spde$n.spde)
  indexs1 <- inla.spde.make.index("s1", spde$n.spde)
  
  
  
  # Match Areal data and point data, Match dppoint and dearea
  locin <- st_join(de1, de2, left = F)
  
  ### PAULA: Keep both points and areas are inside study region. Should join by de2 or boundaryregion?
  # If TRUE, predict in points
  if(dp1ToF){
    locin_pred <- st_join(dp1, de2, left = FALSE)
  }
  
  # If TRUE, predict in areas. Construct prediction surface
  if(dp2ToF){
    # Bounding box
    bb <- unname(attributes(st_geometry(boundaryregion))$bbox)
    # Grid
    x <- seq(bb[1] - 1, bb[3] + 1, length.out = 50)
    y <- seq(bb[2] - 1, bb[4] + 1, length.out = 50)
    coop <- expand.grid(x, y)
    coop_sf <- sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = st_crs(de1))
    # Keep points inside
    dpcontsurface <- coop_sf %>% st_join(boundaryregion, left = FALSE) %>% st_join(de2, left = FALSE)
  }
  
  
  
  # Create projection matrices A
  # Matrix A estimation
  Ae <- inla.spde.make.A(mesh, loc = as.matrix(st_coordinates(locin[, 1]))[, c(1, 2)])
  # Matrix A points prediction
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(locin_pred)[, c(1,2)]))}
  # Matrix A areas prediction (need to predict in a continuous surface first, dense points)
  if(dp2ToF){Ap2 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dpcontsurface[, 1])))}
  
  # Create stk.full
  # btilde0(x) = b0 + b0(x) (Ap1 + indexs1$s1)
  # btilde1(x) = b1 + b1(x) (Ap1 * locin_pred$avalue + indexs$s)
  # stack estimation
  
  stk.p1 = NULL
  stk.p2 = NULL
  stk.e <- inla.stack(tag = "est", data = list(y = locin$pvalue),
                      A = list(1,1,1,1,Ae * locin$avalue, Ae),
                      effects = list(data.frame(Intercept = rep(1, nrow(locin)), X = locin$avalue),temp_mn = de1$temp_mn ,
                                     precptn = de1$precptn, rod_scl = de1$rod_scl, s = indexs$s, s1 = indexs1$s1))
  
  
  
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA),
                                  A = list(1,1,1,1, Ap1 * locin_pred$avalue, Ap1),
                                  effects = list(data.frame(Intercept = rep(1, nrow(locin_pred)), X = locin_pred$avalue),temp_mn = dp1$temp_mn,
                                                 precptn = dp1$precptn, rod_scl = dp1$rod_scl,  s = indexs$s, s1 = indexs1$s1))}
  # stack prediction areas
  
  # construct stack full with the data we have
  stk.full <-  do.call(inla.stack, list(stk.e, stk.p1, stk.p2)[c( T,dp1ToF,dp2ToF)])
  
  
  # Specify formula downscaling
  formula <- y ~ 0 + Intercept + X + temp_mn + precptn + rod_scl + f(s,model = spde)+ f(s1, model = spde)
  
  # Call inla()
  res <- inla(formula, data = inla.stack.data(stk.full), 
              control.compute=list(dic=TRUE),
              control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full)))
  
  dic <- res$dic$dic
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions_p(stk.full, res, "pred1", locin_pred)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions_a(stk.full, res, "pred2", dpcontsurface, dp2)}
  
  return(list(dp1, dp2,dic))
}

fnPredictDownUK_3 <- function(depoint, dearea, dppoint = NULL, dparea = NULL, boundaryregion,
                              mesh = NULL, priorspdesigma = NULL, priorspderange = NULL){
  
  
  # Use 1 for points and 2 for areas
  # datasets estimation
  de1 <- depoint
  de2 <- dearea
  # datasets prediction
  dp1 <- dppoint
  dp2 <- dparea
  
  # Set avalue and pvalue instead of value
  de1$pvalue <- de1$value
  de2$avalue <- de2$value
  dp1$pvalue <- dp1$value
  dp2$avalue <- dp2$value  
  
  # Check inputs
  
  # check function
  # TO DO : check crs of input is crsproj
  fnCheckInputsDown(de1, de2, dp1, dp2, boundaryregion)
  
  # Logical values indicating what datasets I have
  dp1ToF <- !is.null(dp1)
  dp2ToF <- !is.null(dp2)
  
  
  # Mesh
  if (is.null(mesh) == T) {
    mesh <- fnCreateMesh(de1, boundaryregion)
  }
  
  # Create spde and index
  if(!is.null(priorspdesigma) & !is.null(priorspderange)){
    fnCheckPrior(priorspdesigma, priorspderange)
    spde <- inla.spde2.pcmatern(mesh = mesh, prior.range = priorspderange, prior.sigma = priorspdesigma)
  }else{
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = T)
  }
  indexs <- inla.spde.make.index("s", spde$n.spde)
  indexs1 <- inla.spde.make.index("s1", spde$n.spde)
  
  
  
  # Match Areal data and point data, Match dppoint and dearea
  locin <- st_join(de1, de2, left = F)
  
  ### PAULA: Keep both points and areas are inside study region. Should join by de2 or boundaryregion?
  # If TRUE, predict in points
  if(dp1ToF){
    locin_pred <- st_join(dp1, de2, left = FALSE)
  }
  
  # If TRUE, predict in areas. Construct prediction surface
  if(dp2ToF){
    # Bounding box
    bb <- unname(attributes(st_geometry(boundaryregion))$bbox)
    # Grid
    x <- seq(bb[1] - 1, bb[3] + 1, length.out = 50)
    y <- seq(bb[2] - 1, bb[4] + 1, length.out = 50)
    coop <- expand.grid(x, y)
    coop_sf <- sf::st_as_sf(coop, coords = c('Var1','Var2'), crs = st_crs(de1))
    # Keep points inside
    dpcontsurface <- coop_sf %>% st_join(boundaryregion, left = FALSE) %>% st_join(de2, left = FALSE)
  }
  
  
  
  # Create projection matrices A
  # Matrix A estimation
  Ae <- inla.spde.make.A(mesh, loc = as.matrix(st_coordinates(locin[, 1]))[, c(1, 2)])
  # Matrix A points prediction
  if(dp1ToF){Ap1 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(locin_pred)[, c(1,2)]))}
  # Matrix A areas prediction (need to predict in a continuous surface first, dense points)
  if(dp2ToF){Ap2 <- inla.spde.make.A(mesh = mesh, loc = as.matrix(st_coordinates(dpcontsurface[, 1])))}
  
  # Create stk.full
  # btilde0(x) = b0 + b0(x) (Ap1 + indexs1$s1)
  # btilde1(x) = b1 + b1(x) (Ap1 * locin_pred$avalue + indexs$s)
  # stack estimation
  
  stk.p1 = NULL
  stk.p2 = NULL
  stk.e <- inla.stack(tag = "est", data = list(y = locin$pvalue),
                      A = list(1,1, Ae * locin$avalue, Ae),
                      effects = list(data.frame(Intercept = rep(1, nrow(locin)), X = locin$avalue),temp_mn = de1$temp_mn ,
                                      s = indexs$s, s1 = indexs1$s1))
  
  
  
  if(dp1ToF){stk.p1 <- inla.stack(tag = "pred1", data = list(y = NA),
                                  A = list(1,1, Ap1 * locin_pred$avalue, Ap1),
                                  effects = list(data.frame(Intercept = rep(1, nrow(locin_pred)), X = locin_pred$avalue),temp_mn = dp1$temp_mn,
                                                  s = indexs$s, s1 = indexs1$s1))}
  # stack prediction areas
  
  # construct stack full with the data we have
  stk.full <-  do.call(inla.stack, list(stk.e, stk.p1, stk.p2)[c( T,dp1ToF,dp2ToF)])
  
  
  # Specify formula downscaling
  formula <- y ~ 0 + Intercept + X + temp_mn  + f(s,model = spde)+ f(s1, model = spde)
  
  # Call inla()
  res <- inla(formula, data = inla.stack.data(stk.full), 
              control.compute=list(dic=TRUE),
              control.predictor = list(compute = TRUE, A = inla.stack.A(stk.full)))
  
  dic <- res$dic$dic
  # Retrieve predictions points
  # Predictions points
  if(dp1ToF){dp1 <- fnRetrievePredictions_p(stk.full, res, "pred1", locin_pred)}
  # Predictions areas
  if(dp2ToF){dp2 <- fnRetrievePredictions_a(stk.full, res, "pred2", dpcontsurface, dp2)}
  
  return(list(dp1, dp2,dic))
}
