if(input$analysis == "Consumer data"){
  df.raw <- convertToFactors(df.raw, c(input$Consumers, input$Products,
                                       input$Consumerfact))

  withProgress( message = "Calculating, please wait",
                detail = "This may take a few moments...",{ 
                  res <- tryCatch({conjoint(structure = input$struct, data = df.raw, 
                                            response = input$Response, 
                                            fixed = list(Product=input$Products, 
                                                         Consumer=input$Consumerchar), 
                                            random = input$Consumers, 
                                            facs = c(input$Consumerfact, input$Consumers,
                                                     input$Products), 
                                            corr = FALSE, alpha.random = 0.1, 
                                            alpha.fixed = 0.05)}, 
                                  error = function(e) { NULL })     
                })
}else{
  df.raw <- convertToFactors(df.raw, c(input$Assessors, input$Products, 
                                       input$Replications))
  
  
  withProgress(message = "Calculating, please wait",
               detail = "This may take a few moments...", {

                 
                 res <- tryCatch({sensmixed(input$Attributes,
                                            prod_effects=input$Products,
                                            replication = input$Replications,
                                            assessor=input$Assessors,
                                            data=df.raw,
                                            MAM = as.logical(input$MAM),
                                            product_structure = as.numeric(input$struct),
                                            error_structure = input$errstruct,
                                            control = list(
                                              calc_post_hoc = as.logical(input$calc_post_hoc),
                                              alpha.random = as.numeric(input$alpharand),
                                              alpha.fixed = as.numeric(input$alphafixed),
                                              reduce.random = as.logical(input$simplerr),
                                              keep.effs = c(unlist(strsplit(
                                                input$keep," ")),
                                                paste(paste(input$Products,
                                                          collapse = ":"),
                                                    input$Assessors, sep = ":"),
                                                input$Assessors),
                                              parallel = FALSE,
                                              MAM_mult_scaling = as.logical(input$multMAM),
                                              MAM_oneway_rand = FALSE))},
                                 error = function(e) { NULL })
                 if(as.logical(input$MAM) == TRUE){
                   res.MAM <- tryCatch({sensmixed(input$Attributes,
                                                  prod_effects=input$Products, 
                                                  replication = input$Replications,
                                                  assessor=input$Assessors, 
                                                  data=df.raw, 
                                                  control=
                                                    list(MAM_balanced = 
                                                           as.logical(input$MAM))
                                                 )}, 
                                       error = function(e) { NULL })
                 }
                 else
                   res.MAM <- NULL
                 
                 res$MAMan <- res.MAM
                   
                 
               })
}