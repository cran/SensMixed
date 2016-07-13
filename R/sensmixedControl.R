## DOC: ../man/sensmixedControl.Rd
sensmixedControl <- function(MAM_mult_scaling = FALSE, MAM_oneway_rand = FALSE,
                             MAM_balanced = FALSE, MAM_adjusted = FALSE, 
                             MAM_alpha_conditional = 1,
                             calc_post_hoc = FALSE, parallel = FALSE, 
                             reduce.random=TRUE, alpha.random = 0.1, 
                             alpha.fixed = 0.05, interact.symbol = ":", 
                             keep.effs = NULL)
  {
    structure(namedList(mamControl =
                          namedList(MAM_mult_scaling, 
                                    MAM_oneway_rand,
                                    MAM_balanced, 
                                    MAM_adjusted, 
                                    MAM_alpha_conditional),
                        calc_post_hoc, parallel, 
                        reduce.random, alpha.random, 
                        alpha.fixed, interact.symbol, 
                        keep.effs),
              class = c("sensmixedControl"))
  }