
require(SensMixed)
load(system.file("testdata","sensBO.RData",package="SensMixed"))


testBO <- FALSE

if(testBO)
  result_MAM <- sensmixed(names(originalBOdadta)[8:10],
                       prod_effects=c("Track", "Car", "SPL"),
                       replication="RepetitionFixed", assessor="Participant", 
                      product_structure=3, data=originalBOdadta, 
                      parallel = FALSE, error_structure = "3-WAY")


