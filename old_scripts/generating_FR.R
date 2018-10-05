# Generating Semantic Analysis
data_loc = "/Users/mgavin/Desktop/geodata/"

for (f in 1:5) {
  print(f)
  load(paste("~/projects/geosemantics/inst/extdata/ch1_kd",f,".rda", sep=""))
  load(paste("~/projects/geosemantics/inst/extdata/ch1_pd",f,".rda", sep=""))
  

  
  pk = pd %*% t(kd)
  
  FR = apply(pk, 1, function(x) return(x / sum(x)))
  
  FR = t(FR)
  
  save(FR, file = paste(data_loc, "FR", f, ".rda", sep=""))
  
}
