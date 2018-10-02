# Geosemantic Analysis Script
PERIOD = 5
library(Matrix)

#### Load data ####
print("...loading data")
data_loc = "/data/userdata/mgavin/geodata/"
load(paste(data_loc,"ch1_kd",PERIOD,".rda", sep=""))
load(paste(data_loc,"ch1_pd",PERIOD,".rda", sep=""))
load(paste(data_loc,"ch1_gazetteer.rda", sep=""))

#### Define functions ####
cosine_similarity = function(x, y) {
  x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
}

place_composition = function(mat, places) {
  places = places[places %in% rownames(mat)]
  if (length(places) > 1) {
    mat = mat[places,]
    mat = as.matrix(mat)
    vec = colSums(mat)
  } else {
    vec = mat[places,]
  }
  return(vec)
}

place_join = function(places, gaz, mode = "all") {
  places = c(places,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
  places = c(places,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$SUBJECT %in% places))])
  places = c(places,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
  places = unique(places)
  
  if (mode == "all") {
    places = c(places, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
    places = c(places, gaz$SUBJECT[intersect(grep("is in", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
    places = c(places,gaz$SUBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$OBJECT %in% places))])
    places = c(places,gaz$OBJECT[intersect(grep("same", gaz$PREDICATE), which(gaz$SUBJECT %in% places))])
    places = unique(places)
  }
  return(places)
}

ppmi = function(mat) {
  total = sum(mat, na.rm = T)
  pcol = apply(mat, 2, function(x) sum(x, na.rm = T) / total)
  prow = apply(mat, 1, function(x) sum(x, na.rm = T) / total)
  pmat = as.matrix(prow) %*% t(pcol)
  mat = mat / total
  mat = mat / pmat
  mat = apply(mat, 1, log)
  mat = t(mat)
  mat[mat < 0] = 0
  mat[is.na(mat)] = 0
  return(mat)
}



#### Initializing Data ####
stops = c("now","called","soule","plain","losse","reading","deal","burning",
          "coun","ile","key","elo","purgatory","nice","wormes","cave","arg",
          "evi","trin","greg","tenet","iland","countrie","hils","dean",
          "continent","pool","thorn","farewel","dive","aur","vvord","split",
          "comb","ain","agag","aver","barking","batter","bein","boot","brest",
          "chimera","cistern","crouch","dole","dina","earne","flint","ford",
          "grays","gulph","hil","hag","hon","horne","inn","julia","leo","ls",
          "martha","marche","maii","novum","novem","noto","pike",
          "pest","orna","pontius","portae","riuer","pris","py","rick",
          "saynt","scole","sug","tos","vend","zea","zele","park","nied",
          "sour","tende","sever","rian","rase","mur","lade","dou","compendium",
          "calvary","cae","caves","capitol","breste","bres","cowes","chang",
          "lade","landes","lewes","lessen","lete","milton","jacobites","acci",
          "vvorld","andrews","airy","chappel","abus","abo","archangel",
          "assur","vils","lyon","dai","cher","widen","curi","pau","arras",
          "tring","mayn","nea","burg","area","dur","foy","lure","pene","dents",
          "britains","englands","wick","craven","dour","dort","sera","mora",
          "dictator","shap","sall","dium","brower","rodrigo","januarius",
          "valladolid","flour","cleri","chu","peak","stella","fens","chard",
          "dukedom","wells")
pd = pd[setdiff(rownames(pd), stops),]
totals = rowSums(pd)
toponyms = names(sort(totals[totals > 1], decreasing = T))


print("...getting keyword stats")
word_totals = colSums(kd)
normed_kd = apply(kd, 1, function(x, sums) {x / sums}, sums = word_totals)
normed_kd = t(normed_kd)
avgs = apply(normed_kd, 1, mean, na.rm = T)
devs = apply(normed_kd, 1, sd, na.rm = T)


print("...initializing S and Z")
Z = matrix(0, nrow(pd), nrow(kd))
rownames(Z) = rownames(pd)
colnames(Z) = rownames(kd)

S = matrix(0, nrow(pd), nrow(kd))
rownames(S) = rownames(pd)
colnames(S) = rownames(kd)

sname = paste("S",PERIOD,".rda", sep = "")
zname = paste("Z",PERIOD,".rda", sep = "")

print("...calculating KK")
KK = kd %*% t(kd)
KK = ppmi(KK)
# dims = 50
# lsa = irlba::irlba(KK, nv = dims)
# KK = lsa$u %*% diag(lsa$d)
# rownames(KK) = rownames(kd)

#### Performing Semantic Analysis ####
("...calculating similarities")
for (i in 1:length(toponyms)) {
  place = toponyms[i]
  places = place_join(place, gaz = gaz, mode = "same")
  vec = place_composition(pd, places)
  ord = which(vec > 0)
  
  # Define Z
  word_totals = colSums(kd[,ord])
  normed_kd = apply(kd[,ord], 1, function(x, sums) {x / sums}, sums = word_totals)
  normed_kd = t(normed_kd)
  pavgs = apply(normed_kd, 1, mean, na.rm = T)
  z = (pavgs - avgs) / devs
  Z[place,] = z
  
  # Define kk
  kk = kd[,ord] %*% t(kd[,ord])
  kk = ppmi(kk)
  # lsa = irlba::irlba(kk, nv = dims)
  # kk = lsa$u %*% diag(lsa$d)
  # rownames(kk) = rownames(kd)
  
  sims = c()
  for (s in 1:nrow(kk)) {
    sim = cosine_similarity(kk[s,], KK[s,])
    sims = c(sims, sim)
  }
  names(sims) = rownames(kk)
  
  S[place,] = sims
  
  if (i %% 10 == 0) {
    print(i)
    print("...saving progress")
    save(S, file = paste(data_loc,sname, sep=""))
    save(Z, file = paste(data_loc,zname, sep=""))
  }
}

