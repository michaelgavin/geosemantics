load("~/projects/geosemantics/inst/extdata/ch1_kd2.rda")
load("~/projects/geosemantics/inst/extdata/ch1_pd2.rda")

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
          "lade","landes","lewes","lessen","lete","milton","dictator","cher",
          "abus","sall","abo","area","dukedom","tros","lyon","divina")

pd = pd[setdiff(rownames(pd), stops),]
totals = rowSums(pd)
toponyms = names(sort(totals, decreasing = T)[1:5000])
pd = pd[toponyms,]
totals = totals[toponyms]

pk = pd %*% t(kd)
pk = ppmi(pk)
dims = 300
lsa = irlba::irlba(pk, nv = dims)
p = lsa$u %*% diag(lsa$d)
rownames(p) = toponyms
w = lsa$v %*% diag(lsa$d)
rownames(w) = rownames(kd)

zmat = pk
avgs = apply(pk, 1, mean)
devs = apply(pk, 1, sd)
for (i in 1:nrow(zmat)) {
  print(i)
  zmat[i,] = sqrt((zmat[i,] - avgs) / devs)
}
zmat[is.na(zmat)] = 0

definitive = function(place_mat, word_mat, z_mat, place) {
  s = similarity(mat = word_mat, term = place_mat[place,], output = "all")
  z = z_mat[place,]
  results = sort(s * z, decreasing = T)[1:12]
  return(results)
}

distinctiveness = function(place_mat, word_mat, z_mat, place) {
  
}

similarity = function(mat, term, output = "top") {
  cos_sim = function(x, y) {
    x %*% y/(sqrt(x %*% x) * sqrt(y %*% y))
  }
  if (length(term) == ncol(mat)) {
    vec = term
  } else {
    vec = mat[term,]
  }
  results = unlist(apply(mat, 1, cos_sim, vec))
  if (output == "top") {
    results = sort(results, decreasing = T)[1:12]
  }
  return(results)
}

# for (f in 1:5) {
#   load(paste("~/projects/geosemantics/inst/extdata/ch1_pd",f,".rda", sep=""))
#   pd = pd[setdiff(rownames(pd), stops),]
#   totals = rowSums(pd)
#   ord = which(totals > 249)
#   n = length(ord)
#   pd = pd[ord,]
#   for (i in 1:nrow(pd)) {
#     print(paste("period",f,"place",i, "of", n))
#     place = rownames(pd)[i]
#     places = place_join(place, gaz = gaz, mode = "same")
#     vec = place_composition(pd, places)
#     ord = which(vec > 0)
#     if (length(ord) > 250) {
#       toponyms = c(toponyms, place)
#     }
#   }
# }
# toponyms = sort(unique(toponyms))

# Now select words for analysis and prepare PK matrices
words = c("indians","ocean","trade","merchant","crime","england",
          "citizens","natives","people","inhabitants","men","women",
          "man","woman","child","master","servant","wife","daughter",
          "son","commodities","ships","houses","city","cities","village",
          "villages","clean","wash","serve","indies","pleasure",
          "love","loue","forgiveness","beauty","hate","anger","patience",
          "compassion","misery","bitter","hunger","meat","drink","eat", 
          "boat", "ship","horse","foot")


# This builds a blank PK for every word
pk_loc = "~/projects/PK/2/"
for (w in 1:length(words)) {
  word = words[w]
  pk = Matrix(0, length(toponyms), nrow(kd))
  rownames(pk) = toponyms
  colnames(pk) = rownames(kd)
  save(pk, file = paste(pk_loc, word, ".rda", sep = ""))
}



# Build a blank WK for each word
wk_loc = "~/projects/WK/2/"
for (p in 1:length(toponyms)) {
  place = toponyms[p]
  wk = Matrix(0, length(words), nrow(kd))
  rownames(wk) = words
  colnames(wk) = rownames(kd)
  save(wk, file = paste(wk_loc, place, ".rda", sep = ""))
}



# Now build KK matrix for each place, save data
results = list()

print("...calculating global semantic data for the period")
kd_world = ppmi(kd)
dims = 50
lsa = irlba::irlba(kd_world, nv = dims)
kd_world = lsa$u %*% diag(lsa$d)
rownames(kd_world) = rownames(kd)
kk_world = kd_world %*% t(kd_world)
for (p in 1:length(toponyms)) {
  results[[p]] = list()
  print(p)
  place = toponyms[p]
  places = place_join(place, gaz = gaz, mode = "same")
  vec = place_composition(pd, places)
  ord = which(vec > 0) # ord tells which documents mention the toponym
  if (length(ord) > 100) {
    dims = 50
  } else {
    next()
  }
  kdp = ppmi(kd[,ord])
  lsa = irlba::irlba(kdp, nv = dims)
  kds = lsa$u %*% diag(lsa$d)
  rownames(kds) = rownames(kd)
  
  
  # Fill in the WK matrix for the place
  print("...filling in word-keyword matrices with similarity scores")
  load(paste(wk_loc, place, ".rda", sep = ""))
  for (w in 1:length(words)) {
    word = words[w]
    vec = semantic_similarity(kds, word, output = "all")
    wk[word,] = vec
  }
  save(wk, file = paste(wk_loc, place, ".rda", sep = ""))
  
  # Now scroll through the PKs and add the vector
  print("... calculating kk matrix")
  kk = kds %*% t(kds)
  print("...adding place_keyword data")
  for (w in 1:length(words)) {
    word = words[w]
    load(paste(pk_loc, word, ".rda", sep = ""))
    pk[place,] = kk[word,]
    save(pk, file = paste(pk_loc, word, ".rda", sep = ""))
  }
  
  # Now calculate statistics
  print("...calculating statistics")
  sim = cosine_similarity(as.vector(kk_world), as.vector(kk))
  results[[p]]$toponym = place
  results[[p]]$total_similarity = sim
  
  sims = c()
  for (w in 1:nrow(kk)) {
    sim = cosine_similarity(kk_world[w,], kk[w,])
    sims = c(sims, sim)
  }
  names(sims) = rownames(kk)
  results[[p]]$summary = summary(sims)
  results[[p]]$similarities = sims
}

save(results, file = "~/projects/results.rda")

# Note: screwed up, only got 247 results (was using wrong toponyms value)
results = results[which(unlist(lapply(results, function(x) x$toponym)) %in% toponyms)]


setwd("~/projects/PK/2")
filenames = dir()
for (f in 1:length(filenames)) {
  print(f)
  fp = filenames[f]
  load(fp)
  if (f == 1) {
    pk_all = pk
  } else {
    pk_all = pk_all + pk
  }
}

sims = c()
vec2 = as.vector(pk_all)
for (f in 1:length(filenames)) {
  print(f)
  fp = filenames[f]
  load(fp)
  vec1 = as.vector(pk)
  sims = c(sims, cosine_similarity(vec1, vec2))
}
names(sims) = words

# Check "wife" and "england
load("~/projects/geosemantics/inst/extdata/ch1_kd2.rda")
load("~/projects/geosemantics/inst/extdata/ch1_pd2.rda")
pd = pd[toponyms,]
word = "pleasure"
ord = which(kd[word,] > 0)
kd = kd[,ord]
pd = pd[,ord]

m = matrix(0, length(toponyms), nrow(kd))
rownames(m) = toponyms
colnames(m) = rownames(kd)
word_vec = matrix(kd[word,], length(ord), 1)

dk = t(kd)
for (i in 1:length(toponyms)) {
  hits = c()
  print(i)
  place = toponyms[i]
  places = place_join(place, gaz = gaz, mode = "same")
  vec = place_composition(pd, places)
  hits = which(vec > 0)
  vec = word_vec[hits] %*% dk[hits,]
  m[i,] = as.vector(vec)
}
m = ppmi(m)

kk = rbind(kk, m)

places = place_join(place, gaz = gaz, mode = "same")
vec = place_composition(pd, places)
ord = which(vec > 0) # ord tells which documents mention the toponym


pk = pd[,ord] %*% t(kd[,ord])