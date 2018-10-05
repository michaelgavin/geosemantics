# S[is.na(S)] = 0
# Z[is.na(Z)] = 0
# ord = which(rowSums(S) != 0)
# S = S[ord,]
# Z = Z[ord,]

place = "borneo"

z = Z[place,]
s = S[place,]
places = place_join(place, gaz = gaz, mode = "same")
vec = place_composition(pd, places)
ord = which(vec > 0)
w = rep(1, length(s))
w[FR[place,] < mean(FR[place,])] = 0
#sem_dist = sqrt(z^2 * s^2) * w
#keywords = names(sort(sem_dist, decreasing = T)[1:100])

d = 1 - s
work = z * d
keywords = names(sort(z * d * w, decreasing = T))[1:20]
keywords2 = names(sort(z *s * w, decreasing = T))[1:20]
caption = paste("name:", place, "top word", keywords2[1])
plot(d * w, z * w, pch = 20, col = "gray", main = caption)
points(d[keywords2],(z)[keywords2], pch = 20, col = "blue")
points(d[keywords],(z)[keywords], pch = 20, col = "red")
points(d[intersect(keywords2, keywords)],(z)[intersect(keywords2, keywords)], pch = 20, col = "purple")



print("central keywords")
print(intersect(keywords2, keywords))

print("over-represented and typical")
print(setdiff(keywords2, keywords))


print("over-represented, but Atypical")
print(setdiff(keywords, keywords2))


## Get general score for conceptual work
sims = apply(S, 1, mean, na.rm = T)
freqs = apply(pd[which(rownames(pd) %in% rownames(S)),], 1, sum, na.rm = T)

work = freqs * (1 - sims)
