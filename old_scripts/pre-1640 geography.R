library(xml2)
library(Matrix)

xml_loc = "~/Documents/eebo_xml/"
fnames = dir(xml_loc)


cleanup = function(fp) {
  txt = read_xml(fp)
  if (length(grep("headed", fp)) == 1) {
    txt = xml_find_first(txt, ".//TEXT")
  } else {
    txt = xml_find_first(txt,".//d1:text")
  }
  txt = xml_text(txt)
  txt = gsub("ſ", "s", txt)
  txt = gsub("[0-9]", "", txt)
  txt = gsub("vv", "w", txt)
  txt = gsub("'d ", "ed ", txt)
  txt = gsub("'ring ", "ering ", txt)
  txt = tolower(txt)
  txt = unlist(strsplit(txt, split = "\\W"))
  txt = txt[txt!=""]
  return(txt)
}

simple_clean = function(txt) {
  txt = xml_text(txt)
  txt = gsub("ſ", "s", txt)
  #txt = gsub("??", "s", txt)
  txt = gsub("[0-9]", "", txt)
  txt = gsub("vv", "w", txt)
  txt = gsub("'d ", "ed ", txt)
  txt = gsub("'ring ", "ering ", txt)
  txt = tolower(txt)
  txt = gsub(" sea", "sea", txt)
  txt = unlist(strsplit(txt, split = "\\W"))
  txt = txt[txt!=""]
  return(txt)
}


for (f in 1:length(ord_sample)) {
  print(f)
  fp = fnames[ord_sample][f]
  txt = cleanup(fp)
  if (f == ord_sample[1]) {
    vocab = table(txt)
  } else {
    txt = table(txt)
    new_words = setdiff(names(txt), names(vocab))
    vocab = c(vocab,txt[new_words])
    words = intersect(names(txt), names(vocab))
    vocab[words] = vocab[words] + txt[words]
  }
  if (f %% 50 == 0) {
    vocab = vocab[vocab > 2]
  }
}

pnames = intersect(gaz$SUBJECT, names(vocab))
keywords = c(names(vocab)[vocab > 1000], pnames)
keywords = sort(unique(keywords))

m = Matrix(0, nrow = length(keywords), ncol = length(ord),
           dimnames = list(keywords,mdata[ord,1]) )

for (j in 1:ncol(m)) {
  
  print(j)
  
  fp = fnames[ord][j]
  txt = cleanup(fp)
  
  txt = txt[txt %in% rownames(m)]
  
  freqs = table(txt)
  
  vec = rep(0, length(rownames(m)))
  names(vec) = rownames(m)
  
  vec[names(freqs)] = freqs
  
  m[,j] = vec

}

setwd(xml_loc)
pnames = unique(gaz$SUBJECT)
place_mat = matrix(0, length(pnames), 100)
rownames(place_mat) = pnames
for (j in 1:100) {
  print(j)
  fp = fnames[j]
  txt = cleanup(fp)
  txt = txt[txt %in% pnames]
  freqs = table(txt)
  place_mat[names(freqs),j] = freqs
}
save(place_mat, file = "~/Downloads/older_place_mat.rda")



# Working with A17900
setwd(xml_loc)
txt = read_xml("A19700.xml")
regions = xml_find_all(txt, ".//d1:div[@type='region']")
SUBJECT = c()
OBJECT = c()
for (i in 1:length(regions)) {
  region = tolower(xml_attrs(regions[[i]])[1])
  divs = xml_find_all(regions[[i]],".//d1:div")
  divs = divs[2:length(divs)]
  
  for(j in 1:length(divs)) {
    div = divs[[j]]
    header = xml_find_first(div, ".//d1:head")
    header = simple_clean(header)
    header = paste(header, collapse = " ")
    SUBJECT = c(SUBJECT, header)
    OBJECT = c(OBJECT, region)
    rows = xml_find_all(div, ".//d1:row")
    rows = simple_clean(rows)
    SUBJECT = c(SUBJECT, rows)
    OBJECT = c(OBJECT, rep(header, length(rows)))
  }
}
PREDICATE = rep("is in", length(SUBJECT))
df = data.frame(SUBJECT, PREDICATE, OBJECT)
stopwords = c("citie", "an", "de", "of", "pole", "south", "the", "s")
ord = which(SUBJECT %in% stopwords == F)
View(df[ord,])
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTk4NjkwNjE3N119
-->