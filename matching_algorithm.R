# This file demonstrates the matching algorithm used in article "Hospital readmissions among homeless people: a cohort study of linked hospitalisation and mortality data in England for 2,772 homeless inpatients"
# [Publication info]

library(data.table)
library(epitools) # for confidence interval of rate ratio (function 'rateratio')
set.seed(4)

#----------------------
# generate example data
#----------------------

n1 <- 750 # number of homeless patients
n2 <- 5000 # number of housed patients

# age, sex and hospital. homeless patients are younger and more likely to be male

homeless <- data.table(grp = 'homeless', age = exp(rnorm(n1, 3.3, 0.2)), sex = sample(c('m', 'f'), n1, replace = T, prob = c(0.7, 0.3)))
housed <- data.table(grp = 'housed', age = runif(n2, 20, 90), sex = sample(c('m', 'f'), n2, replace = T))
all.pt <- rbind(homeless, housed)
age_labels <- paste0(seq(15, 85, 10), '-', seq(25, 95, 10))
all.pt[, age_group := findInterval(age, seq(15, 95, 10))]
all.pt[, age_group := factor(age_group, 1:8, age_labels)]
all.pt[, hospital := sample(LETTERS[1:5], .N, replace = T, prob = 1:5)]
all.pt[, id := .I] # unique patient ID

# create follow-up data (with higher rates for homeless (2x), older, and female patients)

all.pt[, follow.up := rnorm(.N, 1000, 100)]
age.scale <- all.pt$age - min(all.pt$age)
age.scale <- age.scale / max(age.scale) + 1
predicted_rate <- rnorm(n1+n2, 4, 0.75) * ((all.pt$grp == 'homeless') + 1) * age.scale * ifelse(all.pt$sex == 'f', 2, 1)
all.pt[, readmissions := rpois(.N, predicted_rate)]

#--------------------------
# identify possible matches
#--------------------------

mvars <- c('age_group', 'sex', 'hospital') # variables on which to match
matches <- all.pt[grp == 'housed'][all.pt[grp == 'homeless'], on = mvars, allow.cartesian = T]
setnames(matches, c('id', 'i.id', 'age', 'i.age'), c('housed.id', 'homeless.id', 'housed.age', 'homeless.age'))

#-------------------------
# select 1:1 exact matches
#-------------------------

# prints off the number of matches at intervals of 100 matches. 'mtr' is the dataset of successful matches

mtr <- NULL; nl <- 0
while (nrow(matches) > 0) {
  r <- matches[sample(.N, 1)]
  matches <- matches[(housed.id != r$housed.id) & (homeless.id != r$homeless.id)] # remove the individuals from the potential matches
  mtr <- rbind(mtr, r)
  nl <- nl + 1; if (nl/100 == floor(nl/100)) print(nl)
}
mtr[, pair := .I] # id of pairing for paired analysis
all.pt[, matched := id %in% c(mtr$housed.id, mtr$homeless.id)]
pair.ids <- data.table(id = c(mtr$housed.id, mtr$homeless.id), pair = mtr$pair)
all.pt <- pair.ids[all.pt, on = 'id']

# check balance

check.bal <- dcast(all.pt[matched == T], paste0(paste0(mvars, collapse = '+'), '~grp'), value.var = 'grp', fun.aggregate = length)
table(check.bal$homeless == check.bal$housed) # should be all TRUE

#------------------------------------------
# compare rate ratios (homeless vs. housed)
#------------------------------------------

all.pt[, grp := factor(grp, c('housed', 'homeless'))]

# in the whole dataset, homeless patients have a higher rate of readmission
whole_dataset <- all.pt[, .(n = .N, follow.up = sum(follow.up), readmissions = sum(readmissions), rate = sum(readmissions) / sum(follow.up) * 1000), grp]
rateratio(x = whole_dataset$readmissions[2:1], y = whole_dataset$follow.up[2:1])$measure

# compare with crude IRR from poisson regression, using whole dataset
m1 <- glm(readmissions ~ grp + offset(log(follow.up)), all.pt, family = 'poisson')
exp(cbind(m1$coef, confint(m1))) # same

# in the matched dataset, the difference is wider because housed patients are older and a greater proportion are female
matched_only <- all.pt[matched == T, .(n = .N, follow.up = sum(follow.up), readmissions = sum(readmissions), rate = sum(readmissions) / sum(follow.up) * 1000), grp]
rateratio(x = matched_only$readmissions[2:1], y = matched_only$follow.up[2:1])$measure

# compare with adjusted IRR from poisson regression, using the whole dataset
m2 <- glm(readmissions ~ grp + age + sex + hospital + offset(log(follow.up)), data = all.pt, family = 'poisson')
exp(cbind(m2$coef, confint(m2))) # similar
