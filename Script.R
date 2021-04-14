###### Script Aegean Endemics 08-04-2021 ########
library(tidyverse)
library(Taxonstand)
library(V.PhyloMaker)
library(ape)
library(Hmisc)

setwd("G:/Mon Drive/Aegean endemics/Final")

############### Dataset preparation ########################################
data <- read.csv("AegeanData.csv", sep = ";") %>% 
  mutate(Couple = as.factor(ID), 
         Status = as.factor(Status)) %>%
  select(Couple, Status, TaxonName, SVPetLenMin_mm, SVPetLenMax_mm, SepLenMin_mm, SepLenMax_mm, 
         FrLenMin_mm, FrLenMax_mm, SeedLenMin_mm, SeedLenMax_mm, SVPropLenMin_mm, SVPropLenMax_mm, 
         LeafLenMin_mm, LeafLenMax_mm, LeafWidMin_mm, LeafWidMax_mm, StemLenMin_cm, StemLenMax_cm, 
         MinAlt, MaxAlt, AltRange, StartFloPer, EndFloPer, FloPer) %>%
  # Add average value for traits
  rowwise() %>%
  mutate(PetLen = mean(c(SVPetLenMin_mm, SVPetLenMax_mm), na.rm = T),
         SepLen = mean(c(SepLenMin_mm, SepLenMax_mm), na.rm = T), 
         FrLen = mean(c(FrLenMin_mm, FrLenMax_mm), na.rm = T),
         SeedLen = mean(c(SeedLenMin_mm, SeedLenMax_mm), na.rm = T),
         PropLen = mean(c(SVPropLenMin_mm, SVPropLenMax_mm), na.rm = T),
         StemLen = mean(c(StemLenMin_cm, StemLenMax_cm), na.rm = T),
         LeafLen = mean(c(LeafLenMin_mm,LeafLenMax_mm), na.rm = T),
         LeafWid = mean(c(LeafWidMin_mm, LeafWidMax_mm), na.rm = T)) %>%
  ungroup()

## Keep values if present for both species in a pair
data.fin <- data[,1:3]
for(i in 4:ncol(data)){
    d <- data[,c(1:3,i)] %>%
    na.omit() %>%
    group_by(Couple) %>%
    count() %>%
    filter(n == 2)
  d.tmp <- data[data$Couple %in% d$Couple,]
  data.fin <- data.fin %>%
    full_join(d.tmp[,c(1:3,i)])
}

##### Create phylogeny 
## Phylogeny 
tpl <- TPL(data.fin$TaxonName) # Used to extract family
sp.list <- tpl %>%
  select(Species = Taxon, Genus, Family) %>%# Species names were checked and do not follow the plant list.
  mutate(Species = str_replace_all(Species, " ", "_"))

sp.list$Family[sp.list$Family == "Leguminosae"] <- "Fabaceae"
sp.list$Genus[sp.list$Species == "Anthemis_ammanthus_subsp_ammanthus"] <- "Anthemis"
sp.list$Genus[sp.list$Species == "Anthemis_palestina"] <- "Anthemis"
sp.list$Genus[sp.list$Species == "Onopordum_bracteatum_subsp_bracteatum"] <- "Onopordum"
sp.list$Family[sp.list$Genus == "Onopordum"] <- "Compositae"
sp.list$Family[sp.list$Genus == "Anthemis"] <- "Compositae"
sp.list$Genus[sp.list$Species == "Scorzonera_cana"] <- "Scorzonera"
sp.list$Family[sp.list$Genus == "Scorzonera"] <- "Compositae"
sp.list$Genus[sp.list$Species == "Rorippa_thracica"] <- "Rorippa"
sp.list$Family[sp.list$Genus == "Rorippa"] <- "Brassicaceae"
sp.list$Genus[sp.list$Species == "Arenaria_leptoclados"] <- "Arenaria"
sp.list$Family[sp.list$Genus == "Arenaria"] <- "Caryophyllaceae"
sp.list$Genus[sp.list$Species == "Cerastium_decalvans_subsp_glutinosum"] <- "Cerastium"
sp.list$Family[sp.list$Genus == "Cerastium"] <- "Caryophyllaceae"
sp.list$Genus[sp.list$Species == "Minuartia_attica_subsp_idaea"] <- "Minuartia"
sp.list$Family[sp.list$Genus == "Minuartia"] <- "Caryophyllaceae"
sp.list$Genus[sp.list$Species == "Linum_gyaricum_subsp_icaricum"] <- "Linum"
sp.list$Family[sp.list$Genus == "Linum"] <- "Linaceae"
sp.list$Genus[sp.list$Species == "Limonium_sieberi"] <- "Limonium"
sp.list$Family[sp.list$Genus == "Limonium"] <- "Plumbaginaceae"
sp.list$Genus[sp.list$Species == "Aegilops_biuncialis_subsp_archipelagica"] <- "Aegilops"
sp.list$Genus[sp.list$Species == "Aegilops_biuncialis_subsp_biuncialis"] <- "Aegilops"
sp.list$Family[sp.list$Genus == "Aegilops"] <- "Poaceae"
sp.list$Genus[sp.list$Species == "Cyclamen_rhodium_subsp_rhodium"] <- "Cyclamen"
sp.list$Family[sp.list$Genus == "Cyclamen"] <- "Primulaceae"
sp.list$Family[sp.list$Family == "Compositae"] <- "Asteraceae"

data <- cbind(sp.list, data.fin) %>%
  tibble() %>%
  select(-TaxonName) %>%
  mutate(Species = as.factor(Species), 
         Genus = as.factor(Genus),
         Family = as.factor(Family),
         Couple = as.factor(Couple)) %>%
  select(Couple, everything())
write.csv(data, "data.csv")

# Get phylogenetic tree 
tree <- phylo.maker(data[,2:4])$scenario.3
plot(tree, cex = 0.5)
write.tree(tree, "tree.tre")


############### Load data  #####
data.r <- read.csv("data.csv", sep = ",") %>%
  tibble() %>%
  select(-X) %>%
  mutate(Couple = as.factor(Couple), 
         Species = as.factor(Species), 
         Genus = as.factor(Genus), 
         Family = as.factor(Family), 
         Status = as.factor(Status))
tree <- read.tree("tree.tre")
colSums(!is.na(data.r))/2

## Box-Cox transform data for analyses 
bc <- function(x) bestNormalize::boxcox(x , standardize = TRUE)$x.t #BoxCox tranform
data <- data.r %>%
  mutate(MinAlt = MinAlt +1, AltRange = AltRange +1) %>% ## must be strictly positive to boxcox transform
  mutate_if(is.numeric, bc)

### correlation between traits ##########
data.min.max <- data %>%
  select(SVPetLenMin_mm, SVPetLenMax_mm, SepLenMin_mm,
         SepLenMax_mm, FrLenMin_mm, FrLenMax_mm, SeedLenMin_mm, SeedLenMax_mm,
         SVPropLenMin_mm, SVPropLenMax_mm, LeafLenMin_mm, LeafLenMax_mm, 
         LeafWidMin_mm, LeafWidMax_mm, StemLenMin_cm, StemLenMax_cm, MinAlt, 
         MaxAlt, AltRange, StartFloPer, EndFloPer, FloPer)
cor.min.max <- round(cor(data.min.max, method = "pearson", use = "complete.obs"), 2)

data.cor <- data %>%
  select(PetLen, SepLen, PropLen, StemLen, LeafLen, LeafWid, 
         MinAlt, MaxAlt, StartFloPer, EndFloPer)

CorTr <- rcorr(as.matrix(data.cor), type = "pearson")
#write.csv(round(CorTr$r, 2), "r.csv")
#write.csv(round(CorTr$P, 3), "p.csv")
###################################
######## Models ###################

## Model parameters
delta = 0.99
treedepth = 15
iter = 6000
warmup = 2000
chains = 4
seed = 173

#################### Models ##################################################
###### PetLen
d <- data %>% select(Couple, Species, Status, PetLen) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = PetLen) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE)) # Comparison number of larger E vs W
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
PetLen <- brm(formula = PetLen ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
              data = d, data2 = list(phy.mat = phy.mat), family = student(),
              prior = c(prior(normal(0, 1), class = "Intercept"),
                        prior(exponential(1), class = "sigma")),
              control = list(adapt_delta = delta, max_treedepth = treedepth),
              iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
              file = "output/PetLen")

###### SepLen
d <- data %>% select(Couple, Species, Status, SepLen) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = SepLen) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SepLen <- brm(formula = SepLen ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
              data = d, data2 = list(phy.mat = phy.mat), family = student(),
              prior = c(prior(normal(0, 1), class = "Intercept"),
                        prior(exponential(1), class = "sigma")),
              control = list(adapt_delta = delta, max_treedepth = treedepth),
              iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
              file = "output/SepLen")

###### PropLen
d <- data %>% select(Couple, Species, Status, PropLen) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = PropLen) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
PropLen <- brm(formula = PropLen ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
               data = d, data2 = list(phy.mat = phy.mat), family = student(),
               prior = c(prior(normal(0, 1), class = "Intercept"),
                         prior(exponential(1), class = "sigma")),
               control = list(adapt_delta = delta, max_treedepth = treedepth),
               iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
               file = "output/PropLen")

###### StemLen
d <- data %>% select(Couple, Species, Status, StemLen) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = StemLen) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
StemLen <- brm(formula = StemLen ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
               data = d, data2 = list(phy.mat = phy.mat), family = student(),
               prior = c(prior(normal(0, 1), class = "Intercept"),
                         prior(exponential(1), class = "sigma")),
               control = list(adapt_delta = delta, max_treedepth = treedepth),
               iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
               file = "output/StemLen")

###### LeafLen
d <- data %>% select(Couple, Species, Status, LeafLen) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = LeafLen) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
LeafLen <- brm(formula = LeafLen ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
               data = d, data2 = list(phy.mat = phy.mat), family = student(),
               prior = c(prior(normal(0, 1), class = "Intercept"),
                         prior(exponential(1), class = "sigma")),
               control = list(adapt_delta = delta, max_treedepth = treedepth),
               iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
               file = "output/LeafLen")

###### LeafWid
d <- data %>% select(Couple, Species, Status, LeafWid) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
LeafWid <- brm(formula = LeafWid ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
               data = d, data2 = list(phy.mat = phy.mat), family = student(),
               prior = c(prior(normal(0, 1), class = "Intercept"),
                         prior(exponential(1), class = "sigma")),
               control = list(adapt_delta = delta, max_treedepth = treedepth),
               iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
               file = "output/LeafWid")

###### MinAlt
d <- data %>% select(Couple, Species, Status, MinAlt) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = MinAlt) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = MinAlt) %>% summarise_all(mean)
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
MinAlt <- brm(formula = MinAlt ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
              data = d, data2 = list(phy.mat = phy.mat), family = student(),
              prior = c(prior(normal(0, 1), class = "Intercept"),
                        prior(exponential(1), class = "sigma")),
              control = list(adapt_delta = delta, max_treedepth = treedepth),
              iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
              file = "output/MinAlt")
###### MaxAlt
d <- data %>% select(Couple, Species, Status, MaxAlt) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = MaxAlt) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = MaxAlt) %>% summarise_all(mean)
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
MaxAlt <- brm(formula = MaxAlt ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
              data = d, data2 = list(phy.mat = phy.mat), family = student(),
              prior = c(prior(normal(0, 1), class = "Intercept"),
                        prior(exponential(1), class = "sigma")),
              control = list(adapt_delta = delta, max_treedepth = treedepth),
              iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
              file = "output/MaxAlt")

###### StartFloPer
d <- data %>% select(Couple, Species, Status, StartFloPer) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = StartFloPer) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = StartFloPer) %>% summarise_all(mean)
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
StartFloPer <- brm(formula = StartFloPer ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                   data = d, data2 = list(phy.mat = phy.mat), family = student(),
                   prior = c(prior(normal(0, 1), class = "Intercept"),
                             prior(exponential(1), class = "sigma")),
                   control = list(adapt_delta = delta, max_treedepth = treedepth),
                   iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                   file = "output/StartFloPer")
###### EndFloPer
d <- data %>% select(Couple, Species, Status, EndFloPer) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = EndFloPer) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = EndFloPer) %>% summarise_all(mean)
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
EndFloPer <- brm(formula = EndFloPer ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                 data = d, data2 = list(phy.mat = phy.mat), family = student(),
                 prior = c(prior(normal(0, 1), class = "Intercept"),
                           prior(exponential(1), class = "sigma")),
                 control = list(adapt_delta = delta, max_treedepth = treedepth),
                 iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                 file = "output/EndFloPer")


##### Supplementary raw traits #############################
####### SVPetLenMin_mm
d <- data %>% select(Couple, Species, Status, SVPetLenMin_mm) %>% na.omit()# Keep species for which trait is known
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SVPetLenMin_mm <- brm(formula = SVPetLenMin_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                      data = d, 
                      data2 = list(phy.mat = phy.mat),
                      family = student(),
                      prior = c(prior(normal(0, 1), class = "Intercept"),
                                prior(exponential(1), class = "sigma")),
                      control = list(adapt_delta = delta, max_treedepth = treedepth),
                      iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                      file = "output/SVPetLenMin_mm")
###### SVPetLenMax_mm
d <- data %>% select(Couple, Species, Status, SVPetLenMax_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SVPetLenMax_mm <- brm(formula = SVPetLenMax_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                      data = d, data2 = list(phy.mat = phy.mat), family = student(),
                      prior = c(prior(normal(0, 1), class = "Intercept"),
                                prior(exponential(1), class = "sigma")),
                      control = list(adapt_delta = delta, max_treedepth = treedepth),
                      iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                      file = "output/SVPetLenMax_mm")
###### SepLenMin_mm
d <- data %>% select(Couple, Species, Status, SepLenMin_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SepLenMin_mm <- brm(formula = SepLenMin_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                    data = d, data2 = list(phy.mat = phy.mat), family = student(),
                    prior = c(prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(1), class = "sigma")),
                    control = list(adapt_delta = delta, max_treedepth = treedepth),
                    iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                    file = "output/SepLenMin_mm")
###### SepLenMax_mm
d <- data %>% select(Couple, Species, Status, SepLenMax_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SepLenMax_mm <- brm(formula = SepLenMax_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                    data = d, data2 = list(phy.mat = phy.mat), family = student(),
                    prior = c(prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(1), class = "sigma")),
                    control = list(adapt_delta = delta, max_treedepth = treedepth),
                    iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                    file = "output/SepLenMax_mm")
###### FrLenMin_mm
d <- data %>% select(Couple, Species, Status, FrLenMin_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
FrLenMin_mm <- brm(formula = FrLenMin_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                   data = d, data2 = list(phy.mat = phy.mat), family = student(),
                   prior = c(prior(normal(0, 1), class = "Intercept"),
                             prior(exponential(1), class = "sigma")),
                   control = list(adapt_delta = delta, max_treedepth = treedepth),
                   iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                   file = "output/FrLenMin_mm")
###### FrLenMax_mm
d <- data %>% select(Couple, Species, Status, FrLenMax_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
FrLenMax_mm <- brm(formula = FrLenMax_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                   data = d, data2 = list(phy.mat = phy.mat), family = student(),
                   prior = c(prior(normal(0, 1), class = "Intercept"),
                             prior(exponential(1), class = "sigma")),
                   control = list(adapt_delta = delta, max_treedepth = treedepth),
                   iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                   file = "output/FrLenMax_mm")
###### SeedLenMin_mm
d <- data %>% select(Couple, Species, Status, SeedLenMin_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SeedLenMin_mm <- brm(formula = SeedLenMin_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                     data = d, data2 = list(phy.mat = phy.mat), family = student(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     control = list(adapt_delta = delta, max_treedepth = treedepth),
                     iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                     file = "output/SeedLenMin_mm")
###### SeedLenMax_mm
d <- data %>% select(Couple, Species, Status, SeedLenMax_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SeedLenMax_mm <- brm(formula = SeedLenMax_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                     data = d, data2 = list(phy.mat = phy.mat), family = student(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     control = list(adapt_delta = delta, max_treedepth = treedepth),
                     iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                     file = "output/SeedLenMax_mm")
###### SVPropLenMin_mm
d <- data %>% select(Couple, Species, Status, SVPropLenMin_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SVPropLenMin_mm <- brm(formula = SVPropLenMin_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                       data = d, data2 = list(phy.mat = phy.mat), family = student(),
                       prior = c(prior(normal(0, 1), class = "Intercept"),
                                 prior(exponential(1), class = "sigma")),
                       control = list(adapt_delta = delta, max_treedepth = treedepth),
                       iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                       file = "output/SVPropLenMin_mm")
###### SVPropLenMax_mm
d <- data %>% select(Couple, Species, Status, SVPropLenMax_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SVPropLenMax_mm <- brm(formula = SVPropLenMax_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                       data = d, data2 = list(phy.mat = phy.mat), family = student(),
                       prior = c(prior(normal(0, 1), class = "Intercept"),
                                 prior(exponential(1), class = "sigma")),
                       control = list(adapt_delta = delta, max_treedepth = treedepth),
                       iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                       file = "output/SVPropLenMax_mm")
###### LeafLenMin_mm
d <- data %>% select(Couple, Species, Status, LeafLenMin_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
LeafLenMin_mm <- brm(formula = LeafLenMin_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                     data = d, data2 = list(phy.mat = phy.mat), family = student(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     control = list(adapt_delta = delta, max_treedepth = treedepth),
                     iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                     file = "output/LeafLenMin_mm")
###### LeafLenMax_mm
d <- data %>% select(Couple, Species, Status, LeafLenMax_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
LeafLenMax_mm <- brm(formula = LeafLenMax_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                     data = d, data2 = list(phy.mat = phy.mat), family = student(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     control = list(adapt_delta = delta, max_treedepth = treedepth),
                     iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                     file = "output/LeafLenMax_mm")
###### LeafWidMin_mm
d <- data %>% select(Couple, Species, Status, LeafWidMin_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
LeafWidMin_mm <- brm(formula = LeafWidMin_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                     data = d, data2 = list(phy.mat = phy.mat), family = student(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     control = list(adapt_delta = delta, max_treedepth = treedepth),
                     iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                     file = "output/LeafWidMin_mm")
###### LeafWidMax_mm
d <- data %>% select(Couple, Species, Status, LeafWidMax_mm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
LeafWidMax_mm <- brm(formula = LeafWidMax_mm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                     data = d, data2 = list(phy.mat = phy.mat), family = student(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     control = list(adapt_delta = delta, max_treedepth = treedepth),
                     iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                     file = "output/LeafWidMax_mm")
###### StemLenMin_cm
d <- data %>% select(Couple, Species, Status, StemLenMin_cm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
StemLenMin_cm <- brm(formula = StemLenMin_cm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                     data = d, data2 = list(phy.mat = phy.mat), family = student(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     control = list(adapt_delta = delta, max_treedepth = treedepth),
                     iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                     file = "output/StemLenMin_cm")
###### StemLenMax_cm
d <- data %>% select(Couple, Species, Status, StemLenMax_cm) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
StemLenMax_cm <- brm(formula = StemLenMax_cm ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                     data = d, data2 = list(phy.mat = phy.mat), family = student(),
                     prior = c(prior(normal(0, 1), class = "Intercept"),
                               prior(exponential(1), class = "sigma")),
                     control = list(adapt_delta = delta, max_treedepth = treedepth),
                     iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                     file = "output/StemLenMax_cm")

###### AltRange
d <- data %>% select(Couple, Species, Status, AltRange) %>% na.omit()
d %>% select (-Species) %>% pivot_wider(names_from = Status, values_from = AltRange) %>% 
  mutate(EW = E > W, EeW = E ==W) %>% summarise("E>W" = sum(EW == TRUE), "W>E" = sum(EW == FALSE & EeW == FALSE))
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
AltRange <- brm(formula = AltRange ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                data = d, data2 = list(phy.mat = phy.mat), family = student(),
                prior = c(prior(normal(0, 1), class = "Intercept"),
                          prior(exponential(1), class = "sigma")),
                control = list(adapt_delta = delta, max_treedepth = treedepth),
                iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                file = "output/AltRange")
###### StartFloPer
d <- data %>% select(Couple, Species, Status, StartFloPer) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
StartFloPer <- brm(formula = StartFloPer ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                   data = d, data2 = list(phy.mat = phy.mat), family = student(),
                   prior = c(prior(normal(0, 1), class = "Intercept"),
                             prior(exponential(1), class = "sigma")),
                   control = list(adapt_delta = delta, max_treedepth = treedepth),
                   iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                   file = "output/StartFloPer")
###### EndFloPer
d <- data %>% select(Couple, Species, Status, EndFloPer) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
EndFloPer <- brm(formula = EndFloPer ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
                 data = d, data2 = list(phy.mat = phy.mat), family = student(),
                 prior = c(prior(normal(0, 1), class = "Intercept"),
                           prior(exponential(1), class = "sigma")),
                 control = list(adapt_delta = delta, max_treedepth = treedepth),
                 iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
                 file = "output/EndFloPer")
###### FloPer
d <- data %>% select(Couple, Species, Status, FloPer) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
FloPer <- brm(formula = FloPer ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
              data = d, data2 = list(phy.mat = phy.mat), family = student(),
              prior = c(prior(normal(0, 1), class = "Intercept"),
                        prior(exponential(1), class = "sigma")),
              control = list(adapt_delta = delta, max_treedepth = treedepth),
              iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
              file = "output/FloPer")
###### FrLen
d <- data %>% select(Couple, Species, Status, FrLen) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
FrLen <- brm(formula = FrLen ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
             data = d, data2 = list(phy.mat = phy.mat), family = student(),
             prior = c(prior(normal(0, 1), class = "Intercept"),
                       prior(exponential(1), class = "sigma")),
             control = list(adapt_delta = delta, max_treedepth = treedepth),
             iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
             file = "output/FrLen")
###### SeedLen
d <- data %>% select(Couple, Species, Status, SeedLen) %>% na.omit()
p <- keep.tip(tree, d$Species)
phy.mat <- ape::vcv.phylo(p) # Phylo distance matrix
SeedLen <- brm(formula = SeedLen ~ 1 + (1|Status) + (1|gr(Species, cov = phy.mat)), 
               data = d, data2 = list(phy.mat = phy.mat), family = student(),
               prior = c(prior(normal(0, 1), class = "Intercept"),
                         prior(exponential(1), class = "sigma")),
               control = list(adapt_delta = delta, max_treedepth = treedepth),
               iter = iter, warmup = warmup, chains = chains, cores = chains, seed = 13,
               file = "output/SeedLen")

#######################################
###### Model Results ##################
setwd("G:/Mon Drive/Aegean endemics/Output_Hercules_Student")

PetLenMax <- readRDS("SVPetLenMax_mm.RDS")
PetLenMin <- readRDS("SVPetLenMin_mm.RDS")
SepLenMax <- readRDS("SepLenMax_mm.RDS")
SepLenMin <- readRDS("SepLenMin_mm.RDS")
SeedLenMax <- readRDS("SeedLenMax_mm.RDS")
SeedLenMin <- readRDS("SeedLenMin_mm.RDS")
FrLenMax <- readRDS("FrLenMax_mm.RDS")
FrLenMin <- readRDS("FrLenMin_mm.RDS")
LeafLenMax <- readRDS("LeafLenMax_mm.RDS")
LeafLenMin <- readRDS("LeafLenMin_mm.RDS")
LeafWidMax <- readRDS("LeafWidMax_mm.RDS")
LeafWidMin <- readRDS("LeafWidMin_mm.RDS")
PropLenMax <- readRDS("SVPropLenMax_mm.RDS")
PropLenMin <- readRDS("SVPropLenMin_mm.RDS")
StemLenMax <- readRDS("StemLenMax_cm.RDS")
StemLenMin <- readRDS("StemLenMin_cm.RDS")
StartFlo <- readRDS("StartFloPer.RDS")
EndFlo <- readRDS("EndFloPer.RDS")
FloPer <- readRDS("FloPer.RDS")
ChromNumb <- readRDS("ChromNo.RDS")
AltMin <- readRDS("MinAlt.RDS")
AltMax <- readRDS("MaxAlt.RDS")
AltRange <- readRDS("AltRange.RDS")
## Mean values
PetLen <- readRDS("PetLen.RDS")
SepLen <- readRDS("SepLen.RDS")
SeedLen <- readRDS("SeedLen.RDS")
FrLen <- readRDS("FrLen.RDS")
LeafLen <- readRDS("LeafLen.RDS")
LeafWid <- readRDS("LeafWid.RDS")
PropLen <- readRDS("PropLen.RDS")
StemLen <- readRDS("StemLen.RDS")

############## Extract from posterior ############
PetLenMaxPost <- PetLenMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("PetLenMax"))

PetLenMinPost <- PetLenMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("PetLenMin"))

SepLenMaxPost <- SepLenMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("SepLenMax"))

SepLenMinPost <- PetLenMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("SepLenMin"))

SeedLenMaxPost <- SeedLenMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("SeedLenMax"))

SeedLenMinPost <- SeedLenMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("SeedLenMin"))

FrLenMaxPost <- FrLenMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("FrLenMax"))

FrLenMinPost <- FrLenMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("FrLenMin"))

PropLenMaxPost <- PropLenMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("PropLenMax"))

PropLenMinPost <- PropLenMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("PropLenMin"))

StemLenMaxPost <- StemLenMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("StemLenMax"))

StemLenMinPost <- StemLenMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("StemLenMin"))

LeafLenMaxPost <- LeafLenMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("LeafLenMax"))

LeafLenMinPost <- LeafLenMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("LeafLenMin"))

LeafWidMaxPost <- LeafWidMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("LeafWidMax"))

LeafWidMinPost <- LeafWidMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("LeafWidMin"))

StartFloPost <- StartFlo %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Start Flowering"))

EndFloPost <- EndFlo %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("End Flowering"))

FlorPerPost <- FloPer %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("FloPer"))

ChromNumbPost <- ChromNumb %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("ChromNumb"))

AltMinPost <- AltMin %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Altitude Min"))

AltMaxPost <- AltMax %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Altitude Max"))

AltRangePost <- AltRange %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("AltRange"))

## Mean 
PetLenPost <- PetLen %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Petal Length")) 
SepLenPost <- SepLen %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Sepal Length")) 
SeedLenPost <- SeedLen %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Seed Length")) 
FrLenPost <- FrLen %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Fruit Length"))
LeafLenPost <- LeafLen %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Leaf Length"))
LeafWidPost <- LeafWid %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Leaf Width"))
PropLenPost <- PropLen %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Propagule Length"))
StemLenPost <- StemLen %>%
  spread_draws(r_Status[Endemic,Intercept]) %>%
  compare_levels(r_Status, by = Endemic) %>%
  mutate(phylo = rep("Stem Length"))

res <- rbind(PetLenMaxPost, PetLenMinPost,
             SepLenMaxPost,SepLenMinPost,
             SeedLenMaxPost, SeedLenMinPost,
             FrLenMaxPost, FrLenMinPost,
             PropLenMaxPost, PropLenMinPost,
             StemLenMaxPost, StemLenMinPost,
             LeafLenMaxPost, LeafLenMinPost,
             LeafWidMaxPost, LeafWidMinPost,
             StartFloPost, EndFloPost,
             FlorPerPost, ChromNumbPost,
             AltMinPost, AltMaxPost, AltRangePost,
             PetLenPost, SepLenPost, SeedLenPost, 
             FrLenPost, LeafLenPost, 
             LeafWidPost, PropLenPost, StemLenPost)
res %>%
  group_by(phylo) %>%
  median_hdi(condition_mean = r_Status, .width = 0.85) %>%
  View()

res %>%
  group_by(phylo) %>%
  median_hdi(condition_mean = r_Status, .width = 0.85) %>%
  ggplot(aes(y = phylo, x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = paste0("Trait Widespread - Traits Endemic"), 
       y = "", x = "") +
  theme_minimal()

p <- res %>%
  group_by(phylo)%>% 
  count(r_Status < 0) %>%
  mutate(p = n/16000) %>%
  select(-n) %>%
  pivot_wider(id_cols = phylo, names_from = 'r_Status < 0', values_from = p) %>%
  select(Trait = phylo, EW = 'FALSE', WE = 'TRUE')%>%
  ungroup() %>%
  mutate('E < W' = round(EW, 3), 'W < E' = round(WE, 3)) %>%
  select(-c(EW, WE))
View(p)
write.table(p, "Result_p.txt")