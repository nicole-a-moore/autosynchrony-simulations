## exploring spongy moth expansion data 
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
theme_set(theme_bw())

## make file representing extent of quebec
## read in csd shapefile
csd <- read_sf("data-raw/csd-shapefile/lcsd000a24a_e/lcsd000a24a_e.shp") 

csd <- csd[which(csd$PRNAME == "Quebec / Québec"),]
unique(csd$PRNAME)

## plot
plot(csd)

## combine into 1 shape to define AOI
aoi = st_union(csd)
plot(aoi)
aoi = vect(aoi)

st_bbox(aoi)

## read in historic data on sbw distributuon
## note: levels of defoliation are represented different in this data than in later data 
tord_hist = vect("data-raw/TBE_Historique_1967-1991/Historique_TBE_1967_1991.shp") %>%
  makeValid(.)
names(tord_hist)

# Pour cette pÈriode, les classes de dommages causÈs par la tordeuse des bourgeons de líÈpinette sont les suivantes†:
#   
#   Classes de dÈfoliation 
# 
# 0. DÈfoliation nulle (0†%)
# 
# …valuation ‡ líextÈrieur des aires traitÈes (assessment outside of treated areas)
# 1. DÈfoliation lÈgËre†: de 1 ‡ 35 %
# 2. DÈfoliation modÈrÈe†: de 36 ‡ 70 %
# 3. DÈfoliation grave†: de 71 ‡ 100 %
# 
# …valuation dans les aires traitÈes (assessment inside treated areas)
# 4. DÈfoliation lÈgËre bloc†: de 1 ‡ 35 %
# 5. DÈfoliation grave tÍte verte†: de 36 ‡ 50 %
# 6. DÈfoliation modÈrÈe plus†: de 51 ‡ 70 %
# 7. DÈfoliation grave bloc†: de 71 ‡ 100 %
# 
# Classes de mortalitÈ
# 11 .MortalitÈ†: de 0 ‡ 25†%
# 12. MortalitÈ†: de 26 ‡ 50†%
# 13. MortalitÈ†: de 51 ‡ 75†%
# 14. MortalitÈ†: de 76 ‡ 100†%
# 15. MortalitÈ†: de 100†%


## seems like: 
# 1-7 indicate areas that show defoliation (1-3 have not been treated, 4-7 have)
# 11-15 indicate areas where trees have died because of defoliation 

## untreated areas 
tord_hist %>%
  filter(ANNEE == "1976") %>% filter(NIVEAU %in% 1:3) %>%
  ggplot(aes(colour = NIVEAU, fill = NIVEAU)) +
  geom_spatvector(size = 0.1)

## treated areas
## combine 5-6 to make same as 1:3
tord_hist %>%
  filter(ANNEE == "1976") %>% filter(NIVEAU %in% 4:7) %>%
  ggplot(aes(colour = NIVEAU, fill = NIVEAU)) +
  geom_spatvector(size = 0.1)

## treated areas
## combine 5-6 to make same as 1:3
tord_hist %>%
  filter(ANNEE == "1976") %>% filter(NIVEAU %in% 4:7) %>%
  mutate(NIVEAU = ifelse(NIVEAU %in% c(5,6), 2, ifelse(NIVEAU == 4, 1, ifelse(NIVEAU == 7, 3, NA)))) %>%
  ggplot(aes(colour = NIVEAU, fill = NIVEAU)) +
  geom_spatvector(size = 0.1)

## areas with mortality
tord_hist %>%
  filter(ANNEE == "1976") %>% filter(NIVEAU %in% 11:15) %>%
  ggplot(aes(colour = NIVEAU, fill = NIVEAU)) +
  geom_spatvector(size = 0.1)

tord_hist %>%
  mutate(col = ifelse(NIVEAU %in% 1:3, "1", ifelse(NIVEAU %in% 4:7, "1", "2"))) %>%
  filter(ANNEE == "1991") %>% 
  ggplot(aes(colour = col, fill = col, alpha = 0.05)) +
  geom_spatvector(size = 0.1)

## make column for treated vs. not treated areas 
## reclass categories so defoliation is from 1-3 (combine 5-6, change 4 to 1 and 7 to 3)
## add category for mortality (any amount)
hist = tord_hist %>%
  mutate(treated = ifelse(NIVEAU %in% 4:7, "treated", "untreated")) %>%
  mutate(NIVEAU = ifelse(NIVEAU %in% c(2, 5,6), 2, 
                         ifelse(NIVEAU %in% c(1,4), 1, 
                                ifelse(NIVEAU %in% c(3,7), 3, 
                                       ifelse(NIVEAU %in% c(11:15), 4, 
                                              ifelse(NIVEAU == 0, 0, NA)))))) 

aoi = project(aoi, crs(hist))

hist %>%
  # filter(ANNEE == "1970") %>% 
  filter(NIVEAU != 0) %>%
  ggplot(aes(colour = NIVEAU, fill = NIVEAU)) +
  geom_spatvector(size = 0.1, alpha = 0.1) +
  geom_spatvector(data = aoi, inherit.aes = FALSE, fill = "transparent") +
  scale_color_gradient(low = "black", high = "red", limits = c(1,4)) +
  scale_fill_gradient(low = "black", high = "red", limits = c(1,4)) +
  facet_wrap(~ANNEE)

## read in later data
tord_92_06 = vect("data-raw/TBE_Donnees_1992-2006/TBE_1992_2006.shp") %>% makeValid(.)
tord_07_13 = vect("data-raw/TBE_Donnees_2007-2013/TBE_2007_2013.shp") %>% makeValid(.)
tord_14_pres = vect("data-raw/TBE_Donnees_2014-aujourdhui/TBE_2014_2024.gdb") %>% makeValid(.)
tord_24 = vect("data-raw/TBE_Donnees_annee_en_cours_gdb/TBE_2024.gdb") %>% makeValid(.)

names(tord_92_06)
tord_92_06 %>%
  filter(ANNEE == "2006") %>% 
  ggplot(aes(colour = Ia, fill = Ia)) +
  geom_spatvector(size = 0.1)

names(tord_07_13)
tord_07_13 %>%
  filter(ANNEE == "2013") %>% 
  ggplot(aes(colour = Ia, fill = Ia)) +
  geom_spatvector(size = 0.1)

names(tord_14_pres)
tord_14_pres %>%
  filter(ANNEE == "2014") %>% 
  ggplot(aes(colour = Ia, fill = Ia)) +
  geom_spatvector(size = 0.1)

names(tord_24) 
tord_24 %>%
  ggplot(aes(colour = Ia24, fill = Ia24)) +
  geom_spatvector(size = 0.1)

## merge all together 
hist <- select(hist, ANNEE, NIVEAU)

tord_92_06 <- tord_92_06 %>%
  rename("NIVEAU" = Ia) %>%
  select(ANNEE, NIVEAU)

tord_07_13 <- tord_07_13 %>% 
  select(-NIVEAU) %>%
  rename("NIVEAU" = Ia) %>%
  select(ANNEE, NIVEAU)

tord_14_pres <- tord_14_pres %>%
  rename("NIVEAU" = Ia) %>%
  select(ANNEE, NIVEAU)

tord_24$ANNEE = "2024"
tord_24 <- tord_24 %>%
  select(-NIVEAU) %>%
  rename("NIVEAU" = Ia24) %>%
  select(ANNEE, NIVEAU)

## make sure all share same crs 
crs(hist) == crs(tord_92_06)
crs(hist) == crs(tord_07_13)
crs(hist) == crs(tord_14_pres)
crs(hist) == crs(tord_24)
tord_14_pres = project(tord_14_pres, crs(hist))
tord_24 = project(tord_24, crs(hist))

tord_24$ANNEE = as.integer(tord_24$ANNEE)

tord_all <- bind_spat_rows(hist, tord_92_06, tord_07_13, tord_14_pres, tord_24)

writeVector(tord_all, "outputs/data-processed/tordeuse_allyears.shp", overwrite = TRUE)

rm("hist", "tord_92_06", "tord_07_13", "tord_14_pres", "tord_24")

tord_all <- vect("outputs/data-processed/tordeuse_allyears.shp")
tord_all$ANNEE = as.numeric(as.character(tord_all$ANNEE))

## plot yearly distribution 
for(i in 1:length(unique(tord_all$ANNEE))) {
  
  annee = unique(tord_all$ANNEE)[order(unique(tord_all$ANNEE))][i]

  t = tord_all %>%
    filter(NIVEAU != 0) %>%
    filter(ANNEE == annee) 
  
  if(nrow(t) == 0) {
    plot = t %>%
      ggplot(aes()) +
      geom_spatvector(size = 0.1, alpha = 0.1) +
      geom_spatvector(data = aoi, inherit.aes = FALSE, fill = "transparent") +
      labs(title = annee)
  }
  else {
    plot = t %>%
      ggplot(aes(fill = NIVEAU)) +
      geom_spatvector(size = 0.1, colour = "transparent") +
      geom_spatvector(data = aoi, inherit.aes = FALSE, fill = "transparent") +
      #scale_color_gradient(low = "black", high = "red", limits = c(1,4)) +
      scale_fill_gradient(low = "black", high = "red", limits = c(1,4)) +
      labs(title = annee)
  }
  
  ggsave(plot, path = "outputs/figures/tordeuse", filename = paste0("tordeuse_", annee, ".png"), width = 12, height = 8)
  
  print(annee)
}
## note: data missing from 2019

## covert polygons to points representing centroids of each polygon 
cent = centroids(tord_all, inside = TRUE)

## plot centroids to check 
tord_all %>%
  filter(ANNEE == "1975") %>%
  ggplot(aes(fill = NIVEAU)) +
  geom_spatvector(size = 0.1, colour = "transparent") +
  geom_spatvector(data = aoi, inherit.aes = FALSE, fill = "transparent") +
  #scale_color_gradient(low = "black", high = "red", limits = c(1,4)) +
  scale_fill_gradient(low = "black", high = "red", limits = c(1,4)) +
  geom_spatvector(data = filter(cent, ANNEE == "1975"), size = 0.1)

## extract x-y values
coords = crds(cent)
coords_df <- vect(coords)

## extract attributes from underlying polygons
ext = extract(tord_all, coords_df)


i = 1
while(i <= length(unique(tord_all$ANNEE))) {
  ## split by year
  cur = filter(tord_all, ANNEE == unique(tord_all$ANNEE)[i])
  
  ## split by level
  lev = 1
  while(lev <= length(unique(cur$NIVEAU))) {
    
    cur2 = filter(cur, NIVEAU == unique(cur$NIVEAU)[lev])
    
    if(nrow(cur2) != 0) {
      ## extract x-y values
      coords = crds(centroids(cur2, inside = TRUE))
      coords_df <- vect(coords)
      
      ## extract attributes from underlying polygons
      ext = extract(cur2, coords_df) %>%
        unique()
      ext$x = coords[,1]
      ext$y = coords[,2]
      
      ## save 
      if(i == 1 & lev == 1) {
        ext_all = ext
      } 
      else {
        ext_all = rbind(ext, ext_all)
      }
    }
    lev = lev + 1
  }
  
  i =  i + 1
 print(unique(tord_all$ANNEE)[i])
}

ext_all = ext_all %>% 
  select(-id.y)

## save
write.csv(ext_all, "outputs/data-processed/tordeuse_centroids.csv", row.names = F)

ext_all <- read.csv("outputs/data-processed/tordeuse_centroids.csv")

## make quantile regression plot 
ext_all %>%
  filter(NIVEAU != 0) %>%
  ggplot(aes(x = ANNEE, y = y, colour = NIVEAU)) + 
  geom_point(size = 0.1, position = position_jitter()) +
  scale_color_gradient(low = "black", high = "red", limits = c(1,4)) 

## quantile regression
ext_all %>%
  filter(NIVEAU != 0) %>%
  ggplot(aes(x = ANNEE, y = y, colour = NIVEAU)) + 
  geom_point(size = 0.1, position = position_jitter()) +
  scale_color_gradient(low = "black", high = "red", limits = c(1,4)) +
  geom_quantile(quantiles = c(0.05, 0.5, 0.95), colour = "black")

ggsave(path = "outputs/figures", filename = "tordeuse_outbreak-shift.png", width = 8, height = 4)




