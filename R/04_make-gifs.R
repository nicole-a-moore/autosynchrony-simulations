## make range shift gifs
library(magick)
library(magrittr)
library(tidyverse)

dir = "outputs/figures/range-shift-snapshots"

folders = list.files(dir, full.names = T)[which(str_detect(list.files(dir), "rep"))]
folders = folders[which(!str_detect(folders, "gif"))]

## for each folder of range shift png images
f = 1
while(f <= length(folders)) {
  ## order files by time
  files = list.files(folders[f], full.names = TRUE)
  
  t = str_split_fixed(files, "/t", n = 2)[,2]
  t = str_split_fixed(t, "_", n = 2)[,1]
  files = files[order(as.numeric(t))]
  
  ## make a gif
  frames = files %>% 
    image_read() %>% # reads each path file
    image_join() 
  
  ## make sure none failed
  frames %>% # joins images
    image_apply(function(img) image_background(img, "white")) %>%
    image_animate(fps = 25) %>% # animates, can opt for number of loops
    image_write(paste0(folders[f], ".gif")) # write to current dir
  
  f = f + 1
}

## now make gifs of some range edge shifts 
dir = "outputs/figures/range-edge-shifts/Cardinalis_cardinalis"

files = list.files(dir, full.names = T)[which(str_detect(list.files(dir, full.names = T), "cardinalis\\/BBS\\_"))]

## get rid of 1966-1980
files = files[-1]

## order by year
t = str_split_fixed(files, "BBS_Cardinalis_cardinalis_", n = 2)[,2]
t = str_split_fixed(t, "\\.png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 5) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/range-edge-shifts/Cardinalis_cardinalis/full-edge-shift", ".gif")) # write to current dir


files = list.files(dir, full.names = T)[which(str_detect(list.files(dir, full.names = T), "cardinalis\\/intraspecific_BBS\\_"))]

## get rid of 1966-1980
files = files[-1]

## order by year
t = str_split_fixed(files, "BBS_Cardinalis_cardinalis_", n = 2)[,2]
t = str_split_fixed(t, "\\.png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 5) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/range-edge-shifts/Cardinalis_cardinalis/intraspecific-range-edge-shift", ".gif")) # write to current dir





## common raven
dir = "outputs/figures/range-edge-shifts/Corvus_corax"

files = list.files(dir, full.names = T)[which(str_detect(list.files(dir, full.names = T), "corax\\/BBS\\_"))]

## get rid of 1966-1980
files = files[-1]

## order by year
t = str_split_fixed(files, "BBS_Corvus_corax_", n = 2)[,2]
t = str_split_fixed(t, "\\.png", n = 2)[,1]
files = files[order(as.numeric(t))]

## make a gif
frames = files %>% 
  image_read() %>% # reads each path file
  image_join() 

## make sure none failed
frames %>% # joins images
  image_apply(function(img) image_background(img, "white")) %>%
  image_animate(fps = 5) %>% # animates, can opt for number of loops
  image_write(paste0("outputs/figures/range-edge-shifts/Corvus_corax/full-edge-shift", ".gif")) # write to current dir

