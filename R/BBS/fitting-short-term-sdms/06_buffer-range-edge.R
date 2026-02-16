## make polygon representing buffer around species range edges 
library(sf)

## read in list of species with study divided into longitudinal blocks
bbs <- read.csv("outputs/data-processed/BBS_spatial-filter.csv") 

## get canada and us map for plotting
countries <- necountries::countries(c("Canada", "United States of America"), part = TRUE)
countries <- st_union(countries)
# countries <- vect(countries)

## for each species
sp = 1 
while(sp <= length(unique(bbs$genus_sp))) {
  
  ## get species name 
  cur_sp = unique(bbs$genus_sp)[sp]
  
  ## get occurrences across all years from BBS
  occ = filter(bbs, genus_sp == cur_sp) %>%
    filter(!TotalAbd == 0 & !is.na(TotalAbd)) 
    
  ## make data frame an sf data frame
  occ <- vect(occ, geom = c("Longitude", "Latitude")) %>%
    st_as_sf(.) %>%
    select(geometry) %>%
    unique() 
  
  st_crs(occ) = st_crs(countries)
  
  ## draw buffer around points
  buffer = buffer(vect(occ), width = 100000)
  #plot(buffer)
  
  ## dissolve individual polygons into one
  buffer <- st_as_sf(buffer)
  buffer = st_union(buffer)
  #plot(buffer)
  
  ## smooth over
  buffer = buffer(vect(buffer), width = 100000)
  #plot(buffer)
  #plot(occ, add = T)
  
  ## crop by North America's edge
  buffer_na = crop(buffer, vect(countries))
  
  ## cut out middle
  buffer_inv = buffer(buffer_na, width = -500000)
  #plot(buffer_na)
  #plot(buffer_inv, add = T)
  
  buffer_final = erase(buffer_na, buffer_inv)
  #plot(buffer_final)
  
  ## make polygons representing longitudinal bins
  ## get longitudinal bins
  lon_bins <- bbs %>%
    filter(genus_sp == cur_sp) %>%
    select(Longitude_bin) %>%
    distinct() %>% 
    mutate(min = str_split_fixed(Longitude_bin, ",", 2)[,1],
           max = str_split_fixed(Longitude_bin, ",", 2)[,2]) %>%
    mutate(min = as.numeric(as.character(str_remove(min, "\\("))),
           max = as.numeric(as.character(str_remove(max, "\\]"))), 
           bin = 1:nrow(.)) %>%
    filter(!is.na(Longitude_bin))
  
  polys = c()
  for(i in 1:nrow(lon_bins)) {
    polys = st_as_sf(vect(ext(lon_bins$min[i], lon_bins$max[i], 18, 84))) %>%
      rbind(polys, .)
  }
  #plot(polys)
  st_crs(polys) = st_crs(countries)
  
  ## crop buffer by polygons and split into leading and trailing edge
  buffer_polys = c()
  for(i in 1:nrow(polys)) {
    ## perform crop
    cropped = crop(buffer_final, vect(polys[i,]))
    
    ## if there is a range edge within that polygon  
    if(nrow(cropped) != 0) {
      ## split into leading and trailing edge at mid latitude using bbox
      mid_lat = as.vector(ext(cropped)[3] + ext(cropped)[4]) / 2
      bbox_leading = ext(cropped)
      bbox_leading[3] = mid_lat ## ymin
      bbox_trailing = ext(cropped)
      bbox_trailing[4] = mid_lat ## ymax
      
      leading = crop(cropped, bbox_leading)
      trailing = crop(cropped, bbox_trailing)
      # plot(leading)
      # plot(trailing)
      
      ## combine and turn into sf
      leading = st_as_sf(leading) %>%
        mutate(range_edge = "leading")
      trailing = st_as_sf(trailing) %>%
        mutate(range_edge = "trailing")
      
      cropped = rbind(leading, trailing)
      
      ## save leading and trailing halves, cropped by polygon 
      buffer_polys = cropped %>%
        mutate(polyID = i) %>%
        rbind(buffer_polys, .) 
    }
    
  }
  
  ## plot result 
  plot(buffer_polys)
  
  ## change crs
  st_crs(buffer_polys) = st_crs(countries)
  
  ## save for species
  write_sf(buffer_polys, paste0("outputs/data-processed/range-edge-buffers/", str_replace_all(cur_sp, " ", "_"), "_buffer.shp"))
  
}






