#' plot_deer_usa - plot relative number of samples by state.
#'
#' @param meta_data - a data.frame containing GISAID metadata.
#' @param title - string, an optional title
#'
#' from https://stackoverflow.com/questions/52911837/add-points-to-usmap-with-ggplot-in-r
#' 
plot_deer_usa <- function(meta_data,
                          title = NULL) {
  require(albersusa)
  require(tidyverse)
  require(sp)
  
  df <- meta_data %>% 
    filter(Host == 'Odocoileus virginianus') %>%
    select(Location) %>%
    separate(Location, into = c('Region', 'Country', 'State'), sep = ' / ') %>%
    group_by(State) %>%
    count() %>%
    rename(Count = n)
  
  us <- usa_composite(proj = "aeqd")
  
  states_centers <- as.data.frame(state.center)
  states_centers$name <- state.name
  
  states_centers <- states_centers[!(states_centers$name %in% c("Alaska", "Hawaii")),]
  states_centers <- states_centers %>% filter(name %in% df$State)
  
  coordinates(states_centers) <- ~x+y
  proj4string(states_centers) <- CRS(us_longlat_proj)
  states_centers <- spTransform(states_centers, CRSobj = CRS(us_aeqd_proj))
  states_centers <- as.data.frame(coordinates(states_centers))
  
  us_map <- fortify(us, region="name")
  
  states_centers <- states_centers %>% mutate(samples = df$Count)
  
  p <- ggplot() +
    geom_map(
      data = us_map, map = us_map,
      aes(x = long, y = lat, map_id = id),
      color = "#2b2b2b", size = 0.1, fill = NA
    ) +
    geom_point(
      data = states_centers, aes(x, y, size = samples), color = "steelblue"
    ) +
    coord_equal() + # the points are pre-projected
    ggthemes::theme_map()
  
  if(! is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  print(p)
  
  return(NULL)
}