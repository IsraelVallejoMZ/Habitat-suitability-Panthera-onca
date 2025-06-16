
# Final script for modeling the potential habitat suitability of Panthera onca in Mesoamerica

# Load required libraries
options(timeout = 300)
install.packages("ENMeval")
library(rgbif)
library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(maxnet)
library(ENMeval)
library(geodata)
library(viridis)

# 1. Download occurrence data for Panthera onca
occs_raw <- occ_search(
  scientificName = "Panthera onca",
  hasCoordinate = TRUE,
  limit = 2000
)

# 2. Clean and convert to spatial object
occs <- occs_raw$data %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

occs_sf <- st_as_sf(occs, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
plot(st_geometry(occs_sf), main = "Panthera onca Occurrences")

# 3. Download bioclimatic variables
options(timeout = 600)
bio_current <- worldclim_global(var = "bio", res = 10, path = "data")
bio_future_files <- list.files("C:/Users/israe/OneDrive/Desktop/R/Habitat suitability Panthera onca/data_bio2050", pattern = ".tif$", full.names = TRUE)
bio_future <- rast(bio_future_files)

# 4. Define region of interest and crop layers
ext <- ext(-120, -75, 0, 35)
bio_current_crop <- crop(bio_current, ext)
bio_future_crop <- crop(bio_future, ext)
plot(bio_current_crop[[1]], main = "Bio1 - Current Climate")
plot(bio_future_crop[[1]], main = "Bio1 - Future Climate")

# 5. Apply elevation mask
elev <- geodata::elevation_global(res = 10, path = "data")
elev_crop <- crop(elev, ext)
elev_resample <- resample(elev_crop, bio_current_crop)
elev_mask <- elev_resample < 3500
bio_current_masked <- mask(bio_current_crop, elev_mask)
bio_future_masked  <- mask(bio_future_crop, elev_mask)

# 6. Prepare modeling dataset
occs_vect <- vect(occs_sf)
occs_vals <- terra::extract(bio_current_masked, occs_vect)
occs_vals <- occs_vals[complete.cases(occs_vals), ]
occs_vals$presence <- 1

set.seed(42)
bg_points <- spatSample(bio_current_masked[[1]], size = 10000, method = "random", na.rm = TRUE, as.points = TRUE)
bg_vals <- terra::extract(bio_current_masked, bg_points)
bg_vals <- bg_vals[complete.cases(bg_vals), ]
bg_vals$presence <- 0

all_data <- rbind(occs_vals, bg_vals)
all_data <- all_data[ , !names(all_data) %in% "ID"]
predictors <- all_data[, !(names(all_data) %in% "presence")]
response <- all_data$presence

# 7. Train the model
model <- maxnet(p = response, data = predictors, f = maxnet.formula(response, predictors, classes = "default"))
predicted <- predict(model, predictors, type = "cloglog")

# 8. Evaluate the model
eval <- dismo::evaluate(p = predicted[response == 1], a = predicted[response == 0])
print(eval@auc)

# 9. Project suitability
suit_current <- predict(bio_current_masked, model, type = "cloglog", na.rm = TRUE)
plot(suit_current, main = "Habitat suitability - Present")
names(bio_future_masked) <- names(predictors)
suit_future <- predict(bio_future_masked, model, type = "cloglog", na.rm = TRUE)
plot(suit_future, main = "Habitat suitability - 2041–2060 (SSP5-8.5)")

# 10. Binarize maps and analyze changes
threshold <- eval@t[which.max(eval@TPR + eval@TNR)]
print(threshold)
suit_present_bin <- suit_current >= threshold
suit_future_bin  <- suit_future  >= threshold
change_map <- suit_present_bin + suit_future_bin * 2

change_map <- classify(change_map, rcl = matrix(c(
  0, 0, NA,
  1, 1, 1,
  2, 2, 2,
  3, 3, 3
), ncol = 3, byrow = TRUE))

levels(change_map) <- data.frame(value = c(1, 2, 3), category = c("Loss", "Gain", "Stable"))
plot(change_map, col = c("red", "green", "yellow"), main = "Habitat change 2041–2060")

# 11. Calculate proportions
freq <- freq(change_map)
freq <- freq[!is.na(freq$value), ]
total_pixels <- sum(freq$count)
freq$percent <- round(100 * freq$count / total_pixels, 2)
print(freq)

# 12. Add base map and visualize with ggplot2
world <- geodata::gadm(country = c("Mexico", "Guatemala", "Belize", "Honduras", "El Salvador", "Nicaragua", "Costa Rica", "Panama"), level = 0, path = "data")
world_proj <- project(world, crs(change_map))
world_sf <- sf::st_as_sf(world)

change_df <- as.data.frame(change_map, xy = TRUE, na.rm = TRUE)
colnames(change_df) <- c("x", "y", "class")
change_df$class <- factor(as.numeric(change_df$class), levels = c(1, 2, 3), labels = c("Loss", "Gain", "Stable"))

ggplot() +
  geom_tile(data = change_df, aes(x = x, y = y, fill = class)) +
  geom_sf(data = world_sf, fill = NA, color = "gray70", size = 0.3) +
  scale_fill_manual(name = "Habitat Change",
                    values = c("Loss" = "#D7191C", "Gain" = "#1A9641", "Stable" = "#FFD700")) +
  coord_sf(xlim = c(-120, -75), ylim = c(0, 35), expand = FALSE) +
  theme_void() +
  labs(title = "Habitat change - 2041–2060 (SSP5-8.5)") +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), legend.position = "right")
