###########################
###########################
#### patter-intro.R

#### Aims
# (1) Introduce the `patter` R package

#### Prerequisites
# (1) This code runs on Windows/MacOS
# (patter@dev can run on linux, with modifications to this workflow)

#### Issues
# Please report questions or issues on GitHub:
# https://github.com/edwardlavender/patter/issues
# https://github.com/edwardlavender/patter-workshops/issues


###########################
###########################
#### Set up

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(glue)
library(lubridate)
library(patter)
library(spatstat.explore)

#### Load data
# We will load example datasets from `patter`
# * These data are from flapper skate (_Dipturus intermedius_)
map        <- dat_gebco()
moorings   <- copy(dat_moorings)
detections <- copy(dat_acoustics)
archival   <- copy(dat_archival)

#### Julia set up (~5 s after first run!)
# Connect to Julia
# * Set JULIA options in .Rprofile or .Renviron
# * If you experience Julia errors/warning, try re-running `julia_connect()`
julia_connect()
# Set seed in R and Julia
set_seed()


###########################
###########################
#### Workflow with real data


###########################
#### Define study system

#### Define study area
# `map` is a _planar_ SpatRaster that defines the region within which movements are permitted
# * See `terra::project()` to convert lon/lat to UTM
# * See `terra::rasterize()` to convert shapefile to SpatRaster
terra::plot(map)
# Use `set_map()` to export the map to `Julia`
set_map(map)

#### Define study timeline
# We focus on a short period for this example
# The resolution of the timeline defines the resolution of the simulation
timeline <- seq(as.POSIXct("2016-03-26", tz = "UTC"),
                length.out = 1000, by = "2 mins")
# Reduce the timeline if your computer struggles!
# timeline <- timeline[1:500]
timespan <- lubridate::interval(min(timeline), max(timeline))


###########################
#### Movement model

#### Define state
state      <- "StateXY"

#### Define movement model
# This is a simple movement model for illustrative purposes
# `mobility` is the _maximum_ distance travelled between between two time steps
# In our case, this needs to account for movement capacity, current speed etc.
mobility   <- 1100
model_move <-
  move_xy(dbn_length = glue("truncated(Gamma(1, 250.0), upper = {mobility})"),
          dbn_angle = "Uniform(-pi, pi)")

#### Visualise movement model
sim_path_walk(.map = map,
              .timeline = timeline,
              .state = state,
              .model_move = model_move,
              .n_path = 4L, .one_page = TRUE)

# (Optional) tasks
# * How would you modify this movement model to use different parameters?
# * How would you change the distributions in the model?
# * Why doesn't the movement model ever jump on land?
# * Can the simulated individual leave the study area?
# * What happens if you set:
# *     dbn_angle = VonMises(0.0, 0.25)?
# *     dbn_angle = VonMises(0.0, 1.0)? (Be prepared to kill R!) Why?
# * Can you formulate a suitable model for one of your study species?
# * Can you visualise how distinct movement models differ?


###########################
#### Observation model(s)

#### Process datasets
# `patter` does not process datasets for you!

#### Define datasets _for a selected animal and time period_
# Acoustic receivers & detections
moorings <-
  moorings |>
  as.data.frame() |>
  mutate(receiver_timespan = interval(receiver_start, receiver_end)) |>
  filter(lubridate::int_overlaps(receiver_timespan, timespan)) |>
  mutate(receiver_timespan = NULL) |>
  as.data.table()
det <-
  detections |>
  filter(individual_id == 28L) |>
  filter(timestamp %within% timespan) |>
  as.data.table()
# Depth observations
arc <-
  archival |>
  filter(individual_id == 28L) |>
  filter(timestamp %within% timespan) |>
  as.data.table()
# (optional) Validate acoustic/archival datasets
datalist <- pat_setup_data(.map = map,
                           .acoustics = det,
                           .moorings = moorings,
                           .services = NULL,
                           .archival = arc)

#### Assemble observations (acoustics)
# Assemble a timeline of acoustic observations (0, 1) and model parameters
# * We include ModelObsAcousticLogisTrunc model parameters via `.moorings`
# * We revise the detection range
moorings[, receiver_gamma := 1750]
acc <- assemble_acoustics(.timeline = timeline,
                          .acoustics = det,
                          .moorings = moorings)

#### Assemble archival observations
# Assemble a timeline of archival observations and model parameters
# * Here, we include model parameters for `ModelObsDepthNormalTrunc`
# * We need to include a large uncertainty b/c the bathymetry in this region is complex
# * And the example bathymetry dataset is extremely coarse
arc <- assemble_archival(.timeline = timeline,
                         .archival =
                           arc |>
                           rename(obs = depth) |>
                           mutate(depth_sigma = 100,
                                  depth_deep_eps = 100))

#### Collect observations
yobs <- list(ModelObsAcousticLogisTrunc = acc,
             ModelObsDepthNormalTrunc = arc)

# (Optional) tasks
# * Plot the observation models
# * Plot the observed time series
# * How would you assemble another kind of dataset? (Hint, see `?assemble`)
# * Can you assemble your own datasets for the particle filter?
# * Where can you find out more about the built-in `ModelObs` structures?
# * How would you implement another `ModelObs` structure?
# * Can you implement the `ModelObsDepthUniform` structure for these data?


###########################
#### Model inference

#### Run forward filter (~12 s)
fwd <- pf_filter(.map = map,
                 .timeline = timeline,
                 .state = state,
                 .xinit = NULL, .xinit_pars = list(mobility = mobility),
                 .yobs = yobs,
                 .model_obs = names(yobs),
                 .model_move = model_move,
                 .n_particle = 1e5L)

# (optional) Tasks
# * How would you set the initial location in the filter?
# * What happens if you set an initial location?
# * How does computation time change in the filter if you change the timeline or .n_particle?
# * What happens if you change other arguments, such as .n_record?
# * How can we plot the probability distribution for an individual's location at a selected timestep?
# * What about creating an animation of particle movement (hint: see ?`pf_plot_xy()`)?
# * What checks would you do to check `pf_filter()` has worked as expected?
# * At the moment of detection, how far are particles from receivers? Why?
# * How would you analyse filter diagnostics?
# * Are filter diagnostics acceptable in this example?
# * How could we improve filter diagnostics?
# * If we write the outputs from fwd to file, how much storage space do we need?
# * What does that tell us about the storage requirements of a full-scale analysis?

#### Run backward information filter (~12 s)
bwd <- pf_filter(.map = map,
                 .timeline = timeline,
                 .state = state,
                 .xinit = NULL, .xinit_pars = list(mobility = mobility),
                 .yobs = yobs,
                 .model_obs = names(yobs),
                 .model_move = model_move,
                 .n_particle = 1e5L,
                 .direction = "backward")

# (optional) Tasks
# * How do the convergence properties of the forward filter differ from the backward filter?
# * Why do you think there are differences in this case?
# * How might you solve them?

#### Run two-filter smoother (~1 s)
smo <- pf_smoother_two_filter(.n_particle = 100L)

# (optional) tasks
# * Repeat the tasks from the filter to examine the behaviour of the smoother
# * How does computation time change if you change .n_particle?
# * How does this differ from the smoother?
# * How repeatable are multiple runs of the particle filter/smoother?

#### Map utilisation distribution with/without kernel smoothing (~1 s)
pp      <- par(mfrow = c(1, 2))
ud_pou  <- map_pou(.map = map, .coord = smo$states)
ud_dens <- map_dens(.map = map,
                   .coord = smo$states,
                   .discretise = TRUE,
                   .fterra = TRUE,
                   sigma = bw.diggle)
par(pp)

# (optional) tasks
# * How do we interpret these maps?
# * Are you satisfied with this result? If not, why not? How would you improve it?
# * For `map_pou()`, how do the results change with pixel resolution?
# * For `map_dens()`, what is the pixel resolution used by the estimation routine?
# * Does this differ from the resolution of the map?
# * How does `map_dens()` depend on the `sigma` argument?
# * (Advanced) How would you write a loop to trial different sigma parameters?
#   * Hint: see `?patter::cl_lapply()` and `?terra::wrap()`
#   * Can you boost efficiency by parallelising over chunks?
# * How can you boost the speed of utilisation-distribution estimation?
# * How sensitive are these results to our movement and observation models?
# * Can you create the same maps from the filters?
# * How do they differ from maps from the smoother and why?
# * How would you add core/home ranges to these maps (hint: see `?map_hr`)?
# * How would you quantify residency in a 1 km radius around (709619.4, 6255728)?


###########################
###########################
#### Workflow with simulated data

###########################
#### Study system

# We will use the timeline and map defined above.

# (optional) tasks
# * Experiment with this choice


###########################
#### Movement model

# We use the movement model defined above.

# (optional) tasks
# * Experiment with this choice

#### Simulate a 'true' path
path_coord <- sim_path_walk(.map = map,
                            .timeline = timeline,
                            .state = state,
                            .model_move = model_move)

#### Estimate a 'true' pattern of space use
path_ud <- map_dens(.map = map,
                    .coord = path_coord,
                    .discretise = TRUE,
                    .fterra = TRUE,
                    sigma = bw.diggle)


###########################
#### Simulate observations

#### Define sensor parameters for the ModelObsAcousticLogisTrunc struct
# In this example, we use the real receiver positions
pars_ModelObsAcousticLogisTrunc <-
  moorings |>
  select(sensor_id = "receiver_id",
         "receiver_x", "receiver_y",
         "receiver_alpha", "receiver_beta", "receiver_gamma") |>
  as.data.table()

#### Define sensor parameters for the ModelObsDepthNormalTrunc struct
pars_ModelObsDepthNormalTrunc <-
  data.table(sensor_id = 1L,
             depth_sigma = 100,
             depth_deep_eps = 100)

#### Collect sensor parameters for each data type
pars_ModelObs <- list(ModelObsAcousticLogisTrunc = pars_ModelObsAcousticLogisTrunc,
                      ModelObsDepthNormalTrunc = pars_ModelObsDepthNormalTrunc)

#### Simulate observations
obs <- sim_observations(.timeline = timeline,
                        .model_obs = names(pars_ModelObs),
                        .model_obs_pars = pars_ModelObs)

# (optional) tasks
# * Why is the `sensor_id` field included in all `ModelObs` structures?
# * What does `sensor_id = 1L` mean in `pars_ModelObsDepthNormalTrunc`?
# * For a given `ModelObs` structure, how do we find out the required parameters?
# * How would you simulate observations for a hypothetical receiver array?
# * Can you simulate observations using another `ModelObs` structure?
# * What is the structure of the object returned by `sim_observations()`?
# * Can you plot the simulated observations?
# * How would you check the simulated observations align with the observation models?


###########################
#### Model inference

# (optional) tasks
# * Can you implement the filter and smoother using simulated data?
# * How do reconstructed maps match the simulated data?
# * Why are there differences?
# * Can you create a comparable map using the `coa()` function?
# * How does a map based on COAs differ?
# * How does the map change if you only model one data type?
# * What does this tell us about the relative contribution of the two datasets in this system?
# * If you were planning a study in this system, would you invest in more receivers or more archival tags based on these results? Why?


#### End of code.
###########################
###########################
