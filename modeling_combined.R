#' ---
#' title: "Modeling of whale calling behavior"
#' author: ""
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE,
                      fig.width = 10, fig.height = 6)


#' This code combines univariate and multivariate analyses of whale calls.
#'
#' In RStudio, use File -> Compile Report to generate PDF from this R file.

# Clean the workspace
rm(list = ls())

# Packages ----

# data manipulation
library(chron)
library(dplyr)
library(tidyr)
library(readxl)
library(psych)
library(lubridate)
library(suncalc)

# plots
library(ggplot2)
theme_set(theme_light())
library(patchwork)
library(GGally)
library(gridExtra)

# models
library(flexmix)
library(car)
library(nlme)
library(randomForestSRC)
library(olsrr)
# library(rpart)
# library(rpart.plot)
# library(mgcv)
library(randomForest)
# library(ranger)

NTREE <- 500
MTRY <- 2
NSIZE <- 5

get_decimal_hour <- function(dt) {
    lubridate::hour(dt) + lubridate::minute(dt) / 60 + lubridate::second(dt) / 3600
}

# Get partial dependence data for all predictors
get_pdp_data <- function(model, data, predictors, ...) {
    require(pdp)
    pdp_list <- lapply(predictors, function(pred) {
        pd <- pdp::partial(model,
                           pred.var = pred,
                           train = data,
                           plot = FALSE,
                           ...)
        pd$predictor <- pred
        names(pd)[1] <- "x"
        pd
    })
    do.call(rbind, pdp_list)
}

# Inverse of log10, such as for inv.link = invlog10 in get_pdp_data()
invlog10 <- function(x) {
    10^x
}

# Opposite of in
`%!in%` <- Negate(`%in%`)

# Data ----

BuoyReplacement <- as.Date(c("2022-07-20"))

# Load and format the data
D <- readxl::read_xlsx("./data/NARWData_2024-07-29.xlsx",
                       sheet = 1,
                       col_types = c(rep("guess", 3),
                                     rep("numeric", 3),
                                     rep("guess", 9))) %>%
    # format dates and times
    mutate(Date = as.Date(Date),
           Time_NYC = chron::times(Time_NYC),
           TimeLastVP_NYC = chron::times(TimeLastVP_NYC)) %>%
    mutate(IsVPYesterday = TimeLastVP_NYC > Time_NYC) %>%
    mutate(Time_NYC = as.POSIXct(paste(Date, Time_NYC), tz = "America/New_York"),
           TimeLastVP_NYC = as.POSIXct(paste(Date - IsVPYesterday, TimeLastVP_NYC),
                                       tz = "America/New_York")) %>%
    # calculate time variables
    mutate(TimeSinceLastVP = as.numeric(difftime(Time_NYC, TimeLastVP_NYC, units = "hours")),
           Year = as.numeric(format(Time_NYC, "%Y")),
           Month = as.numeric(format(Time_NYC, "%m")),
           Hour = as.numeric(format(Time_NYC, "%H")),
           HourDecimal = get_decimal_hour(Time_NYC),
           DoY = as.numeric(format(Time_NYC, "%j"))) %>%
    mutate(DoYsin = sin(2 * pi * DoY / 365.25),
           DoYcos = cos(2 * pi * DoY / 365.25),
           Hoursin = sin(2 * pi * HourDecimal / 24),
           Hourcos = cos(2 * pi * HourDecimal / 24)) %>%
    # define seasons based on months
    mutate(Season = case_when(
        Month %in% 9:11 ~ "Fall",
        Month %in% c(12, 1, 2) ~ "Winter",
        Month %in% 3:5 ~ "Spring",
        TRUE ~ "Summer")) %>%
    # order seasons chronologically
    mutate(Season = factor(Season, levels = c("Fall", "Winter", "Spring", "Summer"))) %>%
    # drop unused factor levels
    mutate(Season = droplevels(Season)) %>%
    # rename some variables
    rename(Duration = DeltaTime,
           Duration90 = Dur90,
           SPLvessel = SPL_VesselPassage,
           SPLcalling = SPL_CallingPeriod,
           TimeSinceLVP = TimeSinceLastVP) %>%
    # sort chronologically
    arrange(Time_NYC) %>%
    # calculate inter-call intervals
    mutate(ICI = c(NA, as.numeric(difftime(Time_NYC[-1], Time_NYC[-length(Time_NYC)], units = "secs")))) %>%
    # year of deployment
    mutate(YearDeployment = case_when(
        Date < BuoyReplacement ~ "Year 1",
        Date > BuoyReplacement ~ "Year 2"))

# Set the ICI between years of observations as missing (NA).
# Note that the end of the previous season is March and start of the season is November.
tmp_PreviousMonth <- c(NA, D$Month[-nrow(D)])
YearSwitch <- which(!is.na(tmp_PreviousMonth) & (tmp_PreviousMonth == 3 & D$Month == 11))
D$ICI[YearSwitch] <- NA

# Show summary
summary(D)



# Bouts ----

## Define bouts ----

# Plot ICI distribution
D %>%
    ggplot(aes(x = ICI / 60)) +
    geom_histogram(aes(y = ..count..), fill = "blue", alpha = 0.5, color = NA) +
    scale_x_log10() +
    labs(x = "Inter-call interval (minutes)",
         y = "Frequency")
ggsave("images/ICI_hist.png", width = 6, height = 4, dpi = 600)

# Fit a mixture distribution (2 gamma distributions)
set.seed(123)
ICIs <- na.omit(D$ICI)
fit <- flexmix(ICIs ~ 1,
               model = FLXMRglm(ICIs ~ ., family = "Gamma"),
               k = 2)
summary(fit)
parameters(fit)

# Plot histogram of posterior probabilities
plot(fit, root = FALSE, main = "", eps = 0.1,
     xlab = "Probability", ylab = "Frequency")

png("images/ICI_flexmix.png", width = 6, height = 4, units = "in", res = 600)
plot(fit, root = FALSE, main = "", eps = 0.1,
     xlab = "Probability", ylab = "Frequency")
dev.off()

# In fit@cluster, 1s indicate that the current ICI belongs to the bout,
# and 2s indicate that there's a break between bouts (start of a new bout).
# The ICIs miss the 1st value in each year of deployment (which are NA), so it is n + 2.
# Sometimes the definition of 1 and 2 switches.

# Summarize ICIs (minutes) by cluster
psych::describeBy(ICIs/60, fit@cluster)

# Define bouts
D$Bout <- NA
# Make sure that calls at the YearSwitch are also the breaks
ICIclusters <- fit@cluster
for (i in 1:length(YearSwitch)) {
    ii <- YearSwitch[i]
    ICIclusters <- c(ICIclusters[1:(ii - 2)], 2, ICIclusters[-c(1:(ii - 2))])
}

bi <- 1
D$Bout[1] <- bi
for (i in 1:length(ICIclusters)) {
    # If a break, start a new bout (increase bout index bi)
    if (ICIclusters[i] == 2) {
        bi <- bi + 1
    }
    D$Bout[i + 1] <- bi
}

# Summarize data for bout analysis
DB0 <- D %>%
    mutate(Bout = as.factor(Bout)) %>%
    group_by(Bout) %>%
    summarise(Date = Date[1],
              Year = Year[1],
              Month = Month[1],
              DoY = DoY[1],
              Hour = Hour[1],
              Season = Season[1],
              DoYcos = DoYcos[1],
              DoYsin = DoYsin[1],
              Hoursin = Hoursin[1],
              Hourcos = Hourcos[1],
              HourDecimal = HourDecimal[1],

              ICI = mean(ICI, na.rm = TRUE),
              n_calls = n(),
              BoutDur_s = as.numeric(difftime(tail(Time_NYC, 1), Time_NYC[1], units = "secs")),

              Time_NYC = Time_NYC[1],
              TimeSinceLVP = TimeSinceLVP[1],
              SPLcalling = mean(SPLcalling),
              SPLvessel = mean(SPLvessel),

              FreqMin_first = FreqMin[1],
              FreqMin_last = tail(FreqMin, 1),
              FreqMin = mean(FreqMin),

              FreqMax_first = FreqMax[1],
              FreqMax_last = tail(FreqMax, 1),
              FreqMax = mean(FreqMax),

              FreqDelta_first = FreqDelta[1],
              FreqDelta_last = tail(FreqDelta, 1),
              FreqDelta = mean(FreqDelta),

              Duration_first = Duration[1],
              Duration_last = tail(Duration, 1),
              Duration = mean(Duration),

              Duration90_first = Duration90[1],
              Duration90_last = tail(Duration90, 1),
              Duration90 = mean(Duration90)
    ) %>%
    # convert duration to minutes
    mutate(BoutDur_mins = BoutDur_s / 60) %>%
    mutate(BoutDur_mins_log10 = log10(BoutDur_mins))

# Total number of bouts
max(D$Bout)

# Summary of bout stats
DB0 %>%
    dplyr::select(BoutDur_mins, n_calls) %>%
    psych::describe()

# Count bouts with one call only
db1 <- DB0 %>%
    filter(n_calls == 1)
nrow(db1)

# Remove the short bouts from analysis
DB <- DB0 %>%
    filter(n_calls > 1)

write.csv(DB, file = "dataderived/data_bouts.csv", row.names = FALSE)


## Bout summary ----

# Save numeric summary of the data as a table
summary_DB <- psych::describe(DB) %>%
    as.data.frame() %>%
    dplyr::select(n, min, mean, median, max, sd, range)

# Show summary to match with the table in the paper
summary_DB
write.csv(summary_DB, 'dataderived/summary_DB.csv')


# Summarize like Table 3 in Davis et al. (2023)
# https://doi.org/10.1093/icesjms/fsad174
# but by Year and Month and add columns (Bout_dur_mins_sd)

dtmp <- DB %>%
    group_by(Year, Month) %>%
    summarize(n_bouts = n(),
              BoutDur_mins_mean = mean(BoutDur_mins),
              BoutDur_mins_sd = sd(BoutDur_mins),
              BoutDur_mins_95 = quantile(BoutDur_mins, probs = 0.95)
    )

SumTable <- D %>%
    # remove bouts formed by a single upcall
    filter(!is.element(Bout, db1$Bout)) %>%
    group_by(Year, Month) %>%
    summarize(n_calls = n(),
              ICI_95_min = quantile(ICI[-1], probs = 0.95, na.rm = TRUE) / 60
    ) %>%
    left_join(dtmp)

knitr::kable(SumTable, digits = 2)


# Alternatively, summary by year

dtmp <- DB %>%
    group_by(Year) %>%
    summarize(n_bouts = n(),
              BoutDur_mins_mean = mean(BoutDur_mins),
              BoutDur_mins_sd = sd(BoutDur_mins),
              BoutDur_mins_95 = quantile(BoutDur_mins, probs = 0.95)
    )

SumTable <- D %>%
    # Remove bouts formed by a single upcall
    filter(!is.element(Bout, db1$Bout)) %>%
    group_by(Year) %>%
    summarize(n_calls = n(),
              ICI_95_min = quantile(ICI[-1], probs = 0.95, na.rm = TRUE) / 60
    ) %>%
    left_join(dtmp)

knitr::kable(SumTable, digits = 2)


## Bout plots ----

# Plot bout duration distributions
p1 <- DB %>%
    ggplot(aes(x = BoutDur_mins)) +
    geom_histogram(fill = "blue", alpha = 0.5, color = NA) +
    labs(x = "Bout duration (minutes)",
         y = "Frequency")
p2 <- DB %>%
    ggplot(aes(x = BoutDur_mins)) +
    geom_histogram(fill = "blue", alpha = 0.5, color = NA) +
    scale_x_log10() +
    labs(x = "Bout duration (minutes)",
         y = "Frequency")
p1 + p2 +
    plot_annotation(tag_levels = "A")

ggsave("images/bout_durs.png", width = 8, height = 4, dpi = 600)

# Set the response variables and predictors
RESPONSES <- c("BoutDur_mins_log10")

# PREDICTORS <- c("TimeSinceLVP"
#                 ,"SPLcalling"
#                 ,"SPLvessel"
#                 ,"DoY"
#                 ,"HourDecimal"
# )

PREDICTORS_cycle <- c("TimeSinceLVP"
                      ,"SPLcalling"
                      ,"SPLvessel"
                      ,"DoYcos", "DoYsin"
                      ,"Hoursin", "Hourcos"
)

# Response variables are the first rows/columns in the matrix
DB %>%
    dplyr::select(all_of(c(RESPONSES, PREDICTORS_cycle))) %>%
    GGally::ggpairs() +
    ggplot2::theme(
        axis.text = ggplot2::element_text(size = 8),
        strip.text = ggplot2::element_text(size = 8)
    )

ggsave("images/bout_scatter.png", width = 10, height = 8, dpi = 600)

# Show correlations and p-values
DB %>%
    dplyr::select(all_of(c(RESPONSES, PREDICTORS_cycle))) %>%
    as.matrix() %>%
    Hmisc::rcorr()

# Show 95% confidence intervals for correlations of Response with Predictors
sapply(PREDICTORS_cycle, function(x)
    cor.test(unlist(DB[,RESPONSES]), unlist(DB[,x]), conf.level = 0.95)$conf.int
)


## Bout models ----

# Formula relating the response(s) to explanatory variables
fml <- as.formula(paste(RESPONSES,
                        paste(PREDICTORS_cycle, collapse = " + "),
                        sep = " ~ "))

### Linear regression ----

# Full model (response and all predictors)
M_LR_BoutDur <- lm(fml, data = DB)

# Check collinearity of predictors
ols_coll_diag(M_LR_BoutDur)

# Model selected
m_LR_BoutDur <- ols_step_both_p(M_LR_BoutDur, p_remove = 0.05)
summary(m_LR_BoutDur$model)

# Format the output
tmp <- summary(m_LR_BoutDur$model)
knitr::kable(tmp$coefficients, digits = 3)

# Diagnostics
ols_plot_diagnostics(m_LR_BoutDur$model)
acf(m_LR_BoutDur$model$residuals)

par(mfrow = c(2, 2))
plot(m_LR_BoutDur$model)


with(DB,
     cor.test(SPLvessel, SPLcalling)
     )


### Random forest ----

set.seed(123)
M_RF_BoutDur <- randomForest(fml
                             ,ntree = NTREE, mtry = MTRY, nodesize = NSIZE
                             ,sampsize = ceiling(nrow(DB) * 2/3)
                             ,importance = TRUE
                             ,data = DB)
M_RF_BoutDur
plot(M_RF_BoutDur)

importance_values <- M_RF_BoutDur$importance[,1]
importance_df <- tibble(Variable = names(importance_values),
                        Importance = importance_values) %>%
    arrange(Importance)

ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "",
         y = "Importance")
ggsave("images/importance_BoutDur.png", width = 4, height = 3, dpi = 600)

# Get PDP data
pdp_BoutDur <- get_pdp_data(M_RF_BoutDur, DB, PREDICTORS_cycle)

# Convert the predictor column to factor for faceting
pdp_BoutDur$predictor <- factor(pdp_BoutDur$predictor, levels = PREDICTORS_cycle)

# Plot PDPs with fixed y-axes
ggplot(pdp_BoutDur, aes(x = x, y = yhat)) +
    geom_line() +
    facet_wrap(~predictor, scales = "free_x", nrow = 1) +
    labs(x = "Predictor values",
         y = expression(log[10] * "(Bout duration [minutes])"))
ggsave("images/pdp_BoutDur.png", width = 9, height = 2.5, dpi = 600)



# Occurrence analysis ----

# Save numeric summary of the data as a table
summary_D <- psych::describe(D) %>%
    as.data.frame() %>%
    dplyr::select(n, min, mean, median, max, sd, range)

# Show summary to match with the table in the paper
summary_D
write.csv(summary_D, 'dataderived/summary_D.csv')

## Year-season counts ----

# Count number of upcalls from the year of deployment
table(D$YearDeployment)

# Counts per season
table(D$YearDeployment, D$Season)

# % per season
round(table(D$YearDeployment, D$Season) * 100 / as.vector(table(D$YearDeployment)), 0)


## Seasonal plots ----

# See
# https://stackoverflow.com/questions/15575713/modifying-timezone-of-a-posixct-object-without-changing-the-display
# for help changing the time zone without changing the displayed time
DateMonitoring <- read_xlsx("data/Datetime_Local_NARWMonitoring.xlsx") %>%
    arrange(datetime_local) %>%
    # edit timezone
    mutate(datetime_local = as.POSIXct(as.character(datetime_local),
                                       origin = as.POSIXct("1970-01-01"),
                                       tz = "America/New_York")) %>%
    # round to hours
    mutate(datetime_local = lubridate::floor_date(datetime_local, unit = "hour")) %>%
    # keep only unique
    dplyr::distinct(datetime_local, .keep_all = TRUE)

attr(DateMonitoring$datetime_local, "tzone") <- "UTC"

DateMonitoring_day <- DateMonitoring %>%
    # subtract 5 hours to convert to "UTC - 5"
    mutate(datetime_local = datetime_local - 5*3600) %>%
    mutate(datetime_local = as.Date(datetime_local)) %>%
    dplyr::distinct(datetime_local, .keep_all = TRUE)

AllDays <- tibble(date = seq(from = as.Date(min(DateMonitoring$datetime_local)),
                             to = as.Date(max(DateMonitoring$datetime_local)),
                             by = "1 day"),
                  lat = 38.303,
                  lon = -74.645)

sun <- getSunlightTimes(data = AllDays,
                        tz = "UTC", #"America/New_York",
                        keep = c("sunrise", "sunset")) %>%
    # subtract 5 hours to convert to "UTC - 5"
    mutate(sunrise = sunrise - 5*3600,
           sunset = sunset - 5*3600) %>%
    mutate(sunrise_h = get_decimal_hour(sunrise),
           sunset_h = get_decimal_hour(sunset)) %>%
    arrange(date)

# Create a copy of the main data but switch times to UTC - 5
D_utc <- D
attr(D_utc$Time_NYC, "tzone") <- "UTC"
D_utc <- D_utc %>%
    mutate(Time_NYC = Time_NYC - 5*3600) %>%
    mutate(HourDecimal = get_decimal_hour(Time_NYC),
           Date = as.Date(Time_NYC),
           Hour = as.numeric(format(Time_NYC, "%H")))

# Year 1
DateMonitoring_day_0 <- DateMonitoring_day %>%
    dplyr::filter(datetime_local < BuoyReplacement) %>%
    dplyr::filter(datetime_local < as.Date("2022-04-01")) %>%
    dplyr::filter(datetime_local > as.Date("2021-10-20"))
sun_0 <- sun %>%
    dplyr::filter(date < BuoyReplacement) %>%
    dplyr::filter(date < as.Date("2022-04-01")) %>%
    dplyr::filter(date > as.Date("2021-10-20"))

d_0 <- D_utc %>%
    dplyr::filter(Date < BuoyReplacement)

p1 <- ggplot() +
    geom_rect(data = DateMonitoring_day_0,
              aes(ymin = datetime_local, ymax = datetime_local + days(1),
                  xmin = -Inf, xmax = Inf),
              fill = "gray", alpha = 0.5) +
    geom_point(data = sun_0, pch = 16,
               aes(x = sunrise_h, y = date), color = "orange", size = 1) +
    geom_point(data = sun_0, pch = 16,
               aes(x = sunset_h, y = date), color = "orange", size = 1) +
    geom_point(data = d_0,
               aes(x = HourDecimal, y = Date),
               color = "black", size = 1, pch = 15) +
    labs(x = "Hour of day (UTC - 5)",
         y = "") +
    scale_x_continuous(limits = c(0, 24),
                       breaks = seq(0, 24, by = 6),
                       minor_breaks = seq(0, 24, by = 2),
                       expand = c(0.01, 0.01)) +
    scale_y_date(breaks = as.Date(c("2021-11-01", "2022-01-01",
                                    "2022-04-01",
                                    "2023-01-01", "2023-04-01")),
                 date_labels = "%Y %b",
                 # date_minor_breaks = "1 month"
                 minor_breaks = seq(from = as.Date("2021-11-01"),
                                    to = as.Date("2023-04-01"),
                                    by = "1 month"),
                 expand = c(0.01, 0.01)
    )

# Year 2
DateMonitoring_day_0 <- DateMonitoring_day %>%
    dplyr::filter(datetime_local > BuoyReplacement) %>%
    dplyr::filter(datetime_local > as.Date("2022-10-20"))
sun_0 <- sun %>%
    dplyr::filter(date > BuoyReplacement) %>%
    dplyr::filter(date > as.Date("2022-10-20"))
d_0 <- D_utc %>%
    dplyr::filter(Date > BuoyReplacement)

p2 <- ggplot() +
    geom_rect(data = DateMonitoring_day_0,
              aes(ymin = datetime_local, ymax = datetime_local + days(1),
                  xmin = -Inf, xmax = Inf),
              fill = "gray", alpha = 0.5) +
    geom_point(data = sun_0, pch = 16,
               aes(x = sunrise_h, y = date), color = "orange", size = 1) +
    geom_point(data = sun_0, pch = 16,
               aes(x = sunset_h, y = date), color = "orange", size = 1) +
    geom_point(data = d_0,
               aes(x = HourDecimal, y = Date),
               color = "black", size = 1, pch = 15) +
    labs(x = "Hour of day (UTC - 5)",
         y = "") +
    scale_x_continuous(limits = c(0, 24),
                       breaks = seq(0, 24, by = 6),
                       minor_breaks = seq(0, 24, by = 2),
                       expand = c(0.01, 0.01)) +
    scale_y_date(breaks = as.Date(c("2021-11-01", "2022-01-01",
                                    "2022-11-01",
                                    "2023-01-01", "2023-04-01")),
                 date_labels = "%Y %b",
                 # date_minor_breaks = "1 month"
                 minor_breaks = seq(from = as.Date("2021-11-01"),
                                    to = as.Date("2023-04-01"),
                                    by = "1 month"),
                 expand = c(0.01, 0.01)
    )

p1 + p2 +
    plot_annotation(tag_levels = "A")

ggsave("images/occurrence_diel.png", width = 8, height = 4, dpi = 600)


## Diel analysis ----

# Count calls per date and hour, then combine with all monitoring days
D_dh <- D_utc %>%
    dplyr::group_by(Date, Hour) %>%
    dplyr::summarise(Ncalls = n())
# Check
sum(D_dh$Ncalls) == nrow(D)

Ncalls0 <- expand.grid(Date = DateMonitoring_day$datetime_local,
                       Hour = 0:23) %>%
    as_tibble() %>%
    left_join(D_dh, by = c("Date", "Hour")) %>%
    mutate(HourFac = factor(Hour, levels = 0:23),
           Year = "year 1",
           HourSin = sin(2 * pi * Hour / 24),
           HourCos = cos(2 * pi * Hour / 24)) %>%
    arrange(Date, Hour)
# Those counts of calls that are missing during the monitoring are 0
Ncalls0$Ncalls[is.na(Ncalls0$Ncalls)] <- 0
# Check
sum(Ncalls0$Ncalls) == nrow(D)

# Plot number of calls per hour
# Aggregate calls by hour
Ncalls0_agg <- Ncalls0 %>%
    group_by(Hour) %>%
    summarise(Ncalls = sum(Ncalls))
# Check
sum(Ncalls0_agg$Ncalls) == nrow(D)

p3 <- Ncalls0_agg %>%
    ggplot(aes(x = Hour, y = Ncalls)) +
    geom_bar(stat = "identity") +
    labs(x = "Hour of day (UTC - 5)", y = "Number of calls")

p1 + p2 + p3 +
    plot_annotation(tag_levels = "A")

ggsave("images/occurrence_dielNcalls.png", width = 8, height = 3, dpi = 600)


# Use a model to test the diel pattern
# Model 1: based on categorical hours -- has higher AIC and BIC than model 2 below.
# m_ncalls_gls <- nlme::gls(log(Ncalls + 1) ~ HourFac
#           ,correlation = nlme::corAR1()
#           ,data = Ncalls0)

# Model 2: based on sin+cos transformation of hours
m_ncalls_gls <- nlme::gls(log(Ncalls + 1) ~ HourSin + HourCos
                          ,correlation = nlme::corAR1()
                          ,method = "REML"
                          ,data = Ncalls0)

summary(m_ncalls_gls)
# ?nlme:::residuals.gls
# e <- residuals(m_ncalls_gls, type = "pearson")
e <- residuals(m_ncalls_gls, type = "normalized")

# Residual diagnostics plots
par(mfrow = c(2, 2))
plot.ts(e,
        las = 1,
        ylab = "Residuals")
plot(x = fitted(m_ncalls_gls), y = e,
     las = 1,
     ylab = "Residuals", xlab = "Fitted")
plot(x = Ncalls0$Hour, y = e,
     las = 1,
     ylab = "Residuals", xlab = "Hour of day (UTC - 5)")
acf(e,
    las = 1,
    main = "Autocorrelation of residuals")

# Test significance of diel patterns, using ML since the fixed effects are different
m_ncalls_gls <- nlme::gls(log(Ncalls + 1) ~ HourSin + HourCos
                          ,correlation = nlme::corAR1()
                          ,method = "ML"
                          ,data = Ncalls0)

# Reduced model (without predictors)
m_ncalls_gls_reduced <- gls(log(Ncalls + 1) ~ 1
                            ,correlation = corAR1()
                            ,method = "ML"
                            ,data = Ncalls0)

# Likelihood ratio test
# Use the p-value to determine if the predictors significantly improve the model.
anova(m_ncalls_gls, m_ncalls_gls_reduced)

# 95% confidence interval for chi-squared distribution under the null hypothesis
df <- 2 # df.full - df.reduced
c(qchisq(0.025, df), qchisq(0.975, df))



# Call rates ----

## Calculate call rates ----

# Aggregate to calculate call rates
DC <- D %>%
    group_by(Date, Month, DoY, Hour, Season) %>%
    summarise(CallRate = n(),
              HourDecimal = mean(HourDecimal, na.rm = TRUE),
              TimeSinceLVP = mean(TimeSinceLVP),
              SPLcalling = mean(SPLcalling),
              SPLvessel = mean(SPLvessel),
              DoYcos = mean(DoYcos),
              DoYsin = mean(DoYsin),
              Hoursin = mean(Hoursin),
              Hourcos = mean(Hourcos)) %>%
    ungroup() %>%
    mutate(CallRate_log10 = log10(CallRate))

# Save numeric summary of the data as a table
summary_DC <- psych::describe(DC) %>%
    as.data.frame() %>%
    dplyr::select(n, min, mean, median, max, sd, range)

# Show summary to match with the table in the paper
summary_DC
write.csv(summary_DC, 'dataderived/summary_DC.csv')


## Call rates plots ----

# Plot bout duration distributions
p1 <- DC %>%
    ggplot(aes(x = CallRate)) +
    geom_histogram(fill = "blue", alpha = 0.5, color = NA, binwidth = 2) +
    labs(x = "Call rate (number of calls)",
         y = "Frequency")
p2 <- DC %>%
    ggplot(aes(x = CallRate)) +
    geom_histogram(fill = "blue", alpha = 0.5, color = NA) +
    scale_x_log10() +
    labs(x = "Call rate (number of calls)",
         y = "Frequency")
p1 + p2 +
    plot_annotation(tag_levels = "A")

ggsave("images/CallRate.png", width = 8, height = 4, dpi = 600)


# Set the response variables and predictors
RESPONSES <- c("CallRate_log10")

# Response variables are the first rows/columns in the matrix
DC %>%
    dplyr::select(all_of(c(RESPONSES, PREDICTORS_cycle))) %>%
    GGally::ggpairs() +
    ggplot2::theme(
        axis.text = ggplot2::element_text(size = 8),
        strip.text = ggplot2::element_text(size = 8)
    )

ggsave("images/CallRate_scatter.png", width = 10, height = 8, dpi = 600)

# Show correlations and p-values
DC %>%
    dplyr::select(all_of(c(RESPONSES, PREDICTORS_cycle))) %>%
    as.matrix() %>%
    Hmisc::rcorr()

# Show 95% confidence intervals for correlations of Response with Predictors
sapply(PREDICTORS_cycle, function(x)
    cor.test(unlist(DC[,RESPONSES]), unlist(DC[,x]), conf.level = 0.95)$conf.int
)


p1 <- DC %>%
    ggplot(aes(x = SPLcalling, y = CallRate)) +
    geom_point() +
    scale_y_log10() +
    labs(y = "Call rate (number of calls)",
         x = "SPL of calling period (dB)")
p2 <- DB %>%
    ggplot(aes(x = SPLcalling, y = BoutDur_mins)) +
    geom_point() +
    scale_y_log10() +
    labs(y = "Bout duration (minutes)",
         x = "SPL of calling period (dB)")
p3 <- DB %>%
    ggplot(aes(x = SPLvessel, y = BoutDur_mins)) +
    geom_point() +
    scale_y_log10() +
    labs(y = "Bout duration (minutes)",
         x = "SPL of vessel passage (dB)")
p1 + p2 + p3 +
    plot_annotation(tag_levels = "A")
ggsave("images/CallRateBout.png", width = 8, height = 3, dpi = 600)


## Call rate models ----

# Formula relating the response(s) to explanatory variables
fml <- as.formula(paste(RESPONSES,
                        paste(PREDICTORS_cycle, collapse = " + "),
                        sep = " ~ "))

### Linear regression ----

# Full model (response and all predictors)
M_LR_CallRate <- lm(fml, data = DC)

# Check collinearity of predictors
ols_coll_diag(M_LR_CallRate)

# Model selected
m_LR_CallRate <- ols_step_both_p(M_LR_CallRate, p_remove = 0.05)
summary(m_LR_CallRate$model)

# Format the output
tmp <- summary(m_LR_CallRate$model)
knitr::kable(tmp$coefficients, digits = 3)

# Diagnostics
ols_plot_diagnostics(m_LR_CallRate$model)
acf(m_LR_CallRate$model$residuals)

par(mfrow = c(2, 2))
plot(m_LR_CallRate$model)


### Random forest ----

set.seed(123)
M_RF_CallRate <- randomForest(fml
                              ,ntree = NTREE, mtry = MTRY, nodesize = NSIZE
                              ,sampsize = ceiling(nrow(DC) * 2/3)
                              ,importance = TRUE
                              ,data = DC)
M_RF_CallRate
plot(M_RF_CallRate)

importance_values <- M_RF_CallRate$importance[,1]
importance_df <- tibble(Variable = names(importance_values),
                        Importance = importance_values) %>%
    arrange(Importance)

ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(x = "",
         y = "Importance")
ggsave("images/importance_CallRate.png", width = 4, height = 3, dpi = 600)

# Get PDP data
pdp_CallRate <- get_pdp_data(M_RF_CallRate, DC, PREDICTORS_cycle)

# Convert the predictor column to factor for faceting
pdp_CallRate$predictor <- factor(pdp_CallRate$predictor, levels = PREDICTORS_cycle)

# Plot PDPs with fixed y-axes
ggplot(pdp_CallRate, aes(x = x, y = yhat)) +
    geom_line() +
    facet_wrap(~predictor, scales = "free_x", nrow = 1) +
    labs(x = "Predictor values",
         y = expression(log[10] * "(Number of calls per hour)"))
ggsave("images/pdp_CallRate.png", width = 9, height = 2.5, dpi = 600)



# Call characteristics ----

## Plots ----

# Set the response variables and predictors
RESPONSES <- c("FreqMin", "FreqMax", "FreqDelta", "Duration", "Duration90")

# Response variables are the first rows/columns in the matrix
D %>%
    dplyr::select(all_of(c(RESPONSES, PREDICTORS_cycle))) %>%
    GGally::ggpairs() +
    ggplot2::theme(
        axis.text = ggplot2::element_text(size = 8),
        strip.text = ggplot2::element_text(size = 8)
    )

ggsave("images/CallChar_scatter.png", width = 10, height = 8, dpi = 600)

# Show correlations and p-values
D %>%
    dplyr::select(all_of(c(RESPONSES, PREDICTORS_cycle))) %>%
    as.matrix() %>%
    Hmisc::rcorr()

# Show 95% confidence intervals for correlations of Response with Predictors
tmp <- lapply(RESPONSES, function(r)
    sapply(PREDICTORS_cycle, function(x)
        cor.test(unlist(D[,r]), unlist(D[,x]), conf.level = 0.95)$conf.int
    )
)
names(tmp) <- RESPONSES
tmp


# Plot boxplots of ALL responses per season


ps <- lapply(RESPONSES, function(response) {
    ggplot(D, aes(x = .data$Season, y = .data[[response]])) +
        geom_boxplot(fill = "bisque") +
        labs(x = "Season", y = response)
})

# List sample sizes on one plot
ps[[1]] <- ps[[1]] +
    stat_summary(
        fun.data = function(x) {
            return(data.frame(
                y = min(x) - 0.1,
                label = paste0("n = ", length(x))
            ))
        },
        geom = "text",
        size = 4,
        color = "blue"
    )

p2 <- ggplot(DC, aes(x = Season, y = CallRate)) +
    geom_boxplot(fill = "darkolivegreen3") +
    labs(x = "Season", y = "Call rate") +
    stat_summary(
        fun.data = function(x) {
            return(data.frame(
                y = min(x) - 1,
                label = paste0("n = ", length(x))
            ))
        },
        geom = "text",
        size = 4,
        color = "blue"
    )

p3 <- ggplot(DB, aes(x = Season, y = BoutDur_mins)) +
    geom_boxplot(fill = "cornflowerblue") +
    labs(x = "Season", y = "Bout duration") +
    stat_summary(
        fun.data = function(x) {
            return(data.frame(
                y = min(x) - 1,
                label = paste0("n = ", length(x))
            ))
        },
        geom = "text",
        size = 4,
        color = "blue"
    )

# Combine all ggplot objects using patchwork
combined_plot <- wrap_plots(ps) + plot_spacer() + p2 + p3

# Print the combined plot
print(combined_plot) +
    plot_layout(axes = "collect") +
    plot_annotation(tag_levels = "A")

ggsave("images/Response_boxplots.png", width = 9, height = 7, dpi = 600)


## Multivariate linear regression ----

# Set formula object
mult_form_lin <- as.formula(
    paste("cbind(",
          paste0(RESPONSES, collapse = ", "),
          ") ~",
          paste0(PREDICTORS_cycle, collapse = " + ")
    )
)

# Multivariate model of linear parametric form
mult_mod_lin <- lm(mult_form_lin
                   ,data = D)

# The function summary() gives a separate summary for each response variable
summary(mult_mod_lin)

# Run multivariate tests and check if thy agree
car::Anova(mult_mod_lin, test.statistic = "Wilks")
car::Anova(mult_mod_lin, test.statistic = "Pillai")
car::Anova(mult_mod_lin, test.statistic = "Hotelling-Lawley")
car::Anova(mult_mod_lin, test.statistic = "Roy")



## Multiple (individual) linear regressions ----

# Run model selection and save results in a list
invisible(
    lapply(RESPONSES, function(y) { # y = RESPONSES[1]
        print(y)
        fml <- as.formula(paste(y,
                                paste(PREDICTORS_cycle, collapse = " + "),
                                sep = " ~ "))

        # Full model (response and all predictors)
        Mfull <- lm(fml, data = D)

        # Model selected
        Mreduced <- ols_step_both_p(Mfull, p_remove = 0.05)
        print(summary(Mreduced$model))

        # Format the output
        tmp <- summary(Mreduced$model)
        print(knitr::kable(tmp$coefficients, digits = 3))

        # Diagnostics
        par(mfrow = c(2, 2))
        plot(Mreduced$model)
    })
)

## Random forest ----

# Set formula object
mult_form_rf <- as.formula(
    paste("Multivar(",
          paste0(RESPONSES, collapse = ", "),
          ") ~",
          paste0(PREDICTORS_cycle, collapse = " + ")
    )
)

# Multivariate model of nonlinear nonparametric form
set.seed(123)
mult_mod_rf <- rfsrc(mult_form_rf,
                     ntree = NTREE,
                     nodesize = NSIZE,
                     importance = "permute",
                     splitrule = "mahalanobis",
                     data = D)

# Print out RF results for each response
invisible(
    sapply(RESPONSES, function(i) print(mult_mod_rf, outcome.target = i))
)

# Plot RF results for each response
invisible(
    sapply(RESPONSES, function(i) plot(mult_mod_rf, m.target = i))
)

# Extract standardized variable importance
svimp <- get.mv.vimp(mult_mod_rf, standardize = TRUE) %>%
    t() %>%
    as_tibble()

# Calculate average importance for sorting of the variables
svimp_avg <- apply(svimp, 2, mean) %>%
    sort()
svimp_long <- svimp %>%
    pivot_longer(cols = PREDICTORS_cycle,
                 names_to = "Variable",
                 values_to = "Value") %>%
    mutate(Variable = factor(Variable, levels = names(svimp_avg)))

# Plot importances
ggplot(svimp_long, aes(x = Value, y = Variable)) +
    xlim(0, max(svimp)) +
    geom_boxplot(fill = "bisque") +
    xlab("Standardized importance") +
    ylab("")
ggsave("images/importance_CallChar.png", width = 8, height = 5, dpi = 600)


# Get PDP data
pdp_CallChar <- lapply(PREDICTORS_cycle, function(i) { # i = "TimeSinceLVP"

    # Setup x-grid similar to pdp::partial and pdp:::pred_grid
    grid.resolution <- min(length(unique(mult_mod_rf$xvar[, i, drop = TRUE])), 51)
    xx <- seq(from = min(mult_mod_rf$xvar[, i, drop = TRUE], na.rm = TRUE),
              to = max(mult_mod_rf$xvar[, i, drop = TRUE], na.rm = TRUE),
              length = grid.resolution)

    # Use new partial.values xx instead of the mult_mod_rf$xvar[, i]
    partial.obj <- randomForestSRC::partial(mult_mod_rf,
                                            partial.xvar = i,
                                            partial.values = xx)

    pdta <- lapply(RESPONSES, function(j)
        randomForestSRC::get.partial.plot.data(partial.obj, m.target = j))

    lapply(1:length(RESPONSES), function(j)
        tibble(x = pdta[[j]]$x,
               yhat = pdta[[j]]$yhat,
               predictor = i,
               response = RESPONSES[j])) %>%
        bind_rows()
}) %>%
    bind_rows()

# Rename stuff to include units
pdp_CallChar <- pdp_CallChar %>%
    mutate(predictor = gsub("TimeSinceLVP", "TimeSinceLVP (hours)", predictor)) %>%
    mutate(predictor = gsub("SPLcalling", "SPLcalling (dB)", predictor)) %>%
    mutate(predictor = gsub("SPLvessel", "SPLvessel (dB)", predictor)) %>%
    mutate(response = if_else(grepl("Duration", response),
                               paste0(response, " (s)"),
                               response)) %>%
    mutate(response = gsub("FreqDelta", "FreqDelta (Hz)", response)) %>%
    mutate(response = gsub("FreqMax", "FreqMax (Hz)", response)) %>%
    mutate(response = gsub("FreqMin", "FreqMin (Hz)", response))


# Plot PDPs with fixed y-axes
pdp_CallChar %>%
    dplyr::filter(predictor %in% c("TimeSinceLVP (hours)", "SPLcalling (dB)", "SPLvessel (dB)")) %>%
    ggplot(aes(x = x, y = yhat)) +
    geom_line() +
    facet_grid(response ~ predictor, scales = "free") +
    labs(x = "Predictor values",
         y = "Response values")
ggsave("images/pdp_CallChar1.png", width = 9, height = 8, dpi = 600)


pdp_CallChar %>%
    dplyr::filter(predictor %!in% c("TimeSinceLVP (hours)", "SPLcalling (dB)", "SPLvessel (dB)")) %>%
    ggplot(aes(x = x, y = yhat)) +
    geom_line() +
    facet_grid(response ~ predictor, scales = "free") +
    labs(x = "Predictor values",
         y = "Response values")
ggsave("images/pdp_CallChar2.png", width = 9, height = 8, dpi = 600)



save.image("dataderived/DataImage_modeling_combined.RData")


# Citations ----

citation(package = "car", auto = TRUE)
citation(package = "flexmix")
citation(package = "nlme")
citation(package = "randomForestSRC", auto = TRUE)
citation(package = "randomForest", auto = TRUE)
