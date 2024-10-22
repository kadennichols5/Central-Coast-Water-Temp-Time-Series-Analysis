
## Central Coast Ocean Time Series Analysis

This project analyzes water temperature data from the Santa Barbara Channel, focusing on the years 2020-2022. The goal is to create a forecasting model for future water temperatures and investigate any seasonality or cyclical patterns in the data.

### Project Structure

The analysis is conducted in R and uses various libraries for time series analysis and visualization. The main components of the project are:

1. Data preprocessing
2. Exploratory data analysis
3. Time series decomposition
4. Model selection and fitting
5. Forecasting

### Data Source

Data is collected from NOAA NDBC Buoy LLNR196, located 12 Nautical Miles Southwest of Santa Barbara, CA. The original dataset includes various measurements taken every 10 minutes, but this analysis focuses solely on water temperature.

### Key Features

- Aggregation of data into weekly averages
- Conversion of temperature from Celsius to Fahrenheit
- Time series visualization
- Normality checks
- ACF and PACF analysis
- Differencing and Box-Cox transformations
- Variance comparisons

### Dependencies

The project relies on several R packages, including:
tidyverse (includes ggplot2, dplyr, and more)
readr 
zoo
ggfortify
forecast
lubridate
MASS
tseries


### Results

The analysis provides insights into:

- Seasonal patterns in water temperature
- Autocorrelation structures in the time series
- Effects of various transformations on the data
- Potential models for forecasting

### Future Work

- Implement and compare different forecasting models
- Incorporate additional variables for a more comprehensive analysis
- Extend the analysis to cover a longer time period

### Author

Kaden Nichols

Sources:
[1] [National Data Buoy Center](https://www.ndbc.noaa.gov/station_page.php?station=46053)

