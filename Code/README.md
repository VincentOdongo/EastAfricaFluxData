Main Script
carbon_water_dynamics_East_Africa.m: This is the main script for the project. It orchestrates the data processing, calls various functions for analysis, and generates the outputs. The script integrates data readings, applies statistical analyses, and plots the results.

Supporting Functions
These functions are used by the main script to handle specific tasks such as data loading and processing:

KE_AUS_cropland_biomet.m: Function to read and preprocess biometeorological data from the Australian cropland site.
KE_AUS_flux_cropland.m: Function for processing flux data specific to the cropland in Australia, aiding in the calculation of net ecosystem exchange.
KE_KAP_flux_rangeland.m: Dedicated to processing flux data from the Kapiti rangeland site. It handles specific data transformations necessary for the analysis.
KE_KAP_rangeland_biomet.m: Reads and prepares biometeorological data from the Kapiti rangeland, setting up the data for further analysis by the main script.
