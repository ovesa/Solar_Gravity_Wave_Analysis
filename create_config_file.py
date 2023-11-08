## @oanavesa 08/25/22
## This code is designed to create a config.ini file
## Intended to be used to quickly grab information used in analysis

################################################################


################### Import Necessary Modules ###################
from configparser import ConfigParser

################### Creating the Config File ###################

# Get the configparser object
config_object = ConfigParser()


# Create the sections needed
config_object["04252019"] = {
    "solarx": "95.9425782",
    "solary": "23.6802246",
    "dx": "0.6",
    "dt": "11.88",
    "dx_units": "arcsec/pixel",
    "dt_units": "s",
}

config_object["04262019"] = {
    "solarx": "782.408788800",
    "solary": "56.5665528000",
    "dx": "0.6",
    "dt": "16.4030",
    "dx_units": "arcsec/pixel",
    "dt_units": "s",
    "dt_ca_fe7090": "4.394",
    "dt_ca_k7699": "7.366",
    "dt_ca_fe5434": "11.616",
    "dt_fe7090_k7699": "2.971",
    "dt_fe7090_fe5434": "7.227",
    "dt_k7699_fe5434": "4.25",
}

config_object["05062019"] = {
    "solarx": "-708.744140625",
    "solary": "180.312744141",
    "dx": "0.6",
    "dx_units": "arcsec/pixel",
    "dt_units": "s",
    "dt": "16.386",
    "dt_ca_fe7090": "4.39",
    "dt_ca_k7699": "7.36",
    "dt_ca_fe5434": "11.613",
    "dt_fe7090_k7699": "2.97",
    "dt_fe7090_fe5434": "7.22",
    "dt_k7699_fe5434": "4.25",
}

config_object["05092019"] = {
    "solarx": "170.3333497",
    "solary": "208.9144045",
    "dx": "0.6",
    "dx_units": "arcsec/pixel",
    "dt_units": "s",
    "dt": "16.383",
    "dt_ca_fe7090": "3.822",
    "dt_ca_k7699": "7.361",
    "dt_ca_fe5434": "11.328",
    "dt_fe7090_k7699": "3.540",
    "dt_fe7090_fe5434": "7.506",
    "dt_k7699_fe5434": "3.966",
}

# Line Centers for 09May2019
# Ca II 8542.1104
# Fe I 7090.3521
# K I
# Fe I 5434.4824


# Write the above sections to config.ini file
with open("AGWs_config.ini", "w") as conf:
    config_object.write(conf)
