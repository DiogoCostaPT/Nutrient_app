
'''Created by
        Riya Maria Kurian on May 28th, 2019 '''


from urllib.parse import urlsplit
from urllib import request
import sys
import urllib.request
from math import *
import json

def download_stack(Input_URL):
    data = json.loads(Input_URL)
    
    with open('data.txt', 'w') as Input_URL:
        json.dump(data, Input_URL)

    with open('data.txt') as data_file:
        data_loaded = json.load(data_file)

    a = int(input("Enter the measurement you want: ")) #Measurement wanted
    Measurement_wanted = a-1

    Sampling_location = data_loaded['data']['samples'][Measurement_wanted]['at_sampling_location']
    print("At_Sampling Location : " + str(Sampling_location))

    Color_correction = data_loaded['data']['samples'][Measurement_wanted]['color_correction']
    print("Color Correction = " + Color_correction)

    Color_method = data_loaded['data']['samples'][Measurement_wanted]['color_method']
    print("Color Method = " + str(Color_method))

    Date = data_loaded['data']['samples'][Measurement_wanted]['date']
    print("Date : " + Date)

    Error_message = data_loaded['data']['samples'][Measurement_wanted]['error_message']
    print("Error Message : " + Error_message)

    Label = data_loaded['data']['samples'][Measurement_wanted]['label']
    print("Label : " + Label)

    Latitude = data_loaded['data']['samples'][Measurement_wanted]['latitude']
    print("Latitude : " + Latitude)

    Longitude = data_loaded['data']['samples'][Measurement_wanted]['longitude']
    print("Longitude : " + Longitude)

    Light_condition = data_loaded['data']['samples'][Measurement_wanted]['light_condition']
    print("Light condition = " + str(Light_condition))

    Mass_concentration = data_loaded['data']['samples'][Measurement_wanted]['mass_concentration']
    print("Mass Concentration = " + str(Mass_concentration))

    Mass_concentration_uncorrected = data_loaded['data']['samples'][Measurement_wanted]['mass_concentration_uncorrected']
    print("Mass Concentration Uncorrected = " + str(Mass_concentration_uncorrected))

    Measurement = data_loaded['data']['samples'][Measurement_wanted]['measurement']
    print("Measurements = " + str(Measurement))

    Origin_of_water = data_loaded['data']['samples'][Measurement_wanted]['origin_of_water']
    print("Origin of Water : " + Origin_of_water)

    Reference_concentration = data_loaded['data']['samples'][Measurement_wanted]['reference_concentrations']
    print("Reference Concenration = " + str(Reference_concentration))

    Reference_measurement = data_loaded['data']['samples'][Measurement_wanted]['reference_measurements']
    print("Reference Measurement = " + str(Reference_measurement))

    Temperature = data_loaded['data']['samples'][Measurement_wanted]['temperature'].encode('cp1252', errors='ignore').decode('utf-8')
    print("Temperature = " + Temperature)
   

Test = input("Which test results you want(Nitrate or Phosphate: ")
Nitrate = urllib.request.urlopen("https://gwf-nutrient.usask.ca/api/v1/samples?type=nitrate").read()
Phosphate = urllib.request.urlopen("https://gwf-nutrient.usask.ca/api/v1/samples?type=phosphate").read()
if (Test.lower()== "nitrate"):
    download_stack(Nitrate)
elif (Test.lower() == "phosphate"):
    download_stack(Phosphate)
else:
    print("Invalid Entry")

