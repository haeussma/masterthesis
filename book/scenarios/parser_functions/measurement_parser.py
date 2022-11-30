from typing import Dict, List
import re
import numpy as np
import pandas as pd
from datetime import datetime

from .utility import read_photometer, tab_split, string_to_float, to_seconds


def read_measurement_data(path: str) -> Dict:

    # Read raw data, parse metadata
    df = read_photometer(path)

    metadata_row = path.split("/")[-1]
    date = re.findall(r'\d{4}-\d{2}-\d{2}', metadata_row)[0]
    pH = float(re.findall(r'\d\.\d', metadata_row)[0])
    temperature = float(re.findall(r'\d{2} degrees', metadata_row)[0].split(" ")[0])
    
    df = df.iloc[13:-9].reset_index().drop(columns="index")
    df = df.applymap(tab_split)
    initial_substrates=[0,5,10,15,25,50,75,100,150,200]
    time = []
    array = []

    # Extract measurement data
    for index, row in df.iterrows():
        if index % 6 == 0:
            time.append(row.values[0][0])
            df.loc[index, "##BLOCKS= 2"] = row.values[0][2:]
    
        data = row.values[0]
        array.append([string_to_float(x) for x in data])

    # Cleaning and restructuring of data
    array: np.ndarray = np.array(array)
    array = array[~np.isnan(array)]
    array = array.reshape(11,6,20)
    array = array.swapaxes(0,2)

    #(wavelength, concentration, control, replicates, data)
    array = array.reshape((2,10,2,3,11))
    
    #(wavelength, control, concentration, replicates, data)
    array = array.swapaxes(1,2)

    substrate = array[0][0]
    substrate_control = array[0][1]
    product = array[1][0]
    product_control = array[1][1]

    measurements_product_control = []
    for data, init_substrate in zip(product_control, initial_substrates):
        measurement = {}
        measurement["initial_substrate"] = 0
        reps = []
        for replicate in data:
            rep_dict = {}
            rep_dict["replicate"] = replicate.tolist()
            reps.append(rep_dict)
        measurement["data"] = reps
        
        measurements_product_control.append(measurement)

    measurements_substrate_control = []
    for data, init_substrate in zip(substrate_control, initial_substrates):
        measurement = {}
        measurement["initial_substrate"] = init_substrate
        reps = []
        for replicate in data:
            rep_dict = {}
            rep_dict["replicate"] = replicate.tolist()
            reps.append(rep_dict)
        measurement["data"] = reps
        
        measurements_substrate_control.append(measurement)

    measurements_substrate = []
    for data, init_substrate in zip(substrate, initial_substrates):
        measurement = {}
        measurement["initial_substrate"] = init_substrate
        reps = []
        for replicate in data:
            rep_dict = {}
            rep_dict["replicate"] = replicate.tolist()
            reps.append(rep_dict)
        measurement["data"] = reps
        
        measurements_substrate.append(measurement)

    measurements_product = []
    for data in product:
        measurement = {}
        measurement["initial_substrate"] = 0
        reps = []
        for replicate in data:
            rep_dict = {}
            rep_dict["replicate"] = replicate.tolist()
            reps.append(rep_dict)
        measurement["data"] = reps
        
        measurements_product.append(measurement)

    data_dict = {
        "pH": pH,
        "name": f"ABTS oxidation pH {pH} and {temperature}°C",
        "date": str(datetime(*[int(x) for x in date.split("-")],)),
        "time": [to_seconds(x) for x in time],
        "temperature": temperature,
        "s0": measurements_substrate,
        "s1": measurements_product,
        "s2": measurements_substrate_control,
        "s3": measurements_product_control,
    }
    return data_dict


if __name__ == "__main__":
    import os 
    names = []

    directory_measurement_data = "/Users/maxhaussler/Dropbox/master_thesis/data/sdRDM_ABTS_oxidation/TimeCourseData"
    directory_standrd_data = "StandardData"
    directory_spectrum_data = "SpectrumData"

    # Parse measurement data from photometer output
    names = []
    raw_data_dict = {}
    sorted_list = np.sort(os.listdir(directory_measurement_data))
    for path in sorted_list:
        data = read_measurement_data(f"{directory_measurement_data}/{path}")
        pH = data["pH"]
        temp = data["temperature"]
        raw_data_dict[f"{pH} {temp}"] = data
        names.append(f"{pH} {temp}")