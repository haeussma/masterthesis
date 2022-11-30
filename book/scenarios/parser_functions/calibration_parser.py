from typing import Dict, List

import re
import numpy as np

from .utility import read_photometer, tab_split, string_to_float, blockshaped

from CaliPytion.core.calibration import Calibration
from CaliPytion.core.device import Device
from CaliPytion.core.spectrum import Spectrum
from CaliPytion.core.standard import Standard
from CaliPytion.core.series import Series

def read_calibration_data(
    path_standard: str,
    path_spectrum: str,
    species_id,
    wavelengths,
    concentrations: List[float],
    concentration_unit,
    device_manufacturer: str,
    device_model: str,
    spectrum_reactant_concentration: float,
    ) -> Dict[float, Calibration]:

    # Define parameters for parser
    temperature_unit = 'C'

    # Read and parse
    df_standard = read_photometer(path_standard)
    
    # Format df_standard and extract metadata
    metadata_row = str(df_standard.iloc[-1].values)
    doc_name = path_standard.split("/")[-1]
    date = re.findall(r'\d{4}/\d{2}/\d{2}', metadata_row)[0]
    pH = re.findall(r'\d\.\d\+\d\.\d', metadata_row)[0]
    pH = [float(x) for x in pH.split("+")]
    temperature = float(re.findall(r'\d\d', doc_name)[0].split(" ")[0])

    df_standard = df_standard.loc[13:18]
    df_standard = df_standard.reset_index().drop(columns="index")

    # array manipulation to get measurement data
    standard_absorption = np.array([])
    standard_data = df_standard.values
    for i, row in enumerate(standard_data):
        standard_data = str(row).split("\\t")
        standard_data = [string_to_float(x) for x in standard_data]

        if i == 0:
            standard_data = standard_data[1:]
            standard_absorption = np.append(standard_absorption, standard_data)
        else:
            standard_absorption = np.append(standard_absorption, standard_data)
    
    standard_absorption = standard_absorption[~np.isnan(standard_absorption)]
    standard_absorption = standard_absorption.reshape(6,24)

    # Reshape into (pH, wavelength, replicates, concentrations)
    standard_absorption = blockshaped(standard_absorption, 3, 12).reshape(2,2,3,12)

    # Map pH and corresponding mesurement data into dict
    standard_data_dict = dict(zip(pH, standard_absorption))

    # Spectrum 
    df_spectrum = read_photometer(path_spectrum)

    metadata_row = str(df_spectrum.iloc[-1].values)
    date = re.findall(r'\d{4}/\d{2}/\d{2}', metadata_row)[0]
    pH_spectrum = re.findall(r'\d\.\d\+\d\.\d', metadata_row)[0]
    pH_spectrum = [float(x) for x in pH_spectrum.split("+")]

    if pH_spectrum != pH:
        raise ValueError(
            f"pH of files do not match. Check following files:\n{path_standard}\n{path_spectrum}")

    df_spectrum = df_spectrum.loc[13:87]
    df_spectrum = df_spectrum.reset_index().drop(columns="index")
    df_spectrum = df_spectrum.applymap(tab_split)

    spectrum_data = df_spectrum.values.tolist()
    for i, wavelength in enumerate(spectrum_data):
        spectrum_data[i] = list(filter(None, wavelength[0]))

    data = np.array(spectrum_data)
    spectrum_wavelengths = data[:,0]
    spectrum_absorption = data[:,2:].T

    spectrum_data_dict = dict(zip(pH_spectrum, spectrum_absorption))


    ### Map data to calibration datamodel instance ### 
    # Define photometer meta-data
    device = Device(
        manufacturer=device_manufacturer,
        model=device_model)

    # Create instance of calibration data-model
    # Since each instance of the datamodel only deals with one pH, two instances are created in this case.
    instance_dict = {}
    for pH_value, absorbance_data in standard_data_dict.items():

        # Create instance of data model
        instance = Calibration(
            reactant_id=species_id,
            pH=pH_value,
            date=date,
            device=device,
            temperature=temperature,
            temperature_unit=temperature_unit)

        # Create instances of 'Standard' for both wavelengths
        standard_curves = []
        for i, (wavelength, data) in enumerate(zip(wavelengths, absorbance_data)):
            
            absorption_series = []
            for replicate in absorbance_data[i]:
                absorption_series.append(Series(
                    values=list(replicate)))
            instance.add_to_standard(
                concentration=concentrations,
                wavelength=wavelength,
                concentration_unit=concentration_unit,
                absorption=absorption_series
            )

            standard_curve = Standard(
                    wavelength=wavelength,
                    concentration=concentrations,
                    concentration_unit=concentration_unit
                )
            # Add measurerment data of the standards
            for replicate in data:
                standard_curve.add_to_absorption(list(replicate))

            standard_curves.append(standard_curve)

        # Add absorption spectrum
        spectrum = Spectrum(
            concentration=spectrum_reactant_concentration,
            concentration_unit=concentration_unit,
            wavelength=spectrum_wavelengths.tolist())

        spectrum.add_to_absorption(spectrum_data_dict[pH_value].tolist())

        instance.spectrum = spectrum
            
        # Add sub-classes of the datamodel to the corresponding parent attributes of the datamodel
        

        # Manage the instances in for of a dict
        instance_dict[pH_value] = instance

    return instance_dict

if __name__ == "__main__":
    result = read_calibration_data(
        path_standard="/Users/maxhaussler/Dropbox/master_thesis/data/sdRDM_ABTS_oxidation/StandardData/pH 3.0+3.5 25deg standards.txt",
        path_spectrum="/Users/maxhaussler/Dropbox/master_thesis/data/sdRDM_ABTS_oxidation/SpectrumData/pH 3.0+3.5 25deg scan.txt",
        species_id="s0",
        wavelengths=[340, 420],
        concentrations=[0,5,10,15,25,50,75,100,125,150,175,200],
        concentration_unit="umole / l",
        device_manufacturer="MANUFACTURER",
        device_model="SUPERMODEL",
        spectrum_reactant_concentration=69
    )
