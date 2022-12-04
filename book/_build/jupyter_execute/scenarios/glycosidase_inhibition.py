#!/usr/bin/env python
# coding: utf-8

# # __Scenario B:__<br>$\alpha$-glucosidase inhibition by fucoidan
# 
# Data provided by Chantal Daub (Biochemistry, Rhodes University, Makhanda, South Africa)
# 
# ## Project background
# In this scenario, the inhibitory properties of fucoidan, a polysaccharide found in brown algae, on $\alpha$-glucosidase from *Saccharomyces cerevisiae* was investigated. Fucoidans are actively investigated in the fields of anit-cancer, anti-inflammation, and anti-coagulate, to name a fiew ({cite:t}`li2008fucoidan`). In a previous study ({cite:t}`daub2020fucoidan`), fucoidan from *E. maxima* showed an almost 2-fold lower IC<sub>50</sub> value, compared to acarbose. Thus, fucoidan is a potential antidiabetic drug candidate for the treatment of diabetes mellitus.
# In the following analysis, fucoidan from the brown algae species *Ecklonia maxima*, *Ecklonia radiata*, *Fucus vesiculosus*, S. CYM???????, and *Schimmelmannia elegans* were investigated for their inhibition constant $K_{i}$ for $\alpha$-glucosidase inhibition.
# 
# ## Experimental design 
# 
# Extracted fucoidan from the mentioned brown algea species was applied in two different concentrations to enzyme reactions. Additionally, control reactions without inhibitor were performed. For the enzyme reactions, *p*-nitrophenyl-$\alpha$-D-glucopyranoside (pNPG) was applied as a substrate in a concentration range from 0.1 - 5 mM. Product accumulatin was followed photometrically in a micro titer plate at 405 nm for 20 min.
# 
# ## Data preparation
# 
# ### Imports

# In[1]:


from typing import Dict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import re
import os
import pyenzyme as pe
from CaliPytion.tools.standardcurve import StandardCurve
from EnzymePynetics.tools.parameterestimator import ParameterEstimator


import warnings
warnings.filterwarnings('ignore')

colors = list(mcolors.TABLEAU_COLORS.values())


# ### Experimental data
# Time-course data of all kinetic experiments was collected in an Excel file, whereas the meta data was specified in individual EnzymeML Excel templates for each origin species of fucoidan. In the following cell the experimental data is loaded and blanked to subtract the absorbance contribution of enzyme, buffer, and the respective inhibitor for each measurement. Thereafter, the data is written to EnzymeML documents.

# In[2]:


dataset_path = "../../data/glucosidase_inhibition/experimental_data.xlsx"
template_directory = "../../data/glucosidase_inhibition/EnzymeML_templates"   

inhibitors = sorted(pd.ExcelFile(dataset_path).sheet_names)
initial_substrates = [0.1, 0.25, 0.5, 1, 2.5, 5]


data_dict = {}
for inhibitor in inhibitors:
    df = pd.read_excel(dataset_path, sheet_name=inhibitor).set_index("time")
    blanc_no_inhibitor = df["buffer + enzyme"].mean()
    blank_low_inhibitor = df[df.columns[[1,2]]].values.mean()
    blank_high_inhibitor = df[df.columns[[3,4]]].values.mean()
    df = df.iloc[:,17:]
    keys = sorted(df.columns)
    df_no_inhibitor = df[keys[:12]]
    df_low_inhibitor = df[keys[12:24]]
    df_high_inhibitor = df[keys[24:]]

    df_no_inhibitor = df_no_inhibitor.subtract(blanc_no_inhibitor)
    df_low_inhibitor = df_low_inhibitor.subtract(blank_low_inhibitor)
    df_high_inhibitor = df_high_inhibitor.subtract(blank_high_inhibitor)


    data = []
    data.append(df_no_inhibitor.values.T)
    data.append(df_low_inhibitor.values.T)
    data.append(df_high_inhibitor.values.T)

    data_dict[inhibitor] = np.array(data).reshape(18,2,20)

time = df.index.values

# Parse measurement data to EnzymeML documents
def measurement_data_to_EnzymeML(
    template_path: str,
    measurement_data: np.ndarray,
    species_id: str,
    time: np.ndarray,
    data_unit: str,
    time_unit: str
    ) -> pe.EnzymeMLDocument:

    enzmldoc: pe.EnzymeMLDocument = pe.EnzymeMLDocument.fromTemplate(template_path)

    for IDs, concentration in zip(enzmldoc.measurement_dict.keys(), measurement_data):
        for counter, replicate in enumerate(concentration):
            
            enzmldoc.getMeasurement(IDs).addReplicates(pe.Replicate(
                id=f"Measurement{counter}",
                species_id=species_id,
                data=list(replicate),
                data_unit=data_unit,
                time=list(time),
                time_unit=time_unit), enzmldoc)

    return enzmldoc

enzml_docs = []
datas = list(data_dict.values())
for file, data in zip(sorted(os.listdir(template_directory)), datas):
    enzml_docs.append(measurement_data_to_EnzymeML(
        template_path=f"{template_directory}/{file}",
        measurement_data=data,
        time=time,
        species_id="s1",
        data_unit="mmole / l",
        time_unit="min"
    ))


# ### Concentration calculation
# 
# Standard data of the product was loaded from an excel file, and a standard curve was created. Then, the standard curve was applied to the EnzymeML documents, containing the absorption measurements for concentration calculation.

# In[3]:


path_calibration_data = "../../data/glucosidase_inhibition/p-NP_standard.xlsx"


product_standard = StandardCurve.from_excel(
    path=path_calibration_data,
    reactant_id="s1", 
    sheet_name="csv", 
    wavelength=405, 
    concentration_unit = "mmole / l", 
    cutoff_absorption=2)

product_standard.visualize()

# Apply calibration curves to absorption EnzymeML documents
for enzmldoc in enzml_docs:
    product_standard.apply_to_EnzymeML(enzmldoc, "s1")


# ### Parameter estimation
# 
# Parameter estimation was perfomed with EnzymePynetics. Thereby, each data set was fitted to competitive, uncompetitive, and non-competitive inhibition models.

# In[14]:


kinetics = []
for enzmldoc in enzml_docs:
    result = ParameterEstimator.from_EnzymeML(enzmldoc=enzmldoc,reactant_id="s1", inhibitor_id="s2", measured_species="product")
    result.fit_models()
    kinetics.append(result)
    result.visualize(plot_means=True)
    plt.show()


# ## $k_{cat}$ and $K_{m}$
# 
# 

# In[13]:


# Get kinetic parameters of all datasets
kcat = []
kcat_std = []
Km = []
Km_std = []
corr_kcat_km = []
for result in kinetics:
    params = result.get_parameter_dict()

    kcat.append(params["k_cat"].value)
    kcat_std.append(params["k_cat"].stderr)

    Km.append(params["Km"].value)
    Km_std.append(params["Km"].stderr)


    correlation = params["k_cat"].correl
    if correlation == None:
        corr_kcat_km.append(float("nan"))
    else:
        corr_kcat_km.append(correlation["Km"])


df = pd.DataFrame.from_dict({
    'kcat [1/min]':kcat, 
    'kcat stderr':kcat_std, 
    'Km [mM]':Km, 
    'Km stderr':Km_std, 
    "correlation kcat/Km":corr_kcat_km})

df


# In[5]:


enzmldoc = enzml_docs[0]

del enzmldoc.measurement_dict["m6"]
del enzmldoc.measurement_dict["m7"]
del enzmldoc.measurement_dict["m8"]
del enzmldoc.measurement_dict["m9"]
del enzmldoc.measurement_dict["m10"]
del enzmldoc.measurement_dict["m11"]
del enzmldoc.measurement_dict["m12"]
del enzmldoc.measurement_dict["m13"]
del enzmldoc.measurement_dict["m14"]
del enzmldoc.measurement_dict["m15"]
del enzmldoc.measurement_dict["m16"]
del enzmldoc.measurement_dict["m17"]


# In[6]:


kinetics = ParameterEstimator.from_EnzymeML(
    enzmldoc=enzmldoc,
    reactant_id="s1",
    measured_species="product")

kinetics.fit_models(enzyme_inactivation=True)


# In[7]:


kinetics.fit_models(
    enzyme_inactivation=True, 
    start_time_index=8,
    initial_substrate_concs=[0.1, 0.25, 0.5, 1, 5])


# In[8]:


kinetics.visualize(model_name="irreversible Michaelis Menten", title="a-glucosidase reaction")

