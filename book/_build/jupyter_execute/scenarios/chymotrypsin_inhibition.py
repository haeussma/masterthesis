#!/usr/bin/env python
# coding: utf-8

# # __Scenario A:__<br>Chymotrypsin inhibiton by a designed albumin fusion protein
# 
# Data provided by Marwa Mohamed (Institute of Cell Biology and Immunology, University of Stuttgart, Stuttgart, Germany)
# 
# Assessment of the inhibitory properties of an designed albumin as a chymotrypsin inhibitor
# 
# # Project background
# In this scenaro, the inhibitory effect of a human serum albumin mutant on chymotrypsin was investigated. Thereby, the inhibition constant $K_{i}$ of the HSA wild-type was compared to HSA(M3) mutant. Both HSAs were fusionized into fusion proteins individually through a huFc.  
# Experimental data from the HSA(wt)-huFc and HSA(M3)-huFc originate from two independent experiments. In each experiment the initial substrate concentration was varied. Once with the respective inhibitor, and once without. Therefore, $K_{i}$ and $K_{m}$ were determined indenpendently from each other. Additionally, each individual reaction condition was prepared in duplicates to ensure repeatability.  
# Since 
# 
# ## Experimental design
# 
# In order to assess the effect of the introduced mutations to the HSA(M3) variant, $K_{i}$ of the HSA wild-type was compared to HSA(M3) variant. Therefore, two indevidual experiments were conducted, estimating $K_{i}$ for each experiment individually. Each experiment consisted of enzyme reactions with initial substrate concentrations in the range of 0.25 - 2 mM, with and without the respective HSA.
# 
# ## Data preparation
# 
# ### Imports

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyenzyme as pe
import copy
from EnzymePynetics.tools.parameterestimator import ParameterEstimator
from CaliPytion.tools.standardcurve import StandardCurve

import warnings
warnings.filterwarnings('ignore')


# ### Conversion of absortion signal to concentration
# 
# In order to convert the time-course absorption signal into concenrtation, a standard of p-NA was applied in the range of 0 - 0.3 mM in duplicates. The measured absorbance values were stored as an excel-file. In the following cell, a CaliPytion ```StandardCurve``` is generated directly from the excel-data.

# In[2]:


product_standard = StandardCurve.from_excel(
    path="../../data/chymotrypsin_inhibition/pNA-standard.xlsx",
    reactant_id="s1",
    wavelength=410,
    sheet_name="csv", 
    concentration_unit="mmole / l",
    temperature=30,
    temperature_unit="C")
    
product_standard.visualize()


# Based on the Akaike information criterion (AIC), the relation between concentration and absorption is best described by a quadratic function. Additionally, the visualization of the fit shows 

# ### Load experimental data
# 
# All experimental data was filled in the EnzymeML Excel template. In the next cell, ```EnzymeMLDocument```s are created by reading the excel template and the measured absorption data is calculated based on the ```StandardCurve``` above.

# In[3]:


# Load data from 
chymo_HSAwt = pe.EnzymeMLDocument.fromTemplate("../../data/chymotrypsin_inhibition/chymo_HSAwt.xlsx")
chymo_HSAM3 = pe.EnzymeMLDocument.fromTemplate("../../data/chymotrypsin_inhibition/chymo_HSA(M3).xlsx")

# Apply standard curve to 'EnzymeMLDocument's
chymo_HSAwt = product_standard.apply_to_EnzymeML(chymo_HSAwt, "s1")
chymo_HSAM3 = product_standard.apply_to_EnzymeML(chymo_HSAM3, "s1")

# Visualize the measurement data
chymo_HSAwt.visualize()
plt.show()


# ## Check comparability of the two data sets
# 
# Since the experimental data from the HSA wild-type and the HSA(M3) variant originate from independent experiments, the control reactions without the respective inhibitor were compared. Therefore, kinetic parameters were estimated for each dataset, after deleting measurement in which inhibitor was present.

# In[14]:


# Create copys of the data sets and delete measuremnts with inhibitor present.
wt_control = copy.deepcopy(chymo_HSAwt)
del wt_control.measurement_dict["m4"]
del wt_control.measurement_dict["m5"]
del wt_control.measurement_dict["m6"]
del wt_control.measurement_dict["m7"]

m3_control = copy.deepcopy(chymo_HSAM3)
del m3_control.measurement_dict["m4"]
del m3_control.measurement_dict["m5"]
del m3_control.measurement_dict["m6"]
del m3_control.measurement_dict["m7"]

# Estimate kinetic parameters of the control reactions of the HSA wild-type data set.
kinetics_wt_control = ParameterEstimator.from_EnzymeML(wt_control, "s1", "product")
kinetics_wt_control.fit_models(stop_time_index=-1, display_output=False)
kinetics_wt_control.visualize()
plt.show()


# In[15]:


# Estimate kinetic parameters of the control reactions of the HSA(M3) data set.
kinetics_m3_control = ParameterEstimator.from_EnzymeML(m3_control, "s1", "product")
kinetics_m3_control.fit_models(stop_time_index=-1, display_output=False)
kinetics_m3_control.visualize()


# Statistical analysis between the parameters $k_{cat}$ and $K_{m}$ high correlation above 0.98. This indicates, that the highest applied substrate concentration of 2 mM is lower than the $K_{m}$ of the enzyme under the given experimental conditions. Therefore, $k_{cat}$ and $K_{m}$ cannot be detemined independently. Instead, the catalytic efficiency $\frac{k_{cat}}{K_{m}}$ is used to assess the comparability of the two data sets.  
# In the cell below, 

# In[7]:


kinetics_wt_control.result_dict.drop(columns=["kcat [1/min]", "Km [mmole / l]"])


# In[8]:


kinetics_m3_control.result_dict.drop(columns=["kcat [1/min]", "Km [mmole / l]"])


# 

# ## Parameter estimation for $K_{i}$

# In[9]:


kinetics_HSAwt = ParameterEstimator.from_EnzymeML(chymo_HSAwt, reactant_id="s1", inhibitor_id="s2", measured_species="product")
kinetics_HSAwt.fit_models(initial_substrate_concs=[0.25, 0.5, 1], stop_time_index=-1, start_time_index=5, display_output=False)
kinetics_HSAwt.visualize()


# In[10]:


kinetics_HSAM3 = ParameterEstimator.from_EnzymeML(chymo_HSAM3, reactant_id="s1", inhibitor_id="s3", measured_species="product")
kinetics_HSAM3.fit_models(initial_substrate_concs=[0.25, 0.5, 1], stop_time_index=-1, start_time_index=1, display_output=False)
kinetics_HSAM3.visualize()


# In[11]:


kinetics_HSAwt.result_dict.drop(columns=["kcat [1/min]", "Km [mmole / l]"])


# In[12]:


kinetics_HSAM3.result_dict.drop(columns=["kcat [1/min]", "Km [mmole / l]"])


# In[13]:


fig, axes = plt.subplots(1,2, figsize=(15,5), sharey=True, sharex=True)
for e, (doc, ax) in enumerate(zip([kinetics_HSAwt ,kinetics_HSAM3], axes.flatten())):
    doc.visualize(ax=ax)
    ax.set_ylabel("4-nitroanilin [mM]")
    ax.set_xlabel("time after reaction start [min]")
    ax.set_xticks([5, 10, 15, 20])

handles, labels = ax.get_legend_handles_labels()

fig.legend(handles, labels, loc="lower center", ncol=2, title="initial SGGPpNA  [mM]", bbox_to_anchor=(0.5,-0.15))
plt.show()

