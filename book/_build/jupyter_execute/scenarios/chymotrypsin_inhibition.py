#!/usr/bin/env python
# coding: utf-8

# # __Scenario A:__<br>Chymotrypsin inhibiton by a designed albumin fusion protein
# 
# Data provided by Marwa Mohamed (Institute of Cell Biology and Immunology, University of Stuttgart, Stuttgart, Germany)
# 
# ## Project background
# In this scenaro, the inhibitory effect of a human serum albumin mutant on chymotrypsin was investigated. Thereby, the inhibition constant $K_{i}$ of the HSA wild-type was compared to HSA(M3) mutant. Both HSAs were individually dimerized into fusion proteins through a huFc.  
# Experimental data from the HSA(wt)-huFc and HSA(M3)-huFc originate from two independent experiments. In each experiment the initial substrate concentration was varied. Once with the respective inhibitor, and once without. Therefore, $K_{i}$ and $K_{m}$ were determined indenpendently from each other. Additionally, each individual reaction condition was prepared in duplicates to ensure repeatability.  
# Since 
# 
# ### Experimental design
# 
# In order to assess the effect of the introduced mutations to the HSA(M3) variant, $K_{i}$ of the HSA wild-type was compared to HSA(M3) variant. Enzyme activity was monitored by measuring the product formation of p-Nitroanilin (p-NA) photometrically at 410 nm for 30 min at 30°C. Therefore,  Succinyl-gly-gly-phe-p-nitroanilide (SGGPpNA) was applied as substrate in a concentration range of 0.25 - 2 mM. For concentration calculations, p-NA standard was prepared in the range of 0 - 0.3 mM in duplicates. $K_{i}$ of HSA(wt) and HSA(M3) on chymotrypsin were investigated in independent experiments. Each experiment consisted of enzyme reactions with and without the respective HSA variant. The enzyme reactions contined dfhdfhdfh µM of enzyme and 26.88 µM of the respective HSA variant, if applied.
# 
# ### Data management
# 
# Experimental data and meta data was filled in EnzymeML Excel templates for each of the two inhibition experiments respectively. Calibration data was stored as Excel files. Measurement data was already blanked.
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
from IPython.display import display
from EnzymePynetics.tools.parameterestimator import ParameterEstimator
from CaliPytion.tools.standardcurve import StandardCurve

import warnings
warnings.filterwarnings('ignore')


# ### Concentration calculation
# 
# Product standard data was imported directly from an Excel file. Then, a standard curve was created.

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


# Based on the Akaike information criterion (AIC), the relation between concentration and absorption is best described by a quadratic function. The figure above visualizes the calibration measurements and the fitted calibration model.

# ### Experimantal data
# 
# The EnzymeML documents were loaded and the standard curve was applied to the absorption data for concentration calculation.

# In[3]:


# Load data from 
chymo_HSAwt = pe.EnzymeMLDocument.fromTemplate("../../data/chymotrypsin_inhibition/chymo_HSAwt.xlsx")
chymo_HSAM3 = pe.EnzymeMLDocument.fromTemplate("../../data/chymotrypsin_inhibition/chymo_HSA(M3).xlsx")

# Apply standard curve to 'EnzymeMLDocument'
chymo_HSAwt = product_standard.apply_to_EnzymeML(chymo_HSAwt, "s1")
chymo_HSAM3 = product_standard.apply_to_EnzymeML(chymo_HSAM3, "s1")


# ## Comparability of the experiments
# 
# Since the experimental data from the HSA wild-type and the HSA(M3) variant originate from independent experiments, the control reactions without the respective inhibitor were compared by performing a parameter estimation. Thereby, $k_{cat}$ and $K_{m}$ were highly correlated (corr > 0.98). Hence, catalytic efficiency $\frac{k_{cat}}{K_{m}}$ was used to assess comparability between the data sets. High correlation indicates, that the highest initial substrate concentration is too low, compared to the true $K_{m}$ of the enzyme under the given experimental conditions. In this case, higher substrate concentration were not applied for multiple reasons. On the one hand dimethyl sulfoxide (DMSO) was used as a co-solvent of the substrate, which inhibits enzyme activity {cite:t}`busby1999effect`. Hence, higher initial substrate concentrations would have led to higher enzyme inhibition. On the other hand, high substrate viscosity denied the application of higher concentrations without sacrificing pipetting percision.

# In[4]:


# Create copys of the data sets and delete measuremnts with inhibitor.
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
print("Kinetic parameters of HSA(wt) chymotrypsin control reactions:")
display(kinetics_wt_control.result_dict.drop(columns=["kcat [1/min]", ]))


# Estimate kinetic parameters of the control reactions of the HSA(M3) data set.
kinetics_m3_control = ParameterEstimator.from_EnzymeML(m3_control, "s1", "product")
kinetics_m3_control.fit_models(stop_time_index=-1, display_output=False)
print("\nKinetic parameters of HSA(M3) chymotrypsin control reactions:")
display(kinetics_m3_control.result_dict)


fig, axes = plt.subplots(1,2, figsize=(12.8,4.8), sharey=True, sharex=True)
for e, (doc, ax, title) in enumerate(zip([kinetics_wt_control ,kinetics_m3_control], axes.flatten(), ["chymotrypsin control reactions HSA(wt)", "chymotrypsin control reactions HSA(M3)"])):
    doc.visualize(ax=ax, title=title)
    ax.set_ylabel("4-nitroanilin [mM]")
    ax.set_xlabel("time after reaction start [min]")
    ax.set_xticks([5, 10, 15, 20])

handles, labels = ax.get_legend_handles_labels()

fig.legend(handles, labels, loc="lower center", ncol=4, title="initial SGGPpNA  [mM]", bbox_to_anchor=(0.5,-0.15))
plt.tight_layout()


# _Fig. 1: Measurement data and fitted irreversible Michaelis-Menten model for chymotrypsin reactions without HSA inhibitor._
# 
# The above figure visualizes the control reactions of each HSA inhibition experiment. Each dataset was fitted against the models listed in the above table. 
# Each dataset is best described by the irreversible Michaelis-Menten model in therms of AIC and standard deviation on the estimated parameters. Models with product or substrate inhibition resulted in large uncertanties above 80 % on the parameter estimates. Thus, substrate and product inhibition were ruled out for the given reactions.  $\frac{k_{cat}}{K_{m}}$ was estimated to be 25.353 min<sup>-1</sup>mM<sup>-1</sup> ± 8.57%  for the control reaction of the HSA(wt) data set and 24.776 min<sup>-1</sup>mM<sup>-1</sup> ± 4.35% for the HSA(M3) data set. As a result, the two experiments showed to be comparable, since the catalytic efficiency differs less than 3 % between the two data sets.
# 
# ## Determination and comparison of $K_{i}$
# 
# The data set of chymotrypsin inhibition by HSA(M3) contained negative absorption values for the first measurement point. Presumably from an incorrect blank measurement. Therefore, only measurement data from the second data point (minute 5 and onward) was considered for parameter estimation. Additionally, measuremet
# The tables below show the parameter estimastes for all applied kinetic models.

# In[5]:


# Parameter estimation for HSA(wt) data set
kinetics_HSAwt = ParameterEstimator.from_EnzymeML(chymo_HSAwt, reactant_id="s1", inhibitor_id="s2", measured_species="product")
kinetics_HSAwt.fit_models(initial_substrate_concs=[0.25, 0.5, 1, 2], stop_time_index=-1, start_time_index=5, display_output=False)
print("Kinetic parameters estimates for all models of chymotrypsin inhibition by HSA(wt):")
display(kinetics_HSAwt.result_dict.drop(columns=["kcat [1/min]", "Km [mmole / l]"]))

# Parameter estimation for HSA(M3) data set
kinetics_HSAM3 = ParameterEstimator.from_EnzymeML(chymo_HSAM3, reactant_id="s1", inhibitor_id="s3", measured_species="product")
kinetics_HSAM3.fit_models(initial_substrate_concs=[0.25, 0.5, 1, 2], stop_time_index=-1, start_time_index=1, display_output=False)
print("\nKinetic parameters estimates for all models of chymotrypsin inhibition by HSA(M3):")
display(kinetics_HSAM3.result_dict.drop(columns=["kcat [1/min]", "Km [mmole / l]"]))

# Visualize experimental data and fitted models
fig, axes = plt.subplots(1,2, figsize=(12.8,4.8), sharey=True, sharex=True)
for e, (doc, ax, title) in enumerate(zip([kinetics_HSAwt ,kinetics_HSAM3], axes.flatten(), ["chymotrypsin inhibition by HSA(wt)", "chymotrypsin inhibition by HSA(M3)"])):
    doc.visualize(ax=ax, title=title)
    ax.set_ylabel("4-nitroanilin [mM]")
    ax.set_xlabel("time after reaction start [min]")
    ax.set_xticks([5, 10, 15, 20])

handles, labels = ax.get_legend_handles_labels()

fig.legend(handles, labels, loc="lower center", ncol=2, title="initial SGGPpNA  [mM]", bbox_to_anchor=(0.5,-0.2))
plt.tight_layout()


# _Fig. 1: Measurement data and fitted product inhibition model for chymotrypsin reactions with respective HSA inhibitior._
# 
# Both reaction systems are best described by the competitive inhibition model, which is indicated by the lowest AIC and standard deviation on the estimated parameters. Thereby, a $K_{i}$ of 0.460 mM ± 27.24% was estimated for HSA(wt) and 0.059 mM ± 8.51% for HSA(M3). This resembles a roughly 7-fold increase in affinity of HSA(M3) to the enzyme compared to the HSA(wt). Since the competitive inhibition model describes the data the best, HSA(M3) presumably interactis with the enzyme in the active site region.  
