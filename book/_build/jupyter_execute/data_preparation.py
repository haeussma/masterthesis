#!/usr/bin/env python
# coding: utf-8

# # Test book
# 

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyenzyme as pe


# In[2]:


df = pd.read_excel("Daten von 3 Experimenten.xlsx", sheet_name="control_reactions").set_index("Unnamed: 0")
experiment_ids = df.loc["Experiment"].values
df = df.drop(index="Experiment")
df


# In[3]:


data = df.values.T.reshape(9,2,22)
labels = [0.25, 0.5, 1, 2, 0.25, 0.5, 1, 2, 0.5]
colors = ["orange", "blue", "green", "red", "orange", "blue", "green", "red", "blue", ]
markers = ["o", "o", "o", "o", "d", "d", "d", "d", "x", ]

for d, label, color, marker in zip(data, labels, colors, markers):
    mean =np.mean(d, axis=0)
    std =np.std(d, axis=0)
    time = df.index.values

    plt.scatter(time, mean, color=color, label=f"{label} mM", marker=marker)
plt.legend(title="initial substrate")
plt.title("Comparison of chymotrypsin control reactions across three experiment")
plt.ylabel("absorption")
plt.xlabel("time [min]")
plt.show()


# In[4]:


import pyenzymekinetics as pek


# In[5]:


enzmldoc = pe.EnzymeMLDocument.fromTemplate("chymo_HSAwt and HSA(M3).xlsx")
#enzmldoc.visualize()


# In[6]:


product_standard = pek.StandardCurve.from_excel("../round5_control/pNA-standard.xlsx", sheet_name="csv", concentration_unit="mM", substance_name="p-NP")
#product_standard.visualize()


# In[7]:


enzmldoc = pek.abso_to_conc(enzmldoc, product_standard, reactant_id="s1")


# In[8]:


#kin_two_inhibitiors = pek.EnzymeKinetics.from_EnzymeML(enzmldoc, inhibitor_id="s2", inhibitor2_id="s3")
#kin_two_inhibitiors.fit_models(initial_substrates=[0.25, 0.5, 1], start_time=1)
#kin_two_inhibitiors.visualize_fit()


# In[9]:


#enzmldoc.toFile("", name="two_inhibitors")


# ## Seperate analysis

# In[10]:


chymo_HSAwt = pe.EnzymeMLDocument.fromTemplate("/Users/maxhaussler/Dropbox/master_thesis/data/marwa/final/chymo_HSAwt.xlsx")
chymo_HSAM3 = pe.EnzymeMLDocument.fromTemplate("/Users/maxhaussler/Dropbox/master_thesis/data/marwa/final/chymo_HSA(M3).xlsx")

chymo_HSAwt = pek.abso_to_conc(chymo_HSAwt, product_standard, reactant_id="s1")
chymo_HSAM3 = pek.abso_to_conc(chymo_HSAM3, product_standard, reactant_id="s1")


# In[11]:


kin_HSAwt = pek.EnzymeKinetics.from_EnzymeML(chymo_HSAwt, inhibitor_id="s2")
kin_HSAwt.fit_models(initial_substrates=[0.25, 0.5, 1], start_time=1)
kin_HSAwt.visualize_fit("competitive inhibition")


# In[12]:


kin_HSAM3 = pek.EnzymeKinetics.from_EnzymeML(chymo_HSAM3, inhibitor_id="s3")
kin_HSAM3.fit_models(initial_substrates=[0.25, 0.5, 1], start_time=1)
kin_HSAM3.visualize_fit("competitive inhibition")


# In[13]:


kin_HSAM3.visualize_fit()

