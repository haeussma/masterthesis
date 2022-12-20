# Results

## Workflow for enzyme kinetics parameter estimation

In this work, a Python-based analysis workflow for kinetic parameter estimation of enzyme reactions was developed. The workflow itself consists of the steps (i) data acquisition, (ii) raw data preparation, (iii) concentration calculation, (iv) quality control, (v) modeling, and (vi) saving of results. For precise concentration calculations was well as kinetic parameter estimation the two Python modules CaliPytion and EnzymePynetics were respectively developed over the course of this work.
The two packages were designed to seamlessly interact with EnzymeML documents which were modified using PyEnzyme and thus allowing an continuous data stream from analytical raw data to kinetic parameters.  
Starting from raw data of photometric measurements, which are stored in an EnzymeML document, concentration values are calculated by CaliPytion, whereupon kinetic parameters are determined by EnzymePynetics. Finally, the modeling results are written back to the EnzymeML document via PyEnzyme (Fig. 1).

![Fig. 1](images/concept_workflow.png)
_Fig. 1: Components of the kinetic parameter estimation workflow for enzyme reaction data._

A key aspect of the workflow design was to enable reproducible analysis. Therefore, documentation of raw data treatment prior to modeling is vital. As a result, Jupyter Notebooks were used as a platform, since they allow to write text and code within the same document. In consequence, data analysis as well as the documentation of data analysis are unified within the same document. Ultimately, each data treatment step from raw data to kinetic parameters is comprehensively documented.

xxx#TODO check necessity of packages here

### CaliPytion

CaliPytion was developed to facilitate robust conversion of measurement data to concentration data. Based on an analyte standard, a calibration curve is created. Since linear and non-linear relations between analyte concentration and analytic signal may occur, {cite}`hsu2010effect`, {cite}`martin2017fitting`, linear, polynomial, exponential, and rational calibration equations were implemented. After fitting, the best model is preselected based on Akaike information criterion {cite}`akaike1998information`, offering a metric to compare models with different number of parameters. Thereby, AIC describes the information loss relative to alternative models, while penalizing the use of additional parameters. CaliPytion was designed to work with PyEnzyme. Hence, a `StandardCurve` of CaliPytion can natively by applied to an `EnzymeMLDocument` of PyEnzyme. As a result, absorption data of an EnzymeML document can be converted into concentration data based on the previously generated standard curve.

### EnzymePynetics

EnzymePynetics was developed to enable easy applicable parameter estimation of single substrate enzyme reaction directly from EnzymeML documents. In it's core, the `ParameterEstimator` of EnzymePynetics fits experimental data to different Michaels-Menten type models. Besides the irreversible Michaelis-Menten model, competitive, non-competitive, uncompetitive, and partially competitive inhibition models were implemented. Hence, the inhibition constant of potentially inhibiting substrate or product can be estimated. Furthermore, the inhibition constant of an applied inhibitor can be quantified. Since, the experimental data is fitted against multiple kinetic models, the best fitting model is suggested based on AIC. Hence, a proposal of the enzyme's kinetic mechanism is given based on the best fitting model.
For quality control, visualization of experimental data as well as the fitted kinetic model was a priority during development. Thereby, detection of systematic deviation between model and measurement data is possible. In consequence, potential experimental errors can be identified. Furthermore, the `ParameterEstimator` allows to subset the measurement data by time as well as initial substrate concentration without deleting data. Hence, identified systematic deviations (e.g. lag-phases or measurements with incorrect enzyme concentrations) can be excluded from parameter estimation.  
Like the `StandardCurve` of CaliPytion, the `ParameterEstimator` of EnzymePynetics natively supports `EnzymeMLDocument`s as a data source.

## Scenario-driven workflow development

The development of the workflow was driven by different research scenarios of EnzymeML project partners. Goal of all projects was to reliably determine kinetic parameters of enzyme reactions.
The collaboration consisted of multiple rounds of lab experiments (preformed by the experimental partners) and kinetic modeling of the experimental data. Each round included (i) wet lab experiments, (ii) kinetic modeling, (iii) design of follow-up experiments, and (iv) discussion of results with the respective project partners. Hence, a short feedback loop between lab experiments and modeling-based experimental suggestions was established.  
In parallel, the Python modules CaliPytion and EnzymePynetics were developed. Thereby, individual requirements of each research scenario fostered the implementation of various functionalities, resulting in a feature-rich yet generic workflow.

## FAIR data analysis with Jupyter Notebooks

All data analysis of the project was carried out in Jupyter Notebook by applying the developed workflow for kinetic parameter estimation. Within the next chapters, the latest results of all four experimental scenarios are shown. Thereby, the applicability of the workflow is demonstrated while analyzing and discussing the results of the respective scenario.
Since the analysis was conducted in Jupyter Notebooks, each of the following chapter is an executed Jupyter Notebook. Hence all figures and tables for data visualization are generated at runtime of the analysis.
Thereby, the applicability of Jupyter Notebooks for reproducible and comprehensive computational workflows in an scientific environment is demonstrated.
Every scenario consists of a short description of the project's background and methodology carried out by the project partners in the wet lab. Additionally, data preparation, kinetic modeling steps, as well as the results are shown and project specific results are discussed.
All notebooks can be launched interactively by clicking on the {fa}`rocket`-icon in the upper section of each scenario page. As a result, the analysis itself is transparent, repeatable and compliant with FAIR data principles. Therefore, it is highly recommended to [read this work in its native form](https://haeussma.github.io/masterthesis/welcome.html).