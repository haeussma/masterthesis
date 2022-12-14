# Results

## Python packages

In this work, two Python modules were developed, which form the components of an analysis workflow from raw data to kinetic parameter estimates for enzyme reactions. Thereby, measurement data was initially stored in an EnzymeML document, whereupon the raw data was converted into concentration data with CaliPytion. Thereafter, kinetic parameters were estimated with EnzymePynetics. Finally, the parameter estimates were written back to the EnzymeML document together with information on the respective kinetic model.  
CaliPytion was developed to facilitate robust conversion of raw measurement data to concentration data. Based on a standard of the analyte a calibration curve is created, describing the respective signal to concentration relationship. Since the relation can be linear or non-liner, depending on the utilized detection range of the instrument as well as solvent properties {cite}`hsu2010effect`, {cite}`martin2017fitting`, linear, polynomial, exponential, and rational calibration equations were implemented. Model selection was based on Akaike information criterion {cite}`akaike1998information`, offering a metric to compare models with different number of parameters.

Beginning from measurement data of an analytical device, an EnzymeML domument

## Scenario-driven workflow development

The development of the workflow was driven by different research scenarios of EnzymeML project partners. The goal of all project was to reliably estimate kinetic parameters of enzyme reactions. The collaboration consisted of multiple rounds of lab experiments and kinetic modeling of the respective data. Each round consisted of (i) wet lab experiments, (ii) kinetic modeling, (iii) design of follow-up experiments, and (iv) discussion of results with project partners. Hence, a short feedback loop between lab experiments and modeling-based experimental suggestions was established.  
In parallel, the python modules [CaliPytion](method:calipytion) and [EnzymePynetics](method:enzymepynetics) were developed to provide a workflow from analytical raw data to kinetic parameter estimates. Thereby, the individual requirements of each research scenario fostered the implementation of different features. Hence, a generic workflow for parameter estimation of enzyme kinetics was established. The developed workflow is schematically visualized in figure 1.

![Fig. 1](images/concept_workflow.png)
_Fig. 1: Workflow for kinetic parameter estimation of enzyme reactions._

The workflow was designed around EnzymeML, whereas it's data model serves

- Data model vs. format
- FAIR

The following chapters show the developed workflow applied in different research scenarios. Each each of the following sections is an executable JupyterNotebook. Hence all figures for data visualization were generated at runtime of the analysis.
Each notebook consists of a short description of the project's background and methodology carried out by the project partners in the wet lab. Additionally kinetic modeling steps as well as the results are shown. Lastly, project specific results are discussed.
hg
