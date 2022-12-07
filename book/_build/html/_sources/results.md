# Results

## Project-driven workflow development

The development of the workflow was driven by different research scenarios of EnzymeML project partners. The goal of all project was to reliably estimate kinetic parameters of enzyme reactions. The collaboration consisted of multiple rounds of lab experiments and kinetic modeling of the respective data. Each round consisted of (i) wet lab experiments, (ii) kinetic modeling, (iii) design of follow-up experiments, and (iv) discussion of results with project partners. Hence, a short feedback loop between lab experiments and modeling-based experimental suggestions was established.  
In parallel, the python packages [CaliPytion](method:calipytion) and [EnzymePynetics](method:enzymepynetics) were developed to provide a workflow from analytical raw data to kinetic parameter estimates. Thereby, the individual requirements of each research scenario, fostered the implementation of different features. Hence, a generic workflow for parameter estimation of enzyme kinetics was established. The developed workflow is schematically visualized in figure 1.

![Fig. 1](images/concept_workflow.png)
_Fig. 1: Workflow for kinetic parameter estimation of enzyme reactions._

The workflow was designed around EnzymeML, whereas it's data model serves