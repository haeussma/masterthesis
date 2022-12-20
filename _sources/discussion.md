# Discussion

Availability of raw data is the prerequisite to reproduce results of an data analysis process and additionally reinforces confidence in the validity of scientific findings.
Hence, the analysis workflow was conceptualized starting from photometric measurement data. For that purpose, discussing the output format of the analytical device with the experimental partners proved to be valuable to plan the data acquisition process. Hence, a project specific parser function was established, which transferred the measurement data to the respective EnzymeML containing information on the measurement conditions. Thereby, error-prone and tedious manual copying between files was avoided, ensuring raw data integrity from the very start of the workflow.

## Data integrity through accurate concentration calculation

Another aspect of data integrity is the accurate concentration calculation based on the measured absorption signal. The relation between concentration to absorption signal of an analyte can be linear as well as non-liner. Therefore, wrong assumption of a linear relationship albeit the relation is better described non-linearly, leads to avoidable inaccuracy in concentration calculation {cite}`hsu2010effect` {cite}`martin2017fitting` and thus imprecise kinetic parameter estimates.
Hence, concentration calculation based on extinction coefficients was discouraged, if the underlying relation between absorption and concentration of an analyte was not self-determined for the respective reaction system.
To enable precise concentration calculation of linear and non-linear absorption to concentration relationships, CaliPytion was applied with analyte standards provided by the experimental partners.
All of the analyzed standard curves, except for scenario C ABTS at pH 5, were best described by non-linear calibration equations based on their fit quality determined by AIC. Hence, uncertainties from concentration calculation were reduced and therefore their impact on the kinetic parameters reduced.

## Quality control of experimental data through kinetic modeling

Besides data integrity, high quality data is required for accurate parameter estimates, whereas unnoticed experimental errors unintentionally distort the resulting kinetic parameters. Therefore, quality control of measurement data is an important aspect of every data analysis process. EnzymePynetics allows visualization of measurement data together with the fitted kinetic model. Ultimately, systematic deviations of single measurements were identified by simultaneously visualizing measurement data and the fitted model.
By this method, errors with pipetting, mixing, and temperature equilibration during experimental preparation as well as issues with the experimental design itself were identified in different scenarios.  
Systematic pipetting errors were disclosed by systematic deviations between the measurement data and the fitted model.
In scenario D individual enzyme reactions showed to have deviating substrate concentrations from the assay protocol, whereas modeling without the presumably erroneous measurements resulted in a good fit between data and model as well as low standard deviation on the parameter estimates. Hence the experimental partners were advised to repeat the deviating measurement.  
In two projects the progress-curve of the reactions indicate a lag phase after reaction start. In one case, small pipetting volumes of enzyme to start the reaction led to inhomogeneous mixing and thus increasing reaction rates until the enzyme was distributed evenly. As a measure, the project partners were advised to increase the pipetting volume of the enzyme solution.
In another case, temperature incubation effects in MTP caused by prolonged assay preparation times resulted in an initial lag phase. Due to small reaction volumes and low mass, MTPs have a low heat capacity and thus high susceptibility to temperature change. Hence, project partners were advised to pre-incubate the MTP within the photometer at reaction temperature for 5 min and then start the reaction by adding enzyme as quick as possible.
In another project, modeling results of the estimated $k_{cat}$ and $K_{m}$ were used to iteratively improve the design of the enzyme assay. In addition, the correlation between the aforementioned parameters was used as a measure to ensure that the highest substrate concentration applied was sufficiently high and thus the kinetic parameters could be determined independently.
Thereby, the appropriate enzyme concentration and substrate concentration range were identified through multiple round of lab experimental with subsequent kinetic modeling.

Assessing the data quality through modeling and enhancing assay designs proved to be a strength of the implemented progress curve method on which the parameter estimation of this workflow is based on. In contrast to the predominantly applied initial rates method {cite}`tang2010precise`, the progress curve method offers intrinsic advantages with regard to methods reproducibility.
In initial rate kinetics no consensus on the linear reaction period on which the kinetic parameters are estimated exists. Hence, the linear period is manually determined or only assumed. Ultimately, the resulting parameters are influenced by the choice of linear period. For progress curve analysis the entire dataset of an enzyme assay is used. If the fit statistics reveal that the data is not in accordance with the model, either the assay should be repeated due to an underlying issue or the model should be questioned. Contrarily, initial rates method can always be applied by arbitrarily choosing any linear subset or only the initial two data points of an time course data set. In discussions with project partners this showed to be a common practice although scientifically at least questionable, depending on the circumstances.

**importance of percise parameters and how my workflow helps**

fggf #TODO Importance of precise kinetic parameters for biocatalysis
With upcoming big data technologies like

Cumulative error, which ultimately leads to unreproducible kinetic parameters

- initial rates vs progress-curve analysis

## FAIRness

FAIR science

FAIR guiding principles are not exclusive to data but include methods as well as entire workflows and are a best practice in data steward ship{cite}`wilkinson2016fair`. Thus, FAIR guiding principles were implemented on multiple levels of this thesis. On the data level, all experimental data as well as the modeling results were stored in EnzymeML files, which are compliant with FAIR data principles {cite}`pleiss2021standardized` {cite}`range2022enzymeml`. On a methods level, the workflow and all of its components were designed in a FAIR fashion.  
The Python packages are findable and accessible on PyPI and GitHub, which present the most important distribution platforms for Python code. The software is interoperable with other software. This was achieved by segregating the data model from the functionality of the software. Hence other Python tools can utilize the functionalities of CaliPytion and EnzymePynetics by serving the underlying data models. Thus no modifications on the software functionalities are necessary. Furthermore, the data model and its vocabulary is described in the specifications of [CaliPytion](https://github.com/FAIRChemistry/CaliPytion/blob/main/specifications/CalibrationModel.md) as well as [EnzymePynetics](https://github.com/haeussma/EnzymePynetics/blob/main/specifications/EnzymeKinetics.md). Furthermore, the developed software packages are reusable, since the software can be installed and applied by anyone due to the documentation.  
Many of the FAIR software principles also apply on workflow level. Accessibility is provided by storing the Jupyter Notebooks on GitHub. Interoperability is given by the modular design, which allows to integrate the workflow in other Python-based workflows, which also makes it reusable.

einschränkung longtherm

Furthermore, progress curve analysis is more transparent about the underlying data compared

sharing the raw data and comprehensive description of individual data treatment steps,

- sophisticated methods, otherwise wrong conclusions (italian group)

**Jupyter Notebook and Book**

contemporary

Besides the data analysis, this entire thesis was conceptualized adhering to FAIR data principles . In consequence, this thesis was written was a Jupyter Book, which allow to combine multiple Jupyter Notebooks with text chapters in a structured document. Ultimately, this thesis is findable and accessible on GitHub, interoperable through Binder and the Jupyter Notebook format itself. Hence, making the work of this thesis reusable in a contemporary

Printed version is read twice, digital version unlimited hence optimized for digital experience

- Besides the methodology documentation is bad

- Need for documentation without limiting experimental possibilities -> Jupyter notebooks
- great documentation (FAIR)

## Reproducibility of the workflow for kinetic parameter estimation

Besides FAIR, science should also be fair.

Reproducibility is a relative term which needs reference.

In sum minimizing errors copying method, manual handling errors and (calculation)
**Methods reproducibility**

**Results reproducibility**

An secondary necessity for reproducibility is documentation. Thereby, both laboratory procedures as well as of data analytic procedures need to be documented sufficiently to make experimental

Methods reproducibility is given

- Lab vs computer documentation
- Strenda guidelines are insufficient?? focused on lab, neglecting data analytic aspect

## Results reproducibility of the parameter estimation workflow

## opportunities:

- how precise should parameters be
- randomized order in MTPs
