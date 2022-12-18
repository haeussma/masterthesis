# Discussion

reproducibility:

- Data integrity
- Modeling documentation

Hence, the analysis workflow was conceptualized starting from photometric measurement data. For that purpose, discussing the output format of the analytical device with the experimental partners proved to be valuable to plan the data acquisition process. Hence, a project specific parser function was established, which transferred the measurement data to the respective EnzymeML containing information on the measurement conditions. Thereby, error-prone and tedious manual copying between files was avoided, ensuring raw data integrity.

## Data integrity through accurate concentration calculation

Another aspect of data integrity in the developed workflow is the accurate concentration calculation based on the measured absorption signal. The relation between concentration to absorption signal of an analyte can be linear as well as non-liner. Therefore, wrong assumption of a linear relationship albeit the relation is better described non-linearly, leads to avoidable inaccuracy in concentration calculation {cite}`hsu2010effect` {cite}`martin2017fitting` and thus imprecise kinetic parameter estimates.
Hence, concentration calculation based on extinction coefficients was discouraged, if the underlying relation between absorption and concentration of an analyte was not self-determined for the respective reaction system.
To enable precise concentration calculation of linear and non-linear absorption to concentration relationships, CaliPytion was applied with analyte standards provided by the experimental partners.
All of the analyzed standard curves, except for scenario C ABTS at pH 5, were best described by non-linear calibration equations based on their fit quality determined by AIC. Hence, uncertainties from concentration calculation were reduced and therefore their impact on the kinetic parameters reduced.

## Quality control of experimental data through kinetic modeling

Besides data integrity, high quality data is required for accurate parameter estimates, whereas unnoticed experimental errors unintentionally distort the resulting kinetic parameters. Therefore, quality control of measurement data is an important aspect of every data analysis process. EnzymePynetics allows visualization of measurement data together with the fitted kinetic model. Ultimately, systematic deviations of single measurements were identified by simultaneously visualizing measurement data and the fitted model.
By this method, errors with pipetting, mixing, and temperature equilibration during experimental preparation as well as issues with the experimental design itself were identified in different scenarios.  
Systematic pipetting errors were disclosed by systematic deviations between the measurement data and the fitted model.
In scenario D individual enzyme reactions showed to have deviating substrate concentrations from the assay protocol, whereas modeling without the presumably erroneous measurements resulted in a good fit between data and model as well as low standard deviation on the parameter estimates. Hence the experimental partners were advised to repeat the deviating measurement.  
In two projects the progress-curve of the reactions indicate a lag phase after reaction start. In one case, small pipetting volumes of enzyme to start the reaction led to inhomogeneous mixing and thus increasing reaction rates until the enzyme was distributed evenly. In result, the project partners were advised to increase the pipetting volume of the enzyme solution.
In another case, temperature incubation effects in MTP caused by prolonged assay preparation times resulted in an initial lag phase. Due to small reaction volumes and low mass, MTPs have a low heat capacity and thus high susceptibility to temperature change. Hence, project partners were advised to pre-incubate the MTP within the photometer at reaction temperature for 5 min and then start the reaction by adding enzyme as quick as possible.
In another project, modeling predictions of the estimated $k_{cat}$ and $K_{m}$ together with their correlation were used to iteratively improve the design of the enzyme assay. Thereby, the appropriate enzyme concentration and substrate concentration range were identified through multiple round of lab experimental with subsequent kinetic modeling.

Assessing the data quality through modeling proved to be a strength of fitting integrated rate equations to continuous assay data from photometric measurements. Besides quality control, the progress curve analysis method offer several advantages over the predominantly applied initial rates method {cite}`tang2010precise`, if more than two measurement points exist. A key advantage is reproducibility of the analysis. In initial rate kinetics no consensus on the linear reaction period exists. Hence, the linear period on which the kinetic parameters are estimated on is manually determined or assumed. Ultimately, the resulting parameters are influenced by the choice of linear period. By only evaluating the linear reaction period, only a subset of all available data is utilized for parameter estimation.

This emphasizes the advantages of progress curve analysis over the predominately used initial rates method {cite}`tang2010precise`. Although results from initial rate kinetics can as well be visualized as a double-reciprocal (Lineweaver–Burk) diagram

This is because of both the method and the data which is required for the method.

Kinetic parameters are mostly determined by initial rates method, whereas initial rates of reactions with varying substrate concentration are determined and either . Thereby, rates from a subjectively chosen linear time period of the reaction are chosen

Cumulative error, which ultimately leads to unreproducible kinetic parameters

- Lab vs computer documentation

- initial rates vs progress-curve analysis

## What makes this workflow experimental data for enzyme kinetics

Availability of raw data is the prerequisite to reproduce the results of an data analysis process.

An secondary necessity for reproducibility is documentation. Thereby, both laboratory procedures as well as of data analytic procedures need to be documented sufficiently to make experimental

- Currently: unfair initial rates, lineweaver burk
  sharing the raw data and comprehensive description of individual data treatment steps,

**Findable**

**Accessible**

**Interoperable**

**Reusable**

- Strenda guidelines are insufficient??
- Precision of parameters

**model selection and miss selection**

- kinetic and calibration
- correlation between parameters:
- short experimental time
- low initial substrate concenration

- how precise should parameters be
- randomized order in MTPs
- ENZYME INACTIVATION

  - initial rates not suited to estimate conversion after 24 h.

- sophisticated methods, otherwise wrong conclusions (italian group)

**Jupyter Notebook and Book**

- Besides the methodology documentation is bad

- Need for documentation without limiting experimental possibilities -> Jupyter notebooks
- great documentation (FAIR)

In this project, a robust and FAIR workflow from analytical raw-data to kinetic parameters was established. Thereby, the output of a analytical device was used directly as the input of the modeling pipeline.
