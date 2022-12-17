# Discussion

reproducibility:

- Data integrity
- Modeling documentation

Availability of raw data is the prerequisite for reproducible data analysis. The lack of raw data is therefore one of the reasons for the reproducibility crisis {cite:cts}`miyakawa2020no`. Hence, the analysis workflow was conceptualized starting from photometric measurement data. For that purpose, discussing the output format of the analytical device with the experimental partners proved to be valuable to plan the data acquisition process. Hence, a project specific parser function was established, which transferred the measurement data to the respective EnzymeML containing information on the measurement conditions. Thereby, error-prone and tedious manual copying between files was avoided, ensuring raw data integrity.

## Data integrity through accurate concentration calculation

Another aspect of data integrity in the developed workflow is the accurate concentration calculation based on the measured absorption signal. The relation between concentration to absorption signal of an analyte can be linear as well as non-liner. Therefore, wrong assumption of a linear relationship albeit the relation is better described non-linearly, leads to avoidable inaccuracy in concentration calculation {cite}`hsu2010effect` {cite}`martin2017fitting` and thus imprecise kinetic parameter estimates.
Hence, concentration calculation based on extinction coefficients was discouraged, if the underlying relation between absorption and concentration of an analyte was not self-determined for the respective reaction system.
To enable precise concentration calculation of linear and non-linear absorption to concentration relationships, CaliPytion was applied with analyte standards provided by the experimental partners.
All of the analyzed standard curves, except for scenario C ABTS at pH 5, were best described by non-linear calibration equations based on their fit quality determined by AIC. Hence, uncertainties from concentration calculation were reduced and therefore their impact on the kinetic parameters reduced.

## Quality control of experimental data through kinetic modeling

Besides data integrity, high quality data is required for accurate parameter estimates, whereas unnoticed experimental errors unintentionally distort the resulting kinetic parameters. Therefore, quality control of the measurement data is an important aspect of every data analysis process. EnzymePynetics allows visualization of measurement data together with the fitted kinetic model. Ultimately, systematic deviations of single measurements were identified by simultaneously visualizing measurement data and the fitted model.
By this method, systematic errors with pipetting, mixing, and temperature equilibration during experimental preparation as well as issues with the experimental design itself were identified in different scenarios.  
Systematic pipetting errors were disclosed by systematic deviations between the measurement data and the fitted model.
In scenario D individual enzyme reactions showed to have deviating substrate concentrations from the assay protocol, whereas modeling without the presumably erroneous measurements resulted in a good fit between data and model as well as low standard deviation on the parameter estimates. Hence the experimental partners were advised to repeat the deviating measurement.  
In two projects the progress-curve of the reactions indicate a lag phase after reaction start. In one case, small pipetting volumes of enzyme to start the reaction led to inhomogeneous mixing and thus increasing reaction rates until the enzyme was distributed evenly. In result, the project partners were advised to increase the pipetting volume of the enzyme solution.
In another case, temperature incubation effects in MTP caused by prolonged assay preparation times resulted in an initial lag phase. Due to small reaction volumes and low mass, MTPs have a low heat capacity and thus high susceptibility to temperature change. Hence, project partners were advised to pre-incubate the MTP within the photometer at reaction temperature for 5 min and then start the reaction by adding enzyme as quick as possible.
In another project, modeling predictions of the estimated $k_{cat}$ and $K_{m}$ together with their correlation were used to iteratively improve the design of the enzyme assay. Thereby, appropriate enzyme concentration and substrate concentration range were identified through multiple round of lab experimental with subsequent kinetic modeling.

All of the mentioned experimental issues were disclosed by visualizing measurement data together with fitted rate equations. This emphasizes the advantages of progress curve analysis, which was implemented in EnzymePynetics over the commonly used initial rates method.

This is because of both the method and the data which is required for the method.

Kinetic parameters are mostly determined by initial rates method, whereas initial rates of reactions with varying substrate concentration are determined and either . Thereby, rates from a subjectively chosen linear time period of the reaction are chosen

Cumulative error, which ultimately leads to unreproducible kinetic parameters

- Lab vs computer documentation

sharing the raw data and comprehensive description of individual data treatment steps,

Raw data liked to executable code is the gold standard in computational science {cite}`peng2011reproducible`

- Currently: unfair initial rates, lineweaver burk
- Besides the methodology documentation is bad

## FAIR workflow for enzyme kinetics

**Findable**

**Accessible**

**Interoperable**

**Reusable**

### Reproducible enzyme kinetics

- Strenda guidelines are insufficient??
- Precision of parameters

**errors in data sets and detection of them**

- generals problems in modelling
  - error sources
    - wrong concentration
    - initial rates

**model selection and miss selection**

- kinetic and calibration
- correlation between parameters:
- short experimental time
- low initial substrate concenration
- unpercise data --> systematic deviations between mesaurements

- how precise should parameters be
- initial rates vs progress-curve analysis
- randomized order in MTPs
- ENZYME INACTIVATION
  - initial rates not suited to estimate conversion after 24 h.

**Jupyter Notebook and Book**

- Need for documentation without limiting experimental possibilities -> Jupyter notebooks
- great documentation (FAIR)

In this project, a robust and FAIR workflow from analytical raw-data to kinetic parameters was established. Thereby, the output of a analytical device was used directly as the input of the modeling pipeline.

Despite the robust workflow, after the data is written to an EnzymeML document, manual copying of data into the EnzymeML spreadsheet is prone to human error. Especially, for large datasets with many time points and reactants. Hence, manufacturers should enable direct access to the measurement device, for instance via an API. (Control and laer experiments)
