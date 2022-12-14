# Discussion

Analysis of the four scenarios proved the applicability and advantages of the developed workflow.

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

- General problems in datasets and experimental preparation
  - systematic errors, temperature (gradients and incubation)
  - outliers (pipetting errors,)
  - experimental problems (viscosity)
- generals problems in modelling
  - error sources
    - wrong concentration calculation calibration equation and extincition coefficient
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
