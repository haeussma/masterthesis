# Discussion

Availability of raw data is the prerequisite for reproducible data analysis {cite}`miyakawa2020no`. Hence, the analysis workflow was conceptualized starting from photometric measurement data. For that purpose, discussing the output format of the analytical device with the experimental partners proved to be valuable to plan the data acquisition process. Hence, a project specific parser function was established, which transferred the measurement data to the respective EnzymeML document with information on the measurement conditions. Thereby, error-prone and tedious manual copying between files was avoided, ensuring raw data integrity.  
Since modeling was based on concentration data, thoroughly raw data preparation was crucial not to influence the modeling results by wrong concentration calculations. Besides blanking of the absorption signal, correct concentration calculation is important for accurate parameter estimates. Therefore, in three of the four scenarios standard curves of the analyte under assay conditions were recorded, whereas in scenario D concentrations were calculated via an extinction coefficient. All of the analyzed standard curves, except for scenario C ABTS at pH 5, were best described by non-linear calibration equations based on their fit quality determined by AIC.
Since, experiments were design to take advantage of the full measurement range of the photometer, which includes the non-linear detection range below the upper detection limit, non-linear relations between analyte and instrument response were expected.
This was expected, since

Concentration calculation by extinction coefficient is discouraged, if it is not specifically determined
as well as concentration calculation

Transferring the measurement signal into concentration values constitutes the second step of the workflow.

Cumulative error, which ultimately leads to unreproducible kinetic parameters

initial rates by linear regressin are highly sensitive to the chosen observaation window

- Lab vs computer documentation

sharing the raw data and comprehensive description of individual data treatment steps,
In biocatalysis

Raw data liked to executable code is the gold standard in computational science {cite}`peng2011reproducible`

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
