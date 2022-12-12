# Methods

## 1. EnzymeML

In its core, EnzymeML is a data model, structuring data and metadata of biocatalytic reactions. Thereby, information on reaction conditions, substrate and product measurement data, as well as estimated kinetic parameters are documented {cite}`pleiss2021standardized`. Additionally, the enzyme and the reactants are specified by protein sequence or InChI respectively, enabling clear description of all involved species in an enzyme reaction. EnzymeML is conceptualized as an exchange format based on the Systems Biology Markup Language {cite}`hucka2003systems` and adheres to STRENDA guidelines {cite}`tipton2014standards`, which define minimal experimental reporting standards. Hence, EnzymeML serves as an exchange format for biocatalytic data between experimentalist, modelers and database providers, which is compliant with FAIR data principles {cite}`wilkinson2016fair`.  
EnzymeML documents can be read, edited, and written via the Python API PyEnzyme, providing the possibility to integrate PyEnzyme in Python data analysis workflows.

(method:calipytion)=

## CaliPytion

CaliPytion was developed to provide an easy way to use linear and non linear calibration equations to calculate concentrations of analytical raw data based on a standard measurements.

### Data model for calibration data

The standard curve functionality of CaliPythion is based on a data model, structuring data and metadata of calibration measurements. Thereby, calibration conditions like temperature and pH, and information of the analytical device can be provided. Furthermore, the name as well as an ID can be specified for the analyzed substance. All of this information is stored in the `Calibration` root object. Additionally, the root object can contain a `Standard` and a `Spectrum`. A standard contains measurements of multiple predefined concentrations and a wavelength at which the measurement was conducted, in the case of spectrophotometry. A `Spectrum` can be defined by providing measurement data as well as the respectuve wavelengths at which the data was measured.

### Linear and non-linear fitting of calibration equations

Finds optimal model parameters through non-linear least squares fitting of experimental data to models
Since no true model from the relation between ..., and calibrations ar known to not always be linear

the following equations were implemented

```{math}
:label: linear
A = ax
```

```{math}
:label: quadratic
A = ax^2 + bx
```

```{math}
:label: 3rd_polynominal
A = ax^3 + bx^2 + cx
```

```{math}
:label: poly_e
A = ae^{\frac{b}{x}}
```

```{math}
:label: fractional
A = \frac{ax}{b+x}
```

When a `StandardCurve` is created, the measurement data, which is defined in `Calibration` is used to fit the listed calibration model equations to the data. The model, which represents the relation between analytical signal and concentration the best, is determined based on the lowest Akaike information criterion (AIC). Concentrations are calculated by calculating the root of the fitted calibration model.  
After `StandardCurve` initialization, concentrations can be calculated by calling the `get_concentration()` method. Alternatively, a `StandardCurve` can be directly applied to an `EnzymeMLDocument` by calling the `apply_to_EnzymeML()` method

```python
product_standard = StandardCurve()
enzmldoc_absorption = EnzymeMLDocument()

enzmldoc_concentration = product_standard.apply_to_EnzymeML(
    enzmldoc=enzmldoc_absorption,
    species_id="s1")
```

(method:enzymepynetics)=

## EnzymePynetics

### Overview

EnzymePynetics is a a python package, for kinetic parameter estimation of single-substrate enzyme reactions, which was developed during this thesis. The `ParameterEstimator` of EnzymePynetics estimates the kinetic parameters $k_{cat}$ and $K_{m}$ by fitting time-course measurement data of enzyme reactions to different Michaelis-Menten models. Thereby, the residuals between measurement data and integrated Michaelis-Menten rate equations are minimized through a non-linear least-squares algorithm. Additionally, the inhibition constant $K_{i}$ can be assessed for potential substrate or product inhibition. Furthermore, the inhibition constant of an enzyme inhibitor can be determined.

### Data model

Data models build the backbone of applications, by defining the relations between informations in an hierarchical manner.
EnzymePynetics is based on a data model, resembling the experimental design of an enzyme kinetics assay. Thereby, all relevant data and meta data of an kinetic experiment are ordered in the base object `EnzymeKineticsExperiment`. On the metadata side, the base object consists of the attributes temperature with its respective unit, pH, and the name of the measured substance. Additionally, it can be specified whether the measurement data originates from substrate or product measurements. On the data side, `EnzymeKineticsExperiment` contains one or multiple `Measurements`. Each measurement stores the information of an experimental condition, to which the enzyme was subjected. Therefore, each `Measurement` contains information on the initial substrate concentration, enzyme concentration, and inhibitor concentration, if present, along with the respective concentration units. Each `Measurement` contains the measured data, which itself consist of one or multiple replicates of the respective experimental condition.  
The data model was was generated using [sdRDM](https://github.com/JR-1991/software-driven-rdm), a python tool allowing the creation and versioning of data models.

An extensive documentation of the data model can be accessed in the [specifications](https://github.com/haeussma/EnzymePynetics/blob/main/specifications/EnzymeKinetics.md) of the the software package.

### 3.3 Kinetic models

Besides the irreversible Michaelis-Menten rate equation (Eq. 1) inhibition models for competitive (Eq. 2), uncompetitive (Eq. 3), and non-competitive) inhibition (Eq. 4) were implemented. Thereby, $S$, $E$, and $I$ denote the concentration of substrate, enzyme, and inhibitor, respectively. In therms of kinetic parameters, $k_{cat}$ denotes the turnover number, $K_{m}$ the Michaelis-Menten constant of the substrate, whereas $K_{ic}$ and $K_{iu}$ describe the competitive and uncompetitive inhibition constant, respectively.

```{math}
:label: irreversible_mm
\frac{dS}{dt} = -\frac{k_{cat} * E * S}{K_{m} + S}
```

```{math}
:label: competitive_inhibition
\frac{dS}{dt} = -\frac{k_{cat} * E * S}{K_{m} * (1+\frac{I}{K_{ic}}) + S}
```

```{math}
:label: uncompetitive_inhibition
\frac{dS}{dt} = -\frac{k_{cat} * E * S}{K_{m} * (1+\frac{I}{K_{iu}}) * S}
```

```{math}
:label: noncompetitive_inhibition
\frac{dS}{dt} = -\frac{k_{cat} * E * S}{K_{m} * (1+\frac{I}{K_{ic}}) + (1+\frac{I}{K_{iu}}) + S}
```

By default, $E$ is assumed to be constant throughout the reaction. If required, the kinetic models can be extended by the parameter $K_{inact}$, which describes the time-dependent inactivation rate of the enzyme. The decrease in active enzyme is modeled by Eq. 5:

```{math}
:label: enzyme_inactivation
\frac{dE}{dt} = -K_{inact} * E
```

### 3.4 Initialization

Whereas the `EnzymeKineticsExperiment` object solely serves a
The `ParameterEstimator` harbors the functionalities for parameter estimation. Data can be provided by passing data as an `EnzymeKineticsExperiment` object. Alternatively, an EnzymeML documents can be provided as the data source via [PyEnzyme](https://github.com/EnzymeML/PyEnzyme) software.  
Initially, product concentration is calculated, if the input data is from substrate measurements. Substrate is calculated, if product data was provided. The calculation is based on the initial substrate concentration, which needs to be provided for all measurements. The parameter estimation is then based on substrate data.  
In the background, rough estimates for $k_{cat}$, $K_{m}$, and $K_{i}$ are calculated based on the fastest reaction rate in the provided dataset. Initial parameter estimates are needed as a starting point for the solver, whereas the parameter space is limited to ±1000-fold its estimated value for each parameter.

```python
from EnzymePynetics.tools.parameterestimator import ParameterEstimator
from EnzymePynetics.core.enzymekineticsexperiment import EnzymeKineticsExperiment
import pyenzyme

# Define data
experimental_data = EnzymeKineticsExperiment(
    pH=7,
    reactant_name="test substance",
    ...)

enzymeml_document = pyenzyme.EnzymeMLDocument.fromFile("enzymeML_dataset.omex")

# Initialize the parameter estimator from an 'EnzymeKineticsExperiment' instance
estimator = ParameterEstimator(data=experimental_data)

# ... or directly from an EnzymeML document.
estimator = ParameterEstimator.from_EnzymeML(
    enzmldoc=enzymeml_document,
    reactant_id="s1",
    measured_species="product",
    inhibitor_id="s2")

# Fit experimental data to kinetic models and get the fit report of all models
estimator.fit_models()

# Visualize data with the best fitting model
estimator.visualize()
```

### 3.5 Model selection

### 3.6 Visualization

Akaike information criterion
