# Methods

## 1. EnzymeML

In its core, EnzymeML is a data model, structuring data and metadata of biocatalytic reactions.

An EnzymeML document stores information on the reaction conditions like pH, temperature,
Thereby relevant information on the reaction condition like pH
Therein information on
substrates, and products of enzyme reactions and
THereby information on molecules
Thereby, properties of the reaction vessel, molecules involved in the reaction, as well as reaction properties like stoichiometry, pH, and temperature can be documented. Furthermore, properties of the protein catalyst like sequence, EC-number, UniProt ID can be stored alongside of potential reaction modifiers like inhibitory molecules. {cite}`pleiss2021standardized`

Within the measurements section of the data model, measurement data of one or multiple measurements is stored. Each measurement contains data of identical measurement conditions as well as information about these conditions. Thus, data and meta data from experiments with differing conditions can be stored in a structured way.

The data model consists of descriptions for the reaction vessel,

(method:calipytion)=

## CaliPytion

### Data model for calibration data

### Linear and non-linear fitting of calibration equations

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

### 3.5 Visualization

## 4. Model selection

Akaike information criterion
