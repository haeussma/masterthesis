# Methods

## EnzymeML and PyEnzyme

EnzymeML documents were read, edited, and written via the Python application programming interface PyEnzyme {cite}`pyenzyme_2021`, providing the possibility to integrate data from EnzymeML documents in Python-based data analysis workflows. Thereby, EnzymeML documents were used as a data source.

### Creation of EnzymeML documents

`EnzymeMLDocuments` were created using the EnzymeML Excel template, in combination with the `.fromTemplate()` method of PyEnzyme. Within the spreadsheet, data and metadata of an experiment were filled in the respective sheets. Alternatively, only metadata was entered to the spreadsheet, whereas measurement data was parsed by a custom Python function from the output file of the analytical device to the measurement data section of the EnzymeML document.

### Saving of modeling results

Modeling results from parameter estimation were written back to the EnzymeML document, which contained the measurement data by using PyEnzyme.

## CaliPytion

CaliPytion was developed to provide an easy way to find the best model by which the relationship between analytical signal to analyte concentration of photometric measurements can be determined. The resulting `StandardCurve` is created based on a provided analyte standard. Besides a linear calibration equation, non-linear calibration models were implemented. All implemented calibration equations {eq}`linear`, {eq}`quadratic` {eq}`3rd_polynominal`, {eq}`poly_e` {eq}`fractional` are listed below, whereas $A$ denotes the absorption, $x$ the concentration, and $a, b, c$ the respective parameters of each equation.
During the fitting process, all calibration equations are fitted against the analyte standard by non-linear least-square minimization. Therefore, the Python library lmfit {cite}`newville2016lmfit` was utilized. The best fitting model is then selected based on the lowest Akaike information criterion (AIC). Furthermore the fit of the models can be visually assessed by calling the `.visualize()` method of the `StandardCurve` class. After the calibration curves are fitted, concentration of unknown samples are determined by calculating the root of the parameterized calibration model.

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

(CaliPytion_datamodel)=

### Data model

The `StandardCurve` functionality of CaliPythion is based on a [data model](https://github.com/FAIRChemistry/CaliPytion/blob/main/specifications/CalibrationModel.md), structuring data and metadata of calibration measurements. Thereby, calibration conditions like temperature, pH, and information about the analytical device can be provided. Furthermore, the name as well as an ID can be specified for the analyzed substance. All of this information is stored in the `Calibration` root object. Additionally, the root object can contain a `Standard` and a `Spectrum`.  
A `Standard` contains measurements of multiple predefined concentrations and the wavelength at which the measurement was conducted. Additionally, a `Spectrum` can be defined by providing measurement data as well as the respective wavelengths at which the data was measured. The data model was generated using [sdRDM](https://github.com/JR-1991/software-driven-rdm), a Python tool allowing the creation and versioning of data models.

### Initialization

A `StandardCurve` can be initialized directly from an Excel file by calling `StandardCurve.from_Excel()`
Therefore, the first column must contain concentration data, whereas all following columns contain the respective measurement data. The first row is not considered and should be used to label the data.
Alternatively a `StandardCurve` can be initialized by providing a `Calibration` object, containing the calibration data.

After `StandardCurve` initialization, concentrations can be calculated by calling the `.get_concentration()` method. Alternatively, a `StandardCurve` can be directly applied to an `EnzymeMLDocument` by calling the `.apply_to_EnzymeML()` method.

### Code availability

CaliPytion was published on the Python packaging index ([PyPI](https://pypi.org/project/CaliPytion/)) whereas the source code is accessible on [GitHub](https://github.com/FAIRChemistry/CaliPytion).

(method:enzymepynetics)=

## EnzymePynetics

EnzymePynetics is a Python package for kinetic parameter estimation of single-substrate enzyme reactions. The `ParameterEstimator` of EnzymePynetics estimates the turnover number and Michaelis-Menten constant by fitting time-course measurement data of enzyme reactions to different Michaelis-Menten models. Thereby, the residuals between measurement data and integrated Michaelis-Menten rate equations are minimized through a non-linear least-square algorithm utilizing the lmfit library. Additionally, the inhibition constant $K_{i}$ can be assessed for potential substrate or product inhibition. Furthermore, $K_{i}$ of an enzyme inhibitor apart from substrate or product can be determined for inhibitor studies.

(EnzymePynetics_datamodel)=

### Data model

EnzymePynetics is based on a data model, resembling the experimental design of an enzyme kinetics assay. Thereby, all relevant data and metadata of a kinetic experiment are structured in the `EnzymeKineticsExperiment` base object. Therein, reaction temperature with its respective unit, pH value, and the name of the substance from which the measurement data originates are described. An `EnzymeKineticsExperiment` contains one or multiple `Measurements`, which contain the measurement data. Each measurement stores information of the experimental conditions of an enzyme reaction. Hence, each `Measurement` contains information on the initial substrate concentration, enzyme concentration, and inhibitor concentration, along with the respective concentration units. Each `Measurement` contains measured data, which itself consist of one or multiple replicates of the respective experimental condition.  
The data model was generated using [sdRDM](https://github.com/JR-1991/software-driven-rdm).

Extensive documentation of the data model is accessible in the [specifications](https://github.com/haeussma/EnzymePynetics/blob/main/specifications/EnzymeKinetics.md) of the software package.

### Kinetic models

Besides the irreversible Michaelis-Menten rate equation {eq}`irreversible_mm` inhibition models for competitive {eq}`competitive_inhibition`, uncompetitive {eq}`uncompetitive_inhibition`, non-competitive {eq}`noncompetitive_inhibition` and partially competitive inhibition {eq}`partially_inhibition` were implemented, according to {cite}`bisswanger2017enzyme`. Thereby, $S$, $E$, and $I$ respectively denote the concentration of substrate, enzyme, and inhibitor. In terms of kinetic parameters, $k_{cat}$ denotes the turnover number, $K_{m}$ the Michaelis-Menten constant of the substrate, whereas $K_{ic}$ and $K_{iu}$ respectively describe the competitive and uncompetitive inhibition constant.

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
\frac{dS}{dt} = -\frac{k_{cat} * E * S}{K_{m} * (1+\frac{I}{K_{iu}}) + S}
```

```{math}
:label: noncompetitive_inhibition
\frac{dS}{dt} = -\frac{k_{cat} * E * S}{K_{m} * (1+\frac{I}{K_{ic}}) + (1+\frac{I}{K_{iu}}) + S}
```

```{math}
:label: partially_inhibition
\frac{dS}{dt} = -\frac{k_{cat} * E * S}{K_{m} * \frac{(1+\frac{I}{K_{ic}})}{(1+\frac{I}{K_{iu}})} + S}
```

By default, $E$ is assumed to be constant throughout the reaction. If required, the kinetic models can be extended by the parameter $k_{inact}$, which describes the time-dependent inactivation rate of the enzyme. The decrease in active enzyme concentration is modeled as describe in equation {eq}`enzyme_inactivation`:

```{math}
:label: enzyme_inactivation
\frac{dE}{dt} = -k_{inact} * E
```

Each kinetic model therefore consists of a set of ordinary differential equations, describing the change in concentration of substrate, product and enzyme.

### Run parameter estimation

**Initialization**  
The `EnzymeKineticsExperiment` object solely serves as a data container, whereas
the `ParameterEstimator` harbors the functionalities for parameter estimation. Data can be provided as an `EnzymeKineticsExperiment` object. Alternatively, an `EnzymeMLDocument` document can be provided as the data source via PyEnzyme by calling the `.from_EnzymeML()` method. Within the method, the ID of the measured species needs to be specified according to the nomenclature of the `EnzymeMLDocument`. If the document contains information on an inhibitor, the respective id needs to be provided as well, if its inhibition constant should be estimated.

Initially, missing data for the fitting process is calculated based on the assumption of mass conservation. As such, missing product concentration $P$ is calculated based on the given substrate measurement data $S$ and the specified initial substrate concentration $S_{0}$ for each measurement $t$ {eq}`missing_product`. If product data is provided, missing substrate concentrations are calculated accordingly {eq}`missing_substrate`.

```{math}
:label: missing_product
P_{t} = S_{0} - S_{t}
```

```{math}
:label: missing_substrate
S_{t} = S_{0} - P_{t}
```

**Fitting of models**  
After the `ParameterEstimator` is initialized, all kinetic models are initialized by calling the `.fit_models()` method. Furthermore, it can be specified if enzyme inactivation should be considered for modeling. Depending on whether an inhibitor was specified, the models are initialized accordingly. If no inhibitor was specified, product and substrate inhibition models are initialized besides the irreversible Michaelis model. If inhibitor data was provided, the inhibition models are initialized with concentration data of the inhibitor. Thereafter, kinetic parameters of each model are initialized with estimates based on the highest reaction rate in the dataset. For minimization, the parameter space is limited to ± 1000-fold of the initial parameter estimates.  
After all models are fitted, an overview table is printed which lists all kinetic parameters of each model together with the respective 1σ standard deviation in percent. The table is sorted by ascending AIC.

**Visualization**  
After model fitting, the modeling results can be visualized by calling the `.visualize()` method, which was implemented utilizing Matplotlib {cite}`caswell2020matplotlib`.
Thereby, the measurement data is visualized together with the fitted kinetic model. If the experiment was carried out in replicates, mean values of the measurement data with the corresponding standard deviations are shown. If the experiment data contained reactions with different inhibitor concentrations, individual inhibitor concentrations are denoted by individual markers.  
By default, the best fitting model according to its AIC is visualized. Different models can be visualized by passing the name of the respective model to the function call. Besides the figure, a detailed statistical fit report is printed, if not specified otherwise.

**Data units**  
So far, EnzymePynetics neither converts nor validates units of the measurement data. Hence, all data needs to be provided using the same molar or mass unit. Only the concentration unit of the inhibitor may differ from all other data units.

### Code availability

EnzymePynetics was published on the Python packaging index ([PyPI](https://pypi.org/project/EnzymePynetics/)) whereas the source code is available on [GitHub](https://github.com/haeussma/EnzymePynetics).

## Model comparison

The kinetic model, which best describes the experimental data, was selected based on AIC, standard deviation of parameter estimates, and visual fit between measurement data and fitted model. AIC serves as a statistical metric for information loss, which allows to relatively compare different models with differing number of parameters for a given data set. Hence, AIC can be applied for model selection {cite}`arnold2010uninformative`, {cite}`akaike1998information`. Thereby, models with a lower AIC indicate less information loss.
Since AIC is based on the chi-square statistic, it does not consider the standard deviation of the estimated parameters of a model. Therefore, standard deviation was additionally considered in model selection. Models with low standard deviation were therefore preferred over models with high standard deviation. Furthermore, fit quality was assessed visually by confirming that the model describes the progress curve of the experimental data.

## Jupyter Notebook

All developed tools were deployed in Jupyter Notebooks {cite}`kluyver2016jupyter` for data analysis within a Python environment. A Jupyter Notebook is a digital document containing code cells and text cells. Thereby, executable code with its corresponding output can be supplemented with narrative text, figures, and equations.
Besides its documentation capabilities, notebooks can be uploaded to platforms like GitHub or be deployed with Binder {cite}`ragan2018binder`.
GitHub is a cloud-based Git repository, which enables versioning as well as collaborative working on projects involving any kind of code. Besides GitHub, sharing of Jupiter Notebooks is possible with the open-source project Binder, which allows hosting of notebooks. Thus, Jupyter Notebooks can be shared in an executable form, which is independent from the operating system and the location of the user. Only an internet connection and a browser is required.
As a result, Jupyter Notebooks, which are executable through Binder, form an user-friendly environment for scientific data analysis. Due to the extensive documentation capabilities of Jupyter Notebooks, the analysis workflow is comprehensible and accessible for programmers as well as for non-programmers.

## Jupyter Book

This written thesis, except for the cover page, was entirely conceptualized as a Jupyter Book, which is an open-source tool that allows to incorporate computational material in publication-quality books {cite}`jupyterbook_2020`. Thereby, Jupyter Notebooks can natively be integrated in Jupyter Books, thus making all advantages of the notebook also available in the book.
Jupyter Books allow text formatting with Markdown, which is a simple markdown language. Furthermore, Jupyter Book supports cross-referencing, citations, and numbering of equations.
Together with all data, the final and interactive book was deployed on GitHub Pages whereas the print version was generated by exporting the JupyterBook to LaTeX. The LaTeX version was then adjusted for printing.

All data of this JupyterBook as well as the experimental data is stored on the data repository of the University of Stuttgart (DaRUS) and can be found under the following DOI: [https://doi.org/10.18419/darus-3337](https://doi.org/10.18419/darus-3337).
In order to host the JupyterBook on the internet, the data from DaRUS is mirrored to GitHub. The JupyterBook can be found under the following URL: [https://haeussma.github.io/masterthesis/welcome.html](https://haeussma.github.io/masterthesis/welcome.html).
