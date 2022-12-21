# Conclusion and outlook

In this project, a reproducible and FAIR workflow from analytical raw-data to kinetic parameters was established. Thereby, the output of a analytical device is directly used as the input of the modeling pipeline. By considering possible nonlinear relationships in the conversion of analytic data to concentration data, a source of inaccuracies is avoided. The parameter estimation is based on progress curve analysis. Since various kinetic models were implemented, models can be compared and insight into the kinetic mechanism can be gained. Therefore, the implemented parameter estimator is suited for inhibition studies. After fitting of the models, visualization of the fitted model together with measurement data allows to identify individual systematic deviation between the model and the data. Hence, providing a measure for quality control.
The workflow is based on the EnzymeML format

The developed workflow allows scientist with basic programming knowledge to utilize progress curve analysis of
In combination with the EnzymeML format the workflow allows
By using
The reproducibility of the workflow is mainly attributed to the implemented progress curve analysis method, together with the intrinsic characteristics of the Jupyter Notebook format. Due to the format and the associated Jupyter ecosystem, the analysis can be easily shared.
Thereby the

Cumulative error, which ultimately leads to unreproducible kinetic parameters

In sum minimizing errors copying method, manual handling errors and (calculation)
Analysis of the four scenarios proved the applicability and advantages of the developed workflow.

In combination with jupyternotebook the developed workflow is fair. The EnzymeMLdoument istselfe is not

In future developments, the workflow could be extended to enable parameter estimation of enzyme reactions with more than one substrate, by implementing more complex models. Furthermore, part of the quality control, which is currently carried out by visualizing, could be automated. Thereby, the user would be informed about systematic errors in the data set, which might have resulted from experimental errors. Moreover, tools could be implemented which assist scientist in assay design with focus on the design of minimal experiments. Hence, after an initial analysis round the user could be informed if the given information is sufficient for parameter estimation, or if the data set contains to little information. In result, initial substrate concentrations for the next round of lab experiments could be suggested.
