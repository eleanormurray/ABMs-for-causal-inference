# ABMs-for-causal-inference
The respository 'ABMs for causal inference' archives the SAS code for comparing agent-based models and the parametric g-formula, published in the American Journal of Epidemiology in 2017 [1]. The code is also included as an appendix to the American Journal of Epidemiology publication. 

This repository includes code for simulating data with time-varying treatment-confounder feedback, and estimating the effect of treatment using an agent-based model and the parametric g-formula. This code demonstrates the potential for biases when using agent-based models for causal inference. SAS 9.4 was used for all analyses. 

The code appendix contains the following programs:
1.	RunMacros.sas
This program takes as input the parameters used to simulate the dataset for analyses, the number of subjects in the desired simulated cohort, and the number of bootstrap samples to run, and feeds this information into the macros which perform the analyses.
2.	Macros.sas 
This program contains the macros required to execute the simulation, except the previously published g-formula macro[2] (version: June-2015), restricted cubic spline macro[3,4], and SAS %lst macro[5].

References:

[1] Murray EJ, Robins JM, George R. Seage III, Freedberg KA, Hernan MA. A Comparison of Agent-Based Models and the Parametric G-Formula for Causal Inference. American Journal of Epidemiology. 2017;186(2):131-42 [http://bit.ly/2tKYt75].
[2]	Taubman SL, Young JG, Picciotto S, Logan R, Hernán MA. G-formula SAS macro 2012. Version June-2015. Available from: http://www.hsph.harvard.edu/causal/software/.
[3] Harrell FE. %rcspline macro. Clinical Biostatistics Duke University Medical Center, 1988.
[4] Devlin TF, Weeks BJ, editors. Spline functions for logistic regression modeling. Proceedings of the Eleventh Annual SAS Users Group International; 1986 February 9-12; Atlanta, Georgia: Cary NC: SAS Institute.
[5] Sample 44124: Counting the number of missing and non-missing values for each variable in a data set: SAS Institute Inc;  [11/30/2016]. Available from: http://support.sas.com/kb/44/124.html.


