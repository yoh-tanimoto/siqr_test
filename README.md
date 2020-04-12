# siqr_test
The SIQR model of epidemiology with test and confounding disease.

This is a python code of a modified version of the SIQR model of epidemiology. A fixed number of tests is conducted every day and if one is infected and tested, one goes to Q(uarantined). Tests are conducted on people with symptoms, and there is another confounding disease with the same symptoms with a fixed number of patients.

# trasym
This is a variation of the SIQR model, where there is a compartment "A" of asymptomatic carriers. Apart from the usual SIQR model, a part of infected is asymptomatic and hence is not quarantined, but goes to the compartment A. Those in A can infect further and the secondary infection can be either symptomatic or asymptomatic. There is also a maximum number of tests, so that only a certain number of infected can be quarantined a day.

The original code of the SIR model by Christian Hill is here
https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/

