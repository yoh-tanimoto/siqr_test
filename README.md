# siqr_test
The SIQR model of epidemiology with test and confounding disease.

This is a python code of a modified version of the SIQR model of epidemiology. A fixed number of tests is conducted every day and if one is infected and tested, one goes to Q(uarantined). Tests are conducted on people with symptoms, and there is another confounding disease with the same symptoms with a fixed number of patients.

# trasym
This is a variation of the SIQR model, where there is a compartment "A" of asymptomatic carriers. Apart from the usual SIQR model, a part of infected is asymptomatic and hence is not quarantined, but goes to the compartment A. Those in A can infect further and the secondary infection can be either symptomatic or asymptomatic. There is also a maximum number of tests, so that only a certain number of infected can be quarantined a day.

# tracing
A variation of trasym, where asymptomatic carriers are also detected and quarantined with a certain rate.

# import
A variation of siqr-test, where there is the compoartment of asymptomatic carriers (A). There is a fixed number of confounding disease (Flu) and a fixed number of tests are conducted every day, both to A and to I. If a carrier (symptomatic or not) is tested, s/he goes to Q. R is the number of recovered from I and A and Rq is the recovered from Q (this distinction is made to separate the total number of ever-quarantined (Q + Rq), or total positive).
Interestingly, whether the contagion diverges or not depends both on the number of tests and the initial number of infection.

The differential equation is the following:

    dSdt = -(beta1+beta2) * S * (I+A) / N
    
    dIdt = beta1 * S * (I+A) / N - Test * I / (2 * (I + Flu)) - gamma1 * I
    
    dQdt = Test * I / (2 * (I + Flu)) + Test * I / (2 * (A + Flu)) - gamma1 * Q
    
    dAdt = beta2 * S * (I+A) / N - Test * I / (2 * (A + Flu)) - gamma2 * A
    
    dRdt = gamma1 * I + gamma2 * A
    
    dRqdt = gamma1 * Q
    
where S: susceptible, I: infected (with symptons), Q: quarantined, A: asymptomatic, R: recovered.
The sum I + A contributes to new infection (dIdt and dAdt), which is proportional to the number of contacts with S. As for I (symptomatic), the half of the tests are conducted on each day on the population I + Flu, so (assuming the 100% sensitivity of the test) only Test * I / (2 * (I + Flu)) are quarantined. The testing on A should be modelled on contact-tracing, but currently it is not implemented, so the equation for A is the same as that for I, except parameters.
The plots are made for I+A (real infected, not just observed), Q+Rq (total positive), new positive and the rate of positive among tests.

# seir_ld
A variation of the SEIR model with a lockdown, or a similar sudden change of contact rate.


The original code of the SIR model by Christian Hill is here
https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/
