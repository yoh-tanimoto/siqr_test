# siqr_test
The SIQR model of epidemiology with test and confounding disease.

This is a python code of a modified version of the SIQR model of epidemiology. A fixed number of tests is conducted every day and if one is infected and tested, one goes to Q(uarantined). Tests are conducted on people with symptoms, and there is another confounding disease with the same symptoms with a fixed number of patients.

# siqar
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

# siqar_test

Here the model is written as a difference equation (not a differential equation). Contact-tracing is modelled in such a way that, if some symptomatic is tested positive (and hence quarantined), we assume that there is one asymptomatic in A and this person gets isolated as well (hence giving a partial answer to the problem above). A better way of modelling contact-tracing is welcome.

# seir_ld
A variation of the SEIR model with a lockdown, or a similar sudden change of contact rate.

# trasym
This is a variation of the SIQR model, where there is a compartment "A" of asymptomatic carriers. Apart from the usual SIQR model, a part of infected is asymptomatic and hence is not quarantined, but goes to the compartment A. Those in A can infect further and the secondary infection can be either symptomatic or asymptomatic. There is also a maximum number of tests, so that only a certain number of infected can be quarantined a day.

# tracing
A variation of trasym, where asymptomatic carriers are also detected and quarantined with a certain rate.

# massteststratified

A variation of the SIR model with testing. In this model, the population is split into two groups (N1 and N2), those who get tested and those who do not.
The testing rate is r (r * N tests every day), carried out only to S1 and I1. We assume that the contacts between the groups N1 and N2 are such that if one person in the group N2
have M contacts per day, then mu * M are in N1 and (1-mu) * M are in N2. The differential equation is then given by

    dS1dt = -beta * (S1 * I1 * (1 - mu * N2/N1) / N1 + S1 * I2 * mu / N1)

    dS2dt = -beta * (S2 * I1 *  mu * N2/N1 / N2 + S2 * I2 * (1-mu) / N2)

    dI1dt = beta * (S1 * I1 * (1 - mu * N2/N1) / N1 + S1 * I2 * mu / N1)  - gamma * I1 - s * r * N * I1 / (S1+I1)

    dI2dt = beta * (S2 * I1 *  mu * N2/N1 / N2 + S2 * I2 * (1-mu) / N2)  - gamma * I2

    dR1dt = gamma * I1  + s * r * N * I1 / (S1+I1)

    dR2dt = gamma * I2

The linearized equation about I1, I2 is given by the matrix ((beta * (1 - mu * N2/N1) - gamma - s * r * N/N_1, beta * mu), (beta * mu * N2/N1, beta * (1-mu) - gamma)).
Whether the infection grows exponentially or not is determined by the determinant of this matrix (it grows if the determinant is negative).

The original code of the SIR model by Christian Hill is here
https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/
