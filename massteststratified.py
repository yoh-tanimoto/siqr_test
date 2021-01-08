# the original python code of the SIR model is https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/ by Christian Hill
# This is a variation of the SIR model with mass testing and isolation, with two groups of people who
# get/do not get testing (S1 and I1 get testing, while S2 and I2 do not get testing)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population: N
N = 100000000
# Testing rate, sensitivity, acceptance
r = 0.3
s = 1
a = 0.8

# Initial number of infected, quarantined, and removed individuals, I0 and R0.
I10, I20, R10, R20 = 100, 100, 0, 0
# Everyone else, S0, is susceptible to infection initially. Testing acceptance, a.
S10 = a * (N - I10 - I20 - R10 - R20)
S20 = (1-a)*(N - I10 - I20 - R10- R20)
N1 = a * N
N2 = (1-a) * N

# Contact rate, beta, and mean recovery rate, gamma, (in 1/days), mixing rate, mu.
beta, gamma, mu = 0.3, 0.15, 0.1
# A grid of time points (in days)
t = np.linspace(0, 500, 500)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma, mu):
    S1, S2, I1, I2, R1, R2 = y
    dS1dt = -beta * (S1 * I1 * (1 - mu * a) / N1 + S1 * I2 * mu / N1)
    dS2dt = -beta * (S2 * I1 * mu * a / N2 + S2 * I2 * (1-mu) / N2)
    dI1dt = beta * (S1 * I1 * (1 - mu * a) / N1 + S1 * I2 * mu / N1)  - gamma * I1 - s * r * N * I1 / (S1+I1)
    dI2dt = beta * (S2 * I1 * mu * a / N2 + S2 * I2 * (1-mu) / N2)  - gamma * I2
    dR1dt = gamma * I1  + s * r * N * I1 / (S1+I1)
    dR2dt = gamma * I2
    return dS1dt, dS2dt, dI1dt, dI2dt, dR1dt, dR2dt

# Initial conditions vector
y0 = S10, S20, I10, I20, R10, R20
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, mu))
S1, S2, I1, I2, R1, R2 = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S1/100000000, 'y', alpha=0.5, lw=2, label='S(t): Susceptible, testing')
ax.plot(t, I1/100000000, 'm', alpha=0.5, lw=2, label='I(t): Infected, testing')
ax.plot(t, S2/100000000, 'b', alpha=0.5, lw=2, label='S(t): Susceptible, no testing')
ax.plot(t, I2/100000000, 'r', alpha=0.5, lw=2, label='I(t): Infected, no testing')
ax.plot(t, R1/100000000, 'g', alpha=0.5, lw=2, label='R(t): Recovered/isolated with immunity')
ax.plot(t, R2/100000000, 'k', alpha=0.5, lw=2, label='R(t): Recovered with immunity')
# ax.plot(t, (S2+I2+R2)/100000000, 'k', alpha=0.5, lw=2, label='R(t): Recovered with immunity')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (100000s)')
ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()
#plt.savefig('masstest.png')

