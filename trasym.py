# the original python code of the SIR model is https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/ by Christian Hill

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pylab

# Total population: N
N = 100000000
# Initial number of infected, quarantined, asymptomatic and recovered individuals, I0, Q0, A0 and R0.
I0, Q0, A0, R0 = 100, 0, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - Q0 - A0 - R0
# Contact rate, beta1 (symptomatic) beta2 (asymptomatic), and mean recovery rate, gamma1 (symptomatic) gamma2(asymptomatic) (in 1/days), quarantine rate, delta.
beta1, beta2, gamma1, gamma2, delta = 0.25, 0.25, 0.2, 0.2, 1
# A grid of time points (in days)
t = np.linspace(0, 100, 100)

# The SIQAR model differential equations.
# Symptomatic patients are detedted and quarantined. Asymptomatic carriers remain undetected.
# The maximum tracing capacity is set to 1000

def deriv(y, t, N, beta1, beta2, gamma1, gamma2, delta):
    S, I, Q, A, R = y
    dSdt = -(beta1+beta2) * S * (I+A) / N
    dIdt = beta1 * S * (I+A) / N - gamma1 * I - min(delta * I, 1000)
    dQdt = min(delta * I, 1000) - gamma1 * Q
    dAdt = beta2 * S * (I+A) / N - gamma2 * A
    dRdt = gamma1 * (I+Q) + gamma2 * A
    return dSdt, dIdt, dQdt, dAdt, dRdt

# Initial conditions vector
y0 = S0, I0, Q0, A0, R0
# Integrate the SIQAR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta1, beta2, gamma1, gamma2, delta))
S, I, Q, A, R = ret.T

# Plot the data
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S/100000000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I/100000000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, Q/100000000, 'y', alpha=0.5, lw=2, label='Quarantined')
ax.plot(t, A/100000000, 'k', alpha=0.5, lw=2, label='Asymptomatic')
ax.plot(t, R/100000000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.plot(t, Q/(I+Q+A), 'm', alpha=0.5, lw=2, label='Detection rate')
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
