# the original python code of the SIR model is https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/ by Christian Hill

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population: N
N = 100000000
# Initial number of infected, quarantined, and recovered individuals, I0 and R0.
I0, R0 = 100, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days), waning rate, xi
beta, gamma, xi = 0.3, 0.2, 0.01
# A grid of time points (in days)
t = np.linspace(0, 500, 500)

# The SIR model differential equations.
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    dSdt = -beta * S * I / N + xi * R
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I - xi * R
    return dSdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S/100000000, 'b', alpha=0.5, lw=2, label='S(t): Susceptible')
ax.plot(t, I/100000000, 'r', alpha=0.5, lw=2, label='I(t): Infected')
ax.plot(t, R/100000000, 'g', alpha=0.5, lw=2, label='R(t): Recovered with immunity')
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
#plt.savefig('sirs.png')

