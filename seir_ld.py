# the original python code of the SIR model is https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/ by Christian Hill

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population: N
N = 100000000
# Initial number of exposed, infected, and recovered individuals, E0, I0, and R0.
E0, I0, R0 = 10000, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - E0 - I0 - R0
# mean recovery rate, gamma, (in 1/days), the inverse latent period, sigma.
gamma, sigma = 0.2, 0.5
# A grid of time points (in days)
t = np.linspace(0, 100, 100)

# The time-dependent contact rate, beta, modelling a lockdown or a similar sudden change of contact rate.
def beta(s):
    if s < 50: 
       return 0.5
    else: 
       return 0.10

# The SEIR model differential equations.
def deriv(y, t, N, gamma, sigma):
    S, E, I, R = y
    dSdt = -beta(t) * S * I / N
    dEdt = beta(t) * S * I / N - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, gamma, sigma))
S, E, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t), Q(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S/100000000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/100000000, 'r', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/100000000, 'y', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/100000000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.plot(t, (sigma * E)/10000000, 'm', alpha=0.5, lw=2, label='New positive (x 10)')
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
plt.title("SEIR model with beta = 0.5 (t<50), 0.1 (t>=50)")
plt.show()

