# the original python code of the SIR model is https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/ by Christian Hill

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Total population: N, influenza-like disease Inf
N = 100000000
Inf = 1000
# Initial number of infected, quarantined, asymptomatic and recovered individuals, I0, Q0, A0 and R0
I0, Q0, A0, R0, Rq0= 0, 0, 4000, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - Q0 - A0 - R0 - Rq0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days), quarantine rate, mu,  maximum number of tests per day, Test.
beta1, beta2, gamma1, gamma2, delta1, delta2, Test = 0.25, 0.25, 0.2, 0.2, 0.01, 0.05, 2000
# A grid of time points (in days)
t = np.linspace(0, 40, 40)

# The SIR model differential equations.
def deriv(y, t, N, beta1, beta2, gamma1, gamma2, delta1, delta2, Test):
    S, I, Q, A, R, Rq = y
    dSdt = -(beta1+beta2) * S * (I+A) / N
    dIdt = beta1 * S * (I+A) / N - Test * I / (2 * (I + Inf)) - gamma1 * I
    dQdt = Test * I / (2 * (I + Inf)) + Test * I / (2 * (A + Inf)) - gamma1 * Q
    dAdt = beta2 * S * (I+A) / N - Test * I / (2 * (I + Inf)) - gamma2 * A
    dRdt = gamma1 * I + gamma2 * A
    dRqdt = gamma1 * Q
    return dSdt, dIdt, dQdt, dAdt, dRdt, dRqdt

# Initial conditions vector
y0 = S0, I0, Q0, A0, R0, Rq0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta1, beta2, gamma1, gamma2, delta1, delta2, Test))
S, I, Q, A, R, Rq = ret.T

# Plot the data on three separate curves for S(t), I(t), Q(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, S/100000000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, (I+A)/10000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, (Q+Rq)/100000000,    'm', alpha=0.5, lw=2, label='Positive')
ax.plot(t, (Test * I / (2 * (I + Inf)) + Test * I / (2 * (A + Inf)))/10000, 'y', alpha=0.5, lw=2, label='New positive (x 10000)')
ax.plot(t, (I / (I + Inf)), 'k', alpha=0.5, lw=2, label='Positivity rate')
#ax.plot(t, A/1000000, 'k', alpha=0.5, lw=2, label='Asymptomatic')
#ax.plot(t, R/100000000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
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
#plt.savefig('tracingl3.png')

