# The original code of the SIR model using SciPy odeint by Christian Hill is here https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/

import numpy as np
import matplotlib.pyplot as plt

# Total population: N, influenza-like disease Flu
N = 100000000
Flu = 1000
# Initial number of infected, quarantined, asymptomatic and recovered individuals, I0, Q0, A0 and R0
I0, Q0, A0, R0, Rq0= 1, 0, 0, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - Q0 - A0 - R0 - Rq0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days), quarantine rate, mu,  maximum number of tests per day, Test.
beta, beta2, gamma, gamma2, delta, Test = 0.25, 0.25, 0.2, 0.2, 0.01, 2000
# A grid of time points (in days)
TIME = 100
t = np.linspace(0, TIME, TIME)
S, I, A, Q, Rq, R = np.empty(TIME), np.empty(TIME), np.empty(TIME), np.empty(TIME), np.empty(TIME), np.empty(TIME)

S[0], I[0], Q[0] = S0, I0, Q0

# The SIR model differential equations.

for x in range (0,TIME-1):
  S[x+1] = S[x] - beta * I[x] * S[x] / N
  I[x+1] = I[x] + beta * I[x] * S[x] / N - gamma * I[x]
  R[x+1] = R[x] + gamma * I[x]

#    dSdt = -(beta1+beta2) * S * (I+A) / N
#    dIdt = beta1 * S * (I+A) / N - Test * I / (2 * (I + Flu)) - gamma1 * I
#    dQdt = min(delta *Q, Test / 2) + Test * I / (2 * (A + Flu)) - gamma1 * Q
#    dAdt = beta2 * S * (I+A) / N - min(delta *Q, Test / 2) - gamma2 * A
#    dRdt = gamma1 * I + gamma2 * A
#    dRqdt = gamma1 * Q



# Plot the data on three separate curves for S(t), I(t), Q(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S/100000000, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, I/100000000, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/100000000, 'g', alpha=0.5, lw=2, label='Recovered')


#ax.plot(t, (I+A)/10000, 'r', alpha=0.5, lw=2, label='Infected (x 10000)')
#ax.plot(t, (Q+Rq)/100000000,    'm', alpha=0.5, lw=2, label='Positive')
#ax.plot(t, (Test * I / (2 * (I + Flu)) + Test * I / (2 * (A + Flu)))/10000, 'y', alpha=0.5, lw=2, label='New positive (x 10000)')
#ax.plot(t, (I / (I + Flu)), 'k', alpha=0.5, lw=2, label='Positivity rate')

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
plt.title("SIR model as a difference equation. \nFlu = " + str(Flu) + ", Test = " + str(Test) + "/day, A_0 = " + str(A0))
plt.show()
#plt.savefig('import_3.png')

