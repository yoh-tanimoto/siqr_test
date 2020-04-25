# The original code of the SIR model using SciPy odeint by Christian Hill is here https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/
import numpy as np
import matplotlib.pyplot as plt

# Total population: N, influenza-like disease Flu
N = 100000000
Flu = 1000
# Initial number of infected, quarantined, asymptomatic and recovered individuals, I0, Q0, A0 and R0
I0, Q0, A0, R0, Rq0= 0, 0, 500, 0, 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - Q0 - A0 - R0 - Rq0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days), quarantine rate, mu,  maximum number of tests per day, Test.
beta1, beta2, gamma1, gamma2, delta, Test = 0.25, 0.25, 0.2, 0.2, 1, 700
# A grid of time points (in days)
TIME = 50
t = np.linspace(0, TIME, TIME)
S, I, A, Q, Rq, R = np.empty(TIME), np.empty(TIME), np.empty(TIME), np.empty(TIME), np.empty(TIME), np.empty(TIME)

S[0], I[0], A[0], Q[0], Rq[0], R[0] = S0, I0, A0, Q0, Rq0, R0
Np = 0 # new positive

# The SIQR model differential equations.

for x in range (0,TIME-1):
  S[x+1]  = S[x]  - (beta1+beta2) * (I[x]+A[x]) * S[x] / N
  I[x+1]  = I[x]  + beta1 * (I[x]+A[x]) * S[x] / N - Test * I[x] / (2 * (I[x] + Flu)) - gamma1 * I[x]
  Q[x+1]  = Q[x]  + Test * I[x] / (2 * (I[x] + Flu)) + delta * Np - gamma1 * Q[x]
  A[x+1]  = A[x]  + beta2 * S[x] * (I[x]+A[x]) / N - delta * Np - gamma2 * A[x] # new symptomatic is associated with an asymptomatic contact
  R[x+1]  = R[x]  + gamma1 * I[x] + gamma2 * A[x]
  Rq[x+1] = Rq[x] + gamma1 * Q[x]
  Np = Test * I[x] / (2 * (I[x] + Flu)) # new symptomatic positive of the day before

# Plot the data on three separate curves for S(t), I(t), Q(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
#ax.plot(t, S/100000000, 'b', alpha=0.5, lw=2, label='Susceptible')
#ax.plot(t, I/100000000, 'r', alpha=0.5, lw=2, label='Infected')
#ax.plot(t, R/100000000, 'g', alpha=0.5, lw=2, label='Recovered')


ax.plot(t, (I+A)/10000, 'r', alpha=0.5, lw=2, label='Undetected infected (x 10000)')
ax.plot(t, (Q+Rq)/10000,    'm', alpha=0.5, lw=2, label='Positive confirmed (x 10000)')
ax.plot(t, (Test * I / (2 * (I + Flu)) + Test * I / (2 * (A + Flu)))/10000, 'y', alpha=0.5, lw=2, label='New positive (x 10000)')
ax.plot(t, (I / (I + Flu)), 'k', alpha=0.5, lw=2, label='Positivity rate')

#ax.plot(t, A/1000000, 'k', alpha=0.5, lw=2, label='Asymptomatic')
#ax.plot(t, R/100000000, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number / ' + str(N))
ax.set_ylim(0,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.title("SIAQR model as a difference equation with fixed number of tests. \n A: asymptomatic.\n Flu = " + str(Flu) + ", Test = " + str(Test) + "/day, A_0 = " + str(A0) + ", delta = " + str(delta))
plt.show()
#plt.savefig('siqar_test3.png')

