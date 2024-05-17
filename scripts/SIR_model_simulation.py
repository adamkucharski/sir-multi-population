import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define the SIR model differential equations.
def deriv(y, t, N, beta, gamma, alpha1, alpha2, alpha3):
    S, I, R, S2, I2, R2, S3, I3, R3 = y
    foi1 = beta * (I/N + alpha1*I2/N + alpha3*I3/N)
    foi2 = beta * (I2/N + alpha1*I/N + alpha2*I3/N)
    foi3 = beta * (I3/N + alpha2*I2/N + alpha3*I/N)
    dSdt = -S * foi1
    dIdt = S * foi1 - gamma * I
    dRdt = gamma * I
    dS2dt = -S2 * foi2
    dI2dt = S2 * foi2 - gamma * I2
    dR2dt = gamma * I2
    dS3dt = -S3 * foi3
    dI3dt = S3 * foi3 - gamma * I3
    dR3dt = gamma * I3
    return dSdt, dIdt, dRdt, dS2dt, dI2dt, dR2dt, dS3dt, dI3dt, dR3dt

# Main function to run the SIR model simulation
def run_simulation():
    # Total population, N.
    N = 10000
    # Initial number of infected and recovered individuals, I0 and R0.
    I0, R0 = 1, 0
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0
    # Contact rates, beta, and mean recovery rate, gamma, (in 1/days).
    beta, gamma = 0.4, 1./10 
    # Relative contact rates between populations
    alpha1, alpha2, alpha3 = 0.0001, 0.0001, 0.000
    # A grid of time points (in days)
    t = np.linspace(0, 160, 160)

    # Initial conditions vector
    y0 = S0, I0, R0, S0, I0, R0, S0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, beta, gamma, alpha1, alpha2, alpha3))
    S, I, R, S2, I2, R2, S3, I3, R3 = ret.T

    # Plot the data on three separate curves for S(t), I(t) and R(t)
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
    ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
    ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered with immunity')
    ax.set_xlabel('Time /days')
    ax.set_ylabel('Number (10000s)')
    ax.set_ylim(0,10000)
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.show()

if __name__ == "__main__":
    run_simulation()
