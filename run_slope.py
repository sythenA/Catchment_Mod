
import slope
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from .analytic import slope_analytic as sa


ns = slope.numerical_sol(0.01, 0.1, 100.0, 201, 10.0)
ns.rainfallParse([[0.0, 100.0/3600.0/1000], [1800.0, 100.0/3600.0/1000],
                  [2600.0, 0.0]])
ns.overland()
out_q = ns.q_out
out_q = np.array(out_q)

ana_rise = sa.rising_phase(0.01, 0.1, 100.0, 100.0, 1800.0)
ana_fall = sa.falling_phase(0.01, 0.1, 100.0, 100.0, 1800.0, 2600.0)

ana_out = np.concatenate((ana_rise, ana_fall), axis=0)

ana_plot_c = list()
for j in range(0, 1801, 100):
    ana_plot_c.append(ana_out[j])

for j in range(1802, 2602, 100):
    ana_plot_c.append(ana_out[j])

ana_plot_c = np.array(ana_plot_c)

plt.plot(out_q[:, 0], out_q[:, 1], label='Numeric')
plt.plot(ana_plot_c[:, 0], ana_plot_c[:, 1], 'o', label='Analytic',
         markerfacecolor='none', markeredgecolor='r')
plt.legend(loc=2)
plt.show()
