from channel import channel
import numpy as np
import matplotlib.pyplot as plt
import slope
from .analytic import channel_analytic


def writeFile(fileName, Data):
    f = open(fileName, 'w')
    for i in range(0, len(Data)):
        line = ''
        for j in range(0, len(Data[i])):
            line = line + '{:>15}'.format(Data[i][j])
            line += '    '
        line = line[:-1]
        line += '\n'
        f.write(line)


def toString(mat):
    dataMatrix = list()
    for i in range(0, len(mat)):
        line = list()
        for j in range(0, len(mat[i])):
            line.append(str(mat[i][j]))
        dataMatrix.append(line)

    return dataMatrix


"""
cat_Slope = slope.numerical_sol(0.01, 0.1, 100.0, 201, 5.0)
cat_Slope.rainfallParse([[0.0, 100.0/3600.0/1000], [1800.0, 100.0/3600.0/1000],
                        [5600.0, 0.0/3600.0/1000]])
cat_Slope.overland()
cat_Outflow = toString(cat_Slope.q_out)
writeFile('catchment_Inflow.txt', cat_Outflow)"""

catchment = channel(200.0, 1.0, 101, 0.02, 0.005, 5.0, 100.0,
                    inflow_file='catchment_Inflow.txt')
catchment.rainParse([[0.0, 100.0/3600.0/1000], [1800.0, 100.0/3600.0/1000],
                    [5600.0, 0.0/3600.0/1000]])
catchment.run()

a_ch = channel_analytic.ana_channel(100.0, 1800.0, 5400.0, 5.0, 100., 200.,
                                    0.1, 0.02, 0.01, 0.005, 1.0)
a_ch.run()
a_outQ = np.array(a_ch.outQ)

q_out = np.array(catchment.q_out)

plt.plot(q_out[:, 0], q_out[:, 1], label='Numerical')
plt.plot(a_outQ[:, 0], a_outQ[:, 1], label='Analytical')
plt.xlabel('Time (s)')
plt.ylabel('flowrate ($m^3/s$)')
plt.title('Stream Outflow')
plt.legend(loc=1)
plt.savefig('Stream.png')
plt.show()
