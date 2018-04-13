
import numpy as np
from math import sqrt, ceil


class numerical_sol:
    def __init__(self, slope, n, length, grds, dt, init_cond=None):
        dx = length/(grds - 1)
        self.length = length
        self.dx = dx
        self.grds = grds
        self.dt = dt
        self.alpha = 1.0/n*sqrt(slope)

        x = np.zeros(grds+1)
        for i in range(1, grds+1):
            x[i] = dx*i
        self.x = x

        if init_cond:
            self.init_cond = np.zeros(grds+1)
            self.init_cond[0: len(init_cond)] = init_cond
            self.init_cond[-1] = self.init_cond[-2]
        else:
            self.init_cond = np.zeros([grds+1])
            for i in range(0, len(self.init_cond)):
                self.init_cond[i] = 10**-9

    def rainfallParse(self, rainData, unit='sec'):
        def inBetween(time, rainData):
            for i in range(1, len(rainData)):
                if time > rainData[i-1, 0] and time <= rainData[i, 0]:
                    return rainData[i, 1]
                elif time == rainData[i-1, 0]:
                    return rainData[i-1, 1]

            raise ValueError('Did not find matching time in input rain data, \
input time: %f', time)

        rainData = np.array(rainData)
        startingTime = min(rainData[:, 0])
        if unit == 'sec':
            timeLength = max(rainData[:, 0]) - min(rainData[:, 0])
        elif unit == 'min':
            timeLength = (max(rainData[:, 0]) - min(rainData[:, 0]))*60
            rainData[:, 0] = rainData[:, 0]*60
        elif unit == 'hour':
            timeLength = (max(rainData[:, 0]) - min(rainData[:, 0]))*3600
            rainData[:, 0] = rainData[:, 0]*3600
        if timeLength % 60.0 == 0:
            timeGrids = int(ceil(timeLength/60.0)) + 1
        else:
            timeGrids = int(ceil(timeLength/60.0))

        rainParse = np.zeros([timeGrids, 2])
        for i in range(0, timeGrids):
            rainParse[i, 0] = startingTime + 60.0*i
            rainParse[i, 1] = inBetween(rainParse[i, 0], rainData)

        self.rainData = rainParse

    def overland(self):
        depth0 = self.init_cond  # Initial condition of water stage
        N = self.grds - 1
        dx = self.dx
        alpha = self.alpha

        t0 = min(self.rainData[:, 0])
        te = max(self.rainData[:, 0])

        c_t = t0  # Current time, accumulated in claculation process

        dep_Mat = np.zeros(N)
        h2dx2 = np.zeros(self.grds+1)
        while c_t < te:
            # Solve Ax = B
            # where x = [S_1, ..., S_N]
            # B Maxtrix
            dep_Mat[0] = (depth0[2] + depth0[0] - 2*depth0[1])*6/dx**2
            dep_Mat[-1] = (depth0[-1] + depth0[-3] - 2*depth0[-2])*6/dx**2
            for i in range(1, N-1):
                dep_Mat[i] = 6*(depth0[i+1] + depth0[i-1] - 2*depth0[i])/dx**2

            # print 'dep_Mat'
            # print dep_Mat[-1], depth0[-1], depth0[-2], depth0[-3]

            Sol_Mat = np.zeros([N, N])
            Sol_Mat[0, 0:3] = [4, 0, 0]
            Sol_Mat[N-1, N-3: N] = [0, 0, 4]
            for i in range(1, N-1):
                Sol_Mat[i, i] = 2.0
                Sol_Mat[i, i-1] = 1.0
                Sol_Mat[i, i+1] = 1.0

            Sol_Mat, dep_Mat = self.solveS1(Sol_Mat, dep_Mat)
            Sol_Mat, dep_Mat = self.solveS2(Sol_Mat, dep_Mat)
            # print Sol_Mat, dep_Mat

            h2dx2[1:-1] = dep_Mat[:]
            h2dx2[0] = 2*dep_Mat[0] - dep_Mat[1]
            h2dx2[-1] = 2*dep_Mat[-1] - dep_Mat[-2]

            depth0, c_t = self.depth_iterate1(h2dx2, depth0.copy(), c_t)
            print c_t, alpha*depth0[-1]**(5.0/3)
            print depth0

    # First sweep
    def solveS1(self, Sol_Mat, dep_Mat):
        N = len(dep_Mat)
        Sol_Mat[0, :] = Sol_Mat[0, :]/Sol_Mat[0, 0]
        dep_Mat[0] = dep_Mat[0]/Sol_Mat[0, 0]

        Sol_Mat[N-1, :] = Sol_Mat[N-1, :]/Sol_Mat[N-1, N-1]
        dep_Mat[N-1] = dep_Mat[N-1]/Sol_Mat[N-1, N-1]

        for i in range(1, N-1):
            Sol_Mat[i, :] = Sol_Mat[i, :] - Sol_Mat[i-1, :]
            dep_Mat[i] = dep_Mat[i] - dep_Mat[i-1]
            Sol_Mat[i, :] = Sol_Mat[i, :]/Sol_Mat[i, i]
            dep_Mat[i] = dep_Mat[i]/Sol_Mat[i, i]

        """
        Sol_Mat[N-1, :] = Sol_Mat[N-1, :] - Sol_Mat[N-2, :]
        dep_Mat[N-1] = dep_Mat[N-1] - dep_Mat[N-2]
        Sol_Mat[N-1, :] = Sol_Mat[N-1, :]/Sol_Mat[N-1, N-1]
        dep_Mat[N-1] = dep_Mat[N-1]/Sol_Mat[N-1, N-1]"""

        return Sol_Mat, dep_Mat

    # Second sweep
    def solveS2(self, Sol_Mat, dep_Mat):
        N = len(dep_Mat)
        for i in range(N-2, 0, -1):
            Sol_Mat[i, :] = Sol_Mat[i, :] - Sol_Mat[i, i+1]*Sol_Mat[i+1, :]
            dep_Mat[i] = dep_Mat[i] - Sol_Mat[i, i+1]*dep_Mat[i+1]

        Sol_Mat[0, :] = Sol_Mat[0, :] - Sol_Mat[0, 2]*Sol_Mat[2, :]
        dep_Mat[0] = dep_Mat[0] - Sol_Mat[0, 2]*dep_Mat[2]

        return Sol_Mat, dep_Mat

    def get_rain(self, time):
        rainData = self.rainData
        for i in range(1, len(rainData)):
            if time > rainData[i-1, 0] and time <= rainData[i, 0]:
                return rainData[i, 1]
            else:
                return rainData[i-1, 1]

    def depth_iterate1(self, Sol_Mat, depth0, ct):
        # print depth0
        dt = self.dt
        N = self.grds

        depth1 = depth0.copy()
        for j in range(1, N):
            depth1[j], xl = self.depth_iterate2(Sol_Mat, depth0, j, ct, dt)
        depth1[-1] = depth1[-2]

        ct = ct + dt
        return depth1, ct

    def depth_iterate2(self, Sol_Mat, depth1, i, ct, dt):
        # print "depth1 %f" % depth1[i]
        alpha = self.alpha
        dx = self.dx
        x = self.x
        m = 5.0/3
        Err = 1.0

        o_hp = depth1[i]
        hp = depth1[i]
        hl = depth1[i-1]
        xl = x[i-1]
        while Err > 10**-6:
            A = (Sol_Mat[i+1] - Sol_Mat[i])/6.0/dx
            B = Sol_Mat[i]/2.0
            C = ((depth1[i+1] - depth1[i])/dx -
                 (2*dx*Sol_Mat[i] + dx*Sol_Mat[i+1])/6)
            D = depth1[i]

            hpl = alpha*m*(o_hp**(m-1) + hl**(m-1))/2
            omega = hpl*dt/dx

            # print 'omega: {%f}', omega
            # print 'hpl: {%f}', hpl
            if np.isnan(omega):
                print hp, hl

            hl = A*((1-omega)*dx)**3 + B*((1-omega)*dx)**2 + C*(1-omega)*dx + D
            # print "hl %f" % hl
            if hl < 0:
                print self.get_rain(ct), self.get_rain(ct+dt)
                print A, B, C, D
                print omega, i
                print hl
                print hp
                input()

            hp = hl + dt/2.*(self.get_rain(ct) + self.get_rain(ct+dt))
            xl = x[i] - alpha*m*dt/2*(hp**(m-1) + hl**(m-1))

            Err = abs(o_hp - hp)/o_hp
            # print Err, xl, hp

            o_hp = hp
        """
        if i == 1:
            hp = hl
        """

        return hp, xl
