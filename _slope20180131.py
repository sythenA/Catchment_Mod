
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

        x = np.zeros(grds)
        for i in range(1, grds):
            x[i] = dx*i
        self.x = x

        if init_cond:
            self.init_cond = np.zeros(grds)
            self.init_cond[0: len(init_cond)] = init_cond
        else:
            self.init_cond = np.zeros(grds)
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
        alpha = self.alpha
        self.q_out = list()

        t0 = min(self.rainData[:, 0])
        te = max(self.rainData[:, 0])

        c_t = t0  # Current time, accumulated in claculation process

        while c_t < te:
            depth0, c_t = self.depth_iterate1(depth0.copy(), c_t, self.dt)
            print c_t, alpha*depth0[-1]**(5.0/3)
            self.q_out.append([c_t, alpha*depth0[-1]**(5.0/3)])
            print self.get_rain(c_t)*3600.0*1000.0

    def get_rain(self, time):
        rainData = self.rainData
        for i in range(1, len(rainData)):
            if time > rainData[i-1, 0] and time <= rainData[i, 0]:
                return rainData[i, 1]
            elif time == rainData[i-1, 0]:
                return rainData[i-1, 1]

    def depth_iterate1(self, depth0, ct, dt):
        N = self.grds

        Courant_OK = False
        while Courant_OK != True:
            Courant_OK = self.courant_Check(depth0, ct, dt)
            if not Courant_OK:
                dt = 0.5*dt

        depth1 = np.zeros(len(depth0))
        for j in range(1, N):
            depth1[j], xl = self.depth_iterate2(depth0, j, ct, dt)
        depth1[-1] = depth1[-2]

        ct = ct + dt
        return depth1, ct

    def courant_Check(self, depth0, ct, dt):
        N = self.grds
        dt_Moved = np.zeros(N)
        m = 5.0/3
        for j in range(0, N):
            try:
                dt_Moved[j] = 0.5*self.alpha*m*(
                        depth0[j]**(m-1) +
                              (depth0[j] + self.get_rain(ct+dt))**(m-1))*dt
            except(TypeError):
                return False
        for j in range(0, N):
            if dt_Moved[j] >= 0.8*self.dx:
                return False

        return True

    def depth_iterate2(self, depth1, i, ct, dt):

        try:
            hp, xl1 = self.newton(depth1, i, ct, dt)
        except(ValueError):
            hp, xl1 = self.bisection(depth1, i, ct, dt)

        return hp, xl1

    def newton(self, depth1, i, ct, dt):
        # print 'Use Newton Method'
        dx = self.dx
        m = 5.0/3
        alpha = self.alpha
        hl0 = depth1[i-1] + 0.5*(depth1[i-1] - depth1[i])
        xl0 = 0.5*dx
        Err = 1.0
        while Err > 10**-3:
            hp1 = hl0 + self.get_rain(ct+ 0.5*dt)*dt
            hp2 = hl0 + self.get_rain(ct+dt)*dt
            trip = 1./6*alpha*m*dt*(hl0**(m-1) + 4*hp1**(m-1) + hp2**(m-1))
            f = dx - xl0 - trip
            # f = dx - xl0 - 1./2*alpha*m*(hl0**(m-1) + hp**(m-1))*dt
            dfdx = (-1 - 0.5*alpha*m*(m-1)*dt/dx*(
                depth1[i] - depth1[i-1])*(hp2**(m-2) + hl0**(m-2)))
            xl1 = xl0 - f/dfdx
            hl1 = depth1[i-1] + (xl1/dx)*(depth1[i] - depth1[i-1])

            Err = abs(xl1 - xl0)/xl0
            xl0 = xl1
            hl0 = hl1
            if np.isnan(xl1):
                raise ValueError

        hl1 = depth1[i-1] + (xl1/dx)*(depth1[i] - depth1[i-1])
        hp = hl1 + self.get_rain(ct+dt)*dt

        return hp, xl1

    def bisection(self, depth1, i, ct, dt):
        # print 'Use Bisection Method'
        dx = self.dx
        Err = 1.0
        x1 = dx
        x2 = 0.5*dx
        x3 = 0.0

        while Err > 10**-3:
            f1 = self.function(depth1, i, x1, ct, dt)
            f2 = self.function(depth1, i, x2, ct, dt)
            f3 = self.function(depth1, i, x3, ct, dt)
            if np.isnan(f2):
                input()

            if f1*f2 < 0:
                x3 = x2
                x2 = 0.5*(x1 + x3)
                Err = abs(x2 - x3)/x2

            elif f2*f3 < 0:
                x1 = x2
                x2 = 0.5*(x1 + x3)
                Err = abs(x2 - x1)/x2

        hp = depth1[i-1] + (x2/dx)*(depth1[i] - depth1[i-1])

        return hp, x2

    def function(self, depth1, i, xl, ct, dt):
        dx = self.dx
        alpha = self.alpha
        m = 5.0/3
        hl = depth1[i-1] + (xl/dx)*(depth1[i] - depth1[i-1])
        hp1 = hl + self.get_rain(ct+0.5*dt)*dt
        hp2 = hl + self.get_rain(ct+dt)*dt

        trip = 1./6*alpha*m*dt*(hl**(m-1) + 4*hp1*(m-1) + hp2**(m-1))
        f = dx - xl - trip
        # f = dx - xl - 1./2*alpha*m*(hl**(m-1) + hp**(m-1))*dt

        return f
