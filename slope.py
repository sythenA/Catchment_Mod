
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

        self.currentSitu = list()
        for i in range(0, self.grds):
            self.currentSitu.append([[self.x[i], 10**-10]])
        self.currentSitu[0][0][1] = 0.0

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
        while Courant_OK is not True:
            Courant_OK = self.courant_Check(depth0, ct, dt)
            if not Courant_OK:
                dt = 0.5*dt

        depth1 = np.zeros(len(depth0))
        for j in range(1, N):
            depth1[j], xl = self.depth_iterate2(depth0, j, ct, dt)
        depth1[-1] = depth1[-2]

        self.movingWave(ct, dt, depth1)

        ct = ct + dt
        return depth1, ct

    def courant_Check(self, depth0, ct, dt):
        N = self.grds
        dt_Moved = np.zeros(N)
        m = 5.0/3
        for j in range(0, N):
            try:
                dt_Moved[j] = 0.5*self.alpha*m*(
                    depth0[j]**(m-1) + (depth0[j] +
                                        self.get_rain(ct+dt))**(m-1))*dt
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

    def movingWave(self, ct, dt, depth1):
        currentSitu = self.currentSitu
        alpha = self.alpha
        m = 5.0/3
        for j in range(0, len(currentSitu)-1):
            outRegion = list()
            if currentSitu[j]:
                for k in range(0, len(currentSitu[j])):
                    rain = 0.5*(self.get_rain(ct) + self.get_rain(ct+dt))
                    hl = currentSitu[j][k][1]
                    xl = currentSitu[j][k][0]

                    if rain != 0.:
                        xl = xl + alpha/rain*((hl + rain*dt)**m - hl**m)
                    else:
                        xl = xl + alpha*m*hl**(m-1)*dt

                    currentSitu[j][k][0] = xl
                    currentSitu[j][k][1] = hl + rain*dt

                for g in range(0, len(currentSitu[j])):
                    if currentSitu[j][g][0] >= self.x[j+1]:
                        outRegion.append(g)
                if outRegion:
                    currentSitu[j] = [s for l, s in enumerate(currentSitu[j])
                                      if l not in outRegion]

            currentSitu[j].insert(0, [self.x[j], depth1[j]])

        currentSitu[-1].append([self.x[-1], depth1[-1]])
        currentSitu[-1].pop(0)

        self.currentSitu = currentSitu

    def situToArray(self):
        currentSitu = self.currentSitu
        situArray = list()

        for i in range(0, len(currentSitu)):
            if currentSitu[i]:
                for j in range(0, len(currentSitu[i])):
                    situArray.append(currentSitu[i][j])

        situArray = np.array(situArray)
        return situArray

    def getDHDX(self, x):
        situArray = self.situToArray()
        for i in range(1, len(situArray)):
            if x > situArray[i-1, 0] and x <= situArray[i, 0]:
                return (situArray[i, 1] - situArray[i-1, 1])/(
                    situArray[i, 0] - situArray[i-1, 0])
            elif x == situArray[i-1, 0]:
                return (situArray[i, 1] - situArray[i-2, 1])/(
                    situArray[i, 0] - situArray[i-2, 0])

    def getDepth(self, x):
        situArray = self.situToArray()
        """
        for i in range(1, len(situArray)):
            if x > situArray[i-1, 0] and x <= situArray[i, 0]:
                h = situArray[i-1, 1] + (x - situArray[i-1, 0])*(
                    situArray[i, 1] - situArray[i-1, 1])/(
                    situArray[i, 0] - situArray[i-1, 0])
                return h
            elif x == situArray[i-1, 0]:
                return situArray[i-1, 1]"""
        return np.interp(x, situArray[:, 0], situArray[:, 1])

    def newton(self, depth1, i, ct, dt):
        # print 'Use Newton Method'
        dx = self.dx
        m = 5.0/3
        alpha = self.alpha
        hl0 = depth1[i-1] + 0.5*(depth1[i-1] - depth1[i])
        xl0 = 0.5*dx
        Err = 1.0
        counter = 0
        while Err > 10**-3:
            hp = hl0 + self.get_rain(ct+dt)*dt
            if self.get_rain(ct+dt) != 0.:
                f = dx - xl0 - alpha/self.get_rain(ct+dt)*(hp**m - hl0**m)
                """
                dfdx = -1 - alpha/self.get_rain(ct+dt)*m*hp**(m-1)*(
                    depth1[i] - depth1[i-1])"""
                dfdx = -1 - alpha/self.get_rain(
                    ct+dt)*m*hp**(m-1)*self.getDHDX(self.x[i-1]+xl0)
            else:
                f = dx - xl0 - alpha*m*hl0**(m-1)*dt
                """
                dfdx = -1 - alpha*m*(m-1)*hl0**(m-2)*dt*(
                    depth1[i] - depth1[i-1])"""
                dfdx = -1 - alpha*m*(m-1)*hl0**(m-2)*self.getDHDX(
                    self.x[i-1]+xl0)

            xl1 = xl0 - f/dfdx
            hl1 = self.getDepth(xl1 + self.x[i-1])

            Err = abs(xl1 - xl0)/xl0
            xl0 = xl1
            hl0 = hl1
            counter += 1
            if np.isnan(xl1):
                raise ValueError
            elif counter > 100:
                raise ValueError

        # hl1 = depth1[i-1] + (xl1/dx)*(depth1[i] - depth1[i-1])
        hl1 = self.getDepth(xl1 + self.x[i-1])
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

        hp = self.getDepth(x2 + self.x[i-1])
        hp = hp + self.get_rain(ct+dt)*dt

        return hp, x2

    def function(self, depth1, i, xl, ct, dt):
        dx = self.dx
        alpha = self.alpha
        m = 5.0/3

        hl = self.getDepth(self.x[i-1]+xl)
        hp = hl + self.get_rain(ct+dt)*dt
        if self.get_rain(ct+dt) != 0:
            f = dx - xl - alpha/self.get_rain(ct+dt)*(hp**m - hl**m)
        else:
            f = dx - xl - alpha*m*dt*hl**(m-1)

        return f
