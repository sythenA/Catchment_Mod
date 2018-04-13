
import numpy as np
from math import sqrt, ceil


class channel:
    def __init__(self, length, B, N, n, S, dt, Lc, **kwargs):
        self.alpha = 1./n*sqrt(S)
        self.B = B
        self.grds = N
        self.dx = length/(N-1)
        self.dt = dt
        self.Lc = Lc
        self.Ls = length

        x = np.zeros(N)
        for i in range(1, N):
            x[i] = self.dx*i
        self.x = x

        try:
            inflow_file = kwargs['inflow_file']
            inflow = self.readText(inflow_file)
            self.cat_inflow = inflow
        except(IOError):
            pass

    def rainParse(self, rainData, unit='sec'):
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

    def readText(self, fileName):
        f = open(fileName, 'r')
        data = list()
        for line in f:
            container = list()
            for num in line.split():
                container.append(float(num))
            data.append(container)

        return np.array(data)

    def get_rain(self, time):
        rainData = self.rainData
        for i in range(1, len(rainData)):
            if time > rainData[i-1, 0] and time <= rainData[i, 0]:
                return rainData[i, 1]
            elif time == rainData[i-1, 0]:
                return rainData[i-1, 1]

    def get_catInflow(self, time):
        cat_inflow = self.cat_inflow
        """
        for i in range(1, len(cat_inflow)):
            if time > cat_inflow[i-1, 0] and time <= cat_inflow[i, 0]:
                outflow = cat_inflow[i-1, 1] + (time - cat_inflow[i-1, 0])*(
                    cat_inflow[i, 1] - cat_inflow[i-1, 1])
                return outflow
            elif time == cat_inflow[i-1, 0]:
                return cat_inflow[i-1, 1]"""
        return np.interp(time, cat_inflow[:, 0], cat_inflow[:, 1])

    def courant_Check(self, depth0, ct, dt):
        N = self.grds
        dt_Moved = np.zeros(N)
        m = 5.0/3
        for j in range(0, N):
            try:
                ql = (self.get_catInflow(ct+dt) +
                      self.get_catInflow(ct))/self.B
                hp = depth0[j] + ql*dt
                Rp = self.B*hp/(self.B + 2*hp)
                hl = depth0[j]
                Rl = self.B*hl/(self.B + 2*hl)
                if ql != 0:
                    dt_Moved[j] = self.alpha/ql*(Rp**(m-1)*hp -
                                                 Rl**(m-1)*hl)
                elif ql == 0:
                    dt_Moved[j] = self.alpha*m*Rl**(m-1)*dt

            except(TypeError):
                return False

        for j in range(0, N):
            if dt_Moved[j] >= 0.8*self.dx:
                return False

        return True

    def run(self):
        maxTime = max(self.rainData[:, 0])
        minTime = min(self.rainData[:, 0])
        alpha = self.alpha
        self.q_out = list()

        depth0 = np.zeros(self.grds)

        self.currentSitu = list()
        for i in range(0, self.grds):
            self.currentSitu.append([[self.x[i], 10**-10]])
        self.currentSitu[0][0][1] = 0.0

        ct = minTime  # Current Time
        while ct <= maxTime:
            depth0, ct = self.depth_iterate1(depth0.copy(), ct)
            HRadius = (self.B*depth0[-1]/(self.B + 2*depth0[-1]))
            flowrate = alpha*HRadius**(2./3)*self.B*depth0[-1]
            print ct, flowrate
            self.q_out.append([ct, flowrate])

    def depth_iterate1(self, depth0, ct):
        dt = self.dt
        N = self.grds

        print 'current time: %f' % ct

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
                    ql = (self.get_catInflow(ct+dt) +
                          self.get_catInflow(ct))/self.B
                    hl = currentSitu[j][k][1]
                    xl = currentSitu[j][k][0]
                    Rl = (self.B*hl/(self.B + 2*hl))
                    Rp = (self.B*(hl + ql*dt)/(self.B + 2*(hl + ql*dt)))

                    if ql != 0.:
                        xl = xl + alpha/ql*(Rp**(m-1)*(hl + ql*dt) -
                                            Rl**(m-1)*hl)
                    else:
                        xl = xl + alpha*m*Rl**(m-1)*dt

                    currentSitu[j][k][0] = xl
                    currentSitu[j][k][1] = hl + ql*dt

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
        ql = (self.get_catInflow(ct+dt) + self.get_catInflow(ct))/self.B
        while Err > 10**-5:
            try:
                hp = hl0 + ql*dt
                Rp = (self.B*hp/(self.B + 2*hp))
                Rl0 = (self.B*hl0/(self.B + 2*hl0))
                if ql != 0.:
                    f = dx - xl0 - alpha/ql*(Rp**(m-1)*hp - Rl0**(m-1)*hl0)
                    dfdx = -1 - alpha/ql*m*Rp**(m-1)*self.getDHDX(
                        self.x[i-1] + xl0)
                elif ql == 0.:
                    f = dx - xl0 - alpha*m*Rl0**(m-1)*dt
                    dfdx = -1 - alpha*m*(m-1)*Rl0**(m-2)*dt*self.getDHDX(
                        self.x[i-1] + xl0)
            except(TypeError, ValueError):
                raise ValueError

            xl1 = xl0 - f/dfdx
            hl1 = self.getDepth(self.x[i-1] + xl1)

            Err = abs(xl1 - xl0)/xl0
            xl0 = xl1
            hl0 = hl1
            if np.isnan(xl1):
                raise ValueError

            counter += 1
            if counter > 50:
                raise ValueError

        hp = hl1 + ql*dt

        return hp, xl1

    def bisection(self, depth1, i, ct, dt):
        # print 'Use Bisection Method'
        dx = self.dx
        Err = 1.0
        x1 = dx
        x2 = 0.5*dx
        x3 = 0.0

        while Err > 10**-5:
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

        ql = (self.get_catInflow(ct+dt) + self.get_catInflow(ct))/self.B
        hp = self.getDepth(self.x[i-1] + x2) + ql*dt

        return hp, x2

    def function(self, depth1, i, xl, ct, dt):
        dx = self.dx
        alpha = self.alpha
        m = 5.0/3
        hl = self.getDepth(self.x[i-1] + xl)
        Rl = (self.B*hl/(self.B + 2*hl))

        ql = (self.get_catInflow(ct+dt) + self.get_catInflow(ct))/self.B
        hp = hl + ql*dt
        Rp = (self.B*hp/(self.B + 2*hp))
        if ql != 0.:
            f = dx - xl - alpha/ql*(Rp**(m-1)*hp - Rl**(m-1)*hl)
        elif ql == 0:
            f = dx - xl - alpha*m*Rl**(m-1)*dt

        return f
