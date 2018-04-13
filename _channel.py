
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
        except:
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

    def courant_Check(self, flowrate0, ct, dt):
        N = self.grds
        dt_Moved = np.zeros(N)
        m = 3.0/2
        for j in range(0, N):
            try:
                ql = (self.get_catInflow(ct+dt) + self.get_catInflow(ct))
                if ql != 0:
                    Q0 = flowrate0[j]
                    y0 = self.getH(Q0)
                    A0 = self.B*y0

                    Ap = A0 + ql*dt
                    yp = Ap/self.B

                    dt_Moved[j] = self.alpha/ql*self.B*(yp**m - y0**m)
                elif ql == 0:
                    y0 = self.getH(flowrate0[j])
                    A0 = self.B*y0
                    dt_Moved[j] = self.alpha*m*y0**(m-1)*dt

            except(TypeError):
                return False

        for j in range(0, N):
            if dt_Moved[j] >= 0.9*self.dx:
                return False

        return True

    def getH(self, Qa):
        y0 = (Qa/self.alpha/self.B)**(3./5)
        HRadius = (self.B/(self.B + 2*y0))**(2./3)
        Q0 = self.alpha*HRadius*self.B*y0
        Err = abs(Qa - Q0)

        while Err > 1.0E-8:
            dfdy = (5*self.B + 6*y0)/3/y0/(self.B + 2*y0)
            y1 = y0 - (1. - Qa/Q0)/dfdy

            HRadius = (self.B*y1/(self.B + 2*y1))**(2./3)
            Q1 = self.alpha*HRadius*self.B*y1
            Err = abs(Q1 - Q0)

            y0 = y1
            Q0 = Q1

        return y0

    def run(self):
        max_rain = max(self.rainData[:, 1])
        self.max_Q = 2*max_rain*self.Lc*self.Ls

        maxTime = max(self.rainData[:, 0])
        minTime = min(self.rainData[:, 0])
        self.q_out = list()

        flowrate0 = np.zeros(self.grds)

        self.currentSitu = list()
        for i in range(0, self.grds):
            self.currentSitu.append([[self.x[i], 10**-10]])
        self.currentSitu[0][0][1] = 0.0

        ct = minTime  # Current Time
        while ct < maxTime:
            flowrate0, ct = self.depth_iterate1(flowrate0.copy(), ct)
            print ct, flowrate0[-1]
            self.q_out.append([ct, flowrate0[-1]])

    def depth_iterate1(self, flowRate0, ct):
        dt = self.dt
        N = self.grds

        print 'current time: %f' % ct

        Courant_OK = False

        while Courant_OK is not True:
            Courant_OK = self.courant_Check(flowRate0, ct, dt)
            if not Courant_OK:
                dt = 0.8*dt

        flowRate1 = np.zeros(len(flowRate0))
        for j in range(1, N):
            flowRate1[j], xl = self.depth_iterate2(flowRate0, j, ct, dt)
        self.movingWave(ct, dt, flowRate1)

        ct = ct + dt
        return flowRate1, ct

    def depth_iterate2(self, flowRate1, i, ct, dt):
        try:
            Qp, xl1 = self.newton(flowRate1, i, ct, dt)
        except(ValueError):
            Qp, xl1 = self.bisection(flowRate1, i, ct, dt)

        return Qp, xl1

    def movingWave(self, ct, dt, depth1):
        currentSitu = self.currentSitu
        alpha = self.alpha
        m = 3.0/2
        for j in range(0, len(currentSitu)-1):
            outRegion = list()
            if currentSitu[j]:
                for k in range(0, len(currentSitu[j])):
                    ql = (self.get_catInflow(ct+dt) + self.get_catInflow(ct))
                    Ql = currentSitu[j][k][1]
                    xl = currentSitu[j][k][0]
                    y0 = self.getH(Ql)
                    A0 = self.B*y0

                    if ql != 0.:
                        Ap = A0 + ql*dt
                        yp = Ap/self.B
                        xl = xl + alpha/ql*self.B*(yp**m - y0**m)
                    else:
                        xl = xl + alpha*m*y0**(m-1)*dt

                    cA = A0 + ql*dt
                    currentSitu[j][k][0] = xl
                    currentSitu[j][k][1] = self.fromAtoQ(cA)

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

    def getDADX(self, x):
        situArray = self.situToArray()

        for i in range(1, len(situArray)):
            if x > situArray[i-1, 0] and x <= situArray[i, 0]:
                yp = self.getH(situArray[i, 1])
                yl = self.getH(situArray[i-1, 1])

                return (yp - yl)/(situArray[i, 0] - situArray[i-1, 0])
            elif x == situArray[i-1, 0]:
                yp = self.getH(situArray[i-1, 1])
                yl = self.getH(situArray[i-2, 1])

                return (yp - yl)/(situArray[i, 0] - situArray[i-2, 0])

    def getFlowRate(self, x):
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

    def fromAtoQ(self, Ap):
        yp = Ap/self.B
        HRadius = (Ap/(self.B + 2*yp))
        Qp = self.alpha*HRadius**(2./3)*Ap

        return Qp

    def newton(self, flowRate1, i, ct, dt):
        # print 'Use Newton Method'
        dx = self.dx
        m = 3.0/2
        alpha = self.alpha
        Ql0 = flowRate1[i-1] + 0.5*(flowRate1[i-1] - flowRate1[i])
        xl0 = 0.5*dx
        Err = 1.0

        counter = 0
        ql = (self.get_catInflow(ct+dt) + self.get_catInflow(ct))
        while Err > 10**-5:
            yl0 = self.getH(Ql0)
            Al0 = self.B*yl0

            Ap = Al0 + ql*dt
            yp = Ap/self.B

            try:
                if ql != 0.:
                    f = dx - xl0 - alpha/ql*self.B*(yp**m - yl0**m)
                    dfdx = -1 - alpha/ql*self.B*m*yp**(m-1)*self.getDADX(
                        self.x[i-1] + xl0)
                elif ql == 0.:
                    f = dx - xl0 - alpha*m*yl0**(m-1)*dt
                    dfdx = -1 - alpha*m*(m-1)*yl0**(m-2)*dt*self.getDADX(
                        self.x[i-1] + xl0)
            except(TypeError, ValueError):
                raise ValueError

            xl1 = xl0 - f/dfdx
            Ql1 = self.getFlowRate(self.x[i-1] + xl1)

            Err = abs(xl1 - xl0)/xl0
            xl0 = xl1
            Ql0 = Ql1

            if np.isnan(xl1):
                raise ValueError

            counter += 1
            if counter > 50:
                raise ValueError

        yp = self.getH(Ql1)
        Ap = self.B*yp + ql*dt
        Qp = self.fromAtoQ(Ap)

        return Qp, xl1

    def bisection(self, flowRate1, i, ct, dt):
        # print 'Use Bisection Method'
        dx = self.dx
        Err = 1.0
        x1 = dx
        x2 = 0.5*dx
        x3 = 0.0

        while Err > 10**-5:
            f1 = self.function(flowRate1, i, x1, ct, dt)
            f2 = self.function(flowRate1, i, x2, ct, dt)
            f3 = self.function(flowRate1, i, x3, ct, dt)
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

        ql = (self.get_catInflow(ct+dt) + self.get_catInflow(ct))
        Q0 = self.getFlowRate(self.x[i-1] + x2)
        y0 = self.getH(Q0)
        Ap = self.B*y0 + ql*dt
        Qp = self.fromAtoQ(Ap)

        return Qp, x2

    def function(self, flowrate1, i, xl, ct, dt):
        dx = self.dx
        alpha = self.alpha
        m = 3.0/2
        Ql = self.getFlowRate(self.x[i-1] + xl)
        yl = self.getH(Ql)
        Al = self.B*yl

        ql = (self.get_catInflow(ct+dt) + self.get_catInflow(ct))
        Ap = Al + ql*dt
        yp = Ap/self.B
        if ql != 0.:
            f = dx - xl - alpha/ql*self.B*(yp**m - yl**m)
        elif ql == 0:
            f = dx - xl - alpha*m*yl**(m-1)*dt

        return f
