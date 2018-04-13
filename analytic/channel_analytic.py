
import numpy as np
from math import sqrt


class ana_channel:
    def __init__(self, max_rain, to, tmax, dt, Lc, Ls, nc, ns, Sc, Ss, B):
        self.max_rain = max_rain/3600.0/1000.0
        self.dt = dt
        self.B = B
        alpha_s = 1/ns*sqrt(Ss)
        alpha_c = 1/nc*sqrt(Sc)
        self.alpha_s = alpha_s  # stream
        self.alpha_c = alpha_c  # catchment
        H_max = self.getDepth(2*self.max_rain*Lc*Ls)
        self.A_max = self.B*H_max  # Maximum cross-section in channel
        self.t_max = tmax
        self.Lc = Lc
        self.Ls = Ls
        self.qc = self.max_rain*Lc  # Catchment inflow rate at steady-state
        self.tc = (self.Lc*self.max_rain**(-2./3)/self.alpha_c)**(3./5)
        # tc: Time to steady-state on the slopes
        self.to = to  # rain stops
        # ------------------------------------------------------------------- #
        upper = self.A_max
        lower = 2*self.max_rain*Lc*self.tc
        self.Lambda = upper/lower
        # ------------------------------------------------------------------- #

    def getDepth(self, Q):
        y0 = (Q/self.alpha_s/self.B)**(2./5)
        HRadius = self.B*y0/(self.B + 2*y0)
        Q0 = self.alpha_s*HRadius**(2./3)*self.B*y0

        Err = abs(Q-Q0)
        while Err > 1.0E-8:
            y1 = y0 - (1. - Q/Q0)/((5*self.B + 6*y0)/3/y0/(self.B + 2*y0))
            HRadius = self.B*y1/(self.B + 2*y1)
            Q1 = self.alpha_s*HRadius**(2./3)*self.B*y1
            Err = abs(Q1 - Q0)

            y0 = y1
            Q0 = Q1

        return y0

    def phase1Int(self, t1, t0):
        dt = (t1 - t0)/200.
        intPart = np.zeros(201)

        ct = t0 + dt
        for i in range(1, 200):
            if i % 2 == 0:
                intPart[i] = 2*(1./3*(ct**3 - t0**3))**0.5
            elif i % 2 != 0:
                intPart[i] = 4*(1./3*(ct**3 - t0**3))**0.5
            ct = ct + dt
        intPart[-1] = (1./3*(t1**3 - t0**3))**0.5

        return sum(intPart)*dt/3

    def phase1(self, ct, t1):
        t1 = t1/self.tc
        t = ct/self.tc

        t0 = self.phase1Bisec(ct)
        Area = self.Lambda**-1*(1./3*(t**3 - t0**3))

        return Area

    def phase1Bisec(self, ct):
        def function(t1, t0):
            f = 1. - self.Lambda**-1.5*1.5*self.phase1Int(t1, t0)
            return f
        rct = ct/self.tc
        U_init = rct
        D_init = 0.0

        Err = abs(function(rct, U_init) - function(rct, D_init))
        while Err > 1.0E-6:
            f = function(rct, 0.5*(U_init + D_init))

            if f > 0:
                U_init = 0.5*(U_init + D_init)
            elif f < 0:
                D_init = 0.5*(U_init + D_init)
            elif f == 0:
                return 0.5*(U_init + D_init)

            Err = abs(function(rct, U_init) - function(rct, D_init))

        return 0.5*(U_init + D_init)

    def phase2(self, ct):
        rct = ct/self.tc

        t0 = self.phase2Bisec(rct)
        Area = self.Lambda**-1*(1./3*(1.0 - t0**3) + (rct - 1.0))

        return Area

    def phase2Bisec(self, ct):
        def function(t, t0):
            f = (self.Lambda**1.5 - ((t - 1.0) + 1./3*(1.0 - t0**3))**1.5 +
                 (1./3*(1. - t0**3))**1.5 -
                 1.5*self.phase1Int(1.0, t0))
            return f

        U_init = 1.0
        D_init = 0.0

        Err = abs(function(ct, U_init) - function(ct, D_init))
        while Err > 1.0E-8:
            fun = function(ct, 0.5*(U_init + D_init))
            if fun > 0.:
                U_init = 0.5*(U_init + D_init)
            elif fun < 0.:
                D_init = 0.5*(U_init + D_init)
            elif fun == 0:
                return 0.5*(U_init + D_init)

            Err = abs(function(ct, U_init) - function(ct, D_init))

        return 0.5*(U_init + D_init)

    def integrateG(self, dotT):
        intPart = np.zeros(101)
        dn = dotT/100.

        cn = dn
        for i in range(1, 100):
            if i % 2 != 0:
                intPart[i] = 4.*(2./3*cn**3 + cn - 2./3*(cn**2 + 1)**1.5)
            elif i % 2 == 0:
                intPart[i] = 2.*(2./3*cn**3 + cn - 2./3*(cn**2 + 1)**1.5)
            cn += dn
        intPart[-1] = 2./3*dotT**3 + dotT - 2./3*(dotT**2 + 1)**1.5

        return sum(intPart)*dn/3

    def phase3(self, ct):
        def xToA(x):
            q = 2*self.max_rain*x*self.Ls*self.Lc
            A = self.fromQtoA(q)
            return A/self.A_max

        def fun(u, t0):
            part1 = 2./3*(1.0E-8**2+1)**1.5 - 2./3*1.0E-8**3 - 1.0E-8
            part2 = 2./3*(u**2 + 1)**1.5 - 2./3*u**3 - u
            return self.Lambda**-1*(t0 + part1 - part2)

        ct = (ct - self.to)/self.tc + 1.0
        t0 = self.phase3Bisec(ct-1.0)
        Area = fun(ct-1.0, t0)

        return Area

    def phase3Bisec(self, t):
        def function(t, t0):
            fun = (self.Lambda**1.5 - 1.5*self.phase3Int(t, t0) -
                   t0**1.5)
            return fun

        U_init = 1.0
        D_init = 0.0

        f = function(t, 0.5*(U_init + D_init))
        # print 'f=%f' %f
        Err = abs(function(t, U_init) - function(t, D_init))

        while Err > 1.0E-8:
            f = function(t, 0.5*(U_init + D_init))
            if f > 0:
                D_init = 0.5*(U_init + D_init)
            elif f < 0:
                U_init = 0.5*(U_init + D_init)
            elif f == 0:
                return 0.5*(U_init + D_init)

            Err = abs(function(t, U_init) - function(t, D_init))
            # print t, U_init, D_init, function(t, 0.5*(U_init+D_init))

        return 0.5*(U_init + D_init)

    def phase3Int(self, u, t0):
        def xToA(x):
            q = 2*self.max_rain*x*self.Ls*self.Lc
            A = self.fromQtoA(q)
            return A/self.A_max

        def function(u, t0):
            try:
                part1 = 2./3*(1.0E-8**2 + 1)**1.5 - 2./3*1.0E-8**3 - 1.0E-8
                part2 = 2./3*(u**2 + 1)**1.5 - 2./3*u**3 - u
                return (t0 + part1 - part2)**0.5
            except(ValueError):
                print u, part1, part2, A0*self.Lambda
                input()

        dt = u/200.
        intPart = np.zeros(201)

        ct = dt
        for i in range(0, 200):
            if i % 2 == 0:
                intPart[i] = 2*function(ct, t0)
            elif i % 2 != 0:
                intPart[i] = 4*function(ct, t0)
            ct = ct + dt
        intPart[-1] = function(u, t0)

        return sum(intPart)*dt/3

    def phase4(self, ct, t3):
        def function(u, u0):
            part1 = 2./3*(u0**2 + 1)**1.5 - 2./3*u0**3 - u0
            part2 = 2./3*(u**2 + 1)**1.5 - 2./3*u**3 - u
            return (part1 - part2)*self.Lambda**-1

        ct = (ct - self.to)/self.tc
        t3 = (t3 - self.to)/self.tc
        u0 = self.phase4Bisec(ct, 0.)

        Area = function(ct, u0)
        return Area

    def phase4Bisec(self, u, u_init):
        def function(u, u0):
            return self.Lambda**1.5 - 1.5*self.phase4Int(u, u0)

        U_init = u
        D_init = u_init

        Err = abs(function(u, U_init) - function(u, D_init))

        while Err > 1.0E-8:
            f = function(u, 0.5*(U_init + D_init))
            if f < 0:
                D_init = 0.5*(U_init + D_init)
            elif f > 0:
                U_init = 0.5*(U_init + D_init)
            elif f == 0:
                return 0.5*(U_init + D_init)
            Err = abs(function(u, U_init) - function(u, D_init))

        return 0.5*(U_init + D_init)

    def phase4Int(self, u, u0):
        def function(u, u0):
            try:
                part1 = 2./3*(u0**2 + 1)**1.5 - 2./3*u0**3 - u0
                part2 = 2./3*(u**2 + 1)**1.5 - 2./3*u**3 - u
                return sqrt(part1 - part2)
            except(ValueError):
                print part1, part2, u, u0
                input()

        du = (u - u0)/200
        intPart = np.zeros(201)

        cu = u0 + du
        for i in range(1, 200):
            if i % 2 == 0:
                intPart[i] = 2*function(cu, u0)
            elif i % 2 != 0:
                intPart[i] = 4*function(cu, u0)
            cu = cu + du
        intPart[-1] = function(u, u0)

        return sum(intPart)*du/3.

    def fromQtoA(self, Qa):
        H0 = (Qa/self.B/self.alpha_s)**0.6
        hradius = (self.B*H0/(self.B + 2*H0))
        Q0 = self.alpha_s*hradius**(2./3)*self.B*H0

        Err = abs(Qa - Q0)
        while Err > 1.0E-8:
            H1 = H0 - (1. - Qa/Q0)/(5*self.B + 6*H0)*3*H0*(self.B + 2*H0)
            hradius = (self.B*H1/(self.B + 2*H1))
            Q1 = self.alpha_s*hradius**(2./3)*self.B*H1

            H0 = H1
            Q0 = Q1
            Err = abs(Qa - Q1)
        return H0

    def getQ(self, A):
        # Input Cross-section Area, return flowrate
        H = A/self.B
        HRadius = self.B*H/(self.B + 2.*H)

        return self.alpha_s*HRadius**(2./3)*self.B*H

    def steadyTime(self):
        return self.Lambda + 1.0

    def secInt(self, t):
        def function(u):
            part1 = 2./3*(1.0E-8**2 + 1)**1.5 - 2./3*1.0E-8**3 - 1.0E-8
            part2 = 2./3*(u**2 + 1)**1.5 - 2./3*u**3 - u
            return sqrt(part1 - part2)

        dt = (t-1.0)/200
        intPart = np.zeros(201)

        u = dt
        for i in range(1, 200):
            if i % 2 == 0:
                intPart[i] = 2*function(u)
            elif i % 2 != 0:
                intPart[i] = 4*function(u)
            u = u + dt
        intPart[-1] = function(t - 1.0)

        return sum(intPart)*dt/3

    def secondReach(self):
        def function(t):
            f = self.Lambda**1.5 - 1.5*self.secInt(t)
            return f

        U_init = 2.0
        D_init = 1.001

        Err = abs(function(U_init) - function(D_init))
        while Err > 1.0E-8:
            f = function(0.5*(U_init + D_init))
            if f < 0:
                U_init = 0.5*(U_init + D_init)
            elif f > 0:
                D_init = 0.5*(U_init + D_init)
            elif f == 0:
                return 0.5*(U_init + D_init)

            Err = abs(function(U_init) - function(D_init))

        return 0.5*(U_init + D_init)

    def run(self):
        dt = self.dt
        tmax = self.t_max
        x = np.arange(0., 1.001, 0.01)
        self.x = x
        ct = 0.0 + self.dt  # Current Time
        tc = self.tc  # Time of steady on two slopes
        Amax = self.A_max

        characteristics = list()  # Characteristic lines
        out_Q = list()
        out_Q.append([0., 0.])

        Area = np.zeros(len(x))

        # Initial Characteristics
        for i in range(0, len(x)):
            characteristics.append([x[i], 0.0, 0.0, x[i]])

        t1 = (self.Lambda**1.5*5/sqrt(3))**0.4*self.tc
        # first wave start at x=0 reaches 1.0
        t2 = self.steadyTime()*self.tc
        # Time to reach steady-state under steady rainfall
        t3 = self.secondReach()
        t3 = (t3 - 1.0)*self.tc + self.to

        while ct < tmax:
            if ct <= t1:
                Area = self.Lambda**-1*1./3*(ct/self.tc)**3
                Area = Area*Amax
                out_Q.append([ct, self.getQ(Area)])

            if ct > t1 and ct <= tc:
                Area = self.phase1(ct, t1)*Amax
                out_Q.append([ct, self.getQ(Area)])

            elif ct > tc and ct <= t2:
                Area = self.phase2(ct)
                Area = Area*Amax
                out_Q.append([ct, self.getQ(Area)])

            elif ct > t2 and ct <= self.to:
                out_Q.append([ct, self.getQ(self.A_max)])

            elif ct > self.to and ct <= t3:
                Area = self.phase3(ct)
                Area = Area*self.A_max
                out_Q.append([ct, self.getQ(Area)])
            elif ct > t3:
                Area = self.phase4(ct, t3)
                Area = Area*self.A_max
                out_Q.append([ct, self.getQ(Area)])

            print ct, self.getQ(Area)
            ct = ct + dt  # Add dt to current time

        self.outQ = np.array(out_Q)
