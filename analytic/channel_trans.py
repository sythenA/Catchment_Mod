
import numpy as np
from math import sqrt


class c_transition:
    def __init__(self, i0, i1, Lc, Ls, nc, ns, Sc, Ss, dt, B):
        i0 = i0/1000.0/3600.0
        i1 = i1/1000.0/3600.0
        self.k = i1/i0
        self.dt = dt
        self.q_max = i1*Lc
        self.Q_max = 2*i1*Lc*Ls

        self.B = B
        self.alpha_s = 1./ns*sqrt(Ss)
        self.A_max = self.fromQtoA(self.Q_max)*B
        self.Lc = Lc
        self.Ls = Ls
        self.mc = 2.0
        self.ms = 1.5

        alpha_c = 1./nc*sqrt(Sc)
        self.tc = (i1*Lc/alpha_c)**0.5/i1
        self.Lambda = self.A_max/self.q_max/self.tc

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

    def getMaxTime(self):
        return 1.0 + self.Lambda

    def t4XInt(self, t):
        ms = self.ms
        dt = t/20.
        ct = 0.
        intVal = 0.0

        for i in range(1, 20):
            ct = i*dt
            if i % 2 == 0:
                intVal += 2*ms*self.t4DepInt(ct)**(ms-1)
            else:
                intVal += 4*ms*self.t4DepInt(ct)**(ms-1)

        return dt*intVal/3.

    def t4DepInt(self, t):
        dt = t/20.
        intVal = self.getSideInflow(0)
        for i in range(1, 20):
            ct = i*dt
            if i % 2 != 0:
                intVal += 4*self.getSideInflow(ct)
            else:
                intVal += 2*self.getSideInflow(ct)

        return dt*intVal/3.

    def getT4(self):
        print "getting t4:"

        def function(t):
            ms = self.ms
            Lambda = self.Lambda
            f = Lambda**ms - self.t4XInt(t)
            return f

        Ulim = 1.0
        Llim = 0.0
        Err = abs(Ulim - Llim)
        while Err > 1.0E-6:
            f = function(0.5*(Ulim + Llim))
            if f < 0:
                Ulim = 0.5*(Ulim + Llim)
            elif f > 0:
                Llim = 0.5*(Ulim + Llim)
            elif f == 0:
                return 0.5*(Ulim + Llim)
            Err = abs(Ulim - Llim)
            print Ulim, Llim

        val = (Ulim + Llim)*0.5
        print "t4 = %f" % val
        return val

    def sideInt(self, x0, t):
        mc = self.mc
        dt = t/20.0
        val = ((x0/self.k)**(1./mc))**(mc-1)

        for i in range(1, 21):
            ct = dt*i
            if i % 2 == 0:
                val += 2*((x0/self.k)**(1./mc) + ct)**(mc-1)
            else:
                val += 4*((x0/self.k)**(1./mc) + ct)**(mc-1)

        return dt*val/3.

    def getSideInflow(self, t):
        mc = self.mc

        def function(x0, t):
            f = 1.0 - x0 - self.mc*self.sideInt(x0, t)
            return f

        Ulim = 1.0
        Llim = 0.0
        Err = 1.0

        while Err > 1.0E-6:
            x0 = 0.5*(Ulim + Llim)
            f = function(x0, t)
            if f < 0:
                Ulim = 0.5*(Ulim + Llim)
            elif f > 0:
                Llim = 0.5*(Ulim + Llim)
            elif f == 0:
                break

            Err = abs(Ulim - Llim)

        return ((x0/self.k)**(1./mc) + t)**mc

    def phase1Int(self, xs, t):
        ms = self.ms
        dt = t/20.
        val = ((xs/self.k)**(1/ms) + self.Lambda**-1*self.t4DepInt(t))**(ms-1)

        for i in range(1, 20):
            ct = i*dt
            if i % 2 == 0:
                val += 2*((xs/self.k)**(1/ms) +
                          self.Lambda**-1*self.t4DepInt(ct))**(ms-1)
            else:
                val += 4*((xs/self.k)**(1/ms) +
                          self.Lambda**-1*self.t4DepInt(ct))**(ms-1)

        return dt*val/3.

    def phase1(self, t, init_U):
        def function(xs, t):
            ms = self.ms
            Lambda = self.Lambda
            f = Lambda - Lambda*xs - ms*self.phase1Int(xs, t)
            return f

        ms = self.ms
        Ulim = init_U
        if init_U > 0.6:
            Llim = 0.5
        else:
            Llim = 0.0
        Err = abs(Ulim - Llim)

        while Err > 1.0E-4:
            f = function(0.5*(Ulim + Llim), t)
            if f > 0:
                Llim = 0.5*(Ulim + Llim)
            elif f < 0:
                Ulim = 0.5*(Ulim + Llim)
            elif f == 0:
                return 0.5*(Ulim + Llim)
            Err = abs(Ulim - Llim)
            print Ulim, Llim

        xs = 0.5*(Ulim + Llim)

        A = (xs/self.k)**(1/ms) + self.Lambda**-1*self.t4DepInt(t)

        return A**ms, xs

    def run(self):
        tmax = self.getMaxTime()*self.tc
        self.t4 = self.getT4()

        dt = self.dt
        t = dt
        outQ = list()
        outQ.append([0.0, self.Q_max/self.k])

        xs = 1.0
        while t <= tmax:
            if t <= self.t4*self.tc:
                ct = t/self.tc
                Q, xs = self.phase1(ct, xs)
                outQ.append([t, Q])
                print t, Q
            t += dt
