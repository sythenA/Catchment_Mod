
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

    def t4DepInt(self):
        Ulim = 1.0
        Llim = 0.0

    def getSideInflow(self, t):
        mc = self.mc
        def depthInt(x0, t):
            mc = self.mc
            dt = t/200.0
            box = np.zeros(201)

            for i in range(0, 201):
                ct = dt*i
                box[i] = ((x0/self.k)**(1./mc) + ct)**(mc-1)

            for j in range(1, 201):
                if j % 2 == 0:
                    box[j] = 2*box[j]
                elif j % 2 != 0:
                    box[j] = 4*box[j]

            return dt*sum(box)/3.

        def function(x0, t):
            f = 1.0 - x0 - self.mc*depthInt(x0, t)
            return f

        Ulim = 1.0
        Llim = 0.0
        Err = 1.0

        while Err > 1.0E-8:
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

    def run(self):
        t = 0
        tmax = self.getMaxTime()*self.tc

        dt = self.dt
        outQ = list()
        outQ.append([0.0, self.Q_max/self.k])
