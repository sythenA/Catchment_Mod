
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
        dt = t/50.
        ct = 0.
        intVal = 0.0

        for i in range(1, 50):
            ct = i*dt
            if i % 2 == 0:
                intVal += 2*ms*self.t4DepInt(ct)**(ms-1)
            else:
                intVal += 3*ms*self.t4DepInt(ct)**(ms-1)
        intVal += ms*self.t4DepInt(t)**(ms-1)

        return dt*3*intVal/8.

    def t4DepInt(self, t):
        dt = t/50.
        intVal = self.getSideInflow(0)
        for i in range(1, 50):
            ct = i*dt
            if i % 2 != 0:
                intVal += 3*self.getSideInflow(ct)
            else:
                intVal += 2*self.getSideInflow(ct)
        intVal += self.getSideInflow(t)

        return dt*3*intVal/8.

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
        dt = t/50.0
        val = ((x0/self.k)**(1./mc))**(mc-1)

        for i in range(1, 50):
            ct = dt*i
            if i % 2 == 0:
                val += 2*((x0/self.k)**(1./mc) + ct)**(mc-1)
            else:
                val += 3*((x0/self.k)**(1./mc) + ct)**(mc-1)
        val += ((x0/self.k)**(1./mc) + t)**(mc-1)

        return dt*3*val/8.

    def getSideInflow(self, t):
        mc = self.mc

        def function(x0, t):
            f = 1.0 - x0 - self.mc*self.sideInt(x0, t)
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

    def phase1Int(self, xs, t):
        ms = self.ms
        dt = t/50.
        val = ((xs/self.k)**(1/ms) + self.Lambda**-1*self.t4DepInt(t))**(ms-1)

        for i in range(1, 50):
            ct = i*dt
            if i % 2 == 0:
                val += 2*((xs/self.k)**(1/ms) +
                          self.Lambda**-1*self.t4DepInt(ct))**(ms-1)
            else:
                val += 3*((xs/self.k)**(1/ms) +
                          self.Lambda**-1*self.t4DepInt(ct))**(ms-1)
        val += ((xs/self.k)**(1/ms) + self.Lambda**-1*self.t4DepInt(t))**(ms-1)

        return dt*3*val/8.

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
                break
            Err = abs(Ulim - Llim)
            # print Ulim, Llim

        xs = 0.5*(Ulim + Llim)

        A = (xs/self.k)**(1/ms) + self.Lambda**-1*self.t4DepInt(t)

        return A**ms, xs

    def phase2DepInt(self, t1, t):
        dt = (t - t1)/50.0
        val = self.getSideInflow(t1)

        for i in range(1, 50):
            ct = t1 + dt*i
            if i % 2 == 0:
                val += 2*self.getSideInflow(ct)
            else:
                val += 3*self.getSideInflow(ct)
        val += self.getSideInflow(t)

        # print t, t1, val*dt/3./self.Lambda
        return val*3*dt/8.

    def phase2XInt(self, t1, t):
        dt = (t - t1)/50.0
        ms = self.ms
        val = 0

        for i in range(1, 50):
            ct = t1 + dt*i
            if i % 2 == 0:
                val += 2*self.phase2DepInt(t1, ct)**(ms-1)
            else:
                val += 3*self.phase2DepInt(t1, ct)**(ms-1)
        val += self.phase2DepInt(t1, t)**(ms-1)

        return val*3*dt/8.

    def phase2(self, t, init_t1):
        Ulim = t
        Llim = init_t1

        Lambda = self.Lambda
        ms = self.ms
        Err = abs(Ulim - Llim)

        while Err > 1.0E-6:
            f = Lambda**1.5 - ms*self.phase2XInt(0.5*(Ulim + Llim), t)
            if f > 0:
                Ulim = 0.5*(Ulim + Llim)
            elif f < 0:
                Llim = 0.5*(Ulim + Llim)
            elif f == 0:
                break
            Err = abs(Ulim - Llim)

        t1 = 0.5*(Ulim + Llim)
        A = Lambda**-1*self.phase2DepInt(t1, t)

        return A**ms, t1

    def phase3XInt(self, t, t0):
        dt = (t - 1.0)/50.0
        ms = self.ms

        val = (self.phase2DepInt(t0, 1.0))**(ms-1)
        for i in range(1, 50):
            ct = 1.0 + dt*i
            if i % 2 == 0:
                val += 2*(self.phase2DepInt(t0, 1.0) + (ct-1.0))**(ms-1)
            else:
                val += 3*(self.phase2DepInt(t0, 1.0) + (ct-1.0))**(ms-1)
        val += (self.phase2DepInt(t0, 1.0) + (t-1.0))**(ms-1)

        return val*3*dt/8.

    def phase3(self, t, init_L):
        Ulim = 1.0
        Llim = init_L
        Lambda = self.Lambda
        ms = self.ms

        Err = abs(Ulim - Llim)
        while Err > 1.0E-6:
            f = Lambda**1.5 - ms*(self.phase3XInt(t, 0.5*(Ulim + Llim)) +
                                  self.phase2XInt(0.5*(Ulim + Llim), 1.0))
            if f > 0:
                Ulim = 0.5*(Ulim + Llim)
            elif f < 0:
                Llim = 0.5*(Ulim + Llim)
            elif f == 0:
                break
            Err = abs(Ulim - Llim)
            print f, Ulim, Llim

        t0 = 0.5*(Ulim + Llim)
        A = Lambda**-1*(self.phase2DepInt(t0, 1.0) + (t - 1.0))

        return A**ms, t0

    def run(self):
        tmax = self.getMaxTime()*self.tc
        self.t4 = self.getT4()

        dt = self.dt
        t = dt
        outQ = list()
        outQ.append([0.0, self.Q_max/self.k])

        xs = 1.0
        t1 = 0.0
        t0 = self.t4
        while t <= tmax:
            if t < self.t4*self.tc:
                ct = t/self.tc
                Q, xs = self.phase1(ct, xs)
                outQ.append([t, Q])
                print "current time = %f, flowrate = %f" % (t, Q)
            elif t >= self.t4*self.tc and t <= self.tc:
                ct = t/self.tc
                Q, t1 = self.phase2(ct, t1)
                outQ.append([t, Q])
                print "current time = %f, flowrate = %f" % (t, Q)
            elif t > self.tc:
                ct = t/self.tc
                Q, t0 = self.phase3(ct, t0)
                outQ.append([t, Q])
                print "\n"
                print "current time = %f, flowrate = %f" % (t, Q)
                print "\n"
            t += dt

        outQ = np.array(outQ)
        outQ[:, 1] = outQ[:, 1]*self.Q_max
        self.outQ = outQ
