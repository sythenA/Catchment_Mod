import numpy as np
from math import sqrt


def rising_phase(slope, n, length, intensity, duration):
    intensity = intensity/3600.0/1000.0
    ana_out = list()
    alpha = 1./n*sqrt(slope)
    m = 5.0/3
    tc = (intensity**(1-m)*length/alpha)**(1/m)
    max_flow = alpha*(intensity*tc)**m
    t = 0.0
    while t <= tc:
        ana_out.append([t, alpha*(intensity*float(t))**m])
        t = t + 1.0

        print t

    while t < duration:
        ana_out.append([t, max_flow])
        t = t + 1.0

        print t

    return np.array(ana_out)


def exit_Depth(slope, n, length, intensity, t):
    alpha = 1./n*sqrt(slope)
    xl0 = length - 1.0
    Err = 1.0
    m = 5.0/3

    while Err > 10**-3:
        f = length - xl0 - alpha*m*(intensity*xl0/alpha)**((m-1)/m)*t
        dfdx = -1 - (m-1)*intensity*(intensity*xl0/alpha)**(-1/m)*t

        xl1 = xl0 - f/dfdx

        Err = abs(xl1 - xl0)/xl0

        xl0 = xl1

        print f, dfdx

    return xl1*intensity


def falling_phase(slope, n, length, intensity, duration, max_time):
    t = duration
    ana_out = list()
    intensity = intensity/3600.0/1000

    while t < max_time:
        t = t + 1.0
        flowrate = exit_Depth(slope, n, length, intensity, t-duration)
        ana_out.append([t, flowrate])

        print t

    return np.array(ana_out)
