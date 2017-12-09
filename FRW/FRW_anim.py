from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Some unit conversion factors:
H0 = 68              #Hubble constant H = 70 km/s/Mpc
Mpc = 3.085677581e19   # km
km =  1.0
Gyr = 3.1536e16        #second

H_0 = (H0 * km * Gyr) /  Mpc    #Here H_0 is in the units of 1/Gyr

# If the universe had expanded with a constant rate:
print 1/H_0
# Notice that 13.9 Gyr is a good approximation to known age (13.8 Gyr)


# Friedmann Equation:

def Friedmann(a,m,r,l):
    Omega_0 = m + r + l
    t = (1/ np.sqrt( (r/(a*a)) + (m/a) + (l * (a*a)) + (1-Omega_0)  ) ) /H_0
    #t_0 = (1/( (m/a) + (lam * (a*a)))**(0.5))/gyr # Matter - Lambda Flat Universe Approximation

    return t

# Trapezoid Integration Method:
def Trapezoidal(a,b,m,r,l):
    n = 1000                 # Step number
    deltaX = (b-a)/n         # Step size

    AGE = 0

    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)

    for i in range(n):
        x[i] = a + i*deltaX
        y[i] = Friedmann(x[i],m,r,l)
        z[i] = (deltaX/2) * (2*np.sum(y) - y[0] - y[n-1])

        if (x[i] == 1 or 1 <= x[i] <= 1+deltaX):
            AGE = z[i]

    #result = (deltaX/2) * (2*np.sum(y) - y[0] - y[n-1])

    print 'Age of the universe with m = %5.3f, r = %5.3f, lambda = %5.3f is %5.3f Gyr' %(m,r,l,AGE)

    return x,z,AGE

plt.rc('text', usetex=True)  #In order to use LaTeX
plt.rc('font', family='serif') #In order to use Serif (mathced font with LaTeX)

a,t,age = Trapezoidal(1e-10,10,0.3,0,0.68)
a1,t1,age1 = Trapezoidal(1e-10,10,1,0,-0.3)
a2,t2,age2 = Trapezoidal(1e-10,10,0.6,0.1,0.3)

# Visualize:
def update(num, t, a,t1,a1,t2,a2, line1,line2,line3):
    line1.set_data(t[:num], a[:num])
    line1.axes.axis([0, 80, 0, 5])
    line1.axes.set_xlabel = '$t$ (Gyr)'
    line1.axes.set_ylabel = '$a(t)$'
    line2.set_data(t1[:num], a1[:num])
    line3.set_data(t2[:num], a2[:num])
    return line1,line2,line3

fig = plt.figure(figsize=(8, 6), dpi=150, facecolor='white', edgecolor='w')
ax1 = fig.add_subplot(111)

line1, = ax1.plot(t, a, color='blue', label = '$\Lambda$ - CDM')
line2, = ax1.plot(t1, a1, color='orange', label = '$\Lambda$ Collapse')
line3, = ax1.plot(t2, a2, color='red', label = '$\Omega_m = 0.6, \Omega_\Lambda = 0.3$')

ax1.set_xlabel('$t$ (Gyr)')
ax1.set_ylabel('$a(t)$')
ax1.legend()

ani = animation.FuncAnimation(fig, update, len(a), fargs=[t, a, t1,a1,t2,a2,line1,line2,line3],
                              interval=15, blit=True)

ani.save('FRW1.mp4',extra_args=['-vcodec', 'libx264'])

# Marker on the current ages
ax1.scatter(age,1,color='blue')
ax1.scatter(age1,1,color='orange')
ax1.scatter(age2,1,color='red')
plt.show()
