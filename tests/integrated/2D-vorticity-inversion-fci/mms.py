from __future__ import division
from __future__ import print_function

from boutdata.mms import Metric, sin, cos, Div_par, Grad_par, exprToStr, diff, y, t, DDY, x, z, DDX, DDZ, Delp2, D2DX2
from math import pi

# Length of the y domain
Ly = 2.0 * pi

# Atomic mass number
AA = 1.0
B = 1.0 + 0.0 * x
# metric tensor
metric = Metric()  # Identity

qe = 1.60217663e-19
Me = 9.1093837e-31
e0 = 8.85418781e-12
Pi = pi
q = 1.0
Tnorm = 5
Nnorm = 1e18
Bnorm = 1.0
Omega_ci = qe * Bnorm / (1836.0*Me);
rho_s = 0.00022847

# Define solution in terms of input x,y,z
omega = 0.0001
#n = 1 + 0.1*sin(2 * pi * x) * sin(3 * z)# * sin(t*omega)
phi =  50.0 * sin(2.0 * pi * x) * sin(2 * z) / Tnorm
n = 1.0 + 0.25 * sin(4.0 * pi * x) * sin(4 * z + 2.31312)
Ti = (20.0 + 5.25 * sin(2.0 * pi * x) * sin(2 * z + 0.3132131)) / Tnorm


replace = [(x, metric.x), (z, metric.z * 2.0 * pi ) ]
phi = phi.subs(replace)
n = n.subs(replace)
Ti = Ti.subs(replace)
B = B.subs(replace)

pre_phi = AA / B**2
pre_Pi = AA / (q * B**2)
Pi = n * Ti
Vort = (DDX(pre_phi * DDX(phi)) + DDZ(pre_phi * DDZ(phi)) + DDX(pre_Pi * DDX(Pi)) + DDZ(pre_Pi * DDZ(Pi)))*(rho_s**2)

# Substitute back to get input y coordinates
replace = [(metric.x, x), (metric.z, z / (2.0 * pi) ) ]
Vort = Vort.subs(replace)
Pi = Pi.subs(replace)
n = n.subs(replace)
print("[Vort]")
print("function = " + exprToStr(Vort))
print("bndry_all = free_o2")
print("[Ph+]")
print("function = " + exprToStr(Pi))
print("bndry_all = free_o2")
print("[Nh+]")
print("function = " + exprToStr(n))
print("bndry_all = free_o2")

