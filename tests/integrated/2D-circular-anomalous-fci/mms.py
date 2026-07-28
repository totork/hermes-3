from __future__ import division
from __future__ import print_function

from boutdata.mms import Metric, sin, cos, Div_par, Grad_par, exprToStr, diff, y, t, DDY, x, z, DDX, DDZ, Delp2, D2DX2, D2DZ2
from math import pi

# Length of the y domain
Ly = 2.0 * pi

# Atomic mass number
AA = 1.0

# metric tensor
metric = Metric()  # Identity

qe = 1.60217663e-19
Me = 9.1093837e-31
e0 = 8.85418781e-12
Pi = pi

Tnorm = 5
Nnorm = 1e18
Bnorm = 1.0
Omega_ci = qe * Bnorm / (1836.0*Me);
rho_s = 0.0002284697436697996

#n = 1 + 0.1*sin(2 * pi * x) * sin(3 * z)# * sin(t*omega)
n = 5.0 + 0.1 * sin(2.0 * pi * x )#  * sin(2 * z+ 0.231)


D = (1.0 +  0.0 * x)/(rho_s*rho_s*Omega_ci)


xmin = 0.1
xmax = 0.2
Lx = xmax - xmin
rad = (xmin + x * Lx)/rho_s



#x_new = (metric.x + xmin)/(xmax-xmin)
#z_new = metric.z * 2.0 * pi / (2.0 * pi * x_new)
#replace = [ (x, x_new ), (z, z_new)  ]
replace = [ (x, metric.x / Lx) ]  

n = n.subs(replace)
D = D.subs(replace)

##############################
# Calculate time derivatives

#dndt = D * Delp2(n) * rho_s**2   
dndt = D * ( 1.0 / rad * DDX(rad * (DDX(n)))) * rho_s**2
#############################
# Calculate sources

Sn = diff(n, t) - dndt
# Substitute back to get input y coordinates
#replace = [ ( metric.x , 1.0 / x_new), (metric.z, 1.0 / z_new)]
replace = [ ( metric.x , x*Lx) ]
n = n.subs(replace)
Sn = Sn.subs(replace)

print("[n]")
print("solution = " + exprToStr(n))
print("\nsource = " + exprToStr(Sn))
