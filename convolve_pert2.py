import sys
import numpy

data = numpy.loadtxt('scr/wfirst_cl.dat')

lmax = numpy.shape(data)[0]+1
ncol = numpy.shape(data)[1]-1

for l in range(2,lmax):
  data[l-2,1:] /= l*(l+1)/2/numpy.pi
  data[l-2,1:] *= (l*(l+1)/2.)**2

data2 = numpy.zeros((ncol))

# Cm
Cm = numpy.zeros((lmax-1))
#for l in range(8,1442): Cm[l-2] = 2*numpy.pi/l/(l+1)
Cm[int(sys.argv[1])-2] = 1.
sum = 0.
for l in range(2,lmax): sum += Cm[l-2] * (2*l+1)/(4*numpy.pi)
Cm /= sum

ntheta = 512

for l in range(2,lmax+1):
  s='{:4d}'.format(l)
  data2[:] = 0.
  for lpp in range(2,lmax+1):
    if Cm[lpp-2]>1e-49:
      for itheta in range(ntheta):
        theta = (itheta+.5)/ntheta*numpy.pi
        lpx = l-lpp*numpy.cos(theta)
        lpy = lpp*numpy.sin(theta)
        lp = int(numpy.floor(numpy.sqrt(lpx**2+lpy**2)))
        if lp>=2 and lp<=lmax: data2 += data[lp-2,1:] * Cm[lpp-2] * ((lpx**2-lpy**2)/(lpx**2+lpy**2))**2 * lpp
  data2 *= 1./(2.*numpy.pi)/ntheta
  for i in range(ncol): s+=' {:19.12E}'.format(data2[i])
  print(s)
