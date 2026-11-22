#!/usr/bin/env python3
"""Optional oracle check. Install scipy and mpmath separately; neither is a library dependency.
Pass the absolute path to the built phylogenetics shared library as the first argument.
"""
import sys,ctypes as c,math,random
from pathlib import Path
import scipy.special as s
import mpmath as mp
lib=c.CDLL(str(Path(sys.argv[1]).resolve()))
for n,k in [('pgamma',3),('qgamma',3),('pnorm',3),('qnorm',3),('dweibull',3),('pweibull',3),('qweibull',3),('dlogis',3),('plogis',3),('qlogis',3)]:
 f=getattr(lib,'crpx_'+n);f.argtypes=[c.c_double]*k+[c.c_bool];f.restype=c.c_double
rng=random.Random(321);worst=0;bad=[]
for i in range(5000):
 a=10**rng.uniform(-5,5);x=a*10**rng.uniform(-3,2)
 actual=lib.crpx_pgamma(x,a,1,True);reference=s.gammainc(a,x)
 if not math.isfinite(actual):bad.append((a,x,actual,'nonfinite log CDF'));continue
 if reference>1e-290:
  expected=math.log(reference);error=abs(actual-expected)/max(1,abs(expected));worst=max(error,worst)
  if not math.isfinite(actual) or error>3e-11:bad.append((a,x,actual,expected))
print('5000 CDF grid worst relative log error',worst,'failures',bad[:5]);assert not bad
bad=[]
for i in range(1000):
 a=10**rng.uniform(-2,5);p=10**rng.uniform(-200,-.0001)
 actual=lib.crpx_qgamma(p,a,1,False);expected=s.gammaincinv(a,p)
 if not math.isfinite(actual) or actual<0:bad.append((a,p,actual,expected));continue
 if expected>1e-290 and (not math.isfinite(actual) or abs(actual/expected-1)>5e-10):bad.append((a,p,actual,expected))
print('1000 quantile grid failures',bad[:5]);assert not bad
mp.mp.dps=70
for lp in [-1000,-1e6]:
 q=lib.crpx_qnorm(lp,0,1,True);recovered=mp.log(mp.erfc(-mp.mpf(q)/mp.sqrt(2))/2)
 assert abs(float(recovered)-lp)<max(1,abs(lp))*3e-16
print('high precision log-normal tails passed')
