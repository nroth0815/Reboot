import numpy as np
import time
#import pynbody

class Kern(object): 

  def f(self, ik,jk,lk, iq,jq,lq, eps=1e-16):
  
    k=np.array([ik,jk,lk])
    q=np.array([iq,jq,lq])
    p=k-q  
  
    #value=0.5*(ik*ik+jk*jk+lk*lk)/( (iq*iq+jq*jq+lq*lq)+eps ) * np.dot(q,k-q)/ (np.dot(k-q,k-q) +eps)  
    value=0.5*np.linalg.norm(k)**2./( np.linalg.norm(q)**2. + eps )*np.dot(q,p)/(np.linalg.norm(p)**2 +eps)  
  
    return value
    #print value
  
  
  def h(self, ik,jk,lk, iq,jq,lq, eps=1e-16):
  
    k=np.array([ik,jk,lk])
    q=np.array([iq,jq,lq])      
  
    value=np.dot(k, q)/(np.dot(q,q)+eps)
    return value

  def kernel(self, q1, q2, q3, p1, p2, p3):
    q=np.array([q1,q2,q3])
    p=np.array([p1, p2, p3])
    resu= 5./7. + 2./7. * np.dot(q,p)**2/(np.linalg.norm(q)**2*np.linalg.norm(p)**2+1e-16) + 1./2. * np.dot(q,p)*(1./(np.linalg.norm(q)**2+1e-16) +1./(np.linalg.norm(p)**2+1e-16 ) )
    return resu

  def unq(self, res, k, t=False):

    ind=0
    arr=np.zeros(res*res*res)

    if t == True:
        start=time.time()

    for i1 in range(-res/2+1, res/2+1):
        p1=k[0]-i1
        if p1 > res/2 or p1 < -res/2+1:
            #print i1, p1
            continue
        for i2 in range(-res/2+1, res/2+1):
            p2=k[1]-i2
            if p2 > res/2 or p2 < -res/2+1:
                continue
            for i3 in range(-res/2+1, res/2+1):
                p3=k[2]-i3
                if p3 > res/2 or p3 < -res/2+1:
                    continue

                arr[ind]=5*self.h(k[0], k[1], k[2], i1, i2, i3 ) + 5*self.h(k[0], k[1], k[2],  p1, p2, p3 ) + 4*self.f(k[0], k[1], k[2], i1, i2, i3)
                ind+=1

    if t == True:
        stop=time.time()
        print "Seconds passed:", stop-start

    uq=np.unique(arr)

    return float(len(uq))/res**3

    #print float(len(uq))/res**3
    #return uq


def test2(res, t=False):

    kern=Kern()

    ind=0
    arr=np.zeros(res*res*res)

    for k1 in range(-res/2+1, res/2+1):
        for k2 in range(-res/2+1, res/2+1):
            for k3 in range(0, res/2+1):

                arr[ind]=kern.unq(res, np.array([k1, k2, k3]), t=t)
                ind+=1

    return arr


#0.142857*( 0.000181971*3.85782e-05 -  0.000120421*  0.000113936 ) +  0.714286* (0.0006482857* -0.733568 - 7.11065e-06*0)+ 0.142857*(0.000129346* 3.85782e-05 - -6.87027e-05* -0.000113936)+0.285714*( -0.000141438* -0.000183734- -8.82801e-06*-2.84072e-05)+0.714286 *(-0.000460343* -0.000632829- 0.000426549*7.52807e-05 )+0.285714*(5.13121e-05*-0.000287232 - -0.000221203*4.6093e-05)+0.285714 *(-0.000287232*5.13121e-05- 4.6093e-05*-0.000221203)+0.714286 *(-0.000632829*-0.000460343- 7.52807e-05*0.000426549)+0.285714*(-0.000183734*-0.000141438--2.84072e-05*-8.82801e-06)+0.142857*(3.85782e-05*0.000129346--0.000113936*-6.87027e-05)+0.714286*(-0.733568*0.000648001 -0*7.11065e-06 )+0.142857*( 3.85782e-05*0.000181971 - 0.000113936*0.000120421)

def test(k, q):

 kern=Kern()

 print k, " | ", q

 print (5*kern.h(k[0],k[1],k[2], q[0],q[1],q[2])+5*kern.h(k[0],k[1],k[2], k[0]-q[0],k[1]-q[1],k[2]-q[2])+4*kern.f(k[0], k[1], k[2], q[0], q[1], q[2]))/14.
 print kern.kernel(q[0], q[1], q[2], k[0]-q[0], k[1]-q[1], k[2]-q[2])
 #print kern.kernel(k[0]-q[0], k[1]-q[1], k[2]-q[2], q[0], q[1], q[2])

