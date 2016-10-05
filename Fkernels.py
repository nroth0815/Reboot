import numpy as np
import time
#import pynbody

class Kern(object): 

  def f(self, ik,jk,lk, iq,jq,lq, eps=1e-16):
  
    k=np.array([ik,jk,lk])
    q=np.array([iq,jq,lq])      
  
    #value=0.5*(ik*ik+jk*jk+lk*lk)/( (iq*iq+jq*jq+lq*lq)+eps ) * np.dot(q,k-q)/ (np.dot(k-q,k-q) +eps)  
    value=0.5*np.linalg.norm(k)/( np.linalg.norm(q) + eps ) * np.dot(q,k-q)/ (np.dot(k-q,k-q) +eps)  
  
    return value
    #print value
  
  
  def h(self, ik,jk,lk, iq,jq,lq, eps=1e-16):
  
    k=np.array([ik,jk,lk])
    q=np.array([iq,jq,lq])      
  
    value=np.dot(k, q)/(np.dot(q,q)+eps)
    return value


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


def test(res, t=False):

    kern=Kern()

    ind=0
    arr=np.zeros(res*res*res)

    for k1 in range(-res/2+1, res/2+1):
        for k2 in range(-res/2+1, res/2+1):
            for k3 in range(0, res/2+1):

                arr[ind]=kern.unq(res, np.array([k1, k2, k3]), t=t)
                ind+=1

    return arr






