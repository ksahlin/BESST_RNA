
#import sys
#from scipy.special import erf
import math
from Norm import erf
#from scipy.stats import norm
#from scipy.constants import pi
pi = 3.14159265359
def GapEstimator(mean,sigma,read_length,obs,c1_len,c2_len):
    d_ML=CalcMLvaluesOfdGeneral(mean,sigma,read_length,c1_len,c2_len,obs)
    #print "GAP=",d_ML, 'Data obs: ', obs
    return d_ML

def PreCalcMLvaluesOfdLongContigs(mean,stdDev,readLen):
    def Nom(z,mean,stdDev):
        nom=-(1+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        return nom
    def Denom(d,readLen,mean,stdDev):
        first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*(1+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))
        second=stdDev*2.718**(-((mean-d-2*readLen+1)**2)/(float(2*stdDev**2)))
        denom=first+second
        return denom

    def CalcAofd(d):
        #transform to z ~N(0,1)
        z=((d+2*readLen)-mean)/float(stdDev)    
        nom=Nom(z,mean,stdDev)
        denom=Denom(d,readLen,mean,stdDev)
        val=nom/denom              
        return val
    #we cannot span gaps outside this interval, for fair coverages with -e set to not too low value
    d_upper=int(mean+2*stdDev)
    d_lower=int(-2*stdDev)
    dValuesTable={}
    prev_obs=d_lower
    for d in range(d_lower,d_upper+1):
        Aofd=CalcAofd(d)
        obs=d+Aofd*stdDev**2
        obs=int(round(obs,0))
        dValuesTable[obs]=d
        #fill in missing values of the table here
        if obs - prev_obs > 1:
            n = obs - prev_obs
            for i in range(0,n):
                dValuesTable[prev_obs+i+1]=d
        prev_obs=obs
        
    return dValuesTable

def funcDGeneral(d,mean,stdDev,c1Len,c2Len,readLen):
    c_min=min(c1Len,c2Len)
    c_max=max(c1Len,c2Len)

    def Nominator(d,c_min,c_max,c1Len,c2Len,readLen):
        nomin1=( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) +erf((mean-d-c_max-readLen-1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        nomin2=-(erf((c_min+d+readLen+1-mean)/(2**0.5*float(stdDev)))+erf((mean-d-2*readLen+1)/(2**0.5*float(stdDev))))*(pi/2)**0.5
        nomin=nomin1+nomin2
        return nomin
    
    def Denominator(d,c1Len,c2Len,readLen):
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term2=(c_min-readLen+1)/2.0*(erf((c_max+d+readLen-mean)/((2**0.5)*stdDev))- erf((c_min+d+readLen-mean)/((2**0.5)*stdDev))   )
    
        first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*( erf((c_min+d+readLen-mean)/(2**0.5*float(stdDev))) - erf((d+2*readLen-1-mean)/(2**0.5*float(stdDev)))  )
        second=stdDev*( 2.718**(-( (d+2*readLen-1-mean)**2)/(float(2*stdDev**2))) - 2.718**(-( (c_min+d+readLen-mean)**2)/(float(2*stdDev**2)))) 
        term1=first+second

        first=((pi/2)**0.5)*(c_min+c_max+d-mean+1)*( erf((c_min+c_max+d-mean)/(2**0.5*float(stdDev))) - erf((c_max+readLen+d-mean)/(2**0.5*float(stdDev)))  )
        second=stdDev*( 2.718**(-( (c_min+c_max+d-mean)**2)/(float(2*stdDev**2))) - 2.718**(-( (c_max+readLen+d-mean)**2)/(float(2*stdDev**2))))       
        term3=first+second
        
        denom=term1+term2+term3
        return denom
    
    denominator=Denominator(d,c1Len,c2Len,readLen)
    nominator=Nominator(d,c_min,c_max,c1Len,c2Len,readLen)
    Aofd=nominator/denominator
    func_of_d=d+Aofd*stdDev**2
    return func_of_d

def CalcMLvaluesOfdGeneral(mean,stdDev,readLen,c1Len,c2Len,obs):
    #do binary search among values
    d_upper=int(mean+2*stdDev)
    d_lower=int(-2*stdDev)
    while d_upper-d_lower>1:
        d_ML=(d_upper+d_lower)/2.0
        func_of_d=funcDGeneral(d_ML,mean,stdDev,c1Len,c2Len,readLen)
        if func_of_d>obs:
            d_upper=d_ML
        else:
            d_lower=d_ML

    
    d_ML=(d_upper+d_lower)/2.0
    return int(round(d_ML,0))


