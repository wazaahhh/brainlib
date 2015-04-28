import numpy as np
from scipy import stats
from scipy.interpolate import interp1d

'''Here all functions related to the analysis of EEG signal, 
using methods based on the power spectrum'''


def pSpectrum(vector):
    
    A = np.fft.fft(vector)
    ps = np.abs(A)**2
    ps = ps[:len(ps)/2]
        
    return ps

#ps = pSpectrum(test[500:600])
#loglog(ps)

def entropy(power_spectrum,q):
    q = float(q)
    
    power_spectrum = np.array(power_spectrum)
        
    if not q ==1:
        S = 1/(q-1)*(1-np.sum(power_spectrum**q))
    else:
        S = - np.sum(power_spectrum*np.log2(power_spectrum))
        
    return S

def entropySeries(time,rawData,q,windowSize=1,normalize=True):
    
    
    #t = tSeries['rawTime']
    #raw = tSeries['rawValue']
    
    t = np.array(time)
    raw = np.array(rawData)
    
    intervales = np.arange(np.trunc(min(t)+1),np.trunc(max(t)+1),windowSize)

    T = []
    S = [] 
    spectra = []
    
    for i,ix in enumerate(intervales[:-1]):
        #print i,ix
        
        c = (t>=ix)*(t<intervales[i+1])
        
        vector = raw[c]
        ps = np.array(pSpectrum(vector)) # Power Spectrum
        
        try:
            spectra += [ps]
        except:
            spectra = [ps]

        
        
        #if normalize:
        ps = ps/np.sum(ps)
        
        s = entropy(ps,q) #Compute Entropy
        S.append(s)
        T.append(i)

    S = np.array(S)
    T = np.array(T)
    
    if normalize:
        S = (S-np.mean(S))/np.std(S)
        
    return {'time' : list(intervales[:-1]), 'S': list(S),'spectra':spectra}


def avgPowerSpectrum(pspectra,xOutput):
    '''Computes an average power spectrum
    from multiple power spectra, and returns 
    (interpolated) evaluations for xOutput values '''
    
    if isinstance(pspectra,dict):
        pspectra = pspectra.values()
    
    l = len(pspectra)
    array = np.zeros([l,len(xOutput)])
    S = []
    
    for i,ps in enumerate(pspectra):
        #print i
        s = entropy(ps/np.sum(ps),1)
        S.append(s)
        
        x = np.arange(1,len(ps)+1)
        f = interp1d(x,ps/np.sum(ps))
        try:
            array[i] = f(xOutput)
        except:
            array[i,0]=-1
            continue
        
    index =np.argwhere(array[:,0]==-1)
    array = np.delete(array,index,0)
    S = np.delete(S,index,0)
    return array,S


def pSpectraCharacterize(pSpectrDic,index,plot=False,color='b'):
    '''computes power spectrum characteristics 
    for a dictionary/array of power spectra'''
    
    slope1 = []
    slope2 = []
    S = []
    
    for i in index:
        try:
            dic = pinkNoiseCharacterize(pSpectrDic['spectra'][i])
        except:
            print "error index %s"%i
            continue
        
        slope1 = np.append(slope1,dic['slope1'])
        slope2 = np.append(slope2,dic['slope2'])
        S = np.append(S,pSpectrDic['S'][i])
              
    if plot:
        pl.subplot(131)
        fit = stats.linregress(slope1,S)
        print fit
        pl.plot(slope1,S,'.',c=color)
        pl.plot(slope1,slope1*fit[0]+fit[1],ms='-',c=color)
        pl.xlabel("slope1")
        pl.ylabel("normalized entropy")
        
        pl.subplot(132)
        pl.plot(slope2,S,'.',c=color)
        pl.xlabel("slope2")
        pl.ylabel("normalized entropy")
        
        pl.subplot(133)
        pl.plot(slope1,slope2,'.',c=color)
        pl.xlabel("slope1")
        pl.ylabel("slope2")
        
    return {'S':S,'slope1':slope1,'slope2':slope2}

def pinkNoiseCharacterize(pspectrum,normalize=True,plot=False):
    '''Compute main power spectrum characteristics'''
    if normalize:
        pspectrum = pspectrum/np.sum(pspectrum)
    
    S = entropy(pspectrum,1)
    
    x = np.arange(1,len(pspectrum)+1)
    lx = np.log10(x)
    ly = np.log10(pspectrum)
    
    c1 = (x > 0)*(x < 80)
    c2 = x >= 80
    
    fit1 = stats.linregress(lx[c1],ly[c1])
    fit2 = stats.linregress(lx[c2],ly[c2])
    
    #print fit1
    #print fit2
    
    if plot:
        plot(lx,ly)
        plot(lx[c1],lx[c1]*fit1[0]+fit1[1],'r-')
        plot(lx[c2],lx[c2]*fit2[0]+fit2[1],'r-')
        
    return {'S':S,'slope1':fit1[0],'slope2':fit2[0]}
