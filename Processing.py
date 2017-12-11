from scipy import signal as sp
from scipy import ndimage as nd
import numpy as np

class Processing:
    def lowPasFIRFilter(self,signal,fc,sampleF):
        sampleRate = sampleF
        nyq_rate = sampleRate/2.0
        Wn = fc/nyq_rate 
        order = 50
        a = sp.firwin(order, Wn)
        filtered = sp.filtfilt(a,1.0,signal)
        return filtered

    def highPassFIRFilter(self, signal,fc,sampleF):
        sampleRate = sampleF
        nyq_rate = sampleRate/2.0
        Wn = fc/nyq_rate 
        order = 99
        a = sp.firwin(order, Wn, pass_zero=False)
        filtered = sp.filtfilt(a,1.0,signal)
        return filtered
    
    def NotchFilter(self, signal, fc, sampleF):
        nyqRate = sampleF/2.0
        f1 = (fc-.5)/nyqRate
        f2 = (fc+.5)/nyqRate
        bandwidth = [f1,f2]
        order = 49
        a = sp.firwin(order, bandwidth)
        Filtered = sp.filtfilt(a,1,signal)
        return Filtered

  
    def getACcomponent(self, measure):
        mean = np.mean(measure)
        measure = measure-mean
        return measure
    
    def getDCComponent(self,measure):
        DCcomponent = np.mean(measure)
        return DCcomponent
    
    def calcSpO2(self, acRed,acIR,DCIR,DCRed):
        RR = (acRed/DCRed) / (acIR/DCIR)
        spO2Array = 110-25 * RR
        Spo2Value =int( np.round(np.mean(spO2Array),0))

        return Spo2Value
        
    def Normalize(self, measure):
        abs = np.max(np.abs(measure))
        measureN = measure/abs
        measureN = np.round(measureN,4)
        return measureN
    
    # def delbaselinedrift(self, measure, sampleF):
    #         #200 ms window for the first time
    #     n  = (200* sampleF)/1000
    #     line1 = np.array([])
    #     for k in range(len(measure)):
    #         line0 = np.array([]) 
    #         lim1 = k - n/2
    #         lim2 = k + (n/2) -1
    #         if lim2 > len(measure):
    #             lim2 = len(measure)
    #         if lim1 <= 0:
    #             for k2 in range((n/2)+lim1,lim2):
    #                 line0 = np.append(line0, measure[k])            
    #         else: 
    #             k1 = 0
    #             for k2 in range (lim1,lim2):
    #                 line0 = np.append(line0,measure[k2])    
    #         line0 = np.sort(line0)
    #         mean = nd.median_filter(line0,size = (1,len(line0)))
    #         line1 = np.append(line1,mean)

            
        # #600 ms window fot the second time
        # n = (600 *sampleF)/1000
        # line2 = np.array([])
        # for k in range(len(measure)):
        #     line0 = np.array([]) 
        #     lim1 = k - n/2
        #     lim2 = k + (n/2) -1
        #     if lim2 > len(measure):
        #         lim2 = len(measure)
        #     if lim1 <= 0:
        #         for k2 in range((n/2)+lim1,lim2):
        #            line0 = np.append(line0,line1[k])
        #     else: 
        #         k1 = 0
        #         for k2 in range (lim1,lim2):
        #             line0 =np.append(line0,line1[k2])
            
        #     line0 = np.sort(line0)
        #     mean = nd.median_filter(line0,size = (1,len(line0)))
        #     line2 = np.append(line2,mean)
       
        # measure = measure- line2
    
        # return measure