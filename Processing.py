from scipy import signal as sp
import numpy as np

class Processing:
    def lowPasFIRFilter(self,signal,fc,sampleF):
        sampleRate = sampleF
        nyq_rate = sampleRate/2.0
        Wn = fc/nyq_rate 
        order = 100
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
    
    def ratioOfRatios(self, measureRed,measureIR):
        dcRed= self.getDCComponent(measureRed)
        acRed = self.getACcomponent(measureRed)
        dcIR = self.getDCComponent(measureIR)
        acIR = self.getACcomponent(measureIR)

        RR = (acRed/dcRed) / (acIR/dcIR)
        return RR

    def calcSpO2(self, measureRed, measureIR):
        RR = self.ratioOfRatios(measureRed, measureIR)
        spO2Array = 110-25 * RR
        Spo2Value =int( np.round(np.mean(spO2Array),0))

        return Spo2Value
        
    def Normalize(self, measure):
        abs = np.max(np.abs(measure))
        measureN = measure/abs
        measureN = np.round(measureN,4)
        return measureN
    
    def delbaselinedrift(self, measure, sampleF):
       #200 ms window
       n  = (200* sampleF)/1000
       k = 1
       for x in range(len(measure)):
           zeroline = sp.zeros(1,n)
           lim1 = x - n/2
           lim2 = x -(n/2) -1