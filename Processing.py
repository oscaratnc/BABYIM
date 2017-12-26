from scipy import signal as sp
import numpy as np

class Processing:
    def movMean (self, signal, window):
        WinFilt = np.repeat(1.0,window)/window
        meanfilSignal = np.convolve(signal,WinFilt,'valid')
        return meanfilSignal  
    
    def BPButterFilter(self,signal,flow,fhigh,sampleF,order):
        sampleRate = sampleF
        nyq_rate = sampleRate/2.0
        Wn = [flow/nyq_rate,fhigh/nyq_rate]
        n = order
        [b,a] = sp.butter(n,Wn,'bandpass')
        filtered = sp.filtfilt(b,a,signal)
        return filtered

    def getACcomponent(self, measure):
        mean = np.mean(measure)
        measure = (measure-mean)
        return measure
    
    def getDCComponent(self,measure):
        DCcomponent  = np.mean(measure)
        return DCcomponent
    
    def calcSpO2(self,signalRed, signalIR):
        DCRed = self.getDCComponent(signalRed)
        acRed = self.getACcomponent(signalRed)
        DCIR = self. getDCComponent(signalIR)
        acIR = self.getACcomponent(signalIR)
        
        RR =round(np.mean((acRed/DCRed)/(acIR/DCIR)),4)
        print "RR: ", RR
        spO2Value =  96.545 + 0.616 * RR
        Spo2Value =int(np.round(spO2Value,0))
        return Spo2Value
        
    def Normalize(self, measure):
        abs = np.max(np.abs(measure))
        measureN = measure/abs
        measureN = np.round(measureN,4)
        return measureN



    def heartRateCalc(self, signalIR):
        for i in range(1,len(signalIR)):
            signalIR[i] = signalIR[i]-signalIR[i-1]
        
        # for y in range(len(signalIR)):
        #     max= np.max(signalIR)
        #     if signalIR[y] <.7*max:
        #         signalIR[y]=0

            return signalIR