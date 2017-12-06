import Spo2Sensor as Sp2
import time
from smbus2 import SMBus
import RPi.GPIO as GPIO
import Adafruit_MCP3008
import wiringpi
import Processing as pr
import numpy as np
from scipy import signal as sp
from gpiozero import Button
GPIO.setmode(GPIO.BCM)

class Sensor:
    
    def __init__(self,samplerateSpo2):
        #Definitions for ECG acquisition
        print "begin ECG config"
        CLK  = 11
        MISO = 9
        MOSI = 10
        CS   = 8
        self.mcp = Adafruit_MCP3008.MCP3008(CLK, CS, MISO, MOSI)
        print "ECG config ready"
       

        #Definitions fot Spo2 acquisition
        self.samplerate = samplerateSpo2
        self.Spo2 = Sp2.Spo2Sensor(sampleAvg= 4,sampleRate=self.samplerate)
        print "SpO2 config ready"

        #Array variables to store samples
        self.ecgValues = np.array([])
        self.Red = np.array([])
        self.IR = np.array([])
        self.Spo2Value = 0
        print "config done"


    def getECG(self, numSeconds):
        print "Begin ECG measure"
        sampleRate = 250
        samplePeriod = 1/250
        starTime = wiringpi.millis()

        while wiringpi.millis()-starTime < numSeconds*1000: 
            Ecg = round((self.mcp.read_adc(1)*3.3)/1024,3)
            self.ecgValues = np.append(self.ecgValues,Ecg)
            wiringpi.delayMicroseconds(samplePeriod)
    
    

    def getSpo2(self,numSeconds):
        print "begin SPO2 measure"
        startTime = wiringpi.millis()
        newSample = False
        AFthreshold= 17
        self.Spo2.enableAfull()
        self.Spo2.setFIFOAF(AFthreshold)
        interrupt  = Button(7)

        while wiringpi.millis()-startTime < numSeconds*1000:
            interrupt.when_activated = self.Spo2.sampleAvailable()
            if self.Spo2.newSample == True:
                self.Spo2.readSample()
                self.Spo2.newSample = False
                
        
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

        print "Buffer IR: ", len(self.Spo2.buffer_ir)
        print "Buffer Red: ", len(self.Spo2.buffer_red)

        pro = pr.Processing()

        #get Red and Ir buffers
        self.IR = self.Spo2.buffer_ir
        self.Red = self.Spo2.buffer_red

        #Delay Signal to avoid start overshoot
        #self.Red = pro.delaySignal(self.Red)
        #self.IR = pro.delaySignal(self. IR)

        #Median filter to the signals
        self.IR = sp.medfilt(self.IR)
        self.Red = sp.medfilt(self.Red)

        #low pass filter at 60hz
        #self.IR = pro.NotchFilter(self.IR, 60,self.samplerate)
        #self.Red = pro.NotchFilter(self.Red, 60,self.samplerate)

        #lowpass filter at 6Hz:
        self.IR = pro.lowPasFIRFilter(self.IR, 6,self.samplerate)
        self.Red = pro.lowPasFIRFilter(self.Red, 6,self.samplerate)
        
        #Compute Spo2Value:
        self.Spo2Value = pro.calcSpO2(self.Red,self.IR)
        print "Spo2: ", self.Spo2Value, "%"

        #get AC componente to plot the signal:
        self.Red = pro.getACcomponent(self.Red)
        self.IR = pro.getACcomponent(self.IR)
        #self.IR = pro.Normalize(self.IR)


        
        
        


    

    

    

        
    
    
            
        
    
        

        
        #      # print (wiringpi.millis()-startTime)/1000
        #       reD = Spo2Sensor.getRed()
        #       iR  = Spo2Sensor.getIR()
        #       #print "R: ", reD , "IR: ", iR
        #       self.Red = np.append(self.Red,reD)  
        #       self.IR = np.append(self.IR,iR)
    
        #   #self.Red = Spo2Sensor.lowPasFilter(self.Red,6,samplerate)
        #   #self.Red = Spo2Sensor.removeDC(self.Red)
    
    
        #   #self.IR = Spo2Sensor.lowPasFilter(self.IR,6,samplerate)
        #   #self.IR = Spo2Sensor.removeDC(self.IR)
    
    

        #   print "min IR:", min(self.IR)
        #   print "max IR:", max(self.IR)
        #   print "min RED:", min(self.Red)
        #   print "max RED: ", max(self.Red)

    











