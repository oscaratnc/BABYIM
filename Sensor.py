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
    
    def __init__(self):
        #Array variables to store samples
        self.ecgValues = np.array([])
        self.Red = np.array([])
        self.IR = np.array([])
        self.Spo2Value = 0
       
        
    def configSensors(self,samplerateSpo2):
        #Definitions for ECG acquisition
        print 
        CLK  = 11
        MISO = 9
        MOSI = 10
        CS   = 8
        self.mcp = Adafruit_MCP3008.MCP3008(CLK, CS, MISO, MOSI)
        print "ECG config ready"

        #Definitions fot Spo2 acquisition
        self.samplerateSpo2 = samplerateSpo2
        self.Spo2 = Sp2.Spo2Sensor(sampleAvg= 4,sampleRate=self.samplerateSpo2)
        print "SpO2 config ready"

    def getECG(self, numSeconds):
        print "Begin ECG measure"
        sampleRate = 250
        samplePeriod = (1/250)*1000
        starTime = wiringpi.millis()

        while wiringpi.millis()-starTime < numSeconds*1000: 
            Ecg = round((self.mcp.read_adc(1)*3.3)/1024,3)
            self.ecgValues = np.append(self.ecgValues,Ecg)
            wiringpi.delay(samplePeriod)
            
    def getSpo2read(self,numSeconds):
        print "begin SPO2 measure"
        startTime = wiringpi.millis()
        newSample = False
        AFthreshold= 20
        self.Spo2.enableAfull()
        self.Spo2.setFIFOAF(AFthreshold)
        interrupt  = Button(7)

        while wiringpi.millis()-startTime < numSeconds*1000:
            interrupt.when_activated = self.Spo2.sampleAvailable()
            if self.Spo2.newSample == True:
                self.Spo2.readSample()
                self.Spo2.newSample = False    
            print (wiringpi.millis()-startTime)/1000
          
        print "Spo2 measure ready"
                
        print "processing signal...."
        print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

        print "SF Reading IR: ", len(self.Spo2.buffer_ir)
        print "SF Reading Red: ", len(self.Spo2.buffer_red)
        
        pro = pr.Processing()

        #get Red and Ir buffers
        self.IR = self.Spo2.buffer_ir
        self.Red = self.Spo2.buffer_red
        

        # #Normalize Red and IR signals}
        self.Red = pro.Normalize(self.Red)
        self.IR = pro.Normalize(self.IR)
 
        #Mean filter widnow = 4
        self.IR  = pro.movMean(self.IR,4)
        self.Red = pro.movMean(self.Red,4)

        #Butterword 4th order bandpass filter .5-6Hz
        self.IR = pro.BPButterFilter(self.IR,0.5,6.0,self.samplerateSpo2)
        self.Red = pro.BPButterFilter(self.Red,0.5,6.0,self.samplerateSpo2)
        self.generateDataFile()

    def Spo2Valuecalc(self):
        # #Compute Spo2Value:
        self.Spo2Value = pro.calcSpO2(Red,IR)
        print "Spo2: ", self.Spo2Value, "%"
        return Spo2Value
        
    def generateDataFile(self):
        
        file = open("DataFile_Read.csv", "w")
        
        file.write(str(self.Red.tolist()))
        file.write("\n")
        file.write(str(self.IR.tolist()))
        file.write("\n")
        file.write(str(self.ecgValues.tolist()))
        file.close()

        print "File Created"
        return True


        

        
        
        


    

    

    

        
    
    
            
        
    
        

        
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

    











