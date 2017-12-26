from Sensor import Sensor 
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg

numSecondsECG=2
numSecondsSpO2=5
sampleRate = 400.0
T= 1/sampleRate


Sensors = Sensor()
Sensors.configSensors(sampleRate)

print "Reading ECG Signal...."
Sensors.getECG(numSecondsECG,sampleRate)
print "ECG done"
print "Reading SPo2Value......"
Sensors.getSpo2read(numSecondsSpO2)


Spo2Value = Sensors.Spo2Value
ECG = Sensors.ecgValues
RED  = Sensors.Red
IR = Sensors.IR
HR = Sensors.HR
print HR

fftFILTIR = Sensors.IR_Filtered_FFT
fftFiltRed = Sensors.Red_Filtered_FFT

app = QtGui.QApplication([])

win = pg.GraphicsWindow() 
win.resize(600,500)
win.setWindowTitle("SpO2 Original")
pg.setConfigOptions(antialias= True)


p1 = win.addPlot(title="Spo2  IR  Original")
p1.plot(IR, pen = pg.mkPen(color = 'r', width= 2))

win.nextRow()
p2 = win.addPlot(title = "SpO2 Red Original")
p2.plot(HR,pen = pg.mkPen(color = 'w', width= 2))


win.nextRow()

N = len(fftFILTIR)
xfIR = np.linspace(0.0,30.0)
p3 = win.addPlot(title = "fft IR Filtered")
yfft =  np.abs(fftFILTIR[0:50])
p3.plot(xfIR, yfft,pen = pg.mkPen(color = 'g', width= 2))

win.nextRow()
p4 = win.addPlot(title="FFT RED Filtered")
N = len(fftFiltRed)
xfRed = np.linspace(0.0,30.0)
yfft= np.abs(fftFiltRed[0:50])
p4.plot(xfRed,yfft, pen = pg.mkPen(color = 'b', width= 2))


if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
