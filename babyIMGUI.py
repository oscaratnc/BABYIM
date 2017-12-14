from Sensor import Sensor 
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg

numSecondsECG=2
numSecondsSpO2=3
sampleRateSpo2 = 400
Spo2Value = 0

Sensors = Sensor()
Sensors.configSensors(sampleRateSpo2)

print "Reading ECG Signal...."
Sensors.getECG(numSecondsECG)
print "ECG done"
print "Reading SPo2Value......"
Sensors.getSpo2read(numSecondsSpO2)
Spo2Value = Sensors.Spo2Valuecalc()
print "SPO2 done"


ECG = Sensors.ecgValues
RED  = Sensors.Red
IR = Sensors.IR


app = QtGui.QApplication([])

win = pg.GraphicsWindow()
win.resize(600,800)
win.setWindowTitle("Signals Ploting")

pg.setConfigOptions(antialias= True)
# p1 = win.addPlot(title="ECG")
# p1.plot(ECG, pen = pg.mkPen(color = 'r', width= 2))

win.nextRow()
p2 = win.addPlot(title = "RED LED")
p2.plot(RED, pen = pg.mkPen(color = 'g', width= 2))

if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
