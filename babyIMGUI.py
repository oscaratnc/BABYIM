from Sensor import Sensor 
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg

numSecondsECG=2
numSecondsSpO2=15
sampleRateSpo2 = 200


Sensors = Sensor(sampleRateSpo2)
print "ECG reading"
Sensors.getECG(numSecondsECG)
print "ECG done"
Sensors.getSpo2(numSecondsSpO2)
print "SPO2 done"

#Sensors.generateDataFile()


ECG = Sensors.ecgValues
RED  = Sensors.Red
IR = Sensors.IR
#PPG = RED/IR
#print ECG
#print RED
#print IR

app = QtGui.QApplication([])

win = pg.GraphicsWindow()
win.resize(800,600)
win.setWindowTitle("Signals Ploting")

pg.setConfigOptions(antialias= True)
p1 = win.addPlot(title="ECG")
p1.plot(ECG, pen = pg.mkPen(color = 'r', width= 2))

win.nextRow()
p2 = win.addPlot(title = "RED LED")
p2.plot(RED, pen = pg.mkPen(color = 'g', width= 2))


win.nextRow()
p3= win.addPlot(title = "IR LED")
p3.plot(IR, pen = pg.mkPen(color = 'b', width= 2))

if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()