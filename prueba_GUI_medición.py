import pyqtgraph as pg 
from pyqtgraph import QtGui, QtCore


app = QtGui.Qapplication(sys.argv)
pg.plot(x = [0,1,2,3,4], y = [4,5,9,6])

status = app.exec_()
sys.exit(status)