from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtCore import Qt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class MathTextLabel(QtWidgets.QWidget):

    def __init__(self, mathText, parent=None, **kwargs):
        super(QtWidgets.QWidget, self).__init__(parent, **kwargs)

        l = QVBoxLayout(self)
        l.setContentsMargins(0, 0, 0, 0)

        # r, g, b, a = self.palette().base().color().getRgbF()

        self._figure = Figure()
        self._canvas = FigureCanvas(self._figure)
        l.addWidget(self._canvas)
        self._figure.clear()
        text=self._figure.suptitle(
            mathText,
            x=0.0,
            y=1.0,
            horizontalalignment='left',
            verticalalignment='top',
            size=QtGui.QFont().pointSize()*2
        )
        self._canvas.draw()

        (x0, y0), (x1, y1) = text.get_window_extent().get_points()
        w = x1-x0
        h = y1-y0

        self._figure.set_size_inches(w/80, h/80)
        self.setFixedSize(w, h)

if __name__=='__main__':
    from sys import argv, exit

    class Widget(QtWidgets.QWidget):
        def __init__(self, parent=None, **kwargs):
            super(QtWidgets.QWidget, self).__init__(parent, **kwargs)

            l=QVBoxLayout(self)
            mathText=r'$X_k = \sum_{n=0}^{N-1} x_n . e^{\frac{-i2\pi kn}{N}}$'
            l.addWidget(MathTextLabel(mathText, self), alignment=Qt.AlignHCenter)

    a=QtWidgets.QApplication(argv)
    w=Widget()
    w.show()
    # w.raise_()
    exit(a.exec_())