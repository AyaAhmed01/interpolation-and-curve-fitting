import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
from PyQt5.QtWidgets import QLineEdit, QLabel, QPushButton, QPlainTextEdit, QSlider
from PyQt5 import QtCore, QtGui, uic, QtWidgets
import pyqtgraph as pg
import os
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator



class MyWidget(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        uic.loadUi('gui.ui', self)
        bg = self.palette().window().color()
        cl = (bg.redF(), bg.greenF(), bg.blueF())
        self.fig = Figure(edgecolor=cl, facecolor=cl)
        self.error = Figure(edgecolor=cl, facecolor=cl)
        self.canvas_equation = FigureCanvas(self.fig)
        self.expression.addWidget(self.canvas_equation)
        self.canvas_error_map = FigureCanvas(self.error)
        self.error_map_view.addWidget(self.canvas_error_map)
        self.select_chunk_Box.valueChanged.connect(lambda: self.ch_equation())
        self.actionopen.triggered.connect(lambda: self.Open())
        self.Order_text.textChanged.connect(lambda: self.Order())
        self.NumChunks_text.textChanged.connect(lambda: self.fitting())
        self.extrapolation_text.textChanged.connect(lambda: self.extrapolation())
        self.start_button.clicked.connect(lambda: self.start_errormap())
        self.pens = [pg.mkPen('r'), pg.mkPen('b'), pg.mkPen('g')]



def extrapolation(self):
        self.equ_arr.clear()
        self.NumChunks_text.setText("1")
        self.extra = int(self.extrapolation_text.text())
        order, num_chunks = self.get_values()
        self.size = int(len(self.x_data) * self.extra / 100)   # size taken from gui, 
        self.new_x = self.x_data[:self.size - 1]     # to make them interpolation
        self.new_y = self.y_data[:self.size - 1]  
        self.extrapolated_x = self.x_data[self.size:]  # x to make extrapolation for
        coeffs = poly.polyfit(self.new_x, self.new_y, order)   # 
        x_fitline = np.linspace(self.new_x[0], self.new_x[-1], num=len(self.new_x) * 10)
        x_fitting = np.array(x_fitline)
        x_arr= np.append(x_fitting, self.extrapolated_x)
        y_fitline = poly.polyval(x_arr, coeffs)
        self.maingraph.clear()
        self.maingraph.plot(self.x_data, self.y_data, pen=self.pens[2])
        self.maingraph.plot(x_arr, y_fitline, pen=None, symbol='o')
        for i in range(0, len(self.new_y)):
            self.error_per_chunk += ((self.new_y[i] - y_fitline[i])**2)
        fx = str(round(coeffs[0], 4))
        for i in range(1, len(coeffs)):
            if coeffs[i] > 0:
                const = " + " + str(round(coeffs[i], 4)) + " x^{{{}}}".format(i)
            else:
                const = " - " + str(abs(round(coeffs[i], 4))) + " x^{{{}}}".format(i)
            fx += const
        error = "\n, Percentage Error = " + str(round(self.error_per_chunk, 3)) + " %"
        fx += error
        self.equ_arr.append(fx)
        print(self.equ_arr)
        self.ch_equation()

def main():
    app = QtWidgets.QApplication(sys.argv)
    application = MyWidget()
    application.show()
    app.exec_()

if __name__ == "__main__":
    main()