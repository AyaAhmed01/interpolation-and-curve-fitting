from threading import Thread
import sys
import threading
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
from PyQt5.QtWidgets import QLineEdit, QLabel, QPushButton, QPlainTextEdit, QSlider,QMessageBox
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
        uic.loadUi('components.ui', self)
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
        self.start_button.clicked.connect(lambda: self.start_or_end_errormap())
        self.pens = [pg.mkPen('r'), pg.mkPen('b'), pg.mkPen('g')]

    def Open(self):
        self.file = QtWidgets.QFileDialog.getOpenFileNames(
           self, 'Open only txt or CSV or xls', os.getenv('HOME'))
        data = pd.read_csv(self.file[0][0])
        self.y_data = data.values[:, 1]
        self.x_data = data.values[:, 0]
        self.maingraph.plot(self.x_data, self.y_data, pen=self.pens[2])
        self.Order_text.setText("0")
        self.NumChunks_text.setText("1")
        self.extrapolation_text.setText("100")
        self.fitting()
        self.ch_equation()

    def Order(self):
        if self.extrapolation_text.text() != ("100" or"0"):
            self.NumChunks_text.setText("1")
            self.extrapolation()
        else:
            self.fitting()

    def get_values(self):
        order = int(self.Order_text.text())
        num_chunks = int(self.NumChunks_text.text())
        return order, num_chunks

    x_each_chunk = []
    chunk_index = []
    chunks_points = []
    y_each_chunk = []
    equ_arr = []
    exit_error_map = False

    def fitting(self):
        order, num_chunks = self.get_values()
        chunk_width = (max(self.x_data) - min(self.x_data)) / num_chunks
        self.chunk_index.clear()
        self.equ_arr.clear()
        for i in range(0, num_chunks + 1):
            ind = min(self.x_data) + (i * chunk_width)
            self.chunk_index.append(ind)

        self.maingraph.clear()
        self.maingraph.plot(self.x_data, self.y_data, pen=self.pens[2])
        self.new_chunk = 0
        for j in range(1, len(self.chunk_index) + 1):
            self.x_each_chunk.clear()
            self.y_each_chunk.clear()
            for k in range(self.new_chunk, len(self.x_data)):
                if self.x_data[k] <= self.chunk_index[j]:
                    self.x_each_chunk.append(self.x_data[k])
                    self.y_each_chunk.append(self.y_data[k])
                else:
                    break
            self.new_chunk = k
            self.interpolation(self.x_each_chunk, self.y_each_chunk, order)

    def interpolation(self, x, y, order):
        coeffs = poly.polyfit(x, y, order)
        x_fitline = np.linspace(x[0], x[-1], num=len(x) * 10)
        y_fitline = poly.polyval(x_fitline, coeffs)
        self.error_per_chunk = 0
        self.maingraph.plot(x_fitline, y_fitline, pen=None, symbol='o')
        for i in range(0, len(y)):
            self.error_per_chunk += ((y[i] - y_fitline[i])**2)
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
        return coeffs

    def ch_equation(self):
        order, num_chunks = self.get_values()
        self.select_chunk_Box.setMinimum(1)
        self.select_chunk_Box.setMaximum(num_chunks)
        box_value = self.select_chunk_Box.value()
        eqn = self.equ_arr[box_value-1]
        self.fig.clear()
        self.fig.suptitle("$f(x) = {{{}}}$".format(str(eqn)),
                          x=0.0, y=0.5,
                          horizontalalignment='left',
                          verticalalignment='center')
        self.canvas_equation.draw()

    extrapolated_x = []
    
    def extrapolation(self):
        self.equ_arr.clear()
        self.extra = int(self.extrapolation_text.text())
        order = self.get_values()        
        self.size = int(len(self.x_data) * self.extra / 100)   # size taken from gui, 
        self.new_x = self.x_data[:self.size - 1]                # to make them interpolation
        self.new_y = self.y_data[:self.size - 1]  
        self.extrapolated_x = self.x_data[self.size:]           # x to make extrapolation for
        coeffs = poly.polyfit(self.new_x, self.new_y, order)   # array of coeff. of interpolation as (a + bx + cx^2 + dx^3 +.....) 
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
    
    def start_or_end_errormap(self):
        if self.start_button.text() == "Start":
            self.exit_error_map = False
            x_param = self.x_error_map.currentText()
            y_param = self.y_error_map.currentText()
            if x_param == y_param:  # blank errormap if same attributes on both axis
                msg = QMessageBox()
                msg.setWindowTitle("Error")
                msg.setText("X-axis and Y-axis have the same variable, please change any one of them.")
                msg.exec_()
                return
            print("before starting thread\n")
            map_thread = threading.Thread(target= self.start_errormap)
            map_thread.start()
        else:          # end thread
            self.start_button.setText("Start")
            self.error.clear()
            self.canvas_error_map.draw()
            print("before ending thread\n")
            self.exit_error_map = True 

    def start_errormap(self):
        self.error.clear()
        self.canvas_error_map.draw()
        self.x_param = self.x_error_map.currentText()
        self.y_param = self.y_error_map.currentText()
        print("in start errormap\n")
        self.complete = 0
        self.start_button.setText("cancel")
        
        x_range = int(self.x_max.text())
        y_range = int(self.y_max.text())
        self.step = (100 / (x_range*y_range))
        print("after step\n")
        if self.x_param == 'order':                  # x_param, y_param are taken from user
            if self.y_param == 'number of chunks':    # chunk on y 
                self.overlap_vals = int(self.const_text.text()) / 100.0
                self.chunk_vals = np.arange(1, y_range+1)
                self.order_vals = np.arange(0, x_range+1)
                self.create_errormap(self.order_vals, self.chunk_vals)
            else:    # y is overlapping 
                self.chunk_vals = int(self.const_text.text())
                self.overlap_vals = np.arange(0, y_range+1, 2) / 100.0
                self.order_vals = np.arange(0, x_range+1)
                self.create_errormap(self.order_vals, self.overlap_vals)
            

        if self.x_param == 'number of chunks':
            if self.y_param == 'order':    # order on y 
                self.overlap_vals = int(self.const_text.text()) / 100.0
                self.chunk_vals = np.arange(1, x_range+1)
                self.order_vals = np.arange(0, y_range+1)
                self.create_errormap(self.chunk_vals, self.order_vals)
            else:    # y is  overlapping
                self.order_vals = int(self.const_text.text())
                self.overlap_vals = np.arange(0, y_range+1, 2) / 100.0
                self.chunk_vals = np.arange(1, x_range+1)
                self.create_errormap(self.chunk_vals, self.overlap_vals)
        
        if self.x_param == 'overlapping':
            if self.y_param == 'order':    # order on y 
                self.chunk_vals = int(self.const_text.text())
                self.overlap_vals = np.arange(0, x_range+1, 2) / 100.0
                self.order_vals = np.arange(0, y_range+1)
                self.create_errormap(self.overlap_vals, self.order_vals)
            else:    # y is  chunks
                self.order_vals = int(self.const_text.text()) 
                self.chunk_vals = np.arange(1, y_range+1)
                self.overlap_vals = np.arange(0, x_range+1, 2) / 100.0
                self.create_errormap(self.overlap_vals, self.chunk_vals)
    
    
    def create_errormap(self, x_axis, y_axis):
        self.mat_xaxis = x_axis
        self.mat_yaxis = y_axis
        self.error_array = []
        progress_state = self.calculate_errormap()  # fills the error_array
        if progress_state == False:      # exit error map if cancel button is pushed
            self.progressBar.setValue(0)
            return
        self.error_array = np.array(self.error_array)       # convert to np array to do reshape
        self.error_mat = self.error_array.reshape(len(self.mat_xaxis), len(self.mat_yaxis))
        self.error_matrix = np.flip(self.error_mat, 0) # need flipping to make start index of the error array at the origin of imshow graph
        self.error.clear()
        ax = self.error.add_subplot(111)
        map = ax.imshow(self.error_matrix, cmap='inferno', extent=[min(self.mat_xaxis), max(self.mat_xaxis), min(self.mat_yaxis), max(self.mat_yaxis)],  aspect='auto')    # self.mat_xaxis, self.mat_yaxis,
        divider3 = make_axes_locatable(ax)
        cax3 = divider3.append_axes("right", size="10%")
        plt.colorbar(map, cax=cax3)
        self.canvas_error_map.draw()
        self.progressBar.setValue(100)
        self.start_button.setText("Start")

    def calculate_errormap(self):
        for x_val in self.mat_xaxis:
            for y_val in self.mat_yaxis:
                if self.exit_error_map:    # end loading of error map
                    return False
                print("in calculate errormap\n")
                error_val = self.fill_error_mat(x_val, y_val)
                self.error_array.append(error_val)
                self.complete += self.step
                self.progressBar.setValue(self.complete)
        return True        

    x_chunks_arr = []
    y_chunks_arr = []
    
    
    def fill_error_mat(self, x_val, y_val):
        if isinstance(self.overlap_vals, float):      # when overlap is const
            overlap_val = self.overlap_vals
            if self.x_param == "number of chunks":
                num_chunk = x_val
                order = y_val
            if self.y_param == "number of chunks":
                num_chunk = y_val
                order = x_val
        else:
            if self.x_param == "overlapping":  # then y is chunk or order
                overlap_val = x_val
                if isinstance(self.order_vals, int):  # then order is const
                    order = self.order_vals
                    num_chunk = y_val
                else:  # then chunks is const
                    num_chunk = self.chunk_vals
                    order = y_val
            else:  # then x is chunk or order
                overlap_val = y_val
                if isinstance(self.order_vals, int):  # then order is const
                    order = self.order_vals
                    num_chunk = x_val
                else:  # then chunks is const
                    num_chunk = self.chunk_vals
                    order = x_val
                

        total_error = []
        chunk_width = ((max(self.x_data) - min(self.x_data)) / num_chunk) + (overlap_val / 2)
        self.chunk_index.clear()
        self.chunk_index.append(min(self.x_data))
        for i in range(0, (2 * num_chunk)):
            if (i % 2) == 0:
                ind = self.chunk_index[i] + chunk_width
            else:
                ind = self.chunk_index[i] - overlap_val
            self.chunk_index.append(ind)
        self.chunk_index.pop()
        self.new_chunk = 0
        for j in range(0, len(self.chunk_index)-1):
            self.x_each_chunk.clear()
            self.y_each_chunk.clear()
            if (j % 2) == 0:
                for k in range(self.new_chunk, len(self.x_data)):
                    if self.x_data[k] >= self.chunk_index[j] and self.x_data[k] <= self.chunk_index[j+1]:
                        self.x_each_chunk.append(self.x_data[k])
                        self.y_each_chunk.append(self.y_data[k])
                        if (j+2) != len(self.chunk_index):
                            if self.x_data[k] <= self.chunk_index[j+2]:
                                self.new_chunk = k+1
                    else:
                        break
                self.x_chunks_arr.append(self.x_each_chunk)
                self.y_chunks_arr.append(self.y_each_chunk)

        for i in range(0, len(self.x_chunks_arr)):
            coeffs = poly.polyfit(self.x_chunks_arr[i], self.y_chunks_arr[i], order)
            x_fitline = np.linspace(self.x_chunks_arr[0], self.x_chunks_arr[-1], num=len(self.x_chunks_arr) * 1)
            y_fitline = poly.polyval(x_fitline, coeffs)

        k = 0
        for i in range(len(self.y_chunks_arr)):
            for j in range(len(self.y_chunks_arr[i])):
                total_error.append(((self.y_chunks_arr[i][j] * y_fitline[k]) / self.y_chunks_arr[i][j]) * 100)
                ++k
        return np.median(total_error)
    
def main():
    app = QtWidgets.QApplication(sys.argv)
    application = MyWidget()
    application.show()
    app.exec_()

if __name__ == "__main__":
    main()