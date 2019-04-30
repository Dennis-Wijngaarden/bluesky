"""Plots the SSD in a external matplotlib plot"""
# Import the gloval bluesky objects
from bluesky import stack  #, settings, navdb, traf, sim, scr, tools
from bluesky import traf
from bluesky import sim
import bluesky as bluesky

import threading
import matplotlib.pyplot as plt
import numpy as np

ssd_plotter = None

# Initialisation of plugin
def init_plugin():
    global ssd_plotter
    ssd_plotter = SSD_plotter()

    config = {
        'plugin_name':      'PLOT_SSD',
        'plugin_type':      'sim',
        'update_interval':  1.0,
        'update':           ssd_plotter.update,
    }

    stackfunctions = {
        "PLOTSSD" : [
            'PLOTSSD acid, ON/OFF',
            'acid, onoff',
            ssd_plotter.enable_plot,
            'SSD PLOTTER'
        ]
    }

    return config, stackfunctions

class SSD_plotter(object):
    def __init__(self):
        self.thread = None
        self.plot_elements = {}

    def enable_plot(self, idx, state):
        # Check active state of plot and activate/deactivate
        if traf.SSD_plot_on[idx]:
            if state == True:
                return False, 'SSD plot already enabled'
            else:
                self.plot_elements[traf.id[idx]].stop()
                traf.SSD_plot_on[idx] = False
                return True, 'SSD plot disabled'

        if not traf.SSD_plot_on[idx]:
            if state == True:
                self.plot_elements[traf.id[idx]] = SSD_plot_element(traf.id[idx])
                traf.SSD_plot_on[idx] = True
                return True, 'SSD plot enabled'
            else:
                return False, 'SSD plot not enabled'
    
    def update(self):
        # Check if aircaft which needs to be plotted still exists
        for key in self.plot_elements:
            try:
                key in traf.id
            except:
                self.plot_elements[key].stop()
                self.plot_elements.pop(key)
                    
        # Check if ARV and FRV exist
        try:
            traf.asas.ARV_calc
        except:
            return

        # Update ARV and FRV of plot threads
        for key in self.plot_elements:
            self.plot_elements[key].ARV = traf.asas.ARV[traf.id2idx(key)]
            self.plot_elements[key].FRV = traf.asas.FRV[traf.id2idx(key)]
            self.plot_elements[key].update(traf.id2idx(key))

class SSD_plot_element(object):
    def __init__(self, identity):

        self.ARV = None
        self.FRV = None

        # Plot functions
        plt.ion()
        self.fig, self.ax = plt.subplots()
        self.fig.canvas.set_window_title('SSD ' + identity)
        self.fig.canvas.draw()
        #self.ax.figure.canvas.mpl_connect('close_event', self.stop())

    def update(self, idx):
        # First clear previous data
        self.ax.clear()

        v_own = np.array([traf.gseast[idx], traf.gsnorth[idx]])
        vmin = traf.Vmin[idx]
        vmax = traf.Vmax[idx]
        N_angle = 180

        angles = np.arange(0, 2 * np.pi, 2 * np.pi / N_angle)
        xyc = np.transpose(np.reshape(np.concatenate((np.sin(angles), np.cos(angles))), (2, N_angle)))
        SSD_lst = [list(map(list, np.flipud(xyc * vmax))), list(map(list , xyc * vmin ))]

        #Outer circle: vmax
        SSD_outer = np.array(SSD_lst[0])
        x_SSD_outer = np.append(SSD_outer[:,0] , np.array(SSD_outer[0,0]))
        y_SSD_outer = np.append(SSD_outer[:,1] , np.array(SSD_outer[0,1]))
        #Inner circle: vmin
        SSD_inner = np.array(SSD_lst[1])
        x_SSD_inner = np.append(SSD_inner[:,0] , np.array(SSD_inner[0,0]))
        y_SSD_inner = np.append(SSD_inner[:,1] , np.array(SSD_inner[0,1]))

        self.ax.plot(x_SSD_outer, y_SSD_outer, color = '#404040')
        self.ax.plot(x_SSD_inner, y_SSD_inner, color = '#404040')

        if self.ARV == None:
            return

        # Plot Conflict Zones
        for i in range(len(self.FRV)):
            FRV_element = np.array(self.FRV[i])
            x_FRV = np.append(FRV_element[:,0], np.array(FRV_element[0,0]))
            y_FRV = np.append(FRV_element[:,1], np.array(FRV_element[0,1]))
            self.ax.fill(x_FRV, y_FRV, color = '#FF0000') #red

        # Plot Conflict Free Zones
        for i in range(len(self.ARV)):
            ARV_element = np.array(self.ARV[i])
            x_ARV = np.append(ARV_element[:,0], np.array(ARV_element[0,0]))
            y_ARV = np.append(ARV_element[:,1], np.array(ARV_element[0,1]))
            self.ax.fill(x_ARV, y_ARV, color = '#C0C0C0') #grey

        # make inner circle red
        self.ax.fill(x_SSD_inner, y_SSD_inner, color = '#ffffff') #white

        self.ax.axis('equal')
        self.fig.canvas.draw()

    def stop(self):
        plt.close()