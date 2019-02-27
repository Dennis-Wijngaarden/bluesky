""" Streams airmap postion data of first aircraft of traffic array"""
# Import the global bluesky objects
from bluesky import stack  #, settings, navdb, traf, sim, scr, tools
from bluesky import traf
from bluesky import sim
from bluesky.tools.aero import ft
import bluesky as bs

from nanomsg import Socket, PUB
import airmap_lib.messages_pb2 as messages_pb2
import time

airmap_streamer = None

# Initialisation of plugin
def init_plugin():
    global airmap_streamer
    airmap_streamer = Airmap_streamer()

    config = {
        'plugin_name':      'AIRMAP_STREAMER',
        'plugin_type':      'sim',
        'update_interval':  1.0,
        'update':           airmap_streamer.update,
    }

    stackfunctions = {
        'AIRMAPSTREAM': [
            'AIRMAPSTREAM ON/OFF',
            '[onoff]',
            airmap_streamer.enable,
            'AIRMAP STREAM']
    }

    return config, stackfunctions

class Airmap_streamer(object):
    def __init__(self):
        self.enabled = False
        self.s1 = Socket(PUB)
        self.gps_message = messages_pb2.GPSMessage()

    def update(self):
        if self.enabled:
            self.gps_message.timestamp = 0 #bs.sim.simt
            self.gps_message.lat = int(bs.traf.lat[0] * 10. ** 7) #521234157
            self.gps_message.lon = int(bs.traf.lon[0] * 10. ** 7) #41122334
            self.gps_message.alt_msl = int(bs.traf.alt[0] * ft * 1000.)  #10
            self.gps_message.alt_agl = int(bs.traf.alt[0] * ft * 1000.)  #10
            self.s1.send(self.gps_message.SerializeToString())

    def enable(self, flag=True):
        prev_flag = self.enabled
        self.enabled = flag 
        # Enable socket when switched on
        if (prev_flag == False and flag == True):
            self.s1.bind('ipc:///tmp/gps_position.sock'.encode('utf-8'))
        # Swirch off airmap streamer
        if (prev_flag == True and flag == False):
            self.s1.close()