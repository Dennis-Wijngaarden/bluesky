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
        'AIRMAPSTREAMON': [
            'AIRMAPSTREAMON acid, instance',
            'acid,int',
            airmap_streamer.enable,
            'Switch On Airmap Stream'],
        'AIRMAPSTREAMOFF': [
            'AIRMAPSTREAMOFF acid',
            'acid',
            airmap_streamer.disable,
            'Switch Off Airmap Stream']
    }

    return config, stackfunctions

class Airmap_streamer(object):
    def __init__(self):
        self.sockets = {} # acid : Airmap_stream_socket

    def update(self):
        for key in self.sockets.keys():
            self.sockets[key].update()

    def enable(self, idx, instance):
        acid = traf.id[idx]
        
        # Check if acid already has a open socket
        try:
            self.sockets[acid]
            return False, 'Airmap stream socket already open for this aircraft'
        except:
            self.sockets[acid] = Airmap_stream_socket(acid, instance)
            return True, 'Airmap stream socket opened for aircraft'
  
    def disable(self, idx):
        acid = traf.id[idx]

        try:
            self.sockets[acid].close()
            del self.sockets[acid]
            return True, 'Airmap stream socket succesfully deleted'  
        except:
            return False, 'Airmap stream socket doesnt exist for this aircraft'
    

class Airmap_stream_socket(object):
    def __init__(self, acid, instance):
        self.socket = Socket(PUB)
        self.gps_message = messages_pb2.GPSMessage()
        self.acid = acid
        self.instance = instance

        # Bind socket
        self.socket.bind(('ipc:///tmp/gps_position_' + str(instance)).encode('utf-8'))
    
    def update(self):
        idx = traf.id2idx(self.acid)
        self.gps_message.timestamp = int(time.time() * 1000.)
        self.gps_message.lat = int(bs.traf.lat[idx] * 10.**7)
        self.gps_message.lon = int(bs.traf.lon[idx] * 10. ** 7) 
        self.gps_message.vn = int(bs.traf.gsnorth[idx] * 1000.)
        self.gps_message.ve = int(bs.traf.gseast[idx] * 1000.)
        self.gps_message.vd = int(bs.traf.vs[idx] * -1000.)
        self.gps_message.alt_msl = int(bs.traf.alt[idx] * ft * 1000.) 
        self.gps_message.alt_agl = int(bs.traf.alt[idx] * ft * 1000.) 
        self.socket.send(self.gps_message.SerializeToString())

    def close(self):
        self.socket.close()