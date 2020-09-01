""" BlueSky plugin template. The text you put here will be visible
    in BlueSky as the description of your plugin. """
# Import the global bluesky objects. Uncomment the ones you need
from bluesky import stack  #, settings, navdb, traf, sim, scr, tools
from bluesky.tools.aero import ft, kts

import numpy as np
import threading
import zmq
import utm.common_utm_pb2 as messages_pb2

state_receiver = None

### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():
    global state_receiver
    state_receiver = StateReceiver()

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'UTM_STATE_RECEIVER',

        # The type of this plugin.
        'plugin_type':     'sim',

        # Plugin update interval in seconds.
        'update_interval': 2.5,

        # Update at the end of a timestep.
        'update':          update,
        }

    stackfunctions = {
            # The command name for your function
            'OWNSTATERECEIVE': [
                # A short usage string. This will be printed if you type HELP <name> in the BlueSky console
                'OWNSTATERECEIVE ON/OFF',

                # A list of the argument types your function accepts.
                '[onoff]',

                # The name of your function in this plugin
                state_receiver.enable_disable,

                # a longer help text of your function.
                'OWNSTATE_STREAM']
        }
    return config, stackfunctions

class StateReceiver(object):
    def __init__(self):
        self.thread = None
        self.enabled = False

    def enable_disable(self, enable):
        if self.enabled:
            # if it needs to be disabled
            if not enable:
                self.thread.stop()
                del self.thread
                self.thread = None
                self.enabled = False
                return True, "Disabled"
            else:
                return False, "Already enabled"
        else:
            # it needs to be enabled
            if enable:
                self.thread = StateReceiverThread()
                self.enabled = True
                return True, "Enabled"
            else:
                return False, "Already disabled"

class StateReceiverThread(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
        self.lock = threading.Lock()
        self.stop_event = threading.Event()

        self.spawned = False # Variable to tell if vehicle has already spawned for the first time

        # init socket
        self.context = zmq.Context()
        self.socket = self.context.socket(zmq.SUB)
        self.socket.connect("tcp://127.0.0.1:4243")
        self.socket.setsockopt(zmq.SUBSCRIBE, b'')
        self.state_info_msg = messages_pb2.OwnState()

        self.start()
    
    def run(self):
        while not self.stop_event.is_set():
            try:
                msg = self.socket.recv(flags=zmq.NOBLOCK)
            except zmq.Again:
                continue
            self.state_info_msg.ParseFromString(msg)

            lat = self.state_info_msg.lat / 10000000.
            lon = self.state_info_msg.lon / 10000000.
            alt = self.state_info_msg.alt / 1000. / ft
            v_n = self.state_info_msg.v_n / 1000.
            v_e = self.state_info_msg.v_e / 1000.
            gs = np.sqrt(v_n**2 + v_e**2) / kts
            trk = np.rad2deg(np.arctan2(v_e, v_n))

            if self.spawned:
                stack.stack("MOVE OWN," + str(lat) + "," + str(lon) + "," + str(alt) + "," + str(trk) + "," + str(gs))
            else:
                stack.stack("CRE OWN,DRONE," + str(lat) + "," + str(lon) + "," + str(trk) + "," + str(alt) + "," + str(gs))
                self.spawned = True
            
    def stop(self):
        stack.stack("DEL OWN")
        self.stop_event.set()
        self.socket.close()

def update():
    return