""" BlueSky plugin template. The text you put here will be visible
    in BlueSky as the description of your plugin. """
# Import the global bluesky objects. Uncomment the ones you need
from bluesky import stack, traf  #, settings, navdb, traf, sim, scr, tools
from bluesky.tools.aero import ft, kts
from bluesky.tools import areafilter

import numpy as np
import zmq
import utm.common_utm_pb2 as messages_pb2
import random

### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
traffic_sim = None

def init_plugin():
    global traffic_sim
    traffic_sim = TrafficSim()
    random.seed()

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'UTM_TRAF_SIM',

        # The type of this plugin.
        'plugin_type':     'sim',

        # Plugin update interval in seconds.
        'update_interval': 1.0,

        # Update at the end of a timestep.
        'update':          traffic_sim.update,
        }

    stackfunctions = {
            # The command name for your function
            'START_TRAFSIM': [
                # A short usage string. This will be printed if you type HELP <name> in the BlueSky console
                'START_TRAFSIM lat,lon,lat,lon,n_traf',

                # A list of the argument types your function accepts.
                'latlon, latlon, int',

                # The name of your function in this plugin
                traffic_sim.start,

                # a longer help text of your function.
                'TRAFFIC_SIM'],

            'STOP_TRAFSIM': [
                # A short usage string. This will be printed if you type HELP <name> in the BlueSky console
                'STOP_TRAFSIM',

                # A list of the argument types your function accepts.
                '',

                # The name of your function in this plugin
                traffic_sim.stop,

                # a longer help text of your function.
                'TRAFFIC_SIM'],

        }
    return config, stackfunctions

class TrafficSim(object):
    def __init__(self):
        self.enabled = False

        self.n_traf = 0
        self.lat_top = 0
        self.lon_left = 0
        self.lat_bot = 0
        self.lon_right = 0

        self.v_min = 5.
        self.v_max = 30.
        self.alt_default = 300. #ft
    
    def start(self, lat_top, lon_left, lat_bot, lon_right, n_traf):
        self.enabled = True
        self.lat_top = lat_top
        self.lon_left = lon_left
        self.lat_bot = lat_bot
        self.lon_right = lon_right
        self.n_traf = n_traf
        stack.stack("BOX TRAF," + str(self.lat_top) + "," + str(self.lon_left) + "," + str(self.lat_bot) + "," + str(self.lon_right))

        # create n random vehicles within area
        for i in range(self.n_traf):
            lat = random.uniform(lat_bot, lat_top)
            lon = random.uniform(lon_left, lon_right)
            hdg = random.uniform(0., 360.)
            spd = random.uniform(self.v_min, self.v_max) / kts
            stack.stack("CRE SIM" + str(i) + " DRONE," + str(lat) + "," + str(lon) + "," + str(hdg) + "," + str(self.alt_default) + ","  + str(spd))
    
    def stop(self):
        if self.enabled:
            for i in range(self.n_traf):
                stack.stack("DEL SIM" + str(i))
            stack.stack("DEL TRAF")
            self.enabled = False

    def update(self):
        for i in range(traf.ntraf):
            if "SIM" in traf.id[i]:
                # Get state variables from vehicle
                callsign = traf.id[i]
                lat = traf.lat[i]
                lon = traf.lon[i]
                alt = traf.alt[i]

                # check in area
                inside = areafilter.checkInside("TRAF", lat, lon, alt)

                # Move aircraft if not inside
                if not inside:
                    # Check which vertex of box is violated and move to other side
                    if lat > self.lat_top:
                        traf.move(i, self.lat_bot, lon)
                    elif lat < self.lat_bot:
                        traf.move(i, self.lat_top, lon)
                    elif lon > self.lon_right:
                        traf.move(i, lat, self.lon_left)
                    elif lon < self.lon_left:
                        traf.move(i, lat, self.lon_right)