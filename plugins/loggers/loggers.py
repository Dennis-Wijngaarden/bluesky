""" Bluesky plugin in order to log"""

import numpy as np
from bluesky import traf
from bluesky.tools import datalog, areafilter
from bluesky.tools.aero import nm, kts
from bluesky import settings

logger = None

header_fl_log = \
    '#######################################################\n' + \
    'FLST LOG\n' + \
    'Flight Statistics\n' + \
    '#######################################################\n\n' + \
    'Parameters [Units]:\n' + \
    'Simulation time [s], ' + \
    'Callsign [-], ' + \
    'Lat [deg], ' + \
    'Lon [deg], ' + \
    'Distance Flown [m], ' + \
    'Heading [deg], ' + \
    'Track [deg], ' + \
    'True Airspeed [m/s], ' + \
    'Groundspeed [m/s], ' + \
    'Active waypoint lat [deg], ' + \
    'Active waypoint lon [deg]'

header_conf_log = \
    '#######################################################\n' + \
    'CONF LOG\n' + \
    'Conflict parameters\n' + \
    '#######################################################\n\n' + \
    'Parameters [Units]:\n' + \
    'Simulation time [s], ' + \
    'Ownship [-], ' + \
    'Intruder [-], ' + \
    'qdr [deg], ' + \
    'dist [m], ' + \
    'dcpa [m], ' + \
    'tcpa [s], ' + \
    'tLOS [s]'

header_gf_log = \
    '#######################################################\n' + \
    'GEOFENCE LOG\n' + \
    'GEOFENCE variables\n' + \
    '#######################################################\n\n' + \
    'Parameters [Units]:\n' + \
    'Simulation time [s], ' + \
    'Callsign [-], ' + \
    'In GF [-]'
    


### Initialization function of your plugin. Do not change the name of this
### function, as it is the way BlueSky recognises this file as a plugin.
def init_plugin():

    global logger
    logger = Logger()

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'LOGGERS',

        # The type of this plugin.
        'plugin_type':     'sim',

        # Plugin update interval in seconds.
        'update_interval': 1.0,

        # Update at the end of a timestep.
        'update':          logger.update,

        # Update at the beginning of a timestep.
        #'preupdate':       preupdate
        }

    stackfunctions = {
        # The command name for your function
        'INIT_LOGGERS': [
            # A short usage string. This will be printed if you type HELP <name> in the BlueSky console
            'INIT_LOGGERS TS RS',

            # A list of the argument types your function accepts.
            'int, int',

            # The name of your function in this plugin
            logger.initialize,

            # a longer help text of your function.
            'Print something to the bluesky console based on the flag passed to MYFUN.'],

        # The command name for your function
        'STOP_LOGGERS': [
            # A short usage string. This will be printed if you type HELP <name> in the BlueSky console
            'STOP_LOGGERS',

            # A list of the argument types your function accepts.
            '',

            # The name of your function in this plugin
            logger.stop,

            # a longer help text of your function.
            'Print something to the bluesky console based on the flag passed to MYFUN.']
    }

    return config, stackfunctions

class Logger():
    def __init__(self):
        self.gfl_active = False # Set geofence log non-active
        self.confl_active = False # Set conflict log non-active
        self.TS = None # Test series number
        self.RS = None # Ruleset number
        self.TR = None # Test run/scenario number
        self.initialized = False # Flag if log has been initialized
    
    def update(self):
        # Only update if initialized
        if self.initialized:
            if self.vehicle_at_end_route():
                self.stop()
                return
            self.log_fl()
            self.log_conf()
            try:
                areafilter.areas['GF']
            except:
                pass
            else:
                self.log_gf()


    def initialize(self, TS, RS):
        self.TS = TS
        self.RS = RS
        self.initialized = True

        self.conf_log = datalog.crelog("CONFLOG_TS" + str(TS) + "_RS" + str(RS), None, header_conf_log)
        self.fl_log = datalog.crelog("FLLOG_TS" + str(TS) + "_RS" + str(RS), None, header_fl_log)
        self.gf_log = datalog.crelog("GFLOG_TS" + str(TS) + "_RS" + str(RS), None, header_gf_log)

        self.start()
    
    def start(self):
        self.conf_log.start()
        self.fl_log.start()
        self.gf_log.start()

        self.update()

    def stop(self):
        self.fl_log.reset()
        self.conf_log.reset()
        self.gf_log.reset()
        self.initialized = False

    def log_conf(self):
        for i in range(len(traf.cd.confpairs)):
            vehicle_id = traf.cd.confpairs[i][0]
            vehicle_idx = traf.id2idx(vehicle_id)

            intruder_id = traf.cd.confpairs[i][1]
            intruder_idx = traf.id2idx(intruder_id)

            self.conf_log.log(  traf.id[vehicle_idx],
                                traf.id[intruder_idx],
                                traf.cd.qdr[i],
                                traf.cd.dist[i],
                                traf.cd.dcpa[i],
                                traf.cd.tcpa[i],
                                traf.cd.tLOS[i])
        
    def log_fl(self):
        self.fl_log.log(    traf.id,
                            traf.lat,
                            traf.lon,
                            traf.distflown ,
                            traf.hdg,
                            traf.trk,
                            traf.tas,
                            traf.gs,
                            traf.actwp.lat,
                            traf.actwp.lon)

    def log_gf(self):
        for i in range(traf.ntraf):
            in_gf = areafilter.checkInside('GF', traf.lat[i], traf.lon[i], traf.alt[i])
            self.gf_log.log(    traf.id[i],
                                in_gf)
    
    def vehicle_at_end_route(self):
        for i in range(traf.ntraf):
            if traf.ap.route[i].iactwp < 0:
                return True
        return False