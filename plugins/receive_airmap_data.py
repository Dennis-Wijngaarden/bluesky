"""Receives airmap postion data of nearby aircraft"""
# Import the global bluesky objects
from bluesky import stack  #, settings, navdb, traf, sim, scr, tools
from bluesky import traf
from bluesky import sim
from bluesky.tools.aero import ft
import bluesky as bs

import threading

from nanomsg import Socket, PUB, SUB, SUB_SUBSCRIBE
import airmap_lib.messages_pb2 as messages_pb2
import time

airmap_receiver = None

# Initialisation of plugin
def init_plugin():
    global airmap_receiver
    airmap_receiver = Airmap_receiver()

    config = {
        'plugin_name':      'AIRMAP_RECEIVER',
        'plugin_type':      'sim',
        'update_interval':  1.0,
        'update':           airmap_receiver.update,
    }

    stackfunctions = {
        'AIRMAPRECEIVE': [
            'AIRMAPRECEIVER ON/OFF',
            '[onoff]',
            airmap_receiver.enable,
            'AIRMAP STREAM']
    }

    return config, stackfunctions


class Airmap_receiver(object):
    def __init__(self):
        self.thread = None
        self.enabled = False
        self.traffic = {}

    def update(self):
        time_start = time.time()
        if self.enabled:
            try:
                self.thread.lock.acquire()
                for i in range(len(self.thread.received_traffic_objects)):
                    # Store new or existing id's
                    if (self.thread.received_traffic_objects[i].aircraft_id in self.traffic.keys()):
                        stack.stack("MOVE " +   self.thread.received_traffic_objects[i].aircraft_id + " " + \
                                                str(self.thread.received_traffic_objects[i].lat) + " " + \
                                                str(self.thread.received_traffic_objects[i].lon) + " " + \
                                                str(self.thread.received_traffic_objects[i].alt) + " " + \
                                                str(self.thread.received_traffic_objects[i].heading) + " " + \
                                                str(self.thread.received_traffic_objects[i].groundspeed))
                    else:
                        stack.stack("CRE " +    self.thread.received_traffic_objects[i].aircraft_id + " " + \
                                                "b738 " + \
                                                str(self.thread.received_traffic_objects[i].lat) + " " + \
                                                str(self.thread.received_traffic_objects[i].lon) + " " + \
                                                str(self.thread.received_traffic_objects[i].heading) + " " + \
                                                str(self.thread.received_traffic_objects[i].alt) + " " + \
                                                str(self.thread.received_traffic_objects[i].groundspeed))
                    #print("traffic delay: ", time_start - self.thread.received_traffic_objects[i].recorded_time)
                    self.traffic[self.thread.received_traffic_objects[i].aircraft_id] = self.thread.received_traffic_objects[i]
                
                # cleanup
                pop_keys = []
                for key in self.traffic:
                    if (time.time() - self.traffic[key].time_received) > 5.:
                        pop_keys.append(key)

                for i in range(len(pop_keys)):
                    stack.stack("DEL " + self.traffic[pop_keys[i]].aircraft_id)
                    self.traffic.pop(pop_keys[i])

                # empty buffer
                self.thread.received_traffic_objects = []
            finally:
                self.thread.lock.release()

    def enable(self, flag=True):
        prev_flag = self.enabled
        self.enabled = flag
        # Switch on socket if status turned from off to on
        if (prev_flag == False and flag == True):
            self.thread = Airmap_receiver_thread()
        # Switch off if previously switched on
        if (prev_flag == True and flag == False):
            self.thread.stop()
            self.thread = None

class Airmap_receiver_thread(threading.Thread):
    def __init__(self):
        threading.Thread.__init__(self)
        self.lock = threading.Lock()
        self.stop_event = threading.Event()
        
        self.s1 = Socket(SUB)
        self.s1.connect('ipc:///tmp/traffic_info_0'.encode('utf-8'))
        self.s1.set_string_option(SUB, SUB_SUBSCRIBE, ''.encode('utf-8'))
        self.traffic_info = messages_pb2.TrafficInfo()

        self.received_traffic_objects = []

        self.start()

    def run(self):
        while not self.stop_event.is_set():
            try:
                msg = self.s1.recv()
            except self.s1.recv_timeout():
                continue
            self.traffic_info.ParseFromString(msg)

            try:
                self.lock.acquire()
                
                self.received_traffic_objects.append(ReceivedTraffic(   self.traffic_info.aircraft_id, \
                                                                        self.traffic_info.recorded_time, \
                                                                        self.traffic_info.lat * 10**-7, \
                                                                        self.traffic_info.lon * 10**-7, \
                                                                        self.traffic_info.alt / 1000., \
                                                                        self.traffic_info.groundspeed / 1000., \
                                                                        self.traffic_info.heading))
            finally:
                self.lock.release()

    def stop(self):
        self.stop_event.set()
        self.s1.close()

class ReceivedTraffic(object):
    def __init__(self, aircraft_id, recoreded_time, lat, lon, alt, groundspeed, heading):
        self.aircraft_id    = aircraft_id
        self.recorded_time  = recoreded_time
        self.lat            = lat
        self.lon            = lon
        self.alt            = alt
        self.groundspeed    = groundspeed
        self.heading        = heading
        self.time_received  = time.time()

        #print(aircraft_id, ' ', recoreded_time, ' ', lat, ' ', lon, ' ', alt, ' ', groundspeed, ' ', heading, ' ', self.time_received)
