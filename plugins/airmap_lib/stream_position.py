#!/usr/bin/env python3
from nanomsg import Socket, PUB
import messages_pb2
import time

s1 = Socket(PUB)
s1.bind('ipc:///tmp/gps_position.sock'.encode('utf-8'))

gps_message = messages_pb2.GPSMessage()
gps_message.timestamp = 0
gps_message.lat = 520000000
gps_message.lon = 40000000
gps_message.alt_msl = 10
gps_message.alt_agl = 10

while True:
    s1.send(gps_message.SerializeToString())
    time.sleep(1)

s1.close()

print("Start")
