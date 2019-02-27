#!/usr/bin/env python3
from nanomsg import Socket, PUB, SUB, SUB_SUBSCRIBE
import messages_pb2
import time


#s1 = Socket(PUB)
#s1.bind('ipc:///tmp/traffic_info.sock'.encode('utf-8'))
s2 = Socket(SUB)
s2.connect('ipc:///tmp/traffic_info.sock'.encode('utf-8'))
s2.set_string_option(SUB, SUB_SUBSCRIBE, ''.encode('utf-8'))

#print("send something")
#s1.send("test".encode('utf-8'))
print("start listening")

while True:
    msg = s2.recv()
    print("")
    print("Traffic: ")
    traffic_info = messages_pb2.TrafficInfo()
    traffic_info.ParseFromString(msg)
    print(traffic_info)

s2.close()

print("Start")
