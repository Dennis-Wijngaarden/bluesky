import os
import traffic_generator
import scenario_generator_v2
import route_generator_v2
import geofence_generator_v2

print('Generate traffic')
traffic_generator.generate_traffic()

#print('Generate wind')
#wind_generator.generate_wind()

print('Generate scenario')
scenario_generator_v2.generate_scenario()

print('Generate route')
route_generator_v2.generate_route()

print('Generate geofences')
geofence_generator_v2.generate_geofence()
