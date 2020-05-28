import os
import traffic_generator
import scenario_generator
import wind_generator
import route_generator
import geofence_generator

print('Generate traffic')
traffic_generator.generate_traffic()

#print('Generate wind')
#wind_generator.generate_wind()

print('Generate scenario')
scenario_generator.generate_scenario()

print('Generate route')
route_generator.generate_route()

print('Generate geofences')
geofence_generator.generate_geofences()
