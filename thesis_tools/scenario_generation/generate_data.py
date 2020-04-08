import os
import traffic_generator
import scenario_generator
import wind_generator
import route_generator
import geofence_generator

traffic_generator.generate_traffic()
wind_generator.generate_wind()
scenario_generator.generate_scenario()
route_generator.generate_route()
geofence_generator.generate_geofences()
