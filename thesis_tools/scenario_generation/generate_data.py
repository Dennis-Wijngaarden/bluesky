import os
import traffic_generator
import scenario_generator
import wind_generator
import route_generator

traffic_generator.generate_traffic()
scenario_generator.generate_scenario()
wind_generator.generate_wind()
route_generator.generate_route()
os.system("python thesis_tools/scenario_generation/geofence_generator.py")
