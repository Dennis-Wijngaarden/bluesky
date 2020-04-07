import os
import traffic_generator
import scenario_generator

traffic_generator.generate_traffic()
scenario_generator.generate_scenario()
os.system("python thesis_tools/scenario_generation/wind_generator.py")
os.system("python thesis_tools/scenario_generation/route_generator.py")
os.system("python thesis_tools/scenario_generation/geofence_generator.py")
