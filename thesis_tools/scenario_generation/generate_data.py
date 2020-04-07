import os
import traffic_generator

traffic_generator.generate_traffic()
os.system("python thesis_tools/scenario_generation/scenario_generator.py")
os.system("python thesis_tools/scenario_generation/wind_generator.py")
os.system("python thesis_tools/scenario_generation/route_generator.py")
os.system("python thesis_tools/scenario_generation/geofence_generator.py")
