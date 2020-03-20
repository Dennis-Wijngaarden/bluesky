import parameters
import random
import os

# Initialize randomizer
random.seed()

# Scenario folder
location = "scenario/Thesis/TS1"

# First empty folder with scenarios
filelist = [ f for f in os.listdir(location) if f.endswith(".scn") ]
for f in filelist:
    os.remove(os.path.join(location, f))

batch_file = open(location + "/batch.scn", "w")

for i in range(parameters.N_missions):
    scenario_call = "Thesis/TS1/test" + str(i) + ".scn"
    scenario_file = open(location + "/test" + str(i) + ".scn", "w")
    # write scenario lines
    hdg = random.uniform(0, 360)
    spd = random.uniform(10., 20.)
    scenario_file.write("00:00:00.00>CRE UAV0 UAV_" + str(i) + "_0 " + parameters.ref_position + " " + str(hdg) + " " +\
            str(parameters.ref_alt) + " " + str(spd) + "\n")
    d_psi = random.uniform(0, 360) # [deg]
    cpa = random.uniform(0, 50.) / 1852. # [nm]
    spd = random.uniform(10., 20.) / 1.94384449
    scenario_file.write("00:00:00.00>CRECONFS UAV1 UAV_" + str(i) + "_1 UAV0 " + str(d_psi) + " " + str(cpa) + " " + str(25. + parameters.t_extra) + " 0 0 " + str(spd) + "\n")
    scenario_file.close()
    for j in range(parameters.N_RS):
        # write batch lines
        batch_file.write("00:00:00.00>SCEN test_" + str(i) + "_RS" + str(j + 1) + "\n")
        batch_file.write("00:00:00.00>ASAS ON\n")
        batch_file.write("00:00:00.00>RESO SSDUAV\n")
        batch_file.write("00:00:00.00>PRIORULES ON RS" + str(j + 1) + "\n")
        logname = "TS1_RS" + str(j + 1) + "_" + str(i)
        batch_file.write("00:00:00.00>CRELOG " + logname + " 0.05\n")
        batch_file.write("00:00:00.00>" + logname + " ADD FROM traf id, lat, lon\n")
        batch_file.write("00:00:00.00>" + logname + " ON\n")
        batch_file.write("00:00:00.00>PCALL " + scenario_call + "\n")
        batch_file.write("00:00:00.00>FF\n")
        batch_file.write("00:00:00.00>SCHEDULE 00:05:00 " + logname + " OFF\n")
        batch_file.write("00:00:00.00>SCHEDULE 00:05:00 HOLD\n")

batch_file.close()