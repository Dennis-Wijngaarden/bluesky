#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 08:37:28 2019

@author: dennis
"""
from bluesky.tools import geo
from bluesky.tools.aero import nm
import numpy as np
# Try to import pyclipper
try:
    import pyclipper
except ImportError:
    print("Could not import pyclipper, RESO SSDUAV will not function")

def loaded_pyclipper():
    """ Return true if pyclipper is successfully loaded """
    import sys
    return "pyclipper" in sys.modules

# Start placeholder
def start(asas):
    pass

def detect(asas, traf):
    """ Detect all current conflicts """
    # Construct the SSD
    ConstructSSD(asas, traf)

def resolve(asas, traf):
    """ Resolve all current conflicts """

    # Check if ASAS is ON first!
    if not asas.swasas:
        return

    InitializeSSD(asas, traf.ntraf)

    # Construct the SSD
    ConstructSSD(asas, traf)

    # Get resolved speed-vector
    calculate_resolution(asas, traf)

    # Now assign resolutions to variables in the ASAS class
    # Start with current states, need a copy, otherwise it changes traf!
    asas.trk = np.copy(traf.hdg)
    asas.tas = np.copy(traf.gs)
    # Calculate new track and speed
    # No need to cap the speeds, since SSD implicitly caps
    new_trk  = np.arctan2(asas.asase, asas.asasn) * 180 / np.pi
    new_tas  = np.sqrt(asas.asase ** 2 + asas.asasn ** 2)

    # Sometimes an aircraft is in conflict, but no solutions could be found
    # In that case it is assigned 0 by ASAS, but needs to handled
    asas_cmd = np.logical_and(asas.inconf, new_tas > 0)

    # Assign new track and speed for those that are in conflict
    asas.trk[asas_cmd] = new_trk[asas_cmd]
    asas.tas[asas_cmd] = new_tas[asas_cmd]
    # Not needed as it is a 2D-implementation...
    asas.vs   = traf.vs

def InitializeSSD(asas, ntraf):
    """ Initialize variables for SSD """
    # Need to do it here, since ASAS.reset doesn't know ntraf
    asas.FRV          = [None] * ntraf
    asas.ARV          = [None] * ntraf
    # For calculation purposes
    asas.ARV_calc     = [None] * ntraf
    asas.inrange      = [None] * ntraf
    #asas.inconf       = np.zeros(ntraf, dtype=bool) Commented out as it overwrites already generated asas.inconf
    # Stores resolution vector, also used in visualization
    asas.asasn        = np.zeros(ntraf, dtype=np.float32)
    asas.asase        = np.zeros(ntraf, dtype=np.float32)
    # Area calculation
    asas.FRV_area     = np.zeros(ntraf, dtype=np.float32)
    asas.ARV_area     = np.zeros(ntraf, dtype=np.float32)
    asas.ap_free      = np.ones(ntraf, dtype=bool)

# asas is an object of the ASAS class defined in asas.py
def ConstructSSD(asas, traf):
    """ Calculates the FRV and ARV of the SSD """
    # Parameters
    N_angle = 180                   # [-] Number of points on circle (discretization)
    vmin    = traf.Vmin             # [m/s] Defined in traf.py
    vmax    = traf.Vmax             # [m/s] Defined in traf.py
    hsep    = traf.pzr * nm         # [m] Horizontal separation (5 NM)
    margin  = asas.mar              # [-] Safety margin for evasion
    hsepm   = hsep * margin         # [m] Horizontal separation with safety margin
    alpham  = 0.4999 * np.pi        # [rad] Maximum half-angle for VO

    # Relevant info from traf
    gsnorth = traf.gsnorth
    gseast  = traf.gseast
    ntraf   = traf.ntraf
    hdg     = traf.hdg

    # Local variables, will be put into asas later
    FRV_loc          = [None] * traf.ntraf
    ARV_loc          = [None] * traf.ntraf
    # For calculation purposes
    ARV_calc_loc     = [None] * traf.ntraf
    FRV_area_loc     = np.zeros(traf.ntraf, dtype=np.float32)
    ARV_area_loc     = np.zeros(traf.ntraf, dtype=np.float32)

    # Discretize the circles using points on circle
    angles = np.arange(0, 2 * np.pi, 2 * np.pi / N_angle)
    # Put points of unit-circle in a (180x2)-array (CW)
    xyc = np.transpose(np.reshape(np.concatenate((np.sin(angles), np.cos(angles))), (2, N_angle)))

    # If no traffic
    if ntraf == 0:
        return

    # If only one aircraft
    elif ntraf == 1:
        # Map them into the format ARV wants. Outercircle CCW, innercircle CW
        circle_lst = [list(map(list, np.flipud(xyc * vmax[0]))), list(map(list , xyc * vmin[0]))]
        ARV_loc[0] = circle_lst
        FRV_loc[0] = []
        ARV_calc_loc[0] = ARV_loc[0]
        # Calculate areas and store in asas
        FRV_area_loc[0] = 0
        ARV_area_loc[0] = np.pi * (vmax[0] **2 - vmin[0] ** 2)
        return

    # Construct SSD for each aircraft
    for i in range(ntraf):
        # Map them into the format pyclipper wants. Outercircle CCW, innercircle CW
        circle_tup = (tuple(map(tuple, np.flipud(xyc * vmax[i]))), tuple(map(tuple , xyc * vmin[i])))
        circle_lst = [list(map(list, np.flipud(xyc * vmax[i]))), list(map(list , xyc * vmin[i]))]

        # Calculate SSD only for aircraft in conflict (See formulas appendix)
        if asas.inconf[i]:
            # Check if indiviual ASAS has been switched on
            if not traf.asas_on[i]:
                # Map them into the format ARV wants. Outercircle CCW, innercircle CW
                ARV_loc[i] = circle_lst
                FRV_loc[i] = []
                ARV_calc_loc[i] = ARV_loc[i]
                # Calculate areas and store in asas
                FRV_area_loc[i] = 0
                ARV_area_loc[i] = np.pi * (vmax[i] **2 - vmin[i] ** 2)
                continue

            # Macke a clipper object
            pc = pyclipper.Pyclipper()
            # Add circles (ring-shape) to clipper as subject
            pc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)

            # Data of own aircraft needed to calculate distance and qdr wrt to intruder
            lat1    = traf.lat[i]
            lon1    = traf.lon[i]
            #sinlat1 = np.sin(np.radians(lat1))
            #coslat1 = np.cos(np.radians(lat1))
            r       = rwgs84(lat1)              # Approximation of Earth radius at own aircraft

            # SSD for aircraft i
            for i_other in range(ntraf):
                # Check if the other aircraft is different from the own aircraft
                if i == i_other:
                    continue
                
                # Data of intruder aircraft needed to calculate distance and qdr wrt to own aircraft
                lat2    = traf.lat[i_other]
                lon2    = traf.lon[i_other]
                #sinlat2 = np.sin(np.radians(lat2))
                coslat2 = np.cos(np.radians(lat2))

                # Combined data of own aircraft and intruder to compute relative geometry
                diff_lat = lat2 - lat1
                diff_lon = lon2 - lon1

                sin1 = np.radians(diff_lat)
                sin2 = np.radians(diff_lon)

                # Calculate distance and relative bearing w.r.t. North
                x = r * sin2 * coslat2
                y = r * sin1

                qdr = np.arctan2(x, y)
                dist = np.sqrt(x**2 + y**2)

                # In LOS the VO can't be defines, ast as id dist is on edge
                hsepm_combined = max(hsepm[i], hsepm[i_other])

                if dist < hsepm_combined:
                    dist = hsepm_combined

                # Calculate vertices of Velocity Obstacle (CCW)
                # These are still in relative velocity space, see derivation in appendix
                # Half-angle of the Velocity obstacle [rad]
                # Include safety margin
                alpha = np.arcsin(hsepm[i] / dist)
                # Limit half-angle alpha to 89.982 deg. Ensures that VO can be constructed
                if alpha > alpham:
                    alpha = alpham
                #print("print: ", hsepm[i])
                # Relevant sin/cos/tan
                sinqdr = np.sin(qdr) # East
                cosqdr = np.cos(qdr) # North
                tanalpha = np.tan(alpha)
                cosqdrtanalpha = cosqdr * tanalpha
                sinqdrtanalpha = sinqdr * tanalpha

                # Relevant x1,y1,x2,y2 (x0 and y0 are zero in relative velocity space)
                x1 = (sinqdr + cosqdrtanalpha) * 2 * vmax[i]
                x2 = (sinqdr - cosqdrtanalpha) * 2 * vmax[i]
                y1 = (cosqdr - sinqdrtanalpha) * 2 * vmax[i]
                y2 = (cosqdr + sinqdrtanalpha) * 2 * vmax[i]

                # Get vertices in an x- and y-array of size 3
                x   = np.array([gseast[i_other], x1 + gseast[i_other], x2 + gseast[i_other]])
                y   = np.array([gsnorth[i_other], y1 + gsnorth[i_other], y2 + gsnorth[i_other]])
                xy  = np.transpose(np.reshape(np.concatenate((x, y)), (2, 3)))

                # Add each other aircraft to clipper as clip
                # Scale VO when not in LOS
                # Normally VO shall be added of this other a/c
                VO = pyclipper.scale_to_clipper(tuple(map(tuple, xy)))
                # Add scaled VO to clipper
                pc.AddPath(VO, pyclipper.PT_CLIP, True)

                # Execute clipper command
                FRV = pyclipper.scale_from_clipper(pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
                ARV = pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)

                # Scale back
                ARV = pyclipper.scale_from_clipper(ARV)

                # Check if ARV or FRV is empty
                if len(ARV) == 0:
                    # No aircraft in the vicinity
                    # Map them into the format ARV wants. Outercircle CCW, innercircle CW
                    ARV_loc[i] = []
                    FRV_loc[i] = circle_lst
                    ARV_calc_loc[i] = []
                    # Calculate areas and store in asas
                    FRV_area_loc[i] = np.pi * (vmax[i] **2 - vmin[i] ** 2)
                    ARV_area_loc[i] = 0
                elif len(FRV) == 0:
                    # Should not happen with one a/c or no other a/c in the vicinity.
                    # These are handled earlier. Happens when RotA has removed all
                    # Map them into the format ARV wants. Outercircle CCW, innercircle CW
                    ARV_loc[i] = circle_lst
                    FRV_loc[i] = []
                    ARV_calc_loc[i] = circle_lst
                    # Calculate areas and store in asas
                    FRV_area_loc[i] = 0
                    ARV_area_loc[i] = np.pi * (vmax[i] **2 - vmin[i] ** 2)
                else:
                    # Check multi exteriors, if this layer is not a list, it means it has no exteriors
                    # In that case, make it a list, such that its format is consistent with further code
                    if not type(FRV[0][0]) == list:
                        FRV = [FRV]
                    if not type(ARV[0][0]) == list:
                        ARV = [ARV]
                    # Store in asas
                    FRV_loc[i] = FRV
                    ARV_loc[i] = ARV
                    # Calculate areas and store in asas
                    FRV_area_loc[i] = area(FRV)
                    ARV_area_loc[i] = area(ARV)

                    # Shortest way out prio, so use full SSD (ARV_calc = ARV)
                    ARV_calc = ARV
                    # Update calculatable ARV for resolutions
                    ARV_calc_loc[i] = ARV_calc

    asas.FRV          = FRV_loc
    asas.ARV          = ARV_loc
    asas.ARV_calc     = ARV_calc_loc
    asas.FRV_area     = FRV_area_loc
    asas.ARV_area     = ARV_area_loc
    return

def calculate_resolution(asas, traf):
    """ Calculates closest conflict-free point according to ruleset """
    ARV = asas.ARV_calc
    gsnorth = traf.gsnorth
    gseast  = traf.gseast
    ntraf   = traf.ntraf

    # Loop through SSDs of all aircraft
    for i in range(ntraf):
        # Only those that are in conflict need to resolve
        if asas.inconf[i] and len(ARV[i]) > 0:
            # Loop through all exteriors and append. Afterwards concatenate
            p = []
            q = []
            for j in range(len(ARV[i])):
                p.append(np.array(ARV[i][j]))
                q.append(np.diff(np.row_stack((p[j], p[j][0])), axis=0))
            p = np.concatenate(p)
            q = np.concatenate(q)
            # Calculate squared distance between edges
            l2 = np.sum(q ** 2, axis=1)
            # Catch l2 == 0 (exception)
            same = l2 < 1e-8
            l2[same] = 1.
            # Calc t
            t = np.sum((np.array([gseast[i], gsnorth[i]]) - p) * q, axis=1) / l2
            # Speed of boolean indices only slightly faster (negligible)
            # t must be limited between 0 and 1
            t = np.clip(t, 0., 1.)
            t[same] = 0.
            # Calculate closest point to each edge
            x1 = p[:,0] + t * q[:,0]
            y1 = p[:,1] + t * q[:,1]
            # Get distance squared
            d2 = (x1 - gseast[i]) ** 2 + (y1 - gsnorth[i]) ** 2
            # Sort distance
            ind = np.argsort(d2)
            x1  = x1[ind]
            y1  = y1[ind]

            #print("i = ", i, "d2_min = ", d2[ind][0])

            # Store result in asass
            asas.asase[i] = x1[0]
            asas.asasn[i] = y1[0]

            #print("east: ", x1[0], " north: ", y1[0])
            #print("gseast: ", gseast[i], " gsnorth: ", gsnorth[i])

            # asaseval should be set to True now
            if not asas.asaseval:
                asas.asaseval = True
        
        # Those that are not in conflict will be assigned zeros
        # Or those that have no solutions (full ARV)
        else:
            asas.asase[i] = 0.
            asas.asasn[i] = 0.
 
def rwgs84(latd):
    """ Calculate the earth radius with WGS'84 geoid definition """

    lat     = np.radians(latd)
    a       = 6378137.0       # [m] Major semi-axis WGS-84
    b       = 6356752.314245  # [m] Minor semi-axis WGS-84
    coslat  = np.cos(lat)
    sinlat  = np.sin(lat)
    an      = a * a * coslat
    bn      = b * b * sinlat
    ad      = a * coslat
    bd      = b * sinlat

    anan    = an * an
    bnbn    = bn * bn
    adad    = ad * ad
    bdbd    = bd * bd

    # Calculate radius in meters
    r       = np.sqrt((anan + bnbn) / (adad + bdbd))

    return r

def area(vset):
    """ This function calculates the area of the set of FRV or ARV """
    # Initialize A as it could be calculated iteratively
    A = 0
    # Check multiple exteriors
    if type(vset[0][0]) == list:
        # Calc every exterior separately
        for i in range(len(vset)):
            A += pyclipper.scale_from_clipper(pyclipper.scale_from_clipper(pyclipper.Area(pyclipper.scale_to_clipper(vset[i]))))
    else:
        # Single exterior
        A = pyclipper.scale_from_clipper(pyclipper.scale_from_clipper(pyclipper.Area(pyclipper.scale_to_clipper(vset))))
    return A