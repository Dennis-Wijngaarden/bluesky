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
    InitializeSSD(asas, traf.ntraf)
    ConstructSSD(asas, traf)

def resolve(asas, traf):
    """ Resolve all current conflicts """
    # Check if ASAS is ON first!
    if not asas.swasas:
        return

    #InitializeSSD(asas, traf.ntraf)

    # Construct the SSD
    #ConstructSSD(asas, traf)

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
    gs_ap   = traf.ap.tas
    hdg_ap  = traf.ap.trk
    apnorth = np.cos(hdg_ap / 180 * np.pi) * gs_ap
    apeast  = np.sin(hdg_ap / 180 * np.pi) * gs_ap

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
        if True: #asas.inconf[i]:
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

                geofence_polygons = calc_geofence_polygon(traf, i, i_other, 100, hsepm_combined)
                if geofence_polygons != None:
                    # Convert to form pyclipper wants it
                    VOs = pyclipper.scale_to_clipper(tuple(map(lambda x: tuple(map(lambda y: tuple(y), x)), geofence_polygons)))
                    pc.AddPaths(VOs, pyclipper.PT_CLIP, True)

                if asas.priocode == "RS5":
                    if pyclipper.PointInPolygon(pyclipper.scale_to_clipper((apeast[i],apnorth[i])),VO):
                        asas.ap_free[i] = False

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
    ARV     = asas.ARV_calc
    if asas.priocode == "RS5":
        gsnorth = np.cos(traf.ap.trk / 180 * np.pi) * traf.ap.tas
        gseast  = np.sin(traf.ap.trk / 180 * np.pi) * traf.ap.tas
    else:
        gsnorth = traf.gsnorth
        gseast  = traf.gseast
    ntraf   = traf.ntraf

    # Loop through SSDs of all aircraft
    for i in range(ntraf):
        # Only those that are in conflict need to resolve
        if asas.inconf[i] and len(ARV[i]) > 0:
            # First check if AP-setting is free
            if asas.ap_free[i] and asas.priocode == "RS5":
                asas.asase[i] = gseast[i]
                asas.asasn[i] = gsnorth[i]
            else:
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

                # Check solution (within geofence and if path is free at CPA)
                # Store result in asass
                asas.asase[i] = x1[0]
                asas.asasn[i] = y1[0]

                # asaseval should be set to True now
                if not asas.asaseval:
                    asas.asaseval = True
        
        # Those that are not in conflict will be assigned zeros
        # Or those that have no solutions (full ARV)
        else:
            asas.asase[i] = 0.
            asas.asasn[i] = 0.

def calc_geofence_polygon(traf, i_own, i_int, N, seperation):
    """Calculates the polygon due to a geofence with resolution n"""
    if traf.gs[i_int] < 0.1:
        return None
    # first calculate the relative position of the intruder with respect to the own_aircraft
    r = rwgs84(traf.lat[i_own]) # Earth radius [m] approximated using wgs84

    lat_diff = np.radians(traf.lat[i_int] - traf.lat[i_own])
    lon_diff = np.radians(traf.lon[i_int] - traf.lon[i_own])
    if lon_diff > np.pi:
        lon_diff = lon_diff - 2. * np.pi

    coslat = np.cos(np.radians(traf.lat[i_own]))

    x_int_x = -lon_diff * coslat * r # Relative position of intruder in x (east) direction
    x_int_y = -lat_diff * r          # Relative position of intryder in y (north) direction
    v_int_x = traf.gseast[i_int]        # Intruder velocity [m/s] in east direction
    v_int_y = traf.gsnorth[i_int]       # Intruder velocity [m/s] in north direction

    # Geofence equations
    # Calculate relative positions of the geofence points wrt to the own aircraft
    geofence_xy = []

    for i in range(len(traf.geofence[i_own])):
        lat_diff = np.radians(traf.geofence[i_own][i][0] - traf.lat[i_own])
        lon_diff = np.radians(traf.geofence[i_own][i][1] - traf.lon[i_own])
        if lon_diff > np.pi:
            lon_diff = lon_diff - 2. * np.pi
        geofence_xy.append([lon_diff * coslat * r, lat_diff * r])

    geofence_a = []
    geofence_n = []
    geofence_distv = []
    geofence_dist = []
    geofence_rotation = []
    for i in range(len(geofence_xy)):
        if i == len(geofence_xy) - 1:
            i_next = 0
        else:
            i_next = i + 1
        a, n = line_vec_from_vectors(np.array([geofence_xy[i]]).T, np.array([geofence_xy[i_next]]).T)
        geofence_a.append(a)
        geofence_n.append(n)
        geofence_distv.append(a - np.vdot(a, n) * n)
        geofence_dist.append(np.sqrt(geofence_distv[i][0][0]**2 + geofence_distv[i][1][0]**2))
        geofence_rotation.append(np.arctan2(geofence_distv[i][1][0], geofence_distv[i][0][0]) - 0.5 * np.pi)
    
    # Now loop over geofence and return polygon
    geofence_polygons = []
    for i in range(len(geofence_xy)):
        rotation_geo = geofence_rotation[i]
        v_int_perp = v_int_y * np.cos(rotation_geo) - v_int_x * np.sin(rotation_geo)
        v_int_parallel = v_int_x * np.cos(rotation_geo) + v_int_y * np.sin(rotation_geo)

        d_geo = geofence_dist[i]

        x_int_parallel  = x_int_x * np.cos(rotation_geo) + x_int_y * np.sin(rotation_geo)
        x_int_perp      = x_int_y * np.cos(rotation_geo) - x_int_x * np.sin(rotation_geo)

        phi = 0.5 * np.arctan2(-x_int_parallel, x_int_perp)

        # Create non rotated ellipse
        K_2_parallel = 1. - x_int_perp / d_geo * np.sin(phi)**2 / (np.cos(phi)**2 - np.sin(phi)**2)
        K_parallel = -np.sin(phi) * v_int_perp * (2. + x_int_perp / d_geo)
        c_parallel = -K_parallel / (2. * K_2_parallel)

        K_2_perp = 1. + x_int_perp / d_geo * np.cos(phi)**2 / (np.cos(phi)**2 - np.sin(phi)**2)
        K_perp = -np.cos(phi) * v_int_perp * (2. + x_int_perp / d_geo)
        c_perp = -K_perp / (2. * K_2_perp)

        K = c_parallel**2 * K_2_parallel + c_perp**2 * K_2_perp - v_int_perp**2

        a_2 = K / K_2_parallel
        b_2 = K / K_2_perp

        n_points = N

        # if an ellipse
        if b_2 > 0.:
            if v_int_perp < 0.:
                continue
            a = np.sqrt(a_2)
            b = np.sqrt(b_2)

            a_2_over_b_2 = a_2/b_2

            if a_2_over_b_2 > 200.:
                a_2_over_b_2 = 200.

            angles = np.linspace(-np.pi, np.pi, n_points)
            angles = angles - np.sin(angles * 2.) * a_2_over_b_2 / 400.

            v_res_x_basic = a * np.cos(angles)
            v_res_y_basic = b * np.sin(angles)

        # If an hyperbola
        else:
            a = np.sqrt(a_2)
            b = np.sqrt(-b_2)
            
            x = abs(c_parallel) + 2. * traf.Vmax[i_own]
            
            angle_max = np.log(x / a + np.sqrt(x**2 / a**2 - 1.))
            
            if phi >= 0.:
                angles = np.linspace(angle_max, -angle_max, n_points)
                v_res_x_basic = a * np.cosh(angles)
                v_res_y_basic = b * np.sinh(angles)
            else:
                angles = np.linspace(-angle_max, angle_max, n_points)
                v_res_x_basic = -a * np.cosh(angles)
                v_res_y_basic = b * np.sinh(angles)

        v_res_x_rotated = v_res_x_basic * np.cos(phi + rotation_geo) - v_res_y_basic * np.sin(phi + rotation_geo) 
        v_res_y_rotated = v_res_y_basic * np.cos(phi + rotation_geo) + v_res_x_basic * np.sin(phi + rotation_geo) 
            
        c_x = c_parallel * np.cos(phi + rotation_geo) - c_perp * np.sin(phi + rotation_geo) + v_int_parallel * np.cos(rotation_geo)
        c_y = c_perp * np.cos(phi + rotation_geo) + c_parallel * np.sin(phi + rotation_geo) + v_int_parallel * np.sin(rotation_geo)
            
        v_res_x = v_res_x_rotated + c_x
        v_res_y = v_res_y_rotated + c_y

        geofence_polygons.append(np.array([v_res_x, v_res_y]).T)
        
    return geofence_polygons

def line_vec_from_vectors(vector0, vector1):
    # Vector equation x = a + tn
    a = vector0
    diff = vector1 - vector0
    n = diff / np.sqrt(diff[0][0]**2 + diff[1][0]**2)
    return a, n  

 
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