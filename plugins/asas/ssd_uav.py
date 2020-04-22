''' Conflict resolution based on the SSD algorithm. '''
from bluesky.traffic.asas import ConflictResolution
from bluesky.tools import geo, areafilter
from bluesky.tools.aero import nm, kts
import numpy as np
# Try to import pyclipper
try:
    import pyclipper
except ImportError:
    print("Could not import pyclipper, RESO SSD will not function")
    
def init_plugin():

    # Addtional initilisation code

    # Configuration parameters
    config = {
        # The name of your plugin
        'plugin_name':     'SSD_UAV',

        # The type of this plugin. For now, only simulation plugins are possible.
        'plugin_type':     'sim'
    }

    # init_plugin() should always return these two dicts.
    return config, {}

class SSDUAV(ConflictResolution):
    def setprio(self, flag=None, priocode=''):
        '''Set the prio switch abd the type of prio '''
        if flag is None:
            return True, "PRIORULES [ON/OFF] [PRIOCODE]" + \
                            "\nAvailabe priority codes: " + \
                            "\n     RS1:  Shortest way out" + \
                            "\n     RS2:  Shortest from target heading" + \
                            "\n     RS3:  Heading change" + \
                            "\n     RS4:  Speed change" + \
                            "\nPriority is currently " + ("ON" if self.swprio else "OFF") + \
                            "\nPriority code is currently: " + \
                str(self.priocode)
        options = ["RS1", "RS2", "RS3", "RS4"]
        if priocode not in options:
            return False, "Priority code not Understood. Available Options: " + str(options)
        return super().setprio(flag, priocode)
    
    def loaded_pyclipper():
        """ Return true if pyclipper is succesfully loaded """
        import sys
        return "pyclipper" in sys.modules
    
    def resolve(self, conf, ownship, intruder):
        # Initialize SSD variables with ntraf
        self.initializeSSD(conf, ownship.ntraf)
        
        self.constructSSD(conf, ownship, self.priocode)
        # Get resolved speed-vector
        self.calculate_resolution(conf, ownship)
        
        # Now assign resolutions to variables in the ASAS class
        # Start with current states, need a copy, otherwise it changes traf!
        newtrack = np.copy(ownship.hdg)
        newgs = np.copy(ownship.gs)
        # Calculate new track and speed
        # No need to cap the speeds, since SSD implicitly caps
        new_trk = np.arctan2(conf.asase, conf.asasn) * 180 / np.pi
        new_gs = np.sqrt(conf.asase ** 2 + conf.asasn ** 2)

        # Sometimes an aircraft is in conflict, but no solutions could be found
        # In that case it is assigned 0 by ASAS, but needs to handled
        asas_cmd = np.logical_and(conf.inconf, new_gs > 0)

        # Assign new track and speed for those that are in conflict
        newtrack[asas_cmd] = new_trk[asas_cmd]
        newgs[asas_cmd] = new_gs[asas_cmd]
        # Not needed as it is a 2D-implementation...
        newvs = ownship.vs

        # Cap the velocity
        newgscapped = np.maximum(ownship.perf.vmin, np.minimum(ownship.perf.vmax, newgs))

        alt = ownship.selalt

        return newtrack, newgscapped, newvs, alt
        
    def initializeSSD(self, conf, ntraf):
        # Need to do it here, since ASAS.reset doesn't know ntraf
        conf.FRV = [None] * ntraf
        conf.ARV = [None] * ntraf
        # For calculation purposes
        conf.ARV_calc = [None] * ntraf
        conf.inrange = [None] * ntraf
        # Stores resolution vector, also used in visualization
        conf.asasn = np.zeros(ntraf, dtype=np.float32)
        conf.asase = np.zeros(ntraf, dtype=np.float32)
        # Area calculation
        conf.FRV_area = np.zeros(ntraf, dtype=np.float32)
        conf.ARV_area = np.zeros(ntraf, dtype=np.float32)
        
        conf.ap_free = np.ones(ntraf, dtype=bool)
        
    def constructSSD(self, conf, ownship, priocode="RS1"):
        """ Calculates the FRV and ARV of the SSD """
        # Parameters
        N_angle = 180 # [-] Number of points on circle discretization
        hsep = conf.rpz # [m] Horizontal separation
        margin = self.resofach # [m] Safety margin for avasion (settings.asas_mar)
        hsepm = hsep * margin # [m] Horizontal separation with safety margin
        alpham = 0.4999 * np.pi # [rad] Maximum half-angle for VO
        betalos = np.pi / 4. # [rad] Minimum divertion angle for LoS (45 deg seems optimal)
        adsbmax = 10000. # {m} Maximum ADS-B range
        beta = np.pi / 4. + betalos / 2.
        
        # Relevant info from traf
        gsnorth = ownship.gsnorth
        gseast = ownship.gseast
        lat = ownship.lat
        lon = ownship.lon
        ntraf = ownship.ntraf
        hdg = ownship.hdg
        gs_ap = ownship.ap.tas
        hdg_ap = ownship.ap.trk
        apnorth = np.cos(hdg_ap / 180. * np.pi) * gs_ap
        apeast = np.sin(hdg_ap / 180. * np.pi) * gs_ap
        
        # Local variables, will be put into asas later
        FRV_loc = [None] * ownship.ntraf
        ARV_loc = [None] * ownship.ntraf
        # For calculation purposes
        ARV_calc_loc = [None] * ownship.ntraf
        FRV_area_loc = np.zeros(ownship.ntraf, dtype=np.float32)
        ARV_area_loc = np.zeros(ownship.ntraf, dtype=np.float32)
        
        # Use velocity limits for the ring-shaped part of the SSD
        # Discretize the circles using points on circle
        angles = np.arange(0., 2. * np.pi, 2. * np.pi / N_angle)
        # Put points of unit-circle in a (180x2)-array (Cw)
        xyc = np.transpose(np.reshape(np.concatenate((np.sin(angles), np.cos(angles))), (2, N_angle)))
    
        # If no traffic
        if ntraf == 0:
            return
        
        # Function qdrdist_matrix needs 4 vectors as input (lat1,lon1,lat2,lon2)
        # To be efficient, calculate all qdr and dist in one function call
        # Example with ntraf = 5:   ind1 = [0,0,0,0,1,1,1,2,2,3]
        #                           ind2 = [1,2,3,4,2,3,4,3,4,4]
        # This way the qdrdist is only calculated once between every aircraft
        # To get all combinations, use this function to get the indices
        ind1, ind2 = self.qdrdist_matrix_indices(ntraf)
        # Get absolute bearing [deg] and distance [nm]
        # qdr is defined from [-180,180] deg, w.r.t. North
        [qdr, dist] = geo.qdrdist_matrix(lat[ind1], lon[ind1], lat[ind2], lon[ind2])
        # Put result of function from matrix to ndarray
        qdr = np.reshape(np.array(qdr), np.shape(ind1))
        dist = np.reshape(np.array(dist), np.shape(ind1))
        # SI-units from [deg] to [rad]
        qdr = np.deg2rad(qdr)
        # Get distance from [nm] to [m]
        dist = dist * nm
        
        # In LoS the VO can't be defined, act as if dist is on edge
        dist[dist < hsepm] = hsepm
        
        # Calculate vertices of Velocity Obstacle (CCW)
        # These are still in relative velocity space
        # Half-angle of the Velocity obstacle [rad]
        # Include safety margin
        alpha = np.arcsin(hsepm / dist)
        # Limit half-angle alpha to 89.982 deg. Ensures that VO can be constructed
        alpha[alpha > alpham] = alpham
        # Relevant sin/cos/tan
        sinqdr = np.sin(qdr)
        cosqdr = np.cos(qdr)
        tanalpha = np.tan(alpha)
        cosqdrtanalpha = cosqdr * tanalpha
        sinqdrtanalpha = sinqdr * tanalpha
        
        # Consider every aircraft
        for i in range(ntraf):
            # Calculate SSD only for aircraft in conflict (See formulas appendix)
            if conf.inconf[i]:
                
                vmin = ownship.perf.vmin[i]
                vmax = ownship.perf.vmax[i]

                # wind field around ownship
                vn_wind, ve_wind = ownship.wind.getdata(ownship.lat[i], ownship.lon[i], ownship.alt[i])
                
                # in the first time step, ASAS runs before perf, which means that his value will be zero
                # and the SSD cannot be constructed
                if vmin == vmax == 0:
                    continue
                if (priocode == 'RS3'):
                    v_selected = np.sqrt(gsnorth[i]**2 + gseast[i]**2)
                    circle_tup = (tuple(map(tuple, np.flipud((xyc * (v_selected + 0.001)) + np.array([ve_wind, vn_wind])))), tuple(map(tuple, (xyc * (v_selected - 0.001)) + np.array([ve_wind, vn_wind]))))
                    circle_lst = [list(map(list, np.flipud((xyc * (v_selected + 0.001)) + np.array([ve_wind, vn_wind])))), list(map(list, (xyc * (v_selected - 0.001)) + np.array([ve_wind, vn_wind])))]
                else:
                    vmin = max(vmin, 0.001)
                    # Map them into the format pyclipper wants. Outercircle CCW, innercircle CW
                    circle_tup = (tuple(map(tuple, np.flipud((xyc * vmax) + np.array([ve_wind, vn_wind])))), tuple(map(tuple, (xyc * vmin) + np.array([ve_wind, vn_wind]))))
                    circle_lst = [list(map(list, np.flipud((xyc * vmax) + np.array([ve_wind, vn_wind])))), list(map(list, (xyc * vmin) + np.array([ve_wind, vn_wind])))]
                
                # Relevant x1,y1,x2,y2 (x0 and y0 are zero in relative velocity space)
                x1 = (sinqdr + cosqdrtanalpha) * 20. * vmax
                x2 = (sinqdr - cosqdrtanalpha) * 20. * vmax
                y1 = (cosqdr - sinqdrtanalpha) * 20. * vmax
                y2 = (cosqdr + sinqdrtanalpha) * 20. * vmax
                
                # SSD for aircraft i
                # Get indices that belong to aircraft i
                ind = np.where(np.logical_or(ind1 == i, ind2 == i))[0]
                
                # Check whether there are any aircraft in the vicinity
                if len(ind) == 0:
                    # No aircraft in the vicinity
                    # Map them into the format ARV wants. Outercircle CCW, innercircle CW
                    ARV_loc[i] = circle_lst
                    FRV_loc[i] = []
                    ARV_calc_loc[i] = ARV_loc[i]
                    # Calculate areas and store in asas
                    FRV_area_loc[i] = 0.
                    ARV_area_loc[i] = np.pi * (vmax ** 2 - vmin ** 2)
                else:
                    # The i's of the other aircraft
                    i_other = np.delete(np.arange(0, ntraf), i)
                    # Aircraft that are within ADS-B range
                    ac_adsb = np.where(dist[ind] < adsbmax)[0]
                    # Now account for ADS-B range in indices of other aircraft (i_other)
                    ind = ind[ac_adsb]
                    i_other = i_other[ac_adsb]
                    # Put it in class-object
                    conf.inrange[i] = i_other
                    # VO from 2 to 1 is mirror of 1 to 2. Only 1 to 2 can be constructed in
                    # this manner, so need a correction vector that will mirror the VO
                    fix = np.ones(np.shape(i_other))
                    fix[i_other < i] = -1
                    # Relative bearing [deg] from [-180,180]
                    # (less required conversions than rad in RotA)
                    fix_ang = np.zeros(np.shape(i_other))
                    fix_ang[i_other < i] = 180.
                    
                    # Get vertices in an x- and y-array of size (ntraf-1)*3x1
                    x = np.concatenate((gseast[i_other],
                                        x1[ind] * fix + gseast[i_other],
                                        x2[ind] * fix + gseast[i_other]))
                    y = np.concatenate((gsnorth[i_other],
                                        y1[ind] * fix + gsnorth[i_other],
                                        y2[ind] * fix + gsnorth[i_other]))
                    # Reshape [(ntraf-1)x3] and put arrays in one array [(ntraf-1)x3x2]
                    x = np.transpose(x.reshape(3, np.shape(i_other)[0]))
                    y = np.transpose(y.reshape(3, np.shape(i_other)[0]))
                    xy = np.dstack((x, y))
                    
                    # Make a clipper object
                    pc = pyclipper.Pyclipper()
                    # Add circles (ring-shape) to clipper as subject
                    if priocode == "RS4":
                        hdg_sel = hdg[i] * np.pi / 180
                        xyp = np.array([[np.sin(hdg_sel - 0.0087), np.cos(hdg_sel - 0.0087)],
                                        [0, 0],
                                        [np.sin(hdg_sel + 0.0087), np.cos(hdg_sel + 0.0087)]],
                                       dtype=np.float64)
                        part = pyclipper.scale_to_clipper(tuple(map(tuple, 20. * vmax * xyp)))
                        pc2 = pyclipper.Pyclipper()
                        pc2.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                        pc2.AddPath(part, pyclipper.PT_CLIP, True)
                        vel_cone_scaled = pc2.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO)
                        pc.AddPaths(vel_cone_scaled, pyclipper.PT_SUBJECT, True)
                    else:
                        pc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)

                    # If there is a geofence, calculate relative variables:
                    geofence_defined = False
                    try:
                        areafilter.areas['GF']
                    except:
                        pass
                    else:
                        geofence_defined = True
                        # Load geofence shape
                        geofence = areafilter.areas['GF']
                        # Loop through geofence coordinates
                        coordinates = np.reshape(geofence.coordinates, (int(len(geofence.coordinates) / 2), 2))
                        qdrs_gf = np.array([]) # [deg] in hdg CW
                        dists_gf = np.array([]) # [m]
                        for k in range(len(coordinates)):
                            # Calculate relative qdrs and distances of geofence points w.r.t. ownship
                            qdr_gf, dist_gf = geo.qdrdist(ownship.lat[i], ownship.lon[i], coordinates[k][0], coordinates[k][1])
                            qdrs_gf = np.append(qdrs_gf, qdr_gf)
                            dists_gf = np.append(dists_gf, dist_gf * nm)

                        qdrs_gf_rad = np.deg2rad(qdrs_gf)
                        xs_gf = dists_gf * np.sin(qdrs_gf_rad) # [m] East
                        ys_gf = dists_gf * np.cos(qdrs_gf_rad) # [m] North

                        
                        # Generate data for each geofence segment 0 to 1, 1 to 2, 2 to 3 ..... n to 0.
                        dxs_gf = np.array([])
                        dys_gf = np.array([])
                        for k in range(len(coordinates)):
                            x_from = xs_gf[k]
                            y_from = ys_gf[k]
                            # if last elament (needs to be connected to first element)
                            if k == (len(coordinates) - 1):
                                x_to = xs_gf[0]
                                y_to = ys_gf[0]
                            else:
                                x_to = xs_gf[k + 1]
                                y_to = ys_gf[k + 1]
                            
                            dxs_gf = np.append(dxs_gf, x_to - x_from)
                            dys_gf = np.append(dys_gf, y_to - y_from)

                        # calculate values (phis) of rotation of geofence segments
                        phis_gf = np.arctan2(dys_gf, dxs_gf)
                        x_hats_prime = np.transpose(np.array([np.cos(phis_gf), np.sin(phis_gf)]))
                        y_hats_prime = np.transpose(np.array([-np.sin(phis_gf), np.cos(phis_gf)]))

                    
                    # Add each other other aircraft to clipper as clip
                    for j in range(np.shape(i_other)[0]):
                        # Scale VO when not in LOS
                        if dist[ind[j]] > hsepm:
                            # Normally VO shall be added of this other a/c
                            VO = pyclipper.scale_to_clipper(tuple(map(tuple, xy[j, :, :])))
                        else:
                            # Pair is in LOS, instead of triangular VO, use darttip
                            # Check if bearing should be mirrored
                            if i_other[j] < i:
                                qdr_los = qdr[ind[j]] + np.pi
                            else:
                                qdr_los = qdr[ind[j]]
                            # Length of inner-leg of darttip
                            leg = 20. * vmax / np.cos(beta) * np.array([1, 1, 1, 0])
                            # Angles of darttip
                            angles_los = np.array([qdr_los + 2 * beta, qdr_los, qdr_los - 2 * beta, 0.])
                            # Calculate coordinates (CCW)
                            x_los = leg * np.sin(angles_los)
                            y_los = leg * np.cos(angles_los)
                            # Put in array of correct format
                            xy_los = np.vstack((x_los, y_los)).T
                            # Scale darttip
                            VO = pyclipper.scale_to_clipper(tuple(map(tuple, xy_los)))
                        # Add scaled VO to clipper
                        pc.AddPath(VO, pyclipper.PT_CLIP, True)
                        
                        if priocode == "RS2":
                            if pyclipper.PointInPolygon(pyclipper.scale_to_clipper((apeast[i], apnorth[i])), VO):
                                conf.ap_free[i] = False
                        #elif priocode == "RS4":
                        #    hdg_sel = hdg[i] * np.pi / 180
                        #    xyp = np.array([[np.sin(hdg_sel + 0.0087), np.cos(hdg_sel + 0.0087)],
                        #                    [0, 0],
                        #                    [np.sin(hdg_sel - 0.0087), np.cos(hdg_sel - 0.0087)]],
                        #                   dtype=np.float64)
                        #    part = pyclipper.scale_to_clipper(tuple(map(tuple, 2.1 * vmax * xyp)))
                        #    pc.AddPath(part, pyclipper.PT_SUBJECT, True)

                        # Add Geofence if available
                        if (geofence_defined):
                            # Determine relative distance vector w.r.t. intruder
                            qdr_int, dist_int = geo.qdrdist(ownship.lat[i], ownship.lon[i], ownship.lat[i_other[j]], ownship.lon[i_other[j]])
                            qdr_int_rad = np.deg2rad(qdr_int)
                            x_int = dist_int * nm * np.sin(qdr_int_rad)
                            y_int = dist_int * nm * np.cos(qdr_int_rad)
                            d_int = np.array([x_int, y_int])
                            trk_int = np.deg2rad(ownship.trk[i_other[j]])
                            gs_int = np.deg2rad(ownship.gs[i_other[j]])
                            v_int = np.array([gs_int * np.sin(trk_int), gs_int * np.cos(trk_int)])

                            v_int_dot_y_hats_prime = np.array([])
                            for k in range(len(y_hats_prime)):
                                v_int_dot_y_hats_prime = np.append(v_int_dot_y_hats_prime, np.dot(v_int, y_hats_prime[k]))
                            candidate_gf_segments = np.where(v_int_dot_y_hats_prime < 0)[0]

                            d_int_dot_x_hats_prime = np.array([])
                            d_int_dot_y_hats_prime = np.array([])

                            ds_geo = np.array([]) # [m] Array of distanced w.r.t. geofence of wonship
                            for k in candidate_gf_segments:
                                d_int_dot_x_hats_prime = np.append(d_int_dot_x_hats_prime, np.dot(d_int, x_hats_prime[k]))
                                d_int_dot_y_hats_prime = np.append(d_int_dot_y_hats_prime, np.dot(d_int, y_hats_prime[k]))

                                ds_geo = np.append(ds_geo, -np.dot(np.array([xs_gf[k], ys_gf[k]]), y_hats_prime[k]))
                            
                            phis_prime_gf = 0.5 * np.arctan2(-1. * d_int_dot_x_hats_prime, d_int_dot_y_hats_prime)

                            # Total rotation angle
                            phis_total_gf = phis_gf[candidate_gf_segments] + phis_prime_gf

                            # Secondary axis system primary axes
                            x_hats_2prime = np.transpose(np.array([np.cos(phis_total_gf), np.sin(phis_total_gf)]))
                            y_hats_2prime = np.transpose(np.array([-np.sin(phis_total_gf), np.cos(phis_total_gf)]))
                            
                            # Arrays of dot products
                            d_int_dot_x_hats_2prime = np.array([])
                            d_int_dot_y_hats_2prime = np.array([])
                            v_int_dot_x_hats_2prime = np.array([])
                            v_int_dot_y_hats_2prime = np.array([])
                            d_int_dot_v_int = np.array([])

                            for k in range(len(x_hats_2prime)):
                                d_int_dot_x_hats_2prime = np.append(d_int_dot_x_hats_2prime, np.dot(d_int, x_hats_2prime[k]))
                                d_int_dot_y_hats_2prime = np.append(d_int_dot_y_hats_2prime, np.dot(d_int, y_hats_2prime[k]))

                                v_int_dot_x_hats_2prime = np.append(v_int_dot_x_hats_2prime, np.dot(v_int, x_hats_2prime[k]))
                                v_int_dot_y_hats_2prime = np.append(v_int_dot_y_hats_2prime, np.dot(v_int, y_hats_2prime[k]))

                                d_int_dot_v_int = np.append(d_int_dot_v_int, np.dot(d_int, v_int))

                            # Constants needed to compute geometry of geofence VO's
                            C1s = 1. + np.sin(phis_prime_gf) * d_int_dot_x_hats_2prime / ds_geo
                            C2s = 1. + np.cos(phis_prime_gf) * d_int_dot_y_hats_2prime / ds_geo
                            C3s = -2. * v_int_dot_x_hats_2prime - np.sin(phis_prime_gf) * d_int_dot_v_int / ds_geo
                            C4s = -2. * v_int_dot_y_hats_2prime - np.cos(phis_prime_gf) * d_int_dot_v_int / ds_geo
                            
                            # Center points of geofence VO geometries in double rotated axis system
                            Cxs_2prime = - C3s / (2. * C1s)
                            Cys_2prime = - C4s / (2. * C2s)

                            # semi major axes squared (in case of ellipse)
                            a2s = (- gs_int**2 + C2s * Cys_2prime**2) / C1s + Cxs_2prime**2
                            b2s = (- gs_int**2 + C1s * Cxs_2prime**2) / C2s + Cys_2prime**2

                            # Loop trough a2s and b2s to construct VOs, categorize them 
                            for k in range(len(a2s)):
                                # If Ellipse
                                if (b2s[k] > 0):
                                    ellipse_angles = np.linspace(0., 2. * np.pi, N_angle)
                                    rotated_xs = np.sqrt(a2s[k]) * np.cos(ellipse_angles) + Cxs_2prime[k]
                                    rotated_ys = np.sqrt(b2s[k]) * np.sin(ellipse_angles) + Cys_2prime[k]
                                # If hyperbola
                                else:
                                    tmax = np.log((20. * vmax + np.sqrt(20.**2 * vmax**2 + a2s[k])) / np.sqrt(a2s[k]))
                                    tmin = -tmax
                                    t = np.linspace(tmax, tmin, N_angle)
                                    if (phis_prime_gf[k] > 0):
                                        rotated_xs = -np.sqrt(a2s[k]) * np.cosh(t) + Cxs_2prime[k]
                                    else:
                                        rotated_xs = np.sqrt(a2s[k]) * np.sinh(t) + Cys_2prime[k]
                                non_rotated_xs = rotated_xs * np.cos(phis_total_gf[k]) - rotated_ys * np.sin(phis_total_gf[k])
                                non_rotated_ys = rotated_xs * np.sin(phis_total_gf[k]) + rotated_ys * np.cos(phis_total_gf[k])

                                # Add non roted VO's to clipper
                                xy_gf = []
                                xy_gf.append(non_rotated_xs)
                                xy_gf.append(non_rotated_ys)
                                xy_gf = np.array(xy_gf)
                                xy_gf = np.transpose(xy_gf)
                                xy_gf_tuple = tuple(map(tuple, xy_gf))

                                # Scale VO to clipper
                                VO = pyclipper.scale_to_clipper(xy_gf_tuple)
                                # Add scaled VO to clipper
                                pc.AddPath(VO, pyclipper.PT_CLIP, True)

                    # Execute clipper command
                    FRV = pyclipper.scale_from_clipper(
                        pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))

                    ARV = pyclipper.scale_from_clipper(
                        pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
                    
                    # Check if ARV or FRV is empty
                    if len(ARV) == 0:
                        # Map them into the format ARV wants. Outercircle CCW, innercircle CW
                        ARV_loc[i] = []
                        FRV_loc[i] = circle_lst
                        ARV_calc_loc[i] = []
                        # Calculate areas and store in asas
                        FRV_area_loc[i] = np.pi * (vmax ** 2 - vmin ** 2)
                        ARV_area_loc[i] = 0
                    elif len(FRV) == 0:
                        # Should not happen with one a/c or no other a/c in the vicinity.
                        # These are handled earlier. Happens when RotA has removed all
                        # Map them into the format FRV wants. Outercircle CCW, innercircle CW
                        ARV_loc[i] = circle_lst
                        FRV_loc[i] = []
                        ARV_calc_loc[i] = circle_lst
                        # Calculate areas and store in asas
                        FRV_area_loc[i] = 0
                        ARV_area_loc[i] = np.pi * (vmax ** 2 - vmin ** 2)
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
                        ARV_calc_loc[i] = ARV
                        
                        # Calculate areas and store in asas
                        FRV_area_loc[i] = self.area(FRV)
                        ARV_area_loc[i] = self.area(ARV)
                        
        conf.FRV = FRV_loc
        conf.ARV = ARV_loc
        conf.ARV_calc = ARV_calc_loc
        conf.FRV_area = FRV_area_loc
        conf.ARV_area = ARV_area_loc
        return
    
    def calculate_resolution(self, conf, ownship):
        """ Calculates closest conflict-free point according to ruleset """
        # It's just linalg, however credits to: http://stackoverflow.com/a/1501725
        # Variables
        ARV = conf.ARV_calc
        # Select AP-setting as reference point for closest to target rulesets
        if self.priocode == "RS2":
            gsnorth = np.cos(ownship.ap.trk / 180 * np.pi) * ownship.ap.tas
            gseast = np.sin(ownship.ap.trk / 180 * np.pi) * ownship.ap.tas
        else:
            gsnorth = ownship.gsnorth
            gseast = ownship.gseast
        ntraf = ownship.ntraf
        
        # Loop through SSDs of all aircraft
        for i in range(ntraf):
            # Only those that are in conflict need to resolve
            if conf.inconf[i] and ARV[i] is not None and len(ARV[i]) > 0:
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
                x1 = p[:, 0] + t * q[:, 0]
                y1 = p[:, 1] + t * q[:, 1]
                # Get distance squared
                d2 = (x1 - gseast[i]) ** 2 + (y1 - gsnorth[i]) ** 2
                # Sort distance
                ind = np.argsort(d2)
                x1 = x1[ind]
                y1 = y1[ind]
            
                # Store result in conf
                conf.asase[i] = x1[0]
                conf.asasn[i] = y1[0]
            # Those that are not in conflict will be assigned zeros
            # Or those that have no solutions (full ARV)    
            else:
                conf.asase[i] = 0.
                conf.asasn[i] = 0.
            
    def area(self, vset):
        """ This function calculates the area of the set of FRV or ARV """
        # Initialize A as it could be calculated iteratively
        A = 0
        # Check multiple exteriors
        if type(vset[0][0]) == list:
            # Calc every exterior separately
            for i in range(len(vset)):
                A += pyclipper.scale_from_clipper(
                    pyclipper.scale_from_clipper(pyclipper.Area(pyclipper.scale_to_clipper(vset[i]))))
        else:
            # Single exterior
            A = pyclipper.scale_from_clipper(
                pyclipper.scale_from_clipper(pyclipper.Area(pyclipper.scale_to_clipper(vset))))
        return A
                        
    def qdrdist_matrix_indices(self, ntraf):
        """ This function gives the indices that can be used in the lon/lat-vectors """
        # The indices will be n*(n-1)/2 long
        # Only works for n >= 2, which is logical...
        # This is faster than np.triu_indices :)
        tmp_range = np.arange(ntraf - 1, dtype=np.int32)
        ind1 = np.repeat(tmp_range, (tmp_range + 1)[::-1])
        ind2 = np.ones(ind1.shape[0], dtype=np.int32)
        inds = np.cumsum(tmp_range[1:][::-1] + 1)
        np.put(ind2, inds, np.arange(ntraf * -1 + 3, 1))
        ind2 = np.cumsum(ind2, out=ind2)
        return ind1, ind2