''' Conflict resolution based on the SSD algorithm plotter plugin '''
from bluesky import traf
from bluesky import settings
from bluesky.tools import geo
from bluesky.tools.aero import nm
import numpy as np
import matplotlib
matplotlib.use("Qt5agg")
from matplotlib import pyplot as plt
# Try to import pyclipper
try:
    import pyclipper
except ImportError:
    print("Could not import pyclipper, RESO SSD will not function") 


class Conflict(object):
    def __init__(self):
        self.rpz = settings.asas_pzr * nm
        self.inconf = np.ones(traf.ntraf)

# Define derived class from SSDUAV
class SSDUAVPlotter(object):
    def __init__(self):
        self.priocode = 'RS1'
        self.resofach = settings.asas_mar
        self.resofacv = settings.asas_mar
        self.flight_ids = []
        self.plots = []
        return
    
    def resolve(self, conf, ownship, intruder):
        # Initialize SSD variables with ntraf
        self.initializeSSD(conf, ownship.ntraf)
        
        self.constructSSD(conf, ownship, self.priocode)
        # Get resolved speed-vector
        self.calculate_resolution(conf, ownship)
        return
        
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
                
                # in the first time step, ASAS runs before perf, which means that his value will be zero
                # and the SSD cannot be constructed
                if vmin == vmax == 0:
                    continue
                
                # Map them into the format pyclipper wants. Outercircle CCW, innercircle CW
                vmin = max(0, 0.001)
                circle_tup = (tuple(map(tuple, np.flipud(xyc * vmax))), tuple(map(tuple, xyc * vmin)))
                circle_lst = [list(map(list, np.flipud(xyc * vmax))), list(map(list, xyc * vmin))]
                
                # Relevant x1,y1,x2,y2 (x0 and y0 are zero in relative velocity space)
                x1 = (sinqdr + cosqdrtanalpha) * 2. * vmax
                x2 = (sinqdr - cosqdrtanalpha) * 2. * vmax
                y1 = (cosqdr - sinqdrtanalpha) * 2. * vmax
                y2 = (cosqdr + sinqdrtanalpha) * 2. * vmax
                
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
                    pc.AddPaths(pyclipper.scale_to_clipper(circle_tup), pyclipper.PT_SUBJECT, True)
                    
                    
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
                            leg = 1.1 * vmax / np.cos(beta) * np.array([1, 1, 1, 0])
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
                    
                     # Execute clipper command
                    FRV = pyclipper.scale_from_clipper(pc.Execute(pyclipper.CT_INTERSECTION, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))

                    ARV = pyclipper.scale_from_clipper(pc.Execute(pyclipper.CT_DIFFERENCE, pyclipper.PFT_NONZERO, pyclipper.PFT_NONZERO))
                    
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
        
    def process_stack(self, idx, onoff):
        # Check case on if id does exists
        flight_id = traf.id[idx]
        if (onoff == True):
            if (flight_id in self.flight_ids):
                return
            else:
                self.add_plot(flight_id)
        else:
            if (flight_id in self.flight_ids):
                self.close_plots([self.flight_ids.index(flight_id)])
            else:
                return
        
    def update_plots(self):
        # Check if all flight_ids still exists, else delete plot and from flight_ids list
        close_idxs = []
        for i in range(len(self.flight_ids)):
            if (traf.id2idx == -1):
                close_idxs.append(i)
        self.close_plots(close_idxs)
        
        # Generate ARVs and FRVs
        conf = Conflict()
        self.resolve(conf, traf, traf)
        # Update each plot
        for i in range(len(self.flight_ids)):
            # Plot FRV
            traffic_idx = traf.id2idx(self.flight_ids[i])
            self.plots[i].update(conf.FRV[traffic_idx], conf.ARV[traffic_idx])
    
    def add_plot(self, flight_id):
        self.flight_ids.append(flight_id)
        self.plots.append(SSDPlotElement(flight_id))
    
    def close_plots(self, idxs):
        idxs = idxs[::-1] # from high i to low i
        for i in idxs:
            self.plots[i].close()
            del self.plots[i]
            del self.flight_ids[i]

class SSDPlotElement():
    def __init__(self, flight_id):
        plt.ion()
        self.fig, self.ax = plt.subplots()
        self.ax.axis('equal')
        self.ax.axis('off')
        self.fig.canvas.set_window_title('SSD ' + flight_id)
        self.fig.canvas.draw()
        plt.pause(0.01)
    
    def update(self, FRV, ARV):
        plt.ion()
        self.ax.clear()
        for i in range(len(FRV)):
            FRVi = np.array(FRV[i])
            x_FRVi = np.append(FRVi[:,0] , np.array(FRVi[0,0]))
            y_FRVi = np.append(FRVi[:,1] , np.array(FRVi[0,1]))
            self.ax.fill(x_FRVi, y_FRVi, color = '#FF0000') #red
                         
        for i in range(len(ARV)):
            ARVi = np.array(ARV[i])
            x_ARVi = np.append(ARVi[:,0] , np.array(ARVi[0,0]))
            y_ARVi = np.append(ARVi[:,1] , np.array(ARVi[0,1]))
            self.ax.fill(x_ARVi, y_ARVi, color = '#C0C0C0') #grey
        self.ax.axis('off')
        self.ax.axis('equal')
        self.fig.canvas.draw()
        plt.pause(0.01)
    
    def close(self):
        plt.close(fig=self.fig)
        
    
plotter = SSDUAVPlotter()

def init_plugin():
    # Configuration parameters
    config = {
    
        # The name of your plugin
        'plugin_name':     'SSD_UAV_PLOTTER',

        # The type of this plugin.
        'plugin_type':     'sim',

        # Plugin update interval in seconds.
        'update_interval': 1.0,

        # Update at the end of a timestep.
        'update':          update,

        # Update at the beginning of a timestep.
        #'preupdate':       preupdate
        }
    
    stackfunctions = {
        # The command name for your function
        'PLOTSSD': [
            # A short usage string. This will be printed if you type HELP <name> in the BlueSky console
            'PLOTSSD acid ON/OFF',

            # A list of the argument types your function accepts.
            'acid, onoff',

            # The name of your function in this plugin
            process_stack,

            # a longer help text of your function.
            '']}

    return config, stackfunctions


def update():
    plotter.update_plots()

def process_stack(idx, onoff):
    plotter.process_stack(idx, onoff)

def preupdate():

    return