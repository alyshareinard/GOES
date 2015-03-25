#! /usr/bin/env python2.7
#
#
#
#
#
#
 
 
#import sys
import time
import pdb
#from datetime import date
from datetime import timedelta
from datetime import datetime
from numpy import nonzero
from numpy import greater
from numpy import transpose
from numpy import less
from numpy import mean

"""
Created on Thu Jan 30 17:18:42 2014

Purpose     : Convert quad diode XRS data to determine flare position

Explanation : Algorithm takes in quad diode dark-corrected frames, output from 
                XRS.07 algorithm (information on whether a flare is currently 
                occuring), roll angle (yaw_flip_flag if roll angle isn't available)
                x-y offsets.  Algorithm calculates and subracts
                the background from timesteps preceding the flare, then calculates
                flare position and converts it to heliocentric location

Use         : arcmin_val=main("EVENT_PEAK", flare_start[i], flare_peak[i], rpa_offlimb=True)

Inputs      :   XRS07_status - string indicating the status provided by XRS07 
                            ("EVENT_PEAK", "EVENT_START", "EVENT_RISE"), 
                time_flare_start - time of flare start in python datetime format
                                    or in the form YYMMDDHHMM
                time_current - current time in python datetime format
                                or in the form YYMMDDHHMM

Opt. Inputs : 
              target_num_frames_background - the desired number of background 
                                              frames. Default is 20.
              max_background_delay - the maximum amout of time the algorithm 
                                      will go back when looking for background
                                      frames. Default is 24 hours. 
              delay - the delay (buffer) period before flare start to avoid 
                  including flare data. Default is 5 minutes.  
              rpa_offlimb -- if the flare is offlimb, calculate the radius 
                             and position angle instead of location, set to False 
                             if this is not desired. Default is True
              SXI_rollangle -- needed if EXIS roll angle is not provided. By
                              default set to 0.

Outputs     : Flare location in heliographic coordinates (r and PA for off-limb
                events if desired) and a background flag indicating enough
                background frames were collected (0), fewer than the target number
                of background frames were collected (1), no background frames (2)
                or no data (3)

Opt. Outputs: None

Keywords    : 
    none


Restrictions: None



Written     : A A Reinard, June 30 2014

@author: alyshareinard
"""

def main(XRS07_status, time_flare_start, time_current, target_num_frames_background=20, 
         max_background_delay=timedelta(hours=24), delay=timedelta(minutes=5), rpa_offlimb=True, SXI_rollangle=0):
     
    ################################
    ## adjustable parameters      ##
    ################################
     
    background_flag=0 #0=good, 1=fewer background than desired, 2=no background
    t0=time.time()
    max_rollangle=20 #degrees
    max_xoffset=0.5 #arcmin
    max_yoffset=0.5 #arcmin
    default_value=0 #what roll angle will be set to if it is not available -- if this changes, need to change default of SXIrollangle above
    
    x_factor=1.5
    y_factor=1
    #This is hard coded for SXI, and will need to be changed to match XRS
    xsize=450*5 # 450 pixels * 5 arcseconds/pixel for SXI
    ysize=450*5 # " 

    # yaw_flip_flag
    # EXISIC D164
    # Spacecraft will include in telemetry a three state status item, 
    # at a rate of at least once every 5.0 seconds to indicate the orientation of the spacecraft as one of the following conditions:
    # -Within 5 degrees of upright yaw attitude
    # -Within 5 degrees of inverted yaw attitude
    # -Greater than 5 degrees off upright or inverted yaw attitude            
    # This does not specify what the values will be for each -- here we assign values
    # but this will likely need to be updated
    
    yaw_upright=-1
    yaw_inverted=1
    yaw_off=0

    ##################################
    ## end of adjustable parameters ##
    ##################################


    if XRS07_status != "EVENT_START" and XRS07_status != "EVENT_RISE" and XRS07_status != "EVENT_PEAK":
        return [-99, -99]

    #############################################################################
    ## check the time to see if it's in datetime format, if not convert it     ##  
    #############################################################################
               
    if type(time_current) != datetime:
        try:
            time_current=datetime.strptime(time_current, '%Y%m%d%H%M')
            time_flare_start=datetime.strptime(time_flare_start, '%Y%m%d%H%M')
        except:
            print "Time is not in python datetime format."
            print "YYMMDDHHMM format is also acceptable."
            return [-99, -99]
            
    ##########################        
    ## read in data files   ##
    ##########################

    #read in data file, which contains the quad diode values
    #TODO This section will need to be changed for dealing with dynamic data

    data_file = open("sxi_quadvals.dat", "rb")
    timestep = []
    q1 = []
    q2 = []
    q3 = []
    q4 = []

    for line in data_file:

        temp=(line[0:12])
        temp2=datetime.strptime(temp, '%Y%m%d%H%M')
        timestep.append(temp2)
        q1.append(float(line[14:26]))
        q2.append(float(line[39:52]))
        q3.append(float(line[27:38]))
        q4.append(float(line[53:69]))

    #read in the aux file (contains exptime, rollangle, xoffset, yoffset)
    aux_file = open("sxi_auxfile.dat", "rb")
    timestep_aux=[]
    exptime=[]
    rollangle=[]
    xoffset=[]
    yoffset=[]
    yaw_flip_flag=[]
    fov_eclipse=[]
    fov_offpoint=[]
    for line in aux_file:
        temp=(line[0:12])
        temp2=datetime.strptime(temp, '%Y%m%d%H%M')
        timestep_aux.append(temp2)
        exptime.append(float(line[14:26]))
        rollangle.append(float(line[27:35]))
        xoffset.append(float(line[36:43]))
        yoffset.append(float(line[44:51]))
        yaw_flip_flag.append(float(line[52:59]))
        fov_eclipse.append(float(line[60:67]))
        fov_offpoint.append(float(line[68:75]))


    ###########################################################
    ## check the inputs to be sure the timesteps are aligned ##
    ###########################################################

    for i in range(len(timestep)):
        if timestep[i]!=timestep_aux[i]: 
            print("Timesteps do not match between aux file and data file")
    if timestep != timestep_aux:
        print("Timesteps do not match between aux file and data file")
    


    ###########################
    ## begin data processing ##
    ###########################

    #convert flare_start_time to index into data array
    index_flare_start = min( range( len( timestep ) ), key=lambda i: abs( timestep[i] - time_flare_start) )

    #convert current_time to index into data array
    index_current = min( range( len( timestep ) ), key=lambda i: abs( timestep[i] - time_current ) )

    #check that the required data exists within 1 hour of current time
    if abs(timestep[index_current]-time_current) > timedelta(hours=1):
        return [-99, 3] #background_flag of 3 indicates no data is available

    ############################################
    ## check error flags -- more may be needed ##
    ############################################
    if fov_eclipse[index_current]!=0:
        print "eclipse in progress, results may be affected"
    
    if fov_offpoint[index_current]!=0:
        print "off-point in progress, results may be affected"

    ###########################################################################
    ### determine which values will contribute to the background calculation ##
    ###########################################################################

    #convert delay time to number of frames
    time_end_background = time_flare_start - delay
    index_end_background = min( range( len( timestep ) ), key=lambda i: abs( timestep[i] - time_end_background ) )
    index_begin_background = index_end_background - target_num_frames_background*4
    if (index_end_background == 0):
        index_begin_background = 0
        index_end_background = 0
        num_frames_background = 0
        print("no background frames")
        background_flag=2
    elif (index_begin_background < 0):
        index_begin_background = 0
    #make sure background frames are within max_background_delay (e.g. 1 day) of event start
    if (timestep[index_flare_start] - timestep[index_begin_background] > max_background_delay):
        #if not, choose a new index_begin_background that starts 1 day before the flare start
        index_begin_background = min(range(len(timestep)), key=lambda i: abs(timestep[i] - (time_flare_start - max_background_delay) ) )
    num_frames_background = index_end_background - index_begin_background


    #determine background   
    background_frames=[]
    num_frames_background=0
    for i in range(index_begin_background, index_end_background+1):
        #check that the background frames have lower signal than target frame 
        #and same exposure time, similar x and y offsets and either similar
        #roll angles or the same yaw flip
        #also check for offpointing and eclipse
        #can add other conditions/error flags as needed
        if q1[index_current] > q1[i] and q2[index_current]>q2[i] and \
           q3[index_current]>q3[i] and q4[index_current]>q4[i] and \
           exptime[index_current]==exptime[i] and \
           (yaw_flip_flag[index_current]==yaw_flip_flag[i] or \
           abs(rollangle[index_current]-rollangle[i])<max_rollangle) and \
           abs(xoffset[index_current]-xoffset[i])<max_xoffset and \
           abs(yoffset[index_current]-yoffset[i])<max_yoffset and \
           fov_eclipse[i]==0 and fov_offpoint[i]==0:
            background_frames.append(i)
            num_frames_background+=1
            
    #choose the most recent target number of background frames
    if num_frames_background>target_num_frames_background:
        background_frames=background_frames[num_frames_background-target_num_frames_background:]
        num_frames_background=target_num_frames_background

    if num_frames_background<target_num_frames_background:
        print("background only contains {0} frames").format(num_frames_background)
        background_flag=1
        
    #determine background frames and background roll angle x/y offset
    q1_background_frames=[]
    q2_background_frames=[]
    q3_background_frames=[]
    q4_background_frames=[]
    background_rollangle=[]
    background_xoffset=[]
    background_yoffset=[]
    for i in background_frames:
        q1_background_frames.append(q1[i])
        q2_background_frames.append(q2[i])
        q3_background_frames.append(q3[i])
        q4_background_frames.append(q4[i])
        background_rollangle.append(rollangle[i])
        background_xoffset.append(xoffset[i])
        background_yoffset.append(yoffset[i])
        
    #determine current quad-signal, roll angle and x/y offset    
    q1_current=q1[index_current]
    q2_current=q2[index_current]
    q3_current=q3[index_current]
    q4_current=q4[index_current]


    current_rollangle=rollangle[index_current]
    current_xoffset=xoffset[index_current]
    current_yoffset=yoffset[index_current]

    #if we don't have information about the current roll angle, we need to 
    #look at the yaw_flip_flag        
    reverse=False
    if current_rollangle==default_value:
        #the yaw_flip_flag indicates whether the EXIS N is within 5 degrees from
        #solar north, 5 degrees from solar south or other.  
        if yaw_flip_flag[index_current]==yaw_off:
            print "yaw flip flag indicates EXIS roll angle is more than 5 degrees offset -- calculations may be inaccurate"
        elif yaw_flip_flag[index_current]==yaw_inverted:
            reverse=True

    if num_frames_background >0:
        q1_background=mean(q1_background_frames)
        q2_background=mean(q2_background_frames)
        q3_background=mean(q3_background_frames)
        q4_background=mean(q4_background_frames)
    else:
        print "no background frames"
        q1_background=0
        q2_background=0
        q3_background=0
        q4_background=0
        background_flag=2

    # calculate the background subtracted quad-diod values    
    q1_current_subtracted=q1_current-q1_background
    q2_current_subtracted=q2_current-q2_background
    q3_current_subtracted=q3_current-q3_background
    q4_current_subtracted=q4_current-q4_background


    q_current=q1_current+q2_current+q3_current+q4_current
    q_sum=q1_current_subtracted+q2_current_subtracted+q3_current_subtracted+q4_current_subtracted
    if q_sum == 0 or q_current==0:
        print "Warning: sum of quad diode values is zero", q_sum, q_current

        return [-99, 3]
    
    ###############################
    ## calculate flare position  ##
    ###############################
    
    #These two lines calculate the flare position from quad-diode value
    X_position=-(((q2_current_subtracted+q4_current_subtracted)-(q1_current_subtracted+q3_current_subtracted)))/q_sum
    Y_position=-(((q1_current_subtracted+q2_current_subtracted)-(q3_current_subtracted+q4_current_subtracted)))/q_sum
    
    #if the yaw_flip_flag indicates that EXIS is within 5 degrees of solar S, flip the positions
    if reverse:
        X_position=-X_position
        Y_position=-Y_position
    
    #adjust for x any y offsets, if any
    X_position=X_position - current_xoffset
    Y_position=Y_position - current_yoffset

    #This if statement checks for EXIS roll angle first, this SXI roll angle
    if current_rollangle==default_value:
        #check and see if SXI_rollangle is available -- if it is, use that
        if SXI_rollangle!=default_value:
            current_rollangle=SXI_rollangle
            print "Using SXI roll angle"

    #If roll angle is available, rotate the position
    if current_rollangle!=default_value:     
        xloc_rot = math.cos(current_rollangle)*X_position - math.sin(current_rollangle)*Y_position
        yloc_rot = math.sin(current_rollangle)*X_position + math.cos(current_rollangle)*Y_position
    else:
        print "Roll angle is not available -- flare position accuracy will be limited"
        xloc_rot=X_position
        yloc_rot=Y_position

    #Convert to arcminute based on instrument geometry
    X_arcmin=xloc_rot*xsize/60.
    Y_arcmin=yloc_rot*ysize/60.

    #use subroutine to convert from arcminute to helio coordinates
    helio_or_par=arcmin2hel(X_arcmin, Y_arcmin, date=time_current, rpa_offlimb=rpa_offlimb)


    #arcmin2hel will return position angle and radius for events off the limb
    #print output and if we've reached the peak, save to file. 

    if len(helio_or_par)==1:
        helio_coord=helio_or_par
        print "flare location is: ", helio_coord
        if XRS07_status == "EVENT_PEAK":
            #if we're at the event peak, write inal coordinate to file
            f = open('XRS10_output.txt', 'a')
            output=str(time_current)+ "   "+helio_coord[0]+"\n"
            f.write(output)
            f.close()
            
    elif len(helio_or_par)==2:
        pa=helio_or_par[0]
        r=helio_or_par[1]
        print "PA: ", pa, "r: ", r
        if XRS07_status == "EVENT_PEAK":
            f = open('XRS10_output.txt', 'a')
            output=str(time_current)+"    {0:2.0f}".format(pa[0])+"   {0:2.0f}".format(r[0])+"\n"
            f.write(output)
            f.close()
    
    #append the background flag to flare position and return 
    helio_or_par.append(background_flag)
    return helio_or_par
    
    
    
    
    # -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 17:18:42 2014

converted and adapted from ARCMIN2HEL as in solarsoft (IDL)

Purpose     : Convert arcmin from sun centre to heliographic coords, 
              using the correct (perspective) transformation.

Explanation : Converts an (x,y) position given in arcmins relative to the
              solar disk centre to heliographic coordinates for the date
              supplied (default date = today).  The formula used is
              valid for any viewpoint -- though results can be "funny"
              from viewpoints inside the Sun.

Use         : helio = arcmin2hel(xx,yy,date=date)

Inputs      : xx  -  E/W coordinate in arc minutes relative to sun center
                    (West is positive); can be a vector
              yy  -  S/N coordinate in arc minutes relative to sun center
                     (North is positive); can be a vector

Opt. Inputs : rpa_offlimb -- if the flare is offlimb, calculate the radius 
                             and position angle instead of location

Outputs     : Function returns a 2xN element vector: [lat, long] in
              degrees, where N is number of elements in XX and YY

Opt. Outputs: None

Keywords    : 
    date      -  date/time in python date format
    offlimb  - flag which is true if the coordinates are beyond
                           the solar limb.
    rsun      - solar radius (input in arcsecs)
    radius    - solar radius (output in arcmins)
    p   - If specified, override pb0r for p angle
    b0    - If specified, override pb0r for b angle
    r0    - If specified, override pb0r for solar apparent
    dia:  r0 is the observer's distance from the
         Sun, in solar radii (one may use "zunits" to
         convert from solar radii to, say, kilometers).
   sphere    - If set, return longitude over the whole sphere
        rather than mirror reflecting around the 
        lon=90 meridian.  (Useful for nonzero B angle)


Restrictions: If the supplied coordinates are outside the solar disk, the
              region is projected onto the limb.

Prev. Hist. : Original by J P Wuelser.
                Translated into Python by A. A. Reinard

Written     : CDS/IDL version by C D Pike, RAL, 6 Sept 93
                Python version by A A Reinard, 28 May 2014

@author: alyshareinard

I removed vstruct and stereo from the IDL implementation
"""

import math
import numpy as np
#from datetime import *
#import pdb

def arcmin2hel(xx_in, yy_in, date="", rpa_offlimb=True, error="", radius="",
p="",b0="",r0="",rsun="",sphere=0,_extra="", backside=0, angles=""):
    dummy=[-99., -99.]
    radeg=180./math.pi
    if date=="":
        date=datetime.utcnow()

    if type(xx_in) == list and type(yy_in) == list:
        x_rows=len(xx_in)
        if len(xx_in) != len(yy_in): 
            raise ValueError("Input parameters not compatible") 
            return dummy
    else:
         x_rows=1


    offlimb = [0]*x_rows
    if x_rows == 1:
        offlimb=0.0

#-- check for overriding keywords
    angles=[0.0, 0.0, 0.0]
     
    need_b0=1-is_number(b0)
    need_rad=(1-is_number(r0)) and (1-is_number(radius))

    if need_rad or need_b0:
        angles = pb0r(date)

#-- allow keywords to override individual angles

    if is_number(p): angles[0] = p
    if is_number(b0): angles[1] = b0
    if is_number(r0): angles[2] = radeg*math.atan(1.0/r0)*60. # radians-to-arcmin conversion
    if is_number(rsun): angles[2]=rsun/60.

    b0_r = angles[1]/radeg
    radius = angles[2]
    robs = 1.0/math.tan(radius/radeg/60.)
    if type(xx_in) == list and len(xx_in)>1:
        xxat = [math.tan(float(xi)/60./radeg) for xi in xx_in] #(Convert to radians & tanify)
        yyat = [math.tan(float(yi)/60./radeg) for yi in yy_in] #(Convert to radians & tanify)
    else:
        xxat=[math.tan(float(xx_in)/60./radeg)]
        yyat=[math.tan(float(yy_in)/60./radeg)]

# Convert to cylindrical angular coordinates and azimuth -- makes
# the final transformation easier.  Here, ra is the angle out from
# centerline; phi is the azimuth.  This reduces the problem to 2-D
# geometry in the observer -- Sun-center -- viewpoint plane.


    rat2 = [xi*xi +yi*yi for xi, yi in zip(xxat, yyat)]

    w_rat2 = nonzero(greater(rat2,0))
    w_rat2=w_rat2[0]
    phi=[0]*len(w_rat2)
    if len(w_rat2) > 1:
        for irat2 in w_rat2:
            phi[irat2] = math.atan2(xxat[irat2], yyat[irat2])
    elif len(w_rat2)==1:
        phi=math.atan2(xxat[0], yyat[0])
    max_ra = math.asin(1.0/robs)
    max_rat2 = math.tan(max_ra)*math.tan(max_ra)
    ii = nonzero(greater(rat2, max_rat2)) 
    ii=ii[0]

    if len(ii)== 1:
        rat2 = [max_rat2]
        offlimb = 1
    else:
        for ival in ii:
            rat2[ival] = max_rat2
            offlimb[ival] = 1


#
# Solving for the intersection of the line of sight with the sphere
# gives a z-coordinate (toward the observer) of
#   z = R * (sin(ra))^2 +/- sqrt( Ro^2 - (sin(ra))^2 * R^2 )
# with Ro = the solar radius, ra the angular displacement from disk 
# center, and R the viewpoint distance from Sun center.
#
# We normally want the positive branch, which represents the front
# side of the Sun; but who knows? Someone may want the opposite.
# 
# AAR: leaving this in in case future development wants to give the forecaster 
# the option of saying the flare is on the backside


    ras2 = [0]*len(rat2)
    ras2=[1.0/(1.0+1.0/(val+0.0)) for val in rat2]     

    d1=[(1.0-i) for i in ras2] #force the value to be positive

    d1=[0.0 if i<0.0 else i for i in d1]

    d2=[(1-(robs*robs)*i) for i in ras2]
    d2=[0.0 if i<0.0 else i for i in d2]

    if not backside: 
        x = [ras2[i]*robs + math.sqrt(d1[i]) * math.sqrt(d2[i]) for i in range(len(ras2))]
    else:
       # This branch is for the far side of the sun
        x = [ras2[i]*robs - math.sqrt(d1[i]) * math.sqrt(d2[i]) for i in range(len(ras2))]

    rat2_gt0=list(rat2)
    if [rat2_gt0<0]==True:
        rat2_gt0[rat2_gt0<0]=0
    rr =[ math.sqrt(rat2_gt0[i]) * (robs - x[i]) for i in range(len(x))]

# Now we can convert back to xyz coords and do the 
# helioraphic conversion.  x: towards obs., y: west, z: North

    if len(rr)>1:
        xyz =[x,[math.sin(phii)*rri for phii, rri in zip(phi, rr)],[math.cos(phii)*rri for phii, rri in zip(phi, rr)]]
    if len(rr)==1:
        xyz =[x,[math.sin(phi)*rr[0]],[math.cos(phi)*rr[0]]]
        
#---------------------------------------------------------------------------
#  rotate around y axis to correct for B0 angle (B0: hel. lat. of diskcenter)
#---------------------------------------------------------------------------
    rotmx = [[math.cos(b0_r), 0.0, math.sin(b0_r)], [0.0, 1.0, 0.0], [-math.sin(b0_r), 0.0, math.cos(b0_r)]]

    xyz=transpose(xyz)
    xyz=np.dot(transpose(rotmx),transpose(xyz))
    xyz=transpose(xyz)
    
#---------------------------------------------------------------------------
#  calculate latitude and longitude.
#---------------------------------------------------------------------------

    latitude = [math.asin(xyz[i,2]) for i in range(len(xyz[:,2]))]
    toobig=nonzero(greater(latitude, 89.99/radeg))
    toobig=toobig[0]
    if len(toobig)>0:
        latitude[toobig]=89.99/radeg
    toosmall=nonzero(less(latitude, -89.99/radeg))  
    toosmall=toosmall[0]
    if len(toosmall)>0:
        latitude[toosmall]=-89.99/radeg
    longitude = [math.atan2(xyz[i,1], xyz[i,0]) for i in range(len(xyz[:,0]))] 

#---------------------------------------------------------------------------
#  longitude may be larger than 90 degrees due to nonzero B0: get proper value
#---------------------------------------------------------------------------

    if (1-sphere):
        ii = nonzero(less(xyz[:, 0], 0.0))
        ii=ii[0]
        if len(ii)>=1:
            tmp = xyz[ii, :]
            if len(ii)>1:
                tmp_l=[longitude[iii] for iii in ii]
            else:
                tmp_l = [longitude[ii]]
            jj = nonzero(greater(tmp[:, 1], -0.0000001)) # should be GE 0
            jj=jj[0] 
            if len(jj)>=1 and len(tmp_l)>1:
                tmp_l[jj] = math.pi-tmp_l[jj]
                for index in range(len(ii)):
                    longitude[ii[index]] = tmp_l[index]
            elif len(jj)==1 and len(tmp_l)==1:
                tmp_l = math.pi-tmp_l[0]
                longitude[ii] = tmp_l
            if len(ii)>1:
                for index in range(len(ii)):
                    tmp_l[index] = longitude[ii[index]]
            else:
                tmp_l=longitude[ii]
            jj = nonzero(less(tmp[:, 1], 0.0))
            jj=jj[0]
            if type(tmp_l)==float: tmp_l=[tmp_l]
            if type(jj) == list and len(jj)>=1 and len(tmp_l)>1:
                for index in range(len(jj)):
                    tmp_l[jj[index]] = -math.pi-tmp_l[jj[index]]
#                pdb.set_trace()
                for index in range(len(ii)):
                    longitude[ii[index]] = tmp_l[index]
            elif (type(jj)==float or len(jj)==1) and (len(tmp_l)==1):
                tmp_l = -math.pi-tmp_l[0]
                longitude[ii] = tmp_l      

#---------------------------------------------------------------------------
#  convert to degrees. 
#---------------------------------------------------------------------------

    # before this the routine can handle any 1D array, at this point we assume 
    # one value for longitude and latitude

    longitude=longitude[0]*radeg
    latitude=latitude[0]*radeg

    if offlimb==1 and rpa_offlimb==True:
#    if longitude>90 or latitude>90:
        print "flare is off limb"
        #need to calculate r and PA
        pa=math.atan(-xx_in/yy_in)
        pa=pa*radeg
        if yy_in < 0:
            pa=pa+180.0
        if pa < 0:
            pa=pa+360.0
        if xx_in==0 and yy_in==0:
            pa=-1
        r=math.sqrt(xx_in**2+yy_in**2)/16 #16 arcseconds = solar radius
        print "PA: ", pa, "r: ", r
        return [pa, r]

    if latitude<0:
        helio_coord="S"
    else:
        helio_coord="N"

    helio_coord=helio_coord+"{0:02.0f}".format(abs(latitude))
    
    if longitude<0:
        temp="E"
    else:
        temp="W"

    helio_coord=helio_coord+temp+"{0:02.0f}".format(abs(longitude))
    
    return [helio_coord] 
     
            


def is_number(inarray):
    """
        Test string (or array of strings) to see if it is a number
    """
#SAMPLE CALLING SEQUENCE:
#	out = is_number(inarray)
#	print, is_number(v)
#	print, is_number('xx')
#	out = is_number('66.6e')
#INPUT:
#	inarray - The string(s) to test
#OUTPUT:
#	out	- Boolean array (0 means not a string, 1 means it is)
#METHOD:
#	Use READS and trap on any errors with ON_IOERROR
#HISTORY:
#	Written 29-Oct-97 by M.Morrison
#       Modified 29-Jun-99, Zarro (SM&A/GSFC) - added check for undefined input
#       Modified 29-Sep-00, Zarro (EIT/GSFC) - added check for invalid inputs
#       Modified 21-jun-05, Csillaghy (UAS Switzerland)- added call to is_number2
#       Converted from IDL to Python 28-May-14, Reinard 
    #if it's not one of these options (it's a structure, class or something else) we return 0 
    if not (type(inarray) == list or type(inarray) == float or type(inarray) == int 
    or type(inarray) == str or type(inarray) == long):
        return 0        

    #if it's a list, we go through the items of the list
    if type(inarray)==list:   
        if len(inarray) == 0:  
            return 0
        out=[]
        for i in inarray:
            try:
                float(i)
                out.append(1)
            except ValueError:
                out.append(0)
        return out
    
    #if it's not a list, we just try to make it a float
    try:
        float(inarray)
        return 1
    except ValueError:
        return 0

     
"""
Name        : PB0R

Purpose     : To calculate the solar P, B0 angles and the semi-diameter.

Explanation : Uses semi-rigorous formulae to calculate the solar P (position
               angle of pole) and B0 (latitude of point at disk centre) angles
               and also the semi-diameter of the solar disk at the date/time
               requested.

Use         : IDL> ang = pb0r(date_time,/soho)

Inputs      : date_time  -  the date/time specified in any CDS format

Outputs     : Function returns a 3-element array with
                                   ang(0)  = P  (degrees)
                                   ang(1)  = B0 (degrees)
                                   ang(2)  = R  semi-diameter (arcmin or
                                                               arcsec if
                                                               keyword set)
                                                               
Keywords    : SOHO - if present the semi-diameter returned is as viewed
                     from the SOHO spacecraft
              ARCSEC - returns semi-diameter in arcsec rather than arcmins
              ERROR  - Output keyword containing error message;
                       a null string is returned if no error occurs
              RETAIN - passed to get_orbit to determine whether orbit file
                       is left open or not.
              EARTH  - set internal environment variable to EARTH
                       view
              L0     - L0 value [degrees]
              STEREO = 'A' or 'B' for STEREO Ahead or Behind

Common      : pb0r_common (internal common block)

Category    : Util, coords

Prev. Hist. : Based on Fortran programs by Hohenkerk and Emerson (RGO)

Written     : CDS/IDL version, C D Pike, RAL, 16-May-94

Modified    : Update semi-diameter calculation, CDP, 20-May-94
        Version 3, William Thompson, GSFC, 14 November 1994
               Modified .DAY to .MJD
        Version 4, CDP, 10-Jan-96
               Add SOHO/ARCSEC keywords and make NOW the default.
        Version 5, Liyun Wang, GSFC/ARC, March 12, 1996
               Modified such that point of view can be changed to
               SOHO if the env variable SC_VIEW is set to 1
               (via the call to USE_SOHO_VIEW)
                Added ERROR keyword
        Version 6, Liyun Wang, GSFC/ARC, March 14, 1996
                Replaced call to GET_ORBIT with the IDL call_function
        Version 7, Liyun Wang, GSFC/ARC, March 21, 1996
                Modified such that if no orbit file is found, earth
                view is used
        Version 8, Liyun Wang, GSFC/ARC, April 10, 1996
                Set SC_VIEW to 0 if no orbit files are found
        Version 9, February 6, 1997, Liyun Wang, NASA/GSFC
                Changed call to ANYTIM2JD instead of CDS2JD
        Version 10, April 17 1997, CDP.  Added RETAIN keyword
        Version 11, July 28 1997, DMZ, fixed bug in common block
        Version 12, Nov 17 1997, DMZ, added /EARTH
        Version 13, 26 Jan 1998, William Thompson
                Correct by 1% if no orbit files are found, instead of
                setting SC_VIEW to 0.
                Fix bug involving when to recalculate SOHO positions.
        Version 14, 7 Jan 1999, Zarro (SMA/GSFC) 
                Fixed another bug involving SC_VIEW
                (deprecated EARTH keyword)
        Version 15, 20 Jan 1999, Zarro (SMA/GSFC)
                Added check for GET_ORBIT in !path
        Version 16, 06-Feb-2003, William Thompson, GSFC
                Fixed bug in common block (Previous fix somehow lost?)
        Version 17, 20-Feb-2003, Zarro (EER/GSFC) 
                Added check for IMAGE_TOOL running and removed
                      silly 'goto'
                Modified, 8-Jan-2005, Zarro (L-3Com/GSFC) - added
                      /DEBUG
                Modified, 23-Oct-2007, Zarro (ADNET) 
                     - added SOHO B0 angle
                Modified, 20-Feb-2008, Zarro (ADNET)
                     - added L0 keyword
                Modified, 21-Aug-2008, Zarro (ADNET)
                     - added STEREO by-pass
                Modified, 21-Feb-2009, Zarro (ADNET)
                     - added /VERBOSE
        Version 18, 28-May-2014, Reinard (SWPC/NOAA)
                Translated code to Python
"""
def  pb0r(date_utc, arcsec="N", error="N", earth="Y",debug="N", l0=0.,roll_angle=0.,verbose="N"):

    import datetime
    radeg=180./math.pi
    dtor=math.pi/180.
    
    recal=1
 
#---------------------------------------------------------------------------
#  date supplied?
#---------------------------------------------------------------------------

    if date_utc != 0.0:
#        print type(date_utc)

# require date in python datetime format so as to remove anytim2utc and associated programs
    
        if type(date_utc)!=datetime.datetime and type(date_utc)!=datetime.date:
            print "date is invalid, must be in datetime format"
            return
                          
#---------------------------------------------------------------------------
# TODO does it need to recalculate?
#this was in the original with a common block -- I could keep this and 
#input the old information if that helped with speed. Leaving as is for now
#while I test the speed
#---------------------------------------------------------------------------


    if recal:
        
#---------------------------------------------------------------------------
#  number of Julian days since 2415020.0
#---------------------------------------------------------------------------

        jd = date_to_julian_day(date_utc)

        de=jd-2415020


#---------------------------------------------------------------------------
#  get the longitude of the sun etc.
#---------------------------------------------------------------------------
        temp=sun_pos(de) 
        longmed=temp[0]
        appl=temp[3]
        oblt=temp[4]

#---------------------------------------------------------------------------
#  form aberrated longitude
#---------------------------------------------------------------------------
        lambda1 = longmed - (20.5/3600.0)

#---------------------------------------------------------------------------
#  form longitude of ascending node of sun's equator on ecliptic
#---------------------------------------------------------------------------
        node = 73.666666 + (50.25/3600.0)*( (de/365.25) + 50.0)
        arg = lambda1 - node

#---------------------------------------------------------------------------
#  calculate P, the position angle of the pole
#---------------------------------------------------------------------------
        p = (math.atan(-math.tan(oblt*dtor) * math.cos(appl*dtor)) + 
        math.atan( -0.12722 * math.cos(arg*dtor))) * radeg

#---------------------------------------------------------------------------
#  ... and B0 the tilt of the axis
#---------------------------------------------------------------------------
        b = math.asin( 0.12620 * math.sin(arg*dtor) ) * radeg

#---------------------------------------------------------------------------
#  ... and the semi-diameter
#
#
#  Form the mean anomalies of Venus(MV),Earth(ME),Mars(MM),Jupiter(MJ)
#  and the mean elongation of the Moon from the Sun(D).
#
#---------------------------------------------------------------------------
        t = de/36525.0

        mv = 212.6   + ( (58517.80   * t) % 360.0 )
        me = 358.476 + ( (35999.0498 * t) % 360.0 )
        mm = 319.5   + ( (19139.86   * t) % 360.0 )
        mj = 225.3   + ( ( 3034.69   * t) % 360.0 )
        d = 350.7    + ( (445267.11  * t) % 360.0 )

#---------------------------------------------------------------------------
#  Form the geocentric distance(r) and semi-diameter(sd)
#---------------------------------------------------------------------------
        r = 1.000141 - (0.016748 - 0.0000418*t)*math.cos(me*dtor) 
        - 0.000140 * math.cos(2.0*me*dtor)                       
        + 0.000016 * math.cos((58.3 + 2.0*mv - 2.0*me)*dtor) 
        + 0.000005 * math.cos((209.1 + mv - me)*dtor)            
        + 0.000005 * math.cos((253.8 - 2.0*mm + 2.0*me)*dtor)
        + 0.000016 * math.cos(( 89.5 - mj + me)*dtor)            
        + 0.000009 * math.cos((357.1 - 2.0*mj + 2.0*me)*dtor) 
        + 0.000031 * math.cos(d*dtor)

        sd = (0.2665685/r)*60.0

        output = [p, b, sd] 



        if arcsec.upper()=="Y":
            return [output[0], output[1], output[2]*60.] 
        else:
            return output




    
"""
Name        : SUN_POS
              
Purpose     : Calculate solar ephemeris parameters.
              
Explanation : Allows for planetary and lunar perturbations in the calculation
              of solar longitude at date and various other solar positional
              parameters.
              
Use         : temp= sun_pos(date)
        temp[0] = longitude, temp[1]=ra, temp[2]=dec, temp[3]=app_long, temp[4]=obliq
   
Inputs      : date - fractional number of days since JD 2415020.0 
              
Opt. Inputs : None
              
Outputs     : longitude  -  Longitude of sun for mean equinox of date (degs)
              ra         -  Apparent RA for true equinox of date (degs)
              dec        -  Apparent declination for true equinox of date (degs)
              app_long   -  Apparent longitude (degs)
              obliq      -  True obliquity (degs)
              
Opt. Outputs: All above
              
Keywords    : None

Calls       : None

Common      : None
              
Restrictions: None
              
Side effects: None
              
Category    : Util, coords
              
Prev. Hist. : From Fortran routine by B Emerson (RGO).

Written     : CDS/IDL version by C D Pike, RAL, 17-May-94
              
Modified    : 

Version     : Version 1, 17-May-94
              Version 2, 28-May 2014, A A Reinard (SWPC/NOAA), translated code to Python

"""            
def sun_pos(dd):
    import math
    radeg=180./math.pi
    dtor=math.pi/180.

#  This routine is a truncated version of Newcomb's Sun and
#  is designed to give apparent angular coordinates (T.E.D) to a
#  precision of one second of time

#  form time in Julian centuries from 1900.0
    t = dd/36525.0

#  form sun's mean longitude

    l = (279.696678+((36000.768925*t) % 360.0))*3600.0

#  allow for ellipticity of the orbit (equation of centre)
#  using the Earth's mean anomoly ME

    me = 358.475844 + ((35999.049750*t) % 360.0)
    ellcor  = (6910.1 - 17.2*t)*math.sin(me*dtor) + 72.3*math.sin(2.0*me*dtor)
    l = l + ellcor

# allow for the Venus perturbations using the mean anomaly of Venus MV

    mv = 212.603219 + ((58517.803875*t) % 360.0) 
    vencorr = 4.8 * math.cos((299.1017 + mv - me)*dtor) 
    + 5.5 * math.cos((148.3133 +  2.0 * mv  -  2.0 * me )*dtor)
    + 2.5 * math.cos((315.9433 +  2.0 * mv  -  3.0 * me )*dtor)
    + 1.6 * math.cos((345.2533 +  3.0 * mv  -  4.0 * me )*dtor)
    1.0 * math.cos((318.15   +  3.0 * mv  -  5.0 * me )*dtor)
    l = l + vencorr

#  Allow for the Mars perturbations using the mean anomaly of Mars MM

    mm = 319.529425  +  (( 19139.858500 * t)  %  360.0 )
    marscorr = 2.0 * math.cos((343.8883 -  2.0 * mm  +  2.0 * me)*dtor )
    + 1.8 * math.cos((200.4017 -  2.0 * mm  + me) * dtor)
    l = l + marscorr

# Allow for the Jupiter perturbations using the mean anomaly of
# Jupiter MJ

    mj = 225.328328  +  (( 3034.6920239 * t)  %  360.0 )
    jupcorr = 7.2 * math.cos(( 179.5317 - mj + me )*dtor) 
    + 2.6 * math.cos((263.2167  -  mj ) *dtor) 
    + 2.7 * math.cos(( 87.1450  -  2.0 * mj  +  2.0 * me ) *dtor)
    + 1.6 * math.cos((109.4933  -  2.0 * mj  +  me ) *dtor)
    l = l + jupcorr

# Allow for the Moons perturbations using the mean elongation of
# the Moon from the Sun D

    d = 350.7376814  + (( 445267.11422 * t)  %  360.0 )
    mooncorr  = 6.5 * math.sin(d*dtor)
    l = l + mooncorr

# Allow for long period terms

    longterm  = + 6.4 * math.sin(( 231.19  +  20.20 * t )*dtor)
    l  =    l + longterm
    l  =  ( l + 2592000.0)  %  1296000.0
    longmed = l/3600.0

# Allow for Aberration

    l  =  l - 20.5

# Allow for Nutation using the longitude of the Moons mean node OMEGA

    omega = 259.183275 - (( 1934.142008 * t ) % 360.0 )
    l  =  l - 17.2 * math.sin(omega*dtor)

# Form the True Obliquity

    oblt  = 23.452294 - 0.0130125*t + (9.2*math.cos(omega*dtor))/3600.0

# Form Right Ascension and Declination

    l = l/3600.0
    ra  = math.atan2( math.sin(l*dtor) * math.cos(oblt*dtor), math.cos(l*dtor))  * radeg
    
    if (ra < 0.0):  ra = ra + 360.0

    dec = math.asin(math.sin(l*dtor) * math.sin(oblt*dtor)) * radeg
    
    return((longmed, ra, dec, l, oblt))

 
def date_to_julian_day(my_date):
    """Returns the Julian day number of a date.
        Written by A A Reinard """
    a = (14 - my_date.month)//12
    y = my_date.year + 4800 - a
    m = my_date.month + 12*a - 3
    day_frac=my_date.hour/24.+my_date.minute/24./60. - 0.5 #julian day starts at noon
    return day_frac+my_date.day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045
 


    
