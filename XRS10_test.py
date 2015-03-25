# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 08:57:08 2014

Purpose     : This program is used to test the XRS10.py program using quad-diode data 
derived from the SXI instrument and flare data taken from the NOAA flare catalog.


Explanation : Algorithm runs the XRS10.py algorithm over all the flares in the flare
            catalog. Only a small number have SXI quad-diode proxy data. The
            flare position is calculated for those cases and plots are made comparing
            results from the IDL code and the python code and then comparing
            results from the python code to the actual NOAA flare locations.

Use         : [val, flareloc]=XRS10_test()

Inputs      : None

Opt. Inputs : None
              
Outputs     : arcmin_vals - the values and background flags determined using 
                            the XRS10.py algorithm
              flare_loc - flare locations from the NOAA flare catalog

Opt. Outputs: None

Keywords    : 
    none


Restrictions: None



Written     : A A Reinard, June 30 2014

@author: alyshareinard
"""
#from array import array
import pdb
from XRS10 import *
import time
import matplotlib.pyplot as plt
plt.ion()



def XRS10_test():
    """
    Program to test XRS10.py code -- cycles through flares (in separate file) and has XRS10.py 
    compute the flare location, then compiles the XRS10.py results
    """
    data_file = open("flarefile.dat", "rb")
    flare_start= []
    flare_peak = []
    flare_loc = []
    flare_size = []
    for line in data_file:
        flare_start.append((line[0:12]))
        flare_peak.append((line[13:25]))
        flare_loc.append(line[26:32])
        flare_size.append(line[33:36])
    data_file=open("idl_results.dat", "rb")
    ew_idl=[]
    ns_idl=[]
    for line in data_file:
        ew_idl.append(float(line[0:14]))
        ns_idl.append(float(line[14:28]))
        
    arcmin_vals=[]
    background_flags=[]
    timing=[]
    for i in range(len(flare_start)):
        if(i> -1 and i < 10000):
            if ns_idl[i] > -99:
                print "i", i
                t1=time.time()
                arcmin_val=main("EVENT_PEAK", flare_start[i], flare_peak[i], rpa_offlimb=False)
#                print arcmin_val
                if ns_idl[i]<0:
                    ns="S"
                else:
                    ns="N"
                if ew_idl[i]<0:
                    ew="E"
                else:
                    ew="W"
                print "IDL flare location", ns+"{0:02.0f}".format(abs(ns_idl[i]))+ew+"{0:02.0f}".format(abs(ew_idl[i]))

                print "NOAA flare location", flare_loc[i]
                
                arcmin_vals.append(arcmin_val[0])
                background_flags.append(arcmin_val[1])
                t2=time.time()
                print "time for routine", t2-t1
                timing.append(t2-t1)
            else:
                arcmin_vals.append("NR")
                background_flags.append(2)
    ns_py=[]
    ew_py=[]

    print "time for routine on average", mean(timing)


    #parse the heliographic coordinate value into NS and EW signed values
    for line in arcmin_vals:
      
        if line != "NR":
            ns_val=int(line[1:3])
            try: 
                if line[0]=="S":
                    ns_val=-ns_val
                ns_py.append(ns_val)
            except:
                pdb.set_trace()
            ew_val=int(line[4:6])
            try:                
                if line[3]=="E":
                    ew_val=-ew_val
                ew_py.append(ew_val)
            except:
                pdb.set_trace()
        else:
            ew_py.append(-99)
            ns_py.append(-99)


    # This factor of 1.5 provides a near between the algorithm results and the 
    # NOAA flare locations.  I believe this is needed because the SXI images 
    # are cropped 
    # This correction should not be necessary for GOES-R XRS values, so 
    # it is not included in the XRS10.py file.
    
    ew_factor=1.5

    ew_py_new=[]
    for val in ew_py:
        if val !=-99:
            ew_py_new.append(val*ew_factor)
        else:
            ew_py_new.append(val)
            
    ew_py=ew_py_new
    
    ew_idl_new=[]
    for val in ew_idl:
        if val !=-99:
            ew_idl_new.append(val*ew_factor)
        else:
            ew_idl_new.append(val)
            
    ew_idl=ew_idl_new
    
    # make plots comparing Python and IDL results
    
    plt.figure(1)
    plt.plot(ns_py, ns_idl, 'bo')
    plt.title('All flares, EW locations')  
    plt.xlabel('EW location output from python')
    plt.ylabel('EW location output from IDL')
    plt.axis([-90, 90, -90, 90])
    plt.draw()
    
    plt.figure(2)
    plt.plot(ew_py, ew_idl, 'bo')
    plt.title('All flares, NS locations')  
    plt.xlabel('NS location output from python')
    plt.ylabel('NS location output from IDL')  
    plt.axis([-90, 90, -90, 90])
    plt.draw()
    
    #now get the actual flare location
    
    ns_flare=[]
    ew_flare=[]
    for flare in flare_loc:

        if flare[0]=="S":
            ns_flare.append(-int(flare[1:3]))
        else:
            ns_flare.append(int(flare[1:3]))
    
        if flare[3]=="E":
            ew_flare.append(-int(flare[4:6]))
        else:
            ew_flare.append(int(flare[4:6]))    


    
    #calculate the errors between python and IDL derived locations and the real locations
    err_py=[]
    err_idl=[]
    for i in range(len(ew_py)):        
        real_arcmin=hel2arcmin(ns_flare[i], ew_flare[i], date=flare_peak[i])
        if ew_py[i] != -99:
            py_arcmin=hel2arcmin(ns_py[i], ew_py[i], date=flare_peak[i])
            err_py.append(math.sqrt((py_arcmin[0]-real_arcmin[0])**2+(py_arcmin[1]-real_arcmin[1])**2))
        if ew_idl[i] != -99:
            idl_arcmin=hel2arcmin(ns_idl[i], ew_idl[i], date=flare_peak[i])
            err_idl.append(math.sqrt((idl_arcmin[0]-real_arcmin[0])**2+(idl_arcmin[1]-real_arcmin[1])**2))

    print "error in python values", mean(err_py), std(err_py)
    print "error in idl values", mean(err_idl), std(err_idl)

    #separate C, M and X class flares    
    
    Cflare_indices=[]    
    Mflare_indices=[] 
    Xflare_indices=[] 

    for i in range(len(flare_size)):
        fsize=flare_size[i]
        
        if fsize[0].upper()=="C":
            Cflare_indices.append(i)

        if fsize[0].upper()=="M":
            Mflare_indices.append(i)

        if fsize[0].upper()=="X":
            Xflare_indices.append(i)     
        

    Xerr_py=[]
    Xerr_idl=[]
    Xerr_GB_py=[]
    Xerr_GB_idl=[]
    GB=[]

    for i in Xflare_indices:  

        real_arcmin=hel2arcmin(ns_flare[i], ew_flare[i], date=flare_peak[i])
        if ew_py[i] != -99:
            py_arcmin=hel2arcmin(ns_py[i], ew_py[i], date=flare_peak[i])
            error=math.sqrt((py_arcmin[0]-real_arcmin[0])**2+(py_arcmin[1]-real_arcmin[1])**2)
            Xerr_py.append(error)
            if background_flags[i]==0:
                GB.append(i)
                Xerr_GB_py.append(error)

        if ew_idl[i] != -99:
            idl_arcmin=hel2arcmin(ns_idl[i], ew_idl[i], date=flare_peak[i])
            error=math.sqrt((idl_arcmin[0]-real_arcmin[0])**2+(idl_arcmin[1]-real_arcmin[1])**2)
            Xerr_idl.append(error)
            if background_flags[i]==0: #this assumes that IDL would have the same good background events -- reasonable assumption
                Xerr_GB_idl.append(error)

    Merr_py=[]
    Merr_idl=[]
    Merr_GB_py=[]
    Merr_GB_idl=[]
    for i in Mflare_indices:        
        real_arcmin=hel2arcmin(ns_flare[i], ew_flare[i], date=flare_peak[i])
        if ew_py[i] != -99:
            py_arcmin=hel2arcmin(ns_py[i], ew_py[i], date=flare_peak[i])
            error=math.sqrt((py_arcmin[0]-real_arcmin[0])**2+(py_arcmin[1]-real_arcmin[1])**2)
            Merr_py.append(error)
            if background_flags[i]==0:
                Merr_GB_py.append(error)
        if ew_idl[i] != -99:
            idl_arcmin=hel2arcmin(ns_idl[i], ew_idl[i], date=flare_peak[i])
            error=math.sqrt((idl_arcmin[0]-real_arcmin[0])**2+(idl_arcmin[1]-real_arcmin[1])**2)
            Merr_idl.append(error)
            if background_flags[i]==0:
                Merr_GB_idl.append(error)
                
    Cerr_py=[]
    Cerr_idl=[]        
    Cerr_GB_py=[]
    Cerr_GB_idl=[]

    for i in Cflare_indices:        
        real_arcmin=hel2arcmin(ns_flare[i], ew_flare[i], date=flare_peak[i])
        if ew_py[i] != -99:
            py_arcmin=hel2arcmin(ns_py[i], ew_py[i], date=flare_peak[i])
            error=math.sqrt((py_arcmin[0]-real_arcmin[0])**2+(py_arcmin[1]-real_arcmin[1])**2)
            Cerr_py.append(error)
            if background_flags[i]==0:
                Cerr_GB_py.append(error)
        if ew_idl[i] != -99:
            idl_arcmin=hel2arcmin(ns_idl[i], ew_idl[i], date=flare_peak[i])
            error=math.sqrt((idl_arcmin[0]-real_arcmin[0])**2+(idl_arcmin[1]-real_arcmin[1])**2)
            Cerr_idl.append(error)
            if background_flags[i]==0:
                Cerr_GB_idl.append(error)

    print "Errors: mean/std"
    print "X-class flares -- all events"
    print "Python values: {0:0.2f}".format(mean(Xerr_py)), "/{0:0.2f}".format(std(Xerr_py))
    print "IDL values: {0:0.2f}".format(mean(Xerr_idl)), "/{0:0.2f}".format(std(Xerr_idl))
    
    
    print "X-class flares -- only good backgrounds"
    print "Python values: {0:0.2f}".format(mean(Xerr_GB_py)), "/{0:0.2f}".format(std(Xerr_GB_py))
    print "IDL values: {0:0.2f}".format(mean(Xerr_GB_idl)), "/{0:0.2f}".format(std(Xerr_GB_idl))  

    print ""
    print "M-class flares -- all events"
    print "Python values: {0:0.2f}".format(mean(Merr_py)),  "/{0:0.2f}".format(std(Merr_py))
    print "IDL values: {0:0.2f}".format(mean(Merr_idl)), "/{0:0.2f}".format(std(Merr_idl))
    
    print "M-class flares -- only good backgrounds"    
    print "Python values: {0:0.2f}".format(mean(Merr_GB_py)), "/{0:0.2f}".format(std(Merr_GB_py))
    print "IDL values: {0:0.2f}".format(mean(Merr_GB_idl)), "/{0:0.2f}".format(std(Merr_GB_idl))

    print ""
    print "C-class flares -- all events"                    
    print "Python values: {0:0.2f}".format(mean(Cerr_py)), "/{0:0.2f}".format(std(Cerr_py))
    print "IDL values: {0:0.2f}".format(mean(Cerr_idl)), "/{0:0.2f}".format(std(Cerr_idl))

    print "C-class flares -- only good backgrounds"                    
    print "Python values: {0:0.2f}".format(mean(Cerr_GB_py)), "/{0:0.2f}".format(std(Cerr_GB_py))
    print "IDL values: {0:0.2f}".format(mean(Cerr_GB_idl)), "/{0:0.2f}".format(std(Cerr_GB_idl))

    plt.figure(3)
    Xew_py=[ew_py[val] for val in Xflare_indices]
    Xew_flare=[ew_flare[val] for val in Xflare_indices]

    
    plt.subplot(121)
    plt.plot(Xew_py, Xew_flare, 'bo')
    plt.plot([-90, 90], [-90, 90])
    plt.title('X-class flares EW')    
    plt.xlabel('EW location based on algorithm')
    plt.ylabel('EW location based on NOAA catalog')
    plt.xlim([-90,90])
    plt.ylim([-90,90])
    plt.axis([-90, 90, -90, 90])

    Xns_py=[ns_py[val] for val in Xflare_indices]
    Xns_flare=[ns_flare[val] for val in Xflare_indices]
    plt.subplot(122)
    plt.plot(Xns_py, Xns_flare, 'bo')
    plt.plot([-90, 90], [-90, 90])
    plt.title('X-class flares NS')    
    plt.xlabel('NS location based on algorithm')
    plt.ylabel('NS location based on NOAA catalog')
    plt.xlim([-90,90])
    plt.ylim([-90,90])
    plt.axis([-50, 50, -50, 50])
    plt.draw()
    
    plt.figure(4)
    Mew_py=[ew_py[val] for val in Mflare_indices]
    Mew_flare=[ew_flare[val] for val in Mflare_indices]
    plt.subplot(121)
    plt.plot(Mew_py, Mew_flare, 'bo')
    plt.plot([-90, 90], [-90, 90])
    plt.title('M-class flares EW')  
    plt.xlabel('EW location based on algorithm')
    plt.ylabel('EW location based on NOAA catalog')
    plt.xlim([-90,90])
    plt.ylim([-90,90])
    plt.axis([-90, 90, -90, 90])

    plt.subplot(122)
    Mns_py=[ns_py[val] for val in Mflare_indices]
    Mns_flare=[ns_flare[val] for val in Mflare_indices]

    plt.plot(Mns_py, Mns_flare, 'bo')
    plt.plot([-90, 90], [-90, 90])
    plt.title('M-class flares NS')  
    plt.xlabel('NS location based on algorithm')
    plt.ylabel('NS location based on NOAA catalog')
    plt.xlim([-90,90])
    plt.ylim([-90,90])
    plt.axis([-50, 50, -50, 50])

    plt.draw()
        
    plt.figure(5)
    plt.subplot(121)
    Cew_py=[ew_py[val] for val in Cflare_indices]
    Cew_flare=[ew_flare[val] for val in Cflare_indices]    
    plt.plot(Cew_py, Cew_flare, 'bo')
    plt.plot([-90, 90], [-90, 90])
    plt.title('C-class flares EW')  
    plt.xlabel('EW location based on algorithm')
    plt.ylabel('EW location based on NOAA catalog')
    plt.xlim([-90,90])
    plt.ylim([-90,90])
    plt.axis([-90, 90, -90, 90])
    
    plt.subplot(122)
    Cns_py=[ns_py[val] for val in Cflare_indices]
    Cns_flare=[ns_flare[val] for val in Cflare_indices]    
    plt.plot(Cns_py, Cns_flare, 'bo')
    plt.plot([-90, 90], [-90, 90])
    plt.title('C-class flares NS')  
    plt.xlabel('NS location based on algorithm')
    plt.ylabel('NS location based on NOAA catalog')
    plt.xlim([-90,90])
    plt.ylim([-90,90])
    plt.axis([-50, 50, -50, 50])
    
    plt.draw()



    return arcmin_vals, flare_loc
    
    
    
#+
# Project     : SOHO - CDS     
#                   
# Name        : HEL2ARCMIN()
#               
# Purpose     : Compute position relative to sun centre from heliographic.
#               
# Explanation : Using the input heliographic coordinates of a feature,
#               calculate the position in arcmin relative to the sun centre
#               taking account of the sun's orientation (B0).  The current
#               date is assumed unless specified.  West and  North are 
#               considered positive.  
#
#   You can rely on PB0R to calculate the solar position, 
#   B angle, etc. from SOHO or Earth's vantage points, or 
#   specify complete observer coordinates in the Heliographic 
#   Spherical coordinate system:  B0, L0, R0.  You may also
#   specify the P angle.  
#               
# Use         : IDL> print, hel2armin(ns, ew, date = dat)
#                eg  print,hel2arcmin('S34','E23')
#                or  print,hel2arcmin(-34,-23)
#    or  xy = hel2arcmin(ns,ew,visible,date=dat)
#    
# Inputs      : ns      -  the Heliographic latitude in degrees (can be a
#                          string with N/S first character instead of sign).
#               ew      -  the Heliographic longitude in degrees (can be a
#                          string with E/W first character instead of sign).
#               
# Opt. Inputs : None
#               
# Outputs     : Function returns the (x,y) location in arcmins relative to
#               sun disk centre.
#               
# Opt. Outputs: If mentioned, the VISIBLE parameter gets a boolean
#   array indicating whether each point is in front of
#   theSun.
#               
# Keywords    : date    -  the date to use in the calculation of B0.
#               error   -  Output keyword containing error message;
#                          a null string is returned if no error occurs
#               soho    -  if set uses the SOHO view point rather than
#                          the Earth.  Note this functionality is
#                          duplicated by the system variable SC_VIEW,
#                          which in turn is set by the procedures
#                          USE_EARTH_VIEW or USE_SOHO_VIEW.
#
#   B0  -  The B angle, in degrees
#   P       -  The P angle, in degrees
#   R0      -  The distance of the observer from the Sun, 
#        in solar radii (use "zunits" to convert between
#        solar radii and, say, kilometers)
#   L0      -  The longitude of the observer, relative to Earth,
#        in degrees.
#   
#
# Calls       : PB0R
#               ANYTIM2UTC
#
# Restrictions: None
#
# Side effects: None
#               
# Category    : Utilities, coordinates.
#               
# Prev. Hist. : Yohkoh routine by Hudson/Wuelser.
#
# Written     : CDS version, C D Pike, RAL, 6 Sept 93
#                converted to Python by A A Reinard 30 June 2014
#               
#-            

def hel2arcmin(ns, ew,  date, b0_kw=-99, l0_kw=-99, p_kw=-99, r_kw=-99):
    
    
    if type(date) !=datetime:
        date=datetime.strptime(date, '%Y%m%d%H%M')
    radeg=180./math.pi
    
    if (type(ns) == 'str'):
        n = float(ns[1:4])
        if ns[0].upper() == 'S':
            n = -n
        w = float(ew[1:4])
        if ew[0].upper() == 'E':
            w = -w
    else:
        n = ns
        w = ew

#
#  convert to radians

    lon = w/radeg
    colat = (90. - n)/radeg

#
# get B0 and solar radius, if necessary
#
    if (r_kw == -99 or p_kw == -99 or b0_kw ==-99): 

       angles = pb0r(date)
    else:
       angles = [0, 0, 0]

#    print "angles", angles
#
# Allow keywords to override individual angles
#
    if (p_kw != -99):
        angles[0] = p_kw
    if (b0_kw != -99):
        angles[1] = b0_kw
    if (r_kw != -99):
        angles[2] = radeg*60 * math.atan(1.0/r_kw) # radians-to-arcmin conversion
    sunr = angles[2]                     

#
# vect is the (x,y,z) location of the point for b0 = 0, where x is in the
# direction of Texas, y is west, and z is north. vect1 is rotated by b0. 
#   

#  calculate the result
#
    b0 = angles[1]/radeg
    scl = math.sin(colat)
    ccl = math.cos(colat)
    cb0 = math.cos(b0)
    sb0 = math.sin(b0)
    sl = math.sin(lon)
    cl = math.cos(lon)
       
    answer=[scl*sl*sunr, (-scl * cl * sb0  +  ccl * cb0 ) * sunr]

    return answer



    
    