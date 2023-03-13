#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""

@author: pearlman (Aaron B. Pearlman; aaron.b.pearlman@caltech.edu)
Version 1: 2021.10.20

"""

scriptdir = "/src/bbsearch"

import numpy as np
import glob
import os
import copy
import sys
from subprocess import call, check_call, check_output, Popen

sys.path.append(scriptdir)

import h5py
import filterbank

import click

import logging

LOG_FORMAT = "[%(asctime)s] %(levelname)s "
LOG_FORMAT += "%(module)s::%(funcName)s():l%(lineno)d: "
LOG_FORMAT += "%(message)s"
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


BLOCKSIZE = 1e6

def read_sp_list(sp_fn):
    """ 
    Read a single pulse candidate list. 
    """
    readFile = open(sp_fn, "r")

    # Read in the header from the data file.
    headerString = readFile.readline().strip()

    # Set the delimiter.
    delimiter = " "

    # Read the single pulse list generated by PRESTO.
    dataList = readFile.read().splitlines()

    # Remove empty lines from the list.
    dataList = list(filter(None, dataList))

    numPoints = len(dataList)

    cand_dm = []
    cand_snr = []
    cand_time = []
    cand_sample = []
    cand_boxcarbinwidth = []

    for listIndex in np.arange(0, numPoints, 1):
        line = (" ".join(dataList[listIndex].split())).split(delimiter)

        cand_dm.append(float(line[0]))
        cand_snr.append(float(line[1]))
        cand_time.append(float(line[2]))
        cand_sample.append(float(line[3]))
        cand_boxcarbinwidth.append(float(line[4]))

    readFile.close()
    
    return cand_dm, cand_snr, cand_time, cand_sample, cand_boxcarbinwidth



# RUN NOTES:
# Each filterbank file must have a .fil extension. They get 
# moved to sub-directories, searched, and candidates extracted.
# Don't forget to update scriptdir at the top of the script! 
# This should coincide with the location of the code/repo.
# Example: python baseband_search.py --obsroot sband --nsub 7 --nchansub 32 --edgezap 2 -zapchans 1,2 --dm 87.77 --spthresh 6.0 --width 0.100 --telescope RO --dsn

@click.command()
@click.option("--obsroot", 
              help="Root name of candidate directory (for each filterbank).", 
              type=str)
@click.option("--nsub", 
              help="Number of subbands in baseband data.", 
              type=int)
@click.option("--nchansub", 
              help="Number of channels per subband.", 
              type=int)
@click.option("--edgezap", 
              help="Number of channels to zap at top and bottom of each subband.", 
              type=int)
@click.option("--zapchans", 
              help="Comma separated list of channels to zap. Index 0 corresponds to the highest frequency channel.", 
              type=str, default="")
@click.option("--dm", 
              help="Dispersion measure (for incoherent, inter-channel dedispersion and searching).", 
              type=float)
@click.option("--spthresh", 
              help="S/N threshold to use for single pulse search.", 
              type=float)
@click.option("--width", 
              help="Length of data to extract around each candidate burst (units: seconds).", 
              type=float)
@click.option("--telescope", help="Telescope (Madrid=RO, Goldstone=GS, Canberra=CN).", type=str)
@click.option("--dsn", 
              help="Set a dummy machine_id (999) for the DSN.", 
              is_flag=True)


def pipeline(obsroot, nsub, nchansub, edgezap, zapchans, 
             dm, spthresh, width, telescope, dsn):
    """
    Run search pipeline
    """
    # Get list of filterbank files for the baseband search.
    fil_fn = glob.glob("*.fil")
    fil_fn.sort()

    cand_dir = []
    cwd = os.getcwd()
    
    log.debug("Setting up the data & search directories...")
    
    for i in np.arange(0, len(fil_fn), 1):
        call("mkdir %s%i" % (obsroot, i+1), shell=True)
        
        if (telescope != None and dsn==True):
            fix_cmd = "python %s/fix_sigproc_header.py " %scriptdir +\
                      "--filterbank %s " %fil_fn[i] +\
                      "--telescope %s " %telescope +\
                      "--dsn"
            call(fix_cmd, shell=True)
        elif (telescope != None and dsn==False):
            fix_cmd = "python %s/fix_sigproc_header.py " %scriptdir +\
                      "--filterbank %s " %fil_fn[i] +\
                      "--telescope %s " %telescope +\
            call(fix_cmd, shell=True)
        
        call("mv %s %s%i" % (fil_fn[i], obsroot, i+1), shell=True)
        
        fil_fn[i] = "%s/%s%i/%s" % (cwd, obsroot, i+1, fil_fn[i])
        cand_dir.append("%s/%s%i" % (cwd, obsroot, i+1))
   
    zap_list = []
     
    if (edgezap > 0):
        log.debug("Generating zap channel list. " +\
                  "Edge channels will be zapped in " +\
                  "dedispersed time-series...")
        zap_list_band1 = np.arange(0, nchansub, 1, dtype='int')
        zap_list_band1 = np.concatenate((zap_list_band1[0:edgezap], zap_list_band1[nchansub-edgezap:nchansub]))
    
        # Generate channel idxs to zap. 
        for i in np.arange(0, nsub, 1, dtype='int'):
           zap_add = zap_list_band1 + (i * nchansub)
           zap_list = np.concatenate((zap_list, zap_add))
    
    else:
        log.debug("No edge channels will be zapped...")

    if (zapchans != ""):
        log.debug("Adding channels to zaplist...")
        zapchan_list = [ int(zchan) for zchan in zapchans.split(',') ]
        #print(zapchan_list)
        #print(zap_list)
        zap_list = np.concatenate( (zap_list, zapchan_list) )
    else:
        log.debug("No additional channels will be zapped...")
    
    if len(zap_list):
        Nchans = nsub * nchansub
        zap_chan_str = ""
        zap_list = np.unique(zap_list)[::-1]
       
        # Note that we are flipping the channels to 
        # get the PRESTO ordering scheme (ie, reversed 
        # order to DSN). In PRESTO, channel index 0 
        # corresponds to the lowest frequency channel.  
        # In Aaron's m_fb_zapchan.py tool, channel index 0 
        # corresponds to the highest frequency channel.
        # The user inputs channel indices for zapchans 
        # according to the m_fb_zapchan.py convention,
        # and this is converted to the PRESTO convention.
        for i in np.arange(0, len(zap_list), 1):
            zap_chan_str = zap_chan_str + "%i," %(Nchans - 1 - int(zap_list[i]))
    
        zap_chan_str = zap_chan_str[:-1]
        log.debug("Zaplist: %s" % zap_chan_str)

    else:  pass
 
    # Incoherent dedispersion and single pulse search.
    log.debug("Incoherently dedispersing and single pulse searching...")
    
    dat_fn = []
    sp_fn = []
    
    for i in np.arange(0, len(fil_fn), 1):
       os.chdir(cand_dir[i])
       dedisp_root = (fil_fn[i].split("/")[-1]).split(".fil")[0] + "_DM%.3f" % dm
       
       if len(zap_list):
          zap = "-ignorechan %s" % zap_chan_str
       else:
          zap = ""

       dm_cmd = "prepdata -filterbank %s -dm %.3f -nobary -noclip %s -o %s" % (fil_fn[i], dm, zap, dedisp_root)
       log.debug("De-dispersing...")
       log.debug(dm_cmd)

       call(dm_cmd, shell=True)
       
       dat_fn.append(cand_dir[i] + "/%s.dat" % dedisp_root)
       sp_fn.append(cand_dir[i] + "/%s.singlepulse" % dedisp_root)
       
       # Added the -d 32 option 
       call("python %s/single_pulse_search_w16ms.py %s -t %.3f -m 10.0 -d 32 -b" % (scriptdir, dat_fn[i], spthresh), shell=True)
       
    # Extract filterbank snippets of candidates.
    log.debug("Extracting FRB candidates...")
    for i in np.arange(0, len(sp_fn), 1):
        os.chdir(cand_dir[i])
        
        cand_dm, cand_snr, cand_time, cand_sample, cand_boxcarbinwidth = read_sp_list(sp_fn[i])
        
        inf_fn = glob.glob("*.inf")[0]
        tsamp = check_output("grep Width %s" % inf_fn, shell=True)
        tsamp = tsamp.decode('utf-8')
        tsamp = float(tsamp.split("=  ")[-1])
        
        width_samples = int(np.round(width / tsamp))

        fil_fn_short = fil_fn[i].split("/")[-1]
        
        for j in np.arange(0, len(cand_sample), 1):
            start_sample = cand_sample[j] - int(np.round(width_samples / 2.0))
            
            fb_cand_fn = sp_fn[i].split(".singlepulse")[0] + "_extract_cand%i.fil" % (j+1)
            fb_cand_fn_short = fb_cand_fn.split('/')[-1]

            extract_cmd = "extract %s %i %i > %s" % (fil_fn_short, start_sample, width_samples, fb_cand_fn_short)
            
            log.debug("Extracting candidate...")
            log.debug(extract_cmd)

            call(extract_cmd, shell=True)


if __name__ == "__main__":
    pipeline()
