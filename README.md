# bb_proc
End to end baseband processing for DSN data


## Usage 

    usage: bb_proc.py [-h] -dm DM -nc NCHAN [-m MEMLIM] [-nt NTHREAD] 
                      [-rt RFIDEC] [-snr SNRMIN] [-ezap EDGEZAP] [-nsub NSUB]
                      [-w WIDTH] [-mw MAXWIDTH] [--zerodm] [--badblocks] [-tel TEL]
                      csdir basename outdir
    
    Pipeline to process and search baseband data
    
    positional arguments:
      csdir                 Directory containing *.cs files
      basename              Base name of cs file: {basename}*.cs
      outdir                Output data directory
    
    optional arguments:
      -h, --help            show this help message and exit
      -dm DM, --dm DM       DM for intra-channel coherent dedispersion
      -nc NCHAN, --nchan NCHAN
                            Total number of output channels in filterbank
      -m MEMLIM, --memlim MEMLIM
                            Max memory to use during DADA conversion in GB (def: 16)
      -nt NTHREAD, --nthread NTHREAD
                            Number of threads for digifil processing (def: 1)
      -rt RFIDEC, --rfidec RFIDEC
                            Number of samples to decimate filfile for RFI 
                            diagnostics (def: 512, to skip this: -1)
      -snr SNRMIN, --snrmin SNRMIN
                            Minimum SNR for single pulse search (def: 6.0)
      -ezap EDGEZAP, --edgezap EDGEZAP
                            Subband edge channels to zap (def: 2)
      -nsub NSUB, --nsub NSUB
                            Number of subbands in data (def: 1)
      -w WIDTH, --width WIDTH
                            Output size for candidate snippets in sec (def: 0.1)
      -mw MAXWIDTH, --maxwidth MAXWIDTH
                            Max boxcar width in ms for SP search (def: 10)
      --zerodm              Apply the zero DM filter when dedispersing
      --badblocks           Ignore bad blocks in single pulse search
      -tel TEL, --tel TEL   DSN Telescope Name GS/RO/CN (def: RO)
