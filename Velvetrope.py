#!/usr/bin/env python
"""Velvetrope - an algorithm to find local alignments in sequences

PLEASE READ THE ReadMe.pdf FILE!

__copyright__

__fullUsage__ 

Has been tested on 32/64-bit linux (ubuntu, using apt-get for depends) and Mac OSX (using darwintools for depends)

Depends upon (via python setuptools "easy_install" unless otherwise noted):
    biopython:
        1.50    http://www.biopython.org
    numpy
    scipy
                http://numpy.scipy.org/
    matplotlib
                http://matplotlib.sourceforge.net/
    
    for CUDA need:
    pyCuda, pyTools:
        0.93+   http://mathema.tician.de/software/pycuda
                http://pypi.python.org/pypi/pycuda
    
Non-python-module packages:
    python
        2.6     http://www.python.org/
    swig
                http://www.swig.org/
    gcc
                http://gcc.gnu.org/
    python-dev (only really need the python.h file)
                http://www.python.org/
    
    for CUDA need:
    nvcc:
        3.0     http://developer.nvidia.com/object/cuda_3_0_downloads.html
    
For benchmarker:
    blastdb - Part of BLAST package
    blastp - Part of BLAST package
                http://blast.ncbi.nlm.nih.gov/
    phmmer - Part of HMMer package
                http://hmmer.janelia.org/

Author:
    __author__
    
Version:
    __version__ of __versionDate__
"""
__version__ = '0.16'
__versionDate__ = '19 May 2010'
__usage__ = """Usage: python Velvetrope.py [-options] <inputfile>

where basic options are:
  -h : show brief help on version and usage
"""
__author__ = 'Scott Clark <sc932 at cornell dot edu>'
__description__ = 'a bitwise algorithm for finding local alignments in sequences'
__copyright__ = """
                                   Velvetrope

                                COPYRIGHT NOTICE

Scott Clark Copyright 2009,2010
"""
__fullUsage__ = """# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - #
# Velvetrope - """ + __description__ + """   #
# Velvetrope """ + __version__ + ' (' + __versionDate__ + """)                                                #
# Copyright (C) 2009,2010 """ + __author__ + """               #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - #

Usage: python Velvetrope.py [-options] <inputfile>

Where basic options are:
  -h : show brief help on version and usage
  
Input options:
  -pf <parameter file>     : valid parameter file as described in the manual
                             other options overide file parameters
  -rd <run depth (rgo)>    : rgo for read, generate and output;
                             rg for read, generate only; etc [rgo]
  -rs <reference sequence> : ref seq name (ie 'E coli K12')
                             [defaults to first in file]
  -seq <sequence type>     : sequence type (AA and DNA supported)
                             run DNA3 or DNA6 for more reading frames [AA]
                             
Computing options:
  -cType <type> : Program to use for computing clubs (python,C,CUDA) [C]
  
Output options:
  -o <output file>      : file in which output is stored [VRout+date]
  -odir <output folder> : folder in which all output is stored [/VRout+date]
  -log <logging level>  : logging output level (warning,info,debug) [info]
  
Plotting options (n>0):
  -plotscreen : will send plots to screen as well as file
  -allplots   : make all types of plots
  -CDFplot    : make CDF plot
  -CDFnum <n> : number of CDF plots to make
  -GBGplot    : make Gene By Gene plot
  -Comboplot  : make Combo (one vs all) plot
  -RegOut     : don't show regular expressions of locally aligned areas
  -Latex      : output latex scripts
  -StdOut     : output standard out, like BLAST or HMMer [default on]
  -StatsOut   : output useful statistics on run lengths
  -HTMLout    : output html pages to browse the outputs
  
Algorithm parameters (x>0.0) (n>0):
  -globSig <x> : matches in standard deviations away from the expectation
                 an offset needs to have to trigger the global filter [8.0]
  -locSig <x>  : matches in standard deviations away from the expectation
                 a local vector needs to have vs a window length previous [8.0]
  -locWin <n>  : size of the window in which the local filter searches over, 
                 also the minimum local alignment length [20]
  -locGap <x>  : how close local alignments on the same offset need to be to 
                 be combined (a multiple of window length > 0) [3.5]
  -intron      : turns the intron (repeat) finder and eliminator on
  -intRat <x>  : how much one local alignment needs to overlap another in order
                 to be considered a possible intron (repeat) [0.5]
  -intTot <x>  : how much many local alignments need to overlap another in order
                 to be considered a possible intron (of total) [0.05]
  -regRat <x>  : ratio of test sequences need to have an amino acid for it to be
                 put in the regular expression [0.5]
  -regGap <n>  : length of gap before the next regular expression is started [5]
  -combSpd <n> : speed up the global filter by looking at only every n positions

Please see the user manual for more information on what each parameter means
"""


import logging, sys, os # logging
sys.path.append('src')
import pickle # for saving information generated
from datetime import datetime # time info
import time #timing...
import matplotlib.pylab as plt
# ALSO REQUIRES NUMPY

# Velvetrope files
import PARAMS
from VelvetropeClasses import *
import clubs_output_GBG, clubs_output_CDF, clubs_output_latex
import clubs_reader_genes, clubs_generator, clubs_output_genes, intronFinder
import clubs_output_std,clubs_output_stats, clubs_output_html
import clubgen_pycuda

def RunOutput(setOfAligns,parameters):
    t0 = time.time()
    if parameters['STATS_OUT'] == True:
        clubs_output_stats.generateStatsOut(setOfAligns,parameters)
    if parameters['PLOT_CDF'] == True or parameters['PLOT_GBG'] == True or parameters['PLOT_COMBO'] == True:
        setOfAligns = clubs_output_genes.generatePlotData(setOfAligns,parameters)
    if parameters['STD_OUT'] == True:
        clubs_output_std.generateStdOut(setOfAligns,parameters)
    if parameters['PLOT_LATEX'] == True:
        clubs_output_latex.makeLatexMat(setOfAligns.aligns[0],18,10,10)
    if parameters['PLOT_CDF'] == True:
        clubs_output_CDF.GeneCDFplots(setOfAligns, parameters['PLOT_CDF_NUM'],parameters)
    if parameters['PLOT_GBG'] == True:
        clubs_output_GBG.GeneByGeneOutput(setOfAligns,parameters)
    if parameters['PLOT_COMBO'] == True:
        clubs_output_genes.GenOutput(setOfAligns, parameters)
    #if parameters['REG_EXP_OUT'] == True:
        #clubs_output_genes.RegExpOut(setOfAligns, parameters)
    if parameters['HTML_OUT'] == True:
        clubs_output_html.makeHTML(setOfAligns,parameters)
    logging.info('Generation of output took ' + str(time.time() - t0) + ' seconds.')
    if parameters['PLOT_CDF'] == True or parameters['PLOT_GBG'] == True or parameters['PLOT_COMBO'] == True or parameters['STATS_OUT'] == True:
        if parameters['PLOT_TO_SCREEN'] == True:
            logging.info('Output to screen...')
            plt.show()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print __usage__
        sys.exit(0)
        
    if sys.argv[1] == '--help' or sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--h':
        print __fullUsage__
        sys.exit(0)

    ## parse options
    outdir = ''
    odir = False
    outputFile = ''
    runDepth = 'rgo'
    level_name = ''
    
    for i in range(len(sys.argv[1:])):
        if sys.argv[i][0] == '-':
            if sys.argv[i][1:] == 'log':
                level_name = sys.argv[i+1]
            elif sys.argv[i][1:] == 'odir':
                outdir = sys.argv[i+1]
                odir = True
    
    ## Set up output
    tt = datetime.now()
    dt = tt.timetuple()
    if outdir == '':
        outdir = 'output/VRout_m' + str(dt[1]) + 'd' + str(dt[2]) + 'y' + str(dt[0]) + 'h' + str(dt[3]) + 'min' + str(dt[4])
    #if outdir not in os.listdir('.'):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(outdir + '/images'):
        os.mkdir(outdir + '/images')
        
    if level_name not in PARAMS.LEVELS.keys():
        level_name = 'info'
    level = PARAMS.LEVELS[level_name]
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=outdir + '/logging.out',
                        filemode='w')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(level)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    
    ## Finish parsing file
    
    pfile = False
    for i in range(len(sys.argv[1:])):
        if sys.argv[i][0] == '-':
            if sys.argv[i][1:] == 'pf':
                parameters = PARAMS.readParameters(sys.argv[i+1],True)
                pfile = True
    if pfile == False:
        parameters = PARAMS.readParameters(0,False)
    
    # set up the output directory
    if odir == True:
        parameters['OUTFILE'] = outdir
    if parameters['OUTFILE'] == '':
        parameters['OUTFILE'] = outdir
    
    # read in the rest of the system arguments
    for i in range(1,len(sys.argv)):
        if sys.argv[i][0] == '-':
            if sys.argv[i][1:] == 'rd':
                runDepth = sys.argv[i+1] == 'r' or 'rg' or 'rgo' or 'g' or 'go' or 'o'
                if runDepth == False:
                    logging.error('Unsuported run depth: ' + str(sys.argv[i+1]))
                    print __usage__
            elif sys.argv[i][1:] == 'rs':
                parameters['REF_GENE'] = sys.argv[i+1]
            elif sys.argv[i][1:] == 'seq':
                parameters['SEQ_TYPE'] = sys.argv[i+1]
            elif sys.argv[i][1:] == 'o':
                outputFile = sys.argv[i+1]
            elif sys.argv[i][1:] == 'cType':
                parameters['COMPUTE_TYPE'] = sys.argv[i+1]
            elif sys.argv[i][1:] == 'plotscreen':
                parameters['PLOT_TO_SCREEN'] = True
            elif sys.argv[i][1:] == 'RegOut':
                parameters['REG_EXP_OUT'] = False
            elif sys.argv[i][1:] == 'allplots':
                parameters['PLOT_CDF'] = True
                parameters['PLOT_GBG'] = True
                parameters['PLOT_COMBO'] = True
                parameters['STATS_OUT'] = True
                parameters['REG_EXP_OUT'] = True
                parameters['HTML_OUT'] = True
            elif sys.argv[i][1:] == 'CDFplot':
                parameters['PLOT_CDF'] = True
            elif sys.argv[i][1:] == 'CDFnum':
                parameters['PLOT_CDF_NUM'] = numpy.int(sys.argv[i+1])
            elif sys.argv[i][1:] == 'GBGplot':
                parameters['PLOT_GBG'] = True
            elif sys.argv[i][1:] == 'Comboplot':
                parameters['PLOT_COMBO'] = True
            elif sys.argv[i][1:] == 'Latex':
                parameters['PLOT_LATEX'] = True
            elif sys.argv[i][1:] == 'StdOut':
                parameters['STD_OUT'] = True
            elif sys.argv[i][1:] == 'StatsOut':
                parameters['STATS_OUT'] = True
            elif sys.argv[i][1:] == 'HTMLout':
                parameters['HTML_OUT'] = True
            elif sys.argv[i][1:] == 'globSig':
                parameters['GLOBAL_SIG_LEVEL'] = numpy.float(sys.argv[i+1])
            elif sys.argv[i][1:] == 'locSig':
                parameters['LOCAL_SIG_LEVEL'] = numpy.float(sys.argv[i+1])
            elif sys.argv[i][1:] == 'locWin':
                parameters['LOCAL_WINDOW'] = numpy.int(sys.argv[i+1])
            elif sys.argv[i][1:] == 'locGap':
                parameters['LOCAL_BRIDGE_WIDTH'] = numpy.float(sys.argv[i+1])
            elif sys.argv[i][1:] == 'intron':
                parameters['INTRON_FINDER'] = True
            elif sys.argv[i][1:] == 'intRat':
                parameters['INTRON_OVERLAP'] = numpy.float(sys.argv[i+1])
            elif sys.argv[i][1:] == 'intTot':
                parameters['INTRON_TRUTH_RATIO'] = numpy.float(sys.argv[i+1])
            elif sys.argv[i][1:] == 'regRat':
                parameters['REG_EXP_RATIO'] = numpy.float(sys.argv[i+1])
            elif sys.argv[i][1:] == 'regGap':
                parameters['REG_EXP_GAP'] = numpy.int(sys.argv[i+1])
            elif sys.argv[i][1:] == 'combSpd':
                parameters['COMB_SPEED'] = numpy.int(sys.argv[i+1])
            elif sys.argv[i][1:] == 'pf' or sys.argv[i][1:] == 'log' or sys.argv[i][1:] == 'odir':
                pass
            else:
                logging.error('Could not recognize option : ' + sys.argv[i])
                print __usage__
                sys.exit(1)
            i = i+1
        else:
            parameters['INPUT_FILE'] = sys.argv[i]
            
    if parameters['INPUT_FILE'] == '':
        logging.error('need input file')
        print __usage__
        sys.exit(1)
    
    if parameters['REF_GENE'] == '' and runDepth[0] == 'r':
        pf = open(parameters['INPUT_FILE'], 'r')
        for line in pf:
            if line[0] == '>':
                parameters['REF_GENE'] = line[1:-1]
                break

    
        
    
    if outputFile == '':
        outputFile = 'VelvetropeOutput_m' + str(dt[1]) + 'd' + str(dt[2]) + 'y' + str(dt[0]) + 'h' + str(dt[3]) + 'min' + str(dt[4]) + '__' + str(sys.argv[1])
        
    logging.debug('Starting with parameters ' + str(sys.argv))

    
    
    
    if runDepth[0] == 'r':
        
        [seqOfInt, testSeqs] = clubs_reader_genes.reader(parameters)
        
        if runDepth == 'rg' or runDepth == 'rgo':
            logging.info('Running generator...')
            t0 = time.time()
            setOfAligns = clubs_generator.GenClubs(seqOfInt, testSeqs, parameters)
            logging.info('Generation of clubs took ' + str(time.time() - t0) + ' seconds.')
            if parameters['INTRON_FINDER'] == True:
                setOfAligns = intronFinder.IntronFinder(setOfAligns, parameters)
        else:
            logging.info('Writing output (from read)...')
            # pickle out seqOfInt, testSeqs
            output = open(outputFile + '.vrr.pkl', 'wb')
            pickle.dump([setOfAligns],output)
    
        if runDepth == 'rgo':
            logging.info('Running output...')
            RunOutput(setOfAligns,parameters)
        else:
            logging.info('Writing output (from gen)...')
            # pickle out setOfAligns
            #setOfAligns = clubs_output_genes.generatePlotData(setOfAligns,parameters)
            output = open(outputFile + '.vrg.pkl', 'wb')
            pickle.dump([setOfAligns],output)
    
    elif runDepth[0] == 'g':
        logging.info('Reading in sequences from pickle file: ' + str(inputFile))
        logging.info('Running generator...')
        # pickle in seqOfInt, testSeqs
        pkl_file = open(inputFile,'rb')
        [seqOfInt, testSeqs] = pickle.load(pkl_file)
        t0 = time.time()
        setOfAligns = clubs_generator.GenClubs(seqOfInt, testSeqs, parameters)
        logging.info('Generation of clubs took ' + str(time.time() - t0) + ' seconds.')
        if INTRON_FINDER == True:
            setOfAligns = intronFinder.IntronFinder(setOfAligns, parameters)

        if runDepth == 'go':
            logging.info('Running output...')
            RunOutput(setOfAligns,parameters)
        else:
            logging.info('Writing output (from gen)...')
            #setOfAligns = clubs_output_genes.generatePlotData(setOfAligns,parameters)
            # pickle out setOfAligns
            output = open(outputFile + '.vrg.pkl', 'wb')
            pickle.dump([setOfAligns],output)
            
    elif runDepth[0] == 'o':
        logging.info('Reading in sequences from pickle file: ' + str(inputFile))
        # pickle in setOfAligns
        pkl_file = open(inputFile,'rb')
        [setOfAligns] = pickle.load(pkl_file)
        logging.info('Running output...')
        RunOutput(setOfAligns,parameters)
        
    else:
        logging.info('Was expecting parameters made up of r,g,o. Please see help file')
