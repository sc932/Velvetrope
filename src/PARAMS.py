import logging

LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}
          
##################################
#### Default parameter values ####
##################################

SEQ_TYPE_DEFAULT = 'AA'
GLOBAL_SIG_LEVEL_DEFAULT = 3.0 # 3.0 is defualt
LOCAL_WINDOW_DEFAULT = 20 # 20 is default
LOCAL_SIG_LEVEL_DEFAULT = 3.0 # 3.0 is defualt
LOCAL_BRIDGE_WIDTH_DEFAULT = 3.5 # 3.5 is default
INTRON_FINDER_DEFAULT = False
INTRON_OVERLAP_DEFAULT = 0.5
INTRON_TRUTH_RATIO_DEFAULT = 0.05
PLOT_LATEX_DEFAULT = False
PLOT_CDF_DEFAULT = False
PLOT_CDF_NUM_DEFAULT = 1
PLOT_GBG_DEFAULT = False
PLOT_COMBO_DEFAULT = False
STD_OUT_DEFAULT = True
REG_EXP_OUT_DEFAULT = True
REG_EXP_RATIO_DEFAULT = 0.5
REG_EXP_GAP_DEFAULT = 5
OUTFILE_DEFAULT = ''
PLOT_TO_SCREEN_DEFAULT = False
COMB_SPEED_DEFAULT = 1
STATS_OUT_DEFAULT = False
HTML_OUT_DEFAULT = False
COMPUTE_TYPE_DEFAULT = 'C'

##################################

DNA_READING_FRAMES = 6
DNA_DATA_DEPTHS = [5,8]
AA_READING_FRAMES = 1
AA_DATA_DEPTHS = [0]

AA_TO_5BIT = { # converts AA info into a 5 bit integer
    'G':1,
    'P':2,
    'A':3,
    'V':4,
    'L':5,
    'I':6,
    'M':7,
    'C':8,
    'F':9,
    'Y':10,
    'W':11,
    'H':12,
    'K':13,
    'R':14,
    'Q':15,
    'N':16,
    'E':17,
    'D':18,
    'S':19,
    'T':20,
    '*':21
    }

DNA_TO_8BIT = { #first 5 bits store aa info, last 3 store dna info
    'GGT':1+32, # Glycine, Gly, G
    'GGC':1+64,
    'GGA':1+32+64,
    'GGG':1+128,
    'CCT':2+32, # Proline, Pro, P
    'CCC':2+64,
    'CCA':2+32+64,
    'CCG':2+128,
    'GCT':3+32, # Alanine, Ala, A
    'GCC':3+64,
    'GCA':3+32+64,
    'GCG':3+128,
    'GTT':4+32, # Valine, Val, V
    'GTC':4+64,
    'GTA':4+32+64,
    'GTG':4+128,
    'CTT':5+32, # Leucine, Leu, L
    'CTC':5+64,
    'CTA':5+32+64,
    'CTG':5+128,
    'TTA':5+32+128,
    'TTG':5+64+128,
    'ATT':6+32, # Isoleucine, Ile, I
    'ATC':6+64,
    'ATA':6+32+64,
    'ATG':7+32, # Methionine, Met, M, Start Codon
    'TGT':8+32, # Cysteine, Cys, C
    'TGC':8+64,
    'TTT':9+32, # Phenylalanine, Phe, F
    'TTC':9+64,
    'TAT':10+32, # Tyrosine, Tyr, Y
    'TAC':10+64,
    'TGG':11+32, # Tryptophan, Trp, W
    'CAT':12+32, # Histidine, His, H
    'CAC':12+64,
    'AAA':13+32, # Lysine, Lys, K
    'AAG':13+64,
    'CGT':14+32, # Arginine, Arg, R
    'CGC':14+64,
    'CGA':14+32+64,
    'CGG':14+128,
    'AGA':14+32+128,
    'AGG':14+64+128,
    'CAA':15+32, # Glutamine, Gln, Q
    'CAG':15+64,
    'AAT':16+32, # Asparagine, Asn, N
    'AAC':16+64,
    'GAA':17+32, # Glutamic Acid, Glu, E
    'GAG':17+64,
    'GAT':18+32, # Aspartic Acid, Asp, D
    'GAC':18+64,
    'TCT':19+32, # Serine, Ser, S
    'TCC':19+64,
    'TCA':19+32+64,
    'TCG':19+128,
    'AGT':19+32+128,
    'AGC':19+64+128,
    'ACT':20+32, # Threonine, Thr, T
    'ACC':20+64,
    'ACA':20+32+64,
    'ACG':20+128,
    'TAA':21+32, # Stop codon
    'TAG':21+64,
    'TGA':21+32+64
    }
    
FIVEBIT_TO_AA = {}
for keys in AA_TO_5BIT.keys():
    FIVEBIT_TO_AA[AA_TO_5BIT[keys]] = keys
    
EIGHTBIT_TO_DNA = {}
for keys in DNA_TO_8BIT.keys():
    EIGHTBIT_TO_DNA[DNA_TO_8BIT[keys]] = keys
    
def readParameters(paramFile, isFile):
    params = {
        'INPUT_FILE':'',
        'REF_GENE':'',
        'SEQ_TYPE':SEQ_TYPE_DEFAULT,
        'GLOBAL_SIG_LEVEL':GLOBAL_SIG_LEVEL_DEFAULT,
        'LOCAL_WINDOW':LOCAL_WINDOW_DEFAULT,
        'LOCAL_SIG_LEVEL':LOCAL_SIG_LEVEL_DEFAULT,
        'LOCAL_BRIDGE_WIDTH':LOCAL_BRIDGE_WIDTH_DEFAULT,
        'INTRON_FINDER':INTRON_FINDER_DEFAULT,
        'INTRON_OVERLAP':INTRON_OVERLAP_DEFAULT,
        'INTRON_TRUTH_RATIO':INTRON_TRUTH_RATIO_DEFAULT,
        'PLOT_LATEX':PLOT_LATEX_DEFAULT,
        'PLOT_CDF':PLOT_CDF_DEFAULT,
        'PLOT_CDF_NUM':PLOT_CDF_NUM_DEFAULT,
        'PLOT_GBG':PLOT_GBG_DEFAULT,
        'PLOT_COMBO':PLOT_COMBO_DEFAULT,
        'REG_EXP_OUT':REG_EXP_OUT_DEFAULT,
        'REG_EXP_RATIO':REG_EXP_RATIO_DEFAULT,
        'REG_EXP_GAP':REG_EXP_GAP_DEFAULT,
        'OUTFILE':OUTFILE_DEFAULT,
        'PLOT_TO_SCREEN':PLOT_TO_SCREEN_DEFAULT,
        'COMB_SPEED':COMB_SPEED_DEFAULT,
        'STD_OUT':STD_OUT_DEFAULT,
        'STATS_OUT':STATS_OUT_DEFAULT,
        'HTML_OUT':HTML_OUT_DEFAULT,
        'COMPUTE_TYPE':COMPUTE_TYPE_DEFAULT
        }
        
    if isFile == False:
        return params
    
    logging.debug('Opening parameter file: ' + paramFile)
    pf = open(paramFile, 'r')
    
    for line in pf:
        if line[0] != '%' and len(line) > 1:
            idd = line.split('=')[0].strip()
            val = line.split('=')[1].strip()
            if idd in params.keys():
                logging.debug('Setting ' + idd + ' to ' + val)
                if idd != 'INPUT_FILE' and idd != 'REF_GENE' and idd != 'SEQ_TYPE' and idd != 'LOG_LEVEL' and idd != 'OUTFILE':
                    params[idd] = float(val)
                else:
                    params[idd] = val
            else:
                logging.info('Could not set ' + idd + ' to ' + val)
                
    #if params['REF_GENE'] == '' or params['SEQ_TYPE'] == '':
        #logging.error('Did not set REF_GENE or SEQ_TYPE in parameter file')
        #logging.error('These are required values')
        
    return params