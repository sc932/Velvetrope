from PARAMS import *
from VelvetropeClasses import *
import logging, sys, os # logging

def makeHTML(setOfAligns,parameters):
    logging.info('Making the HTML pages...')
    # make the directory where everything goes
    outdir = parameters['OUTFILE'] + '/html_out'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    homepage = open(outdir + '/VRoutput.html','wb')
    page = makeHomePage(parameters)
    homepage.write(page)
    homepage.close()
    
    homepage = open(parameters['OUTFILE'] + '/VRoutput.html','wb')
    page = makeRedirectPage(outdir)
    homepage.write(page)
    homepage.close()
    
    directoryContents = os.listdir(parameters['OUTFILE'] + '/images')
    
    GBGnum = -1 # for data file
    COMBOnum = 0
    CDFnum = -1 # for data file
    CDFs = []
    for i in range(len(directoryContents)):
        if directoryContents[i][0:3] == 'GBG':
            GBGnum += 1
        elif directoryContents[i][0:5] == 'COMBO':
            COMBOnum += 1
        elif directoryContents[i][0:3] == 'CDF':
            CDFnum += 1
            if directoryContents[i][3:].split('.')[1].strip() == 'png':
                CDFs.append(int(directoryContents[i][3:].split('.')[0]))
                
    CDFs.sort()
    CDFdata = []
    if CDFnum > 0:
        CDFdatafile = open(parameters['OUTFILE'] + '/images/CDFdata.txt','r')
        CDFdatafile.readline()
        for i in range(CDFnum/2):
            line = CDFdatafile.readline()
            data1 = line.split(':')[1].split('!')[0]
            data2 = line.split(':')[1].split('!')[1]
            CDFdata.append([data1,data2])
        CDFdatafile.close()
            
    GBGdata = []
    if GBGnum > 0:
        GBGdatafile = open(parameters['OUTFILE'] + '/images/GBGdata.txt','r')
        GBGdatafile.readline()
        for i in range(GBGnum/2):
            data = []
            line = GBGdatafile.readline()
            rawData = line.split(':')[1].split('!')
            for j in range(len(rawData)):
                data.append(rawData[j].strip())
            GBGdata.append(data)
        GBGdatafile.close()
            
    for i in range(COMBOnum/2):
        combopage = open(outdir + '/combo' + str(i+1) + '.html','wb')
        page = makeComboPage(i,COMBOnum/2,parameters,setOfAligns.seqOfInt.length)
        combopage.write(page)
        combopage.close()
    if COMBOnum < 1:
        combopage = open(outdir + '/combo1.html','wb')
        page = makeComboEmptyPage(parameters)
        combopage.write(page)
        combopage.close()
        
    j = 0
    for i in CDFs:
        cdfpage = open(outdir + '/CDF' + str(i) + '.html','wb')
        page = makeCDFPage(i,CDFs,parameters,CDFdata,j)
        j += 1
        cdfpage.write(page)
        cdfpage.close()
    if CDFnum == 0:
        cdfpage = open(outdir + '/CDF1.html','wb')
        page = makeCDFEmptyPage(parameters)
        cdfpage.write(page)
        cdfpage.close()
        
    for i in range(GBGnum/2):
        gbgpage = open(outdir + '/GBG' + str(i+1) + '.html','wb')
        page = makeGBGPage(i,GBGnum/2,parameters,GBGdata)
        gbgpage.write(page)
        gbgpage.close()
    if GBGnum < 1:
        gbgpage = open(outdir + '/GBG1.html','wb')
        page = makeGBGEmptyPage(parameters)
        gbgpage.write(page)
        gbgpage.close()
        
    statspage = open(outdir + '/statspage.html','wb')
    page = makeStatsPage(parameters)
    statspage.write(page)
    statspage.close()
    
def makeRedirectPage(outdir):
    page = """
    <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
    <html>
    <head>
    <title>Velvet Rope Ouput</title>
    <meta http-equiv="REFRESH" content="0;html_out/VRoutput.html"></HEAD>
    </HTML>"""
    return page
    
def makeStatsPage(parameters):
    page = getSubpageStart(parameters)
    page += '\n<h2>Statistics Page</h2>\n'
    page += '<ul><li>Density (y) vs Run Length (x)</li><li>Run Length Distribution</li></ul>\n'
    page += '<p><a href="../images/STATS_OUT.pdf"><img src="../images/STATS_OUT.png" width="800"></a></p>\n'
    page += '<p><ul>To see a higher resolution image please look at <a href="../images/STATS_OUT.pdf">the PDF</a> in the main directory\n'
    page += getEndOfHome()
    return page
    
        
def makeGBGEmptyPage(parameters):
    page = getCSS(parameters)
    page += '<h2>No Gene-by-gene plots were generated</h2>'
    page += '<p>Try running the with the -GBGplot option enabled i.e.<ul>python Velvetrope.py -GBGplot [other options] [input file]'
    page += getEndOfHome()
    return page
        
def makeCDFEmptyPage(parameters):
    page = getCSS(parameters)
    page += '<h2>No CDF plots were generated</h2>'
    page += '<p>Try running the with the -CDFplot option enabled i.e.<ul>python Velvetrope.py -CDFplot [other options] [input file]</ul>\n'
    page += 'For more CDF plots also set the -CDFnum n option i.e.<ul>python Velvetrope.py -CDFplot -CDFnum 5 [other options] [input file]</ul>\n'
    page += 'For all CDF plots set the -CDFnum 0 option i.e.<ul>python Velvetrope.py -CDFplot -CDFnum 0 [other options] [input file]'
    page += getEndOfHome()
    return page
        
def makeComboEmptyPage(parameters):
    page = getCSS(parameters)
    page += '<h2>No Combo plots were generated</h2>'
    page += '<p>Try running the with the -Comboplot option enabled i.e.<ul>python Velvetrope.py -Comboplot [other options] [input file]'
    page += getEndOfHome()
    return page
        
def makeGBGPage(n,m,parameters,GBGdata):
    page = getSubpageStart(parameters)
    page += '\n<h2>GBG #' + str(n+1) + ' ref gene: ' + parameters['REF_GENE'] + '</h2>\n'
    page += '<p>test seqs:<br><ul>\n'
    for i in range(len(GBGdata[n])):
        page += '<li>' + GBGdata[n][i] + '</li>\n'
    page += '</p></ul>\n'
    page += '<p><a href="../images/GBG' + str(n+1) + '.pdf"><img src="../images/GBG' + str(n+1) + '.png" width="800"></a></p>\n'
    page += '<p>To see a higher resolution image please look at <a href="../images/GBG' + str(n+1) + '.pdf">the PDF</a> in the main directory</p>\n'
    page += '</div>\n<div id="sidebar">\n'
    page += '<h2>Other GBG Plots</h2>\n<ul>\n'
    for i in range(m):
        page += '<li><a href="GBG' + str(i+1) + '.html">Comps ' + str(i*4) + '-' + str((i+1)*4) + '</a></li>\n'
    page += getSubpageEnd()
    return page
        
def makeCDFPage(n,CDFs,parameters,CDFdata,j):
    page = getSubpageStart(parameters)
    page += '\n<h2>CDF #' + str(n) + ' of ' + parameters['REF_GENE'] + ' against ' + CDFdata[j][1] + '</h2>\n'
    page += '<p><a href="../images/CDF' + str(n) + '.pdf"><img src="../images/CDF' + str(n) + '.png" width="800"></a></p>\n'
    page += '<p>To see a higher resolution image please look at <a href="../images/CDF' + str(n) + '.pdf">the PDF</a> in the main directory<br>\n'
    page += 'For more CDF plots also set the -CDFnum n option i.e.<ul>python Velvetrope.py -CDFplot -CDFnum 5 [other options] [input file]</ul>\n'
    page += 'For all CDF plots set the -CDFnum 0 option i.e.<ul>python Velvetrope.py -CDFplot -CDFnum 0 [other options] [input file]</ul></p>\n'
    page += '</div>\n<div id="sidebar">\n'
    page += '<h2>Other CDF Plots</h2>\n<ul>\n'
    j = 0
    for i in CDFs:
        page += '<li><a href="CDF' + str(i) + '.html">CDF #' + str(i) + ': ' + CDFdata[j][1] + '</a></li>\n'
        j += 1
    page += getSubpageEnd()
    return page
    
def makeComboPage(n,m,parameters,maxLen):
    page = getSubpageStart(parameters)
    page += '\n<h2>Comboplot #' + str(n+1) + ' residues ' + str(n*400) + '-' + str(min((n+1)*400,maxLen)) + ' of ' + parameters['REF_GENE'] + '</h2>\n'
    page += '<p><a href="../images/COMBO' + str(n+1) + '.pdf"><img src="../images/COMBO' + str(n+1) + '.png" width="800"></a></p>\n'
    page += '<p>To see a higher resolution image please look at <a href="../images/COMBO' + str(n+1) + '.pdf">the PDF</a> in the main directory</p>\n'
    page += '</div>\n<div id="sidebar">\n'
    page += '<h2>Other Combo Plots</h2>\n<ul>\n'
    for i in range(m):
        page += '<li><a href="combo' + str(i+1) + '.html">Residues ' + str(i*400) + '-' + str(min((i+1)*400,maxLen)) + '</a></li>\n'
    page += getSubpageEnd()
    return page
    
def makeHomePage(parameters):
    page = getCSS(parameters)
    page = page + getMidOfHome()
    page += '<li>' + 'INPUT_FILE' + ' = ' + str(parameters['INPUT_FILE']) + '</li>\n'
    page += '<li>' + 'REF_GENE' + ' = ' + str(parameters['REF_GENE']) + '</li>\n'
    page += '<li>' + 'SEQ_TYPE' + ' = ' + str(parameters['SEQ_TYPE']) + '</li>\n'
    page += '<li>' + 'OUTFILE' + ' = ' + str(parameters['OUTFILE']) + '</li>\n'
    page += '<li>' + 'GLOBAL_SIG_LEVEL' + ' = ' + str(parameters['GLOBAL_SIG_LEVEL']) + '</li>\n'
    page += '<li>' + 'LOCAL_SIG_LEVEL' + ' = ' + str(parameters['LOCAL_SIG_LEVEL']) + '</li>\n'
    page += '<li>' + 'LOCAL_WINDOW' + ' = ' + str(parameters['LOCAL_WINDOW']) + '</li>\n'
    page += '<li>' + 'LOCAL_BRIDGE_WIDTH' + ' = ' + str(parameters['LOCAL_BRIDGE_WIDTH']) + '</li>\n'
    page += '<li>' + 'PLOT_CDF' + ' = ' + str(parameters['PLOT_CDF']) + '</li>\n'
    page += '<li>' + 'PLOT_CDF_NUM' + ' = ' + str(parameters['PLOT_CDF_NUM']) + '</li>\n'
    page += '<li>' + 'PLOT_GBG' + ' = ' + str(parameters['PLOT_GBG']) + '</li>\n'
    page += '<li>' + 'PLOT_COMBO' + ' = ' + str(parameters['PLOT_COMBO']) + '</li>\n'
    page += '<li>' + 'STD_OUT' + ' = ' + str(parameters['STD_OUT']) + '</li>\n'
    page += '<li>' + 'STATS_OUT' + ' = ' + str(parameters['STATS_OUT']) + '</li>\n'
    page += '<li>' + 'REG_EXP_OUT' + ' = ' + str(parameters['REG_EXP_OUT']) + '</li>\n'
    page += '<li>' + 'REG_EXP_RATIO' + ' = ' + str(parameters['REG_EXP_RATIO']) + '</li>\n'
    page += '<li>' + 'REG_EXP_GAP' + ' = ' + str(parameters['REG_EXP_GAP']) + '</li>\n'
    page += '<li>' + 'PLOT_LATEX' + ' = ' + str(parameters['PLOT_LATEX']) + '</li>\n'
    page += '<li>' + 'INTRON_FINDER' + ' = ' + str(parameters['INTRON_FINDER']) + '</li>\n'
    page += '<li>' + 'INTRON_OVERLAP' + ' = ' + str(parameters['INTRON_OVERLAP']) + '</li>\n'
    page += '<li>' + 'INTRON_TRUTH_RATIO' + ' = ' + str(parameters['INTRON_TRUTH_RATIO']) + '</li>\n'
    page += '<li>' + 'COMB_SPEED' + ' = ' + str(parameters['COMB_SPEED']) + '</li>\n'
    page += '<li>' + 'HTML_OUT' + ' = ' + str(parameters['HTML_OUT']) + '</li>\n'
    page += getEndOfHome()
    return page
    
def getSubpageEnd():
    string = """</ul>
    </div>
    <div id="footer">
        <p>Velvetrope (C) 2009, 2010 Scott Clark</p>
    </div>

</div>
</body>
</html>"""
    return string
    
def getSubpageStart(parameters):
    string = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
        "http://www.w3.org/TR/html4/strict.dtd">
<html lang="en">
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <title>Velvetrope Output HTML</title>
    <meta name="description" content="Velvetrope output">
    <meta name="copyright" content="Scott Clark 2010">
    <meta name="author" content="Scott Clark">
    <style type="text/css" media="screen, print, projection">
    body,
    html {
        margin:0;
        padding:0;
        color:#000;
        background:#a7a09a;
    }
    #wrap {
        width:1040px;
        margin:0 auto;
        background:#99c;
    }
    #header {
        padding:5px 10px;
        background:#ddd;
    }
    h1 {
        margin:0;
    }
    #nav {
        padding:5px 10px;
        background:#c99;
    }
    #nav ul {
        margin:0;
        padding:0;
        list-style:none;
    }
    #nav li {
        display:inline;
        margin:0;
        padding:0;
    }
    #main {
        float:right;
        width:800px;
        padding:10px;
        background:#9c9;
    }
    h2 {
        margin:0 0 1em;
    }
    #sidebar {
        float:left;
        width:200px;
        padding:10px;
        background:#99c;
    }
    #footer {
        clear:both;
        padding:5px 10px;
        background:#cc9;
    }
    #footer p {
        margin:0;
    }
    * html #footer {
        height:1px;
    }
    </style>

</head>
<body>
<div id="wrap">
    <div id="header"><h1>Velvetrope output run: """+str(parameters['OUTFILE'])+"""</h1></div>
    <div id="nav">
        <ul>
            <li><a href="VRoutput.html">Home</a></li>
            <li><a href="combo1.html">Comboplots</a></li>
            <li><a href="GBG1.html">Gene-by-gene plots</a></li>
            <li><a href="CDF1.html">CDF plots</a></li>
            <li><a href="statspage.html">Stats plots</a></li>
            <li><a href="../VRstdOut.txt">Standard output</a></li>
            <li><a href="../logging.out">Log file</a></li>
        </ul>
    </div>
    <div id="main">"""
    return string
    
def getMidOfHome():
    string = """        <h2>Choose output to view</h2>
        <ul>
            <li><a href="VRoutput.html">Home</a></li>
            <li><a href="combo1.html">Comboplots</a></li>
            <li><a href="GBG1.html">Gene-by-gene plots</a></li>
            <li><a href="CDF1.html">CDF plots</a></li>
            <li><a href="statspage.html">Stats plots</a></li>
            <li><a href="../VRstdOut.txt">Standard output</a></li>
            <li><a href="../logging.out">Log file</a></li>
        </ul>
        <h2>Parameters run:</h2>
        <p>
        <ul>
        """
    return string
    
def getEndOfHome():
    string = """        </ul></p>
    </div>
    <div id="sidebar">
        <h2>Navigate sequence</h2>
        <ul>
            <li><a href="VRoutput.html">Home</a></li>
        </ul>
    </div>
    <div id="footer">
        <p>Velvetrope (C) 2009, 2010 Scott Clark</p>
    </div>

</div>
</body>
</html>"""
    return string
    
def getCSS(parameters):
    string = """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
        "http://www.w3.org/TR/html4/strict.dtd">
<html lang="en">
<head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
    <title>Velvetrope Output HTML</title>
    <meta name="description" content="Velvetrope output">
    <meta name="copyright" content="Scott Clark 2010">
    <meta name="author" content="Scott Clark">
    <style type="text/css" media="screen, print, projection">
    body,
    html {
        margin:0;
        padding:0;
        color:#000;
        background:#a7a09a;
    }
    #wrap {
        width:1040px;
        margin:0 auto;
        background:#99c;
    }
    #header {
        padding:5px 10px;
        background:#ddd;
    }
    h1 {
        margin:0;
    }
    #nav {
        padding:5px 10px;
        background:#c99;
    }
    #nav ul {
        margin:0;
        padding:0;
        list-style:none;
    }
    #nav li {
        display:inline;
        margin:0;
        padding:0;
    }
    #main {
        float:right;
        width:800px;
        padding:10px;
        background:#9c9;
    }
    h2 {
        margin:0 0 1em;
    }
    #sidebar {
        float:left;
        width:200px;
        padding:10px;
        background:#99c;
    }
    #footer {
        clear:both;
        padding:5px 10px;
        background:#cc9;
    }
    #footer p {
        margin:0;
    }
    * html #footer {
        height:1px;
    }
    </style>

</head>
<body>
<div id="wrap">
    <div id="header"><h1>Velvetrope output run: """+str(parameters['OUTFILE'])+"""</h1></div>
    <div id="nav">
        <ul>
            <li><a href="VRoutput.html">Home</a></li>
            <li><a href="combo1.html">Comboplots</a></li>
            <li><a href="GBG1.html">Gene-by-gene plots</a></li>
            <li><a href="CDF1.html">CDF plots</a></li>
            <li><a href="statspage.html">Stats plots</a></li>
            <li><a href="../VRstdOut.txt">Standard output</a></li>
            <li><a href="../logging.out">Log file</a></li>
        </ul>
    </div>
    <div id="main">"""
    return string