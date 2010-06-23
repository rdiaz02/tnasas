#!/usr/bin/python

## All this code is copyright Ramon Diaz-Uriarte. For security reasons, this is for
## now confidential. No license is granted to copy, distribute, or modify it.
## Once everything is OK, it will be distributed under the GPL.

import sys
import os
import cgi 
import types
import time
import shutil
import string
import signal
import re
import glob
import tarfile

import cgitb
cgitb.enable() ## zz: eliminar for real work?
sys.stderr = sys.stdout ## eliminar?

R_MAX_time = 4 * 3600 ## 4 hours is max duration allowd for any process

## For redirections, from Python Cookbook

def getQualifiedURL(uri = None):
    """ Return a full URL starting with schema, servername and port.
    
    *uri* -- append this server-rooted uri (must start with a slash)
    """
    schema, stdport = ('http', '80')
    host = os.environ.get('HTTP_HOST')
    if not host:
        host = os.environ.get('SERVER_NAME')
        port = os.environ.get('SERVER_PORT', '80')
        if port != stdport: host = host + ":" + port
        
    result = "%s://%s" % (schema, host)
    if uri: result = result + uri
    
    return result

def getScriptname():
    """ Return te scriptname part of the URL."""
    return os.environ.get('SCRIPT_NAME', '')


# def getPathinfo():
#     """ Return the remaining part of the URL. """
#     pathinfo = os.environ.get('PATH_INFO', '')
#     return pathinfo

def clean_for_PaLS(file_in, file_out):
    """ Make sure no file has two consecutive lines that start with '#',
    so there are not lists without genes."""
    f1 = open(file_in, mode = 'r').readlines()
    f2 = open(file_out, mode = 'w')
    maxi = len(f1) - 1
    i = 0
    if len(f1) == 0:
        f2.close()
    else:
        tmp1 = f1[i]
        tmp2 = ' '
        while True:
            if i == maxi:
                break
            tmp2 = f1[i + 1]
            if not tmp1.startswith('#'):
                f2.write(tmp1)
            elif not tmp2.startswith('#'):
                f2.write(tmp1)
            tmp1 = tmp2
            i += 1
    ### make sure last one is written if not a "#"
        if not tmp2.startswith('#'):
            f2.write(tmp2)
        f2.close()


def printPalsURL(newDir,
                 tmpDir,
                 application_url = "http://tnasas.bioinfo.cnio.es",
                 f1 = "Selected.genes.txt",
                 f2 = "Selected.and.CV.selected.txt",
                 s1 = "genes selected in main run",
                 s2 = "genes selected in main run and in CV runs"):
    """ Based on Pomelo II's Send_to_Pals.cgi."""
    f=open(tmpDir + "/idtype")
    idtype = f.read().strip()
    f.close()
    f=open(tmpDir + "/organism")
    organism = f.read().strip()
    f.close()
    if (idtype != "None" and organism != "None"):
        url_org_id = "org=" + organism + "&idtype=" + idtype + "&"
    else:
        url_org_id = ""
    gl_base = application_url + '/tmp/' + newDir + '/'
    gl1 = gl_base + f1
    gl2 = gl_base + f2
    clean_for_PaLS(tmpDir + '/' + f1, tmpDir + '/' + f1)
    clean_for_PaLS(tmpDir + '/' + f2, tmpDir + '/' + f2)
    outstr0 = '<br /> <hr> ' + \
              '<h3> Send results to <a href = "http://pals.bioinfo.cnio.es">' + \
              '<IMG BORDER="0" SRC="../../palsfavicon40.png" align="middle"></a></h3>'
    outstr = outstr0 + \
             '<p> Send set of <a href="http://pals.bioinfo.cnio.es?' + \
             url_org_id + 'datafile=' + gl1 + \
             '">' + s1 + ' to PaLS</a></p>' + \
             '<p> Send set of <a href="http://pals.bioinfo.cnio.es?' + \
             url_org_id + 'datafile=' + gl2 + \
             '">' + s2 + ' to PaLS</a></p>' 
    return(outstr)



def getBaseURL():
    """ Return a fully qualified URL to this script. """
    return getQualifiedURL(getScriptname())


def commonOutput():
    print "Content-type: text/html\n\n"
    print """
    <html>
    <head>
    <title>Tnasas results</title>
    </head>
    <body>
    """
    
## to keep executing myself:
def relaunchCGI():
    print "Content-type: text/html\n\n"
    print """
    <html>
    <head>
    """
    print '<meta http-equiv="Refresh"'
    print 'content="30; URL=' + getBaseURL() + '?newDir=' + newDir + '">'
    print '<title>Tnasas results</title>'
    print '</head> <body>'
    print '<p> This is an autorefreshing page; your results will eventually be displayed here.\n'
    print 'If your browser does not autorefresh, the results will be kept for five days at</p>'
    print '<p><a href="' + getBaseURL() + '?newDir=' + newDir + '">', 'http://tnasas.bioinfo.cnio.es/tmp/'+ newDir + '/results.html</a>.' 
    print '</p> </body> </html>'
    

## Output-generating functions
def printErrorRun():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>Tnasas results </title></head><body>\n")
    outf.write("<h1> ERROR: There was a problem with the R code </h1> \n")
    outf.write("<p>  This could be a bug on our code, or a problem  ")
    outf.write("with your data (that we hadn't tought of). Below is all the output from the execution ")
    outf.write("of the run. Unless it is obvious to you that this is a fault of your data ")
    outf.write("(and that there is no way we could have avoided the crash) ")
    outf.write("please let us know so we can fix the problem. ")
    outf.write("Please send us this URL and the output below</p>")
    ## xx: eliminar, for production, from here to xxx
###     outf.write("<p> This is the output from the R run:<p>")
###     outf.write("<pre>")
###     outf.write(cgi.escape(soFar))
###     outf.write("</pre>")
    ## xxx
    outf.write("<p> This is the results file:<p>")
#    outf.write("<pre>")
    outf.write(resultsFile)
#    outf.write("</pre>")
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


def printOKRun():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">')
    outf.write('\n<html><head>')
    outf.write('\n <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-15">')
##    outf.write('\n <SCRIPT type="text/javascript" SRC="../../aqtree3clickable.js"></SCRIPT> ')
##    outf.write('\n <LINK REL="stylesheet" HREF="../../aqtree3clickable.css"> ')
    outf.write('\n <LINK REL="stylesheet" HREF="../../style1.css"> ')
    outf.write("<html><head><title>Tnasas results </title></head><body>\n")
    outf.write("<h2>Tnasas results </h2>")
    outf.write('<h3>CV error rate vs. number of genes used for classification</h3>')
    outf.write('<IMG BORDER="0" SRC="./predictor_error_rates.png">') 
    outf.write("<br /><br /> <hr>")
    outf.write('<br /><br /><h2> Results <a href="http://tnasas.bioinfo.cnio.es/help/tnasas-help.html#resultstext">(help)</a></h2> \n')
    outf.write(resultsFile)
    ## compress all the results
    allResults = tarfile.open(tmpDir + '/all.results.tar.gz', 'w:gz')
    allResults.add(tmpDir + '/results.txt', 'results.html')
    allResults.add(tmpDir + '/predictor_error_rates.png', 'predictor_error_rates.png')
##    os.system('html2text -width 200 -nobs -o correlationMatrixClusters.txt correlationMatrixCluters.html')
    allResults.close()
    outf.write('<hr> <a href="http://tnasas.bioinfo.cnio.es/tmp/' +
               newDir + '/all.results.tar.gz">Download</a> all figures and text results.')  
    try:
        outf.write(printPalsURL(newDir, tmpDir))
    except:
        None
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


def printRKilled():
    Rresults = open(tmpDir + "/results.txt")
    resultsFile = Rresults.read()
    outf = open(tmpDir + "/pre-results.html", mode = "w")
    outf.write("<html><head><title>Tnasas results </title></head><body>\n")
    outf.write("<h1> ERROR: R process killed </h1> \n")
    outf.write("<p>  The R process lasted longer than the maximum  allowed time, ")
    outf.write(str(R_MAX_time))
    outf.write(" seconds,  and was killed.")
###     outf.write("<p> This is the output from the R run:<p>")
###     outf.write("<pre>")
###     outf.write(cgi.escape(soFar))
###     outf.write("</pre>")
    outf.write("<p> This is the results file:<p>")
    outf.write("<pre>")
    outf.write(cgi.escape(resultsFile))
    outf.write("</pre>")
    outf.write("</body></html>")
    outf.close()
    Rresults.close()
    shutil.copyfile(tmpDir + "/pre-results.html", tmpDir + "/results.html")


    
## Changing to the appropriate directory
    
form = cgi.FieldStorage()
if form.has_key('newDir'):
   value=form['newDir']
   if type(value) is types.ListType:
       commonOutput()
       print "<h1> ERROR </h1>"    
       print "<p> newDir should not be a list. </p>"
       print "<p> Anyone trying to mess with it?</p>"
       print "</body></html>"
       sys.exit()
   else:
       newDir = value.value
else:
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> newDir is empty. </p>"
    print "</body></html>"
    sys.exit()

if re.search(r'[^0-9]', str(newDir)):
## newDir can ONLY contain digits.
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> newDir does not have a valid format. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()
    
redirectLoc = "/tmp/" + newDir
tmpDir = "/http/tnasas/www/tmp/" + newDir

if not os.path.isdir(tmpDir):
    commonOutput()
    print "<h1> ERROR </h1>"    
    print "<p> newDir is not a valid directory. </p>"
    print "<p> Anyone trying to mess with it?</p>"
    print "</body></html>"
    sys.exit()
    

## Were we already done in a previous execution?
## No need to reopen files or check anything else. Return url with results
## and bail out.
if os.path.exists(tmpDir + "/natural.death.pid.txt") or os.path.exists(tmpDir + "/killed.pid.txt"):
    print 'Location: http://tnasas.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'
    sys.exit()

## No, we were not done. Need to examine R output
Rrout = open(tmpDir + "/f1.Rout")
soFar = Rrout.read()
Rrout.close()
finishedOK = soFar.endswith("Normal termination\n")
errorRun = soFar.endswith("Execution halted\n")

if os.path.exists(tmpDir + "/pid.txt"):
    ## do we need to kill an R process?
    if (time.time() - os.path.getmtime(tmpDir + "/pid.txt")) > R_MAX_time:
        try:
            lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
	    os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
                      '; lamhalt -H; lamwipe -H')
        except:
            None
#             os.kill(int(open(tmpDir + "/pid.txt", mode = "r").readline()),
#                 	     signal.SIGKILL)

        printRKilled()
        os.rename(tmpDir + '/pid.txt', tmpDir + '/killed.pid.txt')
        os.remove(tmpDir + '/f1.R')
        try:
            os.system("rm /http/tnasas/www/R.running.procs/R." + newDir + "*")
        except:
            None
        print 'Location: http://tnasas.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'
##                chkmpi = os.system('/http/mpi.log/adhocCheckRmpi.py Tnasas&')
        sys.exit()

if errorRun > 0:
    printErrorRun()
    os.rename(tmpDir + '/pid.txt', tmpDir + '/natural.death.pid.txt')
    os.remove(tmpDir + '/f1.R')
##    chkmpi = os.system('/http/mpi.log/adhocCheckRmpi.py Tnasas&')
    try:
        lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    except:
        None
    try:
        os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
                  '; lamhalt -H; lamwipe -H')
    except:
        None
    try:
        os.system("rm /http/tnasas/www/R.running.procs/R." + newDir + "*")
    except:
        None
    print 'Location: http://tnasas.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'


elif finishedOK > 0:
    ##zz: killing lam seems not to be working from here...
    try:
        lamenv = open(tmpDir + "/lamSuffix", mode = "r").readline()
    except:
        None
    try:
        lamkill = os.system('export LAM_MPI_SESSION_SUFFIX=' + lamenv +
                            '; lamhalt -H; lamwipe -H')
    except:
        None
    printOKRun()
    os.rename(tmpDir + '/pid.txt', tmpDir + '/natural.death.pid.txt')
    os.remove(tmpDir + '/f1.R')
    ##    chkmpi = os.system('/http/mpi.log/adhocCheckRmpi.py Tnasas&')
    try:
        os.system("rm /http/tnasas/www/R.running.procs/R." + newDir  + "*")
    except:
        None
    print 'Location: http://tnasas.bioinfo.cnio.es/tmp/'+ newDir + '/results.html \n\n'

    
else:
    ## we only end up here if: we were not done in a previous run AND no process was overtime 
    ## AND we did not just finish. So we must continue.
    relaunchCGI()
    


