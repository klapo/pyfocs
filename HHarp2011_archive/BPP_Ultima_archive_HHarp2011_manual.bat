@echo off
REM Script to tar and gzip the Ultima data
REM written by Christoph Thomas, COAS, 20 Jul 2011
REM modified by Matthias Zeeman, COAS, 23 Aug 2011
REM   added 1 minute pause
REM   added softcoded paths

IF "%hh%"=="" SET hh=%1
IF "%dd%"=="" SET hh=%2


set yyyy=%date:~-4,4%
set mm=%date:~-10,2%
set dd=%date:~-7,2%
REM set dd=20

set sourcepath=C:\BMM\HHarp2011\temperature\campaign_rev1\
set tmppath=c:\BMM\tmp\
set targetpath=c:\BMM\archive\HHarp2011\
set targetfile=BPP11_HHarp_campaign_rev1

REM Deal with hours between 0 and 9 correctly
REM set hh=%time:~0,2%
REM if "%time:~0,1%"==" " set hh=0%hh:~1,1%
set hh=13
set dd=02
set mm=09

REM setlocal EnableDelayedExpansion
REM for %%A in ( 20 21 22 23 17 ) DO (
REM hh=%%A

echo %hh%

time /T
echo Wait a minute...
REM ping of 5 sec each

REM Create gzipped archve from files channel 1
set channel=wSN
set outfile="%targetpath%\%targetfile%_%channel%_%yyyy%%mm%%dd%-%hh%.tar.gz"
echo Channel %channel% time: %hh%
mkdir "%targetpath%\%channel%\"
echo archiving to %outfile%...
echo Files: %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
bsdtar -czf %outfile% %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml 
echo Delete original files
del /Q "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml"


REM Create gzipped archve from files channel 2
set channel=wNS
set outfile="%targetpath%\%targetfile%_%channel%_%yyyy%%mm%%dd%-%hh%.tar.gz"
echo Channel %channel% time: %hh%
mkdir "%targetpath%\%channel%\"
echo archiving to %outfile%...
echo Files: %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
bsdtar -czf %outfile% %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml 
echo Delete original files
del /Q "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml"


REM Create gzipped archve from files channel 3
set channel=bSN
set outfile="%targetpath%\%targetfile%_%channel%_%yyyy%%mm%%dd%-%hh%.tar.gz"
echo Channel %channel% time: %hh%
mkdir "%targetpath%\%channel%\"
echo archiving to %outfile%...
echo Files: %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
bsdtar -czf %outfile% %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml 
echo Delete original files
del /Q "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml"


REM Create gzipped archve from files channel 4
set channel=bNS
set outfile="%targetpath%\%targetfile%_%channel%_%yyyy%%mm%%dd%-%hh%.tar.gz"
echo Channel %channel% time: %hh%
mkdir "%targetpath%\%channel%\"
echo archiving to %outfile%...
echo Files: %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
bsdtar -czf %outfile% %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml 
echo Delete original files
del /Q "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml"
 


echo Done!
