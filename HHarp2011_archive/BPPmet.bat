@echo off
REM Script to copy and gzip the downloaded data from CR23x weatherstation at BPP
REM written by Christoph Thomas, COAS, 05 Sep 2010

REM Copy current file into new directory
copy c:\Campbellsci\Loggernet\BPP_weatherstation\BPPmet.dat d:\BPPmet\curr\BPPmet.dat
REM Deal with hours between 0 and 9 correctly
set hh=%time:~0,2%
if "%time:~0,1%"==" " set hh=0%hh:~1,1%

REM Archive file as gzip file
copy d:\BPPmet\curr\BPPmet.dat d:\BPPmet\archive\BPPMet_%date:~-4,4%%date:~-10,2%%date:~-7,2%-%hh%%time:~-8,2%.dat
gzip d:\BPPmet\archive\BPPMet_%date:~-4,4%%date:~-10,2%%date:~-7,2%-%hh%%time:~-8,2%.dat

REM Copy gzipped file to ftp directory
copy d:\BPPmet\archive\BPPMet_%date:~-4,4%%date:~-10,2%%date:~-7,2%-%hh%%time:~-8,2%.dat.gz d:\BPPmet\ftpout\

REM Combine with archived data into one file
copy d:\BPPMet\all\BPP_metstation_CR23X.dat+d:\BPPmet\curr\BPPmet.dat d:\BPPmet\all\BPP_metstation_CR23x.dat

echo Done!
