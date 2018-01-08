@echo off
REM Script to tar and gzip the Ultima data
REM written by Christoph Thomas, COAS, 20 Jul 2011

set pp = Mapping 

set yyyy=%date:~-4,4%
set mm=%date:~-10,2%
REM set dd=%date:~-7,2%
set dd=19

REM Deal with hours between 0 and 9 correctly
REM set hh=%time:~0,2%
REM if "%time:~0,1%"==" " set hh=0%hh:~1,1%
set hh=16

echo Channel wSN time: %hh% pp
REM Create gzipped archve from files channel 1
echo copying fies...
copy c:\BMM\pp\temperature\pp\wSN\wSN_%yyyy%%mm%%dd%%hh%*.xml c:\test\
echo archiving...
bsdtar -czf c:\BMM\archive\pp\BPP11_HHarp_%pp%_wSN_%yyyy%%mm%%dd%-%hh%.tar.gz c:\test\*.xml
REM bsdtar -czf c:\BMM\archive\test_waterbath\BPP11_HHarp_testwaterbath_wSN_%yyyy%%mm%%dd%-%hh%.tar.gz c:\BMM\test_waterbath\temperature\test_waterbath\wSN\wSN_%yyyy%%mm%%dd%%hh%*.xml
echo Delete original files
del /Q c:\test\*.xml
REM del /Q c:\BMM\test_waterbath\temperature\test_waterbath\wSN\wSN_%yyyy%%mm%%dd%%hh%*.xml
