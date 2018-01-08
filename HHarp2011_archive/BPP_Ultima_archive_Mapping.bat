@echo off
REM Script to tar and gzip the Ultima data: Special version for daily archives
REM written by Christoph Thomas, COAS, 20 Jul 2011

REM set pp = "Mapping" 

set yyyy=%date:~-4,4%
set mm=%date:~-10,2%
REM set dd=%date:~-7,2%
set dd=20

REM Deal with hours between 0 and 9 correctly
set hh=%time:~0,2%
if "%time:~0,1%"==" " set hh=0%hh:~1,1%
REM set hh=16

echo Channel wSN time: %hh%
REM Create gzipped archve from files channel 1
echo copying fies...
copy c:\BMM\Mapping\temperature\Mapping\wSN\wSN_%yyyy%%mm%%dd%*.xml c:\test\
echo archiving...
bsdtar -czf c:\BMM\archive\Mapping\BPP11_HHarp_Mapping_wSN_%yyyy%%mm%%dd%.tar.gz c:\test\*.xml
REM bsdtar -czf c:\BMM\archive\Mapping\BPP11_HHarp_Mapping_wSN_%yyyy%%mm%%dd%-%hh%.tar.gz c:\BMM\Mapping\temperature\Mapping\wSN\wSN_%yyyy%%mm%%dd%*.xml
echo Delete original files
del /Q c:\test\*.xml
del /Q c:\BMM\Mapping\temperature\Mapping\wSN\wSN_%yyyy%%mm%%dd%*.xml

echo Channel wNS time: %hh%
REM Create gzipped archve from files channel 2
echo copying fies...
copy c:\BMM\Mapping\temperature\Mapping\wNS\wNS_%yyyy%%mm%%dd%*.xml c:\test\
echo archiving...
bsdtar -czf c:\BMM\archive\Mapping\BPP11_HHarp_Mapping_wNS_%yyyy%%mm%%dd%.tar.gz c:\test\*.xml
REM bsdtar -czf c:\BMM\archive\Mapping\BPP11_HHarp_Mapping_wNS_%yyyy%%mm%%dd%-%hh%.tar.gz c:\BMM\Mapping\temperature\Mapping\wNS\wSN_%yyyy%%mm%%dd%*.xml
echo Delete original files
del /Q c:\test\*.xml
del /Q c:\BMM\Mapping\temperature\Mapping\wNS\wNS_%yyyy%%mm%%dd%*.xml

echo Channel bSN time: %hh%
REM Create gzipped archve from files channel 3
echo copying fies...
copy c:\BMM\Mapping\temperature\Mapping\bSN\bSN_%yyyy%%mm%%dd%*.xml c:\test\
echo archiving...
bsdtar -czf c:\BMM\archive\Mapping\BPP11_HHarp_Mapping_bSN_%yyyy%%mm%%dd%.tar.gz c:\test\*.xml
REM bsdtar -czf c:\BMM\archive\Mapping\BPP11_HHarp_Mapping_bSN_%yyyy%%mm%%dd%-%hh%.tar.gz c:\BMM\Mapping\temperature\Mapping\bSN\wSN_%yyyy%%mm%%dd%*.xml
echo Delete original files
del /Q c:\test\*.xml
del /Q c:\BMM\Mapping\temperature\Mapping\bSN\bSN_%yyyy%%mm%%dd%*.xml

echo Channel bNS time: %hh%
REM Create gzipped archve from files channel 4
echo copying fies...
copy c:\BMM\Mapping\temperature\Mapping\bNS\bNS_%yyyy%%mm%%dd%*.xml c:\test\
echo archiving...
bsdtar -czf c:\BMM\archive\Mapping\BPP11_HHarp_Mapping_bNS_%yyyy%%mm%%dd%.tar.gz c:\test\*.xml
REM bsdtar -czf c:\BMM\archive\Mapping\BPP11_HHarp_Mapping_bNS_%yyyy%%mm%%dd%-%hh%.tar.gz c:\BMM\Mapping\temperature\Mapping\bNS\wSN_%yyyy%%mm%%dd%*.xml
echo Delete original files
del /Q c:\test\*.xml
del /Q c:\BMM\Mapping\temperature\Mapping\bNS\bNS_%yyyy%%mm%%dd%*.xml


echo Done!
