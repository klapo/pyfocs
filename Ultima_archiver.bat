@echo off
REM Script to tar and gzip the Ultima data
REM written by Christoph Thomas, COAS, 20 Jul 2011
REM modified by Matthias Zeeman, COAS, 23 Aug 2011
REM   added 1 minute pause
REM   added softcoded paths

REM Paths to search for data
set sourcepath = C:\BMM\HHarp2011\temperature\campaign_rev1\
set targetpath = C:\BMM\archive\HHarp2011\
set targetfile = BPP11_HHarp_campaign_rev1
set channel_1 = channel1
set channel_2 = channel2
set channel_3 = channel3
set channel_4 = channel4

set yyyy=%date:~-4,4%
set mm=%date:~-10,2%
set dd=%date:~-7,2%

REM Deal with hours between 0 and 9 correctly
set hh=%time:~0,2%
if "%time:~0,1%"==" " set hh=0%hh:~1,1%
REM set hh=16

time /T
echo Wait a minute...
REM ping of 5 sec each
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
time /T


REM Create gzipped archve from files channel 1
set outfile = "%targetpath%\%targetfile%_%channel_1%_%yyyy%%mm%%dd%-%hh%.tar.gz"
echo Channel %channel_1% time: %hh%
mkdir "%targetpath%\%channel%\"
REM echo copying fies...
REM xcopy /Q /Y "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml" "c:\tmp\%channel%"
echo files: %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
echo archiving to %outfile%...
REM bsdtar -czf %outfile% "c:\tmp\%channel%\*.xml"
bsdtar -czf %outfile% %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
echo Delete original files
REM del /Q "c:\tmp\%channel%\*.xml"
del /Q "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml"

REM Create gzipped archve from files channel 2
set outfile = "%targetpath%\%targetfile%_%channel_2%_%yyyy%%mm%%dd%-%hh%.tar.gz"
echo Channel %channel_2% time: %hh%
mkdir "%targetpath%\%channel%\"
REM echo copying fies...
REM xcopy /Q /Y "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml" "c:\tmp\%channel%"
echo files: %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
echo archiving to %outfile%...
REM bsdtar -czf %outfile% "c:\tmp\%channel%\*.xml"
bsdtar -czf %outfile% %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
echo Delete original files
REM del /Q "c:\tmp\%channel%\*.xml"
del /Q "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml"

REM Create gzipped archve from files channel 3
set outfile="%targetpath%\%targetfile%_%channel_3%_%yyyy%%mm%%dd%-%hh%.tar.gz"
echo Channel %channel_3% time: %hh%
mkdir "%targetpath%\%channel%\"
REM echo copying fies...
REM xcopy /Q /Y "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml" "c:\tmp\%channel%"
echo files: %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
echo archiving to %outfile%...
REM bsdtar -czf %outfile% "c:\tmp\%channel%\*.xml"
bsdtar -czf %outfile% %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
echo Delete original files
REM del /Q "c:\tmp\%channel%\*.xml"
del /Q "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml"


REM Create gzipped archve from files channel 4
set outfile="%targetpath%\%targetfile%_%channel_4%_%yyyy%%mm%%dd%-%hh%.tar.gz"
echo Channel %channel_4% time: %hh%
mkdir "%targetpath%\%channel%\"
REM echo copying fies...
REM xcopy /Q /Y "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml" "c:\tmp\%channel%"
echo files: %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
echo archiving to %outfile%...
REM bsdtar -czf %outfile% "c:\tmp\%channel%\*.xml"
bsdtar -czf %outfile% %sourcepath%%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml
echo Delete original files
REM del /Q "c:\tmp\%channel%\*.xml"
del /Q "%sourcepath%\%channel%\%channel%_%yyyy%%mm%%dd%%hh%*.xml"

echo Done!

echo Backup files in %targetfile%

call C:HHarp2011_archive\Unison\unison_backup.bat

time /T
echo Wait a couple seconds...
REM ping of 5 sec each
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL

call C:HHarp2011_archive\Unison\unison_backup_mobile.bat

time /T
echo Wait a couple seconds...
REM ping of 5 sec each
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL

call C:HHarp2011_archive\BPP_Ultima_archive_HHarp2011_deleteold.bat

time /T
echo Wait a couple seconds...
REM ping of 5 sec each
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
PING 127.0.0.1 -n 6 >NUL
