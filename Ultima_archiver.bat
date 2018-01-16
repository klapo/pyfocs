@echo off
REM SETLOCAL ENABLEEXTENSIONS
set UNISON=%~dp0Unison\
set path=%UNISON%libs\gtk-runtime-2.16.6.0\Gtk\bin
echo %UNISON%
echo %path%
:: Script to tar and gzip the Ultima data

:: Paths to search for local data
set sourcepath=C:\Users\Karl Lapo\Desktop\
set targetpath=C:\Users\Karl Lapo\Desktop\
set targetfile=archive
:: Root variables for the unison command
set root_local=/Users/Karl Lapo/Desktop/archive/
set root_backup=M:/test/
set root_mobile=E:/test_archive/
set logfile_backup=%root_backup%_logfile.txt
set logfile_mobile=%root_mobile%_logfile.txt
:: Names of the channels
set channel_1=channel 1
set channel_2=channel 2
set channel_3=channel 3
set channel_4=channel 4

:: The below works when the script is run in archiving mode. To run on old
:: data you must specify the dates manually.
set yyyy=%date:~-4,4%
set mm=%date:~-10,2%
set dd=%date:~-7,2%

:: Deal with hours between 0 and 9 correctly
set hh=%time:~0,2%
if "%time:~0,1%"==" " set hh=0%hh:~1,1%

:: Manually specify dates below and uncomment if commented
set yyyy=2017
set hh=15
set dd=13
set mm=12

REM Create gzipped archve from files channel 1
if EXIST %sourcepath%%channel_1%(
  set outfile=%targetpath%%targetfile%\%channel_1%_%yyyy%%mm%%dd%-%hh%.tar.gz
  echo %channel_1% time: %hh%:00 %mm%/%dd%/%yyyy%
  echo files: %sourcepath%%channel_1%\%channel_1%_%yyyy%%mm%%dd%%hh%*
  echo archiving to %outfile%
  bsdtar -czf "%outfile%" "%sourcepath%%channel_1%\%channel_1%_%yyyy%%mm%%dd%%hh%*"
  echo Delete original files
  del /Q "%sourcepath%%channel_1%\%channel_1%_%yyyy%%mm%%dd%%hh%*"
)

echo Done!
echo Backup files in %targetpath%%targetfile%

REM Unison for mobile backup drive
if EXIST %root_mobile% (
  echo Found the mobile backup drive
  echo Attempting unison synchronisation:
  echo  %root_local%  to  %root_mobile%
  "%UNISON%Unison-2.40.61 Text" "%root_local%" "%root_mobile%" -logfile "%logfile_mobile%" -force "%root_local%" -batch -nodeletion "%root_mobile%"
) ELSE (
  echo MOBILE BACKUP DISK NOT FOUND
)


REM Unison for  backup drive
REM if EXIST %root_backup%\NUL(
REM   echo Found the backup drive
REM   echo Attempting unison synchronisation:
REM   echo  %sourcedir%  to  %targetdir%
REM   %~dp0Unison\Unison-2.40.61 Gtk+.exe %root_local% %root_backup% -nodeletion %root_backup% -logfile %logfile_backup% -force %force% -batch
REM ) ELSE (
REM   echo BACKUP DISK NOT PRESENT
REM )

REM time /T
REM echo Wait a couple seconds...
REM REM ping of 5 sec each
REM PING 127.0.0.1 -n 6 >NUL
REM PING 127.0.0.1 -n 6 >NUL
REM PING 127.0.0.1 -n 6 >NUL

REM REM Delete files older than 2 days
REM %~dp0forfiles.exe -p"%targetpath" -m*.gz -d-2 -c"cmd /c del @PATH\@FILE"
REM
REM time /T
REM echo Wait a couple seconds...
REM REM ping of 5 sec each
REM PING 127.0.0.1 -n 6 >NUL
REM PING 127.0.0.1 -n 6 >NUL
REM PING 127.0.0.1 -n 6 >NUL
