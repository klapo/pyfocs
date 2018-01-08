@echo off

REM A Script that tries to start unison from the target USB drive
REM to backup c:\BMM\Archive

set sourcedir=c:\BMM\
set targetdir=h:\BMM\
set targetunison=h:\BMM\programs\Unison\Unison-2.40.61_Text.exe

IF NOT EXIST h:\BMM\programs\Unison\NUL GOTO NODIRPRESENT
  echo Found the backup drive\n
  echo Attempting unison synchronisation:
  echo  %sourcedir%  to  %targetdir%
  %targetunison% HHArchive_backup  -batch

  GOTO EOF

:NODIRPRESENT

echo BACKUP DISK NOT PRESENt!

:EOF