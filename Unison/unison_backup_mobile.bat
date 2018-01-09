@echo off

REM A Script that tries to start unison from the target USB drive
REM to backup c:\BMM\Archive

set sourcedir=c:\BMM\
set targetdir=m:\BMM\
set targetunison=m:\BMM\programs\Unison\Unison-2.40.61_Text.exe

IF NOT EXIST m:\BMM\programs\Unison\NUL GOTO NODIRPRESENT
  echo Found the mobile backup drive
  echo Attempting unison synchronisation:
  echo  %sourcedir%  to  %targetdir%
  %targetunison% HHArchive_backup_mobiledisk  -batch

  GOTO EOF

:NODIRPRESENT

echo BACKUP DISK NOT PRESENT!

:EOF