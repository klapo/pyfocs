@echo off
REM Writtten by Matthias Zeeman
REM This will delete all files older than 1 days in the -p path
forfiles.exe -p"C:\BMM\archive\HHarp2011" -m*.gz -d-1 -c"cmd /c del @PATH\@FILE"
