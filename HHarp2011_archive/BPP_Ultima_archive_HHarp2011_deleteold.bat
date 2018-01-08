@echo off
REM Written by Matthias Zeeman
REM Delete files older than 2 days
forfiles.exe -p"C:\BMM\archive\HHarp2011" -m*.gz -d-2 -c"cmd /c del @PATH\@FILE"
