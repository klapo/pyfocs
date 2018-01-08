@echo off



set bl =gg
setlocal EnableDelayedExpansion
FOR %%A IN (18 19) DO (

set h=%%A
set hh=!h!
echo !h!
echo %bl%!hh!


echo Done!

)

set hh=%2
set dd=%1

echo %hh% %dd%