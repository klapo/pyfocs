SET logfile="C:\Users\Ultima\DarkMix\muppet_archiver\logs\batch.log"
@echo off
@echo Archiving data: %date% %time% >> %logfile%
call C:\Users\Ultima\Anaconda3\Scripts\activate.bat
python "C:\Users\Ultima\silixa_archiver.py" >> %logfile%
@echo finished at %date% %time% >> %logfile%
%* Pause
