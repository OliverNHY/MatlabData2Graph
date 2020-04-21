::filtLine # ffe
::将ffe文件中不包含#符号的行重定向到同名txt文件
@echo off
::cls
setlocal EnableDelayedExpansion
set  Currentdir=%cd%
echo+
echo excute filtLine cutStr fileFormat
echo ____________________________

set str=%1
set FileEx=%2
set /A step=1
set /A countFile=0

if exist .\filtLineV2_txt (
	rd/Q/S .\filtLineV2_txt )
md filtLineV2_txt 
FOR /F "delims=" %%F IN ('dir /b/a-d *.!FileEx!') DO (
set longName=%%F
set "shortName=!longName:.%FileEx%=!"
set "newName=.\filtLineV2_txt\!shortName!.txt"
find /V "!str!" %%F > "!newName!"
echo "%%F" to !newName! 
set /A countFile+=1
)

echo+
echo Deal %countFile% %FileEx% file to filtLineV2_txt.
echo ____________________________
