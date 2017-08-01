set argCount=0
for %%x in (%*) do (
   set /A argCount+=1
)
set USE_OPENMP=ON
set GENERATOR="Visual Studio 14 2015 Win64"
if not "%argCount%"=="0" (
	set GENERATOR=%1
)
if "%2"=="-DUSE_OPENMP" (
	if "%3"=="OFF" (
		set USE_OPENMP=OFF
	)
)

md visual_studio
cd visual_studio
cmake -G %GENERATOR% -DUSE_OPENMP=%USE_OPENMP% ../..
cd ..