rem 作者：彭雷
rem Windows命令行获取当前bat文件所在目录，添加永久系统环境变量的方法
rem 请注意用管理员权限运行该批处理文件，否则会出现find命令无法识别的错误

rem @  关闭单行回显
@echo off rem 从本行开始关闭回显。一般批处理第一行都是这个

ver | find "4.0." > NUL &&  goto win_xp    
ver | find "4.10." > NUL &&  goto win_xp   
ver | find "4.90." > NUL &&  goto win_xp   
ver | find "3.51." > NUL &&  goto win_xp   
ver | find "5.0." > NUL &&  goto win_xp    
ver | find "5.1." > NUL &&  goto win_xp    
ver | find "5.2." > NUL &&  goto win_xp    
ver | find "6.0." > NUL &&  goto win7   
ver | find "6.1." > NUL &&  goto win7    
ver | find "6.2." > NUL &&  goto win7   
ver | find "6.3." > NUL &&  goto win8 

:win_xp 
rem 设置第三方库环境变量THIRD_PARTY
wmic ENVIRONMENT where "name='ALLPROJECT_THIRD_PARTY'" delete
wmic ENVIRONMENT create name="ALLPROJECT_THIRD_PARTY",username="<system>",VariableValue="%~dp03rdPart"
echo %ALLPROJECT_THIRD_PARTY%

rem 设置ALLPROJECT_TRUNK环境变量ALLPROJECT_TRUNK,表示当前的源码主干目录
wmic ENVIRONMENT where "name='ALLPROJECT_TRUNK'" delete 
wmic ENVIRONMENT create name="ALLPROJECT_TRUNK",username="<system>",VariableValue=%~dp0
echo %ALLPROJECT_TRUNK%

rem 设置SDK环境变量ALLPROJECT_SDK_INC、ALLPROJECT_SDK_LIB和ALLPROJECT_SDK_DLL，暂时没用到
rem wmic ENVIRONMENT where "name='ALLPROJECT_SDK_INC'" delete 
rem wmic ENVIRONMENT create name="ALLPROJECT_SDK_INC",username="<system>",VariableValue="%~dp0ALLPROJECT_SDK\include"
rem echo %ALLPROJECT_SDK_INC%

rem wmic ENVIRONMENT where "name='ALLPROJECT_SDK_LIB'" delete 
rem wmic ENVIRONMENT create name="ALLPROJECT_SDK_LIB",username="<system>",VariableValue="%~dp0ALLPROJECT_SDK\lib\"
rem echo %ALLPROJECT_SDK_LIB%

rem wmic ENVIRONMENT where "name='ALLPROJECT_SDK_DLL'" delete 
rem wmic ENVIRONMENT create name="ALLPROJECT_SDK_DLL",username="<system>",VariableValue="%~dp0ALLPROJECT_SDK\dll"
rem echo %ALLPROJECT_SDK_DLL%

rem 设置exe输出路径
wmic ENVIRONMENT where "name='ALLPROJECT_BIN'" delete
wmic ENVIRONMENT create name="ALLPROJECT_BIN",username="<system>",VariableValue="%~dp0OutDir"
echo %ALLPROJECT_BIN%

rem 设置临时文件输出路径
wmic ENVIRONMENT where "name='ALLPROJECT_INTDIR'" delete
wmic ENVIRONMENT create name="ALLPROJECT_INTDIR",username="<system>",VariableValue="%~dp0IntDir"
echo %ALLPROJECT_INTDIR%

goto end 

:win7
@setx ALLPROJECT_THIRD_PARTY "%~dp03rdPart" -m

@setx ALLPROJECT_TRUNK %~dp0 -m

rem @setx ALLPROJECT_SDK_INC "%~dp0ALLPROJECT_SDK\include" -m

rem @setx ALLPROJECT_SDK_LIB "%~dp0ALLPROJECT_SDK\lib" -m

rem @setx ALLPROJECT_SDK_DLL "%~dp0ALLPROJECT_SDK\dll" -m

@setx ALLPROJECT_BIN "%~dp0OutDir" -m

@setx ALLPROJECT_INTDIR "%~dp0IntDir" -m

goto end 

:win8
@setx ALLPROJECT_THIRD_PARTY "%~dp03rdPart" -m

@setx ALLPROJECT_TRUNK %~dp0 -m

rem @setx ALLPROJECT_SDK_INC "%~dp0ALLPROJECT_SDK\include" -m

rem @setx ALLPROJECT_SDK_LIB "%~dp0ALLPROJECT_SDK\lib" -m

rem @setx ALLPROJECT_SDK_DLL "%~dp0ALLPROJECT_SDK\dll" -m

@setx ALLPROJECT_BIN "%~dp0OutDir" -m

@setx ALLPROJECT_INTDIR "%~dp0IntDir" -m

goto end

:end
pause













