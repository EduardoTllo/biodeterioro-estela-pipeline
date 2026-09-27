@echo off
setlocal
echo ========================================================
echo   Compilando tesis con Tectonic LaTeX...
echo ========================================================

set TECTONIC_CMD=tectonic

if exist "%~dp0tectonic.exe" (
    set TECTONIC_CMD="%~dp0tectonic.exe"
) else if exist "%~dp0..\tectonic.exe" (
    set TECTONIC_CMD="%~dp0..\tectonic.exe"
)

%TECTONIC_CMD% "%~dp0main.tex"
if %ERRORLEVEL% equ 0 (
    echo.
    echo ========================================================
    echo   [EXITO] Compilacion completada: main.pdf generado.
    echo ========================================================
) else (
    echo.
    echo ========================================================
    echo   [ERROR] Ocurrio un error durante la compilacion.
    echo ========================================================
)
pause
