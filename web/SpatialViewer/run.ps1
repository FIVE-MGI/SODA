# Start the SpatialViewer tile server (Windows PowerShell; run from anywhere).
#
#   .\run.ps1 -Cache <folder prepare.py wrote into> [-Port 8760] [-Password <pw>] [-NoPassword]
#
# Password: -Password, or $env:SPATIALVIZ_PASSWORD, or a .password file next to server.py (gitignored).
# Python: web\SpatialViewer\.venv, else the repo-root soda-env (or an older .venv), else `python` on PATH.
# On macOS / Linux run the same thing directly:
#   python web/SpatialViewer/server.py --cache <folder> [--port 8760] [--password <pw>]
param(
    [Parameter(Mandatory = $true)][string]$Cache,
    [int]$Port = 8760,
    [string]$Password = "",
    [switch]$NoPassword
)
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$python = "python"
foreach ($venv in @("$here\.venv", "$here\..\..\soda-env", "$here\..\..\.venv")) {
    if (Test-Path "$venv\Scripts\python.exe") { $python = (Resolve-Path "$venv\Scripts\python.exe").Path; break }
}
$pyArgs = @("$here\server.py", "--cache", $Cache, "--port", $Port)
if ($Password -ne "") { $pyArgs += @("--password", $Password) }
if ($NoPassword) { $pyArgs += "--no-password" }
& $python @pyArgs
