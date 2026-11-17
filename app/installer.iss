; Inno Setup script for the SODA desktop app (the same shape as FIVE Atlas's installer.iss).
; Compiled by app\build_release.ps1 when Inno Setup 6 is installed; otherwise only the zip is made.
; By hand:  ISCC.exe /DMyAppVersion=0.1.0 /DStageDir=<release>\SODA /DOutDir=<release> installer.iss
;
; What it does (docs/DESKTOP_APP.md, "The Windows installer"):
;   - installs for the current user only, into %LOCALAPPDATA%\Programs\SODA: no administrator
;     rights, no UAC prompt. /ALLUSERS on the command line installs for everyone into Program Files
;     instead (that one needs admin rights).
;   - a Start-menu entry (SODA, and Uninstall SODA); a desktop shortcut only if ticked.
;   - an uninstaller (Settings -> Apps, or the Start-menu entry) that removes the program folder
;     only: prepared data lives in folders people chose, and the app's logs and numba cache in
;     %LOCALAPPDATA%\SODA, and neither is touched.
;   - silent install for IT:  SODA-<ver>-setup.exe /VERYSILENT /SUPPRESSMSGBOXES /NORESTART
;     (add /DIR="<folder>", /NOICONS for no Start-menu entry, /MERGETASKS=desktopicon for a
;     desktop shortcut, /ALLUSERS for a machine-wide install, /LOG="<file>" for a log);
;     silent uninstall:  "<install folder>\unins000.exe" /VERYSILENT /SUPPRESSMSGBOXES /NORESTART
;
; Nothing machine-specific is compiled in: the paths below are build-time defines (where the
; staged app folder is), and Inno does not store source paths in the setup.exe.

#ifndef MyAppVersion
  #define MyAppVersion "0.0.0-dev"
#endif
; numbers only, for the setup.exe file properties (build_release.ps1 passes 0.1.0 for 0.1.0-dev)
#ifndef MyAppNumVersion
  #define MyAppNumVersion "0.0.0"
#endif
#ifndef StageDir
  #define StageDir "release\SODA"
#endif
#ifndef OutDir
  #define OutDir "release"
#endif

#define MyAppName "SODA"
#define MyAppExeName "SODA.exe"

[Setup]
; The AppId identifies SODA to Windows across versions: a new version installs over the old one
; (same folder, one Settings -> Apps entry). NEVER change it, and never reuse it for another app
; (FIVE Atlas has its own).
AppId={{3B0E6C55-41A7-4F7E-9C1D-5D2B8E0A6F31}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
VersionInfoVersion={#MyAppNumVersion}
VersionInfoProductName={#MyAppName}
VersionInfoProductTextVersion={#MyAppVersion}
VersionInfoDescription={#MyAppName} {#MyAppVersion} setup
DefaultDirName={autopf}\SODA
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
; lets /NOICONS (no Start-menu entry) work in a silent install
AllowNoIcons=yes
OutputDir={#OutDir}
OutputBaseFilename=SODA-{#MyAppVersion}-setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
ArchitecturesAllowed=x64compatible
MinVersion=10.0
; Per-user: no administrator rights, so it installs on a managed lab computer without IT.
; commandline: /ALLUSERS (machine-wide, needs admin) and /CURRENTUSER work; no question is asked.
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=commandline
UninstallDisplayIcon={app}\{#MyAppExeName}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
; unchecked: a silent install adds no desktop shortcut unless asked (/MERGETASKS=desktopicon)
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked

[InstallDelete]
; An update replaces the whole library folder, so no file of the old version is left to be loaded
; by the new one (PyInstaller's _internal is the program's own; nothing of the user's is in it).
; Only when the folder already holds a SODA.exe, so a folder of something else is never emptied.
Type: filesandordirs; Name: "{app}\_internal"; Check: FileExists(ExpandConstant('{app}\{#MyAppExeName}'))

[Files]
Source: "{#StageDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Start {#MyAppName}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Only the program folder (its Python caches can appear in _internal after a run). Prepared data
; lives where the user put it, and the logs and numba cache in %LOCALAPPDATA%\SODA are left behind
; on purpose.
Type: filesandordirs; Name: "{app}\_internal"
