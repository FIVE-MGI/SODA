; Inno Setup script for the SODA desktop app (the same shape as FIVE Atlas's installer.iss).
; Compiled by app\build_release.ps1 when Inno Setup 6 is installed; otherwise only the zip is made.
; By hand:  ISCC.exe /DMyAppVersion=0.1.0 /DStageDir=<release>\SODA /DOutDir=<release> installer.iss

#ifndef MyAppVersion
  #define MyAppVersion "0.0.0-dev"
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
AppId={{3B0E6C55-41A7-4F7E-9C1D-5D2B8E0A6F31}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppVerName={#MyAppName} {#MyAppVersion}
DefaultDirName={autopf}\SODA
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
OutputDir={#OutDir}
OutputBaseFilename=SODA-{#MyAppVersion}-setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
ArchitecturesAllowed=x64compatible
; Per-user: no administrator rights, so it installs on a managed lab computer without IT.
PrivilegesRequired=lowest
UninstallDisplayIcon={app}\{#MyAppExeName}

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; GroupDescription: "Additional shortcuts:"

[Files]
Source: "{#StageDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Start {#MyAppName}"; Flags: nowait postinstall skipifsilent

[UninstallDelete]
; Only the program folder. Prepared data lives where the user put it, and the logs in
; %LOCALAPPDATA%\SODA are left behind on purpose.
Type: filesandordirs; Name: "{app}\_internal"
