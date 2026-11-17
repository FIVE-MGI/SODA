SODA {VERSION}  (Windows 10 / 11, 64-bit)

TO RUN
  Copy this SODA folder to your own computer's hard drive first (not a network
  drive: the app is a few thousand files, and from a network drive it takes
  minutes to start instead of seconds). Then double-click SODA.exe.

  A black console window opens with the address, and your browser opens the
  SODA home page. Keep that window open while you work; closing it stops SODA.

  The first time, Windows may say "Windows protected your PC" because the app
  is not code-signed. Click "More info", then "Run anyway".

WHAT IT DOES
  1. Prepare a folder of datasets: pick the folder that holds one folder per
     dataset, and an empty output folder. SODA works out what each dataset is,
     prepares and exports it into <output folder>\export, writes a log table
     into <output folder>\logs, and checks the export for sample IDs.
  2. Open viewer: the viewer opens in a new browser tab on the export folder.
     The viewer needs an internet connection (it loads its map library from
     the web); your data stays on this computer.
  3. Check a folder for sample IDs: run it on anything before you share it.

  Your data folder is only read, never changed. Nothing is uploaded.

SAMPLE IDS
  Every Prepare ends with the sample ID check. If it says NOT clean, do not
  upload or share the export folder: rename the dataset folders or labels
  (docs/DEIDENTIFICATION.md in the SODA repository) and prepare again.
  If your lab keeps a de-identification key, put it at
     %USERPROFILE%\private\soda\sample_id_key.csv
  (or set SPATIALVIZ_KEY to its path) and SODA uses it.

FOR SCRIPTS
  SODA.exe also works from a command prompt, without the home page:
     SODA.exe prepare --src <data folder> --work <output folder>
     SODA.exe view <export folder>
     SODA.exe check <folder>
     SODA.exe prepare.py detect <folder>     (any SODA script, as with python)
  SODA.exe --help lists everything.

WHERE THINGS GO
  The log of the home page (SODA.log, prepare.log) and a compiled-code cache:
     %LOCALAPPDATA%\SODA
  If SODA will not start, that folder and the console window's text are what
  to send.
