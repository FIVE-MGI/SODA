SODA {VERSION}  (macOS, {ARCH})

WHAT IT RUNS ON
  An Apple Silicon Mac (M1, M2, M3, M4 ...) on macOS {MIN_MACOS} or newer.
  Check: Apple menu -> About This Mac -> "Chip". If it says Intel, this build
  will not open.
  Download the .dmg ON THE MAC ITSELF. Passing it through a Windows computer or
  a shared drive can strip the parts macOS needs and you will get "damaged".

INSTALL
  Double-click the .dmg, drag SODA onto the Applications folder in that window,
  then eject the disk image. Always start it from Applications.

FIRST LAUNCH (only if this copy is not notarised: macOS will refuse it once)
  It will say SODA "is damaged and can't be opened", or "Apple could not verify
  SODA is free of malware". That is Gatekeeper reacting to an app that did not
  come through the App Store, not a fault in the app. The fix, once:
     1. Press Cmd+Space, type  Terminal  and press Return.
     2. Paste this one line exactly and press Return:
            xattr -dr com.apple.quarantine /Applications/SODA.app
        (nothing is printed when it worked)
     3. Open SODA again.
  On macOS 15 Sequoia you can instead open System Settings -> Privacy &
  Security, scroll down to where it says SODA was blocked, and press
  "Open Anyway".

RUNNING IT
  SODA has no window of its own. Open it and, after a few seconds (up to a
  minute the very first time, while macOS checks it), the SODA home page opens
  in your browser at  http://127.0.0.1:8770/app/  (if no tab appears, type that
  address yourself). There is no Dock icon. To stop it, press "Quit SODA" on
  the home page; with no page open it also stops by itself after 30 minutes.

WHAT IT DOES
  1. Prepare a folder of datasets: pick the folder that holds one folder per
     dataset, and an empty output folder. SODA works out what each dataset is,
     prepares and exports it into <output folder>/export, writes a log table
     into <output folder>/logs, and checks the export for sample IDs.
  2. Open viewer: the viewer opens in a new browser tab on the export folder.
     The first time, it needs an internet connection (it loads its map library
     from the web, unless this copy of SODA includes it); after that the browser
     keeps the library. Your data stays on this Mac.
  3. Check a folder for sample IDs: run it on anything before you share it.
  Your data folder is only read, never changed. Nothing is uploaded.

  Data on the lab share: connect to it in Finder first (Go -> Connect to
  Server), then pick the folder under Locations. A typed path must look like
  /Volumes/... (not smb://...).

SAMPLE IDS
  Every Prepare ends with the sample ID check. If it says NOT clean, do not
  upload or share the export folder: rename the dataset folders or labels
  (docs/DEIDENTIFICATION.md in the SODA repository) and prepare again.
  A de-identification key, if your lab keeps one, goes at
     ~/private/soda/sample_id_key.csv

WHERE THINGS GO
  SODA.log, prepare.log and a compiled-code cache live in
     ~/Library/Application Support/SODA
  (Finder: Shift+Cmd+G and paste that path.) If SODA will not start, that log
  is what to send.
