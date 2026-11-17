---
name: Bug report
about: Something does not work as the docs say
title: ""
labels: bug
---

<!--
Before you send this: do not paste sample IDs, patient information, passwords, or a link
that carries a read pass (an Azure "?sv=...&sig=..." or an S3 pre-signed link). Anyone can
read an issue once the repository is public.
-->

**What happened**

<!-- One or two sentences. Paste the exact message, if there was one. -->

**What you expected**

**Steps to make it happen again**

1.
2.
3.

**Where**

- Browser and version (for example Chrome 141, Safari 18, Firefox 142):
- Operating system (for example Windows 11, macOS 15, Ubuntu 24.04; iPhone or Android and its version):
- SODA version or commit (`git log -1 --format=%h`), or the desktop app's version:
- Python version, if it went wrong while preparing (`python --version`):

**Data**

- Data type (for example Visium HD, Xenium, CosMx, MERSCOPE, Stereo-seq, GeoMx DSP, OME-TIFF / QPTIFF, Akoya folder, H&E, `.h5ad`, OME-Zarr, SpatialData):
- Does it also happen with the synthetic demo (`python examples/make_demo.py --out demo/raw`)? yes / no / not tried
- Where the data was opened from: your own computer / a `?data=` link to storage (which: Azure, AWS, Google Cloud, a web server) / a password site / the desktop app

**Log or report**

<!--
Preparing: the last 30 lines of the command's output, or the notebook cell's output.
Viewing: add &diag=1 to the address, press Copy in the panel and paste it here (it never
holds the password or an access token).
-->
