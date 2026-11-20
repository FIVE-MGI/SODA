# Hosting your own data

SODA has two halves:

- **the page**: the files in `web/SpatialViewer/` that a browser runs (`index.html`,
  `viewer.js`, `styles.css`, `lasso_worker.js`, `umapfit.js`, `leiden.js`, `leiden_worker.js`, and your `config.js`);
- **the data**: one folder per dataset, written by `prepare.py export`
  (layout in [DATA_FORMAT.md](DATA_FORMAT.md)). It is plain files: images cut into
  tiles, JSON and small binary arrays. Nothing has to run on the server.

So hosting means putting those two sets of files somewhere a browser can fetch them.

**Looking for the easy way?** [SHARE_YOUR_DATA.md](SHARE_YOUR_DATA.md) covers the two
simplest: upload only the data and share a `?data=` link to the shared SODA page
(https://soda.mgifive.org/), or copy
SODA on GitHub to get your own site. This page is the full reference.

The options below are in order of effort.

| | Who can see it | What you run | Effort |
|---|---|---|---|
| [A. Your own computer](#a-on-your-own-computer) | you | two small Python web servers | minutes |
| [B. Any static web host](#b-any-static-web-host-public-data) | anyone with the link | nothing | an hour |
| [C. Azure, private, with a password](#c-azure-private-data-with-a-password) | people you give the password | one small function | an afternoon |
| [D. AWS (or another cloud), private](#d-aws-s3--cloudfront-or-another-backend) | people you give the password | your own token endpoint | a day, needs cloud experience |
| [E. Live tile server](#e-optional-live-tile-server-serverpy) | people on your network | `server.py` | minutes, but a machine must stay on |

Visitors need an internet connection in every case, also for option A, unless you serve the
two libraries yourself: the page loads its drawing library (deck.gl) and, for the circle tool,
umap-js from `cdn.jsdelivr.net`, or from a `vendor/` folder next to `index.html` after
`python web/SpatialViewer/prepare.py vendor` (copy that folder along with the page). The
viewer needs a browser with WebGL2. It is developed in Chrome on desktop; automated tests
run it in Playwright's Chromium, WebKit (Safari's engine) and Firefox, and the image path was
also tested on an iPhone. The same link works on phones and tablets: the page lays itself
out for the screen. What was tested where: [Browser support](../README.md#browser-support).

Before anything leaves your computer, go through the
[checklist at the end](#before-you-upload-checklist).

**Name your site after SODA.** Please put `soda` in the web address you choose, for
example `https://soda.yourlab.org/` or `https://yourlab.org/soda/`, so visitors know
which software they are using. The SODA label in the bottom corner of the page stays
whatever title you give your instance in `config.js`.

---

## A. On your own computer

This is the way to check your own export before uploading it. From the repository folder,
with the virtual environment active:

```bash
python web/SpatialViewer/static_serve.py export 8765
python -m http.server 8766 --directory web/SpatialViewer
```

Run each in its own terminal window; both keep running until you press Ctrl+C. Then open

```
http://localhost:8766/index.html?static=http://localhost:8765&ds=my_dataset
```

Point the first command at your export folder and set `ds=` to your dataset id. (On
Windows, if `python` is not found, try `py`.)

What the pieces do:

- `static_serve.py <folder> <port>` serves an export folder with the
  `Access-Control-Allow-Origin: *` header the page needs, and `Cache-Control: no-store`
  so a re-export shows up without a hard reload. It answers HTTP Range requests (`206
  Partial Content`), which the page uses to read one gene out of a gene pack; plain
  `python -m http.server` does not, so do not serve the **export** folder with it (the
  page itself is fine). It listens on `127.0.0.1` only, so nobody else on your network
  can reach it.
- `?static=<url>` tells the page to read every dataset from `<url>/<id>/` and to skip
  passwords and tokens entirely. `?ds=<id>` picks the dataset. With both, the dataset
  does not need to be listed in `config.js`, and `config.js` does not need to exist.

## B. Any static web host (public data)

Use this when the data may be public: a published dataset, a teaching set.
Anyone who has, or guesses, the address can download the files.

Copy the page files and the export folders to the host, then tell `config.js` where each
dataset lives, with `tokenEndpoint` empty:

```js
window.SpatialViewerConfig = {
  datasets: {
    my_dataset: { label: "My dataset", base: "data/my_dataset/" },                             // same site
    // other_set: { label: "...", base: "https://data.example.org/viewer/other_set/" },        // elsewhere
  },
  defaultDataset: "my_dataset",
  tokenEndpoint: "",
};
```

`base` is the address of one dataset folder, ending in `/`. A relative `base` such as
`data/my_dataset/` is read relative to the page, so the simplest layout puts the export
next to the page:

```
index.html  viewer.js  styles.css  lasso_worker.js  umapfit.js  leiden.js  leiden_worker.js  config.js  logos/
data/my_dataset/manifest.json
data/my_dataset/tiles/...
```

You do not need to publish the `.py` files, `requirements*.txt` or `run.ps1`. Never
publish the prepare cache (`cache/`, which holds the full-resolution raw image) or
the de-identification key file.

### What the host must do

- **Serve files by exact name.** The viewer never lists folders; it asks for
  `manifest.json`, `tiles/0/3/2/5.png` and so on. Directory listing can (and should) be off.
- **Answer HTTP Range requests.** An export keeps each layer's per-gene data in a few
  packs of about 64 MB (`expr/<key>/genepack/*.pack`, see
  [DATA_FORMAT.md](DATA_FORMAT.md#gene-packs-genepackhash-ipack)) and the page reads one gene
  with `Range: bytes=<from>-<to>`, expecting `206 Partial Content`. Blob storage (Azure, S3,
  Google Cloud Storage), Apache, nginx, IIS, GitHub Pages and Netlify do this for static
  files. Do not let the host or a proxy compress the `.pack` files (`Content-Encoding`): a
  compressed file cannot be read by range. A host that ignores Range sends the whole pack;
  the page then manages with packs up to 160 MB (`packWholeMaxMB` in `config.js`), slowly,
  and refuses larger ones with a message. Use a host that answers Range requests
  (blob storage, S3, GitHub Pages, nginx).
- **Allow cross-origin reads, if the data is on a different site than the page.**
  The data host must answer `GET` with `Access-Control-Allow-Origin: <page origin>`
  (or `*` for public data), allow the `Range` request header in a CORS preflight
  (`Access-Control-Allow-Headers: Range`; newer browsers send a simple `bytes=a-b` range
  without asking, older ones and Safari versions ask first) and expose `Content-Range`
  (`Access-Control-Expose-Headers: Content-Range, Content-Length`), which the page uses to
  check that a ranged answer comes from the pack the manifest lists. Same site: nothing to set.
- **Sensible content types.** The viewer reads `.json` and `.bin` as raw bytes whatever
  the server declares, so `application/octet-stream` for `.bin` is fine. Tiles should come
  as `image/png` or `image/jpeg`. If the host adds `X-Content-Type-Options: nosniff`, make
  sure tile types are right.
- **Cope with many files.** An imaging dataset is mostly image tiles: hundreds to tens of
  thousands of files (tiles grow with channels × image area; see
  [PREPARING_DATA.md](PREPARING_DATA.md#how-big-and-how-long)); its genes add only a few packs
  per layer. (Exports made before the packs add one file per gene and layer: 50,000-70,000
  more for a Visium HD dataset.) Check your host's limits on
  file count, file size (packs are about 64 MB; a single very abundant gene can make one
  larger) and total size before you start.
- **Keep names exact.** File names are case-sensitive on most hosts. Do not let a
  deploy tool rename, minify or "optimise" anything in the data folders.

### Host notes

- **GitHub Pages.** Same site for page and data, so no CORS. GitHub limits a published
  site to 1 GB and the repository should stay small, so this suits only small
  datasets. Git is also slow with many thousands of files.
- **Netlify, Cloudflare Pages and similar.** Same site, no CORS. Check the plan's limits on
  files per deploy and on total size.
- **A lab or departmental web server (Apache, nginx, IIS).** Copy the files into the web
  root. If data and page are on different host names, add the CORS header for the data
  folder. nginx example:

  ```nginx
  location /viewer-data/ {
      add_header Access-Control-Allow-Origin "https://viewer.example.org" always;
      add_header Access-Control-Expose-Headers "Content-Range, Content-Length" always;
      add_header Access-Control-Allow-Headers "Range" always;
      types { application/json json geojson; application/octet-stream bin pack; image/png png; image/jpeg jpg; }
      gzip off;                     # a compressed .pack cannot be read by range
  }
  ```

  nginx answers Range requests for static files by itself (`206`); a preflight `OPTIONS`
  request then gets the headers above.

  For IIS, add MIME maps for `.bin` and `.pack` (`application/octet-stream`) and `.geojson`
  (`application/json`), or IIS refuses to serve those files.
- **Amazon S3 static website / a public S3 bucket.** Upload with
  `aws s3 sync export s3://<bucket>/ ` (it sets content types from the extension).
  Bucket CORS configuration, if the page is elsewhere:

  ```json
  [{ "AllowedOrigins": ["https://viewer.example.org"], "AllowedMethods": ["GET", "HEAD"],
     "AllowedHeaders": ["*"], "ExposeHeaders": ["Content-Range", "Content-Length"], "MaxAgeSeconds": 3600 }]
  ```

- **A public Azure blob container.** See option A in
  [deploy/azure/README.md](../deploy/azure/README.md#a-public-container-for-data-that-may-be-public).

## C. Azure, private data with a password

The data stays in a private Azure blob container. A small function checks one shared
password and gives the browser a read-only token that expires after an hour by default
(`TOKEN_HOURS`, at most 24). This is how
the viewer runs in production. Step-by-step instructions, including CORS, upload with
AzCopy and the Static Web App setup:
**[deploy/azure/README.md](../deploy/azure/README.md)**.

The access token is kept in memory only (not in the address, not in browser storage). A link
made with **Copy link to this view** holds neither the token nor the password: whoever opens it is asked for the password first, then sees
the view.

Without a server, keep the data in a private container and share a link with a read pass
(an Azure SAS; [SHARE_YOUR_DATA.md, Step 5](SHARE_YOUR_DATA.md#step-5-make-the-link)).

## D. AWS S3 + CloudFront, or another backend

The viewer does not care which cloud holds the data. It only needs an endpoint that
behaves like this:

```
POST <tokenEndpoint>   body {"password": "..."}   (Content-Type: application/json)
  200  {"base": "https://.../", "sas": "<query string>", "expires": <ms since 1970>}
  401  wrong password          429  too many attempts (the page shows "wait a few minutes")

GET  <tokenEndpoint>   (no body; the browser sends its cookies)
  200  the same JSON, if a session cookie set by an earlier POST is still valid
  401  otherwise (the page then asks for the password)
```

The page then fetches `<base><dataset id>/<file>?<sas>` for every file, and calls the
endpoint again a few minutes before `expires`. So:

- `base` must end in `/` and contain one folder per dataset id;
- `sas` is appended after a `?` to every path under `base`, unchanged, so it must be one
  query string that is valid for **every** object under `base`;
- the page calls the endpoint with `credentials: "same-origin"`, so cookie renewal only
  works when the endpoint is on the same site as the page. A cross-site endpoint still
  works, but people retype the password each time the token expires.

On AWS, plain S3 pre-signed URLs do not fit, because each one is signed for a single
object. A CloudFront **signed URL with a custom policy** does: its policy may name a
wildcard resource such as `https://dxxxx.cloudfront.net/*`, and its query string
(`Policy=...&Signature=...&Key-Pair-Id=...`) then works for any file under it. An outline:

1. A private S3 bucket holding the export folders; a CloudFront distribution in front of
   it with Origin Access Control, a trusted key group, and a response-headers policy that
   adds CORS for your page's origin.
2. A small function (Lambda behind a function URL or API Gateway, ideally routed through
   the same CloudFront distribution as the page under `/api/*`) that checks the password
   and returns `base = "https://dxxxx.cloudfront.net/"`, `sas = "Policy=...&Signature=...&Key-Pair-Id=..."`,
   and `expires`. Port the throttling and session-cookie logic from
   `deploy/azure/api/token/index.js`.
3. `tokenEndpoint` in `config.js` set to that function's address.

We have not built or tested this AWS setup; treat it as a design, not a recipe. The same
reasoning applies to other clouds: it works wherever the provider can issue one
time-limited query string that is valid for every object under a prefix.

## E. Optional: live tile server (server.py)

`server.py` serves the **prepare cache** directly, with no export step. It makes each tile
on request. Use it to look at a freshly prepared dataset on a lab workstation, or to share
it on your local network.

```bash
python web/SpatialViewer/server.py --cache cache --port 8760 --host 127.0.0.1 --password "<a password>"
```

Then open `http://localhost:8760/SpatialViewer/index.html?api=http://localhost:8760`
and type the password into the page's login box. The server lists every dataset folder in
the cache that has a `manifest.json`.

If you have a `config.js`, it must have `tokenEndpoint: ""` for this mode. With a token
endpoint configured, the page tries to fetch every dataset from the cloud container
instead of from the server. (`useLocalServer: true` in `config.js` lets you leave out
`?api=` when the page is served by `server.py` itself.)

- It refuses to start without a password. Give one with `--password`, the
  `SPATIALVIZ_PASSWORD` environment variable, or a `.password` file next to `server.py`.
  `--no-password` is for testing on your own machine only.
- `--host` defaults to `0.0.0.0`, which accepts connections from other computers.
  Use `--host 127.0.0.1` to keep it to your own machine. The password travels as plain
  HTTP Basic authentication, so put the server behind HTTPS (a reverse proxy) before
  exposing it beyond a trusted network.
- It serves images, gene layers and overlays. The whole-section UMAP, Leiden clusters,
  the circle → UMAP tool and region linking are made by `prepare.py export` and only work
  with exported data (options A-D).
- A page served over `https://` cannot read from an `http://` server (the browser blocks
  mixed content), so `?api=` from a hosted page only works if the server also has HTTPS.

---

## Before you upload: checklist

1. **De-identify and check.** Run the gate on the export and on your config:

   ```bash
   python web/SpatialViewer/deid.py check export/my_dataset web/SpatialViewer/config.js
   ```

   It must end with `[deid] clean`. It finds only sample-ID patterns, so also **read the
   labels yourself**: dataset labels in `config.js` and `manifest.json`, channel names,
   overlay and region names, and the category values of RNA-seq color-by columns. See
   [DEIDENTIFICATION.md](DEIDENTIFICATION.md).
2. **No local paths.** `prepare.py export` drops the `source` block (your file paths) from
   the manifest. Check that nothing else slipped in:

   ```bash
   grep -rl -e "C:\\\\" -e "/Users/" -e "/home/" export/my_dataset --include=*.json
   ```

   Windows PowerShell:

   ```powershell
   Get-ChildItem export\my_dataset -Recurse -Include *.json | Select-String -SimpleMatch -Pattern ':\\', '/Users/', '/home/' | Select-Object -First 5
   ```

   No output is the good result.
3. **Export finished cleanly.** The export log says `[export] verified N channel(s) against
   the pyramid` and has no line starting with `***`. If it does, follow the instruction on
   that line before uploading.
4. **Count files and size**, to check against your host's limits and to know how long
   the upload will take.

   macOS / Linux:

   ```bash
   find export/my_dataset -type f | wc -l
   du -sh export/my_dataset
   ```

   Windows PowerShell:

   ```powershell
   Get-ChildItem export\my_dataset -Recurse -File | Measure-Object -Property Length -Sum
   ```

   (`Count` is the number of files, `Sum` the bytes.)
5. **Test the exact folder locally** with option A (`?static=`), and click through every
   layer, gene panel and the UMAP view, before you upload it.
6. **Upload only the export.** Not the cache folder, not the raw instrument output, and
   never the key file `sample_id_key.csv`.
7. **Permission.** Make sure your data-use agreement, consent and IRB approval allow
   sharing in the form you chose (public link vs password).
