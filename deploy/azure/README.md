# Hosting on Azure with a password

This folder sets up the viewer the way it is run in production:

- the page is an **Azure Static Web App** (free plan is enough),
- the data is in a **private** Azure blob container,
- a small function, `/api/token`, checks one shared password and hands the browser a
  read-only access token (a SAS) that expires after an hour by default (`TOKEN_HOURS`, at
  most 24). The storage key stays on
  the server.

```
browser ── page files ───────────────> Static Web App   (index.html, viewer.js, config.js)
browser ── POST /api/token {password} > token function  (this folder: api/token/index.js)
        <─ {base, sas, expires} ──────
browser ── GET <base><id>/...?<sas> ──> private blob container (your exported datasets)
```

Before you start you need:

- An Azure subscription where you may create a storage account and a Static Web App.
- A GitHub account. The Static Web App deploys from your fork of this repository.
- Exported datasets (`prepare.py export`, see [docs/PREPARING_DATA.md](../../docs/PREPARING_DATA.md))
  that you have checked with `deid.py check` (see [docs/DEIDENTIFICATION.md](../../docs/DEIDENTIFICATION.md)).
- One upload tool: [AzCopy](https://learn.microsoft.com/azure/storage/common/storage-use-azcopy-v10)
  (command line, best for many files) or
  [Azure Storage Explorer](https://azure.microsoft.com/products/storage/storage-explorer/) (graphical).
- Optional: the [Azure CLI](https://learn.microsoft.com/cli/azure/install-azure-cli) (`az`).
  Every step below can also be done by clicking in the Azure portal.

Placeholders used below. Replace them everywhere, including the angle brackets:

| Placeholder | Meaning | Example |
|---|---|---|
| `<rg>` | resource group | `soda-rg` |
| `<account>` | storage account name: 3-24 lowercase letters/digits, globally unique | `mylabsodadata` |
| `<container>` | container name: lowercase letters, digits, hyphens | `soda-data` |
| `<app-name>` | the Static Web App resource name | `soda-mylab` |
| `<site>` | your Static Web App host name, shown in the portal after step 4; if you add a custom domain, please include `soda` (e.g. `soda.yourlab.org`) | `nice-pond-0123.azurestaticapps.net` |
| `<UPLOAD-SAS>` | a write SAS for uploading, made in step 3 | `sv=...&sig=...` |

---

## 1. Create the storage account and a private container

**Portal:** Create a resource → Storage account. Pick your subscription, a resource group,
a name, a region near your users, performance *Standard*, redundancy *LRS*. On the
*Advanced* tab, leave **Allow enabling anonymous access on individual containers** off.
Create it. Then open it → Data storage → Containers → **+ Container**, give it a name,
and keep **Anonymous access level: Private (no anonymous access)**.

**CLI (same thing):**

```bash
az group create --name <rg> --location eastus
az storage account create --name <account> --resource-group <rg> --location eastus \
    --sku Standard_LRS --kind StorageV2 --allow-blob-public-access false --min-tls-version TLS1_2
az storage container create --account-name <account> --name <container> --auth-mode login
```

(`--auth-mode login` needs the *Storage Blob Data Contributor* role on the account. If that
fails, use `--auth-mode key` instead.) On Windows PowerShell, replace the trailing `\`
line continuations with a backtick `` ` ``, or put each command on one line.

## 2. Allow the site to read from the container (CORS)

The page and the data are on different web addresses, so the blob service must say
that your site may read from it. Without this, the browser blocks every request and
the developer console shows a "CORS" error.

**Portal:** storage account → Settings → **Resource sharing (CORS)** → **Blob service** tab.
Add one row:

| Field | Value |
|---|---|
| Allowed origins | `https://<site>` (add your custom domain too if you use one; add `http://localhost:8766` only while testing locally) |
| Allowed methods | `GET`, `HEAD`, `OPTIONS` |
| Allowed headers | `*` |
| Exposed headers | `*` |
| Max age | `3600` |

The viewer sends no custom request headers to storage and only reads the response
body and status, so `*` for headers is more than it needs; it is the simplest value
that works. Do not use `*` for *origins* on a private container.

You only know `<site>` after step 4. Come back and fill it in then.

**CLI:**

```bash
az storage cors add --account-name <account> --services b \
    --methods GET HEAD OPTIONS --origins "https://<site>" \
    --allowed-headers "*" --exposed-headers "*" --max-age 3600
```

## 3. Upload your exported datasets

The export folder looks like `export/<id>/manifest.json`, `export/<id>/tiles/...`.
Each `<id>` folder must end up at the **top of the container**:
`https://<account>.blob.core.windows.net/<container>/<id>/manifest.json`.

### Make a temporary upload SAS

Portal: container → Settings → **Shared access tokens**. Signing method *Account key*,
permissions **Read, Add, Create, Write, List**, expiry tomorrow, allowed protocols
*HTTPS only*. Click **Generate SAS token and URL** and copy the **token** (the part after `?`).
This is a write credential: do not paste it into chat, email, config.js or a commit.

CLI equivalent:

```bash
az storage container generate-sas --account-name <account> --name <container> \
    --permissions racwl --expiry 2030-01-01T00:00Z --https-only --auth-mode key -o tsv
```

(Use a date about one day ahead, not 2030.)

### Upload with AzCopy

An image dataset is tens of thousands of small files. AzCopy handles that well and can
resume an interrupted upload. Run from the repository folder:

```bash
azcopy copy "demo/export/demo_visiumhd" "https://<account>.blob.core.windows.net/<container>?<UPLOAD-SAS>" --recursive
```

That creates `<container>/demo_visiumhd/...`. To upload every dataset in the export
folder at once:

```bash
azcopy copy "demo/export/*" "https://<account>.blob.core.windows.net/<container>?<UPLOAD-SAS>" --recursive
```

Keep the double quotes on Windows too: the SAS contains `&`, which the shell would
otherwise treat as a command separator.

If the upload stops (network drop, laptop sleeps), resume it instead of starting over:

```bash
azcopy jobs list
azcopy jobs resume <job-id> --destination-sas "<UPLOAD-SAS>"
```

Or simply run the same `azcopy copy` again with `--overwrite=false` added: files
already in the container are skipped. After you re-export a dataset, use
`azcopy sync "demo/export/demo_visiumhd" "https://<account>.blob.core.windows.net/<container>/demo_visiumhd?<UPLOAD-SAS>" --recursive`
so changed files are replaced.

AzCopy sets each blob's content type from its file extension (`.json` → `application/json`,
`.png` → `image/png`, `.bin` → `application/octet-stream`), which is what the viewer expects.

### Or upload with Azure Storage Explorer

Sign in, open the storage account → Blob Containers → `<container>`, click **Upload →
Upload Folder** and pick `export/<id>`. Storage Explorer uses AzCopy underneath. For a
first upload of a very large dataset, the command line is easier to resume.

## 4. Create the Static Web App from your fork

1. Fork this repository on GitHub.
2. Azure portal → Create a resource → **Static Web App**. Plan type *Free*.
   Deployment source *GitHub*, then pick your fork and the branch.
3. Build details: Build preset **Custom**, then
   - **App location:** `web/SpatialViewer`
   - **Api location:** `deploy/azure/api`
   - **Output location:** leave empty
4. Create. Azure commits a workflow file to your fork
   (`.github/workflows/azure-static-web-apps-<random>.yml`) and runs it.
5. Edit that workflow file in your fork. Two changes:

   a. Add a step **before** the "Build And Deploy" step, so the site's config file sits
   next to the page (Static Web Apps only reads it from inside the app folder):

   ```yaml
         - name: Add the Static Web Apps config
           run: cp deploy/azure/staticwebapp.config.json web/SpatialViewer/
   ```

   b. In the "Build And Deploy" step's `with:` block, make the locations read:

   ```yaml
             app_location: "web/SpatialViewer"
             api_location: "deploy/azure/api"
             output_location: ""
             skip_app_build: true
   ```

   `skip_app_build: true` uploads the page as it is. The page has no build step, and the
   folder contains a `requirements.txt` that the build system could otherwise take for a
   Python project. The API is still built (its `npm install` runs as usual).

6. Commit. The workflow runs again; its log is under the **Actions** tab of your fork.
   When it is green, the portal's Overview page shows the URL, `https://<site>`.
   Go back to step 2 and put that origin into CORS.

`staticwebapp.config.json` selects Node 20 for the function and returns 404 for the
Python tools, scripts and requirement files that sit next to the page.

## 5. Set the app settings

Portal: Static Web App → Settings → **Environment variables** → *Production*. Add:

| Name | Value | Required |
|---|---|---|
| `VIEWER_PASSWORD` | the shared password people will type | yes |
| `STORAGE_ACCOUNT` | `<account>` | yes |
| `STORAGE_CONTAINER` | `<container>` | yes |
| `STORAGE_KEY` | storage account → Security + networking → **Access keys** → key1 → Key | one of these two |
| `BLOB_SAS` | a read-only (`sp=r`) container SAS, if you may not read account keys | one of these two |
| `TOKEN_HOURS` | lifetime of each token in hours; default 1, at most 24 | no |
| `SESSION_HOURS` | how long a browser may renew without retyping the password; default 1, at most 24 | no |

Click **Apply**. There are no built-in defaults for the account or the container: if a
required setting is missing, `/api/token` answers `500` and names what is missing.

With `STORAGE_KEY`, every visitor gets their own token that expires after
`TOKEN_HOURS`. With `BLOB_SAS` instead, everybody receives that same SAS until its own
expiry date, so keep it read-only and short.

CLI equivalent (the values end up in your shell history, so the portal is safer for the
password and the key):

```bash
az staticwebapp appsettings set --name <app-name> --resource-group <rg> \
    --setting-names VIEWER_PASSWORD="..." STORAGE_ACCOUNT=<account> STORAGE_CONTAINER=<container> STORAGE_KEY="..."
```

## 6. Write config.js

`config.js` names your datasets and your institution. It is listed in `.gitignore` so it
is never committed by accident to the public repository, but the Static Web App deploys
from your fork, so in your fork it must be committed.

macOS / Linux:

```bash
cp web/SpatialViewer/config.example.js web/SpatialViewer/config.js
```

Windows (PowerShell):

```powershell
Copy-Item web\SpatialViewer\config.example.js web\SpatialViewer\config.js
```

Edit it. The parts that matter for this setup:

```js
window.SpatialViewerConfig = {
  datasets: {
    demo_visiumhd: { group: "transcriptomics", label: "Demo section" },   // key = folder name in the container
  },
  defaultDataset: "demo_visiumhd",
  tokenEndpoint: "/api/token",
  // ...the rest as in config.example.js
};
```

Then commit it on purpose and push:

```bash
git add -f web/SpatialViewer/config.js
git commit -m "Site configuration"
git push
```

`config.js` is served to every visitor of the site, before the password. Keep anything
you would not show publicly out of it, including dataset ids and labels that identify
a person. Logos go in `web/SpatialViewer/logos/`, which is also git-ignored: add them
with `git add -f` the same way.

## 7. Test

**The token endpoint.** macOS / Linux:

```bash
curl -i -X POST "https://<site>/api/token" -H "Content-Type: application/json" -d '{"password":"<your password>"}'
```

Windows (PowerShell):

```powershell
Invoke-RestMethod -Method Post -Uri "https://<site>/api/token" -ContentType "application/json" -Body '{"password":"<your password>"}'
```

| Answer | Meaning |
|---|---|
| `200` with `base`, `sas`, `expires` | working |
| `401 wrong password` | the password differs from `VIEWER_PASSWORD` |
| `429` | 30 wrong passwords from your address in 5 minutes; wait 5 minutes |
| `500 server not configured` with a `missing` list | fix those app settings (step 5) |
| `500 could not create an access token` | rare; the function's log (Static Web App → Functions / Application Insights) has the reason |
| `404` | the API was not deployed: check `api_location` in the workflow |

A `200` only proves the password and settings were read. The function cannot tell
whether `STORAGE_KEY` really belongs to `STORAGE_ACCOUNT`; the next test does.

**One data file.** Paste `<base><id>/manifest.json?<sas>` (from that answer) into a
browser. You should see the JSON. `AuthenticationFailed` means the key or SAS does not
match the account; `BlobNotFound` means the folder is not at the top of the container
or the name differs.

**The page.** Open `https://<site>/index.html?ds=<id>`. It asks for the password, then
shows the image. If it stays blank, open the browser's developer console (F12):

- a "CORS" or "Access-Control-Allow-Origin" message → step 2 (origin must match exactly,
  `https://`, no trailing slash);
- `403` on data files → the token expired or the key changed; reload the page;
- `404` on `manifest.json` → the dataset id in `config.js` does not match the folder name.

## Day-to-day

- **Change the password:** edit `VIEWER_PASSWORD` and apply. Every open session stops
  renewing at once (the session cookie is bound to the password); people must type the
  new one. A token already handed out still works until it expires (at most `TOKEN_HOURS`).
- **Cut off all access immediately:** rotate the storage key in the portal (Access keys →
  Rotate key) and put the new key into `STORAGE_KEY`. Every token made with the old key
  stops working.
- **Add a dataset:** export it, run `deid.py check`, upload the `<id>` folder, add one line
  to `config.js`, commit, push.
- **Project description:** the page's "Read more about this project" box loads
  `about.html` from the top of the container (behind the same password). Upload one if
  you want it; without it the box says nothing has been published yet.

## What the function does (for reviewers)

`api/token/index.js`, about 200 lines, one dependency (`@azure/storage-blob`):

- `POST` compares the password with `VIEWER_PASSWORD` in constant time (`crypto.timingSafeEqual`).
- Wrong passwords are slowed by 400 ms each, and an address with 30 failures in 5
  minutes gets `429` for 5 minutes. The counter is in memory, per instance.
- On success it mints a **container-scoped, read-only** SAS (`sp=r`: no list, write or
  delete), HTTPS only, valid from 5 minutes ago (clock skew) for `TOKEN_HOURS`.
- It sets an `HttpOnly; Secure; SameSite=Strict; Path=/api` cookie holding only an
  expiry time and an HMAC. The HMAC key is derived from the storage credential *and* the
  password, so changing either one invalidates every session.
- `GET` with a valid cookie returns a fresh SAS without asking again. That is how a tab
  left open for a long time keeps working.
- Responses are `Cache-Control: no-store`. Errors name missing setting names, never values.

Any other backend can replace it: the viewer only needs
`POST {password} → {base, sas, expires}` and, optionally, cookie-based `GET` renewal
(see [docs/HOSTING.md](../../docs/HOSTING.md)).

---

## Without a server: two simpler options

Both skip the function. Neither can revoke access by changing a password.

### A. Public container (for data that may be public)

1. When creating the storage account, allow anonymous access on containers. Set the
   container's access level to **Blob (anonymous read access for blobs only)**. That
   allows reading a file whose name you know, but not listing the container.
2. Set CORS as in step 2. For public data `*` as the allowed origin is acceptable.
3. Upload as in step 3.
4. In `config.js`, turn the password off and give each dataset its full address:

   ```js
   tokenEndpoint: "",
   datasets: {
     demo_visiumhd: { label: "Demo section",
                      base: "https://<account>.blob.core.windows.net/<container>/demo_visiumhd/" },
   },
   ```

The page itself can then live anywhere: a Static Web App without the API, GitHub Pages,
or your lab's web server (see [docs/HOSTING.md](../../docs/HOSTING.md)).

### B. Encrypted link files (a password, no server)

The data sits in a public-read container under a long random folder name, and the page
holds that address only in encrypted form. The browser decrypts it with the password
(PBKDF2-SHA256, 200,000 rounds, then AES-GCM).

```bash
python web/SpatialViewer/links.py secret
```

prints a random folder name, for example `Qx7...`. Upload each dataset to
`<container>/<random name>/<id>/` in a container with access level *Blob* (as in A).
Then, for each dataset:

```bash
python web/SpatialViewer/links.py encrypt --password "<password>" --url "https://<account>.blob.core.windows.net/<container>/<random name>/demo_visiumhd/" --out web/SpatialViewer/demo_visiumhd.link.txt
```

and in `config.js`:

```js
tokenEndpoint: "",
datasets: {
  demo_visiumhd: { label: "Demo section", linkPath: "demo_visiumhd.link.txt" },
},
```

Commit the `.link.txt` files with the page. Check a link with
`python web/SpatialViewer/links.py decrypt --password "<password>" --file web/SpatialViewer/demo_visiumhd.link.txt`.

Know the limits before you choose this:

- Anyone who has decrypted the address once can keep using it, and share it, for as long
  as the data stays at that address. Changing the password does not lock them out; you
  must move the data to a new random folder and encrypt new link files.
- The link files are public, so a weak password can be guessed offline. Use a long one.
- Decryption needs `https://` (or `http://localhost`).
- `links.py` takes the password on the command line, so it lands in your shell history.
  Clear it afterwards, or run the command in a session that does not keep history.
