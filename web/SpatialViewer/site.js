// The data folder this page opens when nothing else names one. Kept in git, unlike config.js, so a
// site published straight from the repository (GitHub Pages, .github/workflows/pages.yml) has it.
//
//   data: the link of an export folder (the one holding datasets.json, from catalog.py). Here the
//         seven public example datasets (README, Quick start). Empty: the page opens nothing and
//         says how to open a data link. This file is public: only a public-read folder, never a
//         link with a read pass (?sv=…&sig=…).
//
// It is used only when the address has no ?data=, ?static= or ?api=, and config.js names no data,
// datasets, server or password (tokenEndpoint). Any of those wins, so a copy of SODA with its own
// config.js never opens this folder. Only `data` is read from here; every other setting is in
// config.js (see config.example.js).
window.SodaSite = {
  data: "https://fivestore.blob.core.windows.net/soda/examples/",
};
