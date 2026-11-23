// config.js: the settings of YOUR SODA site. This is the one file you edit.
//
// How to edit it on GitHub: open this file in your copy of SODA, click the pencil icon
// ("Edit this file"), make the change, then click "Commit changes". Your site updates by
// itself a minute or two later.
//
// Everything here is public (anyone can read this file on your site), so never put a password,
// a key or a private link in it.
window.SpatialViewerConfig = {

  // ---- your data -------------------------------------------------------------------
  // The web address of your storage folder: the folder you uploaded your SODA export to
  // (the one holding datasets.json). Paste it between the quotes, for example
  //   data: "https://yourstorage.blob.core.windows.net/soda/",
  // While this is empty, your site shows "No data connected yet".
  // Step by step: docs/SHARE_YOUR_DATA.md in the SODA repository.
  data: "",

  // ---- how your site looks -----------------------------------------------------------
  brand: {
    // The name at the top of the page and in the browser tab.
    title: "SODA",

    // Your logos, shown above the title. Upload the image file into the same folder as this file
    // (web/SpatialViewer, "Add file" -> "Upload files" on GitHub), then remove the two slashes (//)
    // in front of the line below and change the file name and the alt text (what a screen reader
    // says instead of the picture). "logo-wide" suits a wide logo, "logo-round" a round one.
    logos: [
      // { src: "my-lab-logo.png", alt: "My lab", className: "logo-wide" },
    ],

    // One line of small print under the controls, e.g. your lab's name. Optional.
    byline: "",
  },

  // More settings (dataset names and groups, colours, limits) are explained in
  // config.example.js, next to this file. Copy only the lines you need into this file.
};
