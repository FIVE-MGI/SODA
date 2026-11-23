// site.js: what the page does when no data is named anywhere. Kept in git and the same on every
// copy of SODA; your own settings (your data folder, title, logos) go in config.js, not here.
//
//   data         an export folder the page opens when neither the link (?data=, ?static=) nor
//                config.js names any data, datasets or password (tokenEndpoint). Empty on purpose:
//                a copy of SODA opens its own folder, named in config.js. This file is public: only
//                a public-read folder could ever go here, never a link with a read pass (?sv=…&sig=…).
//   helpUrl      "how to connect your data", shown when nothing is connected
//   examplesUrl  a SODA page with public example datasets, shown when nothing is connected
//
// showNoData(): the calm "No data connected yet" panel that viewer.js asks for when it has nothing
// to open (instead of an error). It returns true when it drew the panel.
window.SodaSite = {
  data: "",
  helpUrl: "https://github.com/FIVE-MGI/SODA/blob/main/docs/SHARE_YOUR_DATA.md",
  examplesUrl: "https://sidolans01.mgifive.org/SODA-examples/",

  showNoData() {
    if (typeof document === "undefined" || !document.body || document.getElementById("sodaNoData")) return false;
    const css = document.createElement("style");
    css.textContent = `
      .soda-nodata { position: fixed; inset: 0; z-index: 1000; display: flex; align-items: center; justify-content: center;
        padding: 16px; background: var(--bg, #0e1117); color: var(--text, #e6e9ef); overflow-y: auto; }
      .soda-nodata .card { max-width: 30em; width: 100%; padding: 32px 28px; text-align: center; line-height: 1.55;
        background: var(--panel, #161b24); border: 1px solid var(--line, #2a3242); border-radius: 12px; }
      .soda-nodata .mark { font-size: 40px; font-weight: 700; letter-spacing: .14em; margin: 0 0 18px; }
      .soda-nodata h2 { font-size: 19px; font-weight: 600; margin: 0 0 8px; }
      .soda-nodata p { font-size: 14px; color: var(--muted, #8b93a3); margin: 0 0 22px; }
      .soda-nodata a { display: block; margin: 10px auto 0; padding: 9px 14px; max-width: 22em; border-radius: 8px;
        border: 1px solid var(--line, #2a3242); color: var(--accent, #4da3ff); text-decoration: none; font-size: 14px; }
      .soda-nodata a:hover, .soda-nodata a:focus-visible { border-color: var(--accent, #4da3ff); }`;
    document.head.appendChild(css);
    const el = (tag, props, text) => { const e = document.createElement(tag); Object.assign(e, props || {}); if (text) e.textContent = text; return e; };
    const link = (href, text) => el("a", { href, target: "_blank", rel: "noopener" }, text);
    const box = el("div", { id: "sodaNoData", className: "soda-nodata" });
    const card = el("section", { className: "card" });
    card.setAttribute("aria-labelledby", "sodaNoDataTitle");
    card.append(
      el("div", { className: "mark", role: "img", ariaLabel: "SODA" }, "SODA"),
      el("h2", { id: "sodaNoDataTitle" }, "No data connected yet"),
      el("p", {}, "This SODA site is ready, but no data folder has been connected to it yet."),
      link(this.helpUrl, "How to connect your data (for the person who runs this site)"),
      link(this.examplesUrl, "See SODA with example data"));
    box.append(card);
    document.body.appendChild(box);
    return true;
  },
};
