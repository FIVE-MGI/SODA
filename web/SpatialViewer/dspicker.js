// SODA dataset picker: a button that opens a searchable list, in place of a <select>.
// A <select> with a hundred datasets has no search, and every arrow-key press on it opened
// a whole dataset. Here: type to filter on the label and the group heading, the arrow keys
// only move the highlight, Enter or a click opens the dataset, Esc closes the list. Works
// with a keyboard, a mouse and a finger.
//
//   const picker = SodaDatasetPicker.create({ host, labelledBy, onPick });
//   picker.setItems(datasets, groups);   // datasets: [{id, label, group, tags}], groups: {key: heading} in order
//   picker.setValue(id);                 // marks the current dataset; never calls onPick
//
// Filters: a dataset may carry `tags`, short labels such as {tissue: "nerve", condition:
// "injured", slide: "S12"} (config.js `datasets`, or a data folder's datasets.json). When any
// dataset has them, a row of drop-downs sits under the search box, one per tag name, with the
// values found in the data and how many datasets each would leave. They combine with each other
// and with the search, and the search also looks in the tag values. No tags anywhere: no row,
// and the list is exactly what it was. Tags are read as untrusted text (strings or numbers only,
// cut short) and only ever shown as text.
//
// Roles: the button has aria-haspopup="listbox"; the search box is the combobox
// (aria-controls the list, aria-activedescendant the highlighted option); the list is a
// listbox of role="group" sections; the current dataset is the aria-selected option. The
// filters are plain <select>s (labelled), so a keyboard, a screen reader and a phone's own
// picker all work with them; Tab moves from the search box to them, Esc closes the list.
// No dependencies. Styles: dspicker.css.
(function () {
  "use strict";
  let made = 0;
  const PAGE = 10;                                   // PageUp / PageDown step
  const TAG_KEYS = 20, KEY_LEN = 40, VAL_LEN = 80;   // per dataset: at most this many tags, this long
  const FILTERS = 6;                                 // drop-downs shown at most (the first tag names met)
  const sortVals = (typeof Intl !== "undefined" && Intl.Collator)
    ? new Intl.Collator(undefined, { numeric: true, sensitivity: "base" }).compare
    : (a, b) => (a < b ? -1 : a > b ? 1 : 0);

  // {name: value} → [[name, value], ...]: strings (or numbers) only, trimmed, cut short.
  function cleanTags(t) {
    const out = [];
    if (!t || typeof t !== "object" || Array.isArray(t)) return out;
    for (const k of Object.keys(t)) {
      if (out.length >= TAG_KEYS) break;
      const v = t[k];
      if (typeof v !== "string" && !(typeof v === "number" && isFinite(v))) continue;
      const key = k.trim().slice(0, KEY_LEN), value = String(v).trim().slice(0, VAL_LEN);
      if (key && value) out.push([key, value]);
    }
    return out;
  }

  function create(opts) {
    const host = opts.host, onPick = opts.onPick || function () {};
    const n = ++made, uid = (s) => `dsp${n}-${s}`;
    let items = [];          // display order: {id, label, heading, lab, head, tg, tk, gi, el}
    let heads = [];          // one per group shown: {label, el, box}
    let shown = [];          // indices into items that pass the filter, in display order
    let facets = [];         // one per filter drop-down: {key, label, values: [{v, label, opt}], sel, el}
    let fil = null, clr = null;   // the filter row (only when there are facets) and its Clear button
    const kept = new Map();  // filter choices, kept when the list is closed or the items change
    let hot = -1, cur = null, built = false, isOpen = false;

    host.classList.add("dsp");
    host.innerHTML =
      `<button type="button" class="dsp-btn" id="${uid("btn")}" aria-haspopup="listbox" aria-expanded="false" aria-controls="${uid("list")}">` +
        `<span class="dsp-val"></span><i aria-hidden="true"></i></button>` +
      `<div class="dsp-pop" hidden>` +
        `<div class="dsp-top">` +
          `<input type="text" class="dsp-search" id="${uid("search")}" role="combobox" aria-expanded="false" aria-autocomplete="list" ` +
            `aria-controls="${uid("list")}" aria-label="Search datasets" placeholder="Search datasets…" ` +
            `autocomplete="off" autocapitalize="off" autocorrect="off" spellcheck="false" />` +
          `<button type="button" class="dsp-x" aria-label="Close the dataset list">×</button>` +
        `</div>` +
        `<div class="dsp-list" id="${uid("list")}" role="listbox" tabindex="-1"></div>` +
        `<div class="dsp-count" role="status" aria-live="polite"></div>` +
      `</div>`;
    const btn = host.querySelector(".dsp-btn"), val = host.querySelector(".dsp-val"), pop = host.querySelector(".dsp-pop");
    const search = host.querySelector(".dsp-search"), list = host.querySelector(".dsp-list"), count = host.querySelector(".dsp-count");
    if (opts.labelledBy) {
      btn.setAttribute("aria-labelledby", `${opts.labelledBy} ${btn.id}`);
      list.setAttribute("aria-labelledby", opts.labelledBy);
    } else list.setAttribute("aria-label", "Datasets");

    // Same grouping as the <select> had: the headings of `groups` in their order, then
    // "Other" for a dataset whose group is not among them; no groups at all = a flat list.
    function setItems(datasets, groups) {
      groups = groups || {};
      const order = Object.keys(groups).concat([undefined]);
      const bucket = new Map(order.map((g) => [g, []]));
      for (const d of datasets) bucket.get(bucket.has(d.group) ? d.group : undefined).push(d);
      items = []; heads = [];
      for (const g of order) {
        const members = bucket.get(g);
        if (!members.length) continue;
        const flat = g === undefined && order.length === 1;
        const heading = flat ? "" : g === undefined ? "Other" : String(groups[g]);
        const gi = heads.length;
        heads.push({ label: heading, el: null, box: null });
        for (const d of members) {
          const label = String(d.label || d.id), tags = cleanTags(d.tags);
          // tk: tag name → value, both lower case (names and values match whatever their case)
          const tk = new Map(tags.map(([k, v]) => [k.toLowerCase(), v.toLowerCase()]));
          items.push({ id: d.id, label, heading, lab: label.toLowerCase(), head: heading.toLowerCase(),
            tg: tags.map((x) => x[1]).join("\n").toLowerCase(), tk, tags, gi, el: null });
        }
      }
      makeFacets();
      built = false;
      search.placeholder = items.length > 1 ? `Search ${items.length.toLocaleString()} datasets…` : "Search datasets…";
      syncButton();
      if (isOpen) { build(); filter(); }
    }

    // One drop-down per tag name, in the order names first occur, values sorted (S2 before S10).
    // A name every dataset has with the same value would filter nothing: no drop-down for it.
    function makeFacets() {
      const byKey = new Map();
      for (const it of items) for (const [k, v] of it.tags) {
        const key = k.toLowerCase(), val = v.toLowerCase();
        let f = byKey.get(key);
        if (!f) byKey.set(key, (f = { key, label: k, vals: new Map(), have: 0 }));
        if (!f.vals.has(val)) f.vals.set(val, v);
      }
      for (const it of items) for (const f of byKey.values()) if (it.tk.has(f.key)) f.have++;
      facets = [...byKey.values()]
        .filter((f) => f.vals.size > 1 || f.have < items.length)
        .slice(0, FILTERS)
        .map((f) => ({ key: f.key, label: f.label, el: null, sel: "",
          values: [...f.vals].map(([v, label]) => ({ v, label, opt: null })).sort((a, b) => sortVals(a.label, b.label)) }));
      if (fil) { fil.remove(); fil = clr = null; }
      if (!facets.length) return;                    // no tags: the picker is exactly as it was
      fil = document.createElement("div");
      fil.className = "dsp-fil"; fil.setAttribute("role", "group"); fil.setAttribute("aria-label", "Filter datasets");
      facets.forEach((f, fi) => {
        const lab = document.createElement("label");
        lab.className = "dsp-f";
        const name = document.createElement("span");
        name.className = "dsp-fk"; name.textContent = f.label;
        const s = document.createElement("select");
        s.className = "dsp-fs"; s.id = uid(`f${fi}`);
        const any = document.createElement("option");
        any.value = ""; any.textContent = "Any";
        s.appendChild(any);
        f.values.forEach((x, vi) => {
          const o = document.createElement("option");
          o.value = String(vi); o.textContent = x.label;
          x.opt = o; s.appendChild(o);
        });
        const was = kept.get(f.key), at = was ? f.values.findIndex((x) => x.v === was) : -1;
        if (at >= 0) { f.sel = was; s.value = String(at); }
        s.addEventListener("change", () => {
          f.sel = s.value === "" ? "" : f.values[Number(s.value)].v;
          if (f.sel) kept.set(f.key, f.sel); else kept.delete(f.key);
          filter();
        });
        f.el = s;
        lab.appendChild(name); lab.appendChild(s); fil.appendChild(lab);
      });
      clr = document.createElement("button");
      clr.type = "button"; clr.className = "dsp-clr"; clr.textContent = "Clear filters";
      clr.addEventListener("click", () => {
        for (const f of facets) { f.sel = ""; f.el.value = ""; kept.delete(f.key); }
        filter();
        search.focus();
      });
      fil.appendChild(clr);
      fil.addEventListener("keydown", (e) => {
        if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); close(true); }
        // Enter on a drop-down opens the first dataset left, as Enter in the search box does
        else if (e.key === "Enter" && e.target.tagName === "SELECT" && hot >= 0) { e.preventDefault(); choose(hot); }
      });
      pop.insertBefore(fil, list);
    }

    // The option elements are made the first time the list opens, not at page load.
    function build() {
      const frag = document.createDocumentFragment();
      heads.forEach((h, gi) => {
        const box = document.createElement("div");
        box.className = "dsp-g"; box.setAttribute("role", "group");
        if (h.label) {
          const t = document.createElement("div");
          t.className = "dsp-h"; t.id = uid(`g${gi}`); t.textContent = h.label; t.setAttribute("role", "presentation");
          box.setAttribute("aria-labelledby", t.id);
          box.appendChild(t); h.el = t;
        }
        h.box = box; frag.appendChild(box);
      });
      items.forEach((it, i) => {
        const o = document.createElement("div");
        o.className = "dsp-o"; o.id = uid(`o${i}`); o.setAttribute("role", "option");
        o.dataset.i = i; o.textContent = it.label;
        it.el = o; heads[it.gi].box.appendChild(o);
      });
      list.textContent = ""; list.appendChild(frag);
      built = true;
      markCurrent();
    }
    function markCurrent() {
      if (!built) return;
      for (const it of items) {
        const on = it.id === cur;
        it.el.classList.toggle("cur", on);
        it.el.setAttribute("aria-selected", on ? "true" : "false");
      }
    }
    function syncButton() {
      const it = items.find((x) => x.id === cur);
      val.textContent = it ? it.label : cur == null ? "" : String(cur);
      btn.title = val.textContent;
    }

    // Scroll the list (and only the list: scrollIntoView would also move the sidebar and, on a
    // phone, the page) so the option shows; a group's heading sticks to the top, so leave room for it.
    function reveal(i, center) {
      const el = items[i].el, h = heads[items[i].gi].el;
      const top = el.offsetTop, bottom = top + el.offsetHeight, room = h ? h.offsetHeight : 0;
      if (center) list.scrollTop = top - (list.clientHeight - el.offsetHeight) / 2;
      else if (top - room < list.scrollTop) list.scrollTop = top - room;
      else if (bottom > list.scrollTop + list.clientHeight) list.scrollTop = bottom - list.clientHeight;
    }
    function setHot(i, scroll) {
      if (hot >= 0 && items[hot] && items[hot].el) items[hot].el.classList.remove("hot");
      hot = i;
      if (i < 0) { search.removeAttribute("aria-activedescendant"); return; }
      const el = items[i].el;
      el.classList.add("hot");
      search.setAttribute("aria-activedescendant", el.id);
      if (scroll) reveal(i, false);
    }
    // What is typed is looked for as it stands in the label or in the group heading ("cohort 2"
    // lists that group). When nothing holds it whole, every word must be somewhere in the two
    // ("normal coh"). The tag values count as part of the label here ("injured" finds the
    // datasets tagged so). With nothing typed or chosen the highlight starts on the current
    // dataset; with a filter, on the first match.
    // Each value in a drop-down shows how many datasets it would leave, given the search and
    // the other drop-downs.
    function filter() {
      const words = search.value.toLowerCase().split(/\s+/).filter(Boolean), phrase = words.join(" ");
      const has = (it, w) => it.lab.indexOf(w) >= 0 || it.head.indexOf(w) >= 0 || it.tg.indexOf(w) >= 0;
      const whole = (it) => has(it, phrase);
      const each = (it) => words.every((w) => has(it, w));
      const match = !words.length ? () => true : words.length === 1 || items.some(whole) ? whole : each;
      const on = facets.filter((f) => f.sel);
      const tally = facets.map(() => new Map());
      const perGroup = heads.map(() => 0);
      shown = [];
      for (let i = 0; i < items.length; i++) {
        const it = items[i], text = match(it);
        let miss = 0, missed = null;                 // how many chosen filters it fails, and which
        for (const f of on) if (it.tk.get(f.key) !== f.sel) { miss++; missed = f; }
        const ok = text && !miss;
        it.el.hidden = !ok;
        if (ok) { shown.push(i); perGroup[it.gi]++; }
        if (text && miss <= 1) facets.forEach((f, fi) => {
          if (miss && missed !== f) return;
          const v = it.tk.get(f.key);
          if (v !== undefined) tally[fi].set(v, (tally[fi].get(v) || 0) + 1);
        });
      }
      heads.forEach((h, gi) => { h.box.hidden = !perGroup[gi]; });
      facets.forEach((f, fi) => {
        for (const x of f.values) x.opt.textContent = `${x.label} (${(tally[fi].get(x.v) || 0).toLocaleString()})`;
        f.el.classList.toggle("on", !!f.sel);
      });
      if (clr) clr.classList.toggle("on", on.length > 0);
      const narrowed = words.length > 0 || on.length > 0, total = items.length;
      count.textContent = !shown.length ? (on.length ? "No dataset matches these filters." : "No dataset matches.")
        : narrowed ? `${shown.length.toLocaleString()} of ${total.toLocaleString()} datasets`
        : `${total.toLocaleString()} dataset${total === 1 ? "" : "s"}`;
      count.classList.toggle("none", !shown.length);
      const at = narrowed ? -1 : items.findIndex((x) => x.id === cur);
      setHot(shown.length ? (at >= 0 ? at : shown[0]) : -1, false);
      if (narrowed || hot < 0) list.scrollTop = 0;
      else reveal(hot, true);
    }
    function move(step) {
      if (!shown.length) return;
      const at = shown.indexOf(hot);
      const to = at < 0 ? (step > 0 ? 0 : shown.length - 1) : Math.max(0, Math.min(shown.length - 1, at + step));
      setHot(shown[to], true);
    }

    const coarse = () => !!(window.matchMedia && matchMedia("(pointer: coarse)").matches);
    // seed: text to start the filter with (a letter typed on the closed button).
    // focusSearch false: opened by a finger, so the on-screen keyboard stays down until asked for.
    function open(seed, focusSearch) {
      if (!items.length) return;
      if (!built) build();
      isOpen = true;
      pop.hidden = false; host.classList.add("open");
      btn.setAttribute("aria-expanded", "true"); search.setAttribute("aria-expanded", "true");
      search.value = seed || "";
      filter();
      if (focusSearch !== false) search.focus(); else list.focus({ preventScroll: true });
      document.addEventListener("pointerdown", outside, true);
    }
    function close(focusButton) {
      if (!isOpen) return;
      isOpen = false;
      pop.hidden = true; host.classList.remove("open");
      btn.setAttribute("aria-expanded", "false"); search.setAttribute("aria-expanded", "false");
      search.value = ""; setHot(-1, false);
      document.removeEventListener("pointerdown", outside, true);
      if (focusButton) btn.focus({ preventScroll: true });
    }
    function outside(e) { if (!host.contains(e.target)) close(false); }
    // Picking the dataset that is already open only closes the list, as a <select> would.
    function choose(i) {
      const it = items[i]; if (!it) return;
      close(true);
      if (it.id === cur) return;
      setValue(it.id);
      onPick(it.id);
    }
    function setValue(id) { cur = id; syncButton(); markCurrent(); }

    let byFinger = false;
    btn.addEventListener("pointerdown", (e) => { byFinger = e.pointerType === "touch" || (e.pointerType !== "mouse" && coarse()); });
    btn.addEventListener("click", () => {
      const finger = byFinger; byFinger = false;
      if (isOpen) close(true); else open("", !finger);
    });
    btn.addEventListener("keydown", (e) => {
      if (isOpen || e.ctrlKey || e.metaKey || e.altKey) return;
      if (e.key === "ArrowDown" || e.key === "ArrowUp") { e.preventDefault(); open("", true); }
      else if (e.key.length === 1 && e.key !== " ") { e.preventDefault(); open(e.key, true); }   // type-ahead, as on a <select>
    });
    search.addEventListener("input", filter);
    function keys(e) {
      if (e.key === "ArrowDown") { e.preventDefault(); move(1); }
      else if (e.key === "ArrowUp") { e.preventDefault(); move(-1); }
      else if (e.key === "PageDown") { e.preventDefault(); move(PAGE); }
      else if (e.key === "PageUp") { e.preventDefault(); move(-PAGE); }
      else if (e.key === "Enter") { e.preventDefault(); if (hot >= 0) choose(hot); }
      else if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); close(true); }   // only the list: not the page's own Esc
      else if (e.target === list && e.key.length === 1 && !e.ctrlKey && !e.metaKey && !e.altKey) search.focus();   // typing: into the search box
    }
    search.addEventListener("keydown", keys);
    list.addEventListener("keydown", keys);
    // Tab (or Shift+Tab) out of the list closes it; Tab from the search box to the filters does
    // not. A click on the list's own background moves the focus nowhere (no relatedTarget): that
    // stays open, and a click outside is caught by `outside`.
    pop.addEventListener("focusout", (e) => {
      if (isOpen && e.relatedTarget && !pop.contains(e.relatedTarget)) close(false);
    });
    // the mouse must not take the focus out of the search box (touch scrolling is left alone)
    list.addEventListener("mousedown", (e) => { if (document.activeElement === search) e.preventDefault(); });
    list.addEventListener("click", (e) => {
      const o = e.target.closest && e.target.closest(".dsp-o");
      if (o && list.contains(o)) choose(Number(o.dataset.i));
    });
    list.addEventListener("pointermove", (e) => {
      if (e.pointerType !== "mouse") return;
      const o = e.target.closest && e.target.closest(".dsp-o");
      if (o && Number(o.dataset.i) !== hot) setHot(Number(o.dataset.i), false);
    });
    host.querySelector(".dsp-x").addEventListener("click", () => close(true));

    return {
      setItems, setValue, open: (seed) => open(seed || "", true), close: () => close(false),
      value: () => cur, isOpen: () => isOpen, count: () => items.length,
      // what a click on the option does (tests, console)
      pick: (id) => { const i = items.findIndex((x) => x.id === id); if (i >= 0) { if (!built) build(); choose(i); } return i >= 0; },
    };
  }

  window.SodaDatasetPicker = { create };
})();
