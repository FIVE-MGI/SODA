// POST /api/token  {"password": "..."}  -> 200 {base, sas, expires}  | 401 | 429 | 500
// GET  /api/token                       -> the same, authorised by the session cookie
//
// The exported datasets live in a PRIVATE blob container, so a URL on its own
// fetches nothing. This function is the only thing that holds a storage
// credential: it checks the shared password and returns a container-scoped,
// read-only SAS that expires after TOKEN_HOURS. The storage key never reaches
// the browser.
//
//   base     "https://<account>.blob.core.windows.net/<container>/"
//   sas      query string (no leading "?") the viewer appends to every data request
//   expires  milliseconds since 1970 (the viewer renews a few minutes before this)
//
// App settings (Azure portal -> your Static Web App -> Environment variables):
//   VIEWER_PASSWORD    required. The shared password people type.
//   STORAGE_ACCOUNT    required. Storage account name (the part before .blob.core.windows.net).
//   STORAGE_CONTAINER  required. The private container that holds the exported datasets.
//   STORAGE_KEY        an access key of that account. Preferred: a fresh SAS is minted
//                      for every session.
//   BLOB_SAS           a ready-made read-only SAS for the container, used only when
//                      STORAGE_KEY is not set (for when you may not read account keys).
//                      Every visitor then gets the same SAS until it expires.
//   TOKEN_HOURS        optional, default 1, at most 24. Lifetime of each minted SAS.
//   SESSION_HOURS      optional, default 1, at most 24. How long a browser may renew
//                      without retyping the password.
// There are deliberately no defaults for the account or container: a missing
// setting fails loudly instead of pointing at somebody else's storage.
const crypto = require("crypto");
const {
  StorageSharedKeyCredential,
  generateBlobSASQueryParameters,
  ContainerSASPermissions,
  SASProtocol,
} = require("@azure/storage-blob");

const COOKIE = "sv_session";

function hours(name) {
  const raw = process.env[name];
  if (raw === undefined || raw === "") return 1;
  const v = Number(raw);
  if (!Number.isFinite(v) || v <= 0) return 1;
  return Math.min(v, 24);
}

// Everything this function needs, or the list of what is missing / malformed.
// Setting NAMES are reported to the caller; setting VALUES never are.
function settings() {
  const env = process.env;
  const problems = [];
  const account = (env.STORAGE_ACCOUNT || "").trim();
  const container = (env.STORAGE_CONTAINER || "").trim();
  const password = env.VIEWER_PASSWORD || "";
  const key = (env.STORAGE_KEY || "").trim();
  const presetSas = (env.BLOB_SAS || "").replace(/^\?/, "").trim();
  if (!password) problems.push("VIEWER_PASSWORD is not set");
  if (!account) problems.push("STORAGE_ACCOUNT is not set");
  else if (!/^[a-z0-9]{3,24}$/.test(account)) problems.push("STORAGE_ACCOUNT must be 3-24 lowercase letters or digits");
  if (!container) problems.push("STORAGE_CONTAINER is not set");
  else if (!/^[a-z0-9](?:[a-z0-9]|-(?=[a-z0-9])){2,62}$/.test(container))
    problems.push("STORAGE_CONTAINER must be 3-63 lowercase letters, digits or single hyphens");
  if (!key && !presetSas) problems.push("set STORAGE_KEY (preferred) or BLOB_SAS");
  return {
    ok: problems.length === 0, problems, account, container, password, key, presetSas,
    tokenHours: hours("TOKEN_HOURS"), sessionHours: hours("SESSION_HOURS"),
  };
}

// Brute-force throttle. In memory, so it is per running instance and resets on
// restart: enough to make guessing slow, not a substitute for a good password.
// Deliberately generous: a whole lab behind one campus IP or VPN looks like a
// single address here, and a tripped lockout rejects the right password too.
// 30 tries per 5 minutes is still only ~8,600 guesses a day.
const MAX_FAILS = 30;
const LOCKOUT_MS = 5 * 60 * 1000;
const fails = new Map();

function clientIp(req) {
  const fwd = req.headers["x-forwarded-for"] || "";
  return String(fwd).split(",")[0].trim() || "unknown";
}
function locked(ip) {
  const f = fails.get(ip);
  if (!f) return false;
  if (Date.now() > f.until) { fails.delete(ip); return false; }
  return f.count >= MAX_FAILS;
}
function noteFail(ip) {
  const f = fails.get(ip) || { count: 0, until: 0 };
  f.count += 1;
  f.until = Date.now() + LOCKOUT_MS;
  fails.set(ip, f);
}
// Constant-time comparison, so response timing does not reveal how much of a
// guess was right.
function sameString(a, b) {
  const x = Buffer.from(String(a), "utf8");
  const y = Buffer.from(String(b), "utf8");
  if (x.length !== y.length) return false;         // length alone is not secret here
  return crypto.timingSafeEqual(x, y);
}

// ---- session cookie: "<expiry ms>.<hmac>" ------------------------------------
// It carries nothing but its own expiry; it only says "this browser typed the
// password recently", which lets a long-open tab renew its SAS without asking
// again. The HMAC key folds in the storage credential AND the password, so
// changing VIEWER_PASSWORD (or rotating the key) invalidates every outstanding
// session immediately.
function secret(s) {
  return crypto.createHash("sha256")
    .update("sv-session|" + (s.key || s.presetSas) + "|" + s.password)
    .digest();
}
function signSession(s, expiresAt) {
  const mac = crypto.createHmac("sha256", secret(s)).update(String(expiresAt)).digest("base64url");
  return `${expiresAt}.${mac}`;
}
function validSession(s, value) {
  if (!value) return false;
  const dot = String(value).lastIndexOf(".");
  if (dot < 1) return false;
  const expiresAt = String(value).slice(0, dot);
  if (!/^\d+$/.test(expiresAt) || Number(expiresAt) < Date.now()) return false;
  return sameString(String(value), signSession(s, expiresAt));
}
function readCookie(req, name) {
  const raw = req.headers.cookie || "";
  for (const part of raw.split(";")) {
    const eq = part.indexOf("=");
    if (eq > 0 && part.slice(0, eq).trim() === name) return part.slice(eq + 1).trim();
  }
  return null;
}

function mintSas(s) {
  const base = `https://${s.account}.blob.core.windows.net/${s.container}/`;

  // No account key: hand out the preset SAS. The password still gates it; what is
  // lost is the per-session expiry, because every visitor gets the same token
  // until that SAS lapses. Make it read-only (sp=r) and short-lived.
  if (!s.key) {
    const se = new URLSearchParams(s.presetSas).get("se");
    const expires = se ? Date.parse(se) : NaN;
    return {
      base,
      sas: s.presetSas,
      expires: Number.isFinite(expires) ? expires : Date.now() + s.tokenHours * 3600 * 1000,
    };
  }

  const cred = new StorageSharedKeyCredential(s.account, s.key);
  const expiresOn = new Date(Date.now() + s.tokenHours * 3600 * 1000);
  const sas = generateBlobSASQueryParameters(
    {
      containerName: s.container,
      permissions: ContainerSASPermissions.parse("r"),   // read only: no list, write or delete
      startsOn: new Date(Date.now() - 5 * 60 * 1000),    // 5 min back for clock skew
      expiresOn,
      protocol: SASProtocol.Https,
    },
    cred
  ).toString();
  return { base, sas, expires: expiresOn.getTime() };
}

module.exports = async function (context, req) {
  const noStore = { "Content-Type": "application/json", "Cache-Control": "no-store" };
  const ip = clientIp(req);
  const s = settings();

  if (!s.ok) {
    context.log.error("token: not configured - " + s.problems.join("; "));
    context.res = { status: 500, headers: noStore,
                    body: { error: "server not configured", missing: s.problems } };
    return;
  }

  // A tab that already has a session just gets a fresh SAS.
  if (req.method === "GET") {
    if (!validSession(s, readCookie(req, COOKIE))) {
      context.res = { status: 401, headers: noStore, body: { error: "no session" } };
      return;
    }
    try {
      context.res = { status: 200, headers: noStore, body: mintSas(s) };
    } catch (e) {
      context.log.error("token: SAS failed:", e.message);
      context.res = { status: 500, headers: noStore, body: { error: "could not create an access token" } };
    }
    return;
  }

  if (locked(ip)) {
    context.res = { status: 429, headers: noStore, body: { error: "too many attempts - wait a few minutes" } };
    return;
  }

  const body = req.body || {};
  const given = typeof body === "string"
    ? (() => { try { return JSON.parse(body).password; } catch (e) { return ""; } })()
    : body.password;

  if (!given || !sameString(given, s.password)) {
    noteFail(ip);
    await new Promise((r) => setTimeout(r, 400));     // slow down guessing
    context.res = { status: 401, headers: noStore, body: { error: "wrong password" } };
    return;
  }
  fails.delete(ip);

  let out;
  try {
    out = mintSas(s);
  } catch (e) {
    context.log.error("token: SAS failed:", e.message);
    context.res = { status: 500, headers: noStore, body: { error: "could not create an access token" } };
    return;
  }

  const sessionUntil = Date.now() + Math.round(s.sessionHours * 3600 * 1000);
  const cookie = [
    `${COOKIE}=${signSession(s, sessionUntil)}`,
    "HttpOnly",
    "Secure",
    "SameSite=Strict",
    "Path=/api",
    `Max-Age=${Math.round(s.sessionHours * 3600)}`,
  ].join("; ");

  context.res = { status: 200, headers: { ...noStore, "Set-Cookie": cookie }, body: out };
};
