"""Local editor for the pre-registration: edit the text in place and pin comments to it.

    python3 tools/prereg_editor.py [TEMPLATE [RESULTS]]      then open http://localhost:8768

The page shows the draft as rendered, one block per heading and per section body, with every number from the results
file (work/results.json); a number it does not have shows as pending. Saving writes each edited block back into the
template (draft/preregistration_v2.src.md), the file the draft is rendered from: a number whose text the edit left
untouched goes back as its placeholder, and a number changed or added stays as typed, which the page names so that it
can be tied to the pipeline before the next render. Each earlier version is kept in local/prereg_history/, and a save is
refused if the template changed on disk since the page was loaded (reload first). Comments go to
local/prereg_comments.json, mirrored as local/prereg_comments.md. The template's opening comment is not shown. Passages
listed in local/prereg_highlights.json (rendered text, one string each) are marked with a highlighter on the page only:
saving strips the marks, so they never reach the template."""
import collections
import difflib
import hashlib
import html
import json
import re
import sys
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools")); import render as RD; sys.path.pop(0)
DOC, NUMBERS = ROOT / "draft/preregistration_v2.src.md", ROOT / "work/results.json"
NOTES, NOTES_MD, HISTORY = ROOT / "local/prereg_comments.json", ROOT / "local/prereg_comments.md", ROOT / "local/prereg_history"
HIGHLIGHTS = ROOT / "local/prereg_highlights.json"
LEAD = re.compile(r"\A\s*<!--.*?-->\s*", re.S)        # the template's opening comment: kept, not shown
PORT = 8768


def blocks(src):
    """[(kind, level, text)]: a heading line (kind 'h') or the body between headings (kind 'b')."""
    out, body = [], []
    for line in src.split("\n"):
        m = re.match(r"^(#{1,4}) (.*)$", line)
        if m:
            out.append(("b", 0, "\n".join(body).strip())); body = []
            out.append(("h", len(m.group(1)), m.group(2).strip()))
        else:
            body.append(line)
    out.append(("b", 0, "\n".join(body).strip()))
    return out


def join(bl):
    parts = []
    for kind, lev, text in bl:
        if kind == "h":
            parts.append("#" * lev + " " + text)
        elif text:
            parts.append(text)
    return "\n\n".join(parts).rstrip() + "\n"


def render_block(text, reg):
    """The block as the draft shows it, the spans of the numbers' output in it, and each span's placeholder."""
    out, spans, srcs, pos = [], [], [], 0
    for m in RD.PH.finditer(text):
        out.append(text[pos:m.start()]); start = sum(map(len, out))
        try:
            s = RD.fmt(RD.lookup(reg, m.group(1)), m.group(2))
        except (KeyError, IndexError, TypeError, ValueError):
            s = f"[{m.group(1)} pending]"
        out.append(s); spans.append((start, start + len(s))); srcs.append(m.group(0)); pos = m.end()
    out.append(text[pos:])
    return "".join(out), spans, srcs


def unrender(old, spans, srcs, new):
    """The template text of an edited block: each placeholder whose output lies in text the edit left unchanged goes back
    in its place; a number the edit touched stays as typed."""
    eq = [(i1, i2, j1) for tag, i1, i2, j1, _ in difflib.SequenceMatcher(None, old, new, autojunk=False).get_opcodes() if tag == "equal"]
    reps = [(j1 + s - i1, j1 + e - i1, src) for (s, e), src in zip(spans, srcs) for i1, i2, j1 in eq if i1 <= s and e <= i2]
    out, pos = [], 0
    for a, b, src in sorted(reps):
        out += [new[pos:a], src]; pos = b
    return "".join(out + [new[pos:]])


def registry():
    return json.loads(NUMBERS.read_text())["values"] if NUMBERS.exists() else {}


def typed_numbers(old, old_spans, new, new_spans):
    """The numbers in the edited text, outside the placeholders' output, that the text before the edit did not have."""
    def outside(t, spans):
        keep, pos = [], 0
        for a, b in spans:
            keep.append(t[pos:a]); pos = b
        return re.findall(r"\d[\d,.]*\d|\d", " ".join(keep + [t[pos:]]))
    return list((collections.Counter(outside(new, new_spans)) - collections.Counter(outside(old, old_spans))).elements())


def inline(t):
    t = html.escape(t, quote=False)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
    return re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", t)


def render_body(text):
    """Paragraphs, bullet lists (one level of nesting, continuation lines joined), and tables (edited as plain markdown)."""
    out = []
    for para in re.split(r"\n\s*\n", text):
        lines = [ln for ln in para.split("\n") if ln.strip()]
        if not lines:
            continue
        if all(ln.lstrip().startswith("|") for ln in lines):
            out.append(f'<pre class="raw">{html.escape(para.strip())}</pre>'); continue
        if re.match(r"^\d+\. ", lines[0]) and all(re.match(r"^(\d+\. |\s)", ln) for ln in lines):   # a numbered list
            items = []
            for ln in lines:
                m = re.match(r"^(\d+)\. (.*)$", ln)
                if m:
                    items.append([m.group(1), m.group(2).strip()])
                else:
                    items[-1][1] += " " + ln.strip()
            out.append("<ol>" + "".join(f'<li value="{n}">{inline(t)}</li>' for n, t in items) + "</ol>"); continue
        prose, items = [], []
        def flush_prose():
            if prose:
                out.append(f"<p>{inline(' '.join(prose))}</p>"); prose.clear()
        def flush_items():
            if items:
                out.append("<ul>" + "".join(f'<li{" class=\"sub\"" if sub else ""}>{inline(" ".join(t))}</li>' for sub, t in items) + "</ul>"); items.clear()
        for ln in lines:
            m = re.match(r"^(\s*)- (.*)$", ln)
            if m:
                flush_prose(); items.append((len(m.group(1)) > 0, [m.group(2).strip()]))
            elif items and ln.startswith(" "):
                items[-1][1].append(ln.strip())
            else:
                flush_items(); prose.append(ln.strip())
        flush_prose(); flush_items()
    return "".join(out) or "<p></p>"


def mark(h, hl):
    """The page's HTML with each highlighted passage wrapped in <mark>."""
    for x in hl:
        e = html.escape(x, quote=False)
        h = h.replace(e, f"<mark>{e}</mark>")
    return h


def page():
    src = DOC.read_text(); h = hashlib.sha256(src.encode()).hexdigest(); parts = []; reg = registry()
    hl = json.loads(HIGHLIGHTS.read_text()) if HIGHLIGHTS.exists() else []
    for i, (kind, lev, tpl) in enumerate(blocks(src)):
        if i == 0 and kind == "b":
            tpl = LEAD.sub("", tpl)
        text = render_block(tpl, reg)[0]
        if kind == "h":
            parts.append(f'<h{lev} class="slot head" data-slot="{i}" data-head="{html.escape(text)}">{inline(text)}</h{lev}>')
        elif text or i:
            parts.append(f'<div class="slot body" data-slot="{i}">{mark(render_body(text), hl)}</div>')
    return TEMPLATE.replace("{TITLE}", html.escape(DOC.name)).replace("{HASH}", h).replace("{BODY}", "\n".join(parts))


def load_notes():
    return json.loads(NOTES.read_text()) if NOTES.exists() else []


def save_notes(notes):
    NOTES.write_text(json.dumps(notes, indent=1, ensure_ascii=False))
    lines = ["# Comments on the pre-registration draft", "", "Open comments first. Each names the section it is pinned to.", ""]
    for n in sorted(notes, key=lambda n: (n["status"] != "open", n["id"])):
        lines.append(f'## {n["id"]}. [{n["status"]}] in "{n["anchor"]}"')
        if n.get("snippet"):
            lines.append(f'> {n["snippet"]}')
        lines += [n["text"], f'_{n["time"]}_'] + ([f'Reply: {n["reply"]}'] if n.get("reply") else []) + [""]
    NOTES_MD.write_text("\n".join(lines))


def save(edits, h):
    src = DOC.read_text()
    if hashlib.sha256(src.encode()).hexdigest() != h:
        raise ValueError("the file changed on disk since this page was loaded: copy your edits, reload, and save again")
    bl = blocks(src); reg = registry(); typed = []
    for k, text in edits.items():
        kind, lev, tpl = bl[int(k)]
        lead = LEAD.match(tpl).group(0) if int(k) == 0 and kind == "b" and LEAD.match(tpl) else ""
        old, spans, srcs = render_block(tpl[len(lead):], reg)
        new_tpl = unrender(old, spans, srcs, " ".join(text.split()) if kind == "h" else text.strip())
        bl[int(k)] = (kind, lev, (lead.rstrip() + "\n\n" if lead else "") + new_tpl)
        shown, sp, _ = render_block(new_tpl, reg)
        typed += typed_numbers(old, spans, shown, sp)
    new = join(bl)
    HISTORY.mkdir(parents=True, exist_ok=True); stamp = f"{datetime.now():%Y%m%d_%H%M%S}"
    (HISTORY / f"{DOC.stem}_{stamp}.md").write_text(src)
    DOC.write_text(new)
    (HISTORY / f"{DOC.stem}_{stamp}_saved.md").write_text(new)
    return hashlib.sha256(new.encode()).hexdigest(), typed


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def send(self, body, ctype="application/json", code=200):
        body = body.encode() if isinstance(body, str) else body
        self.send_response(code); self.send_header("Content-Type", ctype); self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)

    def do_GET(self):
        if self.path.split("?")[0] in ("/", "/index.html"):
            self.send(page(), "text/html; charset=utf-8")
        elif self.path == "/api/notes":
            self.send(json.dumps(load_notes()))
        else:
            self.send(json.dumps({"error": "not found"}), code=404)

    def do_POST(self):
        data = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
        try:
            if self.path == "/api/save":
                h, typed = save(data["edits"], data["hash"])
                self.send(json.dumps({"ok": True, "saved": len(data["edits"]), "hash": h, "typed": typed}))
            elif self.path == "/api/notes":
                notes = load_notes()
                note = {k: data.get(k) for k in ("anchor", "snippet", "text")}
                note.update(id=max([n["id"] for n in notes] + [0]) + 1, status="open", time=f"{datetime.now():%Y-%m-%d %H:%M}")
                notes.append(note); save_notes(notes); self.send(json.dumps(note))
            elif self.path == "/api/notes/status":
                notes = load_notes()
                for n in notes:
                    if n["id"] == data["id"]:
                        n["status"] = data["status"]
                save_notes(notes); self.send(json.dumps({"ok": True}))
            else:
                self.send(json.dumps({"error": "unknown"}), code=404)
        except Exception as e:                      # shown on the page rather than failing silently
            self.send(json.dumps({"error": str(e)}), code=400)


TEMPLATE = r"""<!doctype html><html lang="en-GB"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{TITLE}</title>
<style>
  :root { --ink: #1d232a; --muted: #5d6670; --paper: #fbfbf8; --rule: #dfe3e6; --blue: #2b5c8a; --amber: #d89b1f; --pin: #a4471c; }
  body { margin: 0; background: var(--paper); color: var(--ink); font: 18px/1.62 Georgia, "Iowan Old Style", serif; }
  main { max-width: 44rem; margin: 0 auto; padding: 3rem 1.25rem 8rem; }
  h1 { font-size: 2rem; line-height: 1.2; margin: 0 0 1.5rem; } h2 { font-size: 1.5rem; margin: 2.6rem 0 .6rem; }
  h3 { font-size: 1.2rem; margin: 2rem 0 .4rem; } h4 { font-size: 1rem; margin: 1.5rem 0 .3rem; }
  p { margin: 0 0 1rem; } ul, ol { margin: 0 0 1rem; padding-left: 1.6rem; } li { margin-bottom: .35rem; } li.sub { margin-left: 1.4rem; list-style: circle; }
  code { font-size: .85em; background: #eef0f1; padding: 0 .2em; border-radius: 3px; }
  pre.raw { font: 14px/1.5 ui-monospace, Menlo, monospace; background: #f1f3f4; padding: .8rem; border-radius: 6px; overflow-x: auto; white-space: pre; }
  mark { background: linear-gradient(transparent 12%, #ffe66d 12%, #ffe66d 88%, transparent 88%); color: inherit; padding: 0 1px; }
  body.nohl mark { background: none; }
  .slot { outline: 1px dashed transparent; outline-offset: 4px; border-radius: 2px; }
  body.ed-text .slot:hover { outline-color: var(--blue); }
  body.ed-text .slot:focus { outline: 2px solid var(--blue); background: rgba(43,92,138,.04); }
  .slot.dirty { outline: 2px solid var(--amber) !important; }
  body.ed-note .slot { cursor: crosshair; } body.ed-note .ed-hover { outline: 2px solid var(--pin) !important; }
  #bar { position: fixed; right: 16px; bottom: 16px; display: flex; gap: 6px; align-items: center; background: #1d232a; color: #fff;
    padding: 8px; border-radius: 10px; font: 14px/1.2 system-ui, sans-serif; box-shadow: 0 6px 24px rgba(0,0,0,.25); z-index: 10; }
  #bar button { font: inherit; color: #fff; background: #333c45; border: 0; border-radius: 6px; padding: 7px 11px; cursor: pointer; }
  #bar button[aria-pressed="true"] { background: var(--blue); } #bar #save.ready { background: var(--amber); color: #1d232a; }
  #bar button:focus-visible { outline: 2px solid #fff; }
  #status { padding: 0 6px; color: #c3cbd2; min-width: 7rem; }
  .pin { position: absolute; width: 24px; height: 24px; margin: -12px 0 0 -12px; border-radius: 50%; background: var(--pin); color: #fff;
    font: 600 12px/24px system-ui, sans-serif; text-align: center; box-shadow: 0 0 0 2px #fff; cursor: pointer; z-index: 9; }
  .pin.resolved { background: #8b949c; opacity: .6; }
  .pop, #list { position: absolute; z-index: 11; background: #fff; color: var(--ink); border-radius: 10px; box-shadow: 0 8px 30px rgba(0,0,0,.25);
    padding: 12px; font: 14px/1.45 system-ui, sans-serif; }
  .pop { width: 320px; } .pop textarea { width: 100%; box-sizing: border-box; min-height: 100px; font: inherit; padding: 8px; border: 1px solid #cfd6dc; border-radius: 6px; }
  .pop .where, #list .a { font-size: 12px; color: var(--muted); margin-bottom: 6px; }
  .pop .row { display: flex; gap: 6px; justify-content: flex-end; margin-top: 8px; }
  .pop button, #list button { font: inherit; border: 0; border-radius: 6px; padding: 5px 11px; cursor: pointer; background: #e7ebee; }
  .pop button.primary { background: var(--blue); color: #fff; }
  #list { position: fixed; right: 16px; bottom: 72px; width: min(380px, calc(100vw - 32px)); max-height: 60vh; overflow: auto; }
  #list .n { border-bottom: 1px solid var(--rule); padding: 8px 0; } #list .n.resolved { opacity: .5; }
  @media (max-width: 600px) { body { font-size: 17px; } #bar { left: 16px; flex-wrap: wrap; } }
</style></head><body>
<main id="doc" data-hash="{HASH}">
{BODY}
</main>
<div id="bar" role="toolbar" aria-label="Editor">
  <button id="t" aria-pressed="true" title="Click any text to edit it">Edit text</button>
  <button id="c" aria-pressed="false" title="Click a section to pin a comment to it">Comment</button>
  <button id="save" title="Save text edits (Cmd+S)">Save</button>
  <button id="show" title="Show all comments">Comments</button>
  <button id="hl" aria-pressed="true" title="Show or hide the highlighted additions">Highlights</button>
  <span id="status" aria-live="polite">Ready</span>
</div>
<script>
(() => {
  const $ = s => document.querySelector(s), $$ = s => [...document.querySelectorAll(s)];
  const status = t => { $("#status").textContent = t; };
  const dirty = new Set(); let mode = "text", notes = [], hash = $("#doc").dataset.hash;
  const md = h => h.replace(/<code>(.*?)<\/code>/g, "`$1`").replace(/<(strong|b)>(.*?)<\/\1>/g, "**$2**").replace(/<(em|i)>(.*?)<\/\1>/g, "*$2*")
    .replace(/<[^>]+>/g, "").replace(/&nbsp;/g, " ").replace(/&lt;/g, "<").replace(/&gt;/g, ">").replace(/&amp;/g, "&").replace(/\s+/g, " ").trim();
  const paras = h => h.replace(/<br\s*\/?>/gi, "\n\n").replace(/<\/(div|p)>/gi, "\n\n").split(/\n\s*\n/).map(md).filter(Boolean);
  function text(el) {
    if (el.classList.contains("head")) return md(el.innerHTML);
    const out = [];
    for (const c of el.children) {
      if (c.tagName === "UL") out.push([...c.children].map(li => (li.classList.contains("sub") ? "   - " : "- ") + md(li.innerHTML)).join("\n"));
      else if (c.tagName === "OL") out.push([...c.children].map((li, k) => (li.value || k + 1) + ". " + md(li.innerHTML)).join("\n"));
      else if (c.tagName === "PRE") out.push(c.innerText.trim());
      else out.push(...paras(c.innerHTML));
    }
    return out.join("\n\n");
  }
  function setMode(m) {
    mode = m; document.body.classList.toggle("ed-text", m === "text"); document.body.classList.toggle("ed-note", m === "note");
    $("#t").setAttribute("aria-pressed", m === "text"); $("#c").setAttribute("aria-pressed", m === "note");
    $$(".slot").forEach(s => s.contentEditable = m === "text" ? "true" : "false");
    status(m === "text" ? "Click text to edit" : "Click a section to comment");
  }
  document.addEventListener("input", e => { const s = e.target.closest && e.target.closest(".slot");
    if (s) { dirty.add(s.dataset.slot); s.classList.add("dirty"); $("#save").classList.add("ready"); status(dirty.size + " unsaved"); } });
  async function save() {
    if (!dirty.size) { status("Nothing to save"); return; }
    const edits = {}; for (const k of dirty) edits[k] = text($(`.slot[data-slot="${k}"]`));
    status("Saving…");
    const r = await fetch("/api/save", { method: "POST", body: JSON.stringify({ edits, hash }) }).then(r => r.json());
    if (r.error) { status("Not saved"); alertBox(r.error); return; }
    hash = r.hash; dirty.clear(); $$(".slot.dirty").forEach(s => s.classList.remove("dirty")); $("#save").classList.remove("ready");
    status("Saved " + r.saved);
    if (r.typed && r.typed.length) alertBox("Saved. Numbers typed by hand, which Claude will tie to the pipeline before the next render: " + r.typed.join(", "));
  }
  function alertBox(msg) { const d = document.createElement("div"); d.className = "pop"; d.style.position = "fixed"; d.style.right = "16px"; d.style.bottom = "72px";
    d.innerHTML = `<div>${msg.replace(/</g, "&lt;")}</div><div class="row"><button>OK</button></div>`; d.querySelector("button").onclick = () => d.remove(); document.body.appendChild(d); }
  $("#save").onclick = save;
  document.addEventListener("keydown", e => { if ((e.metaKey || e.ctrlKey) && e.key === "s") { e.preventDefault(); save(); } });
  addEventListener("beforeunload", e => { if (dirty.size) { e.preventDefault(); e.returnValue = ""; } });
  $("#t").onclick = () => setMode("text"); $("#c").onclick = () => setMode(mode === "note" ? "text" : "note");
  const section = el => { let n = el; while (n && !(n.classList && n.classList.contains("head"))) n = n.previousElementSibling; return n ? n.dataset.head : "top"; };
  function anchorOf(el) { const s = el.closest(".slot"); if (!s) return null;
    const head = s.classList.contains("head") ? s.dataset.head : section(s);
    const sel = getSelection().toString().trim();
    return { node: s, anchor: head, snippet: (sel || s.textContent.trim()).slice(0, 160) }; }
  const locate = n => { const h = $$(".head").find(x => x.dataset.head === n.anchor); return h ? (h.nextElementSibling && h.nextElementSibling.classList.contains("body") ? h.nextElementSibling : h) : null; };
  function pins() { $$(".pin").forEach(p => p.remove());
    const count = {}; notes.forEach(n => { const el = locate(n); if (!el) return; const r = el.getBoundingClientRect(), k = n.anchor; count[k] = (count[k] || 0) + 1;
      const p = document.createElement("div"); p.className = "pin" + (n.status === "open" ? "" : " resolved"); p.textContent = n.id; p.title = n.text;
      p.style.left = (scrollX + r.right + 18) + "px"; p.style.top = (scrollY + r.top + 12 + 28 * (count[k] - 1)) + "px";
      p.onclick = e => { e.stopPropagation(); list(n.id); }; document.body.appendChild(p); }); }
  let hov = null;
  document.addEventListener("mousemove", e => { if (mode !== "note") return; const a = anchorOf(e.target), n = a && a.node;
    if (hov !== n) { hov && hov.classList.remove("ed-hover"); hov = n; n && n.classList.add("ed-hover"); } });
  document.addEventListener("mouseup", e => {
    if (mode !== "note" || e.target.closest("#bar, .pop, #list, .pin")) return;
    const a = anchorOf(e.target); if (!a) return;
    document.querySelector(".pop") && document.querySelector(".pop").remove();
    const pop = document.createElement("div"); pop.className = "pop";
    pop.innerHTML = `<div class="where">In “${a.anchor.replace(/</g, "&lt;")}”: ${a.snippet.slice(0, 90).replace(/</g, "&lt;")}…</div><textarea placeholder="Your comment for Claude"></textarea><div class="row"><button data-x>Cancel</button><button class="primary" data-ok>Add comment</button></div>`;
    pop.style.left = Math.min(e.pageX + 8, scrollX + innerWidth - 340) + "px"; pop.style.top = (e.pageY + 8) + "px";
    document.body.appendChild(pop); pop.querySelector("textarea").focus();
    pop.querySelector("[data-x]").onclick = () => pop.remove();
    pop.querySelector("[data-ok]").onclick = async () => { const t = pop.querySelector("textarea").value.trim(); if (!t) return;
      const n = await fetch("/api/notes", { method: "POST", body: JSON.stringify({ anchor: a.anchor, snippet: a.snippet, text: t }) }).then(r => r.json());
      notes.push(n); pop.remove(); pins(); status("Comment " + n.id + " added"); };
  });
  function list(focus) { let box = $("#list"); if (box && focus === undefined) { box.remove(); return; } box && box.remove();
    box = document.createElement("div"); box.id = "list"; box.innerHTML = notes.length ? "" : "<p>No comments yet. Click Comment, then click a section (select words first to quote them).</p>";
    [...notes].sort((a, b) => (a.status !== "open") - (b.status !== "open") || a.id - b.id).forEach(n => {
      const d = document.createElement("div"); d.className = "n" + (n.status === "open" ? "" : " resolved");
      d.innerHTML = `<b>${n.id}.</b> ${n.text.replace(/</g, "&lt;")}` + (n.reply ? `<div class="a">Reply: ${n.reply.replace(/</g, "&lt;")}</div>` : "") + `<div class="a">${n.anchor.replace(/</g, "&lt;")}, ${n.status}</div>`;
      const go = document.createElement("button"); go.textContent = "Show"; go.onclick = () => { const el = locate(n); el && el.scrollIntoView({ block: "center" }); setTimeout(pins, 400); };
      const st = document.createElement("button"); st.textContent = n.status === "open" ? "Resolve" : "Reopen";
      st.onclick = async () => { n.status = n.status === "open" ? "resolved" : "open"; await fetch("/api/notes/status", { method: "POST", body: JSON.stringify({ id: n.id, status: n.status }) }); pins(); list(n.id); };
      d.append(go, " ", st); box.appendChild(d); if (n.id === focus) setTimeout(() => d.scrollIntoView({ block: "nearest" }), 0); });
    document.body.appendChild(box); }
  $("#show").onclick = () => list();
  $("#hl").onclick = () => { const off = document.body.classList.toggle("nohl"); $("#hl").setAttribute("aria-pressed", !off); };
  addEventListener("resize", pins);
  fetch("/api/notes").then(r => r.json()).then(n => { notes = n; setTimeout(pins, 300); });
  setMode("text");
})();
</script></body></html>"""

if __name__ == "__main__":
    if len(sys.argv) > 1:
        DOC = ROOT / sys.argv[1]
    if len(sys.argv) > 2:
        NUMBERS = ROOT / sys.argv[2]
    srv = ThreadingHTTPServer(("127.0.0.1", PORT), Handler)
    print(f"Editor running at http://localhost:{PORT} (edits go to {DOC.relative_to(ROOT)}, numbers from "
          f"{NUMBERS.relative_to(ROOT) if NUMBERS.is_relative_to(ROOT) else NUMBERS}, comments to local/prereg_comments.md)")
    srv.serve_forever()
