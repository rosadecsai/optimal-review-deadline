// Descargador de ICLR 2018-2026 desde OpenReview, para ejecutarse en el contexto
// de una página de openreview.net con la sesión iniciada (la API exige superar
// una verificación a los clientes anónimos). Guarda un JSON compacto por año en
// la carpeta de Descargas: iclr_YYYY_notes.json
//
// Registros compactos (sin texto): por informe se guardan las fechas (cdate,
// tcdate, tmdate, mdate), la firma anónima, el número de envío, los campos
// cortos (puntuaciones, confianza, recomendación...) y, para los campos de texto
// largos, [caracteres, palabras, líneas, viñetas, interrogaciones].

(function () {
  const S = (window.__ICLR = window.__ICLR || {});
  S.status = 'init';
  S.years = S.years || {};
  S.errors = S.errors || [];
  S.log = S.log || [];

  const YEARS = [
    { y: 2018, api: 'v1', prefix: 'ICLR.cc/2018/Conference/-/Paper.*', subs: ['Blind_Submission', 'Withdrawn_Submission', 'Desk_Rejected_Submission'] },
    { y: 2019, api: 'v1', prefix: 'ICLR.cc/2019/Conference/-/Paper.*', subs: ['Blind_Submission', 'Withdrawn_Submission', 'Desk_Rejected_Submission'] },
    { y: 2020, api: 'v1', prefix: 'ICLR.cc/2020/Conference/Paper.*', subs: ['Blind_Submission', 'Withdrawn_Submission', 'Desk_Rejected_Submission'] },
    { y: 2021, api: 'v1', prefix: 'ICLR.cc/2021/Conference/Paper.*', subs: ['Blind_Submission', 'Withdrawn_Submission', 'Desk_Rejected_Submission'] },
    { y: 2022, api: 'v1', prefix: 'ICLR.cc/2022/Conference/Paper.*', subs: ['Blind_Submission', 'Withdrawn_Submission', 'Desk_Rejected_Submission'] },
    { y: 2023, api: 'v1', prefix: 'ICLR.cc/2023/Conference/Paper.*', subs: ['Blind_Submission', 'Withdrawn_Submission', 'Desk_Rejected_Submission'] },
    { y: 2024, api: 'v2' },
    { y: 2025, api: 'v2' },
    { y: 2026, api: 'v2' },
  ];
  S.YEARS = YEARS;

  const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
  const log = (m) => { S.log.push(new Date().toISOString().slice(11, 19) + ' ' + m); if (S.log.length > 400) S.log.shift(); };

  async function getJSON(url, tries = 8) {
    for (let k = 0; k < tries; k++) {
      try {
        const r = await fetch(url, { credentials: 'include' });
        if (r.status === 200) return await r.json();
        if (r.status === 429 || r.status >= 500) { log('HTTP ' + r.status + ' retry ' + k); await sleep(3000 * (k + 1) * (k + 1)); continue; }
        const t = await r.text();
        throw new Error('HTTP ' + r.status + ' ' + t.slice(0, 200));
      } catch (e) {
        if (k === tries - 1) throw e;
        log('fetch error ' + String(e).slice(0, 80) + ' retry ' + k);
        await sleep(3000 * (k + 1));
      }
    }
  }

  function textMetrics(s) {
    const lines = s.split(/\r?\n/).map((x) => x.trim()).filter(Boolean);
    const words = s.split(/\s+/).filter(Boolean).length;
    const bullets = lines.filter((x) => /^([-*•]|\d+[.)])\s/.test(x)).length;
    const q = (s.match(/\?/g) || []).length;
    return [s.length, words, lines.length, bullets, q];
  }

  function compactContent(content) {
    const Sh = {}, L = {};
    for (const k in content || {}) {
      if (k.startsWith('_')) continue;
      let v = content[k];
      if (v && typeof v === 'object' && !Array.isArray(v) && 'value' in v) v = v.value;
      if (v == null) continue;
      if (typeof v === 'string') {
        if (v.length <= 120) Sh[k] = v; else L[k] = textMetrics(v);
      } else if (typeof v === 'number' || typeof v === 'boolean') {
        Sh[k] = v;
      } else if (Array.isArray(v)) {
        const js = JSON.stringify(v);
        if (js.length <= 300) Sh[k] = v; else Sh[k] = '[array:' + v.length + ']';
      }
    }
    return { S: Sh, L };
  }

  const KEEP_FULL = /Official_Review|Meta_Review|Meta_review|Decision|Review$|Meta$/;

  function suffix(inv) { return inv.split('/-/').pop().split('/').pop(); }

  function compactNote(n, api) {
    const invs = api === 'v2' ? (n.invitations || []) : [n.invitation];
    const inv0 = invs.find((x) => KEEP_FULL.test(x)) || invs[0] || '';
    const num = (inv0.match(/(?:Paper|Submission)(\d+)/) || [])[1];
    const rec = {
      i: n.id, f: n.forum, r: n.replyto, v: suffix(inv0), n: num ? +num : undefined,
      s: (n.signatures || []).map((x) => x.split('/').pop()).join('|'),
      c: n.cdate, tc: n.tcdate, tm: n.tmdate,
    };
    if (api === 'v2') { rec.md = n.mdate; rec.od = n.odate; if (invs.length > 1) rec.vs = invs.map(suffix); }
    if (n.ddate) rec.dd = n.ddate;
    if (KEEP_FULL.test(inv0)) {
      const cc = compactContent(n.content);
      rec.S = cc.S; rec.L = cc.L;
    } else {
      // comentarios y otros: solo longitud total del texto
      let len = 0;
      for (const k in n.content || {}) { let v = n.content[k]; if (v && typeof v === 'object' && 'value' in v) v = v.value; if (typeof v === 'string') len += v.length; }
      rec.cl = len;
    }
    return rec;
  }

  function compactSubmission(n, api) {
    const invs = api === 'v2' ? (n.invitations || []) : [n.invitation];
    const c = n.content || {};
    const val = (k) => { let v = c[k]; if (v && typeof v === 'object' && 'value' in v) v = v.value; return v; };
    return {
      i: n.id, n: n.number, c: n.cdate, tc: n.tcdate, tm: n.tmdate, md: n.mdate, od: n.odate,
      vs: invs.map(suffix), venue: val('venue'), venueid: val('venueid'),
      area: val('primary_area') || val('Please_choose_the_closest_area_that_your_submission_falls_into'),
      nrep: n.details && n.details.directReplies ? n.details.directReplies.length : undefined,
    };
  }

  async function gzipText(txt) {
    const cs = new CompressionStream('gzip');
    const w = cs.writable.getWriter();
    w.write(new TextEncoder().encode(txt)); w.close();
    return await new Response(cs.readable).blob();
  }
  function idb() { return new Promise((res, rej) => { const q = indexedDB.open('iclr_dl', 1); q.onupgradeneeded = () => { q.result.createObjectStore('years'); }; q.onsuccess = () => res(q.result); q.onerror = () => rej(q.error); }); }
  async function idbPut(key, val) { const db = await idb(); return new Promise((res, rej) => { const tx = db.transaction('years', 'readwrite'); tx.objectStore('years').put(val, key); tx.oncomplete = () => res(true); tx.onerror = () => rej(tx.error); }); }
  async function idbGet(key) { const db = await idb(); return new Promise((res, rej) => { const tx = db.transaction('years', 'readonly'); const r = tx.objectStore('years').get(key); r.onsuccess = () => res(r.result); r.onerror = () => rej(r.error); }); }
  async function idbKeys() { const db = await idb(); return new Promise((res, rej) => { const tx = db.transaction('years', 'readonly'); const r = tx.objectStore('years').getAllKeys(); r.onsuccess = () => res(r.result); r.onerror = () => rej(r.error); }); }
  S.idbKeys = idbKeys;
  function triggerDownload(name, blob) {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 120000);
  }
  // Guarda el JSON comprimido (gzip) en IndexedDB y lo descarga como .json.gz
  async function saveJSON(name, obj) {
    const txt = JSON.stringify(obj);
    const gz = await gzipText(txt);
    await idbPut(name, gz);
    triggerDownload(name + '.gz', gz);
    return gz.size;
  }
  S.redownload = async function (name) { const gz = await idbGet(name); if (!gz) return 'no data for ' + name; triggerDownload(name + '.gz', gz); return gz.size; };
  S.saveJSON = saveJSON;

  async function downloadYearV1(cfg, st) {
    const notes = [], subs = [];
    const base = 'https://api.openreview.net/notes?';
    for (const inv of cfg.subs) {
      let offset = 0;
      for (;;) {
        const j = await getJSON(base + 'invitation=' + encodeURIComponent('ICLR.cc/' + cfg.y + '/Conference/-/' + inv) + '&limit=1000&offset=' + offset);
        const arr = j.notes || [];
        for (const n of arr) subs.push(compactSubmission(n, 'v1'));
        st.subs = subs.length; st.pages++;
        if (arr.length < 1000) break;
        offset += 1000; await sleep(250);
      }
    }
    let offset = 0, total = null;
    for (;;) {
      const j = await getJSON(base + 'invitation=' + encodeURIComponent(cfg.prefix) + '&limit=1000&offset=' + offset);
      if (total == null) { total = j.count; st.total = total; }
      const arr = j.notes || [];
      for (const n of arr) notes.push(compactNote(n, 'v1'));
      st.notes = notes.length; st.pages++;
      st.reviews = notes.filter((x) => x.v === 'Official_Review').length;
      if (arr.length < 1000) break;
      offset += 1000; await sleep(250);
    }
    return { subs, notes };
  }

  async function downloadYearV2(cfg, st) {
    const notes = [], subs = [];
    const base = 'https://api2.openreview.net/notes?';
    const LIM = 200;
    let offset = 0;
    for (;;) {
      const j = await getJSON(base + 'invitation=' + encodeURIComponent('ICLR.cc/' + cfg.y + '/Conference/-/Submission') + '&details=directReplies&limit=' + LIM + '&offset=' + offset);
      if (st.total == null) st.total = j.count;
      const arr = j.notes || [];
      for (const n of arr) {
        subs.push(compactSubmission(n, 'v2'));
        for (const rep of (n.details && n.details.directReplies) || []) notes.push(compactNote(rep, 'v2'));
      }
      st.subs = subs.length; st.notes = notes.length; st.pages++;
      st.reviews = notes.filter((x) => x.v === 'Official_Review').length;
      if (arr.length < LIM) break;
      offset += LIM; await sleep(250);
    }
    return { subs, notes };
  }

  async function downloadYear(cfg) {
    const st = (S.years[cfg.y] = { status: 'running', pages: 0, subs: 0, notes: 0, reviews: 0, total: null, started: Date.now() });
    try {
      const res = cfg.api === 'v1' ? await downloadYearV1(cfg, st) : await downloadYearV2(cfg, st);
      const out = { year: cfg.y, api: cfg.api, fetched_at: new Date().toISOString(), config: cfg, n_subs: res.subs.length, n_notes: res.notes.length, subs: res.subs, notes: res.notes };
      st.bytes = await saveJSON('iclr_' + cfg.y + '_notes.json', out);
      await sleep(15000); // dejar que la descarga se complete antes de seguir
      st.status = 'saved'; st.finished = Date.now();
      log('year ' + cfg.y + ' saved: ' + res.notes.length + ' notes, ' + res.subs.length + ' subs, ' + st.bytes + ' bytes');
    } catch (e) {
      st.status = 'error'; st.error = String(e); S.errors.push(cfg.y + ': ' + String(e)); log('year ' + cfg.y + ' ERROR ' + String(e).slice(0, 200));
    }
  }

  S.runAll = async function (years) {
    S.status = 'running';
    const list = YEARS.filter((c) => !years || years.includes(c.y));
    for (const cfg of list) await downloadYear(cfg);
    S.status = 'done';
    log('all done');
  };
  S.summary = function () {
    const o = { status: S.status, errors: S.errors, years: {} };
    for (const y in S.years) { const s = S.years[y]; o.years[y] = { st: s.status, pages: s.pages, subs: s.subs, notes: s.notes, reviews: s.reviews, total: s.total, bytes: s.bytes, mins: ((s.finished || Date.now()) - s.started) / 60000 }; }
    o.lastLog = S.log.slice(-4);
    return o;
  };
  S.status = 'ready';
})();
