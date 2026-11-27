/*
 * Review due dates of TMLR (Transactions on Machine Learning Research) on OpenReview.
 *
 * Paste in the browser console of a tab opened on https://openreview.net while signed in
 * (the API answers anonymous clients with a verification challenge). Ten requests, one per second.
 * It lists the per-paper review invitations TMLR/Paper<n>/-/Review, expired ones included, and saves
 * tmlr_review_duedates.csv with: number, due_ms (due date, milliseconds since the epoch, UTC),
 * end_of_day (1 when the due date is 23:59:59 UTC, the value the journal workflow sets; 0 when it was set by hand).
 *
 * Read on 27 November 2026: 5,651 invitations, papers 1 to 6806 (put under review up to December 2025).
 * From January 2026 the per-paper invitations are no longer listed for signed-in users outside the journal.
 * No reviewer identity is requested or stored.
 */
(async () => {
  const api = 'https://api2.openreview.net';
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  const rows = []; let off = 0;
  for (;;) {
    const r = await fetch(api + '/invitations?invitation=' + encodeURIComponent('TMLR/-/Review') + '&expired=true&limit=1000&offset=' + off, { credentials: 'include' });
    if (!r.ok) throw new Error('http ' + r.status);
    const iv = (await r.json()).invitations || [];
    for (const i of iv) {
      const m = /^TMLR\/Paper(\d+)\/-\/Review$/.exec(i.id);
      if (!m || i.duedate == null) continue;
      const t = i.duedate % 86400000;
      const eod = (t === 86399999 || t === 86399000) ? 1 : 0;
      // end-of-day dates are stored as the last millisecond of the day; dates set by hand are kept to the minute
      const due = eod ? (Math.floor(i.duedate / 86400000) + 1) * 86400000 - 1 : Math.floor(i.duedate / 60000) * 60000;
      rows.push([+m[1], due, eod]);
    }
    if (iv.length < 1000) break; off += 1000; await sleep(1200);
  }
  rows.sort((a, b) => a[0] - b[0]);
  const csv = 'number,due_ms,end_of_day\n' + rows.map(r => r.join(',')).join('\n') + '\n';
  const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv' })); a.download = 'tmlr_review_duedates.csv';
  document.body.appendChild(a); a.click(); a.remove();
  console.log(rows.length + ' due dates');
})();
