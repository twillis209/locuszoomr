/* Exercises LZR.layout()'s label-anchor clamping headlessly.
 *
 * Invoked by tests/testthat/test-layout-clamp.R, which passes the path to
 * inst/js/genetrack-relayout.js as argv[2] (resolved in R via system.file(),
 * so it works in both the source tree and an installed package).
 *
 * LZR.measure() uses a canvas 2d context, which node does not have, so it is
 * replaced with a deterministic monospace-style stub: every character is
 * exactly 6px wide at 14px font. That keeps the assertions about clamping
 * geometry independent of real font metrics.
 */
var fs = require('fs');
var path = require('path');

var jsPath = process.argv[2];
if (!jsPath) {
  process.stderr.write('usage: node layout-clamp.js <path to genetrack-relayout.js>\n');
  process.exit(1);
}
require(jsPath);
var LZR = global.LZR;

/* Deterministic width: 6px per character, scaled by font size / 14. */
LZR.measure = function (text, fontPx) {
  return text.length * 6 * (fontPx / 14);
};

var CASES = JSON.parse(fs.readFileSync(
  path.join(__dirname, 'fixtures', 'layout-clamp-cases.json'), 'utf8'));

var out = CASES.map(function (c) {
  var P = {
    tx: c.tx.map(function (g, i) {
      return {
        id: 'G' + i, name: g.name, label: g.name + '&#8594;',
        start: g.start, end: g.end, strand: '+', priority: i + 1,
        hover: g.name
      };
    }),
    ex: [],
    cfg: {
      gapFrac: 0, fontSizePx: 14, maxrows: c.maxrows || 20,
      textPos: 'top', italics: false,
      geneCol: '#000', exonCol: '#000', exonBorder: '#000'
    }
  };
  /* Blank-label genes keep an empty label so the no-label branch is reachable. */
  P.tx.forEach(function (g, i) { if (c.tx[i].name === '') g.label = ''; });

  var res = LZR.layout(P, c.x0, c.x1, c.pxPerData);
  return {
    name: c.name,
    anchors: res.rows.map(function (r) {
      return {gene: P.tx[r.i].name, mid: r.mid, row: r.row};
    })
  };
});

process.stdout.write(JSON.stringify(out));
