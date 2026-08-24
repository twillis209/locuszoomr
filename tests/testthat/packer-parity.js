/* Runs the JS packer over the shared fixtures and prints JSON to stdout.
 * Invoked by tests/testthat/test-packer-parity.R, which passes the path to
 * genetrack-relayout.js as argv[2] (resolved in R via system.file() so it
 * works both in the source tree and in an installed package, where inst/js
 * is flattened to js/). */
var fs = require('fs');
var path = require('path');

var jsPath = process.argv[2];
require(jsPath);

var cases = JSON.parse(fs.readFileSync(
  path.join(__dirname, 'fixtures', 'packer-cases.json'), 'utf8'));

var out = cases.map(function (c) {
  var items = c.min.map(function (m, i) { return {min: m, max: c.max[i]}; });
  return {name: c.name, rows: global.LZR.packRows(items)};
});

process.stdout.write(JSON.stringify(out));
