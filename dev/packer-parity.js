/* Runs the JS packer over the shared fixtures and prints JSON to stdout.
 * Invoked by tests/testthat/test-packer-parity.R. */
var fs = require('fs');
var path = require('path');

var root = path.join(__dirname, '..');
require(path.join(root, 'inst', 'js', 'genetrack-relayout.js'));

var cases = JSON.parse(fs.readFileSync(
  path.join(root, 'tests', 'testthat', 'fixtures', 'packer-cases.json'), 'utf8'));

var out = cases.map(function (c) {
  var items = c.min.map(function (m, i) { return {min: m, max: c.max[i]}; });
  return {name: c.name, rows: global.LZR.packRows(items)};
});

process.stdout.write(JSON.stringify(out));
