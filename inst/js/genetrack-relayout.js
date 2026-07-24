/* Client-side gene track re-layout for locuszoomr.
 *
 * Injected into plotly widgets by add_genetrack_relayout(). Re-packs the gene
 * annotation panel against the visible window whenever the user zooms or pans.
 * Plain ES5 on purpose: this file ships as package source with no build step.
 */
(function (global) {
  'use strict';

  var LZR = global.LZR || {};

  /* Greedy first-fit packing. items: [{min, max}] in priority order.
   * Returns an array of 1-based row numbers, mirroring R's pack_rows(). */
  LZR.packRows = function (items) {
    var rows = [];
    var out = new Array(items.length);
    for (var i = 0; i < items.length; i++) {
      var it = items[i];
      var placed = false;
      for (var j = 0; j < rows.length && !placed; j++) {
        var clash = false;
        for (var k = 0; k < rows[j].length; k++) {
          var h = rows[j][k];
          if (it.min < h.max && it.max > h.min) { clash = true; break; }
        }
        if (!clash) { rows[j].push(it); out[i] = j + 1; placed = true; }
      }
      if (!placed) { rows.push([it]); out[i] = rows.length; }
    }
    return out;
  };

  global.LZR = LZR;
})(typeof window !== 'undefined' ? window : global);
