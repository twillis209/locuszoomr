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

  /* Convert a trace axis id ("x", "x2") to its layout key ("xaxis", "xaxis2"). */
  LZR.axisKey = function (id) {
    return id.charAt(0) + 'axis' + id.slice(1);
  };

  var measureCtx = null;
  LZR.measure = function (text, fontPx) {
    if (measureCtx === null) {
      measureCtx = document.createElement('canvas').getContext('2d');
    }
    measureCtx.font = fontPx + 'px "Open Sans", verdana, arial, sans-serif';
    return measureCtx.measureText(text).width;
  };

  /* Pack the genes visible in [x0, x1] into rows.
   * Returns {rows: [...], nrow: n, hidden: n} where rows entries carry the
   * gene index, its assigned row, and its label anchor. */
  LZR.layout = function (P, x0, x1, pxPerData) {
    var span = x1 - x0;
    var gap = span * P.cfg.gapFrac;
    var margin = span * 0.1;
    var lo = x0 - margin, hi = x1 + margin;

    var items = [];
    for (var i = 0; i < P.tx.length; i++) {
      var g = P.tx[i];
      if (g.end < lo || g.start > hi) continue;
      /* mirror R: strwidth(paste0("--", gene_name)) */
      var halfw = 0;
      if (g.label !== '') {
        halfw = LZR.measure('--' + g.name, P.cfg.fontSizePx) / pxPerData / 2;
      }
      var mid = (g.start + g.end) / 2;
      items.push({
        i: i, mid: mid, priority: g.priority,
        min: Math.min(g.start, g.end, mid - halfw) - gap / 2,
        max: Math.max(g.start, g.end, mid + halfw) + gap / 2
      });
    }
    items.sort(function (a, b) { return a.priority - b.priority; });

    var assigned = LZR.packRows(items);
    var nrow = 0;
    for (var k = 0; k < assigned.length; k++) {
      items[k].row = assigned[k];
      if (assigned[k] > nrow) nrow = assigned[k];
    }

    var kept = [], hidden = 0;
    for (var m = 0; m < items.length; m++) {
      if (items[m].row <= P.cfg.maxrows) kept.push(items[m]); else hidden++;
    }
    return {rows: kept, nrow: nrow, hidden: hidden};
  };

  LZR.attach = function (el, P) {
    var gd = el;
    var dead = false;

    function fail(err) {
      dead = true;
      if (global.console) global.console.warn('locuszoomr re-layout disabled:', err);
    }

    /* Everything below can throw synchronously (e.g. a malformed payload),
     * and an uncaught throw here would escape the onRender hook entirely,
     * abort remaining htmlwidgets afterRender hooks on this element, and
     * leave no re-layout listener attached. Fail closed instead: warn once
     * and leave the widget exactly as R rendered it. */
    try {
      var Plotly = global.Plotly;
      var xkey = LZR.axisKey(P.idx.xaxis);
      var yref = P.idx.yaxis;

      /* Shapes that are not ours, captured once. Every update rebuilds the
       * array as base.concat(exons), so stored indices never go stale.
       * shapeIdx serialises as a bare number, not a length-1 array, when the
       * gene panel has exactly one shape (jsonlite auto_unbox); normalise to
       * an array before relying on .indexOf(). */
      var idxs = [].concat(P.idx.shapeIdx || []);
      var all = (gd.layout.shapes || []);
      var base = [];
      for (var i = 0; i < all.length; i++) {
        if (idxs.indexOf(i) === -1) base.push(all[i]);
      }
      var baseAnn = (gd.layout.annotations || []).slice();

      var busy = false, timer = null, pending = false, lastRange = null;

      function schedule() {
        if (timer) global.clearTimeout(timer);
        timer = global.setTimeout(function () {
          try { apply(); } catch (err) { fail(err); }
        }, 100);
      }

      function apply() {
        var ax = gd._fullLayout[xkey];
        var rng = ax.range;
        /* No-op guard: dragmode toggles, resizes and legend clicks all fire
         * plotly_relayout without changing the x range. Skip the re-pack (and
         * the label-suppression recompute that goes with it) when the range
         * is unchanged, so those events are visually silent. */
        if (lastRange !== null && rng[0] === lastRange[0] && rng[1] === lastRange[1]) {
          return;
        }
        lastRange = [rng[0], rng[1]];
        var pxPerData = ax._length / (rng[1] - rng[0]);
        var res = LZR.layout(P, rng[0], rng[1], pxPerData);

        var lx = [], ly = [], lt = [], tx = [], ty = [], tt = [];
        var shapes = base.slice();
        for (var i2 = 0; i2 < res.rows.length; i2++) {
          var it = res.rows[i2], g = P.tx[it.i], y = -it.row;
          lx.push(g.start, g.end, null);
          ly.push(y, y, null);
          lt.push(g.hover, g.hover, null);
          if (g.label !== '') {
            tx.push(it.mid); ty.push(y + 0.35); tt.push(g.label);
          }
        }

        if (P.cfg.showExons === false) {
          /* Mirror R/genetrack_ly.R:129-134: one rect per gene, not per exon. */
          for (var g2 = 0; g2 < res.rows.length; g2++) {
            var itg = res.rows[g2], gg = P.tx[itg.i], rr = itg.row;
            shapes.push({
              type: 'rect', fillcolor: P.cfg.geneCol,
              line: {color: P.cfg.exonBorder, width: 1},
              x0: gg.start, x1: gg.end, xref: P.idx.xaxis,
              y0: -rr - 0.15, y1: -rr + 0.15, yref: yref
            });
          }
        } else {
          for (var j = 0; j < P.ex.length; j++) {
            var e = P.ex[j];
            var r = null;
            for (var k = 0; k < res.rows.length; k++) {
              if (res.rows[k].i === e.gene_idx) { r = res.rows[k].row; break; }
            }
            if (r === null) continue;
            shapes.push({
              type: 'rect', fillcolor: P.cfg.exonCol,
              line: {color: P.cfg.exonBorder, width: 0.5},
              x0: e.start, x1: e.end, xref: P.idx.xaxis,
              y0: -r - 0.15, y1: -r + 0.15, yref: yref
            });
          }
        }

        var ann = baseAnn.slice();
        if (res.rows.length === 0) {
          ann.push(LZR.note('No genes in view', P));
        } else if (res.hidden > 0) {
          ann.push(LZR.note(res.hidden + ' genes not shown — zoom in', P));
        }

        busy = true;
        Plotly.restyle(gd, {x: [lx], y: [ly], text: [lt]}, [P.idx.lineTrace])
          .then(function () {
            return Plotly.restyle(gd, {x: [tx], y: [ty], text: [tt]},
                                  [P.idx.labelTrace]);
          })
          .then(function () {
            return Plotly.relayout(gd, {shapes: shapes, annotations: ann});
          })
          .then(function () {
            busy = false;
            /* A relayout arrived while this apply() was in flight; the
             * three-restyle chain takes long enough (150-300ms on a dense
             * locus) that a dragmode="pan" user can easily pan again before
             * it settles. Re-run once more against the now-current range
             * instead of leaving the panel packed for the stale window. */
            if (pending) { pending = false; schedule(); }
          })
          .catch(function (err) { busy = false; fail(err); });
      }

      gd.on('plotly_relayout', function () {
        if (dead) return;
        if (busy) { pending = true; return; }
        schedule();
      });
    } catch (err) {
      fail(err);
    }
  };

  LZR.note = function (text, P) {
    return {
      text: text, showarrow: false,
      xref: P.idx.xaxis + ' domain', yref: P.idx.yaxis + ' domain',
      x: 1, y: 0, xanchor: 'right', yanchor: 'bottom',
      font: {size: P.cfg.fontSizePx * 0.9, color: '#888888'}
    };
  };

  global.LZR = LZR;
})(typeof window !== 'undefined' ? window : global);
