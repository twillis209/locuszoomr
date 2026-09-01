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

      /* Edge-label clamping, mirroring mapRow() in R/genetracks.R. A gene
       * straddling the viewport edge would otherwise anchor its label at its
       * own midpoint, which can sit off-screen, so the bar renders but the
       * name does not. Pull the anchor just inside the boundary instead --
       * but only when the gene is wide enough to hold the whole label, else
       * the text would float beyond the gene it belongs to.
       *
       * The conditions and the sequential left-then-right order match R (a
       * gene wider than the view can satisfy both, and the right-hand test
       * sees the already-clamped anchor), but the boundary does NOT. R
       * clamps to xlim widened by 4%, which works on a graphics device
       * because text may render into the plot margin. Plotly clips hard at
       * the axis range, so an anchor placed at `x0 - 0.04 * span + halfw`
       * is still off-screen whenever halfw is smaller than that overshoot,
       * and only the tail of the label (the strand arrow) shows. Clamp to
       * the visible window itself so the whole label lands inside it.
       *
       * Unlabelled genes are skipped: halfw is 0 for them, which would
       * satisfy the conditions trivially and shift the packing footprint
       * for no visible gain. */
      if (halfw > 0) {
        var gwFull = halfw * 2;
        if ((mid - halfw) < x0 && (x0 + gwFull) < g.end) mid = x0 + halfw;
        if ((mid + halfw) > x1 && (x1 - gwFull) > g.start) mid = x1 - halfw;
      }
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

      /* Shapes that are not ours, captured once; every update rebuilds the
       * array as base.concat(exons). `== null`, not `||`: 0 is a valid
       * shape index, and jsonlite may deliver shapeIdx as a bare number. */
      var idxs = [].concat(P.idx.shapeIdx == null ? [] : P.idx.shapeIdx);
      var all = (gd.layout.shapes || []);
      var base = [];
      for (var i = 0; i < all.length; i++) {
        if (idxs.indexOf(i) === -1) base.push(all[i]);
      }
      var baseAnn = (gd.layout.annotations || []).slice();

      var busy = false, timer = null, pending = false, lastKey = null;

      /* Recognises our own {shapes, annotations} relayout, which would
       * otherwise re-enter the listener and re-trigger apply() forever.
       * This is the structural loop-breaker; the no-op guard in apply() is
       * a cosmetic optimisation and not a safe substitute. */
      function isSelfUpdate(upd) {
        if (!upd) return false;
        var keys = Object.keys(upd);
        if (keys.length === 0) return false;
        for (var k = 0; k < keys.length; k++) {
          if (keys[k] !== 'shapes' && keys[k] !== 'annotations') return false;
        }
        return true;
      }

      function schedule() {
        if (timer) global.clearTimeout(timer);
        timer = global.setTimeout(function () {
          try { apply(); } catch (err) { fail(err); }
        }, 100);
      }

      function apply() {
        var ax = gd._fullLayout[xkey];
        var rng = ax.range;
        var len = ax._length;
        /* Hidden-container guard: a display:none ancestor (tab switch)
         * can emit 'plotly_relayout' while the axis has zero or unset
         * pixel length. Neither throws, so catch() never sees it: len 0
         * gives every gene an infinite footprint (one gene per row), len
         * undefined gives pxPerData NaN (all genes on row 1). Bail before
         * `busy` is set so it cannot be left stuck. */
        if (!(len > 0) || !isFinite(rng[0]) || !isFinite(rng[1])) return;
        /* No-op guard for dragmode toggles and legend clicks. `len` is in
         * the key because a window resize changes `_length` but not the
         * range, and pxPerData depends on it. */
        if (lastKey !== null &&
            rng[0] === lastKey[0] && rng[1] === lastKey[1] && len === lastKey[2]) {
          return;
        }
        lastKey = [rng[0], rng[1], len];
        var pxPerData = len / (rng[1] - rng[0]);
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
        /* One Plotly.update(), not chained restyle/relayout calls: each call
         * repaints the whole figure, so one pass means one redraw. The
         * trace-update arrays are positional: 0 = lineTrace, 1 = labelTrace. */
        Plotly.update(gd,
                      {x: [lx, tx], y: [ly, ty], text: [lt, tt]},
                      {shapes: shapes, annotations: ann},
                      [P.idx.lineTrace, P.idx.labelTrace])
          .then(function () {
            busy = false;
            /* A genuine relayout can arrive while apply() was in flight;
             * re-run against the now-current range rather than leaving the
             * panel packed for a stale window. No loop risk: Plotly.update()
             * emits no 'plotly_relayout', so `pending` is only ever set by a
             * real user event. */
            if (pending) { pending = false; schedule(); }
          })
          .catch(function (err) { busy = false; fail(err); });
      }

      gd.on('plotly_relayout', function (upd) {
        if (dead) return;
        if (isSelfUpdate(upd)) return;
        if (busy) { pending = true; return; }
        schedule();
      });

      /* Re-pack once on the initial render: R packs against a device-based
       * width heuristic, not real canvas text metrics, so without this the
       * widget arrives packed one way and silently reorganises on the first
       * relayout. Inside a hidden container there is no laid-out axis yet
       * and requestAnimationFrame does not fire, so bound the retries; the
       * relayout listener corrects the packing once the container shows. */
      var packTries = 0;
      var initialPack = function () {
        if (dead) return;
        var ax0 = gd._fullLayout && gd._fullLayout[xkey];
        if (!ax0 || !(ax0._length > 0)) {
          if (packTries++ < 60 && global.requestAnimationFrame) {
            global.requestAnimationFrame(initialPack);
          }
          return;
        }
        try { apply(); } catch (err) { fail(err); }
      };
      initialPack();
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
