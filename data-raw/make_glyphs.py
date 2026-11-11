"""Trace the 20 amino-acid letters (+ '-') of a bold monospace font into flattened
outline polygons: inst/extdata/glyphs_dejavu_sans_mono_bold.tsv.

Columns: letter, contour (1..), point order, x, y (font units / unitsPerEm,
baseline at y = 0). Holes are separate contours; draw with the even-odd rule
(ggplot2 geom_polygon `subgroup`).

DejaVu Sans Mono Bold: Bitstream Vera Fonts Copyright (c) 2003 Bitstream, Inc.;
DejaVu changes are in the public domain. See inst/COPYRIGHTS.

Usage: python data-raw/make_glyphs.py /usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf inst/extdata/glyphs_dejavu_sans_mono_bold.tsv
"""
import sys
from fontTools.ttLib import TTFont
from fontTools.pens.basePen import BasePen

N = 8  # points per curve segment

class FlatPen(BasePen):
    def __init__(self, gs):
        super().__init__(gs)
        self.contours = []
        self.cur = None
    def _moveTo(self, p):
        self.cur = [p]
    def _lineTo(self, p):
        self.cur.append(p)
    def _curveToOne(self, p1, p2, p3):
        p0 = self.cur[-1]
        for i in range(1, N + 1):
            t = i / N
            x = (1-t)**3*p0[0] + 3*(1-t)**2*t*p1[0] + 3*(1-t)*t**2*p2[0] + t**3*p3[0]
            y = (1-t)**3*p0[1] + 3*(1-t)**2*t*p1[1] + 3*(1-t)*t**2*p2[1] + t**3*p3[1]
            self.cur.append((x, y))
    def _qCurveToOne(self, p1, p2):
        p0 = self.cur[-1]
        for i in range(1, N + 1):
            t = i / N
            x = (1-t)**2*p0[0] + 2*(1-t)*t*p1[0] + t**2*p2[0]
            y = (1-t)**2*p0[1] + 2*(1-t)*t*p1[1] + t**2*p2[1]
            self.cur.append((x, y))
    def _closePath(self):
        if self.cur:
            if self.cur[0] == self.cur[-1]:
                self.cur.pop()
            self.contours.append(self.cur)
        self.cur = None
    _endPath = _closePath

font = TTFont(sys.argv[1])
gs = font.getGlyphSet()
cmap = font.getBestCmap()
upem = font["head"].unitsPerEm
with open(sys.argv[2], "w") as fh:
    fh.write("letter\tcontour\torder\tx\ty\n")
    for ch in "ACDEFGHIKLMNPQRSTVWY-":
        pen = FlatPen(gs)
        gs[cmap[ord(ch)]].draw(pen)
        for ci, c in enumerate(pen.contours, 1):
            for pi, (x, y) in enumerate(c, 1):
                fh.write(f"{ch}\t{ci}\t{pi}\t{x/upem:.4f}\t{y/upem:.4f}\n")
print("ok")
