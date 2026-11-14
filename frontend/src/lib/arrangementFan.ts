/**
 * The fan under the focal gene — one curve from the gene down to each arrangement it offers.
 *
 * ⭐ **Pure arithmetic over rectangles**, like `popoverPlacement`: jsdom computes no layout, so the
 * caller measures and this decides. `app.js::drawArrFan` is the reference, number for number.
 *
 * ⚠ **The gene moves under a panel that does not.** The track scrolls sideways; the switcher below it
 * does not. So the fan is recomputed from LIVE rectangles on every scroll and resize, and when the
 * focal gene has scrolled out of the scroller the fan draws nothing at all — a line that points at
 * nothing is worse than no line.
 */

import type { Rect } from "./popoverPlacement";

export interface FanRect extends Rect {
  readonly right: number;
}

export interface FanCurve {
  /** SVG path data, in the wrapper's own coordinates. */
  readonly d: string;
  /** The drawn arrangement's curve — `.fan-line.on`, in the accent. */
  readonly isOn: boolean;
}

export interface Fan {
  readonly width: number;
  readonly height: number;
  readonly curves: readonly FanCurve[];
  /** Where the curves meet, under the focal gene — the apex the locus's own count is written at. */
  readonly apexX: number | null;
}

/** What a fan with nothing to join looks like: zero height, so it takes no room. */
export const NO_FAN: Fan = { width: 0, height: 0, curves: [], apexX: null };

/** The published page's floor, so a fan to options sitting right under the gene is still visible. */
const MIN_HEIGHT_PX = 8;
/** Inside this many pixels of the scroller's edge the gene counts as scrolled away. */
const EDGE_PX = 4;

/**
 * The curves from the focal block to each option.
 *
 * ⚠ A cubic that LEAVES the gene vertically and ARRIVES vertically: a straight line to an option far
 * to one side reads as a diagonal crossing the track, which is exactly the confusion the fan exists to
 * remove (`app.js:1712`). The control points sit at 55 % of the drop, as published.
 */
export function arrangementFan(
  focal: FanRect,
  wrapper: FanRect,
  options: readonly { readonly rect: FanRect; readonly isOn: boolean }[],
  scroller: FanRect | null,
  apexDrop = 0,
): Fan {
  if (options.length === 0) return NO_FAN;
  if (scroller !== null && (focal.right < scroller.left + EDGE_PX || focal.left > scroller.right - EDGE_PX)) {
    return NO_FAN;
  }
  const x0 = focal.left + focal.width / 2 - wrapper.left;
  // ⭐ `apexDrop` leaves room above the apex for the locus's own "n genes · x % of genomes" (David,
  // 2026-11-14): the whole above the parts the options divide it into.
  const y0 = apexDrop;
  const height = Math.max(MIN_HEIGHT_PX, ...options.map((option) => option.rect.top - wrapper.top));
  const curves = options.map(({ rect, isOn }) => {
    const x1 = rect.left + rect.width / 2 - wrapper.left;
    const y1 = rect.top - wrapper.top;
    const mid = y0 + (y1 - y0) * 0.55;
    return { d: `M${x0} ${y0} C${x0} ${mid} ${x1} ${mid} ${x1} ${y1}`, isOn };
  });
  return { width: wrapper.width, height, curves, apexX: x0 };
}
