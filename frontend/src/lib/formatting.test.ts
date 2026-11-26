import { describe, expect, it } from "vitest";

import { percentKeepingSide, sharePercent } from "./formatting";

describe("⛔ percentKeepingSide — never on the wrong side of the line its sentence draws", () => {
  it("whole percent where whole percent is honest", () => {
    expect(percentKeepingSide(0.82)).toBe("82%");
    expect(percentKeepingSide(0.96, 0.95)).toBe("96%");
    expect(percentKeepingSide(0)).toBe("0%");
    expect(percentKeepingSide(1)).toBe("100%");
  });

  it("⛔ keeps decimals until the printed value sits on the TRUE side of a strict 95 % line", () => {
    expect(percentKeepingSide(0.9496, 0.95)).toBe("94.96%"); // not shared — "95%" would read as a tie
    expect(percentKeepingSide(0.9504, 0.95)).toBe("95.04%"); // shared — "95%" would read as not
    expect(percentKeepingSide(0.95, 0.95)).toBe("95%"); // exactly on it: 19 of 20 genes
    // …and without a line to respect, the same values round as usual.
    expect(percentKeepingSide(0.9496)).toBe("95%");
  });

  it("⛔ never prints 100 % short of the whole, nor 0 % above nothing — kp-71k 4008's Wzi share", () => {
    expect(percentKeepingSide(69_825 / 69_884, 0.95)).toBe("99.9%");
    expect(percentKeepingSide(59 / 69_884)).toBe("0.1%"); // 0.084 % — a decimal, never "0%"
    expect(percentKeepingSide(2 / 69_884)).toBe("0.003%");
  });
});

describe("sharePercent", () => {
  it("rounds to whole percent in the middle, as the published page did", () => {
    expect(sharePercent(0.68)).toBe("68%");
    expect(sharePercent(0.2)).toBe("20%");
  });

  it("keeps one decimal below 10 %, so a real listed row never reads 0 %", () => {
    expect(sharePercent(0.004)).toBe("0.4%");
    expect(sharePercent(0.046)).toBe("4.6%");
  });

  it("⛔ never prints 100 % for a share short of the whole — 71,136 of 71,435 genomes is 99.6 %", () => {
    expect(sharePercent(71136 / 71435)).toBe("99.6%");
    // A share so close to 1 that one decimal would still round up reads 99.9 %, never "100.0 %".
    expect(sharePercent(71434 / 71435)).toBe("99.9%");
    expect(sharePercent(0.994)).toBe("99%");
    expect(sharePercent(1)).toBe("100%");
  });
});
