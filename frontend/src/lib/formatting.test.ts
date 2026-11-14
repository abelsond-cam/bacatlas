import { describe, expect, it } from "vitest";

import { sharePercent } from "./formatting";

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
