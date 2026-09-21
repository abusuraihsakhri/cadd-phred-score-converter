"use strict";
require("../app.js");
const { convertPhred, rankContext } = globalThis.CaddMath;
const result = convertPhred(20, 8_600_000_000);
if (Math.abs(result.topPercent - 1) > 1e-12) throw new Error("topPercent mismatch");
if (Math.abs(result.percentile - 99) > 1e-12) throw new Error("percentile mismatch");
if (result.rankFromTop !== 86_000_000) throw new Error("rank mismatch");
if (!rankContext(30).includes("0.1%")) throw new Error("context mismatch");
let failed = false;
try { convertPhred(-1); } catch { failed = true; }
if (!failed) throw new Error("negative PHRED should fail");
console.log("web smoke: PASS");
