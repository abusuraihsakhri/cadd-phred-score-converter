(() => {
  "use strict";

  function assertFiniteNonNegative(value, label) {
    if (!Number.isFinite(value) || value < 0) throw new Error(`${label} must be a finite value ≥ 0.`);
  }

  function convertPhred(phred, referenceSize = 8_600_000_000) {
    assertFiniteNonNegative(phred, "PHRED score");
    if (!Number.isSafeInteger(referenceSize) || referenceSize <= 0) {
      throw new Error("Reference variant count must be a positive integer.");
    }
    const tailFraction = 10 ** (-phred / 10);
    const topPercent = tailFraction * 100;
    const percentile = 100 * (1 - tailFraction);
    const rankFromTop = Math.max(1, Math.min(referenceSize, Math.ceil(referenceSize * tailFraction)));
    return { phred, referenceSize, tailFraction, topPercent, percentile, rankFromTop };
  }

  function rankContext(phred) {
    if (phred >= 40) return "approximately within the top 0.01% of possible substitutions";
    if (phred >= 30) return "approximately within the top 0.1% of possible substitutions";
    if (phred >= 20) return "approximately within the top 1% of possible substitutions";
    if (phred >= 10) return "approximately within the top 10% of possible substitutions";
    return "below the top 10% of possible substitutions";
  }

  function formatNumber(value, maxDigits = 6) {
    if (value !== 0 && Math.abs(value) < 0.0001) return value.toExponential(3);
    return new Intl.NumberFormat(undefined, { maximumFractionDigits: maxDigits }).format(value);
  }

  globalThis.CaddMath = { convertPhred, rankContext };

  if (typeof document === "undefined") return;

  const form = document.getElementById("converterForm");
  const phredInput = document.getElementById("phred");
  const referenceInput = document.getElementById("referenceSize");
  const error = document.getElementById("error");
  const themeToggle = document.getElementById("themeToggle");

  function render() {
    try {
      const phred = Number(phredInput.value);
      const referenceSize = Number(referenceInput.value);
      const result = convertPhred(phred, referenceSize);
      document.getElementById("scoreBadge").textContent = `PHRED ${formatNumber(result.phred, 3)}`;
      document.getElementById("tailFraction").textContent = formatNumber(result.tailFraction, 8);
      document.getElementById("topPercent").textContent = `${formatNumber(result.topPercent, 6)}%`;
      document.getElementById("percentile").textContent = `${formatNumber(result.percentile, 6)}%`;
      document.getElementById("rank").textContent = new Intl.NumberFormat().format(result.rankFromTop);
      document.getElementById("context").textContent = `CADD PHRED ${formatNumber(result.phred, 3)} is ${rankContext(result.phred)}.`;
      error.hidden = true;
    } catch (err) {
      error.textContent = err instanceof Error ? err.message : "Invalid input.";
      error.hidden = false;
    }
  }

  form.addEventListener("submit", (event) => { event.preventDefault(); render(); });

  const storedTheme = localStorage.getItem("theme");
  const initialTheme = storedTheme || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
  document.documentElement.dataset.theme = initialTheme;
  themeToggle.setAttribute("aria-pressed", String(initialTheme === "dark"));
  themeToggle.addEventListener("click", () => {
    const next = document.documentElement.dataset.theme === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    localStorage.setItem("theme", next);
    themeToggle.setAttribute("aria-pressed", String(next === "dark"));
  });

  render();
})();
