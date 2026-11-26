/**
 * Drive the built site in a real browser and report what actually happens.
 *
 *   node tools/browser-check.mjs <baseUrl> <screenshotDir>
 *
 * Checks the things a unit test cannot: that the page renders, that the model and the
 * static indexes actually load over HTTP, that a design run produces a populated
 * results table, and that nothing throws or 404s on the way. Also captures light,
 * dark and narrow-viewport screenshots, and flags horizontal overflow at phone width.
 *
 * Exits non-zero if any check fails.
 */

import { chromium } from "playwright";

const baseUrl = process.argv[2] ?? "http://127.0.0.1:8099/";
const shotDir = process.argv[3] ?? "/tmp/splashrna-shots";

const failures = [];
const note = (ok, label, detail = "") => {
  console.log(`${ok ? "PASS" : "FAIL"}  ${label}${detail ? ` — ${detail}` : ""}`);
  if (!ok) failures.push(label);
};
const skip = (label, why) => console.log(`SKIP  ${label} — ${why}`);

/** The large indexes are optional deployments; CI builds without them. */
async function deployed(path) {
  try {
    return (await fetch(new URL(path, baseUrl))).ok;
  } catch {
    return false;
  }
}

/**
 * The page probes for the optional indexes by fetching their manifests, so a 404 on
 * those two paths is the documented way of discovering that a deployment does not have
 * them. Everything else counts as a problem.
 */
const EXPECTED_404 = [/\/offtarget\/manifest\.json$/, /\/genes\/manifest\.json$/];
const isCapabilityProbe = (url) => EXPECTED_404.some((re) => re.test(url));

/** Chromium's console message for a failed subresource carries no URL of its own. */
const GENERIC_LOAD_FAILURE = /^Failed to load resource/;

/** Attach console/network listeners that record anything alarming. */
function watch(page) {
  const consoleErrors = [];
  const pageErrors = [];
  const failedRequests = [];
  const probe404s = [];

  page.on("console", (msg) => {
    if (msg.type() === "error") consoleErrors.push(msg.text());
  });
  page.on("pageerror", (err) => pageErrors.push(String(err)));
  page.on("requestfailed", (req) => {
    if (!isCapabilityProbe(req.url())) {
      failedRequests.push(`${req.url()} (${req.failure()?.errorText})`);
    }
  });
  page.on("response", (res) => {
    if (res.status() < 400) return;
    if (isCapabilityProbe(res.url())) probe404s.push(res.url());
    else failedRequests.push(`${res.url()} → ${res.status()}`);
  });

  return {
    get failedRequests() {
      return failedRequests;
    },
    get pageErrors() {
      return pageErrors;
    },
    /** Console errors, less the generic ones explained by a capability probe. */
    get console() {
      let allowance = probe404s.length;
      return consoleErrors.filter((text) => {
        if (allowance > 0 && GENERIC_LOAD_FAILURE.test(text)) {
          allowance--;
          return false;
        }
        return true;
      });
    },
    get probes() {
      return probe404s;
    },
  };
}

/** The academic-use gate blocks designing until it is satisfied. */
async function satisfyGate(page) {
  await page.check("#academic");
  await page.fill("#email", "someone@example.edu");
}

async function runDesign(page, fasta, { timeout = 120000 } = {}) {
  await satisfyGate(page);
  await page.fill("#fasta", fasta);
  await page.click("#run");
  await page.waitForFunction(
    () => {
      const status = document.getElementById("status");
      return status && (/^Done/.test(status.textContent) || status.classList.contains("error"));
    },
    { timeout },
  );
  return page.textContent("#status");
}

const browser = await chromium.launch();

try {
  // --- 1. Core path: a pasted sequence, no network dependencies -------------------
  {
    const context = await browser.newContext({ viewport: { width: 1280, height: 900 } });
    const page = await context.newPage();
    const problems = watch(page);

    await page.goto(baseUrl, { waitUntil: "networkidle" });
    note((await page.title()).includes("SplashRNA"), "page loads", await page.title());
    await page.screenshot({ path: `${shotDir}/01-initial.png`, fullPage: true });

    // A 600 nt sequence, enough to exercise scoring properly.
    let seq = "";
    let seed = 12345;
    for (let i = 0; i < 600; i++) {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      seq += "ACGT"[(seed >>> 16) & 3];
    }

    const t0 = Date.now();
    const status = await runDesign(page, `> pasted\n${seq}`);
    const elapsed = Date.now() - t0;
    note(/^Done/.test(status), "sequence design completes", status);

    const rows = await page.locator(".result tbody tr").count();
    note(rows === 5, "results table has 5 rows", `got ${rows}`);
    note(elapsed < 30000, "completes in reasonable time", `${(elapsed / 1000).toFixed(1)}s`);

    const guide = await page
      .locator(".result tbody tr:first-child td.seq")
      .first()
      .textContent();
    note(/^[ACGT]{22}$/.test(guide ?? ""), "first guide is a 22-mer", guide ?? "");

    const score = await page
      .locator(".result tbody tr:first-child td.score")
      .first()
      .textContent();
    note(/^-?\d+\.\d{3}$/.test(score ?? ""), "score is rendered", score ?? "");

    const csvDisabled = await page.locator("#csv").isDisabled();
    note(!csvDisabled, "CSV button enabled after a run");

    // The gate must actually block, or it is decoration.
    await page.uncheck("#academic");
    await page.fill("#email", "");
    await page.click("#run");
    await page.waitForFunction(
      () => document.getElementById("status").classList.contains("error"),
      { timeout: 10000 },
    ).catch(() => {});
    const blockedMsg = await page.textContent("#status");
    note(
      /academic user/i.test(blockedMsg ?? ""),
      "academic-use gate blocks designing",
      (blockedMsg ?? "").slice(0, 60),
    );
    await satisfyGate(page);

    await page.screenshot({ path: `${shotDir}/02-results.png`, fullPage: true });

    note(problems.pageErrors.length === 0, "no uncaught JS errors", problems.pageErrors.join(" | "));
    note(
      problems.failedRequests.length === 0,
      "no failed requests",
      problems.failedRequests.slice(0, 3).join(" | "),
    );
    note(problems.console.length === 0, "no console errors", problems.console.slice(0, 3).join(" | "));
    await context.close();
  }

  // --- 2. Entrez gene ID, served by the local static gene index -------------------
  const hasGeneIndex = await deployed("genes/manifest.json");
  const hasOffTarget = await deployed("offtarget/manifest.json");

  if (!hasGeneIndex) {
    skip("Entrez gene index checks", "public/genes is not built (see tools/build_genes.py)");
  } else {
    const context = await browser.newContext({ viewport: { width: 1280, height: 900 } });
    const page = await context.newPage();
    const problems = watch(page);
    await page.goto(baseUrl, { waitUntil: "networkidle" });

    const status = await runDesign(page, "> TP53\n7157");
    note(/^Done/.test(status), "Entrez gene design completes", status);

    const rows = await page.locator(".result tbody tr").count();
    note(rows === 5, "Entrez results populated", `${rows} rows`);

    const offTarget = await page
      .locator(".result tbody tr:first-child td:nth-child(6)")
      .first()
      .textContent();
    if (hasOffTarget) {
      note(
        offTarget === "unique" || /human/.test(offTarget ?? ""),
        "off-target column populated",
        offTarget ?? "",
      );
    } else {
      note(offTarget === "—", "off-target column shows 'not consulted'", offTarget ?? "");
    }

    // NCBI is not consulted when the local gene index answers.
    const hitNcbi = problems.failedRequests.some((r) => r.includes("eutils"));
    note(!hitNcbi, "no NCBI fallback needed");

    await page.screenshot({ path: `${shotDir}/03-entrez.png`, fullPage: true });
    note(problems.pageErrors.length === 0, "no JS errors on Entrez path", problems.pageErrors.join(" | "));
    await context.close();
  }

  // --- 3. Dark mode ----------------------------------------------------------------
  {
    const context = await browser.newContext({
      viewport: { width: 1280, height: 900 },
      colorScheme: "dark",
    });
    const page = await context.newPage();
    await page.goto(baseUrl, { waitUntil: "networkidle" });
    const bg = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
    const isDark = (() => {
      const m = /rgb\((\d+), (\d+), (\d+)\)/.exec(bg);
      return m && Number(m[1]) + Number(m[2]) + Number(m[3]) < 250;
    })();
    note(!!isDark, "dark mode applies a dark background", bg);
    await page.screenshot({ path: `${shotDir}/04-dark.png`, fullPage: true });
    await context.close();
  }

  // --- 4. Phone width --------------------------------------------------------------
  {
    const context = await browser.newContext({ viewport: { width: 390, height: 844 } });
    const page = await context.newPage();
    await page.goto(baseUrl, { waitUntil: "networkidle" });
    await runDesign(page, "> short\n" + "ACGT".repeat(40));

    const overflow = await page.evaluate(
      () => document.documentElement.scrollWidth - document.documentElement.clientWidth,
    );
    note(overflow <= 1, "no horizontal overflow at 390px", `${overflow}px`);
    await page.screenshot({ path: `${shotDir}/05-mobile.png`, fullPage: true });
    await context.close();
  }
} finally {
  await browser.close();
}

console.log(
  failures.length === 0
    ? "\nAll browser checks passed."
    : `\n${failures.length} check(s) FAILED: ${failures.join(", ")}`,
);
process.exit(failures.length === 0 ? 0 : 1);
