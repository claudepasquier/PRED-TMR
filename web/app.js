/**
 * PRED-TMR web portal (phase B).
 *
 * Loads the Emscripten module (createPredTMR, exposed by predtmr.js), then
 * parses a pasted sequence (raw or FASTA/SwissProt/PDB via the shared C++
 * parser) and renders the predicted transmembrane segments.
 *
 * No protein logic lives here: everything is delegated to the validated C++
 * kernel running inside WebAssembly.
 */

"use strict";

/** Residues rendered per line in the annotated sequence view. */
const RESIDUES_PER_LINE = 60;

/** @typedef {{start: number, end: number, score: number}} Segment */
/** @typedef {{name: string, sequence: string, length: number, format: string}} ParsedSequence */
/** @typedef {{segments: Segment[], perResidue: number[], length: number}} Prediction */
/** @typedef {{predict: (sequence: string) => Prediction, parseSequence: (text: string, format: string) => ParsedSequence}} PredTMRModule */
/** @typedef {(options?: Record<string, unknown>) => Promise<PredTMRModule>} PredTMRFactory */

/**
 * Built-in example inputs.
 * @type {Record<string, {label: string, text: string}>}
 */
const EXAMPLES = {
    bacteriorhodopsin: {
        label: "Bacteriorhodopsin",
        text: [
            ">sp|P02945|BACR_HALSA Bacteriorhodopsin (7 transmembrane helices)",
            "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD",
        ].join("\n"),
    },
    "synthetic-tm": {
        label: "Synthetic TM",
        text: [">synthetic single-pass", "MNSVPLVLAAVALLALAAALAAAAPGGGAAA"].join("\n"),
    },
    "no-tm": {
        label: "Soluble protein",
        text: [
            ">sp|P0CG48|UBC_HUMAN Polyubiquitin-C (soluble, no transmembrane segment)",
            "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
        ].join("\n"),
    },
};

/** @type {PredTMRModule | null} */
let module = null;

/** @returns {Record<string, HTMLElement>} */
function elements() {
    const ids = [
        "sequence-input",
        "predict-btn",
        "clear-btn",
        "status",
        "results",
        "empty-state",
        "summary",
        "sequence-view",
        "segments-table",
    ];
    /** @type {Record<string, HTMLElement>} */
    const found = {};
    for (const id of ids) {
        const el = document.getElementById(id);
        if (!el) {
            throw new Error(`Missing element #${id}`);
        }
        found[id] = el;
    }
    return found;
}

/**
 * Sets the status line text and style.
 * @param {HTMLElement} statusEl
 * @param {string} message
 * @param {"idle"|"loading"|"ok"|"error"} tone
 * @returns {void}
 */
function setStatus(statusEl, message, tone) {
    statusEl.textContent = message;
    statusEl.className = `status ${tone === "idle" ? "" : tone}`.trim();
}

/**
 * Ensures the trimmed input is non-empty.
 * @param {string} value
 * @returns {string}
 */
function requireNonEmpty(value) {
    const trimmed = value.trim();
    if (trimmed.length === 0) {
        throw new Error("Please paste a protein sequence first.");
    }
    return trimmed;
}

/**
 * Parses pasted text into a clean residue string using the C++ parser.
 * @param {string} text
 * @returns {ParsedSequence}
 */
function parseInput(text) {
    if (!module) {
        throw new Error("Prediction engine is not ready yet.");
    }
    const parsed = module.parseSequence(requireNonEmpty(text), "auto");
    if (!parsed || typeof parsed.sequence !== "string" || parsed.sequence.length === 0) {
        throw new Error("No residues could be parsed from the input.");
    }
    return parsed;
}

/**
 * Formats a numeric score with fixed precision.
 * @param {number} value
 * @returns {string}
 */
function formatScore(value) {
    return Number.isFinite(value) ? value.toFixed(2) : "—";
}

/**
 * Renders the summary metrics block.
 * @param {HTMLElement} target
 * @param {ParsedSequence} parsed
 * @param {Prediction} prediction
 * @returns {void}
 */
function renderSummary(target, parsed, prediction) {
    const covered = prediction.perResidue.reduce((sum, flag) => sum + (flag ? 1 : 0), 0);
    const name = parsed.name && parsed.name.length > 0 ? parsed.name : "sequence";
    const metrics = [
        ["Name", name],
        ["Format", parsed.format],
        ["Length", String(parsed.length)],
        ["Segments", String(prediction.segments.length)],
        ["Residues in TM", String(covered)],
    ];

    target.replaceChildren();
    for (const [label, value] of metrics) {
        const span = document.createElement("span");
        span.className = "metric";
        span.textContent = `${label}: `;
        const strong = document.createElement("strong");
        strong.textContent = value;
        span.appendChild(strong);
        target.appendChild(span);
    }
}

/**
 * Renders the annotated sequence, highlighting transmembrane residues.
 * @param {HTMLElement} target
 * @param {string} sequence
 * @param {number[]} perResidue
 * @returns {void}
 */
function renderSequence(target, sequence, perResidue) {
    target.replaceChildren();

    for (let offset = 0; offset < sequence.length; offset += RESIDUES_PER_LINE) {
        const line = document.createElement("div");
        line.className = "seq-line";

        const position = document.createElement("span");
        position.className = "seq-pos";
        position.textContent = String(offset + 1);
        line.appendChild(position);

        const chunk = document.createElement("span");
        chunk.className = "seq-chunk";

        const end = Math.min(offset + RESIDUES_PER_LINE, sequence.length);
        for (let index = offset; index < end; index += 1) {
            const residue = document.createElement("span");
            residue.className = perResidue[index] ? "residue tm" : "residue";
            residue.textContent = sequence[index];
            residue.title = `${index + 1}: ${sequence[index]}${perResidue[index] ? " (TM)" : ""}`;
            chunk.appendChild(residue);
        }

        line.appendChild(chunk);
        target.appendChild(line);
    }
}

/**
 * Renders the segments table body.
 * @param {HTMLTableElement} table
 * @param {Segment[]} segments
 * @returns {void}
 */
function renderSegments(table, segments) {
    const body = table.querySelector("tbody");
    if (!body) {
        return;
    }
    body.replaceChildren();

    if (segments.length === 0) {
        const row = document.createElement("tr");
        const cell = document.createElement("td");
        cell.colSpan = 5;
        cell.className = "no-segments";
        cell.textContent = "No transmembrane segment predicted.";
        row.appendChild(cell);
        body.appendChild(row);
        return;
    }

    segments.forEach((segment, index) => {
        const row = document.createElement("tr");
        const values = [
            String(index + 1),
            String(segment.start),
            String(segment.end),
            String(segment.end - segment.start + 1),
            formatScore(segment.score),
        ];
        for (const value of values) {
            const cell = document.createElement("td");
            cell.textContent = value;
            row.appendChild(cell);
        }
        body.appendChild(row);
    });
}

/**
 * Full prediction workflow for the current input.
 * @returns {void}
 */
function runPrediction() {
    const el = elements();
    setStatus(el.status, "Predicting…", "loading");

    try {
        const parsed = parseInput(/** @type {HTMLTextAreaElement} */(el["sequence-input"]).value);
        if (!module) {
            throw new Error("Prediction engine is not ready yet.");
        }
        const prediction = module.predict(parsed.sequence);

        renderSummary(el.summary, parsed, prediction);
        renderSequence(el["sequence-view"], parsed.sequence, prediction.perResidue);
        renderSegments(/** @type {HTMLTableElement} */(el["segments-table"]), prediction.segments);

        el["empty-state"].hidden = true;
        el.results.hidden = false;
        setStatus(
            el.status,
            `Done — ${prediction.segments.length} segment(s) found in ${parsed.length} residues.`,
            "ok",
        );
    } catch (error) {
        setStatus(el.status, `Error: ${error instanceof Error ? error.message : String(error)}`, "error");
    }
}

/**
 * Wires up all DOM event handlers.
 * @param {PredTMRModule} loadedModule
 * @returns {void}
 */
function bindUi(loadedModule) {
    module = loadedModule;
    const el = elements();
    const input = /** @type {HTMLTextAreaElement} */ (el["sequence-input"]);
    const predictBtn = /** @type {HTMLButtonElement} */ (el["predict-btn"]);
    const clearBtn = el["clear-btn"];

    predictBtn.disabled = false;
    setStatus(el.status, "Ready. Paste a sequence and press Predict.", "ok");

    predictBtn.addEventListener("click", runPrediction);

    clearBtn.addEventListener("click", () => {
        input.value = "";
        el.results.hidden = true;
        el["empty-state"].hidden = false;
        setStatus(el.status, "Cleared. Paste a sequence and press Predict.", "idle");
        input.focus();
    });

    document.querySelectorAll(".example-btn").forEach((button) => {
        button.addEventListener("click", () => {
            const key = button.getAttribute("data-example");
            const example = key ? EXAMPLES[key] : undefined;
            if (example) {
                input.value = example.text;
                setStatus(el.status, `Loaded example: ${example.label}.`, "idle");
            }
        });
    });
}

/**
 * Boots the portal: loads the WASM module and binds the UI.
 * @returns {Promise<void>}
 */
async function boot() {
    const el = elements();
    setStatus(el.status, "Loading prediction engine…", "loading");

    const factory = /** @type {PredTMRFactory | undefined} */ (window.createPredTMR);
    if (typeof factory !== "function") {
        setStatus(
            el.status,
            "Error: predtmr.js not found. Build the WebAssembly module first (scripts/build-wasm.sh).",
            "error",
        );
        return;
    }

    try {
        const loaded = await factory();
        bindUi(loaded);
    } catch (error) {
        setStatus(
            el.status,
            `Error: failed to load the WebAssembly module — ${error instanceof Error ? error.message : String(error)}`,
            "error",
        );
    }
}

document.addEventListener("DOMContentLoaded", () => {
    void boot();
});
