/* ===========================================================
   FTIR Web App — Client-side logic
   =========================================================== */

// Composition oxide keys (must match backend DEFAULTS)
const OXIDE_KEYS = [
    "SiO2", "TiO2", "Al2O3", "FeOT", "MnO",
    "MgO", "CaO", "Na2O", "K2O"
];

/* -----------------------------------------------------------
   1. Load defaults from server on page load
   ----------------------------------------------------------- */
document.addEventListener("DOMContentLoaded", async () => {
    try {
        const resp = await fetch("/api/defaults");
        const d = await resp.json();

        // Composition
        for (const ox of OXIDE_KEYS) {
            const el = document.getElementById("comp_" + ox);
            if (el && d.composition && d.composition[ox] !== undefined) {
                el.value = d.composition[ox];
            }
        }

        // Thickness
        setVal("thickness_um", d.thickness_um);
        setVal("thickness_unc_um", d.thickness_unc_um);
        setVal("olivine_fo", d.olivine_fo);

        // Advanced
        setVal("fe3_ratio", d.fe3_ratio);
        setVal("epsilon_oh", d.epsilon_oh);
        setVal("epsilon_h2omol", d.epsilon_h2omol);
        setVal("epsilon_co3", d.epsilon_co3);

        // CO3 model
        if (d.co3_model) {
            document.getElementById("co3_model").value = d.co3_model;
        }

        // Fringe correction
        if (d.fringe_correction) {
            document.getElementById("fringe_correction").value = d.fringe_correction;
        }

        // Thickness mode
        if (d.thickness_mode) {
            document.getElementById("thickness_mode").value = d.thickness_mode;
            switchThicknessMode(d.thickness_mode);
        }

        // Composition-dependent epsilon checkbox
        const cb = document.getElementById("use_composition_epsilon");
        if (cb) {
            cb.checked = !!d.use_composition_epsilon;
            toggleEpsilonFields();
        }
    } catch (e) {
        console.warn("Could not load defaults:", e);
    }
});

function setVal(id, val) {
    const el = document.getElementById(id);
    if (el && val !== undefined && val !== null) el.value = val;
}

/* -----------------------------------------------------------
   App mode switching (Single / Batch)
   ----------------------------------------------------------- */
function switchAppMode(mode) {
    const singlePanel = document.getElementById("single-panel");
    const batchPanel = document.getElementById("batch-panel");
    const btnSingle = document.getElementById("btn-mode-single");
    const btnBatch = document.getElementById("btn-mode-batch");
    const btnRun = document.getElementById("btn-run");
    const btnBatchRun = document.getElementById("btn-batch-run");
    const resultsSection = document.getElementById("results-section");

    if (mode === "batch") {
        singlePanel.style.display = "none";
        batchPanel.style.display = "block";
        btnSingle.classList.remove("active");
        btnBatch.classList.add("active");
        btnRun.style.display = "none";
        btnBatchRun.style.display = "inline-block";
        resultsSection.classList.remove("visible");
    } else {
        singlePanel.style.display = "block";
        batchPanel.style.display = "none";
        btnSingle.classList.add("active");
        btnBatch.classList.remove("active");
        btnRun.style.display = "inline-block";
        btnBatchRun.style.display = "none";
    }
}

function switchBatchParamMode(mode) {
    document.getElementById("batch-csv-panel").style.display =
        mode === "csv" ? "block" : "none";
}

function downloadParamTemplate() {
    const cols = ["filename","thickness_um","thickness_unc_um","olivine_fo",
                  "SiO2","TiO2","Al2O3","FeOT","MnO","MgO","CaO","Na2O","K2O",
                  "co3_model","fringe_correction","fe3_ratio",
                  "use_composition_epsilon"];
    const example = ["sample1.txt","27.4","0.5","85",
                     "49.2","3.18","12.82","14.41","0.22","5.64","9.37","3.22","0.84",
                     "PCA baseline","auto","0.15","FALSE"];
    const csv = cols.join(",") + "\n" + example.join(",") + "\n";
    const blob = new Blob([csv], { type: "text/csv;charset=utf-8;" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "batch_params_template.csv";
    link.click();
    URL.revokeObjectURL(link.href);
}

async function submitBatch() {
    const btn = document.getElementById("btn-batch-run");
    const status = document.getElementById("status");
    const resultsSection = document.getElementById("results-section");

    const batchFiles = document.getElementById("batch_files").files;
    if (!batchFiles || batchFiles.length === 0) {
        showStatus("Please select spectrum files.", "error");
        return;
    }

    const fd = new FormData();

    // Add all spectrum files
    for (let i = 0; i < batchFiles.length; i++) {
        fd.append("batch_files", batchFiles[i]);
    }

    // Add reflectance files (for Mode C, matched by filename)
    const refFiles = document.getElementById("batch_ref_files").files;
    if (refFiles && refFiles.length > 0) {
        for (let i = 0; i < refFiles.length; i++) {
            fd.append("batch_ref_files", refFiles[i]);
        }
    }

    // Parameter mode
    const paramMode = document.getElementById("batch_param_mode").value;
    fd.append("batch_param_mode", paramMode);

    if (paramMode === "csv") {
        const csvFile = document.getElementById("batch_param_csv").files[0];
        if (csvFile) fd.append("batch_param_csv", csvFile);
    }

    // Always send current form parameters (used as unified or defaults)
    const comp = {};
    for (const ox of OXIDE_KEYS) {
        const el = document.getElementById("comp_" + ox);
        if (el && el.value !== "") comp[ox] = parseFloat(el.value);
    }
    fd.append("composition", JSON.stringify(comp));

    const mode = document.getElementById("thickness_mode").value;
    fd.append("thickness_mode", mode);
    fd.append("thickness_um", document.getElementById("thickness_um").value);
    fd.append("thickness_unc_um", document.getElementById("thickness_unc_um").value || "0");
    fd.append("olivine_fo", document.getElementById("olivine_fo").value);

    const foA = document.getElementById("olivine_fo_a");
    if (foA && foA.value) fd.append("olivine_fo_a", foA.value);

    if (mode === "A") {
        const measurements = [];
        document.querySelectorAll("#fringe-table tbody tr").forEach(tr => {
            const hi = tr.querySelector(".fringe-hi").value;
            const lo = tr.querySelector(".fringe-lo").value;
            const n = tr.querySelector(".fringe-n").value;
            if (hi && lo && n) measurements.push([parseFloat(hi), parseFloat(lo), parseInt(n)]);
        });
        fd.append("thickness_measurements", JSON.stringify(measurements));
    }

    fd.append("co3_model", document.getElementById("co3_model").value);
    fd.append("fringe_correction", document.getElementById("fringe_correction").value);
    fd.append("fe3_ratio", document.getElementById("fe3_ratio").value);
    fd.append("use_composition_epsilon",
        document.getElementById("use_composition_epsilon").checked ? "true" : "false");
    if (!document.getElementById("use_composition_epsilon").checked) {
        fd.append("epsilon_oh", document.getElementById("epsilon_oh").value);
        fd.append("epsilon_h2omol", document.getElementById("epsilon_h2omol").value);
        fd.append("epsilon_co3", document.getElementById("epsilon_co3").value);
    }

    // Send
    btn.disabled = true;
    resultsSection.classList.remove("visible");
    showStatus('<span class="spinner"></span> Processing ' + batchFiles.length + ' files...', "loading");

    try {
        const resp = await fetch("/api/batch", { method: "POST", body: fd });

        if (!resp.ok) {
            const data = await resp.json();
            showStatus("Error: " + (data.error || "Unknown error"), "error");
            btn.disabled = false;
            return;
        }

        // Download ZIP
        const blob = await resp.blob();
        const link = document.createElement("a");
        link.href = URL.createObjectURL(blob);
        link.download = "ftir_batch_results.zip";
        link.click();
        URL.revokeObjectURL(link.href);

        showStatus("Batch complete! " + batchFiles.length + " files processed. ZIP downloaded.", "loading");
        status.style.color = "#27ae60";

    } catch (e) {
        showStatus("Network error: " + e.message, "error");
    } finally {
        btn.disabled = false;
    }
}

/* -----------------------------------------------------------
   Changelog
   ----------------------------------------------------------- */
let changelogLoaded = false;

function toggleChangelog() {
    const panel = document.getElementById("changelog-panel");
    const visible = panel.style.display !== "none";
    panel.style.display = visible ? "none" : "block";
    if (!visible && !changelogLoaded) loadChangelog();
}

async function loadChangelog() {
    try {
        const resp = await fetch("/api/changelog");
        const entries = await resp.json();
        const body = document.getElementById("changelog-body");
        if (!entries.length) {
            body.innerHTML = "<p>No changelog entries.</p>";
            return;
        }
        body.innerHTML = entries.map(e =>
            '<div class="log-entry">' +
            '<span class="log-date">' + e.date + '</span>' +
            '<ul class="log-items">' +
            e.items.map(item => '<li>' + item + '</li>').join('') +
            '</ul></div>'
        ).join('');
        changelogLoaded = true;
    } catch (err) {
        document.getElementById("changelog-body").innerHTML = "Failed to load changelog.";
    }
}

/* -----------------------------------------------------------
   2. Thickness mode switching
   ----------------------------------------------------------- */
document.getElementById("thickness_mode").addEventListener("change", function () {
    switchThicknessMode(this.value);
});

function switchThicknessMode(mode) {
    document.querySelectorAll(".mode-panel").forEach(p => p.classList.remove("active"));
    const panel = document.getElementById("panel-" + mode);
    if (panel) panel.classList.add("active");
}

/* -----------------------------------------------------------
   3. Mode A — dynamic fringe rows
   ----------------------------------------------------------- */
function addFringeRow() {
    const tbody = document.querySelector("#fringe-table tbody");
    const tr = document.createElement("tr");
    tr.innerHTML = `
        <td><input type="number" class="fringe-hi" step="any"></td>
        <td><input type="number" class="fringe-lo" step="any"></td>
        <td><input type="number" class="fringe-n" step="1" min="1"></td>
        <td><button type="button" class="btn-sm btn-remove" onclick="removeFringeRow(this)">Remove</button></td>
    `;
    tbody.appendChild(tr);
}

function removeFringeRow(btn) {
    const tbody = document.querySelector("#fringe-table tbody");
    if (tbody.rows.length > 1) {
        btn.closest("tr").remove();
    }
}

/* -----------------------------------------------------------
   4. Toggle epsilon fields based on checkbox
   ----------------------------------------------------------- */
document.getElementById("use_composition_epsilon").addEventListener("change", toggleEpsilonFields);

function toggleEpsilonFields() {
    const checked = document.getElementById("use_composition_epsilon").checked;
    const fields = document.getElementById("epsilon-fields");
    fields.style.display = checked ? "none" : "block";
}

/* -----------------------------------------------------------
   5. Submit analysis
   ----------------------------------------------------------- */
async function submitAnalysis() {
    const btn = document.getElementById("btn-run");
    const status = document.getElementById("status");
    const resultsSection = document.getElementById("results-section");

    // --- Validate transmission file ---
    const transFile = document.getElementById("transmission_file").files[0];
    if (!transFile) {
        showStatus("Please select a transmission spectrum file.", "error");
        return;
    }

    // --- Build FormData ---
    const fd = new FormData();
    fd.append("transmission_file", transFile);

    const refFile = document.getElementById("reflectance_file").files[0];
    if (refFile) fd.append("reflectance_file", refFile);

    // Composition
    const comp = {};
    for (const ox of OXIDE_KEYS) {
        const el = document.getElementById("comp_" + ox);
        if (el && el.value !== "") {
            comp[ox] = parseFloat(el.value);
        }
    }
    fd.append("composition", JSON.stringify(comp));

    // Thickness mode
    const mode = document.getElementById("thickness_mode").value;
    fd.append("thickness_mode", mode);

    if (mode === "A") {
        const measurements = [];
        document.querySelectorAll("#fringe-table tbody tr").forEach(tr => {
            const hi = tr.querySelector(".fringe-hi").value;
            const lo = tr.querySelector(".fringe-lo").value;
            const n = tr.querySelector(".fringe-n").value;
            if (hi && lo && n) {
                measurements.push([parseFloat(hi), parseFloat(lo), parseInt(n)]);
            }
        });
        fd.append("thickness_measurements", JSON.stringify(measurements));
        const foA = document.getElementById("olivine_fo_a").value;
        if (!foA) {
            showStatus("Mode A requires Olivine Fo#.", "error");
            return;
        }
        fd.append("olivine_fo_a", foA);
    } else if (mode === "B") {
        fd.append("thickness_um", document.getElementById("thickness_um").value);
        fd.append("thickness_unc_um", document.getElementById("thickness_unc_um").value || "0");
    } else if (mode === "C") {
        fd.append("olivine_fo", document.getElementById("olivine_fo").value);
    }

    // CO3 model & fringe correction
    fd.append("co3_model", document.getElementById("co3_model").value);
    fd.append("fringe_correction", document.getElementById("fringe_correction").value);

    // Advanced
    fd.append("fe3_ratio", document.getElementById("fe3_ratio").value);
    fd.append("use_composition_epsilon",
        document.getElementById("use_composition_epsilon").checked ? "true" : "false");

    if (!document.getElementById("use_composition_epsilon").checked) {
        fd.append("epsilon_oh", document.getElementById("epsilon_oh").value);
        fd.append("epsilon_h2omol", document.getElementById("epsilon_h2omol").value);
        fd.append("epsilon_co3", document.getElementById("epsilon_co3").value);
    }

    // --- Send request ---
    btn.disabled = true;
    resultsSection.classList.remove("visible");
    showStatus('<span class="spinner"></span> Running analysis...', "loading");

    try {
        const resp = await fetch("/api/process", { method: "POST", body: fd });
        const data = await resp.json();

        if (!resp.ok || data.error) {
            showStatus("Error: " + (data.error || "Unknown error"), "error");
            btn.disabled = false;
            return;
        }

        status.style.display = "none";
        renderResults(data);
        resultsSection.classList.add("visible");

    } catch (e) {
        showStatus("Network error: " + e.message, "error");
    } finally {
        btn.disabled = false;
    }
}

function showStatus(html, cls) {
    const el = document.getElementById("status");
    el.innerHTML = html;
    el.className = cls;
}

/* -----------------------------------------------------------
   6. Render results
   ----------------------------------------------------------- */
let lastResultData = null;  // store for download

function renderResults(data) {
    lastResultData = data;
    // Physical properties
    fillTable("table-physical", buildPhysicalRows(data));

    // H2O peak
    fillTable("table-h2o", buildH2ORows(data));

    // CO3 fit
    fillTable("table-co3", buildCO3Rows(data));

    // Concentrations
    fillTable("table-conc", buildConcRows(data));

    // Figures
    renderFigures(data.figures || {});

    // Interactive charts
    renderCO3Interactive(data);
    renderH2OInteractive(data);
}

function fillTable(tableId, rows) {
    const tbody = document.querySelector("#" + tableId + " tbody");
    tbody.innerHTML = "";
    for (const [label, value] of rows) {
        const tr = document.createElement("tr");
        tr.innerHTML = "<th>" + label + "</th><td>" + value + "</td>";
        tbody.appendChild(tr);
    }
}

function fmt(val, digits) {
    if (val === null || val === undefined) return "--";
    if (typeof val === "string") return val;
    return Number(val).toFixed(digits !== undefined ? digits : 4);
}

function buildPhysicalRows(data) {
    const rows = [];
    const t = data.thickness || {};
    const tVal = t.thickness_um !== undefined ? t.thickness_um : t.average_um;
    rows.push(["Thickness (um)", fmt(tVal, 2) + " +/- " + fmt(t.stdev_um, 2)]);
    if (t.method) rows.push(["Method", t.method]);
    if (t.n_maxima_used !== undefined) {
        rows.push(["Fringes used", t.n_maxima_used + " maxima, " + t.n_minima_used + " minima"]);
    }

    const d = data.density || {};
    if (d.density_gcc !== undefined) {
        rows.push(["Density (g/cc)", fmt(d.density_gcc, 4)]);
    } else if (d.density_g_per_cc !== undefined) {
        rows.push(["Density (g/cc)", fmt(d.density_g_per_cc, 4)]);
    }

    const ri = data.refractive_index || {};
    if (ri.n !== undefined) {
        rows.push(["Refractive index (n)", fmt(ri.n, 4)]);
    }

    if (data.epsilon) {
        const eps = data.epsilon;
        if (eps.epsilon_oh !== undefined) rows.push(["Epsilon OH", fmt(eps.epsilon_oh, 2)]);
        if (eps.epsilon_h2omol !== undefined) rows.push(["Epsilon H2Omol", fmt(eps.epsilon_h2omol, 2)]);
        if (eps.epsilon_co3 !== undefined) rows.push(["Epsilon CO3", fmt(eps.epsilon_co3, 2)]);
    }

    return rows;
}

function buildH2ORows(data) {
    const rows = [];
    const h = data.h2o_peak || {};
    if (h.peak_absorbance !== undefined) rows.push(["Peak absorbance", fmt(h.peak_absorbance, 5)]);
    if (h.peak_wavenumber !== undefined) rows.push(["Peak wavenumber (cm-1)", fmt(h.peak_wavenumber, 1)]);
    if (h.peak_absorbance_unc !== undefined) rows.push(["Absorbance uncertainty", fmt(h.peak_absorbance_unc, 5)]);
    return rows;
}

function buildCO3Rows(data) {
    const rows = [];
    const c = data.co3_fit || {};
    if (c.model !== undefined) rows.push(["Model", c.model]);
    if (c.co3_absorbance !== undefined) rows.push(["CO3 absorbance", fmt(c.co3_absorbance, 5)]);
    if (c.co3_total_unc !== undefined) rows.push(["CO3 uncertainty", fmt(c.co3_total_unc, 5)]);
    if (c.h2omol_1630_absorbance !== undefined) rows.push(["H2Omol 1630 abs", fmt(c.h2omol_1630_absorbance, 5)]);
    if (c.r_squared !== undefined) rows.push(["R-squared", fmt(c.r_squared, 5)]);
    if (c.quality_flag !== undefined) rows.push(["Quality", c.quality_flag]);
    if (c.shift_1430 !== undefined) rows.push(["Shift 1430 (cm-1)", fmt(c.shift_1430, 2)]);
    if (c.shift_1515 !== undefined) rows.push(["Shift 1515 (cm-1)", fmt(c.shift_1515, 2)]);
    if (c.fringe_detected) {
        rows.push(["Fringe detected", "YES (amp=" + fmt(c.fringe_amplitude, 4) +
            ", period=" + fmt(c.fringe_period, 0) + " cm-1, ratio=" + fmt(c.fringe_ratio, 2) + ")"]);
    }
    if (c.fringe_corrected) {
        rows.push(["Fringe corrected", "YES (period=" + fmt(c.fringe_period_est, 0) + " cm-1)"]);
    }
    return rows;
}

function buildConcRows(data) {
    const rows = [];
    const c = data.concentration || {};
    if (c.total_h2o_wt_pct !== undefined) {
        let s = fmt(c.total_h2o_wt_pct, 4) + " wt%";
        if (c.total_h2o_unc_wt_pct !== undefined) s += " +/- " + fmt(c.total_h2o_unc_wt_pct, 4);
        rows.push(["Total H2O", s]);
    }
    if (c.h2omol_1630_wt_pct !== undefined) {
        rows.push(["H2Omol (1630)", fmt(c.h2omol_1630_wt_pct, 4) + " wt%"]);
    }
    if (c.co2_ppm !== undefined) {
        let s;
        if (c.co2_below_detection) {
            s = "Below detection limit (" + fmt(c.co2_detection_limit_ppm, 1) + " ppm)";
        } else {
            s = fmt(c.co2_ppm, 1) + " ppm";
            if (c.co2_unc_ppm !== undefined) s += " +/- " + fmt(c.co2_unc_ppm, 1);
            if (c.co2_detection_limit_ppm !== undefined) s += "  (DL: " + fmt(c.co2_detection_limit_ppm, 1) + " ppm)";
        }
        rows.push(["CO2", s]);
    }
    if (c.density_h2o_iterations !== undefined) {
        let s = c.density_h2o_iterations + " iterations";
        if (c.density_h2o_converged) s += " (converged)";
        else s += " (not converged)";
        rows.push(["Density-H2O iteration", s]);
    }
    return rows;
}

/* -----------------------------------------------------------
   7. Render figures
   ----------------------------------------------------------- */
function renderFigures(figures) {
    const grid = document.getElementById("figures-grid");
    grid.innerHTML = "";

    const titles = {
        "h2o_baseline": "H2O Baseline Fit",
        "co3_fit": "CO3 Doublet Fit",
        "fringes": "Reflectance Fringes (Mode C)"
    };

    for (const [key, dataUri] of Object.entries(figures)) {
        if (!dataUri) continue;
        const wrapper = document.createElement("div");
        wrapper.className = "figure-wrapper";
        wrapper.innerHTML =
            "<h3>" + (titles[key] || key) + "</h3>" +
            '<img src="' + dataUri + '" alt="' + key + '">';
        grid.appendChild(wrapper);
    }

    // Hide card if no figures
    document.getElementById("card-figures").style.display =
        grid.children.length > 0 ? "block" : "none";
}

/* -----------------------------------------------------------
   8. Download results as CSV
   ----------------------------------------------------------- */
function downloadCSV() {
    if (!lastResultData) return;
    const d = lastResultData;
    const headers = [];
    const values = [];

    function add(label, val) {
        headers.push(label);
        values.push(val === null || val === undefined ? "" : val);
    }

    // Physical
    const t = d.thickness || {};
    const tVal = t.thickness_um !== undefined ? t.thickness_um : t.average_um;
    add("Thickness (um)", tVal);
    add("Thickness unc (um)", t.stdev_um);
    if (t.method) add("Thickness method", t.method);

    const dens = d.density || {};
    add("Density (g/cc)", dens.density_gcc || dens.density_g_per_cc);

    const ri = d.refractive_index || {};
    if (ri.n !== undefined) add("Refractive index", ri.n);

    if (d.epsilon) {
        if (d.epsilon.epsilon_oh !== undefined) add("Epsilon OH", d.epsilon.epsilon_oh);
        if (d.epsilon.epsilon_h2omol !== undefined) add("Epsilon H2Omol", d.epsilon.epsilon_h2omol);
        if (d.epsilon.epsilon_co3 !== undefined) add("Epsilon CO3", d.epsilon.epsilon_co3);
    }

    // H2O
    const h = d.h2o_peak || {};
    if (h.peak_absorbance !== undefined) add("H2O peak abs", h.peak_absorbance);
    if (h.peak_wavenumber !== undefined) add("H2O peak wn (cm-1)", h.peak_wavenumber);
    if (h.peak_absorbance_unc !== undefined) add("H2O abs unc", h.peak_absorbance_unc);

    // CO3
    const c3 = d.co3_fit || {};
    if (c3.model !== undefined) add("CO3 model", c3.model);
    if (c3.co3_absorbance !== undefined) add("CO3 abs", c3.co3_absorbance);
    if (c3.co3_total_unc !== undefined) add("CO3 unc", c3.co3_total_unc);
    if (c3.h2omol_1630_absorbance !== undefined) add("H2Omol 1630 abs", c3.h2omol_1630_absorbance);
    if (c3.r_squared !== undefined) add("R-squared", c3.r_squared);
    if (c3.quality_flag !== undefined) add("Quality", c3.quality_flag);
    if (c3.shift_1430 !== undefined) add("Shift 1430 (cm-1)", c3.shift_1430);
    if (c3.shift_1515 !== undefined) add("Shift 1515 (cm-1)", c3.shift_1515);
    if (c3.fringe_detected !== undefined) add("Fringe detected", c3.fringe_detected);
    if (c3.fringe_corrected !== undefined) add("Fringe corrected", c3.fringe_corrected);
    if (c3.doublet_ratio !== undefined) add("Doublet ratio", c3.doublet_ratio);

    // Concentrations
    const cc = d.concentration || {};
    if (cc.total_h2o_wt_pct !== undefined) add("Total H2O (wt%)", cc.total_h2o_wt_pct);
    if (cc.total_h2o_unc_wt_pct !== undefined) add("H2O unc (wt%)", cc.total_h2o_unc_wt_pct);
    if (cc.h2omol_1630_wt_pct !== undefined) add("H2Omol 1630 (wt%)", cc.h2omol_1630_wt_pct);
    if (cc.co2_ppm !== undefined) add("CO2 (ppm)", cc.co2_below_detection ? "BDL" : cc.co2_ppm);
    if (cc.co2_unc_ppm !== undefined) add("CO2 unc (ppm)", cc.co2_unc_ppm);
    if (cc.co2_detection_limit_ppm !== undefined) add("CO2 DL (ppm)", cc.co2_detection_limit_ppm);
    if (cc.co2_below_detection !== undefined) add("CO2 below DL", cc.co2_below_detection);
    if (cc.density_h2o_iterations !== undefined) add("Density-H2O iterations", cc.density_h2o_iterations);
    if (cc.density_h2o_converged !== undefined) add("Density-H2O converged", cc.density_h2o_converged);

    // Build CSV: row 1 = headers, row 2 = values
    const esc = v => {
        const s = String(v);
        return s.includes(",") ? '"' + s + '"' : s;
    };
    const csvContent = headers.map(esc).join(",") + "\n" + values.map(esc).join(",");

    // Trigger download
    const sampleName = d.sample_name || "ftir_results";
    const blob = new Blob([csvContent], { type: "text/csv;charset=utf-8;" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = sampleName + "_results.csv";
    link.click();
    URL.revokeObjectURL(link.href);
}

/* -----------------------------------------------------------
   9. Download figures as ZIP
   ----------------------------------------------------------- */
function downloadFigures() {
    if (!lastResultData || !lastResultData.figures) return;
    const figures = lastResultData.figures;
    const entries = Object.entries(figures).filter(([k, v]) => v);
    if (entries.length === 0) return;

    // If only one figure, download directly as PNG
    if (entries.length === 1) {
        const [key, dataUri] = entries[0];
        const sampleName = lastResultData.sample_name || "ftir";
        const link = document.createElement("a");
        link.href = dataUri;
        link.download = sampleName + "_" + key + ".png";
        link.click();
        return;
    }

    // Multiple figures: use JSZip if available, otherwise download individually
    if (typeof JSZip !== "undefined") {
        const zip = new JSZip();
        const sampleName = lastResultData.sample_name || "ftir";
        for (const [key, dataUri] of entries) {
            const base64 = dataUri.split(",")[1];
            zip.file(sampleName + "_" + key + ".png", base64, { base64: true });
        }
        zip.generateAsync({ type: "blob" }).then(blob => {
            const link = document.createElement("a");
            link.href = URL.createObjectURL(blob);
            link.download = sampleName + "_figures.zip";
            link.click();
            URL.revokeObjectURL(link.href);
        });
    } else {
        // Fallback: download each figure separately
        const sampleName = lastResultData.sample_name || "ftir";
        for (const [key, dataUri] of entries) {
            const link = document.createElement("a");
            link.href = dataUri;
            link.download = sampleName + "_" + key + ".png";
            link.click();
        }
    }
}

/* -----------------------------------------------------------
   10. Interactive CO3 Fit Chart (Plotly)
   ----------------------------------------------------------- */
function renderCO3Interactive(data) {
    const card = document.getElementById("card-co3-interactive");
    const c = data.co3_fit || {};
    const wn = c.wavenumber;
    const obs = c.data;
    const fit = c.fitted;
    const co3sub = c.co3_subtracted;
    const bkg = c.background;
    const baseline = c.baseline;
    const fitRange = c.fit_range;
    const fitWn = c.fit_wavenumber;
    const resid = c.residuals;

    if (!wn || !obs || !fit) {
        card.style.display = "none";
        return;
    }
    card.style.display = "block";

    const model = c.model || "";
    const r2 = c.r_squared;
    const sampleName = data.sample_name || "";
    const shift1430 = c.shift_1430;
    const shift1515 = c.shift_1515;

    // --- Panel 1: Data + Fit + Components ---
    const trData = {
        x: wn, y: obs, mode: "lines", name: "Sample data",
        line: { color: "black", width: 1 }
    };
    const trFit = {
        x: wn, y: fit, mode: "lines", name: "Fitted model",
        line: { color: "red", width: 1.5 }
    };
    const traces1 = [trData, trFit];

    // fixed/taylor: show background (fit without CO3)
    // pca_shift: show PCA baseline
    if (model === "pca_shift") {
        if (baseline) {
            traces1.push({
                x: wn, y: baseline, mode: "lines", name: "PCA baseline",
                line: { color: "green", width: 1, dash: "dash" }
            });
        }
    } else {
        if (bkg) {
            traces1.push({
                x: wn, y: bkg, mode: "lines", name: "Background (excl. CO3)",
                line: { color: "green", width: 1, dash: "dash" }
            });
        }
    }

    // --- Panel 2: CO3 by subtraction (1300-1900 range) ---
    const traces2 = [];
    if (co3sub) {
        // Filter to 1300-1900 cm-1
        const idx2 = [];
        for (let i = 0; i < wn.length; i++) {
            if (wn[i] >= 1300 && wn[i] <= 1900) idx2.push(i);
        }
        const wn2 = idx2.map(i => wn[i]);
        const co3s2 = idx2.map(i => co3sub[i]);
        traces2.push({
            x: wn2, y: co3s2, mode: "lines", name: "CO3 by subtraction",
            line: { color: "purple", width: 1 },
            fill: "tozeroy", fillcolor: "rgba(147,112,219,0.2)"
        });
        traces2.push({
            x: wn2, y: wn2.map(() => 0), mode: "lines", showlegend: false,
            line: { color: "gray", width: 0.5, dash: "dash" }
        });
        // Mark peaks
        if (c.co3_peak_1515 !== undefined) {
            const peakWn1515 = findPeakWn(wn, co3sub, 1496, 1532);
            if (peakWn1515) {
                traces2.push({
                    x: [peakWn1515.wn], y: [peakWn1515.val], mode: "markers+text",
                    name: "1515 peak", marker: { color: "red", size: 8, symbol: "diamond" },
                    text: [peakWn1515.val.toFixed(4) + " @ " + peakWn1515.wn.toFixed(0)],
                    textposition: "top right", textfont: { size: 10 }
                });
            }
        }
        if (c.co3_peak_1430 !== undefined) {
            const peakWn1430 = findPeakWn(wn, co3sub, 1410, 1460);
            if (peakWn1430) {
                traces2.push({
                    x: [peakWn1430.wn], y: [peakWn1430.val], mode: "markers+text",
                    name: "1430 peak", marker: { color: "red", size: 8, symbol: "diamond" },
                    text: [peakWn1430.val.toFixed(4) + " @ " + peakWn1430.wn.toFixed(0)],
                    textposition: "top left", textfont: { size: 10 }
                });
            }
        }
    }

    // --- Panel 3: Residuals ---
    const traces3 = [];
    if (fitWn && resid) {
        traces3.push({
            x: fitWn, y: resid, mode: "lines", name: "Residual (fit - data)",
            line: { color: "black", width: 0.8 },
            fill: "tozeroy", fillcolor: "rgba(250,128,114,0.3)"
        });
        traces3.push({
            x: fitWn, y: fitWn.map(() => 0), mode: "lines", showlegend: false,
            line: { color: "gray", width: 0.5, dash: "dash" }
        });
    }

    // --- Title ---
    let titleExtra = "";
    if (shift1430 !== null && shift1430 !== undefined) {
        titleExtra = "  |  d1430=" + (shift1430 >= 0 ? "+" : "") + shift1430.toFixed(1) +
                     "  d1515=" + (shift1515 >= 0 ? "+" : "") + shift1515.toFixed(1) + " cm-1";
    }
    const modelDisplay = {
        "fixed": "Simple lstsq",
        "taylor": "lstsq with shifted peak position",
        "pca_shift": "PCA baseline"
    }[model] || model;
    const mainTitle = "CO3 Carbonate Fit [" + modelDisplay + "] - " + sampleName +
                      "<br>R2 = " + (r2 !== undefined ? r2.toFixed(6) : "--") + titleExtra;

    // Determine x-axis range for Panel 1 (match static matplotlib: 1200-2000)
    let p1xLo = 1200, p1xHi = 2000;

    // --- Build Plotly figure with subplots ---
    const allTraces = [];
    traces1.forEach(t => { t.xaxis = "x"; t.yaxis = "y"; allTraces.push(t); });
    traces2.forEach(t => { t.xaxis = "x2"; t.yaxis = "y2"; allTraces.push(t); });
    traces3.forEach(t => { t.xaxis = "x3"; t.yaxis = "y3"; allTraces.push(t); });

    const layout = {
        grid: { rows: 3, columns: 1, subplots: [["xy"], ["x2y2"], ["x3y3"]], roworder: "top to bottom" },
        height: 750,
        margin: { t: 60, b: 50, l: 60, r: 30 },
        title: { text: mainTitle, font: { size: 13 } },
        showlegend: true,
        legend: { font: { size: 10 }, x: 0, y: 1, bgcolor: "rgba(255,255,255,0.7)" },
        // Panel 1
        xaxis: { autorange: "reversed", range: [p1xHi, p1xLo], title: "", showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        yaxis: { title: "Absorbance", domain: [0.58, 1.0], range: [0, 3], showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        // Panel 2
        xaxis2: { autorange: "reversed", range: [1900, 1300], title: "", showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        yaxis2: { title: "Absorbance (subtracted)", domain: [0.27, 0.53], showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        // Panel 3
        xaxis3: { autorange: "reversed", title: "Wavenumber (cm-1)", showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        yaxis3: { title: "Residual", domain: [0.0, 0.22], showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        plot_bgcolor: "#fff",
        paper_bgcolor: "#fff",
    };

    const config = {
        responsive: true,
        displayModeBar: true,
        modeBarButtonsToRemove: ["lasso2d", "select2d"],
        toImageButtonOptions: {
            format: "png", filename: (sampleName || "co3_fit") + "_interactive",
            height: 900, width: 1200, scale: 2
        }
    };

    Plotly.newPlot("co3-plotly-chart", allTraces, layout, config);
}

function findPeakWn(wn, arr, lo, hi) {
    let bestVal = -Infinity, bestWn = null;
    for (let i = 0; i < wn.length; i++) {
        if (wn[i] >= lo && wn[i] <= hi && arr[i] > bestVal) {
            bestVal = arr[i];
            bestWn = wn[i];
        }
    }
    return bestWn !== null ? { wn: bestWn, val: bestVal } : null;
}

/* -----------------------------------------------------------
   11. Download CO3 Fit Data as CSV
   ----------------------------------------------------------- */
function downloadCO3FitCSV() {
    if (!lastResultData) return;
    const c = lastResultData.co3_fit || {};
    const wn = c.wavenumber;
    const obs = c.data;
    const fit = c.fitted;
    const co3sub = c.co3_subtracted;
    const bkg = c.background;

    if (!wn || !obs) return;

    const lines = ["wavenumber_cm-1,data,fitted,background,co3_subtracted"];
    for (let i = 0; i < wn.length; i++) {
        lines.push([
            wn[i].toFixed(4),
            obs[i] !== undefined ? obs[i].toFixed(6) : "",
            fit ? fit[i].toFixed(6) : "",
            bkg ? bkg[i].toFixed(6) : "",
            co3sub ? co3sub[i].toFixed(6) : ""
        ].join(","));
    }

    // Append residuals section
    const fitWn = c.fit_wavenumber;
    const resid = c.residuals;
    if (fitWn && resid) {
        lines.push("");
        lines.push("# Residuals in fitting range");
        lines.push("wavenumber_cm-1,residual");
        for (let i = 0; i < fitWn.length; i++) {
            lines.push(fitWn[i].toFixed(4) + "," + resid[i].toFixed(6));
        }
    }

    const sampleName = lastResultData.sample_name || "ftir";
    const blob = new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8;" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = sampleName + "_co3_fit_data.csv";
    link.click();
    URL.revokeObjectURL(link.href);
}

/* -----------------------------------------------------------
   12. Interactive H2O Baseline Chart (Plotly)
   ----------------------------------------------------------- */
function renderH2OInteractive(data) {
    const card = document.getElementById("card-h2o-interactive");
    const h = data.h2o_peak || {};
    const wn = h.wavenumber;
    const raw = h.raw_absorbance;

    if (!wn || !raw) {
        card.style.display = "none";
        return;
    }
    card.style.display = "block";

    // Store spectrum data for recalculation
    window._h2oData = { wn, raw, sampleName: data.sample_name || "" };
    window._h2oClickNext = "low"; // which anchor the next click sets

    // Set anchor input values from backend defaults
    const lowRange = h.baseline_low_range || [2200, 2400];
    const highRange = h.baseline_high_range || [3700, 4000];
    document.getElementById("h2o_anchor_low_min").value = lowRange[0];
    document.getElementById("h2o_anchor_low_max").value = lowRange[1];
    document.getElementById("h2o_anchor_high_min").value = highRange[0];
    document.getElementById("h2o_anchor_high_max").value = highRange[1];

    // Default point-mode values (midpoints of ranges)
    document.getElementById("h2o_point_low").value = ((lowRange[0] + lowRange[1]) / 2).toFixed(0);
    document.getElementById("h2o_point_high").value = ((highRange[0] + highRange[1]) / 2).toFixed(0);

    // Reset to range mode
    document.getElementById("h2o_anchor_mode").value = "range";
    switchAnchorMode();

    // Initial plot with server-computed data
    _h2oPlotWithAnchors(wn, raw, h.baseline, h.baseline_corrected,
                        h.peak_height, h.peak_wavenumber || 3550,
                        lowRange, highRange, null, data.sample_name || "");
    document.getElementById("h2o-peak-display").textContent =
        "Peak height: " + (h.peak_height !== undefined ? h.peak_height.toFixed(5) : "--") +
        " @ " + (h.peak_wavenumber || 3550).toFixed(0) + " cm-1";
}

// anchorPts: null for range mode, or [{wn, ab}, {wn, ab}] for point mode
function _h2oPlotWithAnchors(wn, raw, baseline, corrected, peakHeight, peakWn,
                              lowRange, highRange, anchorPts, sampleName) {
    // Panel 1: raw spectrum + baseline
    const traces = [
        { x: wn, y: raw, mode: "lines", name: "Raw spectrum",
          line: { color: "black", width: 1 }, xaxis: "x", yaxis: "y" },
        { x: wn, y: baseline, mode: "lines", name: "Linear baseline",
          line: { color: "red", width: 1.5, dash: "dash" }, xaxis: "x", yaxis: "y" },
    ];

    // Anchor points (red dots)
    if (anchorPts) {
        // Point mode: explicit anchor positions
        traces.push({
            x: [anchorPts[0].wn, anchorPts[1].wn],
            y: [anchorPts[0].ab, anchorPts[1].ab],
            mode: "markers", name: "Anchor points",
            marker: { color: "red", size: 9, symbol: "circle" },
            xaxis: "x", yaxis: "y"
        });
    } else if (lowRange && highRange) {
        // Range mode: compute averaged anchor positions
        let loWnSum = 0, loAbSum = 0, loN = 0;
        let hiWnSum = 0, hiAbSum = 0, hiN = 0;
        for (let i = 0; i < wn.length; i++) {
            if (wn[i] >= lowRange[0] && wn[i] <= lowRange[1]) {
                loWnSum += wn[i]; loAbSum += raw[i]; loN++;
            }
            if (wn[i] >= highRange[0] && wn[i] <= highRange[1]) {
                hiWnSum += wn[i]; hiAbSum += raw[i]; hiN++;
            }
        }
        if (loN > 0 && hiN > 0) {
            traces.push({
                x: [loWnSum / loN, hiWnSum / hiN],
                y: [loAbSum / loN, hiAbSum / hiN],
                mode: "markers", name: "Anchor points",
                marker: { color: "red", size: 9, symbol: "circle" },
                xaxis: "x", yaxis: "y"
            });
        }
    }

    // Panel 2: corrected spectrum + peak marker
    traces.push({
        x: wn, y: corrected, mode: "lines", name: "Baseline-corrected",
        line: { color: "steelblue", width: 1 },
        fill: "tozeroy", fillcolor: "rgba(70,130,180,0.2)",
        xaxis: "x2", yaxis: "y2"
    });
    traces.push({
        x: wn, y: wn.map(() => 0), mode: "lines", showlegend: false,
        line: { color: "gray", width: 0.5, dash: "dash" },
        xaxis: "x2", yaxis: "y2"
    });

    if (peakHeight !== undefined && peakWn !== undefined) {
        traces.push({
            x: [peakWn], y: [peakHeight], mode: "markers+text",
            name: "Peak @ " + peakWn.toFixed(0),
            marker: { color: "red", size: 10, symbol: "triangle-down" },
            text: [peakHeight.toFixed(4) + " @ " + peakWn.toFixed(0)],
            textposition: "top right", textfont: { size: 11 },
            xaxis: "x2", yaxis: "y2"
        });
    }

    // Shapes: anchor region shading (range mode only)
    const shapes = [];
    if (!anchorPts && lowRange && highRange) {
        shapes.push(
            { type: "rect", xref: "x", yref: "paper", x0: lowRange[0], x1: lowRange[1],
              y0: 0.55, y1: 1.0, fillcolor: "green", opacity: 0.12, line: { width: 0 } },
            { type: "rect", xref: "x", yref: "paper", x0: highRange[0], x1: highRange[1],
              y0: 0.55, y1: 1.0, fillcolor: "orange", opacity: 0.12, line: { width: 0 } }
        );
    }

    const layout = {
        grid: { rows: 2, columns: 1, subplots: [["xy"], ["x2y2"]], roworder: "top to bottom" },
        height: 550,
        margin: { t: 50, b: 50, l: 60, r: 30 },
        title: { text: "H2O Baseline - " + sampleName, font: { size: 13 } },
        showlegend: true,
        legend: { font: { size: 10 }, x: 0, y: 1, bgcolor: "rgba(255,255,255,0.7)" },
        xaxis: { autorange: "reversed", showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        yaxis: { title: "Absorbance", domain: [0.55, 1.0], showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        xaxis2: { autorange: "reversed", title: "Wavenumber (cm-1)", showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        yaxis2: { title: "Corrected Absorbance", domain: [0.0, 0.48], showgrid: true, gridcolor: "rgba(0,0,0,0.1)" },
        shapes: shapes,
        plot_bgcolor: "#fff",
        paper_bgcolor: "#fff",
    };

    const config = {
        responsive: true, displayModeBar: true,
        modeBarButtonsToRemove: ["lasso2d", "select2d"],
        toImageButtonOptions: {
            format: "png", filename: (sampleName || "h2o_baseline") + "_interactive",
            height: 700, width: 1200, scale: 2
        }
    };

    Plotly.newPlot("h2o-plotly-chart", traces, layout, config);

    // Attach click handler for point-pick mode
    const chartDiv = document.getElementById("h2o-plotly-chart");
    chartDiv.removeAllListeners("plotly_click");
    chartDiv.on("plotly_click", function(eventData) {
        if (document.getElementById("h2o_anchor_mode").value !== "point") return;
        if (!eventData.points || eventData.points.length === 0) return;
        // Only respond to clicks on Panel 1 (yaxis "y", not "y2")
        const pt = eventData.points[0];
        if (pt.yaxis && pt.yaxis._id !== "y") return;

        const clickedWn = pt.x;
        if (window._h2oClickNext === "low") {
            document.getElementById("h2o_point_low").value = clickedWn.toFixed(1);
            window._h2oClickNext = "high";
        } else {
            document.getElementById("h2o_point_high").value = clickedWn.toFixed(1);
            window._h2oClickNext = "low";
        }
        recalcH2OPoint();
    });
}

/* -----------------------------------------------------------
   13. Client-side H2O baseline recalculation
   ----------------------------------------------------------- */

// Shared: compute baseline from two anchor points, update plot + display
function _h2oApplyBaseline(loWn, loAb, hiWn, hiAb, anchorPts, lowRange, highRange) {
    if (!window._h2oData) return;
    const { wn, raw, sampleName } = window._h2oData;

    const slope = (hiAb - loAb) / (hiWn - loWn);
    const intercept = loAb - slope * loWn;
    const baseline = wn.map(w => slope * w + intercept);
    const corrected = raw.map((a, i) => a - baseline[i]);

    const peakWn = 3550.0;
    const peakHeight = _interpAt(wn, corrected, peakWn);

    document.getElementById("h2o-peak-display").textContent =
        "Peak height: " + peakHeight.toFixed(5) + " @ " + peakWn.toFixed(0) + " cm-1";
    document.getElementById("h2o-anchor-status").textContent = "";

    _h2oPlotWithAnchors(wn, raw, baseline, corrected, peakHeight, peakWn,
                        lowRange, highRange, anchorPts, sampleName);
}

// --- Range mode recalculation ---
function recalcH2OBaseline() {
    if (!window._h2oData) return;
    const { wn, raw } = window._h2oData;

    const lowMin = parseFloat(document.getElementById("h2o_anchor_low_min").value);
    const lowMax = parseFloat(document.getElementById("h2o_anchor_low_max").value);
    const highMin = parseFloat(document.getElementById("h2o_anchor_high_min").value);
    const highMax = parseFloat(document.getElementById("h2o_anchor_high_max").value);

    if (isNaN(lowMin) || isNaN(lowMax) || isNaN(highMin) || isNaN(highMax)) return;
    if (lowMin >= lowMax || highMin >= highMax) {
        document.getElementById("h2o-anchor-status").textContent = "Invalid range: min must be less than max.";
        return;
    }
    if (lowMax >= highMin) {
        document.getElementById("h2o-anchor-status").textContent = "Anchor ranges must not overlap.";
        return;
    }

    let loWnSum = 0, loAbSum = 0, loN = 0;
    for (let i = 0; i < wn.length; i++) {
        if (wn[i] >= lowMin && wn[i] <= lowMax) {
            loWnSum += wn[i]; loAbSum += raw[i]; loN++;
        }
    }
    if (loN < 3) {
        document.getElementById("h2o-anchor-status").textContent =
            "Not enough points in low anchor range (" + loN + " found, need >= 3).";
        return;
    }

    let hiWnSum = 0, hiAbSum = 0, hiN = 0;
    for (let i = 0; i < wn.length; i++) {
        if (wn[i] >= highMin && wn[i] <= highMax) {
            hiWnSum += wn[i]; hiAbSum += raw[i]; hiN++;
        }
    }
    if (hiN < 3) {
        document.getElementById("h2o-anchor-status").textContent =
            "Not enough points in high anchor range (" + hiN + " found, need >= 3).";
        return;
    }

    _h2oApplyBaseline(loWnSum / loN, loAbSum / loN, hiWnSum / hiN, hiAbSum / hiN,
                      null, [lowMin, lowMax], [highMin, highMax]);
}

// --- Point mode recalculation ---
function recalcH2OPoint() {
    if (!window._h2oData) return;
    const { wn, raw } = window._h2oData;

    const loWn = parseFloat(document.getElementById("h2o_point_low").value);
    const hiWn = parseFloat(document.getElementById("h2o_point_high").value);
    if (isNaN(loWn) || isNaN(hiWn)) return;
    if (loWn >= hiWn) {
        document.getElementById("h2o-anchor-status").textContent = "Low point must be less than high point.";
        return;
    }

    // Interpolate absorbance at selected wavenumbers
    const loAb = _interpAt(wn, raw, loWn);
    const hiAb = _interpAt(wn, raw, hiWn);

    _h2oApplyBaseline(loWn, loAb, hiWn, hiAb,
                      [{ wn: loWn, ab: loAb }, { wn: hiWn, ab: hiAb }],
                      null, null);
}

// --- Mode switching ---
function switchAnchorMode() {
    const mode = document.getElementById("h2o_anchor_mode").value;
    document.getElementById("h2o-range-controls").style.display = mode === "range" ? "flex" : "none";
    document.getElementById("h2o-point-controls").style.display = mode === "point" ? "flex" : "none";
    window._h2oClickNext = "low";
    if (mode === "range") recalcH2OBaseline();
    else recalcH2OPoint();
}

function _interpAt(xArr, yArr, xTarget) {
    const pairs = xArr.map((x, i) => [x, yArr[i]]).sort((a, b) => a[0] - b[0]);
    if (xTarget <= pairs[0][0]) return pairs[0][1];
    if (xTarget >= pairs[pairs.length - 1][0]) return pairs[pairs.length - 1][1];
    for (let i = 0; i < pairs.length - 1; i++) {
        if (pairs[i][0] <= xTarget && pairs[i + 1][0] >= xTarget) {
            const t = (xTarget - pairs[i][0]) / (pairs[i + 1][0] - pairs[i][0]);
            return pairs[i][1] + t * (pairs[i + 1][1] - pairs[i][1]);
        }
    }
    return pairs[0][1];
}

function resetH2OAnchors() {
    if (!lastResultData) return;
    const h = lastResultData.h2o_peak || {};
    const mode = document.getElementById("h2o_anchor_mode").value;

    const lowRange = h.baseline_low_range || [2200, 2400];
    const highRange = h.baseline_high_range || [3700, 4000];
    document.getElementById("h2o_anchor_low_min").value = lowRange[0];
    document.getElementById("h2o_anchor_low_max").value = lowRange[1];
    document.getElementById("h2o_anchor_high_min").value = highRange[0];
    document.getElementById("h2o_anchor_high_max").value = highRange[1];
    document.getElementById("h2o_point_low").value = ((lowRange[0] + lowRange[1]) / 2).toFixed(0);
    document.getElementById("h2o_point_high").value = ((highRange[0] + highRange[1]) / 2).toFixed(0);
    window._h2oClickNext = "low";

    if (mode === "range") recalcH2OBaseline();
    else recalcH2OPoint();
}

// Wire up anchor input event listeners
["h2o_anchor_low_min", "h2o_anchor_low_max",
 "h2o_anchor_high_min", "h2o_anchor_high_max"].forEach(function(id) {
    document.getElementById(id).addEventListener("change", recalcH2OBaseline);
});
["h2o_point_low", "h2o_point_high"].forEach(function(id) {
    document.getElementById(id).addEventListener("change", recalcH2OPoint);
});

/* -----------------------------------------------------------
   14. Download H2O Baseline Data as CSV
   ----------------------------------------------------------- */
function downloadH2OBaselineCSV() {
    if (!lastResultData) return;
    const h = lastResultData.h2o_peak || {};
    const wn = h.wavenumber;
    const raw = h.raw_absorbance;
    const bl = h.baseline;
    const corr = h.baseline_corrected;
    if (!wn || !raw) return;

    const lines = ["wavenumber_cm-1,raw_absorbance,baseline,corrected"];
    for (let i = 0; i < wn.length; i++) {
        lines.push([
            wn[i].toFixed(4),
            raw[i].toFixed(6),
            bl ? bl[i].toFixed(6) : "",
            corr ? corr[i].toFixed(6) : ""
        ].join(","));
    }

    const sampleName = lastResultData.sample_name || "ftir";
    const blob = new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8;" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = sampleName + "_h2o_baseline_data.csv";
    link.click();
    URL.revokeObjectURL(link.href);
}
