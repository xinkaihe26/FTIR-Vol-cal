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
