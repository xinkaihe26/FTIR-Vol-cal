"""
FTIR Web Application — Flask backend.

Usage:
    cd FTIR/web
    python app.py
    Open http://localhost:5000 in a browser.
"""

import sys
import os
import json
import base64
import shutil
import tempfile
from pathlib import Path

# Add parent directory so we can import ftir_tools
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")

from flask import Flask, render_template, request, jsonify
from ftir_tools import process_sample

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024  # 16 MB

# ---------------------------------------------------------------------------
# Default parameter values
# ---------------------------------------------------------------------------
DEFAULTS = {
    "composition": {
        "SiO2": 49.2, "TiO2": 3.18, "Al2O3": 12.82,
        "FeOT": 14.41, "MnO": 0.22, "MgO": 5.64,
        "CaO": 9.37, "Na2O": 3.22, "K2O": 0.84,
    },
    "fe3_ratio": 0.15,
    "co3_model": "fixed",
    "thickness_mode": "B",
    "thickness_um": 27.4,
    "thickness_unc_um": 0.5,
    "olivine_fo": 85,
    "epsilon_oh": 62.3,
    "epsilon_h2omol": 25.0,
    "epsilon_co3": 375.0,
    "use_composition_epsilon": False,
    "fringe_correction": "auto",
}


# ---------------------------------------------------------------------------
# Helper: encode PNG file as base64 data URI
# ---------------------------------------------------------------------------
_CO3_MODEL_MAP = {
    "Simple lstsq": "fixed",
    "lstsq with shifted peak position": "taylor",
    "PCA baseline": "pca_shift",
    # Also accept internal names directly
    "fixed": "fixed",
    "taylor": "taylor",
    "pca_shift": "pca_shift",
}


def _resolve_co3_model(name: str) -> str:
    """Map display name or internal name to internal co3_model value."""
    return _CO3_MODEL_MAP.get(name.strip(), name.strip())


def _encode_figure(fig_path):
    """Read a PNG file and return a base64 data URI string."""
    if fig_path and Path(fig_path).exists():
        data = Path(fig_path).read_bytes()
        b64 = base64.b64encode(data).decode("ascii")
        return f"data:image/png;base64,{b64}"
    return None


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------
@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/defaults", methods=["GET"])
def get_defaults():
    return jsonify(DEFAULTS)


@app.route("/api/changelog", methods=["GET"])
def get_changelog():
    changelog_path = Path(__file__).resolve().parent / "changelog.json"
    if changelog_path.exists():
        return jsonify(json.loads(changelog_path.read_text(encoding="utf-8")))
    return jsonify([])


@app.route("/api/process", methods=["POST"])
def run_analysis():
    tmp_dir = tempfile.mkdtemp(prefix="ftir_web_")
    try:
        # --- Save uploaded files ---
        trans_file = request.files.get("transmission_file")
        if not trans_file or trans_file.filename == "":
            return jsonify({"error": "Transmission spectrum file is required."}), 400
        trans_path = Path(tmp_dir) / "transmission.txt"
        trans_file.save(str(trans_path))

        ref_path = None
        ref_file = request.files.get("reflectance_file")
        if ref_file and ref_file.filename != "":
            ref_path = Path(tmp_dir) / "reflectance.txt"
            ref_file.save(str(ref_path))

        # --- Parse composition ---
        comp_str = request.form.get("composition", "{}")
        try:
            composition = json.loads(comp_str)
        except json.JSONDecodeError:
            return jsonify({"error": "Invalid composition JSON."}), 400
        if not composition:
            return jsonify({"error": "Glass composition is required."}), 400

        # --- Build process_sample kwargs ---
        kwargs = {
            "spectrum_csv": str(trans_path),
            "composition_dict": composition,
            "save_figures": True,
            "output_dir": tmp_dir,
        }

        # Thickness mode
        mode = request.form.get("thickness_mode", "B")
        if mode == "A":
            meas_str = request.form.get("thickness_measurements", "[]")
            try:
                meas = json.loads(meas_str)
                kwargs["thickness_measurements"] = [tuple(m) for m in meas]
            except (json.JSONDecodeError, TypeError):
                return jsonify({"error": "Invalid thickness measurements."}), 400
            fo_a = request.form.get("olivine_fo_a")
            if fo_a:
                kwargs["olivine_fo_a"] = float(fo_a)
        elif mode == "B":
            try:
                kwargs["thickness_um"] = float(request.form.get("thickness_um", 27.4))
                unc = request.form.get("thickness_unc_um", "0")
                kwargs["thickness_unc_um"] = float(unc) if unc else 0.0
            except ValueError:
                return jsonify({"error": "Invalid thickness value."}), 400
        elif mode == "C":
            if not ref_path:
                return jsonify({
                    "error": "Mode C requires a reflectance spectrum file."
                }), 400
            try:
                kwargs["reflectance_csv"] = str(ref_path)
                kwargs["olivine_fo"] = float(
                    request.form.get("olivine_fo", 85))
            except ValueError:
                return jsonify({"error": "Invalid olivine Fo value."}), 400
        else:
            return jsonify({"error": f"Unknown thickness mode: {mode}"}), 400

        # Optional parameters
        fe3 = request.form.get("fe3_ratio")
        if fe3:
            kwargs["fe3_ratio"] = float(fe3)

        co3_model = request.form.get("co3_model", "fixed")
        kwargs["co3_model"] = co3_model

        fringe_correction = request.form.get("fringe_correction", "auto")
        kwargs["fringe_correction"] = fringe_correction

        use_comp_eps = request.form.get("use_composition_epsilon") == "true"
        kwargs["use_composition_epsilon"] = use_comp_eps

        if not use_comp_eps:
            for key in ("epsilon_oh", "epsilon_h2omol", "epsilon_co3"):
                val = request.form.get(key)
                if val:
                    kwargs[key] = float(val)

        # --- Run analysis ---
        result = process_sample(**kwargs)

        # --- Build response ---
        resp = {
            "sample_name": result.get("sample_name", ""),
            "thickness": _clean_dict(result.get("thickness", {})),
            "density": _clean_dict(result.get("density", {})),
            "refractive_index": _clean_dict(result.get("refractive_index", {})),
            "h2o_peak": _clean_dict(result.get("h2o_peak", {})),
            "co3_fit": _clean_dict(result.get("co3_fit", {})),
            "concentration": _clean_dict(result.get("concentration", {})),
        }
        if "epsilon" in result:
            resp["epsilon"] = _clean_dict(result["epsilon"])

        # Encode figures
        figures = {}
        h2o_fig = result.get("h2o_peak", {}).get("figure_path")
        if h2o_fig:
            figures["h2o_baseline"] = _encode_figure(h2o_fig)
        co3_fig = result.get("co3_fit", {}).get("figure_path")
        if co3_fig:
            figures["co3_fit"] = _encode_figure(co3_fig)
        # Mode C fringe figure
        t_fig = result.get("thickness", {}).get("figure_path")
        if t_fig:
            figures["fringes"] = _encode_figure(t_fig)
        resp["figures"] = figures

        return jsonify(resp)

    except Exception as e:
        return jsonify({"error": str(e)}), 500

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


@app.route("/api/batch", methods=["POST"])
def run_batch():
    """Batch process multiple spectra. Returns a ZIP file."""
    import csv
    import io
    import zipfile

    tmp_dir = tempfile.mkdtemp(prefix="ftir_batch_")
    try:
        # --- Collect uploaded spectrum files ---
        files = request.files.getlist("batch_files")
        if not files or all(f.filename == "" for f in files):
            return jsonify({"error": "No spectrum files uploaded."}), 400

        # Save all transmission spectrum files to tmp_dir/trans/
        trans_dir = Path(tmp_dir) / "trans"
        trans_dir.mkdir()
        file_paths = {}
        for f in files:
            if f.filename == "":
                continue
            fpath = trans_dir / f.filename
            f.save(str(fpath))
            file_paths[f.filename] = str(fpath)

        # Save reflectance spectrum files to tmp_dir/ref/ (for Mode C)
        ref_dir = Path(tmp_dir) / "ref"
        ref_dir.mkdir()
        ref_paths = {}
        ref_files = request.files.getlist("batch_ref_files")
        for f in ref_files:
            if f.filename == "":
                continue
            fpath = ref_dir / f.filename
            f.save(str(fpath))
            ref_paths[f.filename] = str(fpath)

        # --- Parse parameters ---
        param_mode = request.form.get("batch_param_mode", "unified")

        # Shared parameters (used for unified mode, or as defaults for csv mode)
        comp_str = request.form.get("composition", "{}")
        try:
            shared_comp = json.loads(comp_str)
        except json.JSONDecodeError:
            shared_comp = {}

        shared_params = {
            "thickness_mode": request.form.get("thickness_mode", "B"),
            "thickness_um": request.form.get("thickness_um", "27.4"),
            "thickness_unc_um": request.form.get("thickness_unc_um", "0"),
            "olivine_fo_a": request.form.get("olivine_fo_a", ""),
            "olivine_fo": request.form.get("olivine_fo", "85"),
            "co3_model": request.form.get("co3_model", "fixed"),
            "fringe_correction": request.form.get("fringe_correction", "auto"),
            "fe3_ratio": request.form.get("fe3_ratio", "0.15"),
            "use_composition_epsilon": request.form.get(
                "use_composition_epsilon") == "true",
            "epsilon_oh": request.form.get("epsilon_oh", "62.3"),
            "epsilon_h2omol": request.form.get("epsilon_h2omol", "25.0"),
            "epsilon_co3": request.form.get("epsilon_co3", "375.0"),
        }

        # --- Per-sample parameters from CSV (if provided) ---
        per_sample_params = {}
        if param_mode == "csv":
            param_file = request.files.get("batch_param_csv")
            if param_file and param_file.filename:
                text = param_file.read().decode("utf-8-sig")
                reader = csv.DictReader(io.StringIO(text))
                for row in reader:
                    fname = row.get("filename", "").strip()
                    if fname:
                        per_sample_params[fname] = row

        # --- Process each file ---
        out_dir = Path(tmp_dir) / "results"
        out_dir.mkdir()

        summary_rows = []
        errors = []

        for filename, fpath in sorted(file_paths.items()):
            try:
                # Build kwargs for this sample
                kwargs = {
                    "spectrum_csv": fpath,
                    "save_figures": True,
                    "output_dir": str(out_dir),
                }

                # Determine parameters (per-sample CSV overrides unified)
                if filename in per_sample_params:
                    p = per_sample_params[filename]
                    # Composition from CSV columns
                    comp = {}
                    for ox in ["SiO2", "TiO2", "Al2O3", "FeOT", "MnO",
                               "MgO", "CaO", "Na2O", "K2O"]:
                        val = p.get(ox, "").strip()
                        if val:
                            comp[ox] = float(val)
                    kwargs["composition_dict"] = comp if comp else shared_comp

                    # Thickness: CSV thickness_um, or reflectance (Mode C)
                    t_um = p.get("thickness_um", "").strip()
                    if t_um:
                        kwargs["thickness_um"] = float(t_um)
                        t_unc = p.get("thickness_unc_um", "0").strip()
                        kwargs["thickness_unc_um"] = float(t_unc) if t_unc else 0.0
                    elif filename in ref_paths:
                        # No thickness_um but reflectance file exists → Mode C
                        kwargs["reflectance_csv"] = ref_paths[filename]
                        fo = p.get("olivine_fo", "").strip()
                        kwargs["olivine_fo"] = float(fo) if fo else float(
                            shared_params["olivine_fo"])
                    else:
                        kwargs["thickness_um"] = float(
                            shared_params["thickness_um"])
                        kwargs["thickness_unc_um"] = float(
                            shared_params["thickness_unc_um"])

                    # Olivine Fo for Mode A
                    fo_a = p.get("olivine_fo", "").strip()
                    if fo_a:
                        kwargs["olivine_fo_a"] = float(fo_a)

                    # Other params: use per-sample if provided, else shared
                    raw_model = p.get(
                        "co3_model", shared_params["co3_model"]).strip() or \
                        shared_params["co3_model"]
                    kwargs["co3_model"] = _resolve_co3_model(raw_model)
                    kwargs["fringe_correction"] = p.get(
                        "fringe_correction",
                        shared_params["fringe_correction"]).strip() or \
                        shared_params["fringe_correction"]
                    fe3 = p.get("fe3_ratio", "").strip()
                    kwargs["fe3_ratio"] = float(fe3) if fe3 else float(
                        shared_params["fe3_ratio"])

                    # Per-sample use_composition_epsilon
                    uce = p.get("use_composition_epsilon", "").strip().upper()
                    if uce in ("TRUE", "1", "YES"):
                        kwargs["use_composition_epsilon"] = True
                    elif uce in ("FALSE", "0", "NO"):
                        kwargs["use_composition_epsilon"] = False
                else:
                    # Unified mode: use shared params
                    kwargs["composition_dict"] = shared_comp
                    mode = shared_params["thickness_mode"]
                    if mode == "B":
                        kwargs["thickness_um"] = float(
                            shared_params["thickness_um"])
                        kwargs["thickness_unc_um"] = float(
                            shared_params["thickness_unc_um"])
                    elif mode == "A":
                        # Mode A requires measurements from form
                        meas_str = request.form.get(
                            "thickness_measurements", "[]")
                        try:
                            meas = json.loads(meas_str)
                            kwargs["thickness_measurements"] = [
                                tuple(m) for m in meas]
                        except (json.JSONDecodeError, TypeError):
                            kwargs["thickness_um"] = float(
                                shared_params["thickness_um"])
                            kwargs["thickness_unc_um"] = 0.0
                        fo_a = shared_params.get("olivine_fo_a", "")
                        if fo_a:
                            kwargs["olivine_fo_a"] = float(fo_a)
                    elif mode == "C":
                        # Mode C: match reflectance file by filename
                        if filename in ref_paths:
                            kwargs["reflectance_csv"] = ref_paths[filename]
                        else:
                            raise ValueError(
                                f"No matching reflectance file for "
                                f"'{filename}'. Upload a reflectance file "
                                f"with the same filename.")
                        kwargs["olivine_fo"] = float(
                            shared_params["olivine_fo"])
                    kwargs["co3_model"] = shared_params["co3_model"]
                    kwargs["fringe_correction"] = shared_params[
                        "fringe_correction"]
                    kwargs["fe3_ratio"] = float(shared_params["fe3_ratio"])

                kwargs["use_composition_epsilon"] = shared_params[
                    "use_composition_epsilon"]
                if not shared_params["use_composition_epsilon"]:
                    kwargs["epsilon_oh"] = float(shared_params["epsilon_oh"])
                    kwargs["epsilon_h2omol"] = float(
                        shared_params["epsilon_h2omol"])
                    kwargs["epsilon_co3"] = float(shared_params["epsilon_co3"])

                # Run
                result = process_sample(**kwargs)

                # Collect summary row
                t = result.get("thickness", {})
                d = result.get("density", {})
                h = result.get("h2o_peak", {})
                c3 = result.get("co3_fit", {})
                cc = result.get("concentration", {})
                t_val = t.get("thickness_um", t.get("average_um", ""))

                summary_rows.append({
                    "filename": filename,
                    "thickness_um": t_val,
                    "thickness_unc_um": t.get("stdev_um", ""),
                    "density_gcc": d.get("density_gcc", ""),
                    "H2O_peak_abs": h.get("peak_height", ""),
                    "H2O_peak_wn": h.get("peak_wavenumber", ""),
                    "H2O_abs_unc": h.get("peak_height_std", ""),
                    "CO3_model": c3.get("model", ""),
                    "CO3_abs": c3.get("co3_absorbance", ""),
                    "CO3_unc": c3.get("co3_total_unc", ""),
                    "R_squared": c3.get("r_squared", ""),
                    "quality": c3.get("quality_flag", ""),
                    "shift_1430": c3.get("shift_1430", ""),
                    "shift_1515": c3.get("shift_1515", ""),
                    "doublet_ratio": c3.get("doublet_ratio", ""),
                    "fringe_detected": c3.get("fringe_detected", ""),
                    "fringe_corrected": c3.get("fringe_corrected", ""),
                    "total_H2O_wt%": cc.get("total_h2o_wt_pct", ""),
                    "H2O_unc_wt%": cc.get("total_h2o_unc_wt_pct", ""),
                    "H2Omol_1630_wt%": cc.get("h2omol_1630_wt_pct", ""),
                    "CO2_ppm": "BDL" if cc.get("co2_below_detection")
                               else cc.get("co2_ppm", ""),
                    "CO2_unc_ppm": cc.get("co2_unc_ppm", ""),
                    "CO2_DL_ppm": cc.get("co2_detection_limit_ppm", ""),
                    "CO2_below_DL": cc.get("co2_below_detection", ""),
                    "density_iterations": cc.get(
                        "density_h2o_iterations", ""),
                })

            except Exception as e:
                errors.append({"filename": filename, "error": str(e)})
                summary_rows.append({
                    "filename": filename,
                    "total_H2O_wt%": f"ERROR: {e}",
                })

        # --- Build ZIP ---
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            # Summary CSV
            if summary_rows:
                cols = list(summary_rows[0].keys())
                csv_buf = io.StringIO()
                writer = csv.DictWriter(csv_buf, fieldnames=cols)
                writer.writeheader()
                for row in summary_rows:
                    writer.writerow(row)
                zf.writestr("batch_results.csv", csv_buf.getvalue())

            # Figures
            for fig_file in out_dir.glob("*.png"):
                zf.write(str(fig_file), f"figures/{fig_file.name}")

            # Error log
            if errors:
                err_lines = [f"{e['filename']}: {e['error']}" for e in errors]
                zf.writestr("errors.txt", "\n".join(err_lines))

        zip_buffer.seek(0)

        from flask import send_file
        return send_file(
            zip_buffer,
            mimetype="application/zip",
            as_attachment=True,
            download_name="ftir_batch_results.zip",
        )

    except Exception as e:
        return jsonify({"error": str(e)}), 500

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _clean_dict(d):
    """Remove non-serializable items (like numpy arrays) from result dict."""
    import numpy as np
    clean = {}
    for k, v in d.items():
        if isinstance(v, dict):
            clean[k] = _clean_dict(v)
        elif isinstance(v, np.ndarray):
            clean[k] = v.tolist()
        elif isinstance(v, (np.floating, np.integer)):
            clean[k] = float(v)
        elif isinstance(v, (list, tuple)):
            clean[k] = [
                float(x) if isinstance(x, (np.floating, np.integer)) else x
                for x in v
            ]
        elif v is None or isinstance(v, (int, float, str, bool)):
            clean[k] = v
        else:
            clean[k] = str(v)
    return clean


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    debug = os.environ.get("FLASK_DEBUG", "1") == "1"
    print(f"Starting FTIR Web App on port {port}...")
    if debug:
        print("Open http://localhost:5000 in your browser.")
    app.run(debug=debug, host="0.0.0.0", port=port)
