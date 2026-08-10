"""
GIWAXS Easy Batch PNG -> NPZ Converter
=======================================

Purpose
-------
Convert one or many rendered GIWAXS/GIXS PNG plots into NPZ files containing:

    intensity : 2-D approximate intensity reconstructed from the PNG colormap
    qr        : 1-D horizontal reciprocal-space axis
    qz        : 1-D vertical reciprocal-space axis
    mask      : 2-D validity mask

* Select several PNG files at once.
* Enter the q-axis limits only once for the whole batch.
* The program automatically detects the colored q-space panel.
* The program automatically identifies common Matplotlib colormaps.
* Your last settings are remembered for the next run.
* No repeated console questions are required.

First run
---------
Run this file with a normal Python 3 interpreter (for example from PyCharm or
by double-clicking/running it from a terminal). Missing pip-installable packages
are installed automatically into that same Python environment.
"""

from __future__ import annotations

# Standard-library imports used for startup, settings, and file handling.
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from typing import Iterable


# =============================================================================
# FIRST-RUN DEPENDENCY SETUP
# =============================================================================
# The converter is intended to be copy/paste friendly.  If a normal Python
# package is missing, install it into the same interpreter that is running this
# script.  This avoids the common PyCharm problem where `pip` installs into a
# different Python environment.
_REQUIRED_PYTHON_PACKAGES = {
    "numpy": "numpy>=1.24",
    "PIL": "Pillow>=10.0",
    "matplotlib": "matplotlib>=3.7",
    "scipy": "scipy>=1.10",
}


def _ensure_python_packages() -> None:
    """Install missing third-party packages into the active Python interpreter."""
    missing = [
        requirement
        for module, requirement in _REQUIRED_PYTHON_PACKAGES.items()
        if importlib.util.find_spec(module) is None
    ]
    if not missing:
        return

    print("Installing required packages once:", ", ".join(missing), flush=True)

    # Most Python installations already include pip.  If not, try Python's
    # built-in ensurepip module before installing the required packages.
    if importlib.util.find_spec("pip") is None:
        subprocess.check_call([sys.executable, "-m", "ensurepip", "--upgrade"])

    try:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", *missing]
        )
    except subprocess.CalledProcessError as exc:
        command = f'"{sys.executable}" -m pip install ' + " ".join(missing)
        raise SystemExit(
            "Automatic dependency installation failed.\n"
            "Check your internet connection and Python package access, then run:\n"
            f"    {command}\n"
            f"Original installation error: {exc}"
        ) from exc


_ensure_python_packages()

# Third-party numerical/image packages.
import numpy as np
from PIL import Image
from matplotlib import colormaps
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.ndimage import binary_closing, binary_fill_holes, label
from scipy.spatial import cKDTree

# Tkinter is included with the standard python.org Windows/macOS installers.
# Some Linux distributions package it separately, so provide a useful message
# if that system component is missing instead of failing with an obscure import.
try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except ImportError as exc:
    raise SystemExit(
        "Tkinter is required for this graphical converter.\n"
        "Windows/macOS: install Python from python.org with Tcl/Tk support.\n"
        "Ubuntu/Debian Linux: sudo apt install python3-tk\n"
        f"Original import error: {exc}"
    ) from exc


APP_TITLE = "GIWAXS Easy Batch PNG → NPZ Converter"
SETTINGS_PATH = Path.home() / ".giwaxs_png_npz_converter_settings.json"
COMMON_COLORMAPS = (
    "jet",
    "turbo",
    "viridis",
    "plasma",
    "inferno",
    "magma",
    "cividis",
    "rainbow",
)

# These defaults match the example Duke stitched WAXS image supplied in the chat.
# They are editable in the GUI and are remembered after conversion.
DEFAULT_SETTINGS = {
    "qr_min": -1.0,
    "qr_max": 2.2,
    "qz_min": -0.1,
    "qz_max": 2.72,
    "auto_crop": True,
    "auto_colormap": True,
    "colormap": "jet",
    "preview_first": False,
}


def load_settings() -> dict:
    settings = dict(DEFAULT_SETTINGS)
    try:
        if SETTINGS_PATH.exists():
            loaded = json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                settings.update(loaded)
    except Exception:
        pass
    return settings


def save_settings(settings: dict) -> None:
    try:
        SETTINGS_PATH.write_text(
            json.dumps(settings, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    except Exception:
        # Conversion should not fail merely because settings could not be saved.
        pass


def normalize_rgb(array: np.ndarray) -> np.ndarray:
    """Return an RGB float array in the range 0..1."""
    array = np.asarray(array)
    if array.ndim == 2:
        array = np.repeat(array[..., None], 3, axis=2)
    if array.ndim != 3 or array.shape[2] not in (3, 4):
        raise ValueError(
            f"Expected grayscale, RGB, or RGBA PNG; received shape {array.shape}."
        )

    rgb = array[..., :3].astype(np.float64)
    if np.issubdtype(array.dtype, np.integer):
        rgb /= float(np.iinfo(array.dtype).max)
    else:
        finite = np.isfinite(rgb)
        if finite.any() and float(np.nanmax(rgb[finite])) > 1.5:
            rgb /= 255.0
    return np.clip(rgb, 0.0, 1.0)


def detect_colored_panel(rgb: np.ndarray) -> tuple[int, int, int, int]:
    """Detect the main colored q-space panel while rejecting labels/colorbars."""
    height, width, _ = rgb.shape
    saturation = np.max(rgb, axis=2) - np.min(rgb, axis=2)
    brightness = np.mean(rgb, axis=2)

    candidate = (
        (saturation > 0.055)
        & (brightness > 0.015)
        & (brightness < 0.985)
    )
    candidate = binary_closing(
        candidate,
        iterations=max(1, min(height, width) // 500),
    )
    candidate = binary_fill_holes(candidate)

    components, count = label(candidate)
    best: tuple[int, int, int, int] | None = None
    best_score = -np.inf

    for component_id in range(1, count + 1):
        ys, xs = np.nonzero(components == component_id)
        if xs.size < 50:
            continue

        x0, x1 = int(xs.min()), int(xs.max()) + 1
        y0, y1 = int(ys.min()), int(ys.max()) + 1
        box_w, box_h = x1 - x0, y1 - y0
        area = box_w * box_h
        aspect = box_w / max(box_h, 1)

        # Reject common narrow colorbar shapes and tiny components.
        if box_w < 0.10 * width or box_h < 0.18 * height:
            continue
        if aspect < 0.18 or aspect > 8.0:
            continue

        fill = xs.size / max(area, 1)
        center_penalty = abs((x0 + x1) / 2 - width / 2) / max(width, 1)
        score = area * (0.45 + fill) * (1.0 - 0.25 * center_penalty)
        if score > best_score:
            best_score = score
            best = (x0, y0, x1, y1)

    if best is None:
        return (0, 0, width, height)

    x0, y0, x1, y1 = best
    pad_x = max(1, int(0.005 * width))
    pad_y = max(1, int(0.005 * height))
    return (
        max(0, x0 - pad_x),
        max(0, y0 - pad_y),
        min(width, x1 + pad_x),
        min(height, y1 + pad_y),
    )


def choose_representative_pixels(rgb: np.ndarray, limit: int = 60_000) -> np.ndarray:
    """Select colored panel pixels for automatic colormap identification."""
    pixels = rgb.reshape(-1, 3)
    saturation = np.max(pixels, axis=1) - np.min(pixels, axis=1)
    brightness = np.mean(pixels, axis=1)
    keep = (
        np.isfinite(pixels).all(axis=1)
        & (saturation > 0.04)
        & (brightness > 0.01)
        & (brightness < 0.99)
    )
    pixels = pixels[keep]
    if pixels.size == 0:
        pixels = rgb.reshape(-1, 3)
    if len(pixels) > limit:
        step = max(1, len(pixels) // limit)
        pixels = pixels[::step][:limit]
    return pixels


def detect_colormap(
    rgb: np.ndarray,
    candidates: Iterable[str] = COMMON_COLORMAPS,
    palette_samples: int = 2048,
) -> tuple[str, dict[str, float]]:
    """Return the closest common Matplotlib colormap and diagnostic scores."""
    pixels = choose_representative_pixels(rgb)
    levels = np.linspace(0.0, 1.0, palette_samples)
    scores: dict[str, float] = {}

    for name in candidates:
        palette = np.asarray(colormaps[name](levels))[:, :3]
        distances, _ = cKDTree(palette).query(pixels, k=1)
        distances = np.asarray(distances, dtype=float)
        # Use both typical and high-percentile mismatch so a few artifacts do
        # not dominate while an incorrect map is still strongly penalized.
        score = float(np.median(distances) + 0.35 * np.quantile(distances, 0.90))
        scores[name] = score

    best = min(scores, key=scores.get)
    return best, scores


def invert_colormap(
    rgb: np.ndarray,
    cmap_name: str,
    samples: int = 4096,
) -> tuple[np.ndarray, np.ndarray]:
    """Recover a normalized scalar image using nearest colormap color."""
    if cmap_name not in colormaps:
        raise ValueError(f"Unknown Matplotlib colormap: {cmap_name}")

    levels = np.linspace(0.0, 1.0, samples)
    palette = np.asarray(colormaps[cmap_name](levels))[:, :3]
    tree = cKDTree(palette)
    distances, indices = tree.query(rgb.reshape(-1, 3), k=1)
    intensity = levels[np.asarray(indices, dtype=int)].reshape(rgb.shape[:2])
    error = np.asarray(distances, dtype=float).reshape(rgb.shape[:2])
    return intensity.astype(np.float32), error.astype(np.float32)


def validate_axes(qr_min: float, qr_max: float, qz_min: float, qz_max: float) -> None:
    values = np.asarray([qr_min, qr_max, qz_min, qz_max], dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("All q-axis limits must be finite numbers.")
    if qr_min >= qr_max:
        raise ValueError("qr minimum must be smaller than qr maximum.")
    if qz_min >= qz_max:
        raise ValueError("qz minimum must be smaller than qz maximum.")


def convert_one(
    png_path: Path,
    qr_min: float,
    qr_max: float,
    qz_min: float,
    qz_max: float,
    auto_crop: bool,
    auto_colormap: bool,
    selected_colormap: str,
) -> tuple[Path, dict]:
    """Convert one PNG and return the output path plus metadata."""
    png_path = png_path.expanduser().resolve()
    if not png_path.exists():
        raise FileNotFoundError(f"File not found: {png_path}")
    if png_path.suffix.lower() != ".png":
        raise ValueError(f"Not a PNG file: {png_path.name}")

    with Image.open(png_path) as image:
        rgb = normalize_rgb(np.asarray(image))

    crop = detect_colored_panel(rgb) if auto_crop else (0, 0, rgb.shape[1], rgb.shape[0])
    x0, y0, x1, y1 = crop
    panel = rgb[y0:y1, x0:x1]
    if panel.shape[0] < 4 or panel.shape[1] < 4:
        raise ValueError(f"Detected panel is too small for {png_path.name}: {crop}")

    if auto_colormap:
        cmap_name, cmap_scores = detect_colormap(panel)
    else:
        cmap_name = selected_colormap
        cmap_scores = {}

    intensity, color_error = invert_colormap(panel, cmap_name)
    height, width = intensity.shape
    qr = np.linspace(qr_min, qr_max, width, dtype=np.float64)
    # Top image row corresponds to the maximum displayed qz.
    qz = np.linspace(qz_max, qz_min, height, dtype=np.float64)

    # Reject pixels that do not reasonably match the selected colormap.
    mask = (
        np.isfinite(intensity)
        & np.isfinite(color_error)
        & (color_error <= 0.10)
    )

    output = png_path.with_name(png_path.stem + "_qspace_axes.npz")
    np.savez_compressed(
        output,
        intensity=intensity,
        qr=qr,
        qz=qz,
        mask=mask,
        color_error=color_error,
        source_png=np.array(str(png_path)),
        crop_xyxy=np.asarray(crop, dtype=np.int32),
        colormap=np.array(cmap_name),
        qr_limits=np.asarray([qr_min, qr_max], dtype=np.float64),
        qz_limits=np.asarray([qz_min, qz_max], dtype=np.float64),
        conversion_note=np.array(
            "Intensity reconstructed from rendered PNG colors; axes generated "
            "from the shared batch limits entered in the converter GUI."
        ),
    )

    # Verify the exact structure expected by the indexing application.
    with np.load(output, allow_pickle=False) as saved:
        if saved["intensity"].shape != (saved["qz"].size, saved["qr"].size):
            raise RuntimeError("Saved NPZ failed the intensity/axis shape check.")
        if saved["mask"].shape != saved["intensity"].shape:
            raise RuntimeError("Saved NPZ failed the mask shape check.")

    metadata = {
        "crop": crop,
        "colormap": cmap_name,
        "valid_percent": float(mask.mean() * 100.0),
        "shape": intensity.shape,
        "scores": cmap_scores,
    }
    return output, metadata


class ConverterApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("780x650")
        self.minsize(700, 560)

        settings = load_settings()
        self.selected_files: list[Path] = []

        self.qr_min_var = tk.StringVar(value=str(settings["qr_min"]))
        self.qr_max_var = tk.StringVar(value=str(settings["qr_max"]))
        self.qz_min_var = tk.StringVar(value=str(settings["qz_min"]))
        self.qz_max_var = tk.StringVar(value=str(settings["qz_max"]))
        self.auto_crop_var = tk.BooleanVar(value=bool(settings["auto_crop"]))
        self.auto_colormap_var = tk.BooleanVar(value=bool(settings["auto_colormap"]))
        self.colormap_var = tk.StringVar(value=str(settings["colormap"]))
        self.preview_first_var = tk.BooleanVar(value=bool(settings["preview_first"]))
        self.file_count_var = tk.StringVar(value="No PNG files selected")
        self.status_var = tk.StringVar(value="Select one or more PNG files to begin.")

        self._build_ui()

    def _build_ui(self) -> None:
        outer = ttk.Frame(self, padding=14)
        outer.pack(fill="both", expand=True)

        title = ttk.Label(
            outer,
            text="Easy batch conversion: select files, enter axes once, convert all",
            font=("TkDefaultFont", 12, "bold"),
        )
        title.pack(anchor="w", pady=(0, 10))

        file_frame = ttk.LabelFrame(outer, text="1. Select PNG files", padding=10)
        file_frame.pack(fill="x", pady=(0, 10))

        buttons = ttk.Frame(file_frame)
        buttons.pack(fill="x")
        ttk.Button(buttons, text="Select PNG files…", command=self.select_files).pack(
            side="left", padx=(0, 8)
        )
        ttk.Button(buttons, text="Select an entire folder…", command=self.select_folder).pack(
            side="left"
        )
        ttk.Label(file_frame, textvariable=self.file_count_var).pack(anchor="w", pady=(8, 0))

        axes_frame = ttk.LabelFrame(
            outer,
            text="2. Reciprocal-space limits — entered once for every selected file",
            padding=10,
        )
        axes_frame.pack(fill="x", pady=(0, 10))

        grid = ttk.Frame(axes_frame)
        grid.pack(anchor="w")
        entries = (
            ("qr minimum (Å⁻¹)", self.qr_min_var),
            ("qr maximum (Å⁻¹)", self.qr_max_var),
            ("qz minimum (Å⁻¹)", self.qz_min_var),
            ("qz maximum (Å⁻¹)", self.qz_max_var),
        )
        for row, (label_text, variable) in enumerate(entries):
            ttk.Label(grid, text=label_text + ":").grid(
                row=row // 2,
                column=(row % 2) * 2,
                sticky="e",
                padx=(0 if row % 2 == 0 else 20, 6),
                pady=4,
            )
            ttk.Entry(grid, textvariable=variable, width=12).grid(
                row=row // 2,
                column=(row % 2) * 2 + 1,
                sticky="w",
                pady=4,
            )

        note = ttk.Label(
            axes_frame,
            text=(
                "These values are remembered after conversion. For a series with the same axes, "
                "you normally enter them only once."
            ),
            wraplength=720,
        )
        note.pack(anchor="w", pady=(6, 0))

        options = ttk.LabelFrame(outer, text="3. Automatic options", padding=10)
        options.pack(fill="x", pady=(0, 10))

        ttk.Checkbutton(
            options,
            text="Automatically detect and crop the colored q-space panel",
            variable=self.auto_crop_var,
        ).pack(anchor="w")
        ttk.Checkbutton(
            options,
            text="Automatically identify the PNG colormap",
            variable=self.auto_colormap_var,
            command=self._update_colormap_state,
        ).pack(anchor="w", pady=(4, 0))

        cmap_row = ttk.Frame(options)
        cmap_row.pack(anchor="w", pady=(6, 0))
        ttk.Label(cmap_row, text="Manual colormap:").pack(side="left", padx=(20, 6))
        self.cmap_combo = ttk.Combobox(
            cmap_row,
            textvariable=self.colormap_var,
            values=COMMON_COLORMAPS,
            state="readonly",
            width=14,
        )
        self.cmap_combo.pack(side="left")

        ttk.Checkbutton(
            options,
            text="Preview the first detected crop before converting",
            variable=self.preview_first_var,
        ).pack(anchor="w", pady=(6, 0))
        self._update_colormap_state()

        action_frame = ttk.Frame(outer)
        action_frame.pack(fill="x", pady=(2, 10))
        self.convert_button = ttk.Button(
            action_frame,
            text="Convert all selected PNG files",
            command=self.convert_selected,
        )
        self.convert_button.pack(side="left")
        ttk.Label(action_frame, textvariable=self.status_var).pack(
            side="left", padx=(12, 0)
        )

        log_frame = ttk.LabelFrame(outer, text="Conversion log", padding=6)
        log_frame.pack(fill="both", expand=True)
        self.log = tk.Text(log_frame, height=14, wrap="word")
        scrollbar = ttk.Scrollbar(log_frame, orient="vertical", command=self.log.yview)
        self.log.configure(yscrollcommand=scrollbar.set)
        self.log.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

    def _update_colormap_state(self) -> None:
        self.cmap_combo.configure(
            state="disabled" if self.auto_colormap_var.get() else "readonly"
        )

    def _append_log(self, message: str) -> None:
        self.log.insert("end", message.rstrip() + "\n")
        self.log.see("end")
        self.update_idletasks()

    def select_files(self) -> None:
        names = filedialog.askopenfilenames(
            title="Select GIWAXS PNG files",
            filetypes=[("PNG images", "*.png"), ("All files", "*.*")],
        )
        if names:
            self.selected_files = [Path(name) for name in names]
            self.file_count_var.set(f"{len(self.selected_files)} PNG file(s) selected")
            self.status_var.set("Ready to convert.")

    def select_folder(self) -> None:
        folder = filedialog.askdirectory(title="Select a folder containing PNG files")
        if folder:
            files = sorted(Path(folder).glob("*.png"))
            self.selected_files = files
            self.file_count_var.set(f"{len(files)} PNG file(s) selected from folder")
            self.status_var.set("Ready to convert." if files else "No PNG files found.")

    def _read_axis_values(self) -> tuple[float, float, float, float]:
        try:
            values = (
                float(self.qr_min_var.get()),
                float(self.qr_max_var.get()),
                float(self.qz_min_var.get()),
                float(self.qz_max_var.get()),
            )
        except ValueError as exc:
            raise ValueError("The four q-axis limits must be valid numbers.") from exc
        validate_axes(*values)
        return values

    def _preview(self, path: Path) -> bool:
        with Image.open(path) as image:
            rgb = normalize_rgb(np.asarray(image))
        crop = detect_colored_panel(rgb) if self.auto_crop_var.get() else (
            0,
            0,
            rgb.shape[1],
            rgb.shape[0],
        )

        window = tk.Toplevel(self)
        window.title("Preview first detected crop")
        window.geometry("1000x520")
        decision = {"accept": False}

        figure = Figure(figsize=(10, 4.5), dpi=100)
        left = figure.add_subplot(1, 2, 1)
        right = figure.add_subplot(1, 2, 2)
        left.imshow(rgb, origin="upper")
        x0, y0, x1, y1 = crop
        from matplotlib.patches import Rectangle

        left.add_patch(
            Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                fill=False,
                linewidth=2.0,
                edgecolor="red",
            )
        )
        left.set_title("Detected panel")
        right.imshow(rgb[y0:y1, x0:x1], origin="upper")
        right.set_title("Converted region")
        left.set_axis_off()
        right.set_axis_off()
        figure.tight_layout()

        canvas = FigureCanvasTkAgg(figure, master=window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        row = ttk.Frame(window, padding=8)
        row.pack(fill="x")

        def accept() -> None:
            decision["accept"] = True
            window.destroy()

        ttk.Button(row, text="Use this crop", command=accept).pack(side="left")
        ttk.Button(row, text="Cancel conversion", command=window.destroy).pack(
            side="left", padx=(8, 0)
        )
        window.transient(self)
        window.grab_set()
        self.wait_window(window)
        return bool(decision["accept"])

    def convert_selected(self) -> None:
        if not self.selected_files:
            messagebox.showwarning(APP_TITLE, "Select at least one PNG file first.")
            return

        try:
            qr_min, qr_max, qz_min, qz_max = self._read_axis_values()
        except Exception as exc:
            messagebox.showerror(APP_TITLE, str(exc))
            return

        if self.preview_first_var.get():
            try:
                if not self._preview(self.selected_files[0]):
                    self.status_var.set("Conversion cancelled.")
                    return
            except Exception as exc:
                messagebox.showerror(APP_TITLE, f"Could not preview the first image:\n{exc}")
                return

        settings = {
            "qr_min": qr_min,
            "qr_max": qr_max,
            "qz_min": qz_min,
            "qz_max": qz_max,
            "auto_crop": self.auto_crop_var.get(),
            "auto_colormap": self.auto_colormap_var.get(),
            "colormap": self.colormap_var.get(),
            "preview_first": self.preview_first_var.get(),
        }
        save_settings(settings)

        self.convert_button.configure(state="disabled")
        self.log.delete("1.0", "end")
        successes = 0
        failures: list[str] = []

        try:
            for index, path in enumerate(self.selected_files, start=1):
                self.status_var.set(f"Converting {index} of {len(self.selected_files)}…")
                self._append_log(f"[{index}/{len(self.selected_files)}] {path.name}")
                try:
                    output, metadata = convert_one(
                        png_path=path,
                        qr_min=qr_min,
                        qr_max=qr_max,
                        qz_min=qz_min,
                        qz_max=qz_max,
                        auto_crop=self.auto_crop_var.get(),
                        auto_colormap=self.auto_colormap_var.get(),
                        selected_colormap=self.colormap_var.get(),
                    )
                    successes += 1
                    self._append_log(
                        "  ✓ "
                        + output.name
                        + f" | crop={metadata['crop']}"
                        + f" | cmap={metadata['colormap']}"
                        + f" | valid={metadata['valid_percent']:.1f}%"
                    )
                except Exception as exc:
                    failures.append(f"{path.name}: {exc}")
                    self._append_log(f"  ✗ Error: {exc}")
        finally:
            self.convert_button.configure(state="normal")

        self.status_var.set(
            f"Finished: {successes} converted, {len(failures)} failed."
        )
        if failures:
            messagebox.showwarning(
                APP_TITLE,
                f"Converted {successes} file(s).\n\n"
                f"{len(failures)} file(s) failed. See the conversion log for details.",
            )
        else:
            messagebox.showinfo(
                APP_TITLE,
                f"Successfully converted {successes} PNG file(s).\n\n"
                "Each output ends with _qspace_axes.npz and contains explicit "
                "intensity, qr, qz, and mask arrays.",
            )


def _running_inside_ipython() -> bool:
    """Return True when the file is being executed by IPython/Jupyter."""
    try:
        return get_ipython() is not None
    except NameError:
        return False


def main():
    """Launch the converter in scripts or Jupyter notebooks.

    Normal Python / PyCharm
    -----------------------
    Tkinter owns the event loop, so ``mainloop()`` is used normally.

    Jupyter / IPython
    -----------------
    IPython is asked to integrate Tk's event loop with the notebook.  This lets
    the GUI stay responsive without permanently blocking the notebook cell.
    The returned ``app`` object is intentionally kept alive by the global below.
    """
    app = ConverterApp()

    if _running_inside_ipython():
        shell = get_ipython()
        try:
            # Equivalent to running `%gui tk` in a notebook cell.
            shell.run_line_magic("gui", "tk")
        except Exception as exc:
            print(
                "Jupyter could not enable Tk GUI integration automatically.\n"
                "If the window does not appear, run `%gui tk` in a cell and then "
                "run this file again.\n"
                f"Details: {exc}"
            )

        app.update_idletasks()
        print(
            "GIWAXS PNG -> NPZ converter launched. "
            "The converter should appear as a separate desktop window."
        )
        return app

    app.mainloop()
    return app


# Keep a global reference so Jupyter does not garbage-collect the GUI window.
_CONVERTER_APP = None


if __name__ == "__main__":
    _CONVERTER_APP = main()
