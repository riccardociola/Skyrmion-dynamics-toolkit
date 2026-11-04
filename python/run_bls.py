"""Run isotropic gamma-point BLS, calculating Wolfram modes only when needed.

From the repository root, for an existing manual export:
    python python/run_bls.py --modes results/modes/bloch_gamma.json

To automatically compute/cache modes for a saved equilibrium:
    python python/run_bls.py --equilibrium NAME.txt

Other options: --force, --use-cache, --timeout-seconds 900.
"""

import argparse
from pathlib import Path
import shutil
import subprocess
import time

from bls import load_magnetic_modes, first_order_intensities

ROOT = Path(__file__).resolve().parents[1]
EQUILIBRIA = ROOT / "results" / "equilibria"
MODE_EXPORTS = ROOT / "results" / "modes"
WOLFRAM_SCRIPT = ROOT / "wolfram" / "export_bls_modes.wls"

# Only changes affecting this export invalidate it. Other Wolfram modules,
# such as plotting or MCA dynamics, are deliberately excluded.
DYNAMICS_SOURCES = (
    "Geometry.wl",
    "FourierBasis.wl",
    "GroundState.wl",
    "DynamicsProjection.wl",
    "DynamicsMagnetic.wl",
    "ExportMagneticModes.wl",
)


def resolve_path(path_text: str, folder: Path) -> Path:
    path = Path(path_text).expanduser()
    if path.is_absolute():
        return path.resolve()
    # Allow either a bare filename (in the expected folder) or a
    # repository-relative path such as results/modes/bloch_gamma.json.
    relative_to_root = ROOT / path
    if relative_to_root.is_file():
        return relative_to_root.resolve()
    return (folder / path).resolve()


def export_is_current(equilibrium: Path, output: Path, *, use_cache: bool) -> bool:
    if not output.is_file():
        return False
    if use_cache:
        return True
    sources = [equilibrium, WOLFRAM_SCRIPT]
    sources.extend(ROOT / "wolfram" / source for source in DYNAMICS_SOURCES)
    sources = [path for path in sources if path.is_file()]
    return output.stat().st_mtime >= max(path.stat().st_mtime for path in sources)


def make_export_if_needed(
    equilibrium: Path,
    output: Path,
    *,
    force: bool = False,
    use_cache: bool = False,
    timeout_seconds: float | None = 900,
) -> None:
    if not equilibrium.is_file():
        raise FileNotFoundError(f"Equilibrium file not found: {equilibrium}")
    if not WOLFRAM_SCRIPT.is_file():
        raise FileNotFoundError(f"Missing Wolfram script: {WOLFRAM_SCRIPT}")

    if not force and export_is_current(equilibrium, output, use_cache=use_cache):
        print(f"Using previously exported modes: {output}")
        return

    executable = shutil.which("wolframscript")
    if executable is None:
        raise RuntimeError(
            "WolframScript is not on PATH. Run "
            "`wolframscript -local -code '2+2'` to test your installation."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    progress_file = Path(str(output) + ".progress.txt")
    console_file = Path(str(output) + ".wolfram.log")
    progress_file.unlink(missing_ok=True)

    # Put -local before -file so neither file argument can be mistaken
    # for the optional Wolfram kernel path accepted by -local.
    command = [
        executable, "-local", "-file", str(WOLFRAM_SCRIPT),
        str(equilibrium), str(output),
    ]
    print("Starting local Wolfram kernel...", flush=True)
    print(f"Wolfram log: {console_file}", flush=True)
    print(f"Progress file: {progress_file}", flush=True)

    started = time.monotonic()
    last_progress = None
    last_status_time = time.monotonic()
    current_stage = "Kernel startup"
    with console_file.open("w", encoding="utf-8") as logfile:
        process = subprocess.Popen(
            command, cwd=ROOT, stdout=logfile, stderr=subprocess.STDOUT,
        )
        try:
            while process.poll() is None:
                if progress_file.is_file():
                    # Mathematica writes this file after each stage, independently
                    # of WolframScript's stdout buffering.
                    line = progress_file.read_text(encoding="utf-8", errors="replace").strip()
                    if line and line != last_progress:
                        print(f"[Wolfram] {line}", flush=True)
                        last_progress = line
                        current_stage = line.split(" | ", 1)[-1]
                        last_status_time = time.monotonic()
                if time.monotonic() - last_status_time > 60:
                    print(f"[Wolfram] Still at: {current_stage} (see {console_file})", flush=True)
                    last_status_time = time.monotonic()
                if timeout_seconds is not None and (time.monotonic() - started) > timeout_seconds:
                    raise TimeoutError(f"Wolfram exceeded {timeout_seconds:g} seconds")
                time.sleep(0.5)
        except (KeyboardInterrupt, TimeoutError):
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            raise
        return_code = process.wait()

    if progress_file.is_file():
        line = progress_file.read_text(encoding="utf-8", errors="replace").strip()
        if line and line != last_progress:
            print(f"[Wolfram] {line}", flush=True)
    if return_code != 0:
        tail = console_file.read_text(encoding="utf-8", errors="replace")[-5000:]
        raise RuntimeError(
            f"WolframScript exited with code {return_code}. "
            f"Full log: {console_file}\nLast log output:\n{tail}"
        )
    if not output.is_file():
        raise RuntimeError(
            f"Wolfram finished, but output is missing: {output}\n"
            f"See: {console_file}"
        )
    print(f"Modes exported: {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--equilibrium", help="Saved isotropic equilibrium filename or path."
    )
    source.add_argument(
        "--modes", help="Use an existing JSON mode export (no Mathematica launched)."
    )
    parser.add_argument("--force", action="store_true", help="Recompute even if cache is current")
    parser.add_argument("--use-cache", action="store_true", help="Use existing export regardless of source timestamps")
    parser.add_argument("--timeout-seconds", type=float, default=900, help="Wolfram time limit (default: 900 seconds)")
    args = parser.parse_args()

    if args.modes:
        output = resolve_path(args.modes, MODE_EXPORTS)
    else:
        equilibrium = resolve_path(args.equilibrium, EQUILIBRIA)
        output = MODE_EXPORTS / f"{equilibrium.stem}_gamma.json"
        make_export_if_needed(
            equilibrium, output, force=args.force, use_cache=args.use_cache,
            timeout_seconds=args.timeout_seconds,
        )

    if not output.is_file():
        raise FileNotFoundError(f"Exported modes not found: {output}")
    data = load_magnetic_modes(output)
    weights = first_order_intensities(
        data, polarization_in=[0, 1, 0], polarization_out=[1, 0, 0]
    )
    print("\nFirst-order BLS weights at k=0:")
    for i, (frequency, weight) in enumerate(zip(data.frequencies, weights)):
        print(f"Mode {i}: frequency={frequency:.5f}, intensity={weight:.8g}")
    print("\nBdG norms:", data.bdg_norms)
    print("Legacy norms:", data.legacy_norms)


if __name__ == "__main__":
    main()