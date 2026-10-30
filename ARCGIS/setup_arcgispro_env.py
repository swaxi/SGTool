"""Create an ArcGIS Pro Python environment for SGTool, Tomofast-x and GeoSubSampler.

The same script is in the ARCGIS folder of all three repositories
(SGTool, tomofast_x_q and GeoSubSampler); running it once sets up one
environment that all three toolboxes can use.

What it does:
  1. Finds ArcGIS Pro and its conda.
  2. Clones the default environment (arcgispro-py3, which cannot be changed)
     into %LOCALAPPDATA%\\ESRI\\conda\\envs\\<name>.
  3. Installs the packages that are missing, from Esri's channels (esri,
     defaults), so arcpy's own packages are not replaced.
  4. Installs pyvista and triangle, which are not on those channels, with pip,
     with every package already in the environment held at its version.
  5. Checks every package imports, then makes the clone the active
     environment for ArcGIS Pro (restart Pro to use it).

Usually run from ArcGIS Pro with the "Set Up Python Environment" tool in each
toolbox. It can also be run from a Windows Command Prompt with ArcGIS Pro's Python:
    "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" setup_arcgispro_env.py
Options:
    --name NAME     environment name (default arcgispro-py3-geotools)
    --no-activate   do not make the environment active in ArcGIS Pro
    --recreate      delete the environment first if it already exists
Running it again on an existing environment just adds what is missing.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

DEFAULT_ENV = "arcgispro-py3-geotools"
CHANNELS = ["esri", "defaults"]

# (conda package, import name, tools that need it)
CONDA_PACKAGES = [
    ("scipy", "scipy", "all"),
    ("matplotlib", "matplotlib", "all"),
    ("pandas", "pandas", "Tomofast-x, GeoSubSampler"),
    ("scikit-learn", "sklearn", "SGTool, GeoSubSampler"),
    ("scikit-image", "skimage", "SGTool"),
    ("pywavelets", "pywt", "SGTool (WTMM)"),
    ("networkx", "networkx", "SGTool, GeoSubSampler"),
    ("shapely", "shapely", "SGTool, GeoSubSampler"),
    ("geopandas", "geopandas", "SGTool, GeoSubSampler"),
    ("pyproj", "pyproj", "SGTool, Tomofast-x"),
    ("fiona", "fiona", "GeoSubSampler"),
    ("defusedxml", "defusedxml", "SGTool"),
]
PIP_PACKAGES = [
    ("pyvista", "pyvista", "Tomofast-x (Visualise)"),
    ("triangle", "triangle", "GeoSubSampler (dyke triangulation)"),
]

# no console window flashing up for each command when run from ArcGIS Pro
_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)


class SetupError(Exception):
    pass


def find_pro():
    """ArcGIS Pro install folder, from the registry or from this Python's path."""
    try:
        import winreg
        for hive in (winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER):
            try:
                with winreg.OpenKey(hive, r"SOFTWARE\ESRI\ArcGISPro") as key:
                    folder = winreg.QueryValueEx(key, "InstallDir")[0]
                if os.path.isdir(folder):
                    return folder
            except OSError:
                pass
    except ImportError:
        pass
    # ...\ArcGIS\Pro\bin\ArcGISPro.exe or ...\ArcGIS\Pro\bin\Python\envs\<env>\python.exe
    parts = os.path.normpath(sys.executable).split(os.sep)
    lower = [p.lower() for p in parts]
    if "bin" in lower:
        return os.sep.join(parts[:lower.index("bin")])
    return None


def _child_env(prefix):
    """Environment variables for the programs this script runs on the new
    environment: run from a toolbox, it inherits those of the ArcGIS Pro
    environment in use (PYTHONHOME, PYTHONPATH, its folders on PATH ...), which
    would make the new environment's python.exe load the other environment's
    packages, so pip installed into it and the check looked in it."""
    env = dict(os.environ)
    for key in ("PYTHONHOME", "PYTHONPATH", "PYTHONSTARTUP", "PYTHONUSERBASE",
                "PYTHONEXECUTABLE", "VIRTUAL_ENV", "CONDA_PREFIX", "CONDA_DEFAULT_ENV",
                "CONDA_PYTHON_EXE", "CONDA_SHLVL", "CONDA_PROMPT_MODIFIER"):
        env.pop(key, None)
    env["PYTHONNOUSERSITE"] = "1"  # nor packages from the user's own site folder
    running = os.path.normcase(os.path.normpath(sys.prefix))
    keep = [d for d in env.get("PATH", "").split(os.pathsep)
            if d and not os.path.normcase(os.path.normpath(d)).startswith(running)]
    first = [prefix, os.path.join(prefix, "Library", "bin"), os.path.join(prefix, "Scripts")]
    env["PATH"] = os.pathsep.join(first + keep)
    return env


class _Setup:
    def __init__(self, log, env=None):
        self.log = log
        self.env = env  # environment variables for the programs run

    def run(self, cmd):
        """Run cmd, passing its output to the log line by line; returns the exit code."""
        self.log("> " + subprocess.list2cmdline(cmd))
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, text=True, errors="replace",
                                creationflags=_NO_WINDOW, env=self.env)
        for line in proc.stdout:
            if line.strip():
                self.log(line.rstrip())
        return proc.wait()

    def missing_modules(self, python, packages):
        """The (package, module, used_by) entries whose module cannot be imported."""
        code = ("import importlib, sys\n"
                "for m in sys.argv[1:]:\n"
                "    try:\n"
                "        importlib.import_module(m)\n"
                "    except Exception:\n"
                "        print(m)\n")
        result = subprocess.run([python, "-c", code] + [m for _, m, _ in packages],
                                capture_output=True, text=True, creationflags=_NO_WINDOW,
                                env=self.env)
        bad = set(result.stdout.split())
        return [p for p in packages if p[1] in bad]

    def clone(self, conda, source, dest):
        self.log(f"Cloning {source}\n     to {dest}\n(this can take several minutes)")
        cmd = [conda, "create", "--clone", source, "--prefix", dest, "--yes", "--quiet"]
        if self.run(cmd + ["--pinned"]) == 0:
            return
        # older/newer conda builds without Esri's --pinned flag
        if os.path.isdir(dest):
            shutil.rmtree(dest, ignore_errors=True)
        if self.run(cmd) != 0:
            raise SetupError("cloning the ArcGIS Pro environment failed (see the conda messages above).")

    def conda_install(self, conda, dest, packages):
        if not packages:
            return []
        base = [conda, "install", "--prefix", dest, "--yes", "--quiet", "--override-channels"]
        for channel in CHANNELS:
            base += ["-c", channel]
        self.log("Installing with conda: " + ", ".join(p for p, _, _ in packages))
        if self.run(base + [p for p, _, _ in packages]) == 0:
            return []
        self.log("Installing them together failed; trying one at a time.")
        return [p for p in packages if self.run(base + [p[0]]) != 0]

    def pip_install(self, python, packages):
        if not packages:
            return []
        names = {p for p, _, _ in packages}
        frozen = subprocess.run([python, "-m", "pip", "list", "--format=freeze"],
                                capture_output=True, text=True, creationflags=_NO_WINDOW,
                                env=self.env).stdout
        pins = [line for line in frozen.splitlines()
                if "==" in line and line.split("==")[0].lower() not in names]
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write("\n".join(pins))
            constraints = f.name
        try:
            self.log("Installing with pip: " + ", ".join(sorted(names)))
            cmd = [python, "-m", "pip", "install", "--progress-bar", "off",
                   "--constraint", constraints]
            if self.run(cmd + [p for p, _, _ in packages]) == 0:
                return []
            return [p for p in packages if self.run(cmd + [p[0]]) != 0]
        finally:
            os.remove(constraints)


def setup(name=DEFAULT_ENV, activate=True, recreate=False, log=print):
    """Create or update the environment. Returns the list of packages still missing;
    raises SetupError if the environment could not be made at all."""
    if os.name != "nt":
        raise SetupError("ArcGIS Pro only runs on Windows.")
    pro = find_pro()
    if not pro:
        raise SetupError("ArcGIS Pro was not found.")
    scripts = os.path.join(pro, "bin", "Python", "Scripts")
    conda = os.path.join(scripts, "conda.exe")
    source = os.path.join(pro, "bin", "Python", "envs", "arcgispro-py3")
    if not os.path.isfile(conda) or not os.path.isdir(source):
        raise SetupError(f"conda or the arcgispro-py3 environment is missing under {pro}")

    envs = os.path.join(os.environ["LOCALAPPDATA"], "ESRI", "conda", "envs")
    dest = os.path.join(envs, name)
    python = os.path.join(dest, "python.exe")
    log(f"ArcGIS Pro:  {pro}\nEnvironment: {dest}")

    if os.path.normcase(os.path.normpath(sys.prefix)) == os.path.normcase(dest):
        raise SetupError(
            f"'{name}' is the environment ArcGIS Pro is using now, and packages cannot be added to "
            "it while it is in use. In Project > Package Manager switch back to arcgispro-py3, "
            "restart ArcGIS Pro and run this again (or choose a different environment name).")
    worker = _Setup(log, _child_env(dest))
    if recreate and os.path.isdir(dest):
        log("Removing the existing environment")
        if worker.run([conda, "env", "remove", "--prefix", dest, "--yes", "--quiet"]) != 0:
            raise SetupError("could not remove the existing environment.")
    if os.path.isfile(python):
        log("The environment already exists; adding anything that is missing.")
    else:
        os.makedirs(envs, exist_ok=True)
        worker.clone(conda, source, dest)

    failed = worker.conda_install(conda, dest, worker.missing_modules(python, CONDA_PACKAGES))
    # anything conda could not supply, try pip
    worker.pip_install(python, failed + worker.missing_modules(python, PIP_PACKAGES))

    log("Checking the packages:")
    still_missing = worker.missing_modules(python, CONDA_PACKAGES + PIP_PACKAGES)
    for entry in CONDA_PACKAGES + PIP_PACKAGES:
        log(f"  {'MISSING' if entry in still_missing else 'ok     '}  {entry[0]:<14} {entry[2]}")

    if activate:
        log("Making it the active ArcGIS Pro environment")
        if worker.run(["cmd", "/c", os.path.join(scripts, "proswap.bat"), "-p", dest]) != 0:
            log(f"proswap failed: in ArcGIS Pro go to Project > Package Manager, "
                f"open the Environment Manager and activate '{name}'.")

    if still_missing:
        log("Finished, but some packages are missing (see above); the tools that need them "
            "will not run. Run this again, or add them in ArcGIS Pro's Package Manager.")
    else:
        log("Done. Restart ArcGIS Pro to use the environment.")
    return [p for p, _, _ in still_missing]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--name", default=DEFAULT_ENV, help="environment name")
    parser.add_argument("--no-activate", action="store_true",
                        help="do not make the environment active in ArcGIS Pro")
    parser.add_argument("--recreate", action="store_true",
                        help="delete the environment first if it already exists")
    args = parser.parse_args()
    try:
        missing = setup(args.name, not args.no_activate, args.recreate,
                        log=lambda msg: print(msg, flush=True))
    except SetupError as exc:
        print(f"\nERROR: {exc}")
        sys.exit(1)
    sys.exit(1 if missing else 0)


if __name__ == "__main__":
    main()
