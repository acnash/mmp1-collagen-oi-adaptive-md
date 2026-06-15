#!/usr/bin/env python3
"""
adaptive_mmp1_unwinding_dual_worker_pbc_corrected.py

Fault-tolerant adaptive OpenMM simulation runner for MMP1-collagen unwinding.

This version corrects the important periodic-boundary-condition issue:
  - The adaptive opening score now uses minimum-image distances from OpenMM Context.getState(enforcePeriodicBox=True)
    with the current periodic box vectors.
  - Selected frames are therefore ranked using PBC-aware inter-chain CA distances.
  - This prevents periodic imaging artefacts from being selected as "unwound" collagen.

It also retains the dual-worker scheduler corrections:
  - The controller does not launch the same worker twice when --concurrent_workers > 1.
  - Active worker directories are excluded from relaunch selection.
  - JSON, NPZ, and checkpoint writes use unique temporary filenames.
  - Worker stdout/stderr logs are overwritten on each launch.
  - On checkpoint resume, existing production.xtc is archived and a fresh trajectory segment is written.

Typical fresh command:
  conda run -n openmm python scripts/adaptive_mmp1_unwinding_dual_worker.py \
    --input_dir generated_mutants/salted_150mM_NaCl/collagen_G978S \
    --system_variant G978S \
    --input_gro NPT_eq_collagen_G978S_150mM_NaCl.gro \
    --platform CUDA \
    --concurrent_workers 4 \
    --mode auto \
    --target_generation 50

Typical resume command:
  conda run -n openmm python scripts/adaptive_mmp1_unwinding_dual_worker.py \
    --resume generated_mutants/salted_150mM_NaCl/collagen_G978S/adaptive_run_G978S_YYYYMMDD_HHMMSS \
    --platform CUDA \
    --concurrent_workers 4 \
    --mode auto \
    --target_generation 50
"""

import argparse
import json
import math
import os
import re
import signal
import subprocess
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, Platform, unit
from openmm.app import (
    GromacsGroFile,
    GromacsTopFile,
    HBonds,
    PME,
    Simulation,
    StateDataReporter,
    XTCReporter,
)


GRO_ATOM_OFFSET = 1
STOP_REQUESTED = False

SYSTEM_VARIANTS = {
    "wild_type": {
        "label": "wild_type",
        "mutated_chain": None,
        "mutation_residue": None,
        "wild_type_residue": None,
        "mutant_residue": None,
    },
    "G978S": {
        "label": "G978S",
        "mutated_chain": "chain_1",
        "mutation_residue": 978,
        "wild_type_residue": "GLY",
        "mutant_residue": "SER",
    },
    "G984C": {
        "label": "G984C",
        "mutated_chain": "chain_1",
        "mutation_residue": 984,
        "wild_type_residue": "GLY",
        "mutant_residue": "CYS",
    },
    "G987R": {
        "label": "G987R",
        "mutated_chain": "chain_1",
        "mutation_residue": 987,
        "wild_type_residue": "GLY",
        "mutant_residue": "ARG",
    },
}

DEFAULT_CONFIG = {
    "input_gro": None,
    "input_top": "system.top",
    "system_variant": "wild_type",
    "run_label": None,
    "collagen_residue_start": 977,
    "collagen_residue_end": 987,
    "mutated_chain": None,
    "temperature_k": 310.15,
    "pressure_bar": 1.0,
    "timestep_ps": 0.002,
    "nonbonded_cutoff_nm": 1.0,
    "nvt_ps": 50.0,
    "npt_equil_ps": 200.0,
    "production_ns": 20.0,
    "trajectory_interval_ps": 10.0,
    "score_interval_ps": 10.0,
    "checkpoint_interval_ps": 200.0,
    "completion_fraction": 0.95,
    "workers_per_generation": 10,
    "top_frames_to_keep": 2,
    "platform": "CPU",
    "minimize_max_iterations": 500,
    "native_retry_limit": 1,
    "catastrophic_failure_fraction": 0.20,
    "slurm_walltime_buffer_min": 20.0,
    "master_seed": 90210,
    "collagen_windows_gro_1based": None,
    "scissile_region": None,
}

TERMINAL_SUCCESS_STATUSES = {"completed", "completed_approximate"}
TERMINAL_FAILURE_STATUSES = {"failed_archived"}
SKIP_STATUSES = TERMINAL_SUCCESS_STATUSES | TERMINAL_FAILURE_STATUSES


def handle_signal(signum, frame):
    global STOP_REQUESTED
    STOP_REQUESTED = True


signal.signal(signal.SIGTERM, handle_signal)
signal.signal(signal.SIGINT, handle_signal)


def now_iso() -> str:
    return datetime.now().isoformat()


def make_timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def write_json(path: Path, data: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    with open(tmp, "w") as handle:
        json.dump(data, handle, indent=2)
    tmp.replace(path)


def read_json(path: Path) -> dict:
    try:
        with open(path, "r") as handle:
            return json.load(handle)
    except Exception:
        return {}


def gro_to_openmm_index(gro_index: int) -> int:
    return int(gro_index) - GRO_ATOM_OFFSET


def steps_from_ps(ps: float, timestep_ps: float) -> int:
    return int(round(float(ps) / float(timestep_ps)))


def steps_from_ns(ns: float, timestep_ps: float) -> int:
    return int(round(float(ns) * 1000.0 / float(timestep_ps)))


def positions_to_numpy_nm(positions):
    return np.array([[p.x, p.y, p.z] for p in positions.value_in_unit(unit.nanometer)], dtype=float)


def box_vectors_to_numpy_nm(box_vectors):
    return np.array([[v.x, v.y, v.z] for v in box_vectors.value_in_unit(unit.nanometer)], dtype=float)


def numpy_to_positions(array_nm):
    return unit.Quantity(np.asarray(array_nm, dtype=float), unit.nanometer)


def numpy_to_box_vectors(array_nm):
    from openmm import Vec3
    array_nm = np.asarray(array_nm, dtype=float)
    return tuple(Vec3(*array_nm[i]) * unit.nanometer for i in range(3))


def save_state_npz(path: Path, positions, box_vectors) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    pos_nm = positions_to_numpy_nm(positions)
    box_nm = box_vectors_to_numpy_nm(box_vectors)

    tmp = path.with_name(f"{path.stem}.tmp.{os.getpid()}.{time.time_ns()}.npz")
    np.savez_compressed(tmp, positions_nm=pos_nm, box_vectors_nm=box_nm)

    if not tmp.exists():
        alt = Path(str(tmp) + ".npz")
        if alt.exists():
            tmp = alt

    tmp.replace(path)


def load_state_npz(path: Path):
    data = np.load(path)
    return numpy_to_positions(data["positions_nm"]), numpy_to_box_vectors(data["box_vectors_nm"])


def minimum_image_displacement_nm(delta_nm: np.ndarray, box_vectors_nm: np.ndarray) -> np.ndarray:
    """
    Return minimum-image displacement for a triclinic periodic box.

    delta_nm is p_i - p_j in nm.
    box_vectors_nm is a 3x3 matrix with OpenMM box vectors as rows.
    """
    box = np.asarray(box_vectors_nm, dtype=float)
    inv_box = np.linalg.inv(box.T)
    fractional = inv_box @ np.asarray(delta_nm, dtype=float)
    fractional -= np.round(fractional)
    return box.T @ fractional


def pbc_distance_nm(p1_nm: np.ndarray, p2_nm: np.ndarray, box_vectors_nm: np.ndarray) -> float:
    disp = minimum_image_displacement_nm(np.asarray(p1_nm) - np.asarray(p2_nm), box_vectors_nm)
    return float(np.linalg.norm(disp))


def write_gro_like_snapshot(path: Path, topology, positions, box_vectors, title: str = "Best frame") -> None:
    pos_nm = positions_to_numpy_nm(positions)
    atoms = list(topology.atoms())

    with open(path, "w") as handle:
        handle.write(f"{title}\n")
        handle.write(f"{len(atoms):5d}\n")

        for atom in atoms:
            resid = atom.residue.id
            try:
                resid_int = int(resid)
            except Exception:
                resid_int = atom.residue.index + 1

            resname = atom.residue.name[:5]
            atomname = atom.name[:5]
            atomid = atom.index + 1
            x, y, z = pos_nm[atom.index]
            handle.write(f"{resid_int % 100000:5d}{resname:<5}{atomname:>5}{atomid % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")

        box_nm = box_vectors_to_numpy_nm(box_vectors)
        handle.write(f"{box_nm[0,0]:10.5f}{box_nm[1,1]:10.5f}{box_nm[2,2]:10.5f}\n")


def enable_periodic_bonded_forces(system):
    periodic_force_classes = {
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "RBTorsionForce",
        "CMAPTorsionForce",
        "CustomBondForce",
        "CustomAngleForce",
        "CustomTorsionForce",
    }

    enabled = []

    for force in system.getForces():
        force_name = force.__class__.__name__
        if force_name in periodic_force_classes and hasattr(force, "setUsesPeriodicBoundaryConditions"):
            force.setUsesPeriodicBoundaryConditions(True)
            enabled.append(force_name)

    return enabled


def auto_detect_gro_file(input_dir: Path) -> str:
    gro_files = sorted(Path(input_dir).glob("*.gro"))
    if not gro_files:
        raise FileNotFoundError(f"No .gro file found in {input_dir}. Use --input_gro to specify one.")
    if len(gro_files) > 1:
        names = ", ".join(p.name for p in gro_files)
        preferred = [p for p in gro_files if "150mM_NaCl" in p.name]
        if len(preferred) == 1:
            return preferred[0].name
        raise RuntimeError(f"Multiple .gro files found in {input_dir}: {names}. Use --input_gro.")
    return gro_files[0].name


def infer_system_variant(input_dir: Path) -> str:
    input_dir = Path(input_dir)
    for variant in ["G978S", "G984C", "G987R"]:
        if (input_dir / f"collagen_{variant}.itp").exists() or variant in input_dir.name:
            return variant
    return "wild_type"


def normalize_system_config(input_dir: Path, config: dict) -> dict:
    input_dir = Path(input_dir)
    config = dict(config)

    variant = config.get("system_variant")
    if variant in {None, "auto"}:
        variant = infer_system_variant(input_dir)
    if variant not in SYSTEM_VARIANTS:
        raise ValueError(f"Unknown system_variant {variant!r}. Choose one of: {', '.join(SYSTEM_VARIANTS)}")

    config["system_variant"] = variant
    variant_meta = SYSTEM_VARIANTS[variant]
    if config.get("mutated_chain") is None:
        config["mutated_chain"] = variant_meta["mutated_chain"]
    if config.get("run_label") is None:
        config["run_label"] = variant_meta["label"]
    if config.get("input_gro") in {None, "", "auto"}:
        config["input_gro"] = auto_detect_gro_file(input_dir)

    return config


def build_topology_and_system(input_dir: Path, config: dict):
    input_dir = Path(input_dir)
    config = normalize_system_config(input_dir, config)
    gro_path = input_dir / config["input_gro"]
    top_path = input_dir / config["input_top"]

    gro = GromacsGroFile(str(gro_path))

    top = GromacsTopFile(
        str(top_path),
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir=str(input_dir),
    )

    system = top.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=float(config["nonbonded_cutoff_nm"]) * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
    )

    enable_periodic_bonded_forces(system)
    return gro, top, system


def residue_id_as_int(residue) -> Optional[int]:
    try:
        return int(residue.id)
    except Exception:
        try:
            return int(residue.index) + 1
        except Exception:
            return None


def expected_residue_name(config: dict, chain_name: str, resid: int) -> str:
    variant = config.get("system_variant", "wild_type")
    meta = SYSTEM_VARIANTS.get(variant, SYSTEM_VARIANTS["wild_type"])
    mutated_chain = config.get("mutated_chain") or meta.get("mutated_chain")
    if (
        meta.get("mutation_residue") == resid
        and mutated_chain == chain_name
    ):
        return str(meta["mutant_residue"])
    return "GLY" if resid in {978, 981, 984, 987, 990} else ""


def get_ca_indices(topology, config: dict) -> Dict[str, List[int]]:
    residue_start = int(config.get("collagen_residue_start", 977))
    residue_end = int(config.get("collagen_residue_end", 987))
    expected_resids = list(range(residue_start, residue_end + 1))
    expected_len = len(expected_resids)
    if expected_len <= 0:
        raise ValueError("collagen_residue_end must be >= collagen_residue_start")

    ca_indices = {}
    atoms = list(topology.atoms())

    windows = config.get("collagen_windows_gro_1based")
    if windows:
        for chain_name, gro_range in windows.items():
            start = gro_to_openmm_index(gro_range[0])
            end = gro_to_openmm_index(gro_range[1])
            selected = []

            for atom in atoms:
                if start <= atom.index <= end and atom.name == "CA":
                    selected.append(atom.index)

            selected = sorted(selected)

            if len(selected) != expected_len:
                raise ValueError(
                    f"{chain_name} should contain {expected_len} CA atoms for residues "
                    f"{residue_start}-{residue_end}, but found {len(selected)}: {selected}"
                )

            ca_indices[chain_name] = selected

        return ca_indices

    collagen_ca_records: List[Tuple[int, int, str, int]] = []
    for atom in atoms:
        if atom.name != "CA":
            continue
        resid = residue_id_as_int(atom.residue)
        if resid is None:
            continue
        if residue_start <= resid <= residue_end:
            collagen_ca_records.append((atom.index, resid, atom.residue.name, atom.residue.index))

    collagen_ca_records = sorted(collagen_ca_records, key=lambda x: x[0])
    required_total = 3 * expected_len
    if len(collagen_ca_records) != required_total:
        raise ValueError(
            f"Expected {required_total} collagen CA atoms for three chains across residues "
            f"{residue_start}-{residue_end}, but found {len(collagen_ca_records)}: "
            f"{[(idx, resid, name) for idx, resid, name, _ in collagen_ca_records]}"
        )

    validation = {}
    for chain_i in range(3):
        chain_name = f"chain_{chain_i + 1}"
        chunk = collagen_ca_records[chain_i * expected_len:(chain_i + 1) * expected_len]
        selected = [idx for idx, _, _, _ in chunk]
        resids = [resid for _, resid, _, _ in chunk]
        resnames = [name for _, _, name, _ in chunk]

        if resids != expected_resids:
            raise ValueError(f"{chain_name} CA residues should be {expected_resids}, but found {resids}")

        validation[chain_name] = dict(zip(resids, resnames))
        ca_indices[chain_name] = selected

    validate_variant_residue_names(config, validation)

    return ca_indices


def validate_variant_residue_names(config: dict, validation: Dict[str, Dict[int, str]]) -> None:
    variant = config.get("system_variant", "wild_type")
    if variant not in SYSTEM_VARIANTS:
        raise ValueError(f"Unknown system_variant {variant!r}")

    meta = SYSTEM_VARIANTS[variant]
    mutation_residue = meta.get("mutation_residue")
    mutated_chain = config.get("mutated_chain") or meta.get("mutated_chain")

    for chain_name, residue_names in validation.items():
        if mutation_residue is None:
            for resid in [978, 984, 987]:
                observed = residue_names.get(resid)
                if observed != "GLY":
                    raise ValueError(f"wild_type expected {chain_name} residue {resid} to be GLY, found {observed}")
            continue

        observed = residue_names.get(int(mutation_residue))
        expected = meta["mutant_residue"] if chain_name == mutated_chain else meta["wild_type_residue"]
        if observed != expected:
            raise ValueError(
                f"{variant} expected {chain_name} residue {mutation_residue} to be {expected}, found {observed}"
            )


def opening_score_from_positions_nm(positions_nm: np.ndarray, ca_indices: Dict[str, List[int]], box_vectors_nm: np.ndarray) -> float:
    """
    PBC-aware opening score in nm.

    Uses minimum-image distances between corresponding CA atoms across the three collagen chains.
    This is the adaptive selection score, so it must not be fooled by periodic imaging.
    """
    values = []

    for i in range(11):
        p1 = positions_nm[ca_indices["chain_1"][i]]
        p2 = positions_nm[ca_indices["chain_2"][i]]
        p3 = positions_nm[ca_indices["chain_3"][i]]

        values.append(pbc_distance_nm(p1, p2, box_vectors_nm))
        values.append(pbc_distance_nm(p1, p3, box_vectors_nm))
        values.append(pbc_distance_nm(p2, p3, box_vectors_nm))

    return float(np.mean(values))


class OpeningScoreReporter:
    def __init__(self, file_path: Path, report_interval: int, ca_indices: Dict[str, List[int]], best_npz_path: Path, best_gro_path: Path, topology, append: bool = True):
        self.file_path = Path(file_path)
        self.reportInterval = int(report_interval)
        self.ca_indices = ca_indices
        self.best_npz_path = Path(best_npz_path)
        self.best_gro_path = Path(best_gro_path)
        self.topology = topology
        self.best_score = -1.0e30
        self.best_step = None

        if append and self.file_path.exists():
            try:
                old = np.genfromtxt(self.file_path, delimiter=",", names=True)
                if old.size > 0 and "opening_score_angstrom" in old.dtype.names:
                    arr = np.atleast_1d(old)
                    idx = int(np.nanargmax(arr["opening_score_angstrom"]))
                    self.best_score = float(arr["opening_score_angstrom"][idx]) / 10.0
                    self.best_step = int(arr["step"][idx])
            except Exception:
                pass
            self._handle = open(self.file_path, "a")
        else:
            self._handle = open(self.file_path, "w")
            self._handle.write("step,time_ps,opening_score_nm,opening_score_angstrom,metric_uses_pbc\n")
            self._handle.flush()

    def describeNextReport(self, simulation):
        steps = self.reportInterval - simulation.currentStep % self.reportInterval
        return (steps, True, False, False, True, None)

    def report(self, simulation, state):
        positions = state.getPositions()
        box_vectors = state.getPeriodicBoxVectors()
        pos_nm = positions_to_numpy_nm(positions)
        box_nm = box_vectors_to_numpy_nm(box_vectors)
        score_nm = opening_score_from_positions_nm(pos_nm, self.ca_indices, box_nm)
        time_ps = state.getTime().value_in_unit(unit.picosecond)

        self._handle.write(f"{simulation.currentStep},{time_ps:.6f},{score_nm:.8f},{score_nm * 10.0:.8f},True\n")
        self._handle.flush()

        if score_nm > self.best_score:
            self.best_score = score_nm
            self.best_step = simulation.currentStep
            save_state_npz(self.best_npz_path, positions, box_vectors)
            write_gro_like_snapshot(self.best_gro_path, self.topology, positions, box_vectors, title=f"Best PBC-corrected opening frame step {simulation.currentStep}")

    def close(self):
        try:
            self._handle.close()
        except Exception:
            pass

    def __del__(self):
        self.close()


def save_openmm_checkpoint_atomic(simulation: Simulation, checkpoint_path: Path) -> bool:
    checkpoint_path = Path(checkpoint_path)
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = checkpoint_path.with_name(f"{checkpoint_path.name}.tmp.{os.getpid()}.{time.time_ns()}")
    try:
        simulation.saveCheckpoint(str(tmp))
        tmp.replace(checkpoint_path)
        return True
    except Exception:
        try:
            if tmp.exists():
                tmp.unlink()
        except Exception:
            pass
        return False


def create_simulation(topology, system, config: dict, seed: int, add_barostat: bool):
    temperature = float(config["temperature_k"]) * unit.kelvin
    timestep = float(config["timestep_ps"]) * unit.picosecond

    if add_barostat:
        pressure = float(config["pressure_bar"]) * unit.bar
        system.addForce(MonteCarloBarostat(pressure, temperature, 25))

    integrator = LangevinMiddleIntegrator(temperature, 1.0 / unit.picosecond, timestep)
    integrator.setRandomNumberSeed(int(seed))

    platform_name = config.get("platform", "CPU")
    platform = Platform.getPlatformByName(platform_name)

    properties = {}
    if platform_name == "CUDA":
        properties["Precision"] = "mixed"

    return Simulation(topology, system, integrator, platform, properties)


def deterministic_seed(config: dict, generation_index: int, worker_index: int, attempt: int, salt: int = 0) -> int:
    master = int(config.get("master_seed", 90210))
    return int(master + generation_index * 100000 + worker_index * 1000 + attempt * 100 + salt)


def initialise_context(input_dir: str, overrides: dict) -> dict:
    input_dir = Path(input_dir).resolve()
    config = dict(DEFAULT_CONFIG)
    config.update({k: v for k, v in overrides.items() if v is not None})
    config = normalize_system_config(input_dir, config)

    run_label = str(config.get("run_label") or config.get("system_variant") or "system")
    run_label = re.sub(r"[^A-Za-z0-9_.-]+", "_", run_label)
    run_dir = input_dir / f"adaptive_run_{run_label}_{make_timestamp()}"
    run_dir.mkdir(parents=True, exist_ok=False)

    context = {
        "created": now_iso(),
        "input_dir": str(input_dir),
        "run_dir": str(run_dir),
        "config": config,
        "next_generation": 0,
        "completed_generations": [],
        "selected_frames": [],
        "lineage": [],
        "notes": [
            "Backwards-compatible adaptive OpenMM MMP1-collagen unwinding run.",
            "Opening score uses minimum-image PBC-corrected CA distances.",
            "Manual mode runs one generation and stops.",
            "Auto mode runs until target_generation has completed.",
            "Periodic bonded forces are enabled for cyclic/periodic collagen topology.",
            "Production checkpoints are rolling and overwritten to conserve space.",
            f"System variant: {config['system_variant']}",
            f"Input GRO: {config['input_gro']}",
        ],
    }

    write_json(run_dir / "context.json", context)
    return context


def apply_runtime_overrides(context: dict, overrides: dict) -> dict:
    config = context["config"]

    for key, value in overrides.items():
        if value is not None:
            config[key] = value

    config = normalize_system_config(Path(context["input_dir"]), config)
    context["config"] = config

    context["last_runtime_overrides"] = {k: v for k, v in overrides.items() if v is not None}
    context["last_updated"] = now_iso()
    write_json(Path(context["run_dir"]) / "context.json", context)
    return context


def load_context(run_dir: Path) -> dict:
    context_path = Path(run_dir) / "context.json"

    if not context_path.exists():
        raise FileNotFoundError(f"Could not find context.json in {run_dir}")

    context = read_json(context_path)

    if "config" not in context:
        context["config"] = {}

    merged_config = dict(DEFAULT_CONFIG)
    merged_config.update(context["config"])
    context["config"] = normalize_system_config(Path(context.get("input_dir", run_dir.parent)), merged_config)

    if "completed_generations" not in context:
        context["completed_generations"] = []

    if "selected_frames" not in context:
        context["selected_frames"] = []

    if "lineage" not in context:
        context["lineage"] = []

    context["last_backward_compatibility_update"] = now_iso()
    write_json(context_path, context)
    return context


def get_worker_dir(run_dir: Path, generation_index: int, worker_index: int) -> Path:
    return Path(run_dir) / f"generation_{generation_index:03d}" / f"worker_{worker_index:03d}"


def get_worker_json(worker_dir: Path) -> Path:
    return Path(worker_dir) / "worker.json"


def read_worker_meta(worker_dir: Path) -> dict:
    return read_json(get_worker_json(worker_dir))


def write_worker_meta(worker_dir: Path, meta: dict) -> None:
    write_json(get_worker_json(worker_dir), meta)


def count_score_frames(worker_dir: Path) -> int:
    csv_path = worker_dir / "opening_scores.csv"
    if not csv_path.exists():
        return 0
    try:
        with open(csv_path, "r") as handle:
            return max(0, sum(1 for _ in handle) - 1)
    except Exception:
        return 0


def worker_completion_fraction(worker_dir: Path, config: dict) -> float:
    target_steps = steps_from_ns(config["production_ns"], config["timestep_ps"])
    score_interval_steps = steps_from_ps(config["score_interval_ps"], config["timestep_ps"])
    expected_frames = max(1, target_steps // score_interval_steps)
    frames = count_score_frames(worker_dir)
    return min(1.0, float(frames) / float(expected_frames))


def is_worker_valid_complete(worker_dir: Path, config: dict) -> bool:
    meta = read_worker_meta(worker_dir)
    status = meta.get("status", "unknown")

    if status in TERMINAL_SUCCESS_STATUSES:
        return True

    if (worker_dir / "best_opening_frame.npz").exists() and (worker_dir / "opening_scores.csv").exists():
        frac = worker_completion_fraction(worker_dir, config)
        if frac >= float(config.get("completion_fraction", 0.95)):
            meta["status"] = "completed_approximate"
            meta["completion_fraction"] = frac
            meta["completed_approximate_at"] = now_iso()
            write_worker_meta(worker_dir, meta)
            return True

    return False


def archive_failed_worker(worker_dir: Path, reason: str) -> Path:
    worker_dir = Path(worker_dir)
    failure_dir = worker_dir / "failure_archive"
    failure_dir.mkdir(parents=True, exist_ok=True)

    meta = read_worker_meta(worker_dir)
    meta["status"] = "failed_archived"
    meta["failure_reason"] = reason
    meta["failed_archived_at"] = now_iso()
    write_worker_meta(worker_dir, meta)

    with open(failure_dir / "failure_reason.txt", "w") as handle:
        handle.write(reason + "\n")

    return failure_dir


def existing_generation_worker_dirs(run_dir: Path, generation_index: int) -> List[Path]:
    generation_dir = Path(run_dir) / f"generation_{generation_index:03d}"
    if not generation_dir.exists():
        return []
    return sorted([p for p in generation_dir.glob("worker_*") if p.is_dir()])


def get_successful_worker_dirs(run_dir: Path, generation_index: int, config: dict) -> List[Path]:
    out = []
    for worker_dir in existing_generation_worker_dirs(run_dir, generation_index):
        if is_worker_valid_complete(worker_dir, config):
            out.append(worker_dir)
    return sorted(out)


def get_archived_failed_worker_dirs(run_dir: Path, generation_index: int) -> List[Path]:
    out = []
    for worker_dir in existing_generation_worker_dirs(run_dir, generation_index):
        meta = read_worker_meta(worker_dir)
        if meta.get("status") == "failed_archived":
            out.append(worker_dir)
    return sorted(out)


def next_worker_index(run_dir: Path, generation_index: int) -> int:
    max_idx = -1
    for worker_dir in existing_generation_worker_dirs(run_dir, generation_index):
        try:
            idx = int(worker_dir.name.split("_")[-1])
            max_idx = max(max_idx, idx)
        except Exception:
            pass
    return max_idx + 1


def get_start_state_pool(context: dict, generation_index: int) -> List[dict]:
    config = context["config"]

    if generation_index == 0:
        return [{"source": "initial", "state_npz": None, "source_generation": None, "source_worker": None} for _ in range(int(config["workers_per_generation"]))]

    if not context.get("selected_frames"):
        raise RuntimeError(f"No selected frames are available for generation {generation_index}")

    previous_selected = context["selected_frames"][-1]["frames"]
    starts = []

    workers_per_source = int(config["workers_per_generation"]) // len(previous_selected)
    remainder = int(config["workers_per_generation"]) % len(previous_selected)

    for source_i, frame in enumerate(previous_selected):
        n = workers_per_source + (1 if source_i < remainder else 0)
        for _ in range(n):
            starts.append(
                {
                    "source": frame.get("worker_dir"),
                    "state_npz": frame.get("best_state_npz"),
                    "source_generation": context["selected_frames"][-1].get("generation"),
                    "source_worker": frame.get("worker_dir"),
                    "source_rank": frame.get("rank"),
                    "source_best_step": frame.get("best_step"),
                    "source_best_score_angstrom": frame.get("best_score_angstrom"),
                }
            )

    return starts


def choose_start_state_for_new_worker(context: dict, generation_index: int, worker_index: int) -> dict:
    pool = get_start_state_pool(context, generation_index)
    if not pool:
        raise RuntimeError("Empty start-state pool")
    return pool[worker_index % len(pool)]


def prepare_worker_if_needed(context: dict, generation_index: int, active_worker_dirs: Optional[Set[Path]] = None) -> Optional[Path]:
    active_worker_dirs = {Path(p).resolve() for p in (active_worker_dirs or set())}

    run_dir = Path(context["run_dir"])
    config = context["config"]
    target_valid = int(config["workers_per_generation"])

    successful = get_successful_worker_dirs(run_dir, generation_index, config)

    if len(successful) >= target_valid:
        return None

    archived_failures = get_archived_failed_worker_dirs(run_dir, generation_index)
    fail_fraction = len(archived_failures) / max(1, target_valid)

    if fail_fraction >= float(config.get("catastrophic_failure_fraction", 0.20)):
        generation_dir = run_dir / f"generation_{generation_index:03d}"
        generation_meta = read_json(generation_dir / "generation.json")
        generation_meta["status"] = "catastrophic_failure_pause"
        generation_meta["failed_archived_workers"] = [str(p) for p in archived_failures]
        generation_meta["message"] = "Catastrophic failure threshold reached. The run is restartable; inspect failed workers and resume to generate replacements."
        generation_meta["updated"] = now_iso()
        write_json(generation_dir / "generation.json", generation_meta)
        raise RuntimeError("Catastrophic failure threshold reached. Campaign paused but restartable.")

    for worker_dir in existing_generation_worker_dirs(run_dir, generation_index):
        worker_dir_resolved = worker_dir.resolve()

        if worker_dir_resolved in active_worker_dirs:
            continue

        meta = read_worker_meta(worker_dir)
        status = meta.get("status", "unknown")

        if status in SKIP_STATUSES:
            continue

        if is_worker_valid_complete(worker_dir, config):
            continue

        return worker_dir

    if active_worker_dirs:
        return None

    worker_index = next_worker_index(run_dir, generation_index)
    worker_dir = get_worker_dir(run_dir, generation_index, worker_index)
    worker_dir.mkdir(parents=True, exist_ok=True)

    start_state = choose_start_state_for_new_worker(context, generation_index, worker_index)

    meta = {
        "worker_index": worker_index,
        "generation": generation_index,
        "status": "queued",
        "created": now_iso(),
        "attempt": 0,
        "max_retries": int(config.get("native_retry_limit", 1)),
        "start_state": start_state,
        "is_replacement": worker_index >= target_valid,
        "replacement_for_worker": None,
        "lineage": start_state,
    }

    write_worker_meta(worker_dir, meta)
    return worker_dir


def prepare_generation(context: dict, generation_index: int) -> None:
    run_dir = Path(context["run_dir"])
    generation_dir = run_dir / f"generation_{generation_index:03d}"
    generation_dir.mkdir(parents=True, exist_ok=True)

    generation_json = generation_dir / "generation.json"
    meta = read_json(generation_json)
    if not meta:
        meta = {
            "generation": generation_index,
            "created": now_iso(),
            "status": "running",
            "workers": [],
        }
    else:
        meta["status"] = "running"
        meta["resumed_or_updated"] = now_iso()

    write_json(generation_json, meta)

    config = context["config"]
    target = int(config["workers_per_generation"])
    existing = existing_generation_worker_dirs(run_dir, generation_index)

    if existing:
        return

    for worker_index in range(target):
        worker_dir = get_worker_dir(run_dir, generation_index, worker_index)
        worker_dir.mkdir(parents=True, exist_ok=True)
        start_state = choose_start_state_for_new_worker(context, generation_index, worker_index)

        meta = {
            "worker_index": worker_index,
            "generation": generation_index,
            "status": "queued",
            "created": now_iso(),
            "attempt": 0,
            "max_retries": int(config.get("native_retry_limit", 1)),
            "start_state": start_state,
            "is_replacement": False,
            "replacement_for_worker": None,
            "lineage": start_state,
        }

        write_worker_meta(worker_dir, meta)


def parse_slurm_time_to_seconds(t: str) -> Optional[int]:
    if not t or t in {"N/A", "UNLIMITED", "NOT_SET"}:
        return None
    try:
        if "-" in t:
            days, rest = t.split("-", 1)
            days = int(days)
        else:
            days, rest = 0, t

        parts = rest.split(":")
        if len(parts) == 3:
            h, m, s = map(int, parts)
        elif len(parts) == 2:
            h = 0
            m, s = map(int, parts)
        else:
            return None

        return days * 86400 + h * 3600 + m * 60 + s
    except Exception:
        return None


def slurm_remaining_seconds() -> Optional[int]:
    job_id = os.environ.get("SLURM_JOB_ID")
    if not job_id:
        return None

    try:
        out = subprocess.check_output(["squeue", "-h", "-j", job_id, "-o", "%L"], text=True, stderr=subprocess.DEVNULL).strip()
        return parse_slurm_time_to_seconds(out)
    except Exception:
        return None


def should_stop_launching_for_walltime(config: dict) -> bool:
    remaining = slurm_remaining_seconds()
    if remaining is None:
        return False
    buffer_s = float(config.get("slurm_walltime_buffer_min", 20.0)) * 60.0
    return remaining <= buffer_s


def detect_gpu_ids(explicit_gpu_ids: Optional[str]) -> List[str]:
    env_visible = os.environ.get("CUDA_VISIBLE_DEVICES")

    if env_visible and env_visible.strip():
        ids = [x.strip() for x in env_visible.split(",") if x.strip()]
        if ids:
            return ids

    try:
        out = subprocess.check_output(["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"], text=True, stderr=subprocess.DEVNULL)
        ids = [line.strip() for line in out.splitlines() if line.strip()]
        if ids:
            return ids
    except Exception:
        pass

    if explicit_gpu_ids:
        ids = [x.strip() for x in explicit_gpu_ids.split(",") if x.strip()]
        if ids:
            return ids

    return ["0"]


def launch_worker_subprocess(script_path: Path, run_dir: Path, generation_index: int, worker_dir: Path, gpu_id: Optional[str], platform: str) -> subprocess.Popen:
    env = os.environ.copy()
    if platform == "CUDA" and gpu_id is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    stdout_path = worker_dir / "worker_stdout.log"
    stderr_path = worker_dir / "worker_stderr.log"

    cmd = [
        sys.executable,
        str(script_path),
        "--run_single_worker",
        "--resume",
        str(run_dir),
        "--generation_index",
        str(generation_index),
        "--worker_dir",
        str(worker_dir),
    ]

    with open(stdout_path, "w") as stdout_handle, open(stderr_path, "w") as stderr_handle:
        proc = subprocess.Popen(cmd, stdout=stdout_handle, stderr=stderr_handle, env=env)

    meta = read_worker_meta(worker_dir)
    meta["status"] = "running"
    meta["launched_at"] = now_iso()
    meta["pid"] = proc.pid
    meta["assigned_gpu_id"] = gpu_id
    meta["platform"] = platform
    write_worker_meta(worker_dir, meta)

    return proc


def run_generation_controller(context: dict, args) -> dict:
    run_dir = Path(context["run_dir"])
    config = context["config"]
    generation_index = int(context["next_generation"])
    platform = config.get("platform", "CPU")
    script_path = Path(__file__).resolve()

    prepare_generation(context, generation_index)

    if platform == "CUDA":
        gpu_ids = detect_gpu_ids(args.gpu_ids)
    else:
        gpu_ids = [None]

    concurrent = int(args.concurrent_workers or 1)
    if platform == "CUDA":
        concurrent = min(concurrent, len(gpu_ids))
    else:
        concurrent = max(1, concurrent)

    active: Dict[Path, subprocess.Popen] = {}
    launch_cursor = 0

    while True:
        successful = get_successful_worker_dirs(run_dir, generation_index, config)

        if len(successful) >= int(config["workers_per_generation"]):
            break

        if should_stop_launching_for_walltime(config):
            if not active:
                print("Remaining SLURM walltime is low. No active workers remain; exiting cleanly.", flush=True)
                raise SystemExit(0)
        else:
            while len(active) < concurrent:
                worker_dir = prepare_worker_if_needed(context, generation_index, active_worker_dirs=set(active.keys()))
                if worker_dir is None:
                    break

                worker_dir = worker_dir.resolve()

                if worker_dir in active:
                    break

                meta = read_worker_meta(worker_dir)
                if meta.get("status") in SKIP_STATUSES:
                    break

                gpu_id = gpu_ids[launch_cursor % len(gpu_ids)] if platform == "CUDA" else None
                launch_cursor += 1
                proc = launch_worker_subprocess(script_path, run_dir, generation_index, worker_dir, gpu_id, platform)
                active[worker_dir] = proc
                print(f"Launched {worker_dir.name} on {'GPU ' + str(gpu_id) if gpu_id is not None else platform}", flush=True)

                successful = get_successful_worker_dirs(run_dir, generation_index, config)
                if len(successful) >= int(config["workers_per_generation"]):
                    break

        if should_stop_launching_for_walltime(config) and active:
            print("Remaining SLURM walltime is low. Requesting active workers to checkpoint and terminate.", flush=True)
            for worker_dir, proc in list(active.items()):
                try:
                    proc.terminate()
                except Exception:
                    pass

        time.sleep(10)

        for worker_dir, proc in list(active.items()):
            ret = proc.poll()
            if ret is None:
                continue

            active.pop(worker_dir)

            meta = read_worker_meta(worker_dir)
            status = meta.get("status", "unknown")

            if ret == 0 and is_worker_valid_complete(worker_dir, config):
                print(f"{worker_dir.name} completed.", flush=True)
                continue

            if status in TERMINAL_SUCCESS_STATUSES:
                print(f"{worker_dir.name} completed.", flush=True)
                continue

            attempt = int(meta.get("attempt", 0))
            max_retries = int(meta.get("max_retries", config.get("native_retry_limit", 1)))

            if status == "terminated_walltime":
                print(f"{worker_dir.name} stopped cleanly for walltime; will resume later.", flush=True)
                continue

            if attempt < max_retries:
                meta["attempt"] = attempt + 1
                meta["status"] = "queued"
                meta["retry_reason"] = f"Subprocess returned {ret}; retrying with new velocities/seed where applicable."
                meta["queued_retry_at"] = now_iso()
                write_worker_meta(worker_dir, meta)
                print(f"{worker_dir.name} failed with return code {ret}; queued retry {attempt + 1}/{max_retries}.", flush=True)
            else:
                archive_failed_worker(worker_dir, f"Subprocess returned {ret} after {attempt} retries.")
                print(f"{worker_dir.name} archived as failed. A replacement worker will be created.", flush=True)

    return finalize_generation(context, generation_index)


def finalize_generation(context: dict, generation_index: int) -> dict:
    run_dir = Path(context["run_dir"])
    config = context["config"]
    generation_dir = run_dir / f"generation_{generation_index:03d}"

    _, top, _ = build_topology_and_system(Path(context["input_dir"]), config)
    ca_indices = get_ca_indices(top.topology, config)

    successful = get_successful_worker_dirs(run_dir, generation_index, config)
    successful = successful[: int(config["workers_per_generation"])]

    rows = []

    for worker_dir in successful:
        meta = read_worker_meta(worker_dir)
        best_score = meta.get("best_score_nm", None)
        if best_score is None:
            best_score = read_best_score_from_csv(worker_dir)

        rows.append(
            {
                "worker_dir": str(worker_dir),
                "worker": worker_dir.name,
                "worker_index": meta.get("worker_index"),
                "status": meta.get("status"),
                "best_score_nm": best_score,
                "best_score_angstrom": None if best_score is None else best_score * 10.0,
                "best_step": meta.get("best_step"),
                "best_state_npz": str(worker_dir / "best_opening_frame.npz"),
                "best_gro": str(worker_dir / "best_opening_frame.gro"),
                "trajectory_xtc": str(worker_dir / "production.xtc"),
                "score_csv": str(worker_dir / "opening_scores.csv"),
                "lineage": meta.get("lineage", {}),
                "is_replacement": meta.get("is_replacement", False),
            }
        )

    ranked = sorted([r for r in rows if r["best_score_nm"] is not None], key=lambda x: x["best_score_nm"], reverse=True)
    selected = ranked[: int(config["top_frames_to_keep"])]

    gen_meta = read_json(generation_dir / "generation.json")
    gen_meta.update(
        {
            "generation": generation_index,
            "status": "completed",
            "completed": now_iso(),
            "ca_indices": ca_indices,
            "workers": rows,
            "selected_top_frames": selected,
            "best_score_nm": selected[0]["best_score_nm"] if selected else None,
            "best_score_angstrom": selected[0]["best_score_angstrom"] if selected else None,
            "opening_metric_uses_pbc": True,
        }
    )
    write_json(generation_dir / "generation.json", gen_meta)

    context["completed_generations"].append(
        {
            "generation": generation_index,
            "generation_dir": str(generation_dir),
            "best_score_nm": gen_meta["best_score_nm"],
            "best_score_angstrom": gen_meta["best_score_angstrom"],
            "platform": config["platform"],
            "production_ns": config["production_ns"],
            "workers_per_generation": config["workers_per_generation"],
            "opening_metric_uses_pbc": True,
        }
    )

    context["selected_frames"].append(
        {
            "generation": generation_index,
            "frames": [
                {
                    "rank": i + 1,
                    "worker_dir": frame["worker_dir"],
                    "best_state_npz": frame["best_state_npz"],
                    "best_gro": frame["best_gro"],
                    "best_score_nm": frame["best_score_nm"],
                    "best_score_angstrom": frame["best_score_angstrom"],
                    "best_step": frame["best_step"],
                    "lineage": frame.get("lineage", {}),
                    "is_replacement": frame.get("is_replacement", False),
                    "opening_metric_uses_pbc": True,
                }
                for i, frame in enumerate(selected)
            ],
        }
    )

    for frame in selected:
        context.setdefault("lineage", []).append(
            {
                "generation": generation_index,
                "selected_worker_dir": frame["worker_dir"],
                "selected_best_gro": frame["best_gro"],
                "selected_best_state_npz": frame["best_state_npz"],
                "selected_best_step": frame["best_step"],
                "selected_best_score_angstrom": frame["best_score_angstrom"],
                "parent": frame.get("lineage", {}),
                "opening_metric_uses_pbc": True,
            }
        )

    context["next_generation"] = generation_index + 1
    context["last_updated"] = now_iso()
    context["opening_metric_uses_pbc"] = True
    write_json(run_dir / "context.json", context)
    return gen_meta


def read_best_score_from_csv(worker_dir: Path) -> Optional[float]:
    csv_path = worker_dir / "opening_scores.csv"
    if not csv_path.exists():
        return None
    try:
        arr = np.genfromtxt(csv_path, delimiter=",", names=True)
        if arr.size == 0:
            return None
        arr = np.atleast_1d(arr)
        best_a = float(np.nanmax(arr["opening_score_angstrom"]))
        return best_a / 10.0
    except Exception:
        return None


def run_single_worker_from_args(args) -> None:
    context = load_context(Path(args.resume))
    generation_index = int(args.generation_index)
    worker_dir = Path(args.worker_dir).resolve()
    run_single_worker(context, generation_index, worker_dir)


def run_single_worker(context: dict, generation_index: int, worker_dir: Path) -> None:
    config = context["config"]
    input_dir = Path(context["input_dir"])

    worker_dir.mkdir(parents=True, exist_ok=True)
    meta = read_worker_meta(worker_dir)
    worker_index = int(meta.get("worker_index", worker_dir.name.split("_")[-1]))
    attempt = int(meta.get("attempt", 0))

    seed = deterministic_seed(config, generation_index, worker_index, attempt, salt=0)

    meta.update(
        {
            "status": "running",
            "started_or_resumed": now_iso(),
            "generation": generation_index,
            "worker_index": worker_index,
            "attempt": attempt,
            "seed": seed,
            "platform": config["platform"],
            "production_ns": config["production_ns"],
            "opening_metric_uses_pbc": True,
        }
    )
    write_worker_meta(worker_dir, meta)

    try:
        gro, top, system_nvt = build_topology_and_system(input_dir, config)
        _, _, system_npt = build_topology_and_system(input_dir, config)

        ca_indices = get_ca_indices(top.topology, config)

        nvt_steps = steps_from_ps(config["nvt_ps"], config["timestep_ps"])
        npt_equil_steps = steps_from_ps(config["npt_equil_ps"], config["timestep_ps"])
        production_steps = steps_from_ns(config["production_ns"], config["timestep_ps"])
        trajectory_interval_steps = steps_from_ps(config["trajectory_interval_ps"], config["timestep_ps"])
        score_interval_steps = steps_from_ps(config["score_interval_ps"], config["timestep_ps"])
        checkpoint_interval_steps = steps_from_ps(config["checkpoint_interval_ps"], config["timestep_ps"])
        completion_steps = int(math.ceil(float(config.get("completion_fraction", 0.95)) * production_steps))

        if is_worker_valid_complete(worker_dir, config):
            return

        if not (worker_dir / "after_npt_equil.npz").exists() and not (worker_dir / "production_checkpoint.chk").exists():
            run_min_nvt_npt(
                worker_dir=worker_dir,
                gro=gro,
                top=top,
                system_nvt=system_nvt,
                system_npt=system_npt,
                config=config,
                seed=seed,
                nvt_steps=nvt_steps,
                npt_equil_steps=npt_equil_steps,
                initial_start_state=meta.get("start_state", {}),
            )

        run_or_resume_production(
            worker_dir=worker_dir,
            top=top,
            system_npt=system_npt,
            config=config,
            seed=seed,
            ca_indices=ca_indices,
            production_steps=production_steps,
            completion_steps=completion_steps,
            trajectory_interval_steps=trajectory_interval_steps,
            score_interval_steps=score_interval_steps,
            checkpoint_interval_steps=checkpoint_interval_steps,
        )

    except Exception as exc:
        meta = read_worker_meta(worker_dir)
        meta["status"] = "failed_retryable"
        meta["failed_at"] = now_iso()
        meta["exception"] = str(exc)
        meta["traceback"] = traceback.format_exc()
        write_worker_meta(worker_dir, meta)
        raise


def run_min_nvt_npt(worker_dir: Path, gro, top, system_nvt, system_npt, config: dict, seed: int, nvt_steps: int, npt_equil_steps: int, initial_start_state: dict) -> None:
    if initial_start_state and initial_start_state.get("state_npz"):
        initial_positions, initial_box_vectors = load_state_npz(Path(initial_start_state["state_npz"]))
    else:
        initial_positions = gro.positions
        initial_box_vectors = gro.getPeriodicBoxVectors()

    sim_nvt = create_simulation(top.topology, system_nvt, config, seed, add_barostat=False)
    sim_nvt.context.setPositions(initial_positions)
    sim_nvt.context.setPeriodicBoxVectors(*initial_box_vectors)
    sim_nvt.context.setVelocitiesToTemperature(float(config["temperature_k"]) * unit.kelvin, seed)

    sim_nvt.reporters.append(
        StateDataReporter(
            str(worker_dir / "nvt.log"),
            max(1, nvt_steps // 10),
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
        )
    )

    sim_nvt.minimizeEnergy(maxIterations=int(config["minimize_max_iterations"]))
    state_after_min = sim_nvt.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    save_state_npz(worker_dir / "after_minimization.npz", state_after_min.getPositions(), state_after_min.getPeriodicBoxVectors())

    sim_nvt.step(nvt_steps)
    state_after_nvt = sim_nvt.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    save_state_npz(worker_dir / "after_nvt.npz", state_after_nvt.getPositions(), state_after_nvt.getPeriodicBoxVectors())

    sim_npt = create_simulation(top.topology, system_npt, config, seed + 17, add_barostat=True)
    sim_npt.context.setPositions(state_after_nvt.getPositions())
    sim_npt.context.setPeriodicBoxVectors(*state_after_nvt.getPeriodicBoxVectors())
    sim_npt.context.setVelocitiesToTemperature(float(config["temperature_k"]) * unit.kelvin, seed + 17)

    sim_npt.reporters.append(
        StateDataReporter(
            str(worker_dir / "npt_equil.log"),
            max(1, npt_equil_steps // 10),
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
        )
    )

    sim_npt.step(npt_equil_steps)
    state_after_npt = sim_npt.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    save_state_npz(worker_dir / "after_npt_equil.npz", state_after_npt.getPositions(), state_after_npt.getPeriodicBoxVectors())

    meta = read_worker_meta(worker_dir)
    meta["equilibration_completed"] = True
    meta["after_npt_equil_npz"] = str(worker_dir / "after_npt_equil.npz")
    meta["updated"] = now_iso()
    write_worker_meta(worker_dir, meta)


def archive_existing_xtc_if_resuming(worker_dir: Path, trajectory_path: Path, resumed_from_checkpoint: bool) -> Path:
    if trajectory_path.exists() and resumed_from_checkpoint:
        archive_dir = worker_dir / "previous_xtc_segments"
        archive_dir.mkdir(parents=True, exist_ok=True)
        archived = archive_dir / f"production_before_resume_{make_timestamp()}_{os.getpid()}.xtc"
        trajectory_path.rename(archived)

        meta = read_worker_meta(worker_dir)
        meta.setdefault("archived_xtc_segments", []).append(str(archived))
        meta["xtc_resume_policy"] = "production.xtc archived on checkpoint resume; fresh XTC segment started"
        meta["xtc_archived_at"] = now_iso()
        write_worker_meta(worker_dir, meta)

    return trajectory_path


def run_or_resume_production(
    worker_dir: Path,
    top,
    system_npt,
    config: dict,
    seed: int,
    ca_indices: Dict[str, List[int]],
    production_steps: int,
    completion_steps: int,
    trajectory_interval_steps: int,
    score_interval_steps: int,
    checkpoint_interval_steps: int,
) -> None:
    checkpoint_path = worker_dir / "production_checkpoint.chk"
    trajectory_path = worker_dir / "production.xtc"
    score_csv_path = worker_dir / "opening_scores.csv"
    best_state_path = worker_dir / "best_opening_frame.npz"
    best_gro_path = worker_dir / "best_opening_frame.gro"

    sim_prod = create_simulation(top.topology, system_npt, config, seed + 31, add_barostat=True)
    resumed_from_checkpoint = False

    if checkpoint_path.exists():
        try:
            sim_prod.loadCheckpoint(str(checkpoint_path))
            resumed_from_checkpoint = True
        except Exception as exc:
            meta = read_worker_meta(worker_dir)
            meta["checkpoint_load_failed"] = str(exc)
            meta["checkpoint_load_failed_at"] = now_iso()
            write_worker_meta(worker_dir, meta)

            try:
                checkpoint_path.rename(worker_dir / f"corrupt_production_checkpoint_{make_timestamp()}.chk")
            except Exception:
                pass

            if not (worker_dir / "after_npt_equil.npz").exists():
                raise RuntimeError("Checkpoint was corrupt and after_npt_equil.npz does not exist.")

    if not resumed_from_checkpoint:
        positions, box_vectors = load_state_npz(worker_dir / "after_npt_equil.npz")
        sim_prod.context.setPositions(positions)
        sim_prod.context.setPeriodicBoxVectors(*box_vectors)
        sim_prod.context.setVelocitiesToTemperature(float(config["temperature_k"]) * unit.kelvin, seed + 31)

    archive_existing_xtc_if_resuming(worker_dir, trajectory_path, resumed_from_checkpoint)
    append_xtc = trajectory_path.exists() and not resumed_from_checkpoint

    sim_prod.reporters.append(XTCReporter(str(trajectory_path), trajectory_interval_steps, append=append_xtc))

    sim_prod.reporters.append(
        StateDataReporter(
            str(worker_dir / "production.log"),
            max(1, production_steps // 100),
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            append=(worker_dir / "production.log").exists(),
        )
    )

    score_reporter = OpeningScoreReporter(
        file_path=score_csv_path,
        report_interval=score_interval_steps,
        ca_indices=ca_indices,
        best_npz_path=best_state_path,
        best_gro_path=best_gro_path,
        topology=top.topology,
        append=True,
    )

    sim_prod.reporters.append(score_reporter)

    meta = read_worker_meta(worker_dir)
    meta["status"] = "running_production"
    meta["resumed_from_checkpoint"] = resumed_from_checkpoint
    meta["production_target_steps"] = production_steps
    meta["production_completion_steps"] = completion_steps
    meta["checkpoint_interval_steps"] = checkpoint_interval_steps
    meta["checkpoint_interval_ps"] = config["checkpoint_interval_ps"]
    meta["updated"] = now_iso()
    meta["opening_metric_uses_pbc"] = True
    write_worker_meta(worker_dir, meta)

    while sim_prod.currentStep < production_steps:
        if STOP_REQUESTED:
            save_openmm_checkpoint_atomic(sim_prod, checkpoint_path)
            meta = read_worker_meta(worker_dir)
            meta["status"] = "terminated_walltime"
            meta["terminated_walltime_at"] = now_iso()
            meta["current_step"] = sim_prod.currentStep
            meta["completion_fraction"] = worker_completion_fraction(worker_dir, config)
            write_worker_meta(worker_dir, meta)
            score_reporter.close()
            return

        remaining = production_steps - sim_prod.currentStep
        chunk = min(checkpoint_interval_steps, remaining)

        sim_prod.step(chunk)
        save_openmm_checkpoint_atomic(sim_prod, checkpoint_path)

        state = sim_prod.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
        save_state_npz(worker_dir / "latest_production_state.npz", state.getPositions(), state.getPeriodicBoxVectors())

        meta = read_worker_meta(worker_dir)
        meta["status"] = "checkpointed"
        meta["current_step"] = sim_prod.currentStep
        meta["last_checkpoint_at"] = now_iso()
        meta["completion_fraction"] = worker_completion_fraction(worker_dir, config)
        meta["best_score_nm"] = score_reporter.best_score
        meta["best_score_angstrom"] = score_reporter.best_score * 10.0
        meta["best_step"] = score_reporter.best_step
        meta["opening_metric_uses_pbc"] = True
        write_worker_meta(worker_dir, meta)

    final_state = sim_prod.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    save_state_npz(worker_dir / "final_frame.npz", final_state.getPositions(), final_state.getPeriodicBoxVectors())

    meta = read_worker_meta(worker_dir)
    meta.update(
        {
            "status": "completed",
            "completed": now_iso(),
            "current_step": sim_prod.currentStep,
            "completion_fraction": 1.0,
            "best_score_nm": score_reporter.best_score,
            "best_score_angstrom": score_reporter.best_score * 10.0,
            "best_step": score_reporter.best_step,
            "best_state_npz": str(best_state_path),
            "best_gro": str(best_gro_path),
            "trajectory_xtc": str(trajectory_path),
            "score_csv": str(score_csv_path),
            "checkpoint": str(checkpoint_path),
            "opening_metric_uses_pbc": True,
        }
    )

    write_worker_meta(worker_dir, meta)
    score_reporter.close()


def validate_args(args):
    if args.run_single_worker:
        if args.resume is None or args.generation_index is None or args.worker_dir is None:
            raise SystemExit("--run_single_worker requires --resume, --generation_index, and --worker_dir.")
        return

    if args.resume is None and args.input_dir is None:
        raise SystemExit("Use either --input_dir for a fresh run or --resume for an existing adaptive run directory.")

    if args.resume is not None and args.input_dir is not None:
        raise SystemExit("Use only one of --input_dir or --resume, not both.")

    if args.mode == "auto" and args.target_generation is None:
        raise SystemExit("Auto mode requires --target_generation N.")

    if args.target_generation is not None and args.target_generation < 0:
        raise SystemExit("--target_generation must be 0 or greater.")

    if args.production_ns is not None and args.production_ns <= 0:
        raise SystemExit("--production_ns must be greater than 0.")

    if args.workers_per_generation is not None and args.workers_per_generation < 1:
        raise SystemExit("--workers_per_generation must be at least 1.")

    if args.concurrent_workers < 1:
        raise SystemExit("--concurrent_workers must be at least 1.")


def build_overrides_from_args(args) -> dict:
    overrides = {}
    for name in [
        "input_gro",
        "input_top",
        "system_variant",
        "run_label",
        "collagen_residue_start",
        "collagen_residue_end",
        "mutated_chain",
        "platform",
        "production_ns",
        "workers_per_generation",
        "checkpoint_interval_ps",
        "master_seed",
        "slurm_walltime_buffer_min",
    ]:
        value = getattr(args, name, None)
        if value is not None:
            overrides[name] = value
    return overrides


def main():
    parser = argparse.ArgumentParser(
        description="Adaptive OpenMM MMP1-collagen unwinding simulation runner with PBC-corrected opening-score selection.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--input_dir", type=str, default=None, help="Directory containing system.top, a .gro coordinate file, topology includes, and force-field include files.")
    parser.add_argument("--resume", type=str, default=None, help="Existing adaptive_run_* directory containing context.json.")
    parser.add_argument("--input_gro", type=str, default=None, help="Input GRO filename inside --input_dir. Use 'auto' or omit to auto-detect.")
    parser.add_argument("--input_top", type=str, default=None, help="Input topology filename inside --input_dir.")
    parser.add_argument("--system_variant", type=str, choices=["auto"] + sorted(SYSTEM_VARIANTS.keys()), default=None, help="System variant for mutation-aware residue validation and run-folder naming.")
    parser.add_argument("--run_label", type=str, default=None, help="Optional label inserted into fresh adaptive run directory names.")
    parser.add_argument("--collagen_residue_start", type=int, default=None, help="First collagen residue included in the CA opening-score window.")
    parser.add_argument("--collagen_residue_end", type=int, default=None, help="Last collagen residue included in the CA opening-score window.")
    parser.add_argument("--mutated_chain", type=str, choices=["chain_1", "chain_2", "chain_3"], default=None, help="Mutated collagen chain for mutation-aware validation. Defaults from --system_variant.")
    parser.add_argument("--platform", type=str, choices=["CPU", "CUDA"], default=None, help="OpenMM platform. Use CUDA for NVIDIA GPU execution.")
    parser.add_argument("--mode", type=str, choices=["manual", "auto"], default="manual", help="manual runs one generation and stops. auto runs until --target_generation is completed.")
    parser.add_argument("--target_generation", type=int, default=None, help="In auto mode, run until this generation index has completed.")
    parser.add_argument("--production_ns", type=float, default=None, help="Length of each worker production run in ns.")
    parser.add_argument("--workers_per_generation", type=int, default=None, help="Number of valid completed workers required per generation.")
    parser.add_argument("--concurrent_workers", type=int, default=1, help="Number of workers to run simultaneously. Use 1 for backwards-compatible serial behaviour.")
    parser.add_argument("--gpu_ids", type=str, default=None, help="Comma-separated GPU IDs to use if CUDA auto-detection fails or if you want explicit control, e.g. 0,1.")
    parser.add_argument("--checkpoint_interval_ps", type=float, default=None, help="Rolling production checkpoint interval in ps.")
    parser.add_argument("--master_seed", type=int, default=None, help="Master seed for deterministic worker seed generation.")
    parser.add_argument("--slurm_walltime_buffer_min", type=float, default=None, help="Stop launching new workers when SLURM remaining walltime is below this many minutes.")

    parser.add_argument("--run_single_worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--generation_index", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--worker_dir", type=str, default=None, help=argparse.SUPPRESS)

    args = parser.parse_args()
    validate_args(args)

    if args.run_single_worker:
        run_single_worker_from_args(args)
        return

    overrides = build_overrides_from_args(args)

    if args.resume is not None:
        run_dir = Path(args.resume).resolve()
        context = load_context(run_dir)
        context = apply_runtime_overrides(context, overrides)
    else:
        context = initialise_context(args.input_dir, overrides)

    print(f"Run directory: {context['run_dir']}")
    print(f"Current next generation: {context['next_generation']}")
    print(f"System variant: {context['config']['system_variant']}")
    print(f"Input GRO: {context['config']['input_gro']}")
    print(f"Mode: {args.mode}")
    print(f"Platform: {context['config']['platform']}")
    print(f"Production length per worker: {context['config']['production_ns']} ns")
    print(f"Workers per generation: {context['config']['workers_per_generation']}")
    print(f"Concurrent workers: {args.concurrent_workers}")
    print(f"Checkpoint interval: {context['config']['checkpoint_interval_ps']} ps")
    print("Opening score: PBC-corrected minimum-image CA distance metric")

    if args.mode == "manual":
        print("This invocation will run one generation and then stop.")
        result = run_generation_controller(context, args)
        print("")
        print("Generation complete.")
        print(f"Generation: {result['generation']}")
        print(f"Best PBC-corrected opening score: {result['best_score_angstrom']:.3f} Å")
        print(f"Run directory: {context['run_dir']}")
        print("")
        print("Resume next generation with:")
        print(f"python {Path(__file__).name} --resume {context['run_dir']} --platform {context['config']['platform']} --concurrent_workers {args.concurrent_workers}")

    else:
        print(f"Auto mode target generation: {args.target_generation}")
        while int(context["next_generation"]) <= int(args.target_generation):
            generation_index = int(context["next_generation"])
            print("")
            print(f"Starting generation {generation_index}")
            result = run_generation_controller(context, args)
            context = load_context(Path(context["run_dir"]))
            print(f"Completed generation {result['generation']}")
            print(f"Best PBC-corrected opening score: {result['best_score_angstrom']:.3f} Å")

        print("")
        print("Auto mode complete.")
        print(f"Completed up to generation {int(context['next_generation']) - 1}.")
        print(f"Run directory: {context['run_dir']}")


if __name__ == "__main__":
    main()
