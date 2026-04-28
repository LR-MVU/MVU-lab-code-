"""Reusable mRNA-related P-body metrics for the reorganized pipeline.

This module refactors the legacy PB index workflow into callable functions so a
separate wrapper can provide user input later.
"""

from __future__ import annotations

import os
import re
from math import pi
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError
from scipy import ndimage
from tifffile import imread
from skimage.measure import label, regionprops


def get_filenames(directory: os.PathLike[str] | str, chan: str, filetype: str = ".tif") -> list[str]:
    """Return sorted filenames in ``directory`` matching a channel token."""

    pattern = re.compile(r".*" + re.escape(chan) + r".*" + re.escape(filetype) + r"$")
    return sorted(
        filename
        for filename in os.listdir(directory)
        if re.match(pattern, filename)
    )


def format_results(rows: list[list[float]]) -> np.ndarray:
    """Pad ragged condition results so they can be written as a table."""

    if not any(len(row) for row in rows):
        return np.empty((0, len(rows)), dtype=object)

    max_len = max(len(row) for row in rows)
    padded = [
        list(row) + [None] * (max_len - len(row))
        for row in rows
    ]
    return np.asarray(padded, dtype=object).T


def _as_bool_mask(mask: np.ndarray) -> np.ndarray:
    return mask.astype(bool)


def _project_to_2d(mask: np.ndarray) -> np.ndarray:
    if mask.ndim == 3:
        return np.max(mask, axis=0)
    return mask


def _read_spot_csv(csv_path: Path) -> np.ndarray:
    """Read a spot CSV into an ``(n, 3)`` integer array."""

    try:
        spots = pd.read_csv(csv_path, header=None).to_numpy()
    except EmptyDataError:
        return np.empty((0, 3), dtype=int)

    if spots.size == 0:
        return np.empty((0, 3), dtype=int)

    if spots.shape[1] < 3:
        raise ValueError(f"Spot table has fewer than 3 columns: {csv_path}")

    return spots[:, :3].astype(int)


def load_rna_spots(
    output_dir: os.PathLike[str] | str,
    mask_basename: str,
    rna_chan: str,
    spot_source: str
) -> tuple[np.ndarray | None, Path]:
    """Load RNA spot coordinates from the explicitly selected output file.

    Parameters
    ----------
    output_dir
        Folder containing RNA spot CSV outputs.
    mask_basename
        Base filename derived from the mask/image set.
    rna_chan
        RNA channel token used in the reorganized output names.
    spot_source
        Required RNA spot source. Supported values are ``"declustered"`` and
        ``"raw"``.

    Returns
    -------
    tuple
        ``(spots, expected_path)`` where ``spots`` is an ``(n, 3)`` array when
        the selected file exists, or ``None`` when the selected file is missing.

    Raises
    ------
    ValueError
        If ``spot_source`` is not one of the supported values.
    """

    output_dir = Path(output_dir)

    if spot_source == "declustered":
        spot_path = output_dir / f"{mask_basename}_mRNA{rna_chan}_spots.csv"
    elif spot_source == "raw":
        spot_path = output_dir / f"{mask_basename}_mRNA{rna_chan}_spots_ind.csv"
    else:
        raise ValueError(
            "spot_source must be either 'declustered' or 'raw'. "
            f"Received: {spot_source!r}"
        )

    if not spot_path.exists():
        print(f"Missing RNA spot file, skipping set: {spot_path}")
        return None, spot_path

    return _read_spot_csv(spot_path), spot_path


def build_rna_mask(image_shape: tuple[int, ...], rna_coords: np.ndarray) -> np.ndarray:
    """Convert RNA spot coordinates into a count mask."""

    rna_mask = np.zeros(image_shape, dtype=np.int32)
    for coord in rna_coords:
        y, x = (int(coord[1]), int(coord[2]))
        if (
            0 <= y < image_shape[0]
            and 0 <= x < image_shape[1]
        ):
            rna_mask[y, x] += 1
    return rna_mask


def localized_pbody_index(
    pb_mask: np.ndarray,
    rna_mask: np.ndarray,
    cyt_mask: np.ndarray,
    cell_area: float,
    pbody_dilation_iterations: int = 2,
    local_shell_iterations: int = 4,
) -> dict[str, list[float] | float]:
    """Compute local and whole-cell mRNA-related P-body metrics.

    ``pbody_dilation_iterations`` controls the first XY dilation step applied
    to PB masks for all calculations. ``local_shell_iterations`` controls how
    many additional 1-pixel XY dilation steps are applied locally to define the
    shell around each PB. The local crop window follows the legacy max-axis
    rule.
    """

    if pbody_dilation_iterations < 1:
        raise ValueError("pbody_dilation_iterations must be at least 1.")
    if local_shell_iterations < 1:
        raise ValueError("local_shell_iterations must be at least 1.")

    cyt_mask = cyt_mask.astype(bool)
    pb_binary = pb_mask > 0
    pb_dilation_kernel = np.ones((3, 3), dtype=np.uint8)

    labeled_pb_mask = label(pb_binary, connectivity=1)
    pb_ids = np.unique(labeled_pb_mask)
    pb_ids = pb_ids[pb_ids != 0]

    local_indices: list[float] = []
    partition_coefficients: list[float] = []

    for pb_id in pb_ids:
        pbody = ndimage.binary_dilation(
            labeled_pb_mask == pb_id,
            structure=pb_dilation_kernel,
            iterations=pbody_dilation_iterations,
        )

        local_shell_object = ndimage.binary_dilation(
            pbody,
            structure=pb_dilation_kernel,
            iterations=local_shell_iterations,
        )

        touches_border = (
            np.any(local_shell_object[0, :]) or
            np.any(local_shell_object[-1, :]) or
            np.any(local_shell_object[:, 0]) or
            np.any(local_shell_object[:, -1])
        )

        if touches_border:
            local_indices.append(np.nan)
            partition_coefficients.append(np.nan)
            continue
        
        shell_outside_cyt = np.any(local_shell_object & (~cyt_mask))
        if shell_outside_cyt:
            local_indices.append(np.nan)
            partition_coefficients.append(np.nan)
            continue

        out_mask = ((local_shell_object) & (~pbody))

        n_in = float(np.sum(rna_mask[pbody]))
        n_out = float(np.sum(rna_mask[out_mask]))

        pb_volume = float(np.count_nonzero(pbody))
        cyt_volume = float(np.count_nonzero(out_mask))
        if pb_volume == 0 or cyt_volume == 0:
            local_indices.append(np.nan)
            partition_coefficients.append(np.nan)
            continue

        total_local_rna = n_in + n_out
        if total_local_rna == 0:
            local_indices.append(np.nan)
            partition_coefficients.append(np.nan)
            continue

        local_index = n_in / (total_local_rna * (pb_volume / (pb_volume + cyt_volume)))
        local_indices.append(local_index)

        if n_out == 0:
            partition_coefficients.append(np.inf if n_in > 0 else np.nan)
        else:
            partition_coefficients.append((n_in / pb_volume) / (n_out / cyt_volume))

    avg_local_index = float(np.nanmean(local_indices)) if local_indices else np.nan
    avg_partition_coefficient = (
        float(np.nanmean(partition_coefficients)) if partition_coefficients else np.nan
    )

    dilated_pb = ndimage.binary_dilation(
        pb_binary,
        structure=pb_dilation_kernel,
        iterations=pbody_dilation_iterations,
    )

    dilated_pb = dilated_pb & cyt_mask
    mrna_in_pb_mask = dilated_pb * rna_mask
    num_rna_in_pb = float(np.sum(mrna_in_pb_mask))
    total_rna = float(np.sum(rna_mask))
    total_pb_area = float(np.count_nonzero(dilated_pb))
    pb_fraction_cell = np.nan if total_rna == 0 else num_rna_in_pb / total_rna

    if total_rna == 0 or total_pb_area == 0 or cell_area == 0:
        total_cell_index = np.nan
    else:
        total_cell_index = num_rna_in_pb / (total_rna * (total_pb_area / cell_area))

    return {
        "pb_local_index": local_indices,
        "pb_local_index_cell": [avg_local_index],
        "pb_total_index_cell": [total_cell_index],
        "partition_coefficient": partition_coefficients,
        "partition_coefficient_cell": [avg_partition_coefficient],
        "pb_fraction": [pb_fraction_cell],
        "rna_in_pbody": [num_rna_in_pb],
    }


def pbody_characteristics(
    pb_mask: np.ndarray,
    pixel_res: float,
    diffraction: float,
) -> dict[str, list[float | int]]:
    """Compute per-PB size/radius values and per-cell PB counts."""

    labeled_pb_mask = label(pb_mask > 0, connectivity=1)
    pb_props = regionprops(labeled_pb_mask)

    pb_sizes: list[float] = []
    pb_radii: list[float] = []
    pb_num = 0

    for prop in pb_props:
        radius = (prop.equivalent_diameter_area * pixel_res / 2) - diffraction
        if radius > 0:
            pb_radii.append(radius)
            pb_sizes.append(radius ** 2 * pi)
            pb_num += 1

    return {
        "pb_size": pb_sizes,
        "pb_radius": pb_radii,
        "pb_number": [pb_num],
    }


def _check_parallel_lists(condition: str, file_lists: dict[str, list[str]]) -> None:
    lengths = {name: len(files) for name, files in file_lists.items()}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"File count mismatch for condition '{condition}': {lengths}")


def collect_condition_metrics(
    path_root: os.PathLike[str] | str,
    condition: str,
    spot_source: str,
    pbody_dilation_iterations: int = 2,
    local_shell_iterations: int = 4,
    pbody_chan: str = "640S",
    rna_chan: str = "555S",
    cell_mask_chan: str = "print",
    nuc_mask_chan: str = "nuc",
    segmentation_folder: str = "Pb_Masks",
    spot_folder: str = "output",
    pixel_res: float = 107.5,
    diffraction: float = 0,
    filter: bool  = False,
    PB_index_threshold: float = None,
) -> dict[str, list[float | int]]:
    """Collect mRNA-related P-body metrics for one condition.

    ``spot_source`` selects exactly one RNA spot CSV format:
    ``"declustered"`` or ``"raw"``. If the selected file is missing for a
    given image/cell set, that set is skipped and the rest of the condition
    continues processing. ``pbody_dilation_iterations`` controls the global PB
    dilation used across calculations, while ``local_shell_iterations``
    controls the additional local shell dilation around each PB.
    """

    path_input = Path(path_root) / condition
    path_pbody = path_input / segmentation_folder
    path_output = path_input / spot_folder

    pbody_mask_files = get_filenames(path_pbody, pbody_chan)
    cell_mask_files = get_filenames(path_input, cell_mask_chan)
    nuc_mask_files = get_filenames(path_input, nuc_mask_chan)

    _check_parallel_lists(
        condition,
        {
            "pbody_mask_files": pbody_mask_files,
            "cell_mask_files": cell_mask_files,
            "nuc_mask_files": nuc_mask_files,
        },
    )

    condition_metrics = {
        "pb_size": [],
        "pb_radius": [],
        "pb_number": [],
        "pb_local_index": [],
        "pb_local_index_cell": [],
        "pb_total_index_cell": [],
        "partition_coefficient": [],
        "partition_coefficient_cell": [],
        "pb_fraction": [],
        "rna_in_pbody": [],
        "rna_count_cell": [],
    }

    for index in range(len(pbody_mask_files)):
        pbody_mask_path = path_pbody / pbody_mask_files[index]
        cell_mask_path = path_input / cell_mask_files[index]
        nuc_mask_path = path_input / nuc_mask_files[index]

        pbody_mask = np.squeeze(imread(pbody_mask_path))
        cell_mask = _as_bool_mask(np.squeeze(imread(cell_mask_path)))
        nuc_mask = _as_bool_mask(np.squeeze(imread(nuc_mask_path)))

        cyt_mask = cell_mask & (~nuc_mask)

        pbody_mask = pbody_mask * cyt_mask

        mask_basename = Path(cell_mask_files[index]).stem
        rna_spots, spot_path = load_rna_spots(
            path_output,
            mask_basename,
            rna_chan,
            spot_source=spot_source,
        )
        if rna_spots is None:
            print(
                "Skipping metric calculation for set because the selected RNA "
                f"spot file was not found: {spot_path}"
            )
            continue

        rna_mask = build_rna_mask(pbody_mask.shape, rna_spots)
        rna_mask = rna_mask * cyt_mask

        cell_area = float(np.count_nonzero(cyt_mask))
        metrics = localized_pbody_index(
            pb_mask=pbody_mask,
            rna_mask=rna_mask,
            cyt_mask=cyt_mask,
            cell_area=cell_area,
            pbody_dilation_iterations=pbody_dilation_iterations,
            local_shell_iterations=local_shell_iterations
        )
        pb_metrics = pbody_characteristics(
            pb_mask=pbody_mask,
            pixel_res=pixel_res,
            diffraction=diffraction,
        )
        if not filter:
            for key in (
                "pb_size",
                "pb_radius",
                "pb_number",
            ):
                condition_metrics[key].extend(pb_metrics[key])
            for key in (
                "pb_local_index",
                "pb_local_index_cell",
                "pb_total_index_cell",
                "partition_coefficient",
                "partition_coefficient_cell",
                "pb_fraction",
                "rna_in_pbody",
            ):
                condition_metrics[key].extend(metrics[key])
            condition_metrics["rna_count_cell"].append(float(np.sum(rna_mask)))
            
        if filter and metrics["pb_total_index_cell"][0] > PB_index_threshold:
            for key in (
                "pb_size",
                "pb_radius",
                "pb_number",
            ):
                condition_metrics[key].extend(pb_metrics[key])
            for key in (
                "pb_local_index",
                "pb_local_index_cell",
                "pb_total_index_cell",
                "partition_coefficient",
                "partition_coefficient_cell",
                "pb_fraction",
                "rna_in_pbody",
            ):
                condition_metrics[key].extend(metrics[key])
            condition_metrics["rna_count_cell"].append(float(np.sum(rna_mask)))

    return condition_metrics


def run_pbody_mrna_metrics(
    path_root: os.PathLike[str] | str,
    conditions: list[str],
    spot_source: str,
    pbody_dilation_iterations: int = 2,
    local_shell_iterations: int = 4,
    pbody_chan: str = "640S",
    rna_chan: str = "555S",
    cell_mask_chan: str = "print",
    nuc_mask_chan: str = "nuc",
    segmentation_folder: str = "Pb_Masks",
    spot_folder: str = "output",
    pixel_res: float = 107.5,
    diffraction: float = 0,
    output_name: str | None = None,
    filter: bool = False,
    PB_index_threshold: float = None,
) -> tuple[dict[str, pd.DataFrame], Path]:
    """Run legacy mRNA-related P-body metrics across multiple conditions.

    ``spot_source`` must be ``"declustered"`` or ``"raw"``. Missing selected
    RNA spot files do not stop the run; the affected image/cell sets are
    skipped, while unsupported ``spot_source`` values raise ``ValueError``.
    ``pbody_dilation_iterations`` controls the global PB dilation used across
    calculations, while ``local_shell_iterations`` controls the additional
    local shell dilation around each PB.

    Returns the dataframes and the workbook path written to disk.
    """

    if not conditions:
        raise ValueError("At least one condition is required.")

    if not filter:
        condition_data = {
            condition: collect_condition_metrics(
                path_root=path_root,
                condition=condition,
                spot_source=spot_source,
                pbody_dilation_iterations=pbody_dilation_iterations,
                local_shell_iterations=local_shell_iterations,
                pbody_chan=pbody_chan,
                rna_chan=rna_chan,
                cell_mask_chan=cell_mask_chan,
                nuc_mask_chan=nuc_mask_chan,
                segmentation_folder=segmentation_folder,
                spot_folder=spot_folder,
                pixel_res=pixel_res,
                diffraction=diffraction
            )
            for condition in conditions
        }
    else:
        condition_data = {
            condition: collect_condition_metrics(
                path_root=path_root,
                condition=condition,
                spot_source=spot_source,
                pbody_dilation_iterations=pbody_dilation_iterations,
                local_shell_iterations=local_shell_iterations,
                pbody_chan=pbody_chan,
                rna_chan=rna_chan,
                cell_mask_chan=cell_mask_chan,
                nuc_mask_chan=nuc_mask_chan,
                segmentation_folder=segmentation_folder,
                spot_folder=spot_folder,
                pixel_res=pixel_res,
                diffraction=diffraction,
                filter=filter,
                PB_index_threshold=PB_index_threshold
            )
            for condition in conditions
        }

    sheet_map = {
        "PB Size": "pb_size",
        "PB Radius": "pb_radius",
        "PB Number": "pb_number",
        "PB Local Index": "pb_local_index",
        "PB Local Index Cell": "pb_local_index_cell",
        "PB Total Index Cell": "pb_total_index_cell",
        "Partition Coefficient": "partition_coefficient",
        "Partition Coefficient Cell": "partition_coefficient_cell",
        "PB Fraction": "pb_fraction",
        "RNA In PB": "rna_in_pbody",
        "RNA Count Cell": "rna_count_cell",
    }

    dataframes: dict[str, pd.DataFrame] = {}
    for sheet_name, metric_key in sheet_map.items():
        rows = [condition_data[condition][metric_key] for condition in conditions]
        dataframes[sheet_name] = pd.DataFrame(format_results(rows), columns=conditions)

    output_path = Path(path_root)
    if output_name is None:
        if not filter:
            output_name = f"PBody_mRNA_metrics_{pbody_chan}_{rna_chan}.xlsx"
        else:
            output_name = f"PBody_mRNA_metrics_FILTERED_{pbody_chan}_{rna_chan}.xlsx"
    workbook_path = output_path / output_name

    with pd.ExcelWriter(workbook_path, engine="openpyxl") as writer:
        for sheet_name, dataframe in dataframes.items():
            dataframe.to_excel(writer, sheet_name=sheet_name, index=False)

    return dataframes, workbook_path


