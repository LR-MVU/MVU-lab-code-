#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Standalone PB mask generation from bigFISH input folders."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile as tf
import tifffile as tiff
from skimage.filters import difference_of_gaussians
from skimage.measure import label, regionprops
#import cv2

import bigfish.segmentation as seg
import bigfish.stack as stack


VALID_MASK_SUFFIXES = {'.tif', '.tiff', '.png'}
VALID_IMAGE_SUFFIXES = {'.tif', '.tiff'}


def get_histogramm_highest_varation_value(array, bins=None):
    """Return the histogram value at the largest positive first derivative."""
    stack.check_parameter(array=np.ndarray)

    if bins is None:
        if np.issubdtype(array.dtype, np.floating):
            bins = 100
        else:
            bins = int(array.max() - array.min())
            bins = max(bins, 1)

    count, values = np.histogram(array, bins=bins)
    gradient = np.gradient(count)
    gradient[gradient < 0] = 0
    max_derivative_index = int(np.argmax(np.abs(gradient)))
    return values[max_derivative_index]


def max_project_to_2d(image):
    """Convert a 3D stack to a 2D max projection, or return a 2D image unchanged."""
    image = np.asarray(image)
    if image.ndim == 2:
        return image.astype(np.float32, copy=False)
    if image.ndim == 3:
        return image.max(axis=0).astype(np.float32, copy=False)
    raise ValueError(f'Expected a 2D or 3D image, got shape {image.shape}.')


def load_mask_2d(mask_path):
    """Load a cell or nucleus mask and convert it to a 2D boolean mask."""
    mask = tf.imread(mask_path)
    mask_2d = max_project_to_2d(mask)
    return mask_2d > 0


def build_match_key(stem, token):
    """Split a filename stem around a token and return a single shared match key."""
    parts = stem.split(token)
    if len(parts) != 2:
        return None
    return ''.join(parts)


def collect_condition_files(
    condition_folder,
    pb_channel,
    nucleus_mask_token,
    cell_mask_token,
):
    """Match PB image stacks to their corresponding cell and nucleus masks."""

    pb_images = {}
    nucleus_masks = {}
    cell_masks = {}

    for path in condition_folder.iterdir():
        if not path.is_file():
            continue

        suffix = path.suffix.lower()
        stem = path.stem

        if suffix in VALID_IMAGE_SUFFIXES and pb_channel in stem:
            match_key = build_match_key(stem, pb_channel)
            if match_key is not None:
                pb_images[match_key] = path
            continue

        if suffix in VALID_MASK_SUFFIXES and nucleus_mask_token in stem:
            match_key = build_match_key(stem, nucleus_mask_token)
            if match_key is not None:
                nucleus_masks[match_key] = path
            continue

        if suffix in VALID_MASK_SUFFIXES and cell_mask_token in stem:
            match_key = build_match_key(stem, cell_mask_token)
            if match_key is not None:
                cell_masks[match_key] = path

    matched_files = []
    missing_pairs = []

    for match_key, pb_image_path in sorted(pb_images.items()):
        cell_mask_path = cell_masks.get(match_key)
        nucleus_mask_path = nucleus_masks.get(match_key)

        if cell_mask_path is None or nucleus_mask_path is None:
            missing_pairs.append((pb_image_path.name, cell_mask_path is not None, nucleus_mask_path is not None))
            continue

        matched_files.append(
            {
                'match_key': match_key,
                'pb_image_path': pb_image_path,
                'cell_mask_path': cell_mask_path,
                'nucleus_mask_path': nucleus_mask_path,
            }
        )

    if missing_pairs:
        print(f'Skipping files with missing masks in {condition_folder}:')
        for name, has_cell, has_nucleus in missing_pairs:
            print(f'  - {name} (cell={has_cell}, nucleus={has_nucleus})')

    return matched_files


def compute_circularity(region):
    """Compute 2D circularity from a labeled connected component region."""
    perimeter = region.perimeter
    if perimeter <= 0:
        return 0.0
    return float(4.0 * np.pi * region.area / (perimeter ** 2))


def filter_pb_components(binary_mask, min_area, max_area, min_circularity, original_image, minimum_intensity):
    """Keep connected components within area and circularity bounds."""
    labeled_mask = label(binary_mask, connectivity=1)
    filtered_mask = np.zeros_like(binary_mask, dtype=bool)

    for region in regionprops(labeled_mask, intensity_image=original_image):
        area = region.area
        circularity = compute_circularity(region)

        if area < min_area:
            continue
        if area > max_area:
            continue
        if circularity < min_circularity:
            continue
        if region.intensity_max < minimum_intensity:
            continue

        filtered_mask[labeled_mask == region.label] = True

    return filtered_mask


def get_component_areas(binary_mask):
    """Measure connected-component areas from a 2D binary mask."""
    labeled_mask = label(binary_mask, connectivity=1)
    return [region.area for region in regionprops(labeled_mask)]


def save_visualizations(visualization_output_folder, visualization_records):
    """Write the main 2x2 quality-control figures into the visualization folder."""
    for record in visualization_records:
        file_name = record['file_name']
        original_image = record['original_image']
        unfiltered_mask = record['unfiltered_mask']
        unfiltered_areas = record['unfiltered_areas']
        pb_mask = record['pb_mask']

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(file_name)

        im = axes[0, 0].imshow(original_image)
        fig.colorbar(im, ax=axes[0, 0], fraction=0.046, pad=0.04)
        axes[0, 0].set_title('Original Image')
        axes[0, 0].axis('off')

        axes[0, 1].imshow(unfiltered_mask, cmap='gray')
        axes[0, 1].set_title('Unfiltered PB Mask')
        axes[0, 1].axis('off')

        axes[1, 0].hist(unfiltered_areas, bins=50)
        axes[1, 0].set_title('Unfiltered PB Size Histogram')
        axes[1, 0].set_xlabel('Area (px)')
        axes[1, 0].set_ylabel('Count')

        axes[1, 1].imshow(pb_mask, cmap='gray')
        axes[1, 1].set_title('Filtered PB Mask')
        axes[1, 1].axis('off')

        fig.tight_layout()
        qc_name = Path(file_name).stem + '_QC.png'
        fig.savefig(visualization_output_folder / qc_name, dpi=300, bbox_inches='tight')
        plt.close(fig)


def save_component_visualizations(visualization_output_folder, visualization_records):
    """Write separate pre-filter component figures for soft-threshold and DoG masks."""
    for record in visualization_records:
        file_name = record['file_name']
        soft_mask = record['soft_mask']
        dog_mask = record['dog_mask']

        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        fig.suptitle(file_name)

        axes[0].imshow(soft_mask, cmap='gray')
        axes[0].set_title('Soft Threshold Mask')
        axes[0].axis('off')

        axes[1].imshow(dog_mask, cmap='gray')
        axes[1].set_title('DoG Mask')
        axes[1].axis('off')

        fig.tight_layout()
        component_name = Path(file_name).stem + '_MaskComponents_QC.png'
        fig.savefig(visualization_output_folder / component_name, dpi=300, bbox_inches='tight')
        plt.close(fig)

def save_condition_unfiltered_histogram(condition_folder, visualization_output_folder, visualization_records):
    """Write one pooled unfiltered PB size histogram for the whole condition."""
    all_unfiltered_areas = []
    for record in visualization_records:
        all_unfiltered_areas.extend(record['unfiltered_areas'])

    fig, ax = plt.subplots(figsize=(8, 6))

    if all_unfiltered_areas:
        max_area = int(np.max(all_unfiltered_areas))
        bins = np.arange(0.5, max_area + 1.5, 1)
        ax.hist(all_unfiltered_areas, bins=bins, edgecolor='black', linewidth=0.5)
    else:
        ax.text(
            0.5,
            0.5,
            'No unfiltered PB objects detected',
            ha='center',
            va='center',
            transform=ax.transAxes,
        )

    ax.set_title(f'{condition_folder.name} Unfiltered PB Size Histogram')
    ax.set_xlabel('Area (px)')
    ax.set_ylabel('Count')

    fig.tight_layout()
    fig.savefig(
        visualization_output_folder / f'{condition_folder.name}_Unfiltered_PB_Size_Histogram.png',
        dpi=300,
        bbox_inches='tight',
    )
    plt.close(fig)

def save_pbodies_per_cell_summary(base_folder, counts_by_condition):
    """Save one combined bar+points plot of PBodies per cell across conditions."""
    summary_output_folder = Path(base_folder) / 'Summary_Results'
    summary_output_folder.mkdir(parents=True, exist_ok=True)

    condition_names = list(counts_by_condition.keys())
    x_positions = np.arange(len(condition_names))
    means = [
        np.mean(counts_by_condition[condition]) if counts_by_condition[condition] else 0
        for condition in condition_names
    ]

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x_positions, means, color='lightgray', edgecolor='black', width=0.6)

    rng = np.random.default_rng(42)
    for i, condition in enumerate(condition_names):
        counts = counts_by_condition[condition]
        if not counts:
            continue

        jitter = rng.uniform(-0.12, 0.12, size=len(counts))
        ax.scatter(
            np.full(len(counts), x_positions[i]) + jitter,
            counts,
            color='black',
            s=35,
            alpha=0.8,
            zorder=3,
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(condition_names, rotation=45, ha='right')
    ax.set_ylabel('PBodies per Cell')
    ax.set_title('PBodies per Cell by Condition')
    fig.tight_layout()

    fig.savefig(
        summary_output_folder / 'PBodies_per_Cell_by_Condition.png',
        dpi=300,
        bbox_inches='tight',
    )
    plt.close(fig)

def image_normalization(image, alpha):
    image_pixels = image.ravel()
    total_pixels = image_pixels.shape[0]
    counts, bin_edges = np.histogram(
        image_pixels,
        bins=256,
        range=(0, max(image_pixels))  # fixed range
    )
    min_val = None
    for i in range(len(counts)):
        if counts[i] > total_pixels * 0.1:
            continue
        elif counts[i] > total_pixels/alpha:
            min_val = bin_edges[i]
            break
    
    max_val = None
    for i in range(len(counts)-1,-1, -1):
        if counts[i] > total_pixels * 0.1:
            continue
        elif counts[i] > total_pixels/alpha:
            max_val = bin_edges[i]
            break

    if min_val is None or max_val is None or max_val <= min_val:
        print("no normalization")
        return image.astype(np.float32)
    
    image = image.astype(np.float32)
    image = (image - min_val) / (max_val - min_val)
    image = np.clip(image, 0, 1)
    image = (image * 65535).astype(np.uint16)
    
    return image

def process_condition(
    condition_folder,
    pb_channel,
    nucleus_mask_token,
    cell_mask_token,
    soft_threshold_fraction,
    low_sigma,
    high_sigma,
    dog_threshold_scale,
    min_pb_area,
    max_pb_area,
    min_circularity,
    alpha,
    minimum_intensity
):
    """Generate PB masks for one condition folder."""
    matched_files = collect_condition_files(
        condition_folder=condition_folder,
        pb_channel=pb_channel,
        nucleus_mask_token=nucleus_mask_token,
        cell_mask_token=cell_mask_token,
    )

    mask_output_folder = condition_folder / 'Pb_Masks'
    visualization_output_folder = condition_folder / 'PB_Visualization'
    mask_output_folder.mkdir(parents=True, exist_ok=True)
    visualization_output_folder.mkdir(parents=True, exist_ok=True)

    visualization_records = []
    pb_counts_per_cell = []

    for matched in matched_files:
        pb_image_path = matched['pb_image_path']
        cell_mask_path = matched['cell_mask_path']
        nucleus_mask_path = matched['nucleus_mask_path']

        original_stack = tf.imread(pb_image_path)
        original_image = max_project_to_2d(original_stack)

        cell_mask = load_mask_2d(cell_mask_path)
        nucleus_mask = load_mask_2d(nucleus_mask_path)

        if cell_mask.shape != original_image.shape:
            raise ValueError(
                f'Cell mask shape {cell_mask.shape} does not match PB image shape {original_image.shape} '
                f'for {pb_image_path.name}.'
            )
        if nucleus_mask.shape != original_image.shape:
            raise ValueError(
                f'Nucleus mask shape {nucleus_mask.shape} does not match PB image shape {original_image.shape} '
                f'for {pb_image_path.name}.'
            )

        allowed_mask = cell_mask & ~nucleus_mask

        percentile = np.percentile(original_image[allowed_mask], 50)
        image_masked = original_image.copy()
        image_masked[~allowed_mask] = percentile

        adjusted_image = image_normalization(image_masked, alpha)

        image_max = adjusted_image.max()
        soft_threshold = soft_threshold_fraction * image_max 

        soft_mask = seg.thresholding(adjusted_image, soft_threshold)

        dog_image = difference_of_gaussians(adjusted_image, low_sigma=low_sigma, high_sigma=high_sigma)
        dog_squared = dog_image ** 2
        
        """"
        raw_threshold_dog = get_histogramm_highest_varation_value(dog_squared)
        threshold_dog = raw_threshold_dog * dog_threshold_scale

        if np.issubdtype(dog_squared.dtype, np.floating):
            histogram_bins = 100
        else:
            histogram_bins = max(int(dog_squared.max() - dog_squared.min()), 1)

        dog_values = dog_squared.ravel()
        count, values = np.histogram(dog_values, bins=histogram_bins)
        gradient = np.gradient(count)
        gradient[gradient < 0] = 0
        bin_centers = 0.5 * (values[:-1] + values[1:])

        plot_count = count[4:]
        plot_left_edges = values[4:-1]
        plot_widths = np.diff(values)[4:]
        plot_centers = bin_centers[4:]
        plot_gradient = gradient[4:]

        fig, ax1 = plt.subplots(figsize=(8, 5))
        ax1.bar(
            plot_left_edges,
            plot_count,
            width=plot_widths,
            align='edge',
            color='gray',
            alpha=0.8,
            edgecolor='black',
            linewidth=0.5,
        )
        ax1.axvline(raw_threshold_dog, color='red', linestyle='--', linewidth=2, label='Max positive derivative')
        ax1.axvline(threshold_dog, color='orange', linestyle=':', linewidth=2, label='Applied DoG threshold')
        ax1.set_title(f'DoG Squared Histogram: {pb_image_path.stem}')
        ax1.set_xlabel('DoG squared value')
        ax1.set_ylabel('Count')

        ax2 = ax1.twinx()
        ax2.plot(plot_centers, plot_gradient, color='blue', linewidth=1.5, label='Positive histogram derivative')
        ax2.set_ylabel('Derivative')

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

        fig.tight_layout()
        hist_name = pb_image_path.stem + '_DoG_squared_histogram.png'
        fig.savefig(visualization_output_folder / hist_name, dpi=300, bbox_inches='tight')
        plt.close(fig)
        """       
        #dog_squared_name = pb_image_path.stem + '_DoG_squared.tif'
        #tiff.imwrite(
        #visualization_output_folder / dog_squared_name,
        #dog_squared.astype(np.float32),
        #)

        threshold_dog = get_histogramm_highest_varation_value(dog_squared)
        threshold_dog *= dog_threshold_scale

        dog_mask = seg.thresholding(dog_squared, threshold_dog)

        constrained_soft_mask = np.asarray(soft_mask & allowed_mask, dtype=np.uint8)
        constrained_dog_mask = np.asarray(dog_mask & allowed_mask, dtype=np.uint8)

        constrained_mask = (constrained_dog_mask > 0) | (constrained_soft_mask > 0)
        unfiltered_mask = constrained_mask.astype(np.uint8)
        unfiltered_areas = get_component_areas(constrained_mask)

        filtered_mask = filter_pb_components(
            constrained_mask,
            min_pb_area,
            max_pb_area,
            min_circularity,
            original_image,
            minimum_intensity
        )
        pb_mask = filtered_mask.astype(np.uint8)
        pb_count = len(get_component_areas(filtered_mask))
        pb_counts_per_cell.append(pb_count)

        mask_name = pb_image_path.stem + '_PBmask.tif'
        tiff.imwrite(mask_output_folder / mask_name, pb_mask.astype(np.uint8))

        visualization_records.append(
            {
                'file_name': pb_image_path.name,
                'original_image': adjusted_image,
                'soft_mask': constrained_soft_mask,
                'dog_mask': constrained_dog_mask,
                'unfiltered_mask': unfiltered_mask,
                'unfiltered_areas': unfiltered_areas,
                'pb_mask': pb_mask,
            }
        )

    if not visualization_records:
        raise FileNotFoundError(f'No matching PB image and mask sets were found in {condition_folder}')

    save_visualizations(visualization_output_folder, visualization_records)
    save_component_visualizations(visualization_output_folder, visualization_records)
    save_condition_unfiltered_histogram(condition_folder, visualization_output_folder, visualization_records)
    print(
        f'Saved {len(visualization_records)} PB masks to {mask_output_folder} '
        f'and QC figures to {visualization_output_folder}')
    return pb_counts_per_cell

def process_conditions(
    base_folder,
    conditions,
    pb_channel,
    nucleus_mask_token,
    cell_mask_token,
    soft_threshold_fraction,
    low_sigma,
    high_sigma,
    dog_threshold_scale,
    min_pb_area,
    max_pb_area,
    min_circularity,
    alpha,
    minimum_intensity,
):
    
    """Loop through user-selected conditions and generate PB masks for each."""
    bigfish_root = Path(base_folder)
    counts_by_condition = {}

    for condition in conditions:
        condition_folder = bigfish_root / condition
        if not condition_folder.exists():
            raise FileNotFoundError(f'Condition folder does not exist: {condition_folder}')

        pb_counts_per_cell = process_condition(
            condition_folder=condition_folder,
            pb_channel=pb_channel,
            nucleus_mask_token=nucleus_mask_token,
            cell_mask_token=cell_mask_token,
            soft_threshold_fraction=soft_threshold_fraction,
            low_sigma=low_sigma,
            high_sigma=high_sigma,
            dog_threshold_scale=dog_threshold_scale,
            min_pb_area=min_pb_area,
            max_pb_area=max_pb_area,
            min_circularity=min_circularity,
            alpha=alpha,
            minimum_intensity=minimum_intensity,
        )

        counts_by_condition[condition] = pb_counts_per_cell

    save_pbodies_per_cell_summary(bigfish_root, counts_by_condition)

"""Editable entry point."""
base_folder = Path(r"\\10.67.15.232\lab\Morgane\MV635-...-642_26-03-12_IF NBDY Dox+ - _ Dcpa1a 647 LSM14A 488\bigFISH_In")
conditions = ['Dox+_CT','Dox+_HS','Dox+_2HR', 'Dox+_3HR', 'Dox+_4HR', 'Dox+_6HR', 'Dox+_8HR', 'Dox+_10HR', 'Dox-_CT', 'Dox-_HS', 'Dox-_2HR', 'Dox-_3HR', 'Dox-_4HR', 'Dox-_6HR', 'Dox-_8HR', 'Dox-_10HR']
#conditions = ['Ctrl', 'HS', '2h', '3h', '4h', '6h', '8h', '10h', '12h', '24h']
process_conditions(
    base_folder=base_folder,
    conditions=conditions,
    pb_channel='_640S_',
    nucleus_mask_token='__nuc_',
    cell_mask_token='__print_',
    soft_threshold_fraction=0.95,
    low_sigma=1,
    high_sigma=5,
    dog_threshold_scale=0.45,
    min_pb_area=15,
    max_pb_area=1000,
    min_circularity=0.6,
    alpha=50000,
    minimum_intensity = 1000,
)

