"""
Enhanced Statistics Module

A Python program to compute statistics of Euler deconvolution estimates
both globally and within grid windows.

This enhanced version includes:
1. Original classic function for global statistics
2. New window_stats function for grid-based local statistics

authors: Felipe F. Melo and Valeria C.F. Barbosa, 2019 (original)
Enhanced with window-based analysis
"""

import numpy as np
import os


def classic(est_classic, area_plt, SI_vet, name, path):
    """
    Original function - computes global statistics for the entire grid
    """
    results = []

    for i in range(len(est_classic)):

        estimates = np.stack(
            (
                est_classic[i][:, 0],
                est_classic[i][:, 1],
                est_classic[i][:, 2],
                est_classic[i][:, 3],
            ),
            axis=-1,
        )

        masked = np.ma.array(
            estimates,
            mask=np.repeat(estimates[:, 0] <= area_plt[0], estimates.shape[1]),
        )
        masked = np.ma.array(
            masked, mask=np.repeat(masked[:, 0] >= area_plt[1], estimates.shape[1])
        )
        masked = np.ma.array(
            masked, mask=np.repeat(masked[:, 1] <= area_plt[2], estimates.shape[1])
        )
        masked = np.ma.array(
            masked, mask=np.repeat(masked[:, 1] >= area_plt[3], estimates.shape[1])
        )

        meanx = np.mean(masked[:, 0] / 1000.0)
        meany = np.mean(masked[:, 1] / 1000.0)
        meanz = np.mean(masked[:, 2] / 1000.0)
        results.append([SI_vet[i], meanx, meany, meanz])

    output = np.array([(results[i]) for i in range(0, len(SI_vet))])
    np.savetxt(
        path + "/" + str(name) + ".txt",
        output,
        fmt="%.3f",
        header="SI, mean x, mean y, mean z",
        comments="",
    )
    return


def _window_index(values, origin, size, n_windows):
    """
    Index of the window holding each value, the window i spanning
    origin + i * size <= value < origin + (i + 1) * size (the bounds computed
    exactly as written, so values on a boundary fall where a direct test of
    each window puts them), and whether the value lies in any window.
    """
    finite = np.isfinite(values)
    v = np.where(finite, values, origin)
    i = np.floor((v - origin) / size).astype(np.int64)
    i = np.where(v < origin + i * size, i - 1, i)
    i = np.where(v >= origin + (i + 1) * size, i + 1, i)
    return i, finite & (i >= 0) & (i < n_windows)


def window_stats(
    est_classic,
    area_plt,
    SI_vet,
    name,
    path,
    data_shape,
    window_size,
    detailed_stats=True,
):
    """
    Compute statistics for each window across the grid

    Parameters:
    * est_classic: list of arrays with estimates for each SI
    * area_plt: [south, north, west, east] - area bounds
    * SI_vet: list of structural indices
    * name: base name for output files
    * path: output directory path
    * data_shape: (rows, cols) - shape of original grid
    * window_size: size of analysis windows
    * detailed_stats: if True, saves detailed statistics; if False, saves only means

    Returns:
    * Dictionary with statistics for each SI
    """

    # Calculate grid parameters
    rows, cols = data_shape
    south, north, west, east = area_plt

    # Calculate coordinate ranges
    x_range = east - west
    y_range = north - south

    # Calculate number of windows that fit in each direction
    n_windows_x = cols // window_size
    n_windows_y = rows // window_size

    # Calculate actual window size in coordinate units
    window_size_x = x_range / n_windows_x
    window_size_y = y_range / n_windows_y

    print(f"Grid analysis: {n_windows_y} x {n_windows_x} windows")
    print(f"Window size: {window_size_x:.2f} x {window_size_y:.2f} coordinate units")

    # Process each SI
    all_results = {}

    for si_idx, SI in enumerate(SI_vet):
        estimates = est_classic[si_idx]

        # Assign each estimate to its window once (rather than testing every
        # estimate against every window), then gather the estimates of each
        # window with a sort, keeping windows in row-major order
        ix, ok_x = _window_index(estimates[:, 1], west, window_size_x, n_windows_x)
        iy, ok_y = _window_index(estimates[:, 0], south, window_size_y, n_windows_y)
        ok = ok_x & ok_y
        key = iy[ok] * n_windows_x + ix[ok]
        order = np.argsort(key, kind="stable")
        key = key[order]
        window_estimates = estimates[ok][order][:, :4]

        windows, start, n_points = np.unique(
            key, return_index=True, return_counts=True
        )
        row = windows // n_windows_x
        col = windows % n_windows_x
        # window bounds and centres, from the same expressions as the bounds
        win_west = west + col * window_size_x
        win_east = west + (col + 1) * window_size_x
        win_south = south + row * window_size_y
        win_north = south + (row + 1) * window_size_y
        win_center_x = (win_west + win_east) / 2.0
        win_center_y = (win_south + win_north) / 2.0

        if len(windows):
            means = np.add.reduceat(window_estimates, start, axis=0) / n_points[:, None]
        else:
            means = np.zeros((0, 4))
        mean_x, mean_y, mean_z, mean_base = means.T
        columns = [
            row, col, win_center_x, win_center_y, n_points,
            mean_x, mean_y, mean_z / 1000.0, mean_base,
        ]

        if detailed_stats:
            if len(windows):
                deviations = window_estimates - np.repeat(means, n_points, axis=0)
                stds = np.sqrt(
                    np.add.reduceat(deviations**2, start, axis=0) / n_points[:, None]
                )
                min_z = np.minimum.reduceat(window_estimates[:, 2], start)
                max_z = np.maximum.reduceat(window_estimates[:, 2], start)
            else:
                stds = np.zeros((0, 4))
                min_z = max_z = np.zeros(0)
            # Single point - no standard deviation possible
            stds[n_points == 1] = 0.0
            std_x, std_y, std_z, std_base = stds.T
            columns += [std_x, std_y, std_z / 1000.0, std_base,
                        min_z / 1000.0, max_z / 1000.0]

        # only windows with estimates are kept
        window_results = np.column_stack(columns).astype(float)
        all_results[SI] = window_results

        # Save results for this SI
        if detailed_stats:
            header = (
                "window_row, window_col, center_x, center_y, n_estimates, "
                "mean_y, mean_x, mean_depth_km, mean_base_level, "
                "std_x, std_y, std_depth_km, std_base_level, "
                "min_depth_km, max_depth_km"
            )
        else:
            header = (
                "window_row, window_col, center_x, center_y, n_estimates, "
                "mean_y, mean_x, mean_depth_km, mean_base_level"
            )
        SI_name = SI
        if SI_name == 0.001:
            SI_name = 0

        output_filename = f"{path}/{name}_window_stats_SI_{SI_name}.txt"
        np.savetxt(
            output_filename,
            window_results,
            fmt="%.6f",
            header=header,
            comments="",
            delimiter=",",
        )

        print(f"Saved window statistics for SI={SI} to {output_filename}")

    # Create summary file with statistics across all windows
    summary_results = []
    for si_idx, SI in enumerate(SI_vet):
        window_data = all_results[SI]
        valid_windows = window_data[window_data[:, 4] > 0]  # Windows with estimates

        if len(valid_windows) > 0:
            total_estimates = np.sum(valid_windows[:, 4])
            mean_depth = np.nanmean(valid_windows[:, 7])  # mean depth across windows
            std_depth = np.nanstd(valid_windows[:, 7])  # std of window means
            n_windows_with_data = len(valid_windows)

            summary_results.append(
                [
                    SI,
                    n_windows_with_data,
                    n_windows_x * n_windows_y,
                    total_estimates,
                    mean_depth,
                    std_depth,
                ]
            )
        else:
            summary_results.append(
                [SI, 0, n_windows_x * n_windows_y, 0, np.nan, np.nan]
            )

    summary_output = np.array(summary_results)
    summary_filename = f"{path}/{name}_window_summary.txt"
    np.savetxt(
        summary_filename,
        summary_output,
        fmt="%.6f",
        header="SI, windows_with_data, total_windows, total_estimates, mean_depth_km, std_depth_km",
        comments="",
        delimiter=",",
    )

    print(f"Saved window summary to {summary_filename}")
    return all_results


def enhanced_analysis(
    est_classic, area_plt, SI_vet, name, path, data_shape, window_size
):
    """
    Wrapper function that performs both global and window-based analysis

    Parameters:
    * est_classic: list of arrays with estimates for each SI
    * area_plt: [south, north, west, east] - area bounds
    * SI_vet: list of structural indices
    * name: base name for output files
    * path: output directory path
    * data_shape: (rows, cols) - shape of original grid
    * window_size: size of analysis windows
    """

    print("=== Enhanced Euler Deconvolution Analysis ===")

    # Run original global analysis
    print("\n1. Computing global statistics...")
    classic(est_classic, area_plt, SI_vet, name, path)
    print("Global statistics saved.")

    # Run new window-based analysis
    print("\n2. Computing window-based statistics...")
    window_results = window_stats(
        est_classic,
        area_plt,
        SI_vet,
        name,
        path,
        data_shape,
        window_size,
        detailed_stats=True,
    )

    print("\n=== Analysis Complete ===")
    print(f"Results saved to: {path}")
    print(f"Files created:")
    print(f"  - {name}.txt (global statistics)")
    print(f"  - {name}_window_stats_SI_*.txt (window statistics for each SI)")
    print(f"  - {name}_window_summary.txt (summary across all windows)")

    return window_results
