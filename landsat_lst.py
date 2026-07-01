#!/usr/bin/env python3
"""Download Landsat 8/9 Collection 2 Level-2 imagery from Microsoft Planetary
Computer, derive Land Surface Temperature (LST) for a bounding box, and export
a publication-ready JPEG map.
"""

import argparse
import math
import sys
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import planetary_computer
import pystac_client
import rasterio
from matplotlib.colors import Normalize
from matplotlib_scalebar.scalebar import ScaleBar
from rasterio.transform import array_bounds
from rasterio.warp import transform_bounds
from rasterio.windows import from_bounds
from shapely.geometry import box, shape

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "landsat-c2-l2"

# ST_B10 scale/offset per USGS Landsat Collection 2 Level-2 product guide.
ST_SCALE, ST_OFFSET = 0.00341802, 149.0
ST_QA_SCALE = 0.01

# QA_PIXEL bit flags (fill, dilated cloud, cirrus, cloud, cloud shadow).
BAD_QA_BITS = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4)


@dataclass
class LSTResult:
    lst_celsius: np.ma.MaskedArray
    transform: rasterio.Affine
    crs: rasterio.crs.CRS
    valid_fraction: float
    mean_uncertainty_k: float
    item_id: str
    datetime: str
    cloud_cover: float


def find_best_scene(bbox, max_cloud_cover=1.0, lookback_days=None, platforms=("landsat-8", "landsat-9")):
    """Return the most recent, fully-covering Landsat item under a cloud-cover
    threshold. Falls back to the lowest cloud-cover fully-covering scene if
    nothing qualifies within the threshold."""
    aoi = box(*bbox)
    catalog = pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)

    search_kwargs = dict(
        collections=[COLLECTION],
        bbox=bbox,
        query={
            "platform": {"in": list(platforms)},
            "landsat:collection_category": {"in": ["T1"]},
        },
        limit=100,
    )
    if lookback_days:
        search_kwargs["datetime"] = f"-P{lookback_days}D/.."

    items = list(catalog.search(**search_kwargs).item_collection())
    if not items:
        raise RuntimeError("No Landsat scenes found for this bounding box.")

    candidates = [it for it in items if shape(it.geometry).contains(aoi)]
    if not candidates:
        raise RuntimeError("No single scene fully covers the requested bounding box.")

    clear = [it for it in candidates if it.properties["eo:cloud_cover"] <= max_cloud_cover]
    if clear:
        clear.sort(key=_ts, reverse=True)  # most recent clear-sky scene wins
        return clear[0]

    candidates.sort(key=lambda it: (it.properties["eo:cloud_cover"], -_ts(it)))  # fall back to least cloudy
    return candidates[0]


def _ts(item):
    from datetime import datetime

    return datetime.fromisoformat(item.properties["datetime"].replace("Z", "+00:00")).timestamp()


def compute_lst(item, bbox):
    """Download the ST_B10, QA_PIXEL and ST-uncertainty bands clipped to bbox
    and return cloud-masked LST in Celsius."""
    lwir_href = item.assets["lwir11"].href
    qa_pixel_href = item.assets["qa_pixel"].href
    st_qa_href = item.assets["qa"].href

    with rasterio.open(lwir_href) as src:
        dst_crs = src.crs
        bounds_proj = transform_bounds("EPSG:4326", dst_crs, *bbox)
        window = from_bounds(*bounds_proj, transform=src.transform)
        st_dn = src.read(1, window=window).astype("float64")
        win_transform = src.window_transform(window)

    with rasterio.open(qa_pixel_href) as src:
        qa_pixel = src.read(1, window=from_bounds(*bounds_proj, transform=src.transform))

    with rasterio.open(st_qa_href) as src:
        st_unc_dn = src.read(1, window=from_bounds(*bounds_proj, transform=src.transform)).astype("float64")

    lst_celsius = (st_dn * ST_SCALE + ST_OFFSET) - 273.15
    st_uncertainty_k = st_unc_dn * ST_QA_SCALE

    bad_mask = (qa_pixel & BAD_QA_BITS).astype(bool)
    valid_fraction = 1 - bad_mask.mean()
    lst_masked = np.ma.masked_where(bad_mask, lst_celsius)

    return LSTResult(
        lst_celsius=lst_masked,
        transform=win_transform,
        crs=dst_crs,
        valid_fraction=valid_fraction,
        mean_uncertainty_k=float(np.nanmean(st_uncertainty_k[~bad_mask])),
        item_id=item.id,
        datetime=item.properties["datetime"],
        cloud_cover=item.properties["eo:cloud_cover"],
    )


def plot_lst(result: LSTResult, out_path, region_name="Study Area"):
    """Render a publication-ready LST map and save it as a JPEG."""
    lst = result.lst_celsius.filled(np.nan)
    h, w = lst.shape
    bounds_proj = array_bounds(h, w, result.transform)
    left, bottom, right, top = transform_bounds(result.crs, "EPSG:4326", *bounds_proj)

    fig, ax = plt.subplots(figsize=(10, 8.5), dpi=300)

    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad(color="#dddddd")
    vmin, vmax = np.nanpercentile(lst, [1, 99])
    norm = Normalize(vmin=vmin, vmax=vmax)

    im = ax.imshow(lst, cmap=cmap, norm=norm, extent=[left, right, bottom, top], interpolation="nearest")

    ax.set_xlabel("Longitude", fontsize=11)
    ax.set_ylabel("Latitude", fontsize=11)
    ax.set_title(f"Land Surface Temperature\n{region_name}", fontsize=15, fontweight="bold", pad=12)

    acq_date = result.datetime[:10]
    ax.text(
        0.5, 1.005,
        f"Landsat OLI/TIRS Collection 2 Level-2 (ST_B10) — {acq_date} — Cloud cover: {result.cloud_cover:.2f}%",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=9, style="italic",
    )

    cbar = fig.colorbar(im, ax=ax, orientation="vertical", fraction=0.04, pad=0.03)
    cbar.set_label("Land Surface Temperature (°C)", fontsize=11)
    ax.set_aspect("equal")

    lat_mid = (bottom + top) / 2
    m_per_deg_lon = 111320 * math.cos(math.radians(lat_mid))
    ax.add_artist(ScaleBar(m_per_deg_lon, units="m", location="lower right", box_alpha=0.7))

    ax.annotate(
        "N", xy=(0.96, 0.90), xytext=(0.96, 0.80), xycoords="axes fraction",
        arrowprops=dict(facecolor="black", width=3, headwidth=10, headlength=8),
        ha="center", va="center", fontsize=12, fontweight="bold",
    )

    fig.text(
        0.01, 0.01,
        "Data: USGS/NASA Landsat Collection 2 Level-2, served via Microsoft Planetary Computer\n"
        f"Scene: {result.item_id} | CRS: EPSG:4326 (display) | Cloud-masked (QA_PIXEL)",
        fontsize=6.5, color="#444444",
    )

    plt.tight_layout(rect=[0, 0.02, 1, 1])
    fig.savefig(out_path, dpi=300, format="jpg", pil_kwargs={"quality": 95}, bbox_inches="tight")
    plt.close(fig)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bbox", nargs=4, type=float, required=True, metavar=("MINX", "MINY", "MAXX", "MAXY"),
                   help="Bounding box in WGS84 lon/lat: minx miny maxx maxy")
    p.add_argument("--out", default="lst_output.jpg", help="Output JPEG path")
    p.add_argument("--region-name", default="Study Area", help="Region label shown in the map title")
    p.add_argument("--max-cloud-cover", type=float, default=1.0,
                   help="Preferred cloud-cover threshold (%%) for scene selection (default: 1.0)")
    p.add_argument("--lookback-days", type=int, default=None,
                   help="Only consider scenes from the last N days (default: unlimited)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    bbox = tuple(args.bbox)

    print(f"Searching Planetary Computer STAC for scenes covering {bbox} ...")
    item = find_best_scene(bbox, max_cloud_cover=args.max_cloud_cover, lookback_days=args.lookback_days)
    print(f"Selected scene: {item.id} ({item.properties['datetime'][:10]}, "
          f"cloud cover {item.properties['eo:cloud_cover']:.2f}%)")

    print("Downloading thermal + QA bands and computing LST ...")
    result = compute_lst(item, bbox)
    valid_lst = result.lst_celsius.compressed()
    print(f"Valid (clear-sky) pixel fraction: {result.valid_fraction:.2%}")
    print(f"LST range: {valid_lst.min():.2f} to {valid_lst.max():.2f} °C "
          f"(mean {valid_lst.mean():.2f} °C)")
    print(f"Mean ST uncertainty: {result.mean_uncertainty_k:.3f} K")

    print(f"Rendering map to {args.out} ...")
    plot_lst(result, args.out, region_name=args.region_name)
    print("Done.")


if __name__ == "__main__":
    sys.exit(main())
