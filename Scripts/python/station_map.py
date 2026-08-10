#!/usr/bin/env python3
"""Interactive maps of the study reservoirs and the water quality sampling locations.

Builds three maps and a print-quality PNG of each:

  station_map            both projects, for orientation
  station_map_lanier     Lake Sidney Lanier and Buford Dam
  station_map_allatoona  Lake Allatoona and Allatoona Dam

Each map is written to ``analysis/maps/`` as two interactive ``.html`` files, one
carrying the title and legend panels and ``<name>_plain.html`` without them for
screen capture, and two headless-Chrome stills: ``<name>.png`` on the default
topographic basemap and ``<name>_satellite.png`` on Esri aerial imagery. A third
capture, built separately with enlarged labels and no panels, is written straight to
``report/latex/src/images/`` for inclusion at full text width. Sixteen basemaps are
offered in the layer selector of the HTML, so only the still images are duplicated.
Label colours follow the basemap, inverting to light text with a dark halo over
imagery and dark cartography.

Coordinates are read from ``Data/station_coordinates.csv`` and
``Data/river_label_anchors.csv`` and are not derived here. Their provenance, carried
in those files, is:
  - Forebay profile stations: EPA Water Quality Portal Station service (WGS84,
    GPS-located).
  - Tailrace gages: USGS NWIS site service (NAD83). The coordinate method code is
    ``N`` for 02334430 and ``M`` for 02394000, so the Allatoona gage position is
    map-interpolated rather than surveyed and is the less precise of the two.
  - River labels are anchored on USGS gages named for the river being labelled, so
    each label sits on a surveyed point of that river rather than an estimated one.

No dam structure coordinate is plotted. The two USGS gages are the stations named
"at Buford Dam" and "at Allatoona Dam" and stand in for the dam locations; plotting
a separate dam point would require a coordinate this study has not sourced.

Both continuous tailrace records described in the report are co-located with the
USGS gage at their project (the USACE-SAM sonde at Allatoona and the Corps project
monitor at Buford, each verified against the gage temperature record), so each
tailrace marker represents every tailwater record at that project.

Run with the clearwater conda env:
    /opt/homebrew/Caskroom/miniconda/base/envs/clearwater/bin/python3 \
        Scripts/python/station_map.py
"""
import math
import os
import tempfile
import time
from typing import Any

import folium
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COORDS = os.path.join(REPO, "Data", "station_coordinates.csv")
RIVERS = os.path.join(REPO, "Data", "river_label_anchors.csv")
OUT_DIR = os.path.join(REPO, "analysis", "maps")
IMAGE_DIR = os.path.join(REPO, "report", "latex", "src", "images")

# Width in pixels of the copies placed in the report. The text block is 397.5 pt,
# or 5.52 in, so 1800 px is about 325 dpi at full text width: past the point where
# more pixels show on paper, and a fraction of the size of the full captures.
REPORT_WIDTH = 1800

# Frame shape used for the report captures. The overview is given more height than
# its on-screen counterpart: printed at text width it is a short figure, and the
# enlarged labels need vertical room they do not have in a 16:9 frame.
REPORT_SIZE = {"both": (1700, 1200), "Lanier": (1250, 1200), "Allatoona": (1400, 1100)}

# Target label size in points once the figure is reduced to text width. Nine points
# proved too large for the wide overview, whose labels then dominated the map.
REPORT_LABEL_PT = 7.5

# LaTeX text block width in points, the figure width the report places these captures
# at. This is the same width the comment on REPORT_WIDTH above rounds to 397.5 pt.
TEXT_WIDTH_PT = 397.485

# Base label size in CSS pixels, the size river_label and the station labels use at
# scale 1.0. The report label scale is the ratio that turns this into REPORT_LABEL_PT
# once the capture is reduced to the text block width.
BASE_LABEL_PX = 12.0

# Project colours match the report figures: Buford/Lanier blue, Allatoona orange.
PROJECT_COLOR = {"Lanier": "#4c72b0", "Allatoona": "#dd8452"}
PROJECT_LABEL = {"Lanier": "Lake Sidney Lanier / Buford Dam",
                 "Allatoona": "Lake Allatoona / Allatoona Dam"}
RIVER_COLOR = "#1a5b8a"

TITLE = {
    "both": ("Water quality sampling locations",
             "Lake Sidney Lanier (Buford Dam, Chattahoochee River) and "
             "Lake Allatoona (Allatoona Dam, Etowah River), Georgia"),
    "Lanier": ("Lake Sidney Lanier and Buford Dam",
               "Chattahoochee River, Apalachicola-Chattahoochee-Flint basin. "
               "Fed by the Chattahoochee and Chestatee Rivers."),
    "Allatoona": ("Lake Allatoona and Allatoona Dam",
                  "Etowah River, Alabama-Coosa-Tallapoosa basin. "
                  "Fed by the Etowah and Little Rivers."),
}

# Record detail shown in each popup, taken from the report data inventory (Table 1).
DETAIL = {
    "21GAEPD_WQX-LK_12_4028": [
        ("Operator", "Georgia EPD reservoir monitoring"),
        ("Parameters", "DO, temperature, conductivity (profiles); nutrient and "
                       "eutrophication suite (surface)"),
        ("Period of record", "1990-2025; depth-resolved 2000-2025"),
        ("Casts", "169; profiles to 47 m"),
    ],
    "21GAEPD_WQX-LK_14_4494": [
        ("Operator", "Georgia EPD reservoir monitoring"),
        ("Parameters", "DO, temperature, conductivity (profiles); nutrient and "
                       "eutrophication suite (surface)"),
        ("Period of record", "1990-2025; depth-resolved 2000-2025"),
        ("Casts", "144; profiles to 45 m"),
    ],
    "02334430": [
        ("Operator", "USGS, plus a co-located USACE project monitor"),
        ("Parameters", "DO, temperature, specific conductance"),
        ("Tailwater DO", "USACE project monitor, hourly 2002, 30-min 2004, 15-min "
                         "2005-2008; USGS daily, 2023-present"),
        ("Other record", "Discharge 1942-present; temperature 1975-present"),
    ],
    "02394000": [
        ("Operator", "USGS, plus a co-located USACE-SAM sonde"),
        ("Parameters", "DO, temperature, specific conductance, pH"),
        ("Tailwater DO", "USGS daily, 2005-01 to 2007-01; USACE-SAM sonde, 30-min, "
                         "2012-2013; grab samples 2019-2020"),
        ("Other record", "Discharge 1938-present; temperature 2005-present"),
    ],
}

TYPE_LABEL = {"forebay_profile": "Forebay profile station (source water)",
              "tailrace": "Tailrace station (release and tailwater)"}
TYPE_ICON = {"forebay_profile": "tint", "tailrace": "signal"}

# Basemaps offered in the layer selector, as (display name, xyzservices provider,
# dark background). Every entry was checked against a tile covering the study area
# and returns imagery; Esri OceanBasemap is deliberately absent because it answers
# HTTP 500 for inland tiles. The dark flag drives label colour: light text with a
# dark halo over photography and dark cartography, dark text with a white halo
# otherwise.
BASEMAPS = [
    ("Topographic - Esri NatGeoWorldMap",       "Esri.NatGeoWorldMap",    False),
    ("Topographic - Esri WorldTopoMap",         "Esri.WorldTopoMap",      False),
    ("Topographic - OpenTopoMap",               "OpenTopoMap",            False),
    ("Topographic - USGS USTopo",               "USGS.USTopo",            False),
    ("Satellite - Esri WorldImagery",           "Esri.WorldImagery",      True),
    ("Satellite - USGS USImagery",              "USGS.USImagery",         True),
    ("Satellite with labels - USGS USImageryTopo", "USGS.USImageryTopo",  True),
    ("Street - OpenStreetMap",                  "OpenStreetMap.Mapnik",   False),
    ("Street - Esri WorldStreetMap",            "Esri.WorldStreetMap",    False),
    ("Street - CartoDB Voyager",                "CartoDB.Voyager",        False),
    ("Light - CartoDB Positron",                "CartoDB.Positron",       False),
    ("Light - Esri WorldGrayCanvas",            "Esri.WorldGrayCanvas",   False),
    ("Dark - CartoDB DarkMatter",               "CartoDB.DarkMatter",     True),
    ("Terrain - Esri WorldTerrain",             "Esri.WorldTerrain",      False),
    ("Terrain - Esri WorldShadedRelief",        "Esri.WorldShadedRelief", False),
    ("Terrain - Esri WorldPhysical",            "Esri.WorldPhysical",     False),
]
DEFAULT_BASEMAP = "Topographic - Esri NatGeoWorldMap"
SATELLITE_BASEMAP = "Satellite - Esri WorldImagery"
BASEMAP_DARK = {name: dark for name, _, dark in BASEMAPS}


def tile_layer(name: str, provider_path: str, show: bool) -> folium.TileLayer:
    """Build a folium TileLayer from an xyzservices provider path.

    Several of these services stop at a coarse native zoom (Esri WorldPhysical at 8,
    the Esri terrain services at 13). Passing max_native_zoom lets Leaflet upsample
    their last available tile instead of showing blank tiles when the map is zoomed
    in past that limit.

    Only the active basemap is shown on load. Every layer added to the map
    contributes its attribution to the attribution control, so adding all sixteen at
    once would bury the map under a block of credit text; the inactive ones stay
    registered with the layer control and attach when selected.

    Parameters
    ----------
    name : str
        Display name in the layer selector.
    provider_path : str
        Dotted path into ``xyzservices.providers``, for example
        ``Esri.NatGeoWorldMap``.
    show : bool
        Whether this layer is active when the map loads.

    Returns
    -------
    folium.TileLayer
        Layer with ``max_zoom`` 20 and ``max_native_zoom`` set from the provider's own
        declared maximum, or unset when the provider does not declare one.
    """
    import xyzservices.providers as xp
    provider = xp
    for part in provider_path.split("."):
        provider = provider[part]
    native = provider.get("max_zoom")
    return folium.TileLayer(tiles=provider, name=name, overlay=False, control=True,
                            show=show, max_zoom=20,
                            max_native_zoom=int(native) if native else None)

# Little River is a genuine Allatoona tributary, but its only USGS gage (02392500)
# lies about 15 km upstream of the Little River embayment, so a label anchored there
# reads as detached from the reservoir. It is named in the Allatoona subtitle instead
# of being pinned. Its row is retained in Data/river_label_anchors.csv.
UNLABELLED_RIVERS = {"Little River"}


def haversine_km(a: tuple[float, float], b: tuple[float, float]) -> float:
    """Great-circle distance in kilometres between two (lat, lon) pairs.

    Parameters
    ----------
    a, b : tuple of float
        Latitude and longitude in decimal degrees. The archived coordinates are WGS84
        for the profile stations and NAD83 for the gages; over this study area the two
        datums differ by well under a metre, which is far below the precision of the
        map-interpolated Allatoona gage position, so they are used together without
        transformation.

    Returns
    -------
    float
        Distance in km on a sphere.

    Notes
    -----
    The radius is the IUGG arithmetic mean radius R1 = (2a + b) / 3 for the WGS84
    ellipsoid, 6371.0088 km. A spherical approximation is used because these
    separations are of order 1 km: substituting the local radius of curvature along
    each line moves the two distances this script prints by about 1 m at Lanier and
    3 m at Allatoona, well inside the 0.01 km they are reported to.
    """
    (la1, lo1), (la2, lo2) = a, b
    r = 6371.0088   # km, IUGG mean radius R1 for WGS84
    p1, p2 = math.radians(la1), math.radians(la2)
    dp, dl = math.radians(la2 - la1), math.radians(lo2 - lo1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(h))


def popup_html(row: Any) -> folium.Popup:
    """Build the click-through popup for one station.

    The record detail comes from the DETAIL table above, keyed on station id, and the
    position line from the archived coordinate file, so the popup shows the coordinate
    together with the service it came from and its datum.

    Parameters
    ----------
    row : namedtuple
        One row of ``Data/station_coordinates.csv`` as produced by ``itertuples``.

    Returns
    -------
    folium.Popup
        Popup capped at 380 px wide.
    """
    rows = "".join(
        f"<tr><td style='padding:2px 8px 2px 0;vertical-align:top;white-space:nowrap'>"
        f"<b>{k}</b></td><td style='padding:2px 0'>{v}</td></tr>"
        for k, v in DETAIL.get(row.station_id, [])
    )
    return folium.Popup(folium.Html(
        f"<div style='font-family:Helvetica,Arial,sans-serif;font-size:12px;"
        f"min-width:300px'>"
        f"<div style='font-size:13px;font-weight:700;color:{PROJECT_COLOR[row.project]}'>"
        f"{row.station_id}</div>"
        f"<div style='margin:2px 0 6px 0'>{row.station_name.title()}</div>"
        f"<div style='margin-bottom:6px;font-style:italic'>{TYPE_LABEL[row.station_type]}</div>"
        f"<table style='border-collapse:collapse'>{rows}"
        f"<tr><td style='padding:2px 8px 2px 0;vertical-align:top'><b>Position</b></td>"
        f"<td style='padding:2px 0'>{row.latitude:.5f}, {row.longitude:.5f} "
        f"({row.datum})<br><span style='color:#666'>{row.coord_source}</span></td></tr>"
        f"</table></div>", script=True), max_width=380)


def halo(dark_background: bool) -> str:
    """CSS text shadow that keeps a label legible against its basemap.

    Dark text with a white halo reads best on the OpenStreetMap base; on aerial
    imagery the same label loses contrast, so the polarity is inverted to light text
    with a dark halo, which is the usual cartographic treatment over photography.

    Parameters
    ----------
    dark_background : bool
        True for aerial imagery and dark cartography, as flagged in BASEMAPS.

    Returns
    -------
    str
        A ``text-shadow`` declaration placing a 1.5 px halo on all four diagonals.
    """
    c = "#000" if dark_background else "#fff"
    return (f"text-shadow:-1.5px -1.5px 0 {c}, 1.5px -1.5px 0 {c},"
            f"-1.5px 1.5px 0 {c}, 1.5px 1.5px 0 {c}")


def river_label(row: Any, dark_background: bool = False,
                scale: float = 1.0) -> folium.Marker:
    """A map label for a named river, anchored on a USGS gage on that river.

    Inflow labels sit to the right of their anchor. Release labels sit well below
    and to the left, because their anchor is the tailrace gage and the station pins
    are drawn at that same point.

    Parameters
    ----------
    row : namedtuple
        One row of ``Data/river_label_anchors.csv`` as produced by ``itertuples``,
        carrying the river name, its role, the anchor gage and its coordinate.
    dark_background : bool, optional
        Invert the label colour for imagery and dark cartography.
    scale : float, optional
        Multiplier on every type size and pixel offset. 1.0 is the on-screen size;
        the report captures use the larger scale computed in :func:`main`.

    Returns
    -------
    folium.Marker
        Marker carrying a DivIcon with the label, and a tooltip naming the anchor gage
        so a reader can see the label is pinned to a surveyed point on that river
        rather than to an estimated river position.
    """
    inflow = row.role == "inflow"
    text = (f"{row.river}<br><span style='font-weight:400;"
            f"font-size:{10 * scale:.0f}px'>"
            f"{'inflow' if inflow else 'release'}</span>")
    align = "left" if inflow else "right"
    color = "#eaf4ff" if dark_background else RIVER_COLOR
    html = (f"<div style='font-family:Helvetica,Arial,sans-serif;"
            f"font-size:{12 * scale:.0f}px;"
            f"font-weight:700;color:{color};white-space:nowrap;"
            f"text-align:{align};{halo(dark_background)};line-height:1.15'>"
            f"{text}</div>")
    # Sit the inflow label above its anchor gage: each of these gages is at a town
    # whose name the basemap prints just below the point, and a lower label collides
    # with it.
    anchor = ((-10 * scale, 32 * scale) if inflow
              else (160 * scale, -22 * scale))
    return folium.Marker(
        location=(row.latitude, row.longitude),
        icon=folium.DivIcon(html=html, icon_size=(150 * scale, 30 * scale),
                            icon_anchor=anchor),
        tooltip=f"{row.river} ({row.role}); label anchored on USGS {row.anchor_gage}")


def build_map(scope: str, sta: pd.DataFrame, riv: pd.DataFrame,
              basemap: str = DEFAULT_BASEMAP, chrome: bool = True,
              scale: float = 1.0) -> folium.Map:
    """Build one folium map.

    scope is 'both', 'Lanier', or 'Allatoona'. basemap names the layer that starts
    active; every entry in BASEMAPS is added and remains switchable in the HTML.

    chrome draws the title and legend panels, which are wanted on screen but not in
    the report, where the LaTeX caption says the same things. scale multiplies the
    label type sizes: a figure reduced to text width needs much larger type in the
    capture than the same map does on a monitor.

    Parameters
    ----------
    scope : str
        ``both``, ``Lanier`` or ``Allatoona``. Selects which stations and river labels
        are drawn and which title text is used.
    sta : pandas.DataFrame
        Contents of ``Data/station_coordinates.csv``.
    riv : pandas.DataFrame
        Contents of ``Data/river_label_anchors.csv``, minus UNLABELLED_RIVERS.
    basemap : str, optional
        Display name of the layer that starts active. Must be a key of BASEMAP_DARK.
    chrome : bool, optional
        Draw the title and legend panels.
    scale : float, optional
        Multiplier on the label type sizes and on the frame padding.

    Returns
    -------
    folium.Map
        Map fitted to the bounding box of the drawn stations and river anchors, padded
        in proportion to the type size so no label is clipped at the frame edge.

    Notes
    -----
    The view is fitted to the data rather than set from a hard-coded centre and zoom,
    so adding a station to the coordinate file reframes the map instead of pushing the
    new point out of frame.
    """
    s = sta if scope == "both" else sta[sta.project == scope]
    r = riv if scope == "both" else riv[riv.project == scope]
    r = r[~r.river.isin(UNLABELLED_RIVERS)]
    dark = BASEMAP_DARK[basemap]

    m = folium.Map(tiles=None, control_scale=True)
    # Leaflet activates whichever visible base layer is added last, so the requested
    # one is held back and added after the rest.
    for name, path, _ in BASEMAPS:
        if name != basemap:
            tile_layer(name, path, show=False).add_to(m)
    active = next(p for n, p, _ in BASEMAPS if n == basemap)
    tile_layer(basemap, active, show=True).add_to(m)

    rivers = folium.FeatureGroup(name="River labels", show=True)
    for row in r.itertuples():
        river_label(row, dark_background=dark, scale=scale).add_to(rivers)
    rivers.add_to(m)

    for project, grp in s.groupby("project"):
        layer = folium.FeatureGroup(name=PROJECT_LABEL[project], show=True)
        for row in grp.itertuples():
            color = PROJECT_COLOR[project]
            folium.CircleMarker(
                location=(row.latitude, row.longitude), radius=11,
                color=color, weight=2, fill=True, fill_color=color,
                fill_opacity=0.25).add_to(layer)
            folium.Marker(
                location=(row.latitude, row.longitude),
                tooltip=f"{row.station_id} - {TYPE_LABEL[row.station_type]}",
                popup=popup_html(row),
                icon=folium.Icon(color="blue" if project == "Lanier" else "orange",
                                 icon=TYPE_ICON[row.station_type],
                                 prefix="glyphicon")).add_to(layer)
            # The forebay and tailrace stations sit about 1.2 to 1.4 km apart, so at
            # reservoir scale their pins nearly coincide. Label each one, offsetting
            # the two in opposite directions so both stay readable in the PNG.
            forebay = row.station_type == "forebay_profile"
            short = row.station_id.replace("21GAEPD_WQX-", "")
            lab = ("Forebay" if forebay else "Tailrace") + f" {short}"
            text_color = "#ffffff" if dark else color
            folium.Marker(
                location=(row.latitude, row.longitude),
                icon=folium.DivIcon(
                    html=f"<div style='font-family:Helvetica,Arial,sans-serif;"
                         f"font-size:{11 * scale:.0f}px;font-weight:700;"
                         f"color:{text_color};"
                         f"white-space:nowrap;{halo(dark)}'>{lab}</div>",
                    icon_size=(170 * scale, 16 * scale),
                    icon_anchor=((-12 * scale, 30 * scale) if forebay
                                 else (-12 * scale, -8 * scale)))).add_to(layer)
        layer.add_to(m)

    lats = list(s.latitude) + list(r.latitude)
    lons = list(s.longitude) + list(r.longitude)
    # Labels are anchored beside their point, so the frame is padded in proportion to
    # the type size, which also keeps the pins clear of the title and legend panels.
    pad = int(85 * max(1.0, scale))
    m.fit_bounds([[min(lats), min(lons)], [max(lats), max(lons)]],
                 padding=(pad, pad))

    head, sub = TITLE[scope]
    title = ("<div style='position:fixed;top:10px;left:50px;z-index:9999;max-width:540px;"
             "background:rgba(255,255,255,0.93);padding:8px 12px;border:1px solid #999;"
             "border-radius:4px;font-family:Helvetica,Arial,sans-serif;font-size:13px'>"
             f"<b>{head}</b><br>{sub}</div>")
    projects = s.project.unique()
    swatches = "".join(
        f"<span style='color:{PROJECT_COLOR[p]}'>&#9679;</span> {PROJECT_LABEL[p]}<br>"
        for p in projects)
    # The legend sits bottom right: on all three maps the stations fall in the left
    # or central part of the frame, so the right corner is the free one.
    legend = ("<div style='position:fixed;bottom:22px;right:12px;z-index:9999;"
              "background:rgba(255,255,255,0.93);padding:8px 12px;border:1px solid #999;"
              "border-radius:4px;font-family:Helvetica,Arial,sans-serif;font-size:12px'>"
              "<b>Stations</b><br>" + swatches +
              "<span style='color:#555'>&#128167; forebay profile "
              "&nbsp;&nbsp;&#128246; tailrace</span><br>"
              f"<span style='color:{RIVER_COLOR};font-weight:700'>River</span>"
              "<span style='color:#666'> labels anchored on USGS gages named for "
              "that river.<br>Tailrace markers are the USGS gages at each dam; the "
              "continuous<br>sondes are co-located with them.</span></div>")
    if chrome:
        m.get_root().html.add_child(folium.Element(title))
        m.get_root().html.add_child(folium.Element(legend))
    folium.LayerControl(collapsed=True).add_to(m)
    return m


def export_png(html_path: str, png_path: str, width: int = 1500, height: int = 1050,
               scale: int = 2, wait: float = 7.0, bare: bool = False,
               control_scale: float = 1.0) -> None:
    """Screenshot a saved folium map with headless Chrome.

    Zoom and layer controls are hidden before capture so the image reads as a
    figure rather than a screenshot of a web application. With bare set, the
    attribution is hidden as well, for report figures whose caption carries the
    basemap credit, and the scale bar is enlarged by control_scale so that it stays
    readable once the figure is reduced to text width.

    Parameters
    ----------
    html_path : str
        Saved folium HTML to load, as a local file.
    png_path : str
        Capture destination. Overwritten if it exists.
    width, height : int, optional
        Browser window size in CSS pixels. The captured image is this multiplied by
        ``scale``.
    scale : int, optional
        Device pixel ratio. 2 gives a capture twice the window size in each direction.
    wait : float, optional
        Seconds to wait after load for the basemap tiles to arrive. This is a fixed
        wait, not a readiness check, so a slow network can produce a capture with
        missing tiles; imagery layers are given a longer wait for that reason. Any
        capture should be looked at before it is placed in the report.
    bare : bool, optional
        Also hide the attribution control, for report figures whose caption carries
        the basemap credit.
    control_scale : float, optional
        CSS scale applied to the scale bar, anchored at its bottom left.
    """
    from selenium import webdriver
    from selenium.webdriver.chrome.options import Options

    opts = Options()
    opts.add_argument("--headless=new")
    opts.add_argument(f"--window-size={width},{height}")
    opts.add_argument(f"--force-device-scale-factor={scale}")
    opts.add_argument("--hide-scrollbars")
    driver = webdriver.Chrome(options=opts)
    try:
        driver.get("file://" + os.path.abspath(html_path))
        time.sleep(wait)  # let the basemap tiles finish loading
        hide = ".leaflet-control-zoom,.leaflet-control-layers"
        if bare:
            hide += ",.leaflet-control-attribution"
        driver.execute_script(
            "document.querySelectorAll(arguments[0])"
            ".forEach(function(e){e.style.display='none';});", hide)
        if control_scale != 1.0:
            driver.execute_script(
                "document.querySelectorAll('.leaflet-control-scale')"
                ".forEach(function(e){e.style.transformOrigin='left bottom';"
                "e.style.transform='scale(' + arguments[0] + ')';});", control_scale)
        time.sleep(0.6)
        driver.get_screenshot_as_file(png_path)
    finally:
        driver.quit()


def downsample(src: str, dst: str, width: int) -> None:
    """Write a width-limited copy of a capture for inclusion in the report.

    Parameters
    ----------
    src : str
        Capture to read.
    dst : str
        Destination. May be the same path as ``src``, which is how :func:`main` uses
        it: the full capture is replaced in place by its reduced copy.
    width : int
        Maximum width in pixels. A narrower image is left alone rather than upscaled.

    Notes
    -----
    Lanczos resampling is used because these captures carry fine label type, which
    bilinear reduction softens noticeably at this ratio. Aspect ratio is preserved.
    """
    from PIL import Image
    im = Image.open(src)
    if im.width > width:
        im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
    im.save(dst, optimize=True)


def main() -> None:
    """Build all three maps, capture them, and print the station separations.

    Writes, per scope, an interactive HTML with the title and legend panels, a second
    interactive copy without them for taking screen captures of individual features,
    a default-basemap PNG and a satellite PNG in ``analysis/maps/``, and a separate
    report capture written straight into ``report/latex/src/images/``.

    The report capture is a separate build, not a reduction of the screen one: at text
    width the on-screen type would fall to about 3 pt, so the labels are enlarged by
    ``label_scale`` before capture and the panels are dropped, because the LaTeX
    caption repeats what they say.
    """
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(IMAGE_DIR, exist_ok=True)
    sta = pd.read_csv(COORDS)
    riv = pd.read_csv(RIVERS)

    # Window shape is matched to each map's content so the fitted view wastes little
    # space: the two-project map spans east to west, each single-project map is
    # roughly as tall as it is wide.
    for scope, stem, size in (("both", "station_map", (1700, 900)),
                              ("Lanier", "station_map_lanier", (1250, 1150)),
                              ("Allatoona", "station_map_allatoona", (1400, 1050))):
        # The interactive map is saved once, opening on the default basemap with
        # every other basemap available from the layer control.
        html = os.path.join(OUT_DIR, stem + ".html")
        build_map(scope, sta, riv, basemap=DEFAULT_BASEMAP).save(html)
        print(f"wrote {os.path.relpath(html, REPO)}")

        # A second interactive copy without the title and legend panels, for taking
        # screen captures of individual features: the basemap selector and the zoom
        # control are kept so the view can be moved and switched to imagery, but
        # nothing overlays the map itself.
        plain = os.path.join(OUT_DIR, stem + "_plain.html")
        build_map(scope, sta, riv, basemap=DEFAULT_BASEMAP, chrome=False).save(plain)
        print(f"wrote {os.path.relpath(plain, REPO)}")
        png = os.path.join(OUT_DIR, stem + ".png")
        export_png(html, png, width=size[0], height=size[1])
        print(f"wrote {os.path.relpath(png, REPO)}")
        rsize = REPORT_SIZE[scope]
        # Capture pixels per text-block point, times the point size wanted on paper,
        # divided by the base label size in CSS pixels.
        label_scale = rsize[0] / TEXT_WIDTH_PT * REPORT_LABEL_PT / BASE_LABEL_PX
        report_png = os.path.join(IMAGE_DIR, stem + ".png")
        fd, tmp_r = tempfile.mkstemp(suffix=".html")
        os.close(fd)
        try:
            build_map(scope, sta, riv, basemap=DEFAULT_BASEMAP, chrome=False,
                      scale=label_scale).save(tmp_r)
            export_png(tmp_r, report_png, width=rsize[0], height=rsize[1],
                       bare=True, control_scale=label_scale)
            downsample(report_png, report_png, REPORT_WIDTH)
            print(f"wrote {os.path.relpath(report_png, REPO)}")
        finally:
            os.unlink(tmp_r)

        # The satellite still is captured from a throwaway copy of the same map with
        # the imagery active, so no second interactive file has to be maintained.
        sat_png = os.path.join(OUT_DIR, stem + "_satellite.png")
        fd, tmp = tempfile.mkstemp(suffix=".html")
        os.close(fd)
        try:
            build_map(scope, sta, riv, basemap=SATELLITE_BASEMAP).save(tmp)
            # Imagery tiles are heavier than the OSM raster, so allow more load time.
            export_png(tmp, sat_png, width=size[0], height=size[1], wait=11.0)
            print(f"wrote {os.path.relpath(sat_png, REPO)}")
        finally:
            os.unlink(tmp)

    print("\nStation separations (great-circle, from the archived coordinates):")
    for project, grp in sta.groupby("project"):
        f = grp[grp.station_type == "forebay_profile"].iloc[0]
        t = grp[grp.station_type == "tailrace"].iloc[0]
        d = haversine_km((f.latitude, f.longitude), (t.latitude, t.longitude))
        print(f"  {project:10s} forebay {f.station_id} to tailrace {t.station_id}: "
              f"{d:.2f} km")


if __name__ == "__main__":
    main()
