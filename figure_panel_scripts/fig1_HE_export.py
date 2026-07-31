#!/usr/bin/env python
# Orient each platform's H&E OME-TIFF to its spatial map and save an RGB PNG.
import tifffile, numpy as np, os
from PIL import Image
from scipy import ndimage
Image.MAX_IMAGE_PIXELS = None

OUT = "/path/to/platform_comparison/results/manuscript_cache"
os.makedirs(OUT, exist_ok=True)
MAXDIM = 2600

# name -> (file, level, (rot90ccw_k, flipLR), component)
# component: None = whole tissue; "small"/"big" = that connected fragment only
spec = {
    "atera": ("/path/to/h_e/atera/WTA_Preview_FFPE_Breast_Cancer_he_image.ome.tif", 4, (3, True), "small"),
    "xenium_280": ("/path/to/h_e/xenium_280/Human_Breast_Biomarkers_S1_Bot_he_image.ome.tif", 4, (1, True), None),
    "xenium_5k": ("/path/to/h_e/xenium_5k/Xenium_Prime_Breast_Cancer_FFPE_he_image.ome.tif", 4, (2, False), None),
}

def load_level(f, lvl):
    with tifffile.TiffFile(f) as tif:
        s = tif.series[0]
        levels = s.levels if hasattr(s, "levels") else [s]
        arr = levels[min(lvl, len(levels) - 1)].asarray()
    if arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[0] < arr.shape[-1]:
        arr = np.moveaxis(arr, 0, -1)              # (S,Y,X) -> (Y,X,S)
    return arr[..., :3].astype(np.uint8)

for name, (f, lvl, (k, flip), comp) in spec.items():
    arr = load_level(f, lvl)
    if comp in ("small", "big"):
        # Keep just ONE tissue fragment as a NATURAL rectangular crop (no masking
        # -> no ragged edges). Locate the fragment's bounding box via a saturation
        # mask (stained tissue is saturated; the faint gray slide-edge frame is
        # not), then crop the ORIGINAL pixels to that box.
        af = arr.astype(float)
        sat = (af.max(2) - af.min(2)) / (af.max(2) + 1e-6)
        tissue0 = (sat > 0.10) & (af.max(2) < 245)
        opened = ndimage.binary_opening(tissue0, iterations=2)
        lab, n = ndimage.label(opened)
        sizes = ndimage.sum(np.ones_like(lab), lab, range(1, n + 1))
        order = np.argsort(sizes)[::-1] + 1        # component ids largest->smallest
        cid = order[0] if comp == "big" else order[1]
        sel = lab == cid
        rows = np.where(sel.any(1))[0]; cols = np.where(sel.any(0))[0]
        pad = 15
        r0 = max(rows.min() - pad, 0); r1 = min(rows.max() + pad, arr.shape[0] - 1)
        c0 = max(cols.min() - pad, 0); c1 = min(cols.max() + pad, arr.shape[1] - 1)
        arr = arr[r0:r1 + 1, c0:c1 + 1]            # natural crop, no whitening
    arr = np.rot90(arr, k)
    if flip:
        arr = arr[:, ::-1]
    if comp is None:
        # whole-tissue crop to bounding box (drop surrounding white)
        tissue = arr.mean(2) < 225
        rows = np.where(tissue.any(1))[0]; cols = np.where(tissue.any(0))[0]
        if rows.size and cols.size:
            arr = arr[rows.min():rows.max() + 1, cols.min():cols.max() + 1]
    im = Image.fromarray(arr)
    w, h = im.size
    if max(w, h) > MAXDIM:
        sc = MAXDIM / max(w, h); im = im.resize((int(w * sc), int(h * sc)))
    im.save(os.path.join(OUT, f"he_oriented_{name}.png"))
    print(name, "level", lvl, "component", comp, "-> oriented", im.size)
print("done")
