import rasterio
from rasterio.mask import mask
import geopandas as gpd
import os
import shutil
import numpy as np
from matplotlib import pyplot as plt


def clip_tile(src_folder, res, tile, aoi, dest_folder):
    geojson = gpd.read_file(aoi)
    tile_loc = os.path.join(src_folder, res, tile)
    with rasterio.open(tile_loc) as src:
        raster_crs = src.crs
        if geojson.crs != raster_crs:
            geojson = geojson.to_crs(raster_crs)

        for idx, feature in enumerate(geojson.__geo_interface__["features"]):
            polygon = [feature["geometry"]]
            try:
                out_image, out_transform = mask(src, polygon, crop=True)
                out_meta = src.meta

                out_meta.update({
                    "driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform
                })
                clipped_folder = str(os.path.join(dest_folder, str(idx+1),res))
                os.makedirs(clipped_folder, exist_ok=True)
                clipped_loc = os.path.join(clipped_folder, f"{tile}.tif")
                with rasterio.open(clipped_loc, 'w', **out_meta) as dest:
                    dest.write(out_image)
            except ValueError:
                print(f"Polygon {idx + 1} does not overlap raster {tile}")


def clip(src_folder,aoi,dest_folder):
    if os.path.exists(dest_folder):
        shutil.rmtree(dest_folder)
    os.mkdir(dest_folder)
    for res in os.listdir(src_folder):
        res_folder = os.path.join(src_folder, res)
        for tile in os.listdir(res_folder):
            clip_tile(src_folder, res, tile, aoi, dest_folder)

def get_band_tile(src_folder,band_number):
    for tile in os.listdir(src_folder):
        tokens = tile.split("_")
        if tokens[2] == band_number:
            return tile

def gen_spectral_index(dest_folder, sis_folder, index_name, band_files, index_function, cmap, resolution):
    sis_subfolder = os.path.join(sis_folder, index_name)
    os.makedirs(sis_subfolder, exist_ok=True)
    for site in os.listdir(dest_folder):
        site_folder = os.path.join(dest_folder, site, resolution)
        band_paths = {band: os.path.join(site_folder, get_band_tile(site_folder, band)) for band in band_files}
        with rasterio.open(band_paths[band_files[0]]) as src:
            width, height = src.width, src.height
        bands = {band: rasterio.open(band_paths[band]).read(1).astype(np.float32) for band in band_files}
        np.seterr(divide='ignore', invalid='ignore')
        index = index_function(*[bands[band] for band in band_files])
        index[np.isnan(index)] = 0
        plt.imshow(index, cmap=cmap)
        plt.colorbar(label=index_name)
        plt.axis('off')
        output_path = os.path.join(sis_subfolder, f"{site}.png")
        plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
        plt.close()

def gen_ndvi(dest_folder, sis_folder):
    gen_spectral_index(dest_folder, sis_folder, "NDVI", ["B04", "B08"], lambda b4, b8: (b8 - b4) / (b8 + b4), "RdYlGn", "R10m")

def gen_evi(dest_folder, sis_folder):
    gen_spectral_index(dest_folder, sis_folder, "EVI", ["B02", "B04", "B08"], lambda b2, b4, b8: 2.5 * (b8 - b4) / (b8 + 6 * b4 - 7.5 * b2 + 1), "RdYlGn", "R10m")

def gen_savi(dest_folder, sis_folder, l=0.5):
    gen_spectral_index(dest_folder, sis_folder, "SAVI", ["B04", "B08"], lambda b4, b8: (b8 - b4) / (b8 + b4 + l) * (1 + l), "RdYlGn", "R10m")

def gen_ndwi(dest_folder, sis_folder):
    gen_spectral_index(dest_folder, sis_folder, "NDWI", ["B8A", "B11"], lambda b8, b11: (b8 - b11) / (b8 + b11), "Blues", "R20m")

def gen_pri(dest_folder, sis_folder):
    gen_spectral_index(dest_folder, sis_folder, "PRI", ["B03", "B04"], lambda b3, b4: (b3 - b4) / (b3 + b4), "RdYlBu", "R10m")

def gen_msavi(dest_folder, sis_folder):
    gen_spectral_index(dest_folder, sis_folder, "MSAVI", ["B04", "B08"], lambda b4, b8: (2 * b8 + 1 - ((2 * b8 + 1) ** 2 - 8 * (b8 - b4)) ** 0.5) / 2, "RdYlGn", "R10m")

def gen_ndbi(dest_folder, sis_folder):
    gen_spectral_index(dest_folder, sis_folder, "NDBI", ["B8A", "B11"], lambda b8, b11: (b11 - b8) / (b11 + b8), "Purples", "R20m")

def gen_ndii(dest_folder, sis_folder):
    gen_spectral_index(dest_folder, sis_folder, "NDII", ["B8A", "B11"], lambda b8, b11: (b8 - b11) / (b8 + b11), "YlGn", "R20m")

def gen_chlorophyll_index(dest_folder, sis_folder):
    gen_spectral_index(dest_folder, sis_folder, "Chlorophyll_Index", ["B05", "B8A"], lambda b5, b8a: (b8a / b5) - 1, "Greens", "R20m")



def generate_sis(src_folder,aoi,dest_folder,sis_folder,skip_clip=False):
    if not skip_clip:
        clip(src_folder,aoi,dest_folder)
    gen_ndvi(dest_folder, sis_folder)
    gen_evi(dest_folder, sis_folder)
    gen_savi(dest_folder, sis_folder)
    gen_ndwi(dest_folder, sis_folder)
    gen_pri(dest_folder, sis_folder)
    gen_msavi(dest_folder, sis_folder)
    gen_ndbi(dest_folder, sis_folder)
    gen_ndii(dest_folder, sis_folder)
    gen_chlorophyll_index(dest_folder, sis_folder)


generate_sis("wimmera","all.geojson","out","sis1",True)