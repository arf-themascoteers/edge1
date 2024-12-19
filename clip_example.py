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

def gen_ndvi(dest_folder, sis_folder):
    sis_subfolder = os.path.join(sis_folder, "NDVI")
    os.makedirs(sis_subfolder, exist_ok=True)
    for site in os.listdir(dest_folder):
        site_folder = os.path.join(dest_folder, site,"R10m")
        b4_file = get_band_tile(site_folder,"B04")
        b8_file = get_band_tile(site_folder,"B08")
        b4_path = os.path.join(site_folder, b4_file)
        b8_path = os.path.join(site_folder, b8_file)
        with rasterio.open(b4_path) as b4_src, rasterio.open(b8_path) as b8_src:
            b4 = b4_src.read(1).astype(np.float32)
            b8 = b8_src.read(1).astype(np.float32)
            np.seterr(divide='ignore', invalid='ignore')
            ndvi = (b8 - b4) / (b8 + b4)
            ndvi[np.isnan(ndvi)] = 0
        plt.imshow(ndvi, cmap='RdYlGn')
        plt.colorbar(label='NDVI')
        plt.axis('off')
        output_path = os.path.join(sis_subfolder, f"{site}.tif")
        plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
        plt.close()

def gen_evi(dest_folder, sis_folder):
    sis_subfolder = os.path.join(sis_folder, "EVI")
    os.makedirs(sis_subfolder, exist_ok=True)
    for site in os.listdir(dest_folder):
        site_folder = os.path.join(dest_folder, site,"R10m")
        b2_file = get_band_tile(site_folder,"B02")
        b4_file = get_band_tile(site_folder,"B04")
        b8_file = get_band_tile(site_folder,"B08")
        b2_path = os.path.join(site_folder, b2_file)
        b4_path = os.path.join(site_folder, b4_file)
        b8_path = os.path.join(site_folder, b8_file)
        with rasterio.open(b2_path) as b2_src, rasterio.open(b4_path) as b4_src, rasterio.open(b8_path) as b8_src:
            b2 = b2_src.read(1).astype(np.float32)
            b4 = b4_src.read(1).astype(np.float32)
            b8 = b8_src.read(1).astype(np.float32)
            np.seterr(divide='ignore', invalid='ignore')
            evi = 2.5 * (b8 - b4) / (b8 + 6 * b4 - 7.5 * b2 + 1)
            evi[np.isnan(evi)] = 0
        plt.imshow(evi, cmap='RdYlGn')
        plt.colorbar(label='EVI')
        plt.axis('off')
        output_path = os.path.join(sis_subfolder, f"{site}.tif")
        plt.savefig(output_path, bbox_inches='tight', pad_inches=0)
        plt.close()



def generate_sis(src_folder,aoi,dest_folder,sis_folder,skip_clip=False):
    if not skip_clip:
        clip(src_folder,aoi,dest_folder)
    gen_ndvi(dest_folder,sis_folder)
    gen_evi(dest_folder,sis_folder)


generate_sis("wimmera","all.geojson","out","sis1",True)