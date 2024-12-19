import rasterio
from rasterio.mask import mask
import geopandas as gpd
import os
import shutil


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


clip("wimmera","all.geojson","out")
