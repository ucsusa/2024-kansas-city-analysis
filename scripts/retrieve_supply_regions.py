import geopandas as gpd


if __name__ == "__main__":
    gdf = gpd.read_file(snakemake.input.community)
        
    centroids = (gdf.to_crs(epsg=snakemake.config['projected_crs'])\
                 .centroid.to_crs(epsg=snakemake.config['geographic_crs']))
    gdf['x'] = centroids.x
    gdf['y'] = centroids.y

    gdf['name'] = [snakemake.config['community_name']]

    gdf.to_file("data/spatial_data/supply_regions.shp")
    
    