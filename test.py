import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

data = gpd.read_file('site_1.geojson')
data.plot()
plt.show()
