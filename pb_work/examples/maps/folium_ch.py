"""=================================================================================================


Michael Gribskov     06 April 2021
================================================================================================="""
import folium
import pandas as pd

# --------------------------------------------------------------------------------------------------
#
# --------------------------------------------------------------------------------------------------
url = (
    "https://raw.githubusercontent.com/python-visualization/folium/master/examples/data"
)
state_geo = f"{url}/us-states.json"
state_unemployment = f"{url}/US_Unemployment_Oct2012.csv"
state_data = pd.read_csv(state_unemployment)

# initialize the map and store it in a m object
m = folium.Map(location=[40, -95], zoom_start=4)

folium.Choropleth(
    geo_data=state_geo,
    name="choropleth",
    data=state_data,
    columns=["State", "Unemployment"],
    key_on="feature.id",
    fill_color="BuPu",
    fill_opacity=0.25,
    line_opacity=.1,
    legend_name="Unemployment Rate (%)",
).add_to(m)

folium.LayerControl().add_to(m)

m.save('folium_chloropleth.html')