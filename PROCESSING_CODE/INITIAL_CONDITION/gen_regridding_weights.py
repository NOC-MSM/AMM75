# Script to produce the regridding weights
import xarray as xr
import numpy as np
import xesmf as xe

PATH_DEST_DOMCFG = "domain_cfg.nc"
#PATH_SRC_T_DATA = "/bodc/SOC220065/GO8p7_JRA55_eORCA12/monthly/T/nemo_g8p7ho_1m_200501-200501_grid-T.nc"
#PATH_SRC_TMASK = "mesh_mask_eORCA025-GO6_subset.nc"
PATH_SRC_DOMCFG_SUBSET = "mesh_mask_eORCA025-GO6_subset.nc"
PATH_WEIGHTS_DIR = "WEIGHTS/"


####################################################
# step 1 get the destination grid sorted
####################################################

# Set path to destination domain_cfg
dom_ds = xr.open_dataset( PATH_DEST_DOMCFG ).squeeze()

# create a file that defines the centre and edge points of the T cells.
ds_out = xr.Dataset()
# t-cells are the centres
ds_out["lon"] = dom_ds.glamt[:,:]
ds_out["lat"] = dom_ds.gphit[:,:]

# f-grid are the edges, for the edges we need to expand the 
# f-grid at the lower and left edges
del_lon = np.diff( dom_ds.glamf[0,:2] )
new_lon = (dom_ds.glamf[:,0].values - del_lon[np.newaxis]).squeeze()
glamf_tmp = np.insert( dom_ds.glamf.values, 0, new_lon, axis=1 )
glamf     = np.insert( glamf_tmp, 0, glamf_tmp[0,:], axis=0 )
ds_out["lon_b"] = xr.DataArray( glamf, dims=["y_b", "x_b"] )

del_lat = np.diff( dom_ds.gphif[:2,0] )
new_lat = ( dom_ds.gphif[0,:].values - del_lat[np.newaxis] ).squeeze()
gphif_tmp = np.insert( dom_ds.gphif.values, 0, new_lat, axis=0 )
gphif     = np.insert( gphif_tmp, 0, gphif_tmp[:,0], axis=1 )
ds_out["lat_b"] = xr.DataArray( gphif, dims=["y_b", "x_b"] )

####################################################
# step 2 get source grid and data
####################################################

# load in the global go8 temp and salinity and mask (pre made with)
#mask = xr.where(~np.isnan(ds_src.vosaline), 1, 0)
#ds_src = xr.open_dataset( PATH_SRC_T_DATA )
#ds_src = ds_src[ ["votemper", "vosaline"] ]
#ds_src["mask"] = xr.open_dataset( PATH_SRC_TMASK ).mask

# Load in src global mask t-grid 
#mask = xr.where(~np.isnan(ds_src.vosaline), 1, 0)
ds_src = xr.open_dataset( PATH_SRC_DOMCFG_SUBSET ).squeeze()
ds_src["mask"] = ds_src.tmask

# We need the src domain_cfg, for efficiency we have made a subset of the global domain
#go8_ds = xr.open_dataset( PATH_SRC_DOMCFG_SUBSET ).squeeze()
# the x and y indicies for the subset slice are saved, e.g.
#print( go8_ds.y_slice )
#print( go8_ds.x_slice )

# subset the temp and salinity to a smaller domain 
t_s_src_ds = ds_src.isel( y = slice(1, None), x = slice(1, None) ).reset_coords()
# add the corner points (f-grid) to allow for conservative regridding
t_s_src_ds["lat_b"] = ds_src.gphif.rename({"y":"y_b", "x":"x_b"})
t_s_src_ds["lon_b"] = ds_src.glamf.rename({"y":"y_b", "x":"x_b"})
# add the centre points (t-grid) to allow for conservative regridding

t_s_src_ds["lat"] = t_s_src_ds.gphit
t_s_src_ds["lon"] = t_s_src_ds.glamt

t_s_src_ds = t_s_src_ds.load()
print(t_s_src_ds)

####################################################
# Step 3 create the regridders and save to disk
####################################################

# create the regirdder weights for each depth and save to disk
# Note conservative_normed stops NANS bleeding into the domain 
# so that we don't unnecessarily have NANS on the finer grid near
# the coast. This has to be done level by level because of the mask
for i in np.arange(0,75):
    t_s_src_ds_i = t_s_src_ds.isel(z=i)
    regrid_con_normed_i = xe.Regridder(t_s_src_ds_i, ds_out, method="conservative_normed")
    regrid_con_normed_i.attrs=dict( description="GLOSEA6 -> AMM75 weights, conservative_normed, level " + str(i) )
    regrid_con_normed_i.to_netcdf( PATH_WEIGHTS_DIR + "con_normed_level_" + str(i) + ".nc" )

# Compute the weights files for bilinear interpolation with extrapolation
# This is required to flood the land points so that there are no NANS on the 
# finer grid around the coast
for i in np.arange(0,75):
    t_s_src_ds_i = t_s_src_ds.isel(z=i)
    regrid_bilin_extrap_i = xe.Regridder(
        t_s_src_ds_i, ds_out, 
        method="bilinear", extrap_method="nearest_s2d", ignore_degenerate=True
    )
    regrid_bilin_extrap_i.attrs=dict( 
        description="GLOSEA6 -> AMM75 weights, bilinear with nearest neighbour extrapolation, level " + str(i) )
    regrid_bilin_extrap_i.to_netcdf( PATH_WEIGHTS_DIR +"bilinear_extrap_level_" + str(i) + ".nc" )
