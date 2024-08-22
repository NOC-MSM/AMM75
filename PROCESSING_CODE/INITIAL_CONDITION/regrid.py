# Script to do conservative regridding of Temperature and Salinity
# From z-level source to s-level destination
import xarray as xr
import numpy as np
import xesmf as xe
from xgcm import Grid

####################################################
# step 0 edit paths
# Destination domain_cfg, Source t-grid data file,
# Source t-grid mask, Source domain_cfg (subset),
# Weights directory, Output filename, 
# Description of source data.
####################################################
PATH_DEST_DOMCFG = "domain_cfg.nc"
PATH_SRC_T_DATA = "/gws/nopw/j04/jmmp/MASS/GloSea6/Daily/glosea6_grid_T_19930101.nc"
PATH_SRC_TMASK = "mesh_mask_eORCA025-GO6_subset.nc"
PATH_SRC_DOMCFG_SUBSET = "mesh_mask_eORCA025-GO6_subset.nc"
PATH_WEIGHTS_DIR = "WEIGHTS/"
FN_OUT = "IC_TandS_GLOASEA6_y1993m01d01.nc"
source_description = "Glosea6 - /gws/nopw/j04/jmmp/MASS/GloSea6/Daily/glosea6_grid_T_19930101.nc"


def get_cell_bounds(e3):
# get the current top and bottom bounds of the cell
    dims = e3.shape
    bounds = np.zeros( ( dims[0]+1, dims[1], dims[2]) ) #z,y,x
    bounds[1:,:,:] = e3.cumsum(dim='z').data
    centre = np.diff(bounds, axis=0) / 2 + bounds[:-1,:,:]
    return bounds, centre
    
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
ds_src = xr.open_dataset( PATH_SRC_T_DATA ).squeeze()
ds_src = ds_src[ ["votemper", "vosaline"] ]
ds_src["mask"] = xr.open_dataset( PATH_SRC_TMASK ).squeeze().tmask.rename({"z":"deptht"})

# We need the src domain_cfg, for efficiency we have made a subset of the global domain
src_dom_ds = xr.open_dataset( PATH_SRC_DOMCFG_SUBSET ).squeeze()
# the x and y indicies for the subset slice are saved, e.g.
#print( go8_ds.y_slice )
#print( go8_ds.x_slice )

# subset the temp and salinity to a smaller domain (needs to be one smaller to 
# get f-grid points outside t-grid
t_s_src_ds = ds_src.isel( y = slice(1, None), x = slice(1, None) ).reset_coords()
# add the corner points (f-grid) to allow for conservative regridding
t_s_src_ds["lat_b"] = src_dom_ds.gphif.rename({"y":"y_b", "x":"x_b"})
t_s_src_ds["lon_b"] = src_dom_ds.glamf.rename({"y":"y_b", "x":"x_b"})
t_s_src_ds = t_s_src_ds.rename({"nav_lat":"lat", "nav_lon":"lon"})

t_s_src_ds = t_s_src_ds.load()
    
####################################################
# Step 3 perform the regrid on temp and salinity 
####################################################

# Load in the weights for the conservative regridding and regrid T&S
# then store each level into a DataArray for temp and salinity
tem_conserve_ls = []
sal_conserve_ls = []
for i in np.arange(0,75):
    fn_i = PATH_WEIGHTS_DIR + "con_normed_level_" + str(i) + ".nc"
    t_s_src_ds_i = t_s_src_ds.isel(deptht=i)    
    regridder_con = xe.Regridder(t_s_src_ds_i, ds_out, method="conservative_normed", weights=fn_i, unmapped_to_nan=True)
    tem_conserve_ls.append( regridder_con( t_s_src_ds_i.votemper ) )
    sal_conserve_ls.append( regridder_con( t_s_src_ds_i.vosaline ) )

tem_con = xr.concat( tem_conserve_ls, dim=t_s_src_ds.deptht )
sal_con = xr.concat( sal_conserve_ls, dim=t_s_src_ds.deptht )

# Load in the weights for the bilinear extrapolating regridding and regrid T&S
# then store each level into a DataArray for temp and salinity
tem_bil_extrap_ls = []
sal_bil_extrap_ls = []
for i in np.arange(0,75):
    fn_i = PATH_WEIGHTS_DIR + "bilinear_extrap_level_" + str(i) + ".nc"
    t_s_src_ds_i = t_s_src_ds.isel(deptht=i)    
    regridder_bil_extrap = xe.Regridder(
        t_s_src_ds_i, ds_out, 
        method="bilinear", extrap_method="nearest_s2d", ignore_degenerate=True,
        weights = fn_i
    )
    tem_bil_extrap_ls.append( regridder_bil_extrap( t_s_src_ds_i.votemper ) )
    sal_bil_extrap_ls.append( regridder_bil_extrap( t_s_src_ds_i.vosaline ) )

tem_bil = xr.concat( tem_bil_extrap_ls, dim=t_s_src_ds.deptht )
sal_bil = xr.concat( sal_bil_extrap_ls, dim=t_s_src_ds.deptht )

# Merge the DataArrays so the bilinear is only used at point 
# that are masked in the conservative fields. We do this so there are
# no masked points in the destination grid that shouldn't be
tem_src = xr.where( np.isnan(tem_con), tem_bil, tem_con ) 
sal_src = xr.where( np.isnan(sal_con), sal_bil, sal_con ) 

####################################################
# Step 4 do the vertical interpolation
####################################################

bounds, centre = get_cell_bounds( src_dom_ds.e3t_0 )
bounds_dest, centre_dest = get_cell_bounds( dom_ds.e3t_0 )

ds_src_pre_zinterp = xr.Dataset()
# make temp and sal extensive quantities
e3t_0_1d = src_dom_ds.e3t_0[:,0,0].rename({"z":"deptht"})
ds_src_pre_zinterp["tem_src"] = tem_src * e3t_0_1d
ds_src_pre_zinterp["sal_src"] = sal_src * e3t_0_1d
ds_src_pre_zinterp = ds_src_pre_zinterp.assign_coords(
    {'zc':bounds[:,0,0]})
#tem_src = tem_src.expand_dims({'zc':76}).assign_coords(zc = ('zc', bounds[:,0,0].squeeze()))
grid = Grid(ds_src_pre_zinterp, coords={'Z':{'center':'deptht', 'outer':'zc'}}, periodic=False)
# for mapping to / from coordinates that are indepenent at each grid point
# there is no way to do this without looping
tem_dest_data = np.empty_like( dom_ds.e3t_0 )
sal_dest_data = np.empty_like( dom_ds.e3t_0 )
for x_i in ds_src_pre_zinterp.x:    
    for y_i in ds_src_pre_zinterp.y:
        tem_dest_data[:,y_i,x_i] = grid.transform(
            ds_src_pre_zinterp.tem_src[:,y_i,x_i], 'Z', bounds_dest[:,y_i,x_i], method='conservative' )
        sal_dest_data[:,y_i,x_i] = grid.transform(
            ds_src_pre_zinterp.sal_src[:,y_i,x_i], 'Z', bounds_dest[:,y_i,x_i], method='conservative' )

####################################################
# Step 5 output to disk
####################################################
        
# Create destination dataset for output
dest_ds = xr.Dataset(attrs=dict(description="Regridded from: " + source_description))
dest_ds["votemper"] = xr.DataArray( tem_dest_data / dom_ds.e3t_0, dims=["z","y","x"], 
              coords = dict( lat=(["y","x"], tem_src.lat.data ), lon=(["y","x"], tem_src.lon.data ),
                            deptht=(["z","y","x"], centre_dest)),
              attrs=dict( description="conservative temperature", units="degC"), name="votemper" )

dest_ds["vosaline"] = xr.DataArray( sal_dest_data / dom_ds.e3t_0, dims=["z","y","x"], 
              coords = dict( lat=(["y","x"], sal_src.lat.data ), lon=(["y","x"], sal_src.lon.data ),
                            deptht=(["z","y","x"], centre_dest)),
              attrs=dict( description="absolute salinity", units="PSU"), name="vosaline" )

dest_ds.to_netcdf( FN_OUT )
