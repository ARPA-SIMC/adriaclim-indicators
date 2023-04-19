import xarray as xr
import string
import numpy as np
import numpy
import numpy as np
import os
from xarrayMannKendall import Mann_Kendall_test
import datetime as datetime

############################################################################
'''
#COMPUTING TREND
#Mann Kendall Test

For more information on the Mann-Kendall method please refer to:

Mann, H. B. (1945). Non-parametric tests against trend, Econometrica, 13, 163-171.

Kendall, M. G. (1975). Rank Correlation Methods, 4th edition, Charles Griffin, London.

Yue, S. and Wang, C. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. Water Resources Management, 18(3), 201–218. doi:10.1023/b:warm.0000043140.61082.60

and

Hussain, M. and Mahmud, I. (2019). pyMannKendall: a python package for non parametric Mann Kendall family of trend tests. Journal of Open Source Software, 4(39), 1556. doi:10.21105/joss.01556

Josué Martínez Moreno, & Navid C. Constantinou. (2021, January 23). josuemtzmo/xarrayMannKendall: Mann Kendall significance test implemented in xarray. (Version v.1.0.0). Zenodo. http://doi.org/10.5281/zenodo.4458777
#https://github.com/josuemtzmo/xarrayMannKendall

Init signature:
Mann_Kendall_test(
    DataArray,
    dim='time',
    alpha=0.01,
    MK_modified=False,
    method='linregress',
    coords_name=None,
)
Docstring:     
Compute linear trends and significance using Mann Kendall test.

Parameters
----------
DataArray : xarray.DataArray
    Dataset to analyse.
dim : str
    Coordiante name in which the linear trend will apply ('time').
alpha: float
    Significance level (default = 0.01)
MK_modified: Boolean
    Modified Mann-Kendall using Yue and Wang (2004) method.
    DOI: https://doi.org/10.1023/B:WARM.0000043140.61082.60
method: str
    Method for linear regression: linregress (default) and theilslopes
coords_name: dict
    Coordinates name dict renames coordinates to 'lon','lat'. 
    Example:   
        coords_name={'xu_ocean':'lon','yu_ocean':'lat','t':time}

Same approach used in:         
Josué Martínez-Moreno;Andrew McC. Hogg;Matthew H. England;Navid C. Constantinou;Andrew E. Kiss;Adele K. Morrison; (2021). Global changes in oceanic mesoscale currents over the satellite altimetry record . Nature Climate Change, (), –. doi:10.1038/s41558-021-01006-9     

Trends, significance and uncertainties: Linear trends are calculated using a linear
least-squares regression for spatially integrated time series. 
All the observed trends are assessed using a Theil–Sen estimator, while the statistical
significance uses a modified Mann–Kendall test. This statistical test takes into
account autocorrelations within the time series. Finally, the reported uncertainties
correspond to the standard error using the effective sample
size from the Mann–Kendall test, that is the standard deviation of the time series
divided by the square root of the effective sample size.

'''

############################################################################
# IMPORT DATA#
############################################################################

ind_path='input_path' 
out_path='output_path'

for indicator in os.listdir(ind_path):
	M = os.path.join(ind_path,indicator)
	dataset = xr.open_dataset(M,chunks={'time': 50, 'x': -1, 'y': -1})
	print(M)
	dataset = dataset.T2
	#data=data.sortby('time').sel(time=slice('1993','2020'))	
	trends=Mann_Kendall_test(dataset,'time',MK_modified=False,method="linregress",alpha=0.05,coords_name = {'time':'Times','x':'west_east','y':'south_north'})

	ind_trend = trends.compute()

	ind_trend.attrs['title'] = 'Indicator Trend: the name of the indicator is specified in the file nomenclature'
	ind_trend.attrs['Description'] = 'Trends, significance and uncertainties: Linear trends are calculated using a linear least-squares regression for spatially integrated time series. All the observed trends are assessed using a Theil–Sen estimator, while the statistical significance uses a modified Mann–Kendall test.'
	ind_trend.attrs['Author'] ='Giusy Fedele'
	ind_trend.attrs['Contact'] = 'giusy.fedele@cmcc.it'
	ind_trend.attrs['Created date'] = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
	ind_trend.attrs['adriaclim_model'] = 'WRF'
	ind_trend.attrs['adriaclim_scale'] = 'Adriatic'
	ind_trend.attrs['adriaclim_type'] = 'trend'
	ind_trend.attrs['adriaclim_timeperiod'] = 'yearly' # It is based on yearly data
	ind_trend.attrs['institution'] = 'CMCC'

	######################################################

	ind_trend['trend'].attrs['units'] = "[unit]/yr"
	ind_trend['trend'].attrs['name'] = 'trend'
	ind_trend['trend'].attrs['long_name'] = "Indicator trend over the entire timeseries, Significance level 0.05"
	ind_trend['trend'].attrs['reference'] = "Josué Martínez Moreno, & Navid C. Constantinou. (2021, January 23). josuemtzmo/xarrayMannKendall: Mann Kendall significance test implemented in xarray. (Version v.1.0.0). Zenodo. http://doi.org/10.5281/zenodo.4458777; Hussain, M. and Mahmud, I. (2019). pyMannKendall: a python package for non parametric Mann Kendall family of trend tests. Journal of Open Source Software, 4(39), 1556. doi:10.21105/joss.01556"
	ind_trend['trend'].attrs['missing_value'] = np.nan
	ind_trend['trend'].attrs['valid_min'] = np.nanmin(ind_trend['trend'])
	ind_trend['trend'].attrs['valid_max'] = np.nanmax(ind_trend['trend'])
	ind_trend['trend'].attrs['valid_range'] = [np.nanmin(ind_trend['trend']),np.nanmax(ind_trend['trend'])]

	######################################################

	ind_trend['signif'].attrs['units'] = ""
	ind_trend['signif'].attrs['name'] = 'signif'
	ind_trend['signif'].attrs['long_name'] = "Trends significance"

	ind_trend['signif'].attrs['missing_value'] = np.nan
	ind_trend['signif'].attrs['valid_min'] = np.nanmin(ind_trend['signif'])
	ind_trend['signif'].attrs['valid_max'] = np.nanmax(ind_trend['signif'])
	ind_trend['signif'].attrs['valid_range'] = [np.nanmin(ind_trend['signif']),np.nanmax(ind_trend['signif'])] 

	comp = dict(zlib=True, complevel=5)
	encoding = {var: comp for var in ind_trend.data_vars}
	ind_trend.to_netcdf('trend_'+indicator)

