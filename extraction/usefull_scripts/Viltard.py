#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 12:53:11 2023

Code qui co-localise les données de SEVIRI et les données de la mosaic radar MF à 1 km.
Require: satpy, cartopy, netcdf4

@author: viltard
"""
import os
from satpy.scene import Scene
import datetime
from glob import glob
from pyresample.geometry import AreaDefinition

import eumdac
import datetime
import time
import shutil
import requests

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from pyproj.crs import CRS

import random

from netCDF4 import Dataset

exit(0)
consumer_key = '4lFF8IFB4q9OIptYzESHT0oFoyga'
consumer_secret = 'uEbbU4yvMIBpllWiM6RURyDjG8Ma'
credentials = (consumer_key, consumer_secret)
token = eumdac.AccessToken(credentials)
datastore = eumdac.DataStore(token)
selected_collection = datastore.get_collection('EO:EUM:DAT:MSG:HRSEVIRI')

MFGrid  = '/net/nfs/precipitations/nicoEUMETSAT/grid_france_1536-1km_200906110000-219912310000.nc'

nc_fid   = Dataset(MFGrid, 'r')
rLat     = nc_fid.variables['latitude'][:,:]
rLon     = nc_fid.variables['longitude'][:,:]
nc_fid.close()

#for YYYY in ["2018"]:
#    for Month in [10]: #, range (1, 13, 1):
# 2008 => 2021
#for YYYY in ["2009"]:
    #for Month in [ 9 ]: #, range (1, 13, 1):
for YYYY in ["2015", "2011", "2013", "2021", "2012", "2015","2016", "2017", "2018", "2014", "2019", "2020", "2009" ]:
    for Month in [ 9, 3, 8, 4, 5, 6, 12, 1, 2, 10, 11]: #, range (1, 13, 1):
        if (Month < 10):
           MM = "0"+str(Month)
        else:
           MM = str(Month)
        #for Day in [1, 2]:
        for Day in range(1, 32, 2):
           if Day < 10:
             DD = "0"+str(Day)
           else:
             DD = str(Day)
           # -- on traite des jours entiers mais random
           if (random.random() > 0.5):
               continue
           # -- reading the MF Mosaic
           Quality = 0
           FileList = glob('/net/nfs/precipitations/SEDOO/Mosaiques/LMDO_1536_1km_5min/'+YYYY+'/'+YYYY+MM+'/'+YYYY+MM+DD+'/'+'cumul_france_1536-1km-5min_'+YYYY+MM+DD+'*.nc')
           if (len(FileList) < 1): 
               print("Lapin trouvé..8-(")
               continue

           for iFile,File in enumerate(FileList):
              print(str(iFile)+"/"+str(len(FileList)))
              nc_fid = Dataset(File, 'r')
              LRR          = nc_fid.variables['ACRR'][:,:,:]
              QRR          = nc_fid.variables['QUALITY'][:,:,:]
              LRR          = LRR[0,:,:].astype(float)
              QRR          = QRR[0,:,:].astype(int)
              Time         = nc_fid.variables['time']
              # -- ACCRR en 1/100 mm et pour 5 minutes, donc en mm/hr:
              # -- si LRR est hors limite ou si le QRR est en-dessous du seuil, on met à -99 <=> missing
              LRR[LRR > 65000]   = -99.
              LRR[QRR > 200]     = -98.
              LRRFull            = np.copy(np.transpose((LRR)))
              LRR[QRR < Quality] = -97.
              #-- conversion en de 1/100 mm sur 5 minutes en mm/hr
              LRR[LRR > 0.]        = LRR[LRR > 0.]/100.*12.
              LRRFull[LRRFull > 0] = LRRFull[LRRFull > 0]/100.*12.

              # -- on s�lectionne que des images radars qui ont de la pluie...
              #    si il y a moins de 100 pixels de qualite 80 superieurs a 1 mmm/hr on passe...
              if (np.sum((LRR >= 1.) & (QRR >= 80.)) < 100.): 
                  continue

              Time       = nc_fid.variables['time']
              Dum        = getattr(Time, 'units')
              iDum       = [int(Dum[14:18]), int(Dum[19:21]), int(Dum[22:24]),
                            int(Dum[25:27]), int(Dum[28:30]), int(Dum[31:33])]
              Time       = Time[0]
              tref       = datetime.datetime(iDum[0], iDum[1], iDum[2], iDum[3], iDum[4], iDum[5])
              Mosaic_time = tref+datetime.timedelta(seconds=int(Time))
              print(Mosaic_time)

              nc_fid.close()

              mi       = str(Mosaic_time)
              # -- 20240824: on ne garde que les trucs qui correspondent aux images SEVIRI (toutes les 15')... 
              #    toutefois, la mosaique a -5 minutes could be best correletd wi the SEVIRI 
              #    image car le balayage se fait du nord au sud sur 10'.
              if ((mi[14:16] == "55") or (mi[14:16] == "10") or (mi[14:16] == "25") or (mi[14:16] == "40")):
                 Fileout  =  'sevmos_'+mi[:10]+'_'+mi[11:13]+mi[14:16]+mi[17:19]+'.nc'
                 DumOut = glob('/net/nfs/precipitations/MatthieuM/DATA/'+Fileout)
                 if (len(DumOut) >= 1):
                    print(Fileout," existe deja, on boucle au suivant")
                    continue

                 Lat = np.copy(rLat)
                 Lon = np.copy(rLon)

                 Good  = np.where(LRR > 0., True, False)
                 nkeep = np.sum(Good)
                 LRR   = np.reshape(LRR[Good], nkeep)
                 Lat   = np.reshape(Lat[Good], nkeep)
                 Lon   = np.reshape(Lon[Good], nkeep)

                 start = Mosaic_time+datetime.timedelta(minutes=-7.49)
                 end   = Mosaic_time+datetime.timedelta(minutes=7.50)
                 #print(start)
                 #print(end)

                 try: 
                    products = selected_collection.search(dtstart=start, dtend=end)
                    print("found "+str(len(products)))
                 except:
                    continue  
                
                 dttime = []
                 for product in products:
                     print(str(product))
                     Dum  = str(product)
                     iDum = [int(Dum[24:28]), int(Dum[28:30]), int(Dum[30:32]),
                             int(Dum[32:34]), int(Dum[34:36]), int(Dum[36:38])]
                     tsat = datetime.datetime(iDum[0], iDum[1], iDum[2], iDum[3], iDum[4], iDum[5])
                     dttime.append(np.abs(tsat-Mosaic_time))
              
                 print(dttime)
                 ichoice = np.argmin(dttime)
                  
                 for iprod, product in enumerate(products):
                    if ( iprod != ichoice): continue
                    print(str(product))
                    Fail_Download = False
                    selected_product = datastore.get_product(product_id=str(product),
                                                             collection_id='EO:EUM:DAT:MSG:HRSEVIRI')
                    if (len(selected_product.entries) > 1):
                      for entries in selected_product.entries:
                         if ('.nat' in entries):
                            with selected_product.open(entry=entries) as fsrc, open(fsrc.name, mode='wb') as fdst:
                                try :
                                   shutil.copyfileobj(fsrc, fdst)
                                   time.sleep(3)
                                   print(f'Download of product {fsrc.name} finished.')
                                except:
                                   print(f'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                                   print(f'Download of product {fsrc.name} was unsuccessful !!!!!!!!!')
                                   print(f'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                                   Fail_Download = True
                                   continue

                 if (Fail_Download):
                     continue

                 filenames = [fsrc.name]
                 print('All downloads are finished.')
                 MSG       = Scene(reader="seviri_l1b_native", filenames=filenames)
                 # 'IR_016',
                 # 'IR_039',
                 # 'IR_087',
                 # 'IR_097',
                 # 'IR_108',
                 # 'IR_120',
                 # 'IR_134',
                 # 'VIS006',
                 # 'VIS008',
                 # 'WV_062',
                 # 'WV_073',
                 DSName = ['WV_062','WV_073','IR_087','IR_097','IR_108', 'IR_120','IR_134',]
                 MSG.load(DSName)

                 #%% -- Resampling the SEVIRI data to MF domain
                 mf_france = AreaDefinition.from_cf(MFGrid, variable="grid_mapping", x='X', y='Y')
   
                 local_MSG  = MSG.resample(mf_france)
                 #Masked_MSG = local_MSG.where(LRRFull >= 0)
                 #local_MSG.show('IR_108')
                 proj_crs = CRS.from_epsg(4326)
                 cart_crs = ccrs.CRS(proj_crs)
                 crs = local_MSG[DSName[1]].attrs['area'].to_cartopy_crs()
                 """
                 #%% Plotting SEVIRI Tb and radar data
                 for DS in DSName:
                     fig = plt.figure(dpi=300)
                     ax = plt.axes(projection=crs)
                     ax.coastlines(resolution="50m", color="yellow", linewidth=0.2)
                     gl = ax.gridlines(draw_labels=True)
                     gl.top_labels = False
                     gl.right_labels = False
                     ax.set_global()
                     plt.imshow(local_MSG[DS], transform=crs, extent=crs.bounds, cmap='gist_heat_r')
                     cbar = plt.colorbar()
                     cbar.set_label("[K]")
                     plt.scatter(x=Lon, y=Lat, c=LRR, transform=cart_crs, cmap="gist_ncar", marker='s', s=0.01, alpha=1)
                     cbar = plt.colorbar()
                     cbar.set_label("[mm/hr]")
                     plt.title(DS+" "+str(Mosaic_time))
                     plt.show()
                 """
                 #%% Writing netcdf with extracted image
                 local_MSG.save_datasets(writer='cf', datasets=DSName, filename=Fileout,
                                         exclude_attrs=['raw_metadata'], 
                                         base_dir='/net/nfs/precipitations/nicoEUMETSAT/DATA')

                 Fileout                   = '/net/nfs/precipitations/nicoEUMETSAT/DATA/'+Fileout
                 w_nc_fid                  = Dataset(Fileout, 'a')
                 w_nc_fid.Contact          = "nicolas.viltard@latmos.ipsl.fr"
                 w_nc_fid.description      = "Generated with satpytest_MkV.py"
                 w_nc_fid.MSGFile          = filenames[0]
                 w_nc_fid.MFMosaicFile     = File
                 w_nc_fid.Comment1         = "First attempt..."
                 #w_nc_fid.Comment2         = "    -na"
                 w_nc_fid.history          = "Created " + str(datetime.datetime.now())
                 w_nc_fid.version          = "v0.1"

                 #--------------------            
                 w_nc_var   = w_nc_fid.createVariable('rain_rate', 'f', ('x', 'y'))
                 MinMax = (np.min(LRRFull),np.max(LRRFull))
                 w_nc_var.setncatts({'long_name': "MF mosaic 5' rain accumulation converted to rain rate", 
                                     'mosaic ref. time':str(Mosaic_time),
                                     'units': "mm/hr", 'valid_range':MinMax, 'missing_value':-99.}) #, '_FillValue':-888})
                 w_nc_fid.variables['rain_rate'][:,:] = LRRFull

                 w_nc_var   = w_nc_fid.createVariable('rain_quality', 'i4', ('x', 'y'))
                 MinMax = (np.min(LRRFull),np.max(LRRFull))
                 w_nc_var.setncatts({'long_name': "Quality: 100= good, 0 = bad",
                                     'units': "percent", 'valid_range':(0, 100), 'missing_value':255})
                 w_nc_fid.variables['rain_quality'][:,:] = QRR
   
                 w_nc_fid.close()

                 for File in filenames:
                    os.remove(File)
