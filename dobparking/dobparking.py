import pandas as pd
import geopandas as gpd
import numpy as np
import statsmodels.api as sm
from shapely import wkt



pd.set_option('display.max_columns', None)
path='C:/Users/Yijun Ma/Desktop/D/DOCUMENT/DCP2020/DOBPARKING/'
#path='/home/mayijun/DOBPARKING/'



#mappluto=gpd.read_file(path+'nyc_mappluto_20v1_shp/MapPLUTO.shp')
#mappluto=mappluto[['BBL','ZoneDist1','LandUse','BldgClass','LotArea','BldgArea','GarageArea','NumFloors','geometry']].reset_index(drop=True)
#mappluto=mappluto.to_crs({'init':'epsg:4326'})
#mappluto.to_file(path+'mappluto/mappluto.shp')



#parkinglots=gpd.read_file(path+'planimetrics/parkinglots.shp')
#parkinglots.crs={'init':'epsg:4326'}
#mappluto=gpd.read_file(path+'mappluto/mappluto.shp')
#mappluto.crs={'init':'epsg:4326'}
#parkinglots=gpd.overlay(mappluto,parkinglots,'union')
#parkinglots=parkinglots[(pd.notna(parkinglots['BBL']))&(pd.notna(parkinglots['OBJECTID']))].reset_index(drop=True)
#parkinglots=parkinglots[['BBL','geometry']].reset_index(drop=True)
#parkinglots=parkinglots.dissolve(by='BBL').reset_index(drop=False)
#parkinglots.crs={'init':'epsg:4326'}
#parkinglots=parkinglots.to_crs({'init':'epsg:6539'})
#parkinglots['parkinglots']=[x.area for x in parkinglots['geometry']]
#parkinglots=parkinglots.groupby('BBL',as_index=False).agg({'parkinglots':'sum'}).reset_index(drop=True)
#parkinglots.to_csv(path+'parkinglots.csv',index=False)



#buildinggarage=pd.read_csv(path+'building.csv',dtype=str)
#buildinggarage=buildinggarage[buildinggarage['FEAT_CODE']=='5110'].reset_index(drop=True)
#buildinggarage=buildinggarage[pd.notna(buildinggarage['MPLUTO_BBL'])].reset_index(drop=True)
#buildinggarage['HEIGHTROOF']=pd.to_numeric(buildinggarage['HEIGHTROOF'])
#buildinggarage=buildinggarage[pd.notna(buildinggarage['HEIGHTROOF'])].reset_index(drop=True)
#buildinggarage['ESTFLOOR']=[int(x/10) for x in buildinggarage['HEIGHTROOF']]
#buildinggarage['BUILDINGGARAGE']=pd.to_numeric(buildinggarage['SHAPE_AREA'])*buildinggarage['ESTFLOOR']
#buildinggarage['BBL']=pd.to_numeric(buildinggarage['MPLUTO_BBL'])
#buildinggarage=buildinggarage.groupby('BBL',as_index=False).agg({'BUILDINGGARAGE':'sum'}).reset_index(drop=True)
#buildinggarage=buildinggarage[['BBL','BUILDINGGARAGE']].reset_index(drop=True)
#buildinggarage.to_csv(path+'buildinggarage.csv',index=False)



df=pd.read_csv(path+'translatedforyijun.csv',dtype=str)
df['BIN']=pd.to_numeric(df['BIN'])
df['carnumbersum']=pd.to_numeric(df['carnumbersum'])
building=pd.read_csv(path+'building.csv',dtype=str)
building=building[pd.notna(building['MPLUTO_BBL'])].reset_index(drop=True)
building['BIN']=pd.to_numeric(building['BIN'])
building['BBL']=pd.to_numeric(building['MPLUTO_BBL'])
df=pd.merge(df,building,how='inner',on='BIN')
mappluto=gpd.read_file(path+'mappluto/mappluto.shp')
mappluto.crs={'init':'epsg:4326'}
df=pd.merge(df,mappluto,how='inner',on='BBL')
buildinggarage=pd.read_csv(path+'buildinggarage.csv',dtype=float)
df=pd.merge(df,buildinggarage,how='left',on='BBL')
df['BUILDINGGARAGE']=df['BUILDINGGARAGE'].fillna(0)
parkinglots=pd.read_csv(path+'parkinglots.csv',dtype=float)
df=pd.merge(df,parkinglots,how='left',on='BBL')
df['parkinglots']=df['parkinglots'].fillna(0)
df['bc']=[str(x)[0:1] for x in df['BldgClass']]
df['estgaragearea']=[int(x/200) for x in df['GarageArea']]
df['estbuildinggarage']=[int(x/200) for x in df['BUILDINGGARAGE']]
df['estparkinglots']=[int(x/400) for x in df['parkinglots']]
df['estcar']=df['estgaragearea']+df['estbuildinggarage']+df['estparkinglots']
df.to_csv(path+'df.csv',index=False)



df=pd.read_csv(path+'df.csv')

k=df[(df['carnumbersum']>0)&(df['GarageArea']>0)&(df['BUILDINGGARAGE']==0)&(df['parkinglots']==0)].reset_index(drop=True)
# 720
k['rate']=k['GarageArea']/k['carnumbersum']
k=k.groupby('bc').agg({'rate':'median','GarageArea':'count'})

k=df[(df['carnumbersum']>0)&(df['GarageArea']==0)&(df['BUILDINGGARAGE']>0)&(df['parkinglots']==0)].reset_index(drop=True)
# 4367
k['rate']=k['BUILDINGGARAGE']/k['carnumbersum']
k=k.groupby('bc').agg({'rate':'median','BUILDINGGARAGE':'count'})

k=df[(df['carnumbersum']>0)&(df['GarageArea']==0)&(df['BUILDINGGARAGE']==0)&(df['parkinglots']>0)].reset_index(drop=True)
# 1428
k['rate']=k['parkinglots']/k['carnumbersum']
k=k.groupby('bc').agg({'rate':'median','parkinglots':'count'})



k=df[(df['carnumbersum']>0)&((df['GarageArea']>0)|(df['BUILDINGGARAGE']>0)|(df['parkinglots']>0))].reset_index(drop=True)
X = k[['GarageArea','BUILDINGGARAGE','parkinglots']]
y = k['carnumbersum']
model = sm.RLM(y, X).fit()
model.summary()
k['predict']=model.predict()
k.to_csv(path+'k.csv',index=False)








k=df[(df['carnumbersum']>0)&(df['GarageArea']==0)&(df['BUILDINGGARAGE']==0)&(df['parkinglots']==0)].reset_index(drop=True)
k=k[(k['bc']=='A')|(k['bc']=='B')|(k['BldgClass']=='C0')].reset_index(drop=True)
k.to_csv(path+'k.csv',index=False)





X = df[['GarageArea','BUILDINGGARAGE','parkinglots','LotArea']]
y = df['carnumbersum']
model = sm.OLS(y, X).fit()
model.summary()
df['predict']=model.predict()





k=df[(df['carnumbersum']!=0)&((df['parkinglots']>0)|(df['GarageArea']>0))].reset_index(drop=True)
k.to_csv(path+'k.csv',index=False)



k=df[(df['carnumbersum']!=0)&(df['parkinglots']==0)&(df['GarageArea']==0)].reset_index(drop=True)
k['rate']=k['LotArea']/k['carnumbersum']
k=k.groupby('bc').agg({'rate':'median','parkinglots':'count'})


k=k[(k['bc']=='A')].reset_index(drop=True)
k.rate.describe()














