import geopandas as gpd
import pandas as pd
import numpy as np
import shapely
import re
import datetime
from shapely import wkt
from geosupport import Geosupport



pd.set_option('display.max_columns', None)
#path='C:/Users/Yijun Ma/Desktop/D/DOCUMENT/DCP2019/ONSTPARKING/'
#path='C:/Users/Y_Ma2/Desktop/ONSTPARKING/'
#path='E:/ONSTPARKING/'
#path='D:/ONSTPARKING/'
path='I:/ONSTPARKING/'



def calcpvmtbearing(pvbr):
    latfrom=np.radians(pvbr['latfrom'])
    latto=np.radians(pvbr['latto'])
    difflong=np.radians(pvbr['longto']-pvbr['longfrom'])
    pvbearing=np.degrees(np.arctan2(np.sin(difflong)*np.cos(latfrom),
                                      np.cos(latfrom)*np.sin(latto)-(np.sin(latfrom)*np.cos(latto)*np.cos(difflong))))
    pvbr['pvbearing']=pvbearing
    return pvbr

def calcpvmtdir(pvdir):
    if (pvdir['pvbearing']>=0)&(pvdir['pvbearing']<90):
        pvdir['pvdir']='N/E'
    elif (pvdir['pvbearing']>=90)&(pvdir['pvbearing']<=180):
        pvdir['pvdir']='S/E'
    elif (pvdir['pvbearing']<0)&(pvdir['pvbearing']>=-90):
        pvdir['pvdir']='N/W'
    elif (pvdir['pvbearing']<-90)&(pvdir['pvbearing']>=-180):
        pvdir['pvdir']='S/W'
    return pvdir

def adjhydrant(hd):
    hdbk=pvmtedgemdn[np.isin(pvmtedgemdn['bkfaceid'],hydrantbuffer.loc[hydrantbuffer['unitid']==hd['unitid'],'bkfaceid'])].reset_index(drop=True)
    if len(hdbk)>0:
        hd['bkfaceid']=hdbk.loc[np.argmin([hd['geometry'].distance(x) for x in hdbk['geometry']]),'bkfaceid']
        hd['snapdist']=min([hd['geometry'].distance(x) for x in hdbk['geometry']])
        hd['adjgeom']=shapely.ops.nearest_points(hd['geometry'],hdbk.loc[np.argmin([hd['geometry'].distance(x) for x in hdbk['geometry']]),'geometry'])[1].wkt
    else:
        print(str(hd['unitid'])+' no bkfaceid joined')
    return hd

def adjdotshp(ds):
    dsbk=pvmtedgemdn[np.isin(pvmtedgemdn['bkfaceid'],dotshpbuffer.loc[(dotshpbuffer['order']==ds['order'])&(dotshpbuffer['seq']==ds['seq']),'bkfaceid'])].reset_index(drop=True)
    if len(dsbk)>0:
        ds['bkfaceid']=dsbk.loc[np.argmin([ds['geometry'].distance(x) for x in dsbk['geometry']]),'bkfaceid']
        ds['snapdist']=min([ds['geometry'].distance(x) for x in dsbk['geometry']])
        ds['adjgeom']=shapely.ops.nearest_points(ds['geometry'],dsbk.loc[np.argmin([ds['geometry'].distance(x) for x in dsbk['geometry']]),'geometry'])[1].wkt
    else:
        print(str(ds['order'])+'|'+str(ds['seq'])+' no bkfaceid joined')
    return ds

def sepdescmutcd(dm):
    if pd.isna(dm['Sign_description']):
        dm['desc']=''
        dm['mutcd']=' '.join(dm['SR_Mutcd_Code'].split()).upper()
    elif pd.isna(dm['Unnamed: 7']):
        dm['desc']=' '.join(dm['Sign_description'].split()).upper()
        dm['mutcd']=' '.join(dm['SR_Mutcd_Code'].split()).upper()
    elif pd.isna(dm['Unnamed: 8']):
        dm['desc']=' '.join((dm['Sign_description']+', '+dm['SR_Mutcd_Code']).split()).upper()
        dm['mutcd']=' '.join(dm['Unnamed: 7'].split()).upper()
    else:
        dm['desc']=' '.join((dm['Sign_description']+', '+dm['SR_Mutcd_Code']+', '+dm['Unnamed: 7']).split()).upper()
        dm['mutcd']=' '.join(dm['Unnamed: 8'].split()).upper()
    return dm

def sepdesc(sds):
    if sds['descnum']==0:
        sds['desc1']=' '.join(sds['desc'].split()).upper()
        sds['dir1']=np.where(len(list(re.finditer('ARROW',sds['desc1'])))>0,1,2)
    elif sds['descnum']==1:
        descarrow=[x.start() for x in list(re.finditer('>',sds['desc']))]
        sds['desc1']=' '.join(sds['desc'][:descarrow[0]+1].split()).upper()
        sds['dir1']=np.where(len(list(re.finditer('<',sds['desc1'])))>0,2,1)
        sds['desctail']=' '.join(sds['desc'][descarrow[0]+1:].split()).upper()
    elif sds['descnum']==2:
        descarrow=[x.start() for x in list(re.finditer('>',sds['desc']))]
        sds['desc1']=' '.join(sds['desc'][:descarrow[0]+1].split()).upper()
        sds['dir1']=np.where(len(list(re.finditer('<',sds['desc1'])))>0,2,1)
        sds['desc2']=' '.join(sds['desc'][descarrow[0]+1:descarrow[1]+1].split()).upper()
        sds['dir2']=np.where(len(list(re.finditer('<',sds['desc2'])))>0,2,1)
        sds['desctail']=' '.join(sds['desc'][descarrow[1]+1:].split()).upper()
    elif sds['descnum']==3:
        descarrow=[x.start() for x in list(re.finditer('>',sds['desc']))]
        sds['desc1']=' '.join(sds['desc'][:descarrow[0]+1].split()).upper()
        sds['dir1']=np.where(len(list(re.finditer('<',sds['desc1'])))>0,2,1)
        sds['desc2']=' '.join(sds['desc'][descarrow[0]+1:descarrow[1]+1].split()).upper()
        sds['dir2']=np.where(len(list(re.finditer('<',sds['desc2'])))>0,2,1)
        sds['desc3']=' '.join(sds['desc'][descarrow[1]+1:descarrow[2]+1].split()).upper()
        sds['dir3']=np.where(len(list(re.finditer('<',sds['desc3'])))>0,2,1)
        sds['desctail']=' '.join(sds['desc'][descarrow[2]+1:].split()).upper()
    return sds

def extractdays(ed):
    ed=ed.copy()
    if pd.notna(re.search('ANYTIME',ed['time'])):
        ed[['m','t','w','r','f','s','u']]+=1
        ed['daysflag']+=1  
    if pd.notna(re.search('ALL DAYS',ed['time'])):
        ed[['m','t','w','r','f','s','u']]=1
        ed['daysflag']+=1
    if pd.notna(re.search('HOTEL LOADING ZONE',ed['time'])):
        ed[['m','t','w','r','f','s','u']]+=1
        ed['daysflag']+=1
    if pd.notna(re.search('FIRE ZONE',ed['time'])):
        ed[['m','t','w','r','f','s','u']]+=1
        ed['daysflag']+=1        
    if ed['time'].replace('<','').replace('>','').replace('-','').replace(' ','')=='':
        ed[['m','t','w','r','f','s','u']]+=1
        ed['daysflag']+=1    
    if pd.notna(re.search('SCHOOL DAYS',ed['time'])):
        ed[['m','t','w','r','f','s','u']]+=1
        ed['daysflag']+=1
    if pd.notna(re.search('MONDAY',ed['time'])):
        ed['m']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('MON',ed['time'])):
        ed['m']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('TUESDAY',ed['time'])):
        ed['t']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('TUES',ed['time'])):
        ed['t']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('WEDNESDAY',ed['time'])):
        ed['w']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('WED',ed['time'])):
        ed['w']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('THURSDAY',ed['time'])):
        ed['r']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('THURS',ed['time'])):
        ed['r']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('FRIDAY',ed['time'])):
        ed['f']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('FRI',ed['time'])):
        ed['f']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('SATURDAY',ed['time'])):
        ed['s']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('SAT',ed['time'])):
        ed['s']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('SUNDAY',ed['time'])):
        ed['u']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('SUN',ed['time'])):
        ed['u']+=1
        ed['daysflag']+=1
    if pd.notna(re.search('MONDAY-FRIDAY',ed['time'])):
        ed[['m','t','w','r','f']]+=1
        ed['daysflag']+=1
    if pd.notna(re.search('MON-FRI',ed['time'])):
        ed[['m','t','w','r','f']]+=1
        ed['daysflag']+=1
    if pd.notna(re.search('MON THRU FRI',ed['time'])):
        ed[['m','t','w','r','f']]+=1
        ed['daysflag']+=1
    if pd.notna(re.search('EXCEPT SUNDAY',ed['time'])):
        ed[['m','t','w','r','f','s']]+=1
        ed['u']=0
        ed['daysflag']+=1
    if pd.notna(re.search('EXCEPT SATURDAY',ed['time'])):
        ed[['m','t','w','r','f','u']]+=1
        ed['s']=0
        ed['daysflag']+=1
    if pd.notna(re.search('INCLUDING SUNDAY',ed['time'])):
        ed[['m','t','w','r','f','s','u']]+=1
        ed['daysflag']+=1
    return ed

def extracthours(eh):
    eh=eh.copy()
    if pd.notna(re.search('\d+AM-\d+AM',eh['time'])):
        eh['starthour']=str(int(re.search('\d+AM-\d+AM',eh['time'])[0].split('-')[0].replace('AM',''))).zfill(2)+':00'
        eh['endhour']=str(int(re.search('\d+AM-\d+AM',eh['time'])[0].split('-')[1].replace('AM',''))).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+AM-\d+PM',eh['time'])):
        eh['starthour']=str(int(re.search('\d+AM-\d+PM',eh['time'])[0].split('-')[0].replace('AM',''))).zfill(2)+':00'
        eh['endhour']=str(int(re.search('\d+AM-\d+PM',eh['time'])[0].split('-')[1].replace('PM',''))+12).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+PM-\d+AM',eh['time'])):
        eh['starthour']=str(int(re.search('\d+PM-\d+AM',eh['time'])[0].split('-')[0].replace('PM',''))+12).zfill(2)+':00'
        eh['endhour']=str(int(re.search('\d+PM-\d+AM',eh['time'])[0].split('-')[1].replace('AM',''))).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+PM-\d+PM',eh['time'])):
        eh['starthour']=str(int(re.search('\d+PM-\d+PM',eh['time'])[0].split('-')[0].replace('PM',''))+12).zfill(2)+':00'
        eh['endhour']=str(int(re.search('\d+PM-\d+PM',eh['time'])[0].split('-')[1].replace('PM',''))+12).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+-\d+AM',eh['time'])):
        eh['starthour']=str(int(re.search('\d+-\d+AM',eh['time'])[0].split('-')[0])).zfill(2)+':00'
        eh['endhour']=str(int(re.search('\d+-\d+AM',eh['time'])[0].split('-')[1].replace('AM',''))).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+-\d+PM',eh['time'])):
        eh['starthour']=str(int(re.search('\d+-\d+PM',eh['time'])[0].split('-')[0])+12).zfill(2)+':00'
        eh['endhour']=str(int(re.search('\d+-\d+PM',eh['time'])[0].split('-')[1].replace('PM',''))+12).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+AM-MIDNIGHT',eh['time'])):
        eh['starthour']=str(int(re.search('\d+AM-MIDNIGHT',eh['time'])[0].split('-')[0].replace('AM',''))).zfill(2)+':00'
        eh['endhour']=str(0).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+PM-MIDNIGHT',eh['time'])):
        eh['starthour']=str(int(re.search('\d+PM-MIDNIGHT',eh['time'])[0].split('-')[0].replace('PM',''))+12).zfill(2)+':00'
        eh['endhour']=str(0).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('MIDNIGHT-\d+AM',eh['time'])):
        eh['starthour']=str(0).zfill(2)+':00'
        eh['endhour']=str(int(re.search('MIDNIGHT-\d+AM',eh['time'])[0].split('-')[1].replace('AM',''))).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('MIDNIGHT-\d+PM',eh['time'])):
        eh['starthour']=str(0).zfill(2)+':00'
        eh['endhour']=str(int(re.search('MIDNIGHT-\d+PM',eh['time'])[0].split('-')[1].replace('PM',''))+12).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+AM-NOON',eh['time'])):
        eh['starthour']=str(int(re.search('\d+AM-NOON',eh['time'])[0].split('-')[0].replace('AM',''))).zfill(2)+':00'
        eh['endhour']=str(12).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('\d+PM-NOON',eh['time'])):
        eh['starthour']=str(int(re.search('\d+PM-NOON',eh['time'])[0].split('-')[0].replace('PM',''))+12).zfill(2)+':00'
        eh['endhour']=str(12).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('NOON-\d+AM',eh['time'])):
        eh['starthour']=str(12).zfill(2)+':00'
        eh['endhour']=str(int(re.search('NOON-\d+AM',eh['time'])[0].split('-')[1].replace('AM',''))).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('NOON-\d+PM',eh['time'])):
        eh['starthour']=str(12).zfill(2)+':00'
        eh['endhour']=str(int(re.search('NOON-\d+PM',eh['time'])[0].split('-')[1].replace('PM',''))+12).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('ANYTIME',eh['time'])):
        eh['starthour']=str(0).zfill(2)+':00'
        eh['endhour']=str(0).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('FIRE ZONE',eh['time'])):
        eh['starthour']=str(0).zfill(2)+':00'
        eh['endhour']=str(0).zfill(2)+':00'
        eh['hoursflag']+=1
    if pd.notna(re.search('HOTEL LOADING ZONE',eh['time'])):
        eh['starthour']=str(0).zfill(2)+':00'
        eh['endhour']=str(0).zfill(2)+':00'
        eh['hoursflag']+=1
    if eh['time'].replace('<','').replace('>','').replace('-','').replace(' ','')=='':
        eh['starthour']=str(0).zfill(2)+':00'
        eh['endhour']=str(0).zfill(2)+':00'
        eh['hoursflag']+=1
    return eh

def checkdotshpadjgeom(dsg):
    dsg=dsg.reset_index(drop=True)
    dsgod=list(dsg['order'])[0]
    if len(dsg)>=2:
        latfrom=np.radians(dsg.loc[0,'lat'])
        latto=np.radians(dsg.loc[len(dsg)-1,'lat'])
        difflong=np.radians(dsg.loc[len(dsg)-1,'long']-dsg.loc[0,'long'])
        dsgbearing=np.degrees(np.arctan2(np.sin(difflong)*np.cos(latfrom),
                                         np.cos(latfrom)*np.sin(latto)-(np.sin(latfrom)*np.cos(latto)*np.cos(difflong))))
        dsg=dsg.loc[[0]].reset_index(drop=True)
        dsg.loc[0,'bearing']=dsgbearing
        dsg.loc[0,'bearingdiff']=abs(dsgbearing-dsg.loc[0,'pvbearing'])
        if (dsg.loc[0,'bearingdiff']>90)&(dsg.loc[0,'bearingdiff']<=270):
            dsg.loc[0,'geom']='LINESTRING ('+', '.join(dsg.loc[0,'geometry'].to_wkt().replace('LINESTRING (','').replace(')','').split(', ')[::-1])+')'
            if dsg.loc[0,'pvdir']=='N/E':
                dsg.loc[0,'pvdir']=='S/W'
            elif dsg.loc[0,'pvdir']=='S/E':
                dsg.loc[0,'pvdir']=='N/W'
            elif dsg.loc[0,'pvdir']=='N/W':
                dsg.loc[0,'pvdir']=='S/E'
            elif dsg.loc[0,'pvdir']=='S/W':
                dsg.loc[0,'pvdir']=='N/E'
        elif (dsg.loc[0,'bearingdiff']<=90)|(dsg.loc[0,'bearingdiff']>270):
            dsg.loc[0,'geom']=dsg.loc[0,'geometry'].to_wkt()
        dsg=dsg[['order','bearing','bkfaceid','pvbearing','pvdir','shplen','bearingdiff','geom']].reset_index(drop=True)
    else:
        print(dsgod+' has only 1 sign')
        dsg=pd.DataFrame(columns=['order','bearing','bkfaceid','pvbearing','pvdir','shplen','bearingdiff','geom'])
    return dsg

def sumsign(ss):
    try:
        ssod=list(ss['order'])[0]
        ss=ss.reset_index(drop=True)
        ss['dist2']=np.roll(ss['dist'],-1)
        ss['arrow2']=np.roll(ss['arrow'],-1)
        ss['cl2']=np.roll(ss['cl'],-1)
        ss['blpl2']=np.roll(ss['blpl'],-1)
        ss['t09302']=0
        ss['t20302']=0
        ss['w16002']=0
        ss=ss.loc[:len(ss)-2].reset_index(drop=True)
        if (ss.loc[0,'cl']==1)&(ss.loc[0,'blpl2']==1)&(ss.loc[len(ss)-1,'blpl']==1&(ss.loc[len(ss)-1,'cl2']==1)):
            ss.loc[0,'t09302']=1
            ss.loc[0,'t20302']=1
            ss.loc[0,'w16002']=1
            ss.loc[len(ss)-1,'t09302']=1
            ss.loc[len(ss)-1,'t20302']=1
            ss.loc[len(ss)-1,'w16002']=1
            for i in ss.index[2:-1]:
                if ss.loc[i,'arrow']=='':
                    ss.loc[ss['dist']==ss['dist'].unique()[list(ss['dist'].unique()).index(ss.loc[i,'dist'])-1],'t09302']+=ss.loc[i,'t0930']
                    ss.loc[ss['dist']==ss['dist'].unique()[list(ss['dist'].unique()).index(ss.loc[i,'dist'])-1],'t20302']+=ss.loc[i,'t2030']
                    ss.loc[ss['dist']==ss['dist'].unique()[list(ss['dist'].unique()).index(ss.loc[i,'dist'])-1],'w16002']+=ss.loc[i,'w1600']
                    ss.loc[ss['dist']==ss.loc[i,'dist'],'t09302']+=ss.loc[i,'t0930']
                    ss.loc[ss['dist']==ss.loc[i,'dist'],'t20302']+=ss.loc[i,'t2030']
                    ss.loc[ss['dist']==ss.loc[i,'dist'],'w16002']+=ss.loc[i,'w1600']
                elif ss.loc[i,'arrow'][0] in ss.loc[i,'pvdir']:
                    ss.loc[ss['dist']==ss.loc[i,'dist'],'t09302']+=ss.loc[i,'t0930']
                    ss.loc[ss['dist']==ss.loc[i,'dist'],'t20302']+=ss.loc[i,'t2030']
                    ss.loc[ss['dist']==ss.loc[i,'dist'],'w16002']+=ss.loc[i,'w1600']
                else:
                    ss.loc[ss['dist']==ss['dist'].unique()[list(ss['dist'].unique()).index(ss.loc[i,'dist'])-1],'t09302']+=ss.loc[i,'t0930']
                    ss.loc[ss['dist']==ss['dist'].unique()[list(ss['dist'].unique()).index(ss.loc[i,'dist'])-1],'t20302']+=ss.loc[i,'t2030']
                    ss.loc[ss['dist']==ss['dist'].unique()[list(ss['dist'].unique()).index(ss.loc[i,'dist'])-1],'w16002']+=ss.loc[i,'w1600']
            ss=ss.loc[ss['dist']!=ss['dist2'],['order','dist','dist2','t09302','t20302','w16002']].reset_index(drop=True)
            ss.columns=['order','dist1','dist2','t0930','t2030','w1600']
            ss['dist1']=ss['dist1']/ss.loc[len(ss)-1,'dist2']
            ss['dist2']=ss['dist2']/ss.loc[len(ss)-1,'dist2']
            return ss
    except:
        print(ssod+' ERROR')

def hydrantsign(hs):
    try:
        hsod=list(hs['order'])[0]
        hs=hs.reset_index(drop=True)
        hs['hdt']=0
        hsgm=list(dotshpadjclean.loc[dotshpadjclean['order']==hsod,'geometry'])[0]
        hdt=hydrantadj[hydrantadj['bkfaceid']==list(dotshpadjclean.loc[dotshpadjclean['order']==hsod,'bkfaceid'])[0]].reset_index(drop=True)
        hdt=[hsgm.project(x,normalized=True) for x in hdt['geometry']]
        for i in hdt:
            hdtmin=hs[(hs['dist1']<=i-15/hsgm.length)&(hs['dist2']>=i-15/hsgm.length)]
            hdtmax=hs[(hs['dist1']<=i+15/hsgm.length)&(hs['dist2']>=i+15/hsgm.length)]
            if (len(hdtmin)==0)&(len(hdtmax)==0):
                hs['hdt']=1
            elif len(hdtmin)==0:
                hdtminmax=pd.concat([hs[:hdtmax.index[0]+1],hdtmax],axis=0,ignore_index=True)
                hdtminmax.loc[len(hdtminmax)-2,'dist2']=i+15/hsgm.length
                hdtminmax.loc[len(hdtminmax)-1,'dist1']=i+15/hsgm.length
                hdtminmax.loc[:len(hdtminmax)-2,'hdt']=1
                hs=hs[hdtmax.index[0]+1:].reset_index(drop=True)
                hs=pd.concat([hs,hdtminmax],axis=0,ignore_index=True)
                hs=hs.sort_values(['dist1','dist2'],ascending=True).reset_index(drop=True)
            elif len(hdtmax)==0:
                hdtminmax=pd.concat([hdtmin,hs[hdtmin.index[0]:]],axis=0,ignore_index=True)
                hdtminmax.loc[0,'dist2']=i-15/hsgm.length
                hdtminmax.loc[1,'dist1']=i-15/hsgm.length
                hdtminmax.loc[1:,'hdt']=1
                hs=hs[:hdtmin.index[0]].reset_index(drop=True)
                hs=pd.concat([hs,hdtminmax],axis=0,ignore_index=True)
                hs=hs.sort_values(['dist1','dist2'],ascending=True).reset_index(drop=True)
            elif hdtmax.index-hdtmin.index>=0:
                hdtminmax=pd.concat([hdtmin,hs[hdtmin.index[0]:hdtmax.index[0]+1],hdtmax],axis=0).reset_index(drop=True)
                hdtminmax.loc[0,'dist2']=i-15/hsgm.length
                hdtminmax.loc[1,'dist1']=i-15/hsgm.length
                hdtminmax.loc[len(hdtminmax)-2,'dist2']=i+15/hsgm.length
                hdtminmax.loc[len(hdtminmax)-1,'dist1']=i+15/hsgm.length
                hdtminmax.loc[1:len(hdtminmax)-2,'hdt']=1
                hs=pd.concat([hs[:hdtmin.index[0]],hs[hdtmax.index[0]+1:]],axis=0,ignore_index=True)
                hs=pd.concat([hs,hdtminmax],axis=0,ignore_index=True)
                hs=hs.sort_values(['dist1','dist2'],ascending=True).reset_index(drop=True)
        return hs
    except:
        print(hsod+' ERROR')

def splitgm(sg):
    try:
        sgod=list(sg['order'])[0]
        sg=sg.reset_index(drop=True)
        sggm=list(dotshpadjclean.loc[dotshpadjclean['order']==sgod,'geometry'])[0]
        splitpos=[x for x in list(sg['dist1'])+[list(sg['dist2'])[-1]]]
        splitter=shapely.geometry.MultiPoint([sggm.interpolate(x,normalized=True) for x in splitpos])
        shapesplit=shapely.ops.split(sggm,splitter.buffer(1e-8))
        shapesplit=[shapesplit[x].wkt for x in range(1,len(shapesplit),2)]
        sg['geom']=shapesplit
        return sg
    except:
        print(sgod+' ERROR')

def lionpkbkface(pb):
    if pb['trafficdir']=='W':
        if pb['parkinglane']==1:
            pb['lparking']=0
            pb['rparking']=1
        elif pb['parkinglane'] in [2,4]:
            pb['lparking']=1
            pb['rparking']=1
    elif pb['trafficdir']=='A':
        if pb['parkinglane']==1:
            pb['lparking']=1
            pb['rparking']=0
        elif pb['parkinglane'] in [2,4]:
            pb['lparking']=1
            pb['rparking']=1
    elif pb['trafficdir']=='T':
        if pb['parkinglane']==1:
            pb['lparking']=0.5
            pb['rparking']=0.5
        elif pb['parkinglane'] in [2,4]:
            pb['lparking']=1
            pb['rparking']=1
    return pb

def lionpkhydrant(lh):
    try:
        lhbk=list(lh['bkfaceid'])[0]
        lh=lh.reset_index(drop=True)
        lhgm=list(pvmtedgemdn.loc[pvmtedgemdn['bkfaceid']==lhbk,'geometry'])[0]
        lh=pd.concat([lh,lh,lh],axis=0,ignore_index=True)
        lh.loc[0,'dist1']=0
        lh.loc[0,'dist2']=15/lhgm.length
        lh.loc[1,'dist1']=15/lhgm.length
        lh.loc[1,'dist2']=(lhgm.length-15)/lhgm.length 
        lh.loc[2,'dist1']=(lhgm.length-15)/lhgm.length
        lh.loc[2,'dist2']=1
        lhdt=hydrantadj[hydrantadj['bkfaceid']==lhbk].reset_index(drop=True)
        lhdt=[lhgm.project(x,normalized=True) for x in lhdt['geometry']]
        for i in lhdt:
            lhdtmin=lh[(lh['dist1']<=i-15/lhgm.length)&(lh['dist2']>=i-15/lhgm.length)]
            lhdtmax=lh[(lh['dist1']<=i+15/lhgm.length)&(lh['dist2']>=i+15/lhgm.length)]
            if (len(lhdtmin)==0)&(len(lhdtmax)==0):
                lh['hdt']=1
            elif len(lhdtmin)==0:
                lhdtminmax=pd.concat([lh[:lhdtmax.index[0]+1],lhdtmax],axis=0,ignore_index=True)
                lhdtminmax.loc[len(lhdtminmax)-2,'dist2']=i+15/lhgm.length
                lhdtminmax.loc[len(lhdtminmax)-1,'dist1']=i+15/lhgm.length
                lhdtminmax.loc[:len(lhdtminmax)-2,'hdt']=1
                lh=lh[lhdtmax.index[0]+1:].reset_index(drop=True)
                lh=pd.concat([lh,lhdtminmax],axis=0,ignore_index=True)
                lh=lh.sort_values(['dist1','dist2'],ascending=True).reset_index(drop=True)
            elif len(lhdtmax)==0:
                lhdtminmax=pd.concat([lhdtmin,lh[lhdtmin.index[0]:]],axis=0,ignore_index=True)
                lhdtminmax.loc[0,'dist2']=i-15/lhgm.length
                lhdtminmax.loc[1,'dist1']=i-15/lhgm.length
                lhdtminmax.loc[1:,'hdt']=1
                lh=lh[:lhdtmin.index[0]].reset_index(drop=True)
                lh=pd.concat([lh,lhdtminmax],axis=0,ignore_index=True)
                lh=lh.sort_values(['dist1','dist2'],ascending=True).reset_index(drop=True)
            elif lhdtmax.index-lhdtmin.index>=0:
                lhdtminmax=pd.concat([lhdtmin,lh[lhdtmin.index[0]:lhdtmax.index[0]+1],lhdtmax],axis=0).reset_index(drop=True)
                lhdtminmax.loc[0,'dist2']=i-15/lhgm.length
                lhdtminmax.loc[1,'dist1']=i-15/lhgm.length
                lhdtminmax.loc[len(lhdtminmax)-2,'dist2']=i+15/lhgm.length
                lhdtminmax.loc[len(lhdtminmax)-1,'dist1']=i+15/lhgm.length
                lhdtminmax.loc[1:len(lhdtminmax)-2,'hdt']=1
                lh=pd.concat([lh[:lhdtmin.index[0]],lh[lhdtmax.index[0]+1:]],axis=0,ignore_index=True)
                lh=pd.concat([lh,lhdtminmax],axis=0,ignore_index=True)
                lh=lh.sort_values(['dist1','dist2'],ascending=True).reset_index(drop=True)
        return lh
    except:
        print(str(lhbk)+' ERROR')

def lionpksplitgeom(lsg):
    try:
        lsgbk=list(lsg['bkfaceid'])[0]
        lsg=lsg.reset_index(drop=True)
        lsggm=list(pvmtedgemdn.loc[pvmtedgemdn['bkfaceid']==lsgbk,'geometry'])[0]
        lsplitpos=[x for x in list(lsg['dist1'])+[list(lsg['dist2'])[-1]]]
        lsplitter=shapely.geometry.MultiPoint([lsggm.interpolate(x,normalized=True) for x in lsplitpos])
        lshapesplit=shapely.ops.split(lsggm,lsplitter.buffer(1e-8))
        lshapesplit=[lshapesplit[x].wkt for x in range(1,len(lshapesplit),2)]
        lsg['geom']=lshapesplit
        return lsg
    except:
        print(str(lsgbk)+' ERROR')



## Clean PAVEMENTEDGE and MEDIAN (8 mins)
#start=datetime.datetime.now()
#pvmtedge=gpd.read_file(path+'PLANIMETRIC/PAVEMENTEDGE.shp')
#pvmtedge=pvmtedge.to_crs({'init':'epsg:4326'})
#pvmtedge['bkfaceid']=pd.to_numeric(pvmtedge['BLOCKFACEI'])
#pvmtedge['shplen']=pd.to_numeric(pvmtedge['SHAPE_Leng'])
#pvmtedge=pvmtedge[['bkfaceid','shplen','geometry']].reset_index(drop=True)
#mdn=gpd.read_file(path+'PLANIMETRIC/MEDIAN.shp')
#mdn=mdn.to_crs({'init':'epsg:4326'})
#mdn['medians']=1
#mdn=mdn[['medians','geometry']].reset_index(drop=True)
#pvmtedgemdn=gpd.sjoin(pvmtedge,mdn,how='left',op='intersects')
#pvmtedgemdn['medians']=np.where(pd.isna(pvmtedgemdn['medians']),0,pvmtedgemdn['medians'])
#pvmtedgemdn=pvmtedgemdn.groupby(['bkfaceid','shplen'],as_index=False).agg({'medians':max})
#pvmtedgemdn=pd.merge(pvmtedge,pvmtedgemdn,how='left',on=['bkfaceid','shplen'])
#pvmtedgemdn=pvmtedgemdn[[type(x)==shapely.geometry.linestring.LineString for x in pvmtedgemdn['geometry']]]
#pvmtedgemdn['latfrom']=[x.xy[1][0] for x in pvmtedgemdn['geometry']]
#pvmtedgemdn['longfrom']=[x.xy[0][0] for x in pvmtedgemdn['geometry']]
#pvmtedgemdn['latto']=[x.xy[1][-1] for x in pvmtedgemdn['geometry']]
#pvmtedgemdn['longto']=[x.xy[0][-1] for x in pvmtedgemdn['geometry']]
#pvmtedgemdn['pvbearing']=np.nan
#pvmtedgemdn['pvdir']=''
#pvmtedgemdn=pvmtedgemdn.apply(calcpvmtbearing,axis=1)
#pvmtedgemdn=pvmtedgemdn.apply(calcpvmtdir,axis=1)
#pvmtedgemdn=pvmtedgemdn[['bkfaceid','pvbearing','pvdir','shplen','medians','geometry']].reset_index(drop=True)
#pvmtedgemdn.crs={'init':'epsg:4326'}
#pvmtedgemdn.to_file(path+'PLANIMETRIC/PVMTEDGEMDN.shp')
#print(datetime.datetime.now()-start)

## Clean HYDRANT (90 mins)
#start=datetime.datetime.now()
#hydrant=gpd.read_file(path+'HYDRANT/HYDRANT.shp')
#hydrant.crs={'init':'epsg:4326'}
#hydrant=hydrant.to_crs({'init':'epsg:6539'})
#hydrant=hydrant[['unitid','geometry']].reset_index(drop=True)
#hydrant['orggeom']=[x.wkt for x in hydrant['geometry']]
#hydrant['bkfaceid']=np.nan
#hydrant['snapdist']=np.nan
#hydrant['adjgeom']=''
#pvmtedgemdn=gpd.read_file(path+'PLANIMETRIC/PVMTEDGEMDN.shp')
#pvmtedgemdn.crs={'init':'epsg:4326'}
#pvmtedgemdn=pvmtedgemdn.to_crs({'init':'epsg:6539'})
#hydrantbuffer=gpd.GeoDataFrame(hydrant[['unitid']],geometry=hydrant['geometry'].buffer(50),crs={'init':'epsg:6539'}).reset_index(drop=True)
#hydrantbuffer=gpd.sjoin(hydrantbuffer,pvmtedgemdn,how='left',op='intersects')
#hydrantbuffer=hydrantbuffer[['unitid','bkfaceid']].dropna().drop_duplicates().reset_index(drop=True)
#hydrant=hydrant.apply(adjhydrant,axis=1)
#hydrant=hydrant.drop(['geometry'],axis=1).dropna().drop_duplicates().reset_index(drop=True)
#hydrant=gpd.GeoDataFrame(hydrant,geometry=hydrant['adjgeom'].map(wkt.loads),crs={'init':'epsg:6539'})
#hydrant=hydrant.to_crs({'init':'epsg:4326'})
#hydrant.to_file(path+'HYDRANT/HYDRANTADJ.shp')
#print(datetime.datetime.now()-start)

## Clean DOTREG Shapefile (500 mins)
#start=datetime.datetime.now()
#pvmtedgemdn=gpd.read_file(path+'PLANIMETRIC/PVMTEDGEMDN.shp')
#pvmtedgemdn.crs={'init':'epsg:4326'}
#pvmtedgemdn=pvmtedgemdn.to_crs({'init':'epsg:6539'})
#dotshp=gpd.read_file(path+'DOTREG/Parking_Regulation_Shapefile/Parking_Regulation_Shapefile.shp')
#dotshp.crs={'init':'epsg:4326'}
#dotshp=dotshp.to_crs({'init':'epsg:6539'})
#dotshp['boro']=[' '.join(x.split()).upper() for x in dotshp['SG_KEY_BOR']]
#dotshp['boro']=np.where(dotshp['boro']=='M',1,np.where(dotshp['boro']=='B',2,np.where(dotshp['boro']=='K',3,np.where(dotshp['boro']=='Q',4,np.where(dotshp['boro']=='S',5,0)))))
#dotshp['order']=[' '.join(x.split()).upper() for x in dotshp['SG_ORDER_N']]
#dotshp['seq']=pd.to_numeric(dotshp['SG_SEQNO_N'])
#dotshp['arrow']=[' '.join(x.split()).upper() if pd.notna(x) else '' for x in dotshp['SG_ARROW_D']]
#dotshp['desc']=[' '.join(x.split()).upper() for x in dotshp['SIGNDESC1']]
#dotshp['mutcd']=[' '.join(x.split()).upper() for x in dotshp['SG_MUTCD_C']]
#dotshp=dotshp[['boro','order','seq','arrow','desc','mutcd','geometry']].reset_index(drop=True)
#dotshpbuffer=gpd.GeoDataFrame(dotshp[['order','seq']],geometry=dotshp['geometry'].buffer(50),crs={'init':'epsg:6539'}).reset_index(drop=True)
#dotshpbuffer=gpd.sjoin(dotshpbuffer,pvmtedgemdn,how='left',op='intersects')
#dotshpbuffer=dotshpbuffer[['order','seq','bkfaceid']].dropna().drop_duplicates().reset_index(drop=True)
#dotshp=dotshp.apply(adjdotshp,axis=1)
#dotshp=dotshp.drop(['geometry'],axis=1).dropna().drop_duplicates().reset_index(drop=True)
#dotshp=gpd.GeoDataFrame(dotshp,geometry=dotshp['adjgeom'].map(wkt.loads),crs={'init':'epsg:6539'})
#dotshp=dotshp.to_crs({'init':'epsg:4326'})
#dotshp.to_file(path+'DOTREG/Parking_Regulation_Shapefile/DOTSHPADJ.shp')
#print(datetime.datetime.now()-start)

## Join PVMTEDGEMDN to DOTREG Shapefile (8 mins)
#start=datetime.datetime.now()
#pvmtedgemdn=gpd.read_file(path+'PLANIMETRIC/PVMTEDGEMDN.shp')
#pvmtedgemdn.crs={'init':'epsg:4326'}
#dotshpadj=gpd.read_file(path+'DOTREG/Parking_Regulation_Shapefile/DOTSHPADJ.shp')
#dotshpadj.crs={'init':'epsg:4326'}
#dotshpadj['lat']=[x.y for x in dotshpadj['geometry']]
#dotshpadj['long']=[x.x for x in dotshpadj['geometry']]
#dotshpadj=pd.merge(dotshpadj,pvmtedgemdn,how='left',on='bkfaceid')
#dotshpadjclean=dotshpadj[['order','bkfaceid']].drop_duplicates().drop_duplicates(['order'],keep=False)
#dotshpadjclean=dotshpadjclean.drop_duplicates(['bkfaceid'],keep=False).reset_index(drop=True)
#dotshpadjclean=pd.merge(dotshpadj,dotshpadjclean,how='left',on=['order'])
#dotshpadjclean=dotshpadjclean.loc[pd.notna(dotshpadjclean['bkfaceid_y']),['order','seq','lat','long','bkfaceid_y','pvbearing','pvdir','shplen','geometry_y']].reset_index(drop=True)
#dotshpadjclean=dotshpadjclean.rename(columns={'bkfaceid_y':'bkfaceid','geometry_y':'geometry'})
#dotshpadjclean['bearing']=np.nan
#dotshpadjclean['bearingdiff']=np.nan
#dotshpadjclean['geom']=np.nan
#dotshpadjclean=dotshpadjclean.groupby('order',as_index=False).apply(checkdotshpadjgeom).reset_index(drop=True)
#dotshpadjclean.to_csv(path+'DOTREG/DOTSHPADJCLEAN.csv',index=False)
#print(datetime.datetime.now()-start)

## Clean DOTREG SIGN (2.5 mins)
#start=datetime.datetime.now()
#dotregsign=pd.read_csv(path+'DOTREG/SIGN.csv',encoding='latin_1',dtype=str)
#dotregsignclean=dotregsign.copy()
#dotregsignclean['boro']=[' '.join(x.split()).upper() for x in dotregsignclean['SRP_Boro']]
#dotregsignclean['boro']=np.where(dotregsignclean['boro']=='M',1,np.where(dotregsignclean['boro']=='B',2,np.where(dotregsignclean['boro']=='K',3,np.where(dotregsignclean['boro']=='Q',4,np.where(dotregsignclean['boro']=='S',5,0)))))
#dotregsignclean['order']=[' '.join(x.split()).upper() for x in dotregsignclean['SRP_Order']]
#dotregsignclean['seq']=pd.to_numeric(dotregsignclean['SRP_Seq'])
#dotregsignclean['dist']=pd.to_numeric(dotregsignclean['SR_Distx'])
#dotregsignclean['arrow']=[' '.join(x.split()).upper() if pd.notna(x) else '' for x in dotregsignclean['SR_Arrow']]
#dotregsignclean['desc']=''
#dotregsignclean['mutcd']=''
#dotregsignclean=dotregsignclean.apply(sepdescmutcd,axis=1)
#dotregsignclean=dotregsignclean[['boro','order','seq','dist','arrow','desc','mutcd']].drop_duplicates().reset_index(drop=True)
#dotregsignclean.to_csv(path+'DOTREG/SIGNCLEAN.csv',index=False)
#print(datetime.datetime.now()-start)

## Clean DOTREG SIGN MUTCD (8 mins)
#start=datetime.datetime.now()
#dotregsignclean=pd.read_csv(path+'DOTREG/SIGNCLEAN.csv',dtype=str,converters={'boro':float,'seq':float,'dist':float})
#mutcd=dotregsignclean.groupby(['mutcd','desc'],as_index=False).agg({'order':'count'}).rename(columns={'order':'count'})
#mutcd['mutcd']=[' '.join(x.split()).upper() for x in mutcd['mutcd']]
#mutcd['desc']=[' '.join(x.split()).upper() for x in mutcd['desc']]
#mutcd=mutcd.drop_duplicates().sort_values('count',ascending=False).reset_index(drop=True)
#mutcd['descnum']=[len(list(re.finditer('>',x))) for x in mutcd['desc']]
#mutcd['desc1']=''
#mutcd['dir1']=0
#mutcd['desc2']=''
#mutcd['dir2']=0
#mutcd['desc3']=''
#mutcd['dir3']=0
#mutcd['desctail']=''
#mutcd=mutcd.apply(sepdesc,axis=1)
#mutcd1=mutcd[['mutcd','descnum','desc','desc1','dir1']].rename(columns={'desc':'orgdesc','desc1':'sepdesc','dir1':'dir'})
#mutcd1['descseq']=1
#mutcd2=mutcd[['mutcd','descnum','desc','desc2','dir2']].rename(columns={'desc':'orgdesc','desc2':'sepdesc','dir2':'dir'})
#mutcd2['descseq']=2
#mutcd3=mutcd[['mutcd','descnum','desc','desc3','dir3']].rename(columns={'desc':'orgdesc','desc3':'sepdesc','dir3':'dir'})
#mutcd3['descseq']=3
#mutcd=pd.concat([mutcd1,mutcd2,mutcd3],axis=0).dropna().reset_index(drop=True)
#mutcd=mutcd.loc[mutcd['dir']!=0,['mutcd','orgdesc','descnum','descseq','sepdesc','dir']].reset_index(drop=True)
#mutcd=mutcd[~np.isin(mutcd['mutcd'],[x+'A' for x in mutcd['mutcd']])].reset_index(drop=True)
#mutcd['type']=''
#mutcd['time']=''
#mutcd['m']=0
#mutcd['t']=0
#mutcd['w']=0
#mutcd['r']=0
#mutcd['f']=0
#mutcd['s']=0
#mutcd['u']=0
#mutcd['daysflag']=0
#mutcd['dayscheck']=0
#mutcd['starthour']=''
#mutcd['endhour']=''
#mutcd['hoursflag']=0
#mutcd['hourscheck']=0
#mutcd['manual']=0
#for i in mutcd.index:
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('DESCRIPTION NOT AVAILABLE',mutcd.loc[i,'sepdesc'])))>0,'DNA',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='DNA':
#        mutcd.loc[i,['m','t','w','r','f','s','u']]=0
#        mutcd.loc[i,['starthour','endhour']]=''
#        mutcd.loc[i,['daysflag','hoursflag']]+=1
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('HOUR PARKING',mutcd.loc[i,'sepdesc'])))>0,'HOUR PARKING',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='HOUR PARKING':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('HOUR PARKING',mutcd.loc[i,'sepdesc']))[0].start()+12:].split())
#        mutcd.loc[i,['m','t','w','r','f','s','u']]=0
#        mutcd.loc[i,['starthour','endhour']]=''
#        mutcd.loc[i,['daysflag','hoursflag']]+=1
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('HMP',mutcd.loc[i,'sepdesc'])))>0,'HMP',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='HMP':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('HMP',mutcd.loc[i,'sepdesc']))[0].start()+3:].split())
#        mutcd.loc[i,['m','t','w','r','f','s','u']]=0
#        mutcd.loc[i,['starthour','endhour']]=''
#        mutcd.loc[i,['daysflag','hoursflag']]+=1
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('HOUR METERED PARKING',mutcd.loc[i,'sepdesc'])))>0,'HOUR METERED PARKING',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='HOUR METERED PARKING':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('HOUR METERED PARKING',mutcd.loc[i,'sepdesc']))[0].start()+20:].split())
#        mutcd.loc[i,['m','t','w','r','f','s','u']]=0
#        mutcd.loc[i,['starthour','endhour']]=''
#        mutcd.loc[i,['daysflag','hoursflag']]+=1
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('COMMERCIAL VEHICLES',mutcd.loc[i,'sepdesc'])))>0,'COMMERCIAL',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='COMMERCIAL':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('COMMERCIAL VEHICLES',mutcd.loc[i,'sepdesc']))[0].start()+19:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('TRUCK \(SYMBOL\)',mutcd.loc[i,'sepdesc'])))>0,'TRUCK',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='TRUCK':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('TRUCK \(SYMBOL\)',mutcd.loc[i,'sepdesc']))[0].start()+14:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('TAXI \(SYMBOL\)',mutcd.loc[i,'sepdesc'])))>0,'TAXI',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='TAXI':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('TAXI \(SYMBOL\)',mutcd.loc[i,'sepdesc']))[0].start()+13:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('TAXI HAIL \(SYMBOL\)',mutcd.loc[i,'sepdesc'])))>0,'TAXI HAIL',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='TAXI HAIL':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('TAXI HAIL \(SYMBOL\)',mutcd.loc[i,'sepdesc']))[0].start()+18:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('TAXI HAILING \(SYMBOL\)',mutcd.loc[i,'sepdesc'])))>0,'TAXI HAILING',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='TAXI HAILING':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('TAXI HAILING \(SYMBOL\)',mutcd.loc[i,'sepdesc']))[0].start()+21:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('FHV \(SYMBOL\)',mutcd.loc[i,'sepdesc'])))>0,'FHV',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='FHV':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('FHV \(SYMBOL\)',mutcd.loc[i,'sepdesc']))[0].start()+12:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('STAR \(SYMBOL\)',mutcd.loc[i,'sepdesc'])))>0,'STAR',mutcd.loc[i,'type'])    
#    if mutcd.loc[i,'type']=='STAR':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('STAR \(SYMBOL\)',mutcd.loc[i,'sepdesc']))[0].start()+13:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('CROSS \(SYMBOL\)',mutcd.loc[i,'sepdesc'])))>0,'CROSS',mutcd.loc[i,'type'])    
#    if mutcd.loc[i,'type']=='CROSS':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('CROSS \(SYMBOL\)',mutcd.loc[i,'sepdesc']))[0].start()+14:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('BUS \(SYMBOL\)',mutcd.loc[i,'sepdesc'])))>0,'BUS',mutcd.loc[i,'type'])    
#    if mutcd.loc[i,'type']=='BUS':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('BUS \(SYMBOL\)',mutcd.loc[i,'sepdesc']))[0].start()+12:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('NO STOPPING',mutcd.loc[i,'sepdesc'])))>0,'NO STOPPING',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='NO STOPPING':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('NO STOPPING',mutcd.loc[i,'sepdesc']))[0].start()+11:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('NO STANDING',mutcd.loc[i,'sepdesc'])))>0,'NO STANDING',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='NO STANDING':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('NO STANDING',mutcd.loc[i,'sepdesc']))[0].start()+11:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.loc[i,'type']=np.where(len(list(re.finditer('NO PARKING',mutcd.loc[i,'sepdesc'])))>0,'NO PARKING',mutcd.loc[i,'type'])
#    if mutcd.loc[i,'type']=='NO PARKING':
#        mutcd.loc[i,'time']=' '.join(mutcd.loc[i,'sepdesc'][list(re.finditer('NO PARKING',mutcd.loc[i,'sepdesc']))[0].start()+10:].split())
#        mutcd.loc[i]=extractdays(mutcd.loc[i])
#        mutcd.loc[i]=extracthours(mutcd.loc[i])
#    mutcd.to_csv(path+'DOTREG/MUTCD.csv',index=False)
#print(datetime.datetime.now()-start)

# Manually check and categorize MUTCD

## Create No Parking Tags for MUTCDCHECK (2 mins)
#start=datetime.datetime.now()
#mutcdcheck=pd.read_csv(path+'DOTREG/MUTCDCHECK.csv',dtype=float,converters={'mutcd':str,'orgdesc':str,'sepdesc':str,'type':str,
#                                                                            'time':str,'starthour':str,'endhour':str})
#mutcdtag=mutcdcheck.copy()
#for i in [str(x)+str(y) for x in ['m','t','w','r','f','s','u'] for y in [str(x)[11:13]+str(x)[14:16] for x in pd.date_range('00:00','23:30',freq='30min')]]:
#    mutcdtag[i]=0
#for i in mutcdtag.index:
#    if (mutcdtag.loc[i,'dayscheck']==1)&(mutcdtag.loc[i,'hourscheck']==1):
#        if pd.Timestamp(mutcdtag.loc[i,'starthour'])<pd.Timestamp(mutcdtag.loc[i,'endhour']):
#            npdays=[x for x in ['m','t','w','r','f','s','u'] if mutcdtag.loc[i,x]>0]
#            nphours=[str(x)[11:13]+str(x)[14:16] for x in pd.date_range(mutcdtag.loc[i,'starthour'],mutcdtag.loc[i,'endhour'],freq='30min')][:-1]
#            npdayshours=[str(x)+str(y) for x in npdays for y in nphours]
#            mutcdtag.loc[i,npdayshours]=1
#        elif pd.Timestamp(mutcdtag.loc[i,'starthour'])>=pd.Timestamp(mutcdtag.loc[i,'endhour']):
#            npdays1=[x for x in ['m','t','w','r','f','s','u'] if mutcdtag.loc[i,x]>0]
#            nphours1=[str(x)[11:13]+str(x)[14:16] for x in pd.date_range(mutcdtag.loc[i,'starthour'],'23:30',freq='30min')]
#            npdayshours1=[str(x)+str(y) for x in npdays1 for y in nphours1]
#            npdays2=[['m','t','w','r','f','s','u','m'][list(['m','t','w','r','f','s','u']).index(x)+1] for x in ['m','t','w','r','f','s','u'] if mutcdtag.loc[i,x]>0]
#            nphours2=[str(x)[11:13]+str(x)[14:16] for x in pd.date_range('00:00',mutcdtag.loc[i,'endhour'],freq='30min')][:-1]
#            npdayshours2=[str(x)+str(y) for x in npdays2 for y in nphours2]
#            mutcdtag.loc[i,npdayshours1+npdayshours2]=1
#mutcdtag.to_csv(path+'DOTREG/MUTCDTAG.csv',index=False)
#print(datetime.datetime.now()-start)

# Manually check and categorize MUTCDTAG

# Join and summarize DOTREG SIGN and Adjusted DOTREG Shapefile (30 mins)
start=datetime.datetime.now()
dotregsignclean=pd.read_csv(path+'DOTREG/SIGNCLEAN.csv',dtype=str,converters={'boro':float,'seq':float,'dist':float})
dotshpadjclean=pd.read_csv(path+'DOTREG/DOTSHPADJCLEAN.csv',dtype=float,converters={'order':str,'pvdir':str,'geom':str})
dotregsignsum=pd.merge(dotregsignclean,dotshpadjclean,how='left',on='order')
dotregsignsum=dotregsignsum[pd.notna(dotregsignsum['geom'])]
dotregsignsum=dotregsignsum[['order','seq','dist','arrow','mutcd','pvdir']].sort_values(['order','seq'],ascending=True).reset_index(drop=True)
mutcdtagcheck=pd.read_csv(path+'DOTREG/MUTCDTAGCHECK.csv',dtype=float,converters={'mutcd':str,'orgdesc':str,'sepdesc':str,'type':str,
                                                                                  'time':str,'starthour':str,'endhour':str})
mutcdtagcheck=mutcdtagcheck[['mutcd','t0930','t2030','w1600']].reset_index(drop=True)
mutcdtagchecka=mutcdtagcheck.copy()
mutcdtagchecka['mutcd']=[str(x)+'A' for x in mutcdtagchecka['mutcd']]
mutcdtagcheck=pd.concat([mutcdtagcheck,mutcdtagchecka],axis=0,ignore_index=True)
dotregsignsum=pd.merge(dotregsignsum,mutcdtagcheck,how='left',on='mutcd')
dotregsignsum['arrow']=dotregsignsum['arrow'].fillna('')
dotregsignsum['cl']=np.where(dotregsignsum['mutcd']=='CL',1,0)
dotregsignsum['blpl']=np.where(np.isin(dotregsignsum['mutcd'],['BL','PL']),1,0)
dotregsignsum=dotregsignsum.groupby(['order','dist','arrow','pvdir'],as_index=False).agg({'t0930':'sum','t2030':'sum','w1600':'sum','cl':'sum','blpl':'sum'}).reset_index(drop=True)
dotregsignsum=dotregsignsum.groupby('order',as_index=False).apply(sumsign).reset_index(drop=True)
dotregsignsum.to_csv(path+'DOTREG/SIGNSUM2.csv',index=False)
print(datetime.datetime.now()-start)

# Add HYDRANT to SIGNSUM (15 mins)
start=datetime.datetime.now()
hydrantadj=gpd.read_file(path+'HYDRANT/HYDRANTADJ.shp')
hydrantadj.crs={'init':'epsg:4326'}
hydrantadj=hydrantadj.to_crs({'init':'epsg:6539'})
dotshpadjclean=pd.read_csv(path+'DOTREG/DOTSHPADJCLEAN.csv',dtype=float,converters={'order':str,'pvdir':str,'geom':str})
dotshpadjclean=gpd.GeoDataFrame(dotshpadjclean,geometry=dotshpadjclean['geom'].map(wkt.loads),crs={'init':'epsg:4326'})
dotshpadjclean=dotshpadjclean.to_crs({'init':'epsg:6539'})
dotregsignsum=pd.read_csv(path+'DOTREG/SIGNSUM2.csv',dtype=float,converters={'order':str})
dotregsignsumhdt=dotregsignsum.groupby('order',as_index=False).apply(hydrantsign).reset_index(drop=True)
dotregsignsumhdt.to_csv(path+'DOTREG/SIGNSUMHDT2.csv',index=False)
print(datetime.datetime.now()-start)

# Split geometry (10 mins)
start=datetime.datetime.now()
dotshpadjclean=pd.read_csv(path+'DOTREG/DOTSHPADJCLEAN.csv',dtype=float,converters={'order':str,'pvdir':str,'geom':str})
dotshpadjclean=gpd.GeoDataFrame(dotshpadjclean,geometry=dotshpadjclean['geom'].map(wkt.loads),crs={'init':'epsg:4326'})
dotregsignsumhdt=pd.read_csv(path+'DOTREG/SIGNSUMHDT2.csv',dtype=float,converters={'order':str})
dotregsignsumhdtgm=dotregsignsumhdt.groupby('order',as_index=False).apply(splitgm).reset_index(drop=True)
dotregsignsumhdtgm.to_csv(path+'DOTREG/SIGNSUMHDTGM2.csv',index=False)
dotregsignsumhdtgm=pd.read_csv(path+'DOTREG/SIGNSUMHDTGM2.csv',dtype=float,converters={'order':str,'geom':str})
dotregsignsumhdtgm=gpd.GeoDataFrame(dotregsignsumhdtgm,geometry=dotregsignsumhdtgm['geom'].map(wkt.loads),crs={'init':'epsg:4326'})
dotregsignsumhdtgm.to_file(path+'SIGNSUMHDTGM2.shp')
print(datetime.datetime.now()-start)

# Other parking without regulations (25 mins)
start=datetime.datetime.now()
lion=gpd.read_file(path+'LION/LION.shp')
lionpk=lion[['LBlockFace','RBlockFace','FeatureTyp','SegmentTyp','TrafDir','Number_Par']].drop_duplicates().reset_index(drop=True)
lionpk=lionpk.dropna(axis=0)
lionpk['lbkfaceid']=pd.to_numeric(lionpk['LBlockFace'])
lionpk['rbkfaceid']=pd.to_numeric(lionpk['RBlockFace'])
lionpk['featuretype']=[' '.join(x.split()).upper() for x in lionpk['FeatureTyp']]
lionpk['segmenttype']=[' '.join(x.split()).upper() for x in lionpk['SegmentTyp']]
lionpk['trafficdir']=[' '.join(x.split()).upper() for x in lionpk['TrafDir']]
lionpk['parkinglane']=pd.to_numeric(lionpk['Number_Par'])
lionpk=lionpk[np.isin(lionpk['featuretype'],['0','6','A','C'])].reset_index(drop=True)
lionpk=lionpk[np.isin(lionpk['segmenttype'],['B','R','U','S'])].reset_index(drop=True)
lionpk=lionpk[lionpk['trafficdir']!='P'].reset_index(drop=True)
lionpk=lionpk[['lbkfaceid','rbkfaceid','trafficdir','parkinglane']].drop_duplicates().reset_index(drop=True)
lionpk['lparking']=np.nan
lionpk['rparking']=np.nan
lionpk=lionpk.apply(lionpkbkface,axis=1)
lionpkl=lionpk[['lbkfaceid','lparking']].rename(columns={'lbkfaceid':'bkfaceid','lparking':'parking'})
lionpkr=lionpk[['rbkfaceid','rparking']].rename(columns={'rbkfaceid':'bkfaceid','rparking':'parking'})
lionpk=pd.concat([lionpkl,lionpkr],axis=0,ignore_index=True)
lionpk=lionpk.drop_duplicates().drop_duplicates(['bkfaceid'],keep=False).reset_index(drop=True)
lionpk=lionpk[lionpk['parking']>0].reset_index(drop=True) # Still need to solve the 0.5 parking
dotshpadjclean=pd.read_csv(path+'DOTREG/DOTSHPADJCLEAN.csv',dtype=float,converters={'order':str,'pvdir':str,'geom':str})
dotshpadjcleanbk=dotshpadjclean[['bkfaceid','order']].reset_index(drop=True)
lionpk=pd.merge(lionpk,dotshpadjcleanbk,how='left',on='bkfaceid')
lionpk=lionpk.loc[pd.isna(lionpk['order']),['bkfaceid','parking']].reset_index(drop=True)
lionpk['dist1']=np.nan
lionpk['dist2']=np.nan
lionpk['hdt']=0
pvmtedgemdn=gpd.read_file(path+'PLANIMETRIC/PVMTEDGEMDN.shp')
pvmtedgemdn.crs={'init':'epsg:4326'}
pvmtedgemdn=pvmtedgemdn.to_crs({'init':'epsg:6539'})
hydrantadj=gpd.read_file(path+'HYDRANT/HYDRANTADJ.shp')
hydrantadj.crs={'init':'epsg:4326'}
hydrantadj=hydrantadj.to_crs({'init':'epsg:6539'})
lionpk=lionpk.groupby('bkfaceid',as_index=False).apply(lionpkhydrant).reset_index(drop=True)
pvmtedgemdn=pvmtedgemdn.to_crs({'init':'epsg:4326'})
lionpk=lionpk.groupby('bkfaceid',as_index=False).apply(lionpksplitgeom).reset_index(drop=True)
lionpk=lionpk[['bkfaceid','dist1','dist2','parking','hdt','geom']].reset_index(drop=True)
lionpk.to_csv(path+'LIONPK2.csv',index=False)
lionpk=pd.read_csv(path+'LIONPK2.csv',dtype=float,converters={'geom':str})
lionpk=gpd.GeoDataFrame(lionpk,geometry=lionpk['geom'].map(wkt.loads),crs={'init':'epsg:4326'})
lionpk.to_file(path+'LIONPK2.shp')
print(datetime.datetime.now()-start)

# Other non-parking blockface (0.5 min)
start=datetime.datetime.now()
pvmtedgemdn=gpd.read_file(path+'PLANIMETRIC/PVMTEDGEMDN.shp')
pvmtedgemdn.crs={'init':'epsg:4326'}
dotshpadjclean=pd.read_csv(path+'DOTREG/DOTSHPADJCLEAN.csv',dtype=float,converters={'order':str,'pvdir':str,'geom':str})
lionpk=pd.read_csv(path+'LIONPK2.csv',dtype=float,converters={'geom':str})
lionpk=lionpk[['bkfaceid','parking']].drop_duplicates().reset_index(drop=True)
addnopk=pd.merge(pvmtedgemdn,dotshpadjclean,how='left',on='bkfaceid')
addnopk=addnopk.loc[pd.isna(addnopk['order']),['bkfaceid','geometry']].reset_index(drop=True)
addnopk=pd.merge(addnopk,lionpk,how='left',on='bkfaceid')
addnopk=addnopk.loc[pd.isna(addnopk['parking']),['bkfaceid','geometry']].reset_index(drop=True)
addnopk.to_file(path+'ADDNOPK2.shp')
print(datetime.datetime.now()-start)

# Summarize parking spaces (3.5 mins)
start=datetime.datetime.now()
county=gpd.read_file(path+'BOUNDARY/COUNTY.shp')
county.crs={'init':'epsg:4326'}
county=county.to_crs({'init':'epsg:6539'})
county['geoid']=county['GEOID'].copy()
dotregsignsumhdtgm=gpd.read_file(path+'SIGNSUMHDTGM2.shp')
dotregsignsumhdtgm.crs={'init':'epsg:4326'}
dotregsignsumhdtgm=dotregsignsumhdtgm.to_crs({'init':'epsg:6539'})
dotregsignsumhdtgm['signsumspaces']=[x.length/21 for x in dotregsignsumhdtgm['geometry']]
dotregsignsumhdtgm=dotregsignsumhdtgm[(dotregsignsumhdtgm['w1600']==0)&(dotregsignsumhdtgm['hdt']==0)].reset_index(drop=True)
dotregsignsumhdtgm=gpd.sjoin(dotregsignsumhdtgm,county,how='left',op='intersects')
dotregsignsumhdtgm=dotregsignsumhdtgm.groupby('geoid',as_index=False).agg({'signsumspaces':'sum'}).reset_index(drop=True)
lionpk=gpd.read_file(path+'LIONPK2.shp')
lionpk.crs={'init':'epsg:4326'}
lionpk=lionpk.to_crs({'init':'epsg:6539'})
lionpk['lionpkspaces']=[x.length/21 for x in lionpk['geometry']]
lionpk=lionpk[(lionpk['parking']>0)&(lionpk['hdt']==0)].reset_index(drop=True)
lionpk=gpd.sjoin(lionpk,county,how='left',op='intersects')
lionpk=lionpk.groupby('geoid',as_index=False).agg({'lionpkspaces':'sum'}).reset_index(drop=True)
pkspaces=pd.merge(dotregsignsumhdtgm,lionpk,how='left',on='geoid').reset_index(drop=True)
pkspaces['pkspaces']=pkspaces['signsumspaces']+pkspaces['lionpkspaces']
pkspaces=pkspaces[['geoid','pkspaces']].reset_index(drop=True)
pkspaces.to_csv(path+'PKSPACES.csv',index=False)
print(datetime.datetime.now()-start)











        
# Compare with DOT shapefile (15 mins)
start=datetime.datetime.now()
def splitgm(sg):
    try:
        sgod=list(sg['order'])[0]
        sg=sg.reset_index(drop=True)
        sggm=list(dotregloclionpv.loc[dotregloclionpv['order']==sgod,'geometry'])[0]
        splitpos=[x for x in list(sg['dist1'])[2:len(sg)-1]]
        splitter=shapely.geometry.MultiPoint([sggm.interpolate(x,normalized=True) for x in splitpos])
        sg=sg.loc[2:len(sg)-2,['order','dist1']].reset_index(drop=True)
        sg['geom']=[x.wkt for x in splitter]
        return sg
    except:
        print(sgod+' ERROR')
dotregloclionpv=pd.read_csv(path+'DOTREG/LOCLIONPV.csv',dtype=float,converters={'order':str,'pvdir':str,'geom':str})
dotregloclionpv=gpd.GeoDataFrame(dotregloclionpv,geometry=dotregloclionpv['geom'].map(wkt.loads),crs={'init':'epsg:4326'})
dotregsignsum=pd.read_csv(path+'DOTREG/SIGNSUM.csv',dtype=float,converters={'order':str})
dotregsignsumgm=dotregsignsum.groupby('order',as_index=False).apply(splitgm).reset_index(drop=True)
dotregsignsumgm=gpd.GeoDataFrame(dotregsignsumgm,geometry=dotregsignsumgm['geom'].map(wkt.loads),crs={'init':'epsg:4326'})
dotregsignsumgm=dotregsignsumgm.to_file(path+'SIGNSUMGMTEST.shp')
print(datetime.datetime.now()-start)





start=datetime.datetime.now()
print(datetime.datetime.now()-start)
















dotregsign=pd.read_csv(path+'DOTREG/SIGN.csv',encoding='latin_1',dtype=str)
dotregsign['order']=[' '.join(x.split()).upper() for x in dotregsign['SRP_Order']]
dotregsign['seq']=pd.to_numeric(dotregsign['SRP_Seq'])
dotregloc=pd.read_csv(path+'DOTREG/LOCATION.csv',dtype=str)
dotregloc['order']=[' '.join(x.split()).upper() for x in dotregloc['order_no']]
dotregsignsum=pd.read_csv(path+'DOTREG/SIGNSUM.csv',dtype=float,converters={'order':str})
dotregsignsumgm=gpd.read_file(path+'SIGNSUMGMTEST.shp')
dotshp=gpd.read_file(path+'DOTREG/Parking_Regulation_Shapefile/Parking_Regulation_Shapefile.shp')


# Order+Seq
# Total: 761403
# Loc: 663316 / 98087
# Sum: 471036 / 290367
# Shp: 733537 / 27866
# In Shp Not Loc: 94502


k=pd.merge(dotregsign,dotregloc,how='left',on='order')
k=k[pd.isna(k['boro'])]
k=k.groupby(['order','seq'],as_index=False).agg({'SRP_Boro':'count'})


m=pd.merge(dotregsign,dotregsignsum,how='left',on='order')
m=m.drop_duplicates(['order','seq'],keep='first')
m=m[pd.isna(m['dist1'])]
m=m.groupby(['order','seq'],as_index=False).agg({'SRP_Boro':'count'})

l=pd.merge(dotregsign,dotshp,how='left',on='order')
l=l.drop_duplicates(['order','seq_x'],keep='first')
l=l[pd.notna(l['boro'])]
l=l.groupby(['order','seq_x'],as_index=False).agg({'SRP_Boro':'count'})




n=pd.merge(l,k,how='left',left_on=['order','seq_x'],right_on=['order','seq'])
n=n[pd.notna(n['seq'])]




# Order
# Total: 91104
# Loc:  81122 / 9982
# Sum: 55637 / 35467
# Shp: 86174 / 4930
# Shpadj: 59666


j=dotregsign.order.unique()

k=pd.merge(dotregsign,dotregloc,how='left',on='order')
k=k[pd.notna(k['boro'])]
k=k.groupby(['order'],as_index=False).agg({'SRP_Boro':'count'})


m=pd.merge(dotregsign,dotregsignsum,how='left',on='order')
m=m.drop_duplicates(['order','seq'],keep='first')
m=m[pd.isna(m['dist1'])]
m=m.groupby(['order'],as_index=False).agg({'SRP_Boro':'count'})

l=pd.merge(dotregsign,dotshp,how='left',left_on='order',right_on='SG_ORDER_N')
l=l.drop_duplicates(['order'],keep='first')
l=l[pd.isna(l['SG_KEY_BOR'])]
l=l.groupby(['order'],as_index=False).agg({'SRP_Boro':'count'})




n=pd.merge(l,k,how='left',left_on=['order','seq_x'],right_on=['order','seq'])
n=n[pd.notna(n['seq'])]

p=dotshpadj.groupby(['order']).agg({'boro':'count'})


q=pd.merge(dotshp,dotregsign,how='left',left_on='SG_ORDER_N',right_on='order')
q=q.drop_duplicates(['SG_ORDER_N','SG_SEQNO_N'],keep='first')
q=q[pd.isna(l['SRP_Order'])]































## Geosupport (10 mins)
#start=datetime.datetime.now()
#g = Geosupport()
#dotreglocgeo=pd.read_csv(path+'DOTREG/LOCATION.csv',dtype=str)
#dotreglocgeo['order']=[' '.join(x.split()).upper() for x in dotreglocgeo['order_no']]
#dotreglocgeo['boro']=[' '.join(x.split()).upper() for x in dotreglocgeo['boro']]
#dotreglocgeo['boro']=np.where(dotreglocgeo['boro']=='M',1,np.where(dotreglocgeo['boro']=='B',2,np.where(dotreglocgeo['boro']=='K',3,np.where(dotreglocgeo['boro']=='Q',4,np.where(dotreglocgeo['boro']=='S',5,0)))))
#dotreglocgeo['onstreet']=[' '.join(x.split()).upper() for x in dotreglocgeo['main_st']]
#dotreglocgeo['fromstreet']=[' '.join(x.split()).upper() for x in dotreglocgeo['from_st']]
#dotreglocgeo['tostreet']=[' '.join(x.split()).upper() for x in dotreglocgeo['to_st']]
#dotreglocgeo['sos']=[' '.join(x.split()).upper() for x in dotreglocgeo['sos']]
#dotreglocgeo=dotreglocgeo[['order','boro','onstreet','fromstreet','tostreet','sos']].reset_index(drop=True)
#dotreglocgeo['physicalid']=np.nan
#dotreglocgeo['nodeidfrom']=np.nan
#dotreglocgeo['nodeidto']=np.nan
#dotreglocgeo['trafficdir']=''
#dotreglocgeo['parkinglane']=np.nan
#dotreglocgeo['bkfaceid']=np.nan
#for i in dotreglocgeo.index:
#    try:
#        borocode=str(dotreglocgeo.loc[i,'boro'])
#        streetname1=str(dotreglocgeo.loc[i,'onstreet'])
#        streetname2=str(dotreglocgeo.loc[i,'fromstreet'])
#        streetname3=str(dotreglocgeo.loc[i,'tostreet'])
#        compassdir=str(dotreglocgeo.loc[i,'sos'])
#        result=g['3C']({'Borough Code-1':borocode,'Street Name-1':streetname1,'Street Name-2':streetname2,'Street Name-3':streetname3,'Compass Direction':compassdir},mode='extended')
#        dotreglocgeo.loc[i,'physicalid']=pd.to_numeric(result['Physical ID'])
#        dotreglocgeo.loc[i,'nodeidfrom']=pd.to_numeric(result['From Node'])
#        dotreglocgeo.loc[i,'nodeidto']=pd.to_numeric(result['To Node'])
#        dotreglocgeo.loc[i,'trafficdir']=' '.join(result['Traffic Direction'].split()).upper()
#        dotreglocgeo.loc[i,'parkinglane']=pd.to_numeric(result['Number of Parking Lanes on the Street'])
#        dotreglocgeo.loc[i,'bkfaceid']=pd.to_numeric(result['Blockface ID'])
#    except:
#        print(str(i))
#dotreglocgeo.to_csv(path+'DOTREG/LOCGEO.csv',index=False)
#print(datetime.datetime.now()-start)

