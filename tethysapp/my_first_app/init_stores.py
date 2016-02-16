from .model import engine, SessionMaker, Base, StreamGage
import pandas as pd

def init_stream_gage_db(first_time):
    
    """
    An example persistent store initializer function
    """
    # Create tables
    Base.metadata.create_all(engine)

    # Initial data
    if True:
        # Make session
        
        session = SessionMaker()

        sites = ['11477000']
        
        sites = pd.read_csv('/usr/lib/tethys/src/tethys_apps/tethysapp/my_first_app/public/huc_18.tsv',sep='\t',header=30,index_col=False,skiprows=[31])
        sites = sites[sites['site_tp_cd'] == 'ST']
        
        
        for idx, site in sites.iterrows():
            try:

        #        data = ulmo.usgs.nwis.get_site_data(site)
        #        self.no     = id
        #        self.name   = data['00010:00001']['site']['name']
        #        lat    = float(data['00010:00001']['site']['location']['latitude'])
        #        lon    = float(data['00010:00001']['site']['location']['longitude'])
                lat = site.dec_lat_va
                lon = site.dec_long_va   
                gage = StreamGage(latitude=lat,longitude=lon)
                session.add(gage)
                print idx
            except: pass

        session.commit()