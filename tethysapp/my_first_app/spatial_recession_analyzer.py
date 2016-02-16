import numpy as np  # numerical methods
import matplotlib.pyplot as plt     # plotting stuff
from pandas import *    # analytics package
import os.path  # operating system interface    
import urllib   # download shit 
import seaborn as sns   # overload plotting defaults so the colors are nice
#import fiona    # import/export various file formats
#from shapely.ops import cascaded_union  # find the union of a collection of shapes
#from shapely.geometry import LineString, MultiLineString, Polygon, MultiPolygon, shape # internal line and polygon representation
#from descartes import PolygonPatch  
#from matplotlib.collections import LineCollection, PatchCollection # make lines and polygons into a plot-able format
import pickle # store objects for later use


def main(gauges=['11477000'],startYear='1900',endYear='2015',doAnalysis=True,doPlotAB=False,doPlotTimeSeries=False,doPlotSpatial=False):
    """
    INPUTS:   gauges to be analyzed [Python list of strings]
              startYear [string]
              endYear [string]
              doAnalysis [optional Boolean, default False]
              doPlotXYZ [optional Boolean, default False]
    
    OUTPUTS:  None
    
    Function will check local data folder to see if the data for the desired years have been downloaded before. If yes, proceed. If not, download from USGS website.
    
    If doAnalysis is True, the function will fit an (a,b) to each dry season recession and then decorrelate the point cloud using Bergner-Zouhar.
    
    If doPlotXYZ is True, the function will make a pretty XYZ plot.
    """
    if doAnalysis:
         rec = DataFrame()  # this is an empty pandas data frame object
         for gauge in gauges:
            print('Processing gauge no. ' + gauge)
            try:
                timeSeries = getTimeSeries(gauge,startYear,endYear)  # columns = ['Agency','Site','Time','Discharge','DischargeQualification']
                print('\tGetting time series')
            except:
                print('\tNo daily streamflow from given time period')
                continue
            
            print('\tFitting recessions') 
            new = getRecessions(gauge,timeSeries) # columns = ['Gauge','StartIdx','EndIdx','A','B']

            # in general, this is where quality checks on the recessions go
            if len(new.index) < 10: 
                print('\tToo few quality recessions')
                continue
            
            print('\tResidual correlation between log(a) and b: ' + str(np.corrcoef([np.log(new.A),new.B])[0][1]))

            rec = concat([rec, new])    # append the recessions from the most recently considered catchment to the whole collection

         gauges = rec['Gauge'].drop_duplicates().values.tolist()    # this is just the gauges that actually produced quality recessions
            
         if len(rec.index) == 0: return # if no gauges turned up good recessions, return

         print(rec)
         if doPlotAB: prettyABPlot(rec)
         if doPlotTimeSeries: prettyTimeSeriesPlot(gauge,timeSeries,zip(rec['StartIdx'],rec['EndIdx']))
                 
    if doPlotSpatial: prettySpatialPlot(gauges)

def getTimeSeries(gauge,startYear,endYear):
    """
        INPUTS: gauge [string]
                startYear [string]
                endYear [string]
    
        OUTPUTS: Pandas dataframe with columns ['Agency','Site','Time','Discharge','DischargeQualification']
    
        If the data for the input gauge number between the start and end years is already stored locally, load it.
        Othewise, download it from USGS.
        Note: We only keep accepted (not provisional data).
        Note: There are some borderline useless columsn (e.g., Agency, DischargeQualification).
    
    """
    dateparse = lambda x: datetime.strptime(x, '%Y-%m-%d') # custom function that will be used to make our x axis ticks work right

    filename = 'data/' + gauge  + '_' + startYear + '_' + endYear + '.tsv'; # this is just a convention
    # if the file exists, read it. if not, download it from USGS
    if os.path.isfile(filename):
        df = read_csv(filename,sep='\t',header=26,index_col=False,parse_dates=[2],date_parser=dateparse,skiprows=[27])
    else:
        response = urllib.request.urlopen('http://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no=' + gauge + '&referred_module=sw&period=&begin_date=' + startYear + '-01-01&end_date=' + endYear + '-12-31')
        tsv = response.read().decode("utf8")
        f = open(filename,'w')
        f.write(tsv)
        df =read_csv(filename,sep='\t',header=26,index_col=False,parse_dates=[2],date_parser=dateparse,skiprows=[27])
    
    df.columns = ['Agency','Site','Time','Discharge','DischargeQualification']
    
    df = df[df.DischargeQualification == 'A']   # only considered accepted (not provisional) data points
    df.Discharge = df.Discharge.astype(float)   # make sure discharge is a float, not a string, since we'll have to operate on it numerically
    return df
        
def getRecessions(gauge,timeSeries,minRecessionLength = 10):
    """
        INPUTS: gauge [string]
                timeSeries [numpy array]
                minRecessionLength [optional integer, default = 10]
    
        OUTPUT: pandas data frame with columns [Gauge, StarIdx, EndIdx, A, B]
    
        Note: The definition of 'recessing' could use some work -- easily replaced.
        Note: The method for fitting the recession to an (a,b) pair is general and easily replaced.
    """
    q = timeSeries.Discharge    

    dq = np.append(0,np.diff(q))    # first derivative (sort of -- dt component...)
    ddq = np.append(np.diff(dq),0)  # second derivative (again, sort of...)
    
    timeSeries.IsReceding = np.logical_and(dq < 0, ddq > 0) # simple definition of recession: dec. and concave down
    
    recessions = []
    start = 0
    A = np.array([])
    B = np.array([])
    for i in range(len(timeSeries)):
        if not timeSeries.IsReceding[i]:
            end = i # if we're not receding, set the end to now
            if end - start > minRecessionLength: 
                recessions += [(start,end)]     # if the recession is good, append it (Python note: plus means 'append' here)
                (a,b) = fitRecession(timeSeries.Time[start:end],timeSeries.Discharge[start:end])
                A = np.append(A,a)
                B = np.append(B,b)    
            start = i+1 # assume that the next time step will start a recession
    
    rec = DataFrame()  
    if len(recessions) == 0: return rec  
    
    rec['StartIdx']   = [start for (start,end) in recessions]
    rec['EndIdx']     = [end for (start,end) in recessions]

    # decorrelate the (a,b) point cloud -- this only affects the a values
    A = BergnerZouhar(A,B)    
    rec['A'] = A
    rec['B'] = B
    
    rec['Gauge'] = gauge
    
    return rec
    
def BergnerZouhar(A,B):
    """ INPUTS: A, B [numpy arrays
        OUTPUTS: a [numpy array]
    
        a is the scaled version of A such that log(a) and B are linearly uncorrelated
        Note: I've numerically tested that the correlation is actually removed.
    """
    num = -np.sum((B - np.mean(B))*(np.log(A) - np.mean(np.log(A))))
    den = np.sum((B - np.mean(B))**2)
    q0 = np.exp(num/den)
    return A*q0**(B-1)

def fitRecession(time,discharge):
    """ INPUTS: time [numpy array]
                discharge [numpy array]
    
        OUTPUTS: recession parameter (a,b) pair
    
        This is a very simple fitting procedure, but it can easily be superceded.
    """
    dq = np.diff(discharge)
    p = np.polyfit(np.log(discharge[1:]),np.log(-dq),1)
    return (np.exp(p[1]),p[0]) #(a,b)

def prettyABPlot(rec):
    """
        INPUTS: rec [pandas Dataframe with columns 'A', 'B', and 'Gauge']
    
        Plot the (a,b) point cloud represented by the Pandas data frame rec. Color is based on the gauge number (in order).
    
    """
    gauges = rec['Gauge'].drop_duplicates().values.tolist()
    sns.set_palette(getCommonColorScheme(len(gauges)))
    sns.lmplot('A','B',data=rec,hue='Gauge',fit_reg=False,legend=False)
    plt.gca().set_xscale('log')
    plt.gca().set_xlim([10**-4,10**2])
    plt.xlabel('log(a) [ ]')
    plt.ylabel('b [ ]')
    plt.legend(loc=1)
    
    plt.savefig('a_b.png')  # don't show it, just save it

def prettyTimeSeriesPlot(gauge,df,recessions=[]):
    """
        INPUTS: gauge [string]
                df [pandas data from with columns Time and Discharge]
                rec [optional list of tuples of the form (startIdx, endIdx), default empty]
    
        Plot the time series associated with the input data frame.
        If the optional argument recessions is provided, then each recession is plotted in a different color.
        Note: df should contain only a single gauge's worth of discharge information.
    
    """
    plt.plot(df.Time,df.Discharge,'k')  # plot the whole time series
    for (start,end) in recessions: plt.plot(df.Time[start:end],df.Discharge[start:end],linewidth=2)

    plt.gcf().autofmt_xdate()   # make the xtick's work with dates
    plt.title('USGS gauge No.: ' + gauge) 
    plt.ylabel('Discharge [cfs]')
    plt.savefig(gauge + '.png')  # don't show it, just save it
        
def prettySpatialPlot(gauges):
    """
        INPUTS: gauges [Python list of strings]
    
        OUTPUTS: none
    
        Make a pretty spatial plot showing where the gauges in the input are located.
    
        Right now this plot features a couple things:
        - streams in the HUC in which any gague is located
        - lithology of the region under consideration
    
        We could easily substitute lithology for land use or take it out completely.
    """
    # get the figure ready
    plt.figure()
    ax = plt.gca()

    # read in all the HUCs starting with '18' [California?]
    sites = read_csv('data/huc_18.tsv',sep='\t',header=30,index_col=False,skiprows=[31])

    # plot the state first
    patches = []    
    cali = MultiPolygon([shape(pol['geometry']) for pol in fiona.open('data/states_21basic/states.shp') if pol['properties']['STATE_ABBR'] == 'CA'])
    for p in cali:
        patches.append(PolygonPatch(p, alpha=1, fc='#FFFFFF', ec='#000000', zorder=1))
    ax.add_collection(PatchCollection(patches, match_original=True))

    # now get ready to plot the lithology
    litho = [shape(pol['geometry']) for pol in fiona.open('data/CAgeol_dd/cageol_poly_dd.shp')]
    rockTypes = [pol['properties']['ROCKTYPE1'] for pol in fiona.open('data/CAgeol_dd/cageol_poly_dd.shp')]     # the file contains a couple rock types -- just take the first

    # map a unique number to each rock type; we'll eventually use to assign colors
    colorDict = {}
    for rockType in rockTypes:
        if rockType == 'water': continue
        if rockType not in colorDict: colorDict[rockType] = len(colorDict)

    ## get list of all hydrologic units in which any gauge is located
    hucList = list(sites[sites.site_no.isin(gauges)].huc_cd)

    ## watersheds in the hydro unit
    mp = MultiPolygon([shape(pol['geometry']) for pol in fiona.open('data/CA_HUC12/CA_HUC12.shp') if int(pol['properties']['HUC12'][:8]) in hucList])
    whole = cascaded_union(mp)

    # add each of the watersheds to a plottable data structure, and then plot
    patches = []
    for idx, p in enumerate(mp):
        patches.append(PolygonPatch(p, fc='#FFFFFF',ec='#555555', alpha=0.0, zorder=1))
    ax.add_collection(PatchCollection(patches, match_original=True))

    # set plot bounds
    minx, miny, maxx, maxy = mp.bounds
    w, h = maxx - minx, maxy - miny
    ax.set_xlim(minx - 0.1 * w, maxx + 0.1 * w)
    ax.set_ylim(miny - 0.1 * h, maxy + 0.1 * h)
    ax.set_aspect(1)

    # if we already have a file containing just califronia streams, load it.
    # otherwise, look through the full stream data frame and find the ones in california. save the work so we don't have to do so again.
    filename = 'data/ca_streams.pickle'
    if os.path.isfile(filename):
        streams = pickle.load(open(filename,'rb'))
        streams = MultiLineString(streams)
        streams = streams.intersection(whole)
    else:
        streams = MultiLineString([shape(pol['geometry']) for pol in fiona.open('data/streaml010g.shp_nt00885/streaml010g.shp') if pol['properties']['State'] == 'CA'])
        streams = [LineString([stream.coords[i][:2] for i in range(len(stream.coords))]) for stream in streams]
        f = open(filename,'rb')
        f.write(streams)
        f.close()

    # make a simple Python list from streams (a MultiLineString object)
    c = []
    for stream in streams: c.append(stream)

    # convert the streams into a plottable object and plot it
    c = LineCollection(c,colors=(sns.xkcd_rgb["windows blue"]),linewidths=1,alpha=0.75)
    ax.add_collection(c)        

    
    patches = []
    cm = sns.color_palette("hls", len(colorDict))
    cm = sns.color_palette("PuOr",len(colorDict),desat=0.75)
    for i,p in enumerate(litho):
        if rockTypes[i] == 'water': continue
        patches.append(PolygonPatch(p, alpha=1, fc=cm[colorDict[rockTypes[i]]], ec=cm[colorDict[rockTypes[i]]], zorder=1))        
    ax.add_collection(PatchCollection(patches, match_original=True))       

    # plot all the gauges using the colors from a common color scheme across all plot functions
    cm = getCommonColorScheme(len(gauges))
    for i,gauge in enumerate(gauges):
        lat = sites[sites.site_no == gauge].dec_lat_va
        lon = sites[sites.site_no == gauge].dec_long_va
        plt.plot(lon,lat,'o',color=cm[i%len(cm)],markersize=10)
        plt.text(lon,lat+.01,gauge,fontsize=10)

    # get rid of the axis ticks    
    ax.set_xticks([])
    ax.set_yticks([])
    
    # save the figure
    plt.savefig('spatial.png')

def getCommonColorScheme(n):
    return sns.color_palette("husl", n)

def allGaugesInHUC(hucs):
    """
        INPUT: hucs [Python list strings]
        OUTPUT: gauges [Python list of strings]
    
        Return a list of all stream gauges that are located within any HUC in the input.
    
        (Note that right now this only works for HUCs starting with '18' [Calfornia?])
    
    """
    sites = read_csv('data/huc_18.tsv',sep='\t',header=30,index_col=False,skiprows=[31]) # read the file
    sites = sites[sites.site_tp_cd == 'ST'] # only want streams
    
    gauges = sites[sites.huc_cd.isin(hucs)].sort_values(by='huc_cd').site_no # sort the gauage IDs so that all gagues always appear in the same order in any plot
    gauges = [str(gauge) for gauge in gauges] # type conversion to play well with other functions
    return gauges    

if __name__ == '__main__':
    bigSur = '11143000'
    scotia = '11477000'
    arroyo = '11098000'
    gauges = [bigSur, scotia, arroyo]
    
    gauges = allGaugesInHUC([18010102])
    #gauges = allGaugesInHUC([18010102,18010105])    

    main(gauges,doAnalysis=True,doPlotAB=True)
