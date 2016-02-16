from tethys_sdk.gizmos import DatePicker, MapView, MVLayer, MVView, TextInput, Button, ButtonGroup, LinePlot, ScatterPlot
from django.shortcuts import render
from django.contrib.auth.decorators import login_required
from .model import SessionMaker, StreamGage, TimeSeries, getRecessions

@login_required()
def home(request):
    """
    Controller for the app home page.
    """
    
    context = {}
    
    return render(request, 'my_first_app/home.html', context)
    
def setup(request):
    """
    Controller for the setup page.
    """

    gages  = TextInput(name='gages', display_text='Gage', initial='11477000')      

    start = DatePicker(name='start',
                                            display_text='Start date',
                                            autoclose=True,
                                            format='yyyy-m-d',
                                            start_date='01/01/1910',
                                            initial='2000-01-01')
    stop = DatePicker(name='stop',
                                            display_text='Stop date',
                                            autoclose=True,
                                            format='yyyy-m-d',
                                            start_date='01/01/1910',
                                            initial='2015-01-01')  
                                            
                                                                           

    run_button = Button(display_text='Run',
                        icon='glyphicon glyphicon-play',
                        style='success',
                        submit=True)
                        
    delete_button = Button(display_text='Delete',
                           icon='glyphicon glyphicon-trash',
                           disabled=True,
                           style='danger')
    
                           
    horizontal_buttons = ButtonGroup(buttons=[run_button, delete_button])


    line_plot_view = None
    scatter_plot_view = None
    if request.POST and 'gages' in request.POST:
        gageName = request.POST['gages'].split(',')
        start = request.POST['start']
        stop = request.POST['stop']
        
        ts = TimeSeries(gageName[0],start,stop)
        rec = getRecessions(gageName,ts)
        
        line_plot_view = LinePlot(
        height='500px',
        width='500px',
        engine='highcharts',
        title='Flow Time Series',
        spline=True,
        x_axis_title='Time',
        y_axis_title='Flow',
        y_axis_units='cfs',
        xAxis={
            'type': 'datetime',
            },
        
        series=[{
               'name': gageName,
               'color': '#0066ff',
               'marker': {'enabled': False},
               'data': zip(ts.time,ts.discharge),
               'dateTimeLabelFormats':{'second':'%Y'},
               }]
        )
        
        scatter_plot_view = ScatterPlot(
        height='500px',
        width='500px',
        engine='highcharts',
        title='Recession Parameters',
        spline=True,
        x_axis_title='log(a)',
        y_axis_title='b',
        x_axis_units = '[]',
        y_axis_units = '[]',
        xAxis = {'type':'logarithmic'},
        series=[{
               'name': gageName,
               'color': '#0066ff',
               'data': zip(rec.A,rec.B),
               'dateTimeLabelFormats':{'second':'%Y'},
               }]
        )
                
        
        
        

    context = {'start': start, 'stop':stop, 'gages': gages, 'buttons': horizontal_buttons, 'line_plot_view':line_plot_view, 'scatter_plot_view':scatter_plot_view}

    return render(request, 'my_first_app/setup.html', context)
    
 
    
def results(request):
    session = SessionMaker()
    gages = session.query(StreamGage).all()
    
    features = []
    for gage in gages:
        gage_feature = {
          'type': 'Feature',
          'geometry': {
            'type': 'Point',
            'coordinates': [gage.longitude, gage.latitude]
          }
        }

        features.append(gage_feature)

    geojson_gages = {
      'type': 'FeatureCollection',
      'crs': {
        'type': 'name',
        'properties': {
          'name': 'EPSG:4326'
        }
      },
      'features': features
    }

    # Define layer for Map View
    geojson_layer = MVLayer(source='GeoJSON',
                            options=geojson_gages,
                            legend_title='HUC 18 stream Gages',
                            legend_extent=[-111.74, 40.22, -111.67, 40.25])    
    view_options = MVView(
        projection='EPSG:4326',
        center=[-100, 40],
        zoom=3.5,
        maxZoom=18,
        minZoom=2)

    wms_layer = MVLayer(source='ImageWMS',
                    options={'url': 'http://mrdata.usgs.gov/services/ca?version=1.1.1&amp;service=WMS',
                                'params': {'layers': 'lith-low,lith-high,faults-low,faults-high','style':'default'},
                                'serverType': 'geoserver'},
                    legend_title='USGS Lithology (WMS)')                             
                                            
    map_view_options = MapView(height='600px',width='100%',layers=[geojson_layer],legend=True)
    
    context = {'map_options': map_view_options}
    
    return render(request, 'my_first_app/results.html', context)