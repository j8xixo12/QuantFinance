import warnings
import pandas as pd
from pandas_datareader import data
import cufflinks as cf
import ipywidgets as wd
from ipywidgets import interact, interact_manual
from plotly.offline import iplot
from IPython.display import display
from datetime import datetime, timedelta
from function_box import function_box
import yfinance as yf
yf.pdr_override() # <== that's all it takes :-)

warnings.simplefilter(action='ignore', category=FutureWarning)

def ta_dashboard(asset, indicator, start_date, end_date, 
                 k, N, MACD_Fast_Periods, MACD_Slow_Periods, MACD_Signal_Periods,
                 RSI_Periods, RSI_UPPER, RSI_LOWER,
                 EMA_periods,
                 SMA_periods,
                 CCI_periods, CCI_UPPER, CCI_LOWER,
                 minute_data):
    if minute_data is True:
        df = data.DataReader(asset,
                        interval = '1m', 
                        start = end_date - timedelta(days = 1), 
                        end = end_date)
    else:
        df = data.DataReader(asset, 
                        start = start_date,   
                        end = end_date)

    qf = cf.QuantFig(df, title=f'TA Dashboard - {asset}', 
                     legend='right', name=f'{asset}')
    qf.add_volume()
            
    if 'Bollinger Bands' in indicator: 
        qf.add_bollinger_bands(periods=N, 
                               boll_std=k)
    if 'MACD' in indicator: 
        qf.add_macd(fast_period=MACD_Fast_Periods, 
                    slow_period=MACD_Slow_Periods, 
                    signal_period=MACD_Signal_Periods)
    if 'RSI' in indicator: 
        qf.add_rsi(periods=RSI_Periods, 
                   rsi_upper=RSI_UPPER, 
                   rsi_lower=RSI_LOWER, 
                   showbands=True)
    if 'EMA' in indicator:
        qf.add_ema(EMA_periods)
    if 'SMA' in indicator:
        qf.add_sma(SMA_periods)
    if 'CCI' in indicator:
        qf.add_cci(periods = CCI_periods,
                   cci_upper = CCI_UPPER, 
                   cci_lower = CCI_LOWER)

    fig = qf.figure()
    iplot(fig)    


def dashboard_init(stocks, indicators):

    stocks_selector = wd.Dropdown(
        options=stocks, 
        value=stocks[0], 
        description= 'Stock'
    )

    indicator_selector = wd.SelectMultiple(
        description = 'Indicator',
        options = indicators, 
        value = [indicators[0]]
    )

    start_date_selector = wd.DatePicker(
        description = 'Start Date', 
        value = pd.to_datetime('2018-01-01'), 
        continuous_update = False
    )

    end_date_selector = wd.DatePicker(
        description='End Date', 
        value=datetime.today(), 
        continuous_update = False
    )

    main_selector_label = wd.Label('Main parameters', 
                                layout=wd.Layout(height='45px'))
    
    checkbox = wd.Checkbox(value=False,
                           description='Minute data',
                           disabled=False,
                           indent=False)

    main_selector_box = wd.VBox(children=[main_selector_label,
                                        checkbox,
                                        stocks_selector,
                                        indicator_selector,
                                        start_date_selector,
                                        end_date_selector])

    controls_dict = {'asset':stocks_selector, 
                    'indicator':indicator_selector, 
                    'start_date':start_date_selector, 
                    'end_date':end_date_selector, 
                    'minute_data': checkbox}

    bb_box_obj = function_box('Bollinger Bands', controls_dict, N = [20, 1, 40, 1],
                                                 k = [2, 0, 40, 1])

    macd_box_obj = function_box('MACD', controls_dict, MACD_Fast_Periods = [5, 5, 240, 5],
                                        MACD_Slow_Periods = [5, 5, 240, 5],
                                        MACD_Signal_Periods = [5, 5, 240, 5])

    rsi_box_obj = function_box('RSI', controls_dict, RSI_Periods = [5, 5, 240, 5],
                                      RSI_UPPER = [70, 1, 100, 1],
                                      RSI_LOWER = [30, 1, 100, 1])

    ema_box_obj = function_box('EMA', controls_dict, EMA_periods = [20, 2, 90, 1])

    sma_box_obj = function_box('SMA', controls_dict, SMA_periods = [20, 2, 90, 1])

    cci_box_obj = function_box('CCI', controls_dict, CCI_periods = [20, 2, 90, 1],
                                      CCI_UPPER = [100, 0, 500, 50],
                                      CCI_LOWER = [-100, -500, 0, -50])


    sec_selector_label = wd.Label('Secondary parameters', 
                                layout=wd.Layout(height='45px'))
    blank_label = wd.Label('', layout=wd.Layout(height='45px'))

    sec_box_1 = wd.VBox([sec_selector_label, bb_box_obj.get_instance(), macd_box_obj.get_instance(), cci_box_obj.get_instance()])
    sec_box_2 = wd.VBox([blank_label, rsi_box_obj.get_instance(), ema_box_obj.get_instance(), sma_box_obj.get_instance()])

    secondary_selector_box = wd.HBox([sec_box_1, sec_box_2])

    ui = wd.HBox([main_selector_box, secondary_selector_box])
    out = wd.interactive_output(ta_dashboard, controls_dict)
    return ui, out
