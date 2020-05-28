import warnings
import pandas as pd
from pandas_datareader import data
import cufflinks as cf
import ipywidgets as wd
from ipywidgets import interact, interact_manual
from plotly.offline import iplot
from IPython.display import display
from datetime import datetime
from function_box import function_box

warnings.simplefilter(action='ignore', category=FutureWarning)

def ta_dashboard(asset, indicator, start_date, end_date, 
                 bb_k, bb_n, macd_fast, macd_slow, macd_signal,
                 rsi_periods, rsi_upper, rsi_lower,
                 ema_periods,
                 sma_periods,
                 cci_periods, cci_upper, cci_lower):
    
    df = data.DataReader(asset,
                    'yahoo',  
                    start=start_date, 
                    end=end_date)

    qf = cf.QuantFig(df, title=f'TA Dashboard - {asset}', 
                     legend='right', name=f'{asset}')
    qf.add_volume()
            
    if 'Bollinger Bands' in indicator: 
        qf.add_bollinger_bands(periods=bb_n, 
                               boll_std=bb_k)
    if 'MACD' in indicator: 
        qf.add_macd(fast_period=macd_fast, 
                    slow_period=macd_slow, 
                    signal_period=macd_signal)
    if 'RSI' in indicator: 
        qf.add_rsi(periods=rsi_periods, 
                   rsi_upper=rsi_upper, 
                   rsi_lower=rsi_lower, 
                   showbands=True)
    if 'EMA' in indicator:
        qf.add_ema(ema_periods)
    if 'SMA' in indicator:
        qf.add_sma(sma_periods)
    if 'CCI' in indicator:
        qf.add_cci(periods = cci_periods,
                   cci_upper = cci_upper, 
                   cci_lower = cci_lower)

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

    main_selector_box = wd.VBox(children=[main_selector_label,
                                        stocks_selector,
                                        indicator_selector,
                                        start_date_selector,
                                        end_date_selector])

    bb_box_obj = function_box('Bollinger Bands', N = [20, 1, 40, 1],
                                                 k = [2, 0, 40, 1])

    macd_box_obj = function_box('MACD', Fast_Periods = [5, 5, 240, 5],
                                        Slow_Periods = [5, 5, 240, 5],
                                        Signal_Periods = [5, 5, 240, 5])

    rsi_box_obj = function_box('RSI', RSI_Periods = [5, 5, 240, 5],
                                      RSI_UPPER = [70, 1, 100, 1],
                                      RSI_LOWER = [30, 1, 100, 1])

    ema_box_obj = function_box('EMA', EMA_periods = [20, 2, 90, 1])

    sma_box_obj = function_box('SMA', SMA_periods = [20, 2, 90, 1])

    cci_box_obj = function_box('CCI', CCI_periods = [20, 2, 90, 1],
                                      CCI_UPPER = [100, 0, 500, 50],
                                      CCI_LOWER = [-100, -500, 0, -50])

    sec_selector_label = wd.Label('Secondary parameters', 
                                layout=wd.Layout(height='45px'))
    blank_label = wd.Label('', layout=wd.Layout(height='45px'))

    sec_box_1 = wd.VBox([sec_selector_label, bb_box_obj.get_instance(), macd_box_obj.get_instance(), cci_box_obj.get_instance()])
    sec_box_2 = wd.VBox([blank_label, rsi_box_obj.get_instance(), ema_box_obj.get_instance(), sma_box_obj.get_instance()])

    secondary_selector_box = wd.HBox([sec_box_1, sec_box_2])

    controls_dict = {'asset':stocks_selector, 
                    'indicator':indicator_selector, 
                    'start_date':start_date_selector, 
                    'end_date':end_date_selector, 
                    'bb_k':bb_box_obj.get_part()[2], 
                    'bb_n':bb_box_obj.get_part()[1],
                    'macd_fast': macd_box_obj.get_part()[1], 
                    'macd_slow': macd_box_obj.get_part()[2], 
                    'macd_signal': macd_box_obj.get_part()[3],
                    'rsi_periods': rsi_box_obj.get_part()[1], 
                    'rsi_upper': rsi_box_obj.get_part()[2],
                    'rsi_lower': rsi_box_obj.get_part()[3],
                    'ema_periods': ema_box_obj.get_part()[1],
                    'sma_periods': sma_box_obj.get_part()[1],
                    'cci_periods': cci_box_obj.get_part()[1],
                    'cci_upper': cci_box_obj.get_part()[2],
                    'cci_lower': cci_box_obj.get_part()[3]}

    ui = wd.HBox([main_selector_box, secondary_selector_box])
    out = wd.interactive_output(ta_dashboard, controls_dict)
    return ui, out
