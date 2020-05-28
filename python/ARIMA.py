import pandas as pd
import numpy as np
import itertools    
import warnings
from statsmodels.tsa.statespace.sarimax import SARIMAX
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from pandas_datareader import data
import matplotlib.pyplot as plt

def load_financial_data(underlying_stock, start_date, end_date):
    print('Loading financial data...')
    df = data.DataReader(underlying_stock, 'yahoo', start_date, end_date)
    return df

def arima_grid_search(dataframe, s):
    p = d = q = range(2)
    param_combinations = list(itertools.product(p, d, q))

    lowest_aic, pdq, pdqs = None, None, None
    total_iterations = 0
    for order in param_combinations:    
        for (p, q, d) in param_combinations:
            seasonal_order = (p, q, d, s)
            total_iterations += 1
            model = SARIMAX(data, order=order, 
                seasonal_order=seasonal_order, 
                enforce_stationarity=False,
                enforce_invertibility=False,
                disp=False,
                simple_differencing=False
            )
            model_result = model.fit(maxiter=500, disp=False)
            if not lowest_aic or model_result.aic < lowest_aic:
                lowest_aic = model_result.aic
                ret = model_result
    return ret

symbol = '2330.TW'

stock = load_financial_data(underlying_stock = symbol, start_date = '2018-01-01',
                                end_date = '2019-01-20')

data = stock['Close']
acf_data = plot_acf(data)
plt.show()
pacf_data = plot_pacf(data)
plt.show()
plt.plot(data)
plt.show()
data.index = pd.DatetimeIndex(data.index).to_period('D')
# print(data.index)

# result = arima_grid_search(data, 7)
# print(result.summary())
# # result.plot_diagnostics(figsize=(12, 8))

# n = len(data.index)
# print(n)
# prediction = result.get_prediction(
#     start=n - 7 * 5, 
#     end=n+14
# )
# prediction_ci = prediction.conf_int()
# # plt.figure(figsize=(12, 6))
# n = len(data.index)
# ax = data[:].plot(label='actual')
# prediction_ci.plot(
#     ax=ax, style=['--', '--'],
#     label='predicted/forecasted')

# ci_index = prediction_ci.index
# lower_ci = prediction_ci.iloc[:, 0]
# upper_ci = prediction_ci.iloc[:, 1]

# ax.fill_between(ci_index, lower_ci, upper_ci,
#     color='r', alpha=.1)

# ax.set_xlabel('Time (years)')
# ax.set_ylabel('Prices')

# plt.legend()

# # fig = plt.figure()
# # ax1 = fig.add_subplot(111, ylabel = symbol + ' price in $')
# # stock['Close'].plot(ax = ax1, color = 'k', lw = 2., legend = True)
# # stock['Adj Close'].plot(ax = ax1, color = 'g', lw = 2., legend = True)
# plt.show()