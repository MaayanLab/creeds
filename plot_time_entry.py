import os, sys
import numpy as np
import datetime, MySQLdb
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange

conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
cur = conn.cursor()

cur.execute("""SELECT * FROM geo2enrichr""")
times = []
for row in cur:
	t = row[-2]
	if t is None:
		times.append(datetime.datetime(2015,2,11,15))
	if t is not None:
		times.append(t)
fig, ax = plt.subplots(figsize=(20,6))

ax.plot_date(times, np.arange(len(times)), '-')
ax.axhline(y=2422, color='r')
ax.xaxis.set_major_locator( DayLocator() )
ax.xaxis.set_minor_locator( HourLocator(np.arange(0,25,6)) )
ax.xaxis.set_major_formatter( DateFormatter('%m-%d') )
ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')

ax.set_ylabel('# NASB microtask entries', fontsize=20)
fig.autofmt_xdate()
plt.show()
