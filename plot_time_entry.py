import os, sys
import numpy as np
import datetime, MySQLdb
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 ## Output Type 3 (Type3) or Type 42 (TrueType)
rcParams['font.sans-serif'] = 'Arial'

from matplotlib.dates import MonthLocator, DayLocator, HourLocator, DateFormatter, drange

sys.path.append('../maayanlab_utils')
from plots import COLORS10

def get_times(table):

	conn = MySQLdb.connect(host='localhost',user='root', passwd='',db='maaya0_crowdsourcing')
	cur = conn.cursor()

	cur.execute("""SELECT * FROM %s""" % table)
	times = []
	for row in cur:
		t = row[-2]
		if t is None:
			times.append(datetime.datetime(2015,2,11,15))
		if t is not None:
			times.append(t)
	return times

fig, ax = plt.subplots(figsize=(20,6))
t_genes = get_times('geo2enrichr')
t_dzs = get_times('geo2enrichr_dz')
t_drugs = get_times('geo2enrichr_drug')

ax.plot_date(t_genes, np.arange(len(t_genes)), color=COLORS10[0], label='gene', ls='-', marker=None, lw=2)
ax.plot_date(t_dzs, np.arange(len(t_dzs)), color=COLORS10[1], label='disease', ls='-', marker=None, lw=2)
ax.plot_date(t_drugs, np.arange(len(t_drugs)), color=COLORS10[2], label='drug', ls='-', marker=None, lw=2)

# ax.axhline(y=2422, color='r')
ax.xaxis.set_major_locator( MonthLocator() )
# ax.xaxis.set_minor_locator( HourLocator(np.arange(0,25,6)) )
ax.xaxis.set_major_formatter( DateFormatter('%Y-%m') )
ax.fmt_xdata = DateFormatter('%Y-%m-%d %H:%M:%S')

ax.legend(prop={'size':18})
ax.set_ylabel('# NASB microtask entries', fontsize=20)
fig.autofmt_xdate()
plt.show()
