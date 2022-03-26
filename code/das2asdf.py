# Purpose:	Convert DAS data to ASDF (h5) format
# Author:	Jun ZHU, modified from Obspy tutorial
# Date:		SEP 13 2021
# Email:	Jun__Zhu@outlook.com


import glob

import pyasdf
from obspy import Stream, Trace, UTCDateTime
from obspy.core.trace import Stats
from obspy.core.inventory import (Inventory, Network, Station, Channel, Site)
from obspy.clients.nrl import NRL

import numpy as np
import matplotlib.pyplot as plt

from read_rd import read_PASSCAL_segy as rsgy


def Das2ASDF(root):
	dirlist = glob.glob(root+'*/')
	nTrace = 1250
	fs = 50
	# create an inventory object for the das array
	inv = Inventory(
			networks=[],
			source='Caltech')
	net = Network(code='RDG',
			stations=[],
			description='A DAS array across the Ridgcrest Fault',
			start_date=UTCDateTime(2019, 7, 6))
	sta = Station(
			code='DAS',
			latitude=35.75,
			longitude=-117.5,
			elevation=0,
			creation_date=UTCDateTime(2021, 9, 13),
			site=Site(name="Ridgecrest Fault"))
	for itrace in range(nTrace):
		cha = Channel(
					code="%04d"%itrace,
					location_code="",
					latitude=1.0,
					longitude=2.0,
					elevation=0,
					depth=0,
					sample_rate=100)
		sta.channels.append(cha)
	nrl = NRL() # used for response
	net.stations.append(sta)
	inv.networks.append(net)
	#
	for dir in dirlist:
		eventid = dir.split('/')[-2]
		ds = pyasdf.ASDFDataSet(dir+eventid+'.h5')
		fsegy = glob.glob(dir+'*.segy')[0]
		fstation = glob.glob(dir+'*.stainfo')[0]
		fevent = glob.glob(dir+'*.info')[0]
		data = rsgy(fsegy, fs, nTrace)
		st = Stream()
		with open(fstation, 'r') as f:
			meta = [x.split(' ') for x in f.readlines()]
			id_chan_list = [int(x[2]) for x in meta]
			snr_list = [float(x[3]) for x in meta]
			arrival_list = [float(x[4][:-1]) for x in meta]
		meta = {'id': [], 'arrival': [], 'snr': [], 'quality': []}
		BadSensorDefaultValue = {'snr': None, 'arrival': None}
		for id_chan in range(nTrace):
			if id_chan in id_chan_list:
				idx = id_chan_list.index(id_chan)
				meta['id'].append(id_chan_list[idx])
				meta['arrival'].append(arrival_list[idx])
				meta['snr'].append(snr_list[idx])
				meta['quality'].append('operational')
			else:
				meta['id'].append(id_chan)
				meta['arrival'].append(BadSensorDefaultValue['arrival'])
				meta['snr'].append(BadSensorDefaultValue['snr'])
				meta['quality'].append('absent')
	#	# for the beauty of displayed image, the data is normalized by each
	#	# sensor's maximum absolute amplitude
	#	plt.imshow(np.transpose(data/np.max(np.abs(data), axis=1)[:, None]), vmin=-0.5, vmax=0.5, aspect='auto', cmap='RdBu')
	#	plt.scatter(id_chan_list, [x+500 for x in arrival_list], c='black', s=1)
	#	plt.xlabel('Channel ID')
	#	plt.ylabel('Sample Point')
	#	plt.show()
		for itrace in range(data.shape[0]):
			stats = Stats()
			stats.experiment = 'Ridgecrest'
			stats.sampling_rate = fs
			stats.npts = len(data[itrace])
			stats.network = 'RDG'
			stats.station = 'DAS'
			stats.channel = "%04d"%meta['id'][itrace]
			stats.arrival = meta['arrival'][itrace]
			stats.quality = meta['quality'][itrace]
			stats.instrument = {'DAS': {'interrogator': None, 'fiber':
				None, 'cable': None, 'fiber_length': -999, 'gauge_length':
				-999, 'pulse_repetition_rate': -999, 'laser_wavelength': -999,
				'spatial_averaging': -999, 'deployed_environment':None}}
			tr = Trace(data=data[itrace], header=stats)
			st += tr
		ds.add_waveforms(st, tag='raw')
		ds.add_stationxml(inv)
		# for test
		ds.waveforms.RDG_DAS.raw[-10:].plot()
	return


if __name__ == "__main__":
	root = '../data/'
	Das2ASDF(root)
