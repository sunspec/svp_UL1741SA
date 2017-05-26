"""
Copyright (c) 2017, Sandia National Labs and SunSpec Alliance
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

Neither the names of the Sandia National Labs and SunSpec Alliance nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions can be directed to support@sunspec.org
"""

import time

try:
    import typhoon.api.hil_control_panel as cp
    from typhoon.api.schematic_editor import model
    import typhoon.api.pv_generator as pv
except Exception, e:
    print('Typhoon HIL API not installed. %s' % e)

data_points = [
    'TIME',
    'DC_V',
    'DC_I',
    'AC_VRMS',
    'AC_IRMS',
    'DC_P',
    'AC_S',
    'AC_P',
    'AC_Q',
    'AC_FREQ',
    'AC_PF',
    'TRIG',
    'TRIG_GRID'
]

# To be implemented later
# typhoon_points_asgc_1 = [
#     'time',
#     'V( V_DC3 )', # DC voltage
#     'I( Ipv )',
#     'V( Vrms1 )',
#     'I( Irms1 )',
#     'DC_P',  # calculated
#     'S',
#     'Pdc',
#     'Qdc',
#     'AC_FREQ',
#     'k',
#     'TRIG',
#     'TRIG_GRID'
# ]
#
# typhoon_points_asgc_3 = [
#     'time',
#     'V( V_DC3 )',  # DC voltage
#     'I( Ipv )',
#     'V( Vrms1 )',
#     'V( Vrms2 )',
#     'V( Vrms3 )',
#     'I( Irms1 )',
#     'I( Irms2 )',
#     'I( Irms3 )',
#     'DC_P',  # calculated
#     'S',
#     'Pdc',
#     'Qdc',
#     'AC_FREQ',
#     'k',
#     'TRIG',
#     'TRIG_GRID'
# ]
#
# typhoon_points_map = {
#     'ASGC3': typhoon_points_asgc_3,  # AGF circuit, 3 phase
#     'ASGC1': typhoon_points_asgc_1,  # AGF circuit, single phase
#     'ASGC_Fault': typhoon_points_asgc_fault,  # ride-through circuit
#     'ASGC_UI': typhoon_points_ui   # unintentional islanding circuit
# }

wfm_channels = ['AC_V_1', 'AC_V_2', 'AC_V_3', 'AC_I_1', 'AC_I_2', 'AC_I_3', 'EXT']

wfm_typhoon_channels = {'AC_V_1': 'V( Vrms1 )',
                        'AC_V_2': 'V( Vrms2 )',
                        'AC_V_3': 'V( Vrms3 )',
                        'AC_I_1': 'I( Irms1 )',
                        'AC_I_2': 'I( Irms2 )',
                        'AC_I_3': 'I( Irms3 )',
                        'EXT': 'Trigger',
                        'V( Vrms1 )': 'AC_V_1',
                        'V( Vrms2 )': 'AC_V_2',
                        'V( Vrms3 )': 'AC_V_3',
                        'I( Irms1 )': 'AC_I_1',
                        'I( Irms2 )': 'AC_I_2',
                        'I( Irms3 )': 'AC_I_3',
                        'Trigger': 'EXT'}

event_map = {'Rising_Edge': 'Rising edge', 'Falling_Edge': 'Falling edge'}

class Device(object):

    def __init__(self, params=None):
        self.params = params
        self.dsm_method = self.params.get('dsm_method')
        self.dsm_id = self.params.get('dsm_id')
        self.comp = self.params.get('comp')
        self.file_path = self.params.get('file_path')
        self.data_file = os.path.join(self.file_path, DATA_FILE)
        self.points_file = os.path.join(self.file_path, POINTS_FILE)
        self.wfm_trigger_file = os.path.join(self.file_path, WFM_TRIGGER_FILE)

        self.data_points = data_points
        self.points_map = typhoon_points_map.get(str(self.dsm_id))
        self.points = None
        self.point_indexes = []

        self.rec = {}
        self.recs = []

        self.read_error_count = 0
        self.read_last_error = ''

        # waveform settings
        self.wfm_sample_rate = None
        self.wfm_pre_trigger = None
        self.wfm_post_trigger = None
        self.wfm_trigger_level = None
        self.wfm_trigger_cond = None
        self.wfm_trigger_channel = None
        self.wfm_timeout = None
        self.wfm_channels = None
        self.wfm_capture_name = None
        self.wfm_capture_name_path = r'C:\captured_signals\capture_test.mat'

        self.ts = self.params.get('ts')

        self.numberOfSamples = self.sample_rate*(self.pre_trigger+self.post_trigger)
        self.decimation = 1
        self.captureSettings = None
        self.triggerOffset = (self.pre_trigger/(self.pre_trigger+self.post_trigger))*100.
        self.triggerSettings = None
        self.channelSettings = None

        # regular python list is used for data buffer
        self.capturedDataBuffer = None
        # data containers
        self.time_vector = None
        self.wfm_data = None
        self.signalsNames = None

    def info(self):
        hw = model.get_hw_settings()
        return 'HIL hardware version: %s' % (hw,)

    def open(self):
        pass

    def close(self):
        pass

    def data_read(self):

        v1 = float(cp.read_analog_signal(name='V( Vrms1 )'))
        v2 = float(cp.read_analog_signal(name='V( Vrms2 )'))
        v3 = float(cp.read_analog_signal(name='V( Vrms3 )'))
        i1 = float(cp.read_analog_signal(name='I( Irms1 )'))
        i2 = float(cp.read_analog_signal(name='I( Irms2 )'))
        i3 = float(cp.read_analog_signal(name='I( Irms3 )'))
        p = float(cp.read_analog_signal(name='Pdc'))  # Note this is the AC power (fundamental)
        va = float(cp.read_analog_signal(name='S'))
        q = float(cp.read_analog_signal(name='Qdc'))
        pf = float(cp.read_analog_signal(name='k'))
        # f = cp.frequency

        dc_v = float(cp.read_analog_signal(name='V( V_DC3 )'))
        dc_i = float(cp.read_analog_signal(name='I( Ipv )'))

        datarec = {'TIME': time.time(),
                   'AC_VRMS_1': v1,
                   'AC_IRMS_1': i1,
                   'AC_P_1': p/3.,
                   'AC_S_1': va/3.,
                   'AC_Q_1': q/3.,
                   'AC_PF_1': pf,
                   'AC_FREQ_1': None,
                   'AC_VRMS_2': v2,
                   'AC_IRMS_2': i2,
                   'AC_P_2': p/3.,
                   'AC_S_2': va/3.,
                   'AC_Q_2': q/3.,
                   'AC_PF_2': pf,
                   'AC_FREQ_2': None,
                   'AC_VRMS_3': v3,
                   'AC_IRMS_3': i3,
                   'AC_P_3': p/3.,
                   'AC_S_3': va/3.,
                   'AC_Q_3': q/3.,
                   'AC_PF_3': pf,
                   'AC_FREQ_3': None,
                   'DC_V': dc_v,
                   'DC_I': dc_i,
                   'DC_P': dc_i*dc_v}

        return datarec

    def waveform_config(self, params):
        """
        Configure waveform capture.

        params: Dictionary with following entries:
            'sample_rate' - Sample rate (samples/sec)
            'pre_trigger' - Pre-trigger time (sec)
            'post_trigger' - Post-trigger time (sec)
            'trigger_level' - Trigger level
            'trigger_cond' - Trigger condition - ['Rising_Edge', 'Falling_Edge']
            'trigger_channel' - Trigger channel - ['AC_V_1', 'AC_V_2', 'AC_V_3', 'AC_I_1', 'AC_I_2', 'AC_I_3', 'EXT']
            'timeout' - Timeout (sec)
            'channels' - Channels to capture - ['AC_V_1', 'AC_V_2', 'AC_V_3', 'AC_I_1', 'AC_I_2', 'AC_I_3', 'EXT']
        """
        self.wfm_sample_rate = params.get('sample_rate')
        self.wfm_pre_trigger = params.get('pre_trigger')
        self.wfm_post_trigger = params.get('post_trigger')
        self.wfm_trigger_level = params.get('trigger_level')
        self.wfm_trigger_cond = params.get('trigger_cond')
        self.wfm_trigger_channel = params.get('trigger_channel')
        self.wfm_timeout = params.get('timeout')
        self.wfm_channels = params.get('channels')

        self.captureSettings = [self.decimation, len(self.channels), self.numberOfSamples]

        # triggerType,triggerSource,threshold,edge,triggerOffset
        self.triggerSettings = ["Analog", self.trigger_channel, self.trigger_level,
                                event_map[self.wfm_trigger_cond], self.triggerOffset]
        # triggerType - type of trigger ("Analog" , "Digital" or "Forced")

        # signals for capturing
        # self.channelSettings = ['V( Vrms1 )', 'V( Vrms2 )', 'V( Vrms3 )', 'I( Irms1 )', 'I( Irms2 )', 'I( Irms3 )']
        self.channelSettings = []
        for i in range(self.wfm_channels):
            self.channelSettings.append(wfm_typhoon_channels.get(self.wfm_channels[i]))

        # regular python list is used for data buffer
        self.capturedDataBuffer = None  # reset the data buffer

    def waveform_capture(self, enable=True, sleep=None):
        """
        Enable/disable waveform capture.
        """
        if enable:
            # start capture process and if everything ok continue...
            if hil.start_capture(self.captureSettings,
                                 self.triggerSettings,
                                 self.channelSettings,
                                 dataBuffer=self.capturedDataBuffer,
                                 fileName=self.wfm_capture_name_path):

                # when capturing is finished...
                while waveform_status == 'ACTIVE':
                    pass

                # unpack data from data buffer
                # signalsNames - list with names
                # wfm_data  - 'numpy.ndarray' matrix with data values,
                # time_vector - 'numpy.array' with time data
                self.signalsNames, self.wfm_data, self.time_vector = capturedDataBuffer[0]

    def waveform_status(self):
        # return INACTIVE, ACTIVE, COMPLETE
        if hil.capture_in_progress():
            stat = 'ACTIVE'
        elif self.capturedDataBuffer is None:
            stat = 'INACTIVE'
        else:
            stat = 'COMPLETE'
        return stat

    def waveform_force_trigger(self):
        """
        Create trigger event with provided value.
        """
        self.triggerSettings = ["Forced"]

    def waveform_capture_dataset(self):
        ds = dataset.Dataset()

        if len(self.signalsNames) == len(self.channelSettings):
            ds.points.append('TIME')
            ds.data.append(self.time_vector)
            for i in range(self.channelSettings):
                ds.points.append(wfm_typhoon_channels.get(self.signalsNames[i]))
                # unpack data for appropriate captured signals
                ds.data.append(self.wfm_data[i])  # first row for first signal and so on
        else:
            self.ts.log_error('Number of channels returned from waveform capture is unexpected. '
                              'Expected %s. Got: %s' % (self.channelSettings, self.signalsNames))

        return ds


if __name__ == "__main__":

    pass


