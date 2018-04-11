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

import sys
import os
import traceback
from svpelab import gridsim
from svpelab import loadsim
from svpelab import pvsim
from svpelab import das
from svpelab import der
import script
import result as rslt

def test_pass_fail(pf_act=None, pf_target=None, pf_msa=None):

    # set pf range and if range spans pf=1
    pf_lower = pf_upper = 1.0
    span = False
    if pf_target < 0:
        offset = 1.0 + pf_target
        if offset <= pf_msa:
            pf_upper = 1.0 - (pf_msa - offset)
            span = True
        else:
            pf_upper = pf_target - pf_msa
        pf_lower = pf_target + pf_msa
    elif pf_target < 1.0:
        offset = 1.0 - pf_target
        if offset <= pf_msa:
            pf_lower = -(1.0 - (pf_msa - offset))
            if pf_lower == -1.0:
                pf_lower = 1.0
            span = True
        else:
            pf_lower = pf_target - pf_msa
        pf_upper = pf_target + pf_msa
    elif pf_target == 1.0:
        pf_upper = 1.0 - pf_msa
        pf_lower = -pf_upper
        span = True

    # check if pf in range
    passfail = 'Fail'
    if ((span and ((pf_act == 1.0) or
                   (pf_act < 0 and abs(pf_act) >= abs(pf_lower)) or
                   (pf_act < 1.0 and pf_act >= pf_upper))) or
        ((pf_act < 0 and pf_act >= pf_upper and pf_act <= pf_lower) or
         (pf_act <= 1.0 and pf_act >= pf_lower and pf_act <= pf_upper))):
        passfail = 'Pass'

    return (passfail, pf_lower, pf_upper)


def test_run():

    result = script.RESULT_FAIL
    grid = None
    pv = None
    daq = None
    eut = None
    rs = None

    sc_points = ['PF_TARGET', 'PF_MAX', 'PF_MIN']

    # result params
    result_params = {
        'plot.title': ts.name,
        'plot.x.title': 'Time (secs)',
        'plot.x.points': 'TIME',
        'plot.y.points': 'AC_PF_1, PF_TARGET',
        'plot.y.title': 'Power Factor',
        'plot.y2.points': 'AC_IRMS_1',
        'plot.y2.title': 'Current (A)'
    }

    try:
        p_rated = ts.param_value('eut.p_rated')
        pf_min_ind = ts.param_value('eut.pf_min_ind')
        pf_min_cap = ts.param_value('eut.pf_min_cap')
        pf_settling_time = ts.param_value('eut.pf_settling_time')
        pf_msa = ts.param_value('eut.pf_msa')

        p_low = p_rated * .2
        pf_mid_ind = (-1 + pf_min_ind)/2
        pf_mid_cap = (1 + pf_min_cap)/2

        '''
        2) Set all AC source parameters to the normal operating conditions for the EUT. 
        '''
        # grid simulator is initialized with test parameters and enabled
        grid = gridsim.gridsim_init(ts)

        # pv simulator is initialized with test parameters and enabled
        pv = pvsim.pvsim_init(ts)
        pv.power_set(p_low)
        pv.power_on()

        # initialize data acquisition
        daq = das.das_init(ts, sc_points = sc_points)
        ts.log('DAS device: %s' % daq.info())

        '''
        3) Turn on the EUT. It is permitted to set all L/HVRT limits and abnormal voltage trip parameters to the
        widest range of adjustability possible with the SPF enabled in order not to cross the must trip
        magnitude threshold during the test.
        '''
        # it is assumed the EUT is on
        eut = der.der_init(ts)
        if eut is not None:
            eut.config()

        '''
        4) Select 'Fixed Power Factor' operational mode.
        '''
        # set power levels that are enabled
        power_levels = []
        if ts.param_value('spf.p_100') == 'Enabled':
            power_levels.append((1, '100'))
        if ts.param_value('spf.p_50') == 'Enabled':
            power_levels.append((.5, '50'))
        if ts.param_value('spf.p_20') == 'Enabled':
            power_levels.append((.2, '20'))

        # set target power factors
        pf_targets = []
        if ts.param_value('spf.pf_min_ind') == 'Enabled':
            pf_targets.append(pf_min_ind)
        if ts.param_value('spf.pf_mid_ind') == 'Enabled':
            pf_targets.append(pf_mid_ind)
        if ts.param_value('spf.pf_min_cap') == 'Enabled':
            pf_targets.append(pf_min_cap)
        if ts.param_value('spf.pf_mid_cap') == 'Enabled':
            pf_targets.append(pf_mid_cap)

        n_r = ts.param_value('spf.n_r')

        # open result summary file
        result_summary_filename = 'result_summary.csv'
        result_summary = open(ts.result_file_path(result_summary_filename), 'a+')
        ts.result_file(result_summary_filename)
        result_summary.write('Result, Test Name, Power Level (%), Iteration, PF Actual, PF Target, PF MSA, PF Min Allowed,'
                             'PF Max Allowed, Dataset File\n')

        for pf in pf_targets:
            for power_level in power_levels:
                '''
                5) Set the input source to produce Prated for the EUT.
                '''
                power, power_label = power_level
                pv.power_set(p_rated * power)
                ts.log('*** Setting power level to %s W (rated power * %s)' % ((p_rated * power), power))

                for count in range(1, n_r + 1):
                    ts.log('Starting pass %s' % (count))
                    '''
                    6) Set the EUT power factor to unity. Measure the AC source voltage and EUT current to measure the
                    displacement
                    '''
                    #ts.log('Fixed PF settings: %s' % eut.fixed_pf())
                    ### use pass_fail to get lower, upper
                    passfail, pf_lower, pf_upper = test_pass_fail(pf_act=1.0, pf_target=1.0, pf_msa=pf_msa)
                    daq.sc['PF_TARGET'] = 1.0
                    daq.sc['PF_MAX'] = pf_upper
                    daq.sc['PF_MIN'] = pf_lower
                    if eut is not None:
                        eut.fixed_pf(params={'Ena': True, 'PF': 1.0})
                        pf_setting = eut.fixed_pf()
                        ts.log('PF setting: %s' % (pf_setting))
                    ts.log('Starting data capture for pf = %s' % (1.0))
                    daq.data_capture(True)
                    ts.log('Sampling for %s seconds' % (pf_settling_time * 3))
                    ts.sleep(pf_settling_time * 3)
                    ts.log('Sampling complete')
                    daq.data_capture(False)
                    ds = daq.data_capture_dataset()
                    filename = 'spf_1000_%s_%s.csv' % (str(power_label), str(count))
                    ds.to_csv(ts.result_file_path(filename))
                    ts.result_file(filename)
                    ts.log('Saving data capture %s' % (filename))

                    '''
                    7) Set the EUT power factor to the value in Test 1 of Table SA12.1. Measure the AC source voltage
                    and EUT current to measure the displacement power factor and record all data.
                    '''
                    ### use pass_fail to get lower, upper
                    passfail, pf_lower, pf_upper = test_pass_fail(pf_act=1.0, pf_target=pf, pf_msa=pf_msa)
                    daq.sc['PF_TARGET'] = pf
                    daq.sc['PF_MAX'] = pf_upper
                    daq.sc['PF_MIN'] = pf_lower
                    if eut is not None:
                        parameters = {'Ena': True, 'PF': pf}
                        ts.log('PF set: %s' % (parameters))
                        eut.fixed_pf(params=parameters)
                        ts.log('PF set exit: %s' % (parameters))
                        pf_setting = eut.fixed_pf()
                        ts.log('PF setting: %s' % (pf_setting))
                    ts.log('Starting data capture for pf = %s' % (pf))
                    daq.data_capture(True)
                    ts.log('Sampling for %s seconds' % (pf_settling_time * 3))
                    ts.sleep(pf_settling_time * 3)
                    ts.log('Sampling complete')
                    daq.data_capture(False)
                    ds = daq.data_capture_dataset()
                    testname = 'spf_%s_%s_%s' % (str(pf * 1000), str(power_label), str(count))
                    filename = '%s.csv' % (testname)
                    ds.to_csv(ts.result_file_path(filename))
                    result_params['plot.title'] = testname
                    ts.result_file(filename, params=result_params)
                    ts.log('Saving data capture %s' % (filename))

                    # create result summary entry
                    pf_point = 'AC_PF_1'
                    pf_data = []
                    try:
                        idx = ds.points.index(pf_point)
                        pf_data = ds.data[idx]
                    except ValueError, e:
                        ts.fail('Data point %s not in dataset' % (pf_point))
                    if len(pf_data) <= 0:
                        ts.fail('No data for data point %s' % (pf_point))

                    # verify starting pf
                    pf_act = float(pf_data[0])

                    passfail, pf_lower, pf_upper = test_pass_fail(pf_act=pf_act, pf_target=1.0, pf_msa=pf_msa)
                    result_summary.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n' %
                                         (passfail, ts.config_name(), power * 100, count, pf_act, pf, pf_msa, pf_lower,
                                          pf_upper, filename))

                    # use the last PF measurement for the pass/fail checks
                    pf_act = float(pf_data[-1])

                    passfail, pf_lower, pf_upper = test_pass_fail(pf_act=pf_act, pf_target=pf, pf_msa=pf_msa)
                    result_summary.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n' %
                                         (passfail, ts.config_name(), power * 100, count, pf_act, pf, pf_msa, pf_lower,
                                          pf_upper, filename))
                    '''
                    8) Repeat steps (6) - (8) for two additional times for a total of three repetitions.
                    '''
                '''
                9) Repeat steps (5) - (7) at two additional power levels. One power level shall be a Pmin or 20% of
                Prated and the second at any power level between 33% and 66% of Prated.
                '''
            '''
            10) Repeat Steps (6) - (9) for Tests 2 - 5 in Table SA12.1
            '''

        '''
        11) In the case of bi-directional inverters, repeat Steps (6) - (10) for the active power flow direction
        '''

        result = script.RESULT_COMPLETE

    except script.ScriptFail, e:
        reason = str(e)
        if reason:
            ts.log_error(reason)
    finally:

        pv.power_set(p_rated)

        if grid is not None:
            grid.close()
        if pv is not None:
            pv.close()
        if daq is not None:
            daq.close()
        if eut is not None:
            eut.close()
        if rs is not None:
            rs.close()

        # create result workbook
        file = ts.config_name() + '.xlsx'
        rslt.result_workbook(file, ts.results_dir(), ts.result_dir())
        ts.result_file(file)

    return result

def run(test_script):

    try:
        global ts
        ts = test_script
        rc = 0
        result = script.RESULT_COMPLETE

        ts.log_debug('')
        ts.log_debug('**************  Starting %s  **************' % (ts.config_name()))
        ts.log_debug('Script: %s %s' % (ts.name, ts.info.version))
        ts.log_active_params()

        result = test_run()

        ts.result(result)
        if result == script.RESULT_FAIL:
            rc = 1

    except Exception, e:
        ts.log_error('Test script exception: %s' % traceback.format_exc())
        rc = 1

    sys.exit(rc)

info = script.ScriptInfo(name=os.path.basename(__file__), run=run, version='1.1.0')

info.param_group('spf', label='Test Parameters')
info.param('spf.p_100', label='Power Level 100% Tests', default='Enabled', values=['Disabled', 'Enabled'])
info.param('spf.p_50', label='Power Level 50% Tests', default='Enabled', values=['Disabled', 'Enabled'])
info.param('spf.p_20', label='Power Level 20% Tests', default='Enabled', values=['Disabled', 'Enabled'])
info.param('spf.n_r', label='Number of test repetitions', default=3)

info.param('spf.pf_min_ind', label='Minimum inductive', default='Enabled', values=['Disabled', 'Enabled'])
info.param('spf.pf_mid_ind', label='Mid-range inductive', default='Enabled', values=['Disabled', 'Enabled'])
info.param('spf.pf_min_cap', label='Minimum capacitive', default='Enabled', values=['Disabled', 'Enabled'])
info.param('spf.pf_mid_cap', label='Mid-range capacitive', default='Enabled', values=['Disabled', 'Enabled'])

info.param_group('eut', label='EUT Parameters', glob=True)
info.param('eut.p_rated', label='P_rated', default=3000)
info.param('eut.phases', label='Phases', default='Single Phase', values=['Single Phase', '3-Phase 3-Wire',
                                                                         '3-Phase 4-Wire'])
info.param('eut.pf_min_ind', label='PF_min_ind', default=.850)
info.param('eut.pf_min_cap', label='PF_min_cap', default=-.850)
info.param('eut.pf_settling_time', label='PF Settling Time (secs)', default=1)
info.param('eut.pf_msa', label='PF Manufacturer Stated Accuracy (PF units)', default=5.0)

der.params(info)
das.params(info)
gridsim.params(info)
loadsim.params(info)
pvsim.params(info)

def script_info():
    
    return info


if __name__ == "__main__":

    # stand alone invocation
    config_file = None
    if len(sys.argv) > 1:
        config_file = sys.argv[1]

    params = None

    test_script = script.Script(info=script_info(), config_file=config_file, params=params)
    test_script.log('log it')

    run(test_script)


