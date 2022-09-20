"""
Copyright (c) 2017, Sandia National Labs, SunSpec Alliance and CanmetENERGY
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
from svpelab import pvsim
from svpelab import das
from svpelab import der
from svpelab import hil
from svpelab import result as rslt
import script
import numpy as np
import collections



def p_target(f, f_nom, hz_start, hz_stop, pwr):
    """
    Interpolation function to find the target power (using the FW parameter definition)
    :param value: freq point for the interpolation
    :param f: FW freq points
    :param p: FW power points
    :return: target power
    """
    f_pct = 100.*(f/f_nom)
    hz_start_pct = 100.*(hz_start/f_nom)
    hz_stop_pct = 100.*(hz_stop/f_nom)
    if f_pct < hz_start_pct:
        p_targ = 100.*pwr
    elif f_pct > hz_stop_pct:
        p_targ = 0.
    else:
        p_targ = (pwr - pwr*((f_pct-hz_start_pct)/(hz_stop_pct-hz_start_pct)))*100
    return float(p_targ)


def fw_point_interp(value, f, p):
    """
    Interpolation function to find the target power (using the FW pointwise definition)

    :param value: freq point for the interpolation
    :param f: FW freq points
    :param p: FW power points
    :return: target power
    """
    if value <= f[0]:  # if freq is less than f[0]
        return float(p[0])
    elif value >= f[-1]:  # if freq is greater than f[end]
        return float(p[-1])
    else:
        for i in range(len(f)):
            if f[i] <= value <= f[i+1]:                
               p_value = p[i] - ((p[i] - p[i+1])/(f[i+1] - f[i]) * (value - f[i]))
               return float(p_value)
            else:
                ts.log_warning('Unable to Interpolate FW function. f=%s, p=%s, value=%s' % (f, p, value))
                ts.log_warning('Returning nominal power...')
                return 100.


def p_msa_range(f_value, f_msa, p_msa_pct, pwr, f_nom, fw_mode, f=None, p=None, f_slope_start=None, f_slope_stop=None):
    """
    Determine power target and the min/max p values for pass/fail acceptance based on manufacturer's specified
    accuracies (MSAs).

    :param f_value: measured freq value
    :param f_msa: manufacturer's specified accuracy of freq
    :param p_msa_pct: manufacturer's specified accuracy of power in percentage
    :param pwr: power level 
    :param f_nom: EUT nominal freq
    :param fw_mode: FW F1
    :param f: FW freq points (list)
    :param p: FW power points (list)
    :param f_slope_start: FW F1
    :param f_slope_stop: FW F2
    :return: points for p_target, p_target_min, p_target_max
    """
    if fw_mode == 'Pointwise':
        p = [x * pwr for x in p]
        p_targ = fw_point_interp(f_value, f, p)      # target power for the voltage measurement
        p1 = fw_point_interp(f_value - f_msa, f, p)  # power target from the lower voltage limit
        p2 = fw_point_interp(f_value + f_msa, f, p)  # power target from the upper voltage limit        
    else:  # Parameters
        p_targ = p_target(f_value, f_nom, f_slope_start, f_slope_stop,pwr)      # target power for freq
        p1 = p_target(f_value - f_msa, f_nom, f_slope_start, f_slope_stop,pwr)  # power target from the lower freq
        p2 = p_target(f_value + f_msa, f_nom, f_slope_start, f_slope_stop,pwr)  # power target from the upper freq
    if p1 >= p_targ:
        # if the FW curve has a negative slope
        # add the power MSA to the high side (left point, p1)
        # subtract the power MSA from the low side (right point, p2)
        #
        #                          \ * (f_value - f_msa, p_upper)
        #                           \
        #                            . (f_value - f_msa, p1)
        #                             \
        #                              x (f_value, p_target)
        #                               \
        #                                . (f_value + f_msa, p2)
        #                                 \
        #     (f_value + f_msa, p_lower) * \

        p_upper = round(p1 + p_msa_pct, 1)
        p_lower = round(p2 - p_msa_pct, 1)
        return p_targ, p_lower, p_upper
    else:
        p_lower = round(p1 - p_msa_pct, 1)
        p_upper = round(p2 + p_msa_pct, 1)
        return p_targ, p_lower, p_upper

def power_total(data, phases):
    """
    Sum the EUT power from all phases
    :param data: dataset
    :param phases: number of phases in the EUT
    :return: total EUT power
    """
    if phases == 'Single phase':
        ts.log_debug('Powers are: %s' % (data.get('AC_P_1')))
        power = data.get('AC_P_1')

    elif phases == 'Split phase':
        ts.log_debug('Powers are: %s, %s' % (data.get('AC_P_1'),
                                             data.get('AC_P_2')))
        power = data.get('AC_P_1') + data.get('AC_P_2')

    elif phases == 'Three phase':
        ts.log_debug('Powers are: %s, %s, %s' % (data.get('AC_P_1'),
                                                 data.get('AC_P_2'),
                                                 data.get('AC_P_3')))
        power = data.get('AC_P_1') + data.get('AC_P_2') + data.get('AC_P_3')
    else:
        ts.log_error('Inverter phase parameter not set correctly.')
        raise

    return abs(power)
def f_mean(data, phases):
    """
    Average the EUT frequency from all phases
    :param data: dataset
    :param phases: number of phases in the EUT
    :return: mean EUT frequency
    """
    if phases == 'Single phase':
        freq = data.get('AC_FREQ_1')
    elif phases == 'Split phase':
        freq = (data.get('AC_FREQ_1') + data.get('AC_FREQ_2'))/2
    elif phases == 'Three phase':
        freq = (data.get('AC_FREQ_1') + data.get('AC_FREQ_2') + data.get('AC_FREQ_3'))/3
    else:
        ts.log_error('Inverter phase parameter not set correctly.')
        raise

    return freq
def test_run():

    result = script.RESULT_FAIL
    daq = None
    grid = None
    pv = None
    eut = None
    chil = None
    result_summary = None
    f_steps_dic = collections.OrderedDict()

    try:
        """
        Configuration
        """
        fw_mode = ts.param_value('eut_fw.fw_mode')
        f_nom = ts.param_value('eut_fw.f_nom')
        f_min = ts.param_value('eut_fw.f_min')
        f_max = ts.param_value('eut_fw.f_max')
        MSAHz = ts.param_value('eut_fw.MSAHz')
        MSA_P = ts.param_value('eut_fw.MSAP')
        t_settling = ts.param_value('eut_fw.ts')
        fstart_min = ts.param_value('eut_fw.fstart_min')
        fstart_max = ts.param_value('eut_fw.fstart_max')
        k_pf_min = ts.param_value('eut_fw.k_pf_min')
        k_pf_max = ts.param_value('eut_fw.k_pf_max')
        k_pf = ts.param_value('eut_fw.kpf')  
        phases = ts.param_value('eut_fw.phases')
        p_rated = ts.param_value('eut_fw.p_rated')
        MSA_P_pct = (MSA_P/p_rated)*100
        eff = {
                1.0 : ts.param_value('eut_fw.efficiency_100')/100,
                0.66:ts.param_value('eut_fw.efficiency_66')/100,
                0.33:ts.param_value('eut_fw.efficiency_33')/100,
        }        
        n_points = ts.param_value('fw.n_points')
        irr = ts.param_value('fw.irr')
        n_iterations = ts.param_value('fw.n_iter')
        curves = ts.param_value('fw.curves')
        """
        Equipment Configuration
        """         
        # initialize hardware-in-the-loop environment (if applicable)
        ts.log('Configuring HIL system...')
        chil = hil.hil_init(ts)
        if chil is not None:
            chil.config()
            
        # initialize grid simulator
        grid = gridsim.gridsim_init(ts)
        
        # initialize pv simulator
        pv = pvsim.pvsim_init(ts)
        pv.power_set(p_rated)
        pv.power_on()
        # DAS soft channels
        das_points = {'sc': ('P_TARGET', 'P_TARGET_MIN', 'P_TARGET_MAX','P_ACT','F_TARGET','F_ACT','event')}
        
        # initialize data acquisition system
        daq = das.das_init(ts, sc_points=das_points['sc'])
        daq.sc['P_TARGET'] = 100
        daq.sc['P_TARGET_MIN'] = 100
        daq.sc['P_TARGET_MAX'] = 100
        daq.sc['F_TARGET'] = f_nom
        daq.sc['event'] = 'None'
        """
        EUT Configuration
        """        
        # Configure the EUT communications
        eut = der.der_init(ts)
        
        if eut is not None:
            eut.config()
            ts.log_debug(eut.measurements())
            ts.log_debug('L/HFRT and trip parameters set to the widest range : f_min:{0} Hz, f_max:{1} Hz'.format(f_min,f_max))
            
            # TODO : Need to update FRT parameters with SunSpec Model reference
            eut_response = eut.frt_stay_connected_high(params={'Ena' : True,'ActCrv':0, 'Tms1':3000,'Hz1' : f_max,'Tms2':160,'Hz2' : f_max})
            ts.log_debug('HFRT and trip parameters from EUT : {}'.format(eut_response))
            eut_response = eut.frt_stay_connected_low(params={'Ena' : True,'ActCrv':0, 'Tms1':3000,'Hz1' : f_min,'Tms2':160,'Hz2' : f_min})
            ts.log_debug('LFRT and trip parameters from EUT : {}'.format(eut_response))

        else:
            ts.log_debug('Set L/HFRT and trip parameters to the widest range of adjustability possible.')
        """
        Test Configuration
        """           
        if curves == 'Both':
            fw_curves = [1, 2]
        elif curves == 'Characteristic Curve 1':
            fw_curves = [1]
        else:  # Characteristic Curve 2
            fw_curves = [2]

        if irr == 'All':
            pv_powers = [1., 0.66, 0.33]
        elif irr == '100%':
            pv_powers = [1.]
        elif irr == '66%':
            pv_powers = [0.66]
        ts.log_debug("Power level tested : {}".format(pv_powers))
        # open result summary file
        result_summary_filename = 'result_summary.csv'
        result_summary = open(ts.result_file_path(result_summary_filename), 'a+')
        ts.result_file(result_summary_filename)
        result_summary.write('Result,Test Name,Power Level,Iteration,direction,Freq_target,Freq_actual,Power_target,Power_actual,P_min,P_max,Dataset File\n')
        """
        Test start
        """
        for fw_curve in fw_curves:
            if fw_curve == 1:  # characteristic curve 1
                hz_stop = round(fstart_min + 100./k_pf_max,3)
                hz_start = fstart_min
                kp_pf = k_pf_max
            else:  # characteristic curve 2
                hz_stop = round(fstart_max + 100./k_pf_min,3)
                hz_start = fstart_max
                kp_pf = k_pf_min
                
            f_points = [hz_start, hz_stop]
            p_points = [100, 0]
            
            # TODO : Need to update FW parameters with SunSpec Model reference
            if eut is not None: 
                if fw_mode == 'Parameters':
                    eut.freq_watt_param(params={'Ena': True,'HysEna': False, 'HzStr': hz_start,'WGra': kp_pf})
                else:  # Pointwise
                    eut.freq_watt(params={'ActCrv': 1})
                    f_points = [hz_start, hz_stop]
                    p_points = [100, 0]
                    parameters = {'hz': f_points, 'w': p_points}
                    ts.log_debug(parameters)
                    eut.freq_watt_curve(id=1, params=parameters)
                    eut.freq_watt(params={'Ena': True})
                ts.log_debug('Initial EUT FW settings are %s' % eut.freq_watt())
            if fw_mode == 'Parameters':
                ts.log('Testing FW function with pointwise mode : f_start = {0} Hz, slope = {1} %'.format(hz_start,kp_pf))
            else:
                ts.log('Testing FW function with parameters mode : f_start = {0} Hz, f_stop = {1} Hz'.format(hz_start,hz_stop))
                
            for power in pv_powers:
                # SA14.3.2(d) to (g)
                if hz_stop < f_max:
                    if fw_mode == 'Parameters':
                        hz_stop = round(hz_start + (power*100.)/kp_pf,3)                        
                    f_steps_up =   list(np.linspace(f_min+MSAHz, hz_start-MSAHz, n_points)) + list(np.linspace(hz_start+MSAHz, hz_stop-MSAHz, 2*n_points)) + list(np.linspace(hz_stop+MSAHz, f_max - MSAHz, n_points))
                    f_steps_dic['up'] = np.around(f_steps_up,  decimals=3)
                else :
                    f_steps_up =   list(np.linspace(f_min+MSAHz, hz_start-MSAHz, n_points)) + list(np.linspace(hz_start+MSAHz, f_max - MSAHz, 2*n_points))
                    f_steps_dic['up'] = np.around(f_steps_up,  decimals=3)
                f_steps_down = list(reversed(f_steps_up))
                f_steps_dic['down'] = np.around(f_steps_down,  decimals=3)         
                ts.log('Testing FW function at the following frequency up points %s' % f_steps_dic['up'])
                ts.log('Testing FW function at the following frequency down points %s' % f_steps_dic['down'])
                
                # Set to the power level and apply effiency correction
                pv_power_setting = (p_rated*power)/eff[power]
                ts.log('Set PV simulator power to {} with efficiency at {} %'.format(p_rated*power,eff[power]*100.))
                pv.power_set(pv_power_setting)
                for n_iter in range(n_iterations):
                    dataset_filename = 'FW_%s_pwr_%0.2f_iter_%s.csv' % (fw_curve, power, n_iter + 1)
                    daq.data_capture(True)
                    for direction ,f_steps in f_steps_dic.iteritems():
                        for f_step in f_steps:
                            ts.log('        Recording power at frequency %0.3f Hz for 2*t_settling = %0.1f sec.' %
                                   (f_step, 2*t_settling))
                            p_targ, p_min, p_max = p_msa_range(f_value=f_step,
                                                               f_msa=MSAHz,
                                                               p_msa_pct=MSA_P_pct,
                                                               pwr=power,
                                                               f_nom=f_nom,
                                                               fw_mode=fw_mode,
                                                               f=f_points,
                                                               p=p_points,
                                                               f_slope_start=hz_start,
                                                               f_slope_stop=hz_stop)
                            daq.sc['F_TARGET'] = f_step
                            daq.sc['P_TARGET'] = p_targ
                            daq.sc['P_TARGET_MIN'] = p_min
                            daq.sc['P_TARGET_MAX'] = p_max
                            daq.sc['event'] = 'f_Step_{}'.format(direction)
                            ts.log('        Powers targ, min, max: %s, %s, %s' % (p_targ, p_min, p_max))
                            grid.freq(f_step)
                            ts.sleep(2 * t_settling)
                            daq.data_sample()
                            data = daq.data_capture_read()
                            #This is to run the FW script with manual DAQ
                            try:
                                AC_W_pct = (power_total(data, phases) / p_rated) * 100.
                                F_ACT = f_mean(data, phases)
                                if daq.sc['P_TARGET_MIN'] <= AC_W_pct <= daq.sc['P_TARGET_MAX']:
                                    passfail = 'Pass'
                                else:
                                    passfail = 'Fail'
                            except:
                                AC_W_pct = 'No Data'
                                F_ACT= 'No Data'
                                passfail = 'Fail'
                            daq.sc['P_ACT'] = AC_W_pct
                            daq.sc['F_ACT'] = F_ACT
                            daq.sc['event'] = 'T_settling_done_{}'.format(direction)
                            
                            daq.data_sample()
                            ts.log('        Powers measured: {} %'.format(AC_W_pct))
                            result_summary.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n' %
                                                 (passfail,
                                                  ts.config_name(),
                                                  power*100.,
                                                  n_iter+1,
                                                  direction,
                                                  f_step,
                                                  F_ACT,
                                                  p_targ,
                                                  AC_W_pct,
                                                  daq.sc['P_TARGET_MIN'],
                                                  daq.sc['P_TARGET_MAX'],
                                                  dataset_filename))
                    # Parameter for the plotting
                    result_params = {
                    'plot.title': 'title_name',
                    'plot.x.title': 'Time (sec)',
                    'plot.x.points': 'TIME',
                    'plot.y.points': 'F_TARGET,F_ACT',
                    'plot.y.title': 'Frequency (Hz)',
                    'plot.F_TARGET.point': 'True',
                    'plot.y2.points': 'P_TARGET,P_ACT',                    
                    'plot.P_TARGET.point': 'True',
                    'plot.P_TARGET.min_error': 'P_TARGET_MIN',
                    'plot.P_TARGET.max_error': 'P_TARGET_MAX',
                    }
                    
                    daq.data_capture(False)
                    ds = daq.data_capture_dataset()
                    ts.log('Saving file: %s' % dataset_filename)
                    ds.to_csv(ts.result_file_path(dataset_filename))
                    result_params['plot.title'] = os.path.splitext(dataset_filename)[0]
                    ts.result_file(dataset_filename, params=result_params)

        result = script.RESULT_COMPLETE
    except script.ScriptFail, e:
        reason = str(e)
        if reason:
            ts.log_error(reason)
    finally:
        if daq is not None:
            daq.close()
        if pv is not None:
            if p_rated is not None:
                pv.power_set(p_rated)
            pv.close()
        if grid is not None:
            if f_nom is not None:
                grid.freq(f_nom)
            grid.close()
        if chil is not None:
            chil.close()
        if eut is not None:
            eut.close()
        if result_summary is not None:
            result_summary.close()

        # create result workbook
        xlsxfile = ts.config_name() + '.xlsx'
        rslt.result_workbook(xlsxfile, ts.results_dir(), ts.result_dir())
        ts.result_file(xlsxfile)

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
        
        ts.svp_version(required='1.5.9')
        result = test_run()

        ts.result(result)
        if result == script.RESULT_FAIL:
            rc = 1

    except Exception, e:
        ts.log_error('Test script exception: %s' % traceback.format_exc())
        rc = 1

    sys.exit(rc)

info = script.ScriptInfo(name=os.path.basename(__file__), run=run, version='1.0.0')

der.params(info)
# EUT FW parameters
info.param_group('eut_fw', label='FW Configuration',glob=True)
info.param('eut_fw.fw_mode', label='Freq-Watt Mode', default='Parameters',values=['Parameters', 'Pointwise'],desc='Parameterized FW curve or pointwise linear FW curve?')
info.param('eut_fw.p_rated', label='Output Power Rating (W)', default=34500.)
info.param('eut_fw.f_nom', label='Nominal AC frequency (Hz)', default=60.)
info.param('eut_fw.f_min', label='Min AC frequency (Hz)', default=49.)
info.param('eut_fw.f_max', label='Max AC frequency (Hz)', default=52.)
info.param('eut_fw.MSAHz', label='Manufacturer\'s stated AC frequency accuracy (Hz)', default=0.1)
info.param('eut_fw.MSAP', label='Manufacturer\'s stated power accuracy (W)', default=10.)
info.param('eut_fw.ts', label='Settling time (s)', default=1.)
info.param('eut_fw.fstart_min', label='Min start of frequency droop (Hz)', default=50.1)
info.param('eut_fw.fstart_max', label='Max start of frequency droop (Hz)', default=51.)
info.param('eut_fw.k_pf_min', label='Min slope of frequency droop (%Prated/Hz)', default=20.)
info.param('eut_fw.k_pf_max', label='Max slope of frequency droop (%Prated/Hz)', default=50.)
info.param('eut_fw.phases', label='Phases', default='Single Phase', values=['Single phase', 'Split phase', 'Three phase'])
info.param('eut_fw.efficiency_33', label='CEC Efficiency list for power level = 33% at nominal VDC', default=97.0)
info.param('eut_fw.efficiency_66', label='CEC Efficiency list for power level = 66% at nominal VDC', default=97.1)
info.param('eut_fw.efficiency_100', label='CEC Efficiency list for power level = 100% at nominal VDC', default=96.9)

info.param_group('fw', label='Test Parameters')
info.param('fw.curves', label='Curves to Evaluate', default='Both',values=['Characteristic Curve 1', 'Characteristic Curve 2', 'Both'])
info.param('fw.irr', label='Power Levels', default='All',values=['100%', '66%', '33%', 'All'])
info.param('fw.n_iter', label='Number of iteration for each test', default=3)
info.param('fw.n_points', label='Number of points tested above f_start', default=3)

gridsim.params(info)
pvsim.params(info)
das.params(info)
hil.params(info)

# info.logo('sunspec.gif')

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


