"""
Copyright (c) 2018, Sandia National Labs, SunSpec Alliance and CanmetENERGY
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
from svpelab import hil
from svpelab import gridsim
from svpelab import pvsim
from svpelab import das
from svpelab import der
from svpelab import result as rslt
import numpy as np
import time
import script
import collections


def p_target(v, v_nom, v_slope_start, v_slope_stop, pwr):
    """
    Calculate the target power from the VW curve parameters
    :param v: voltage point
    :param v_nom: nominal EUT voltage
    :param v_slope_start: VW V1
    :param v_slope_stop: VW V2
    :param pwr : power level
    :return: Expected VW power
    """
    v_pct = 100.*(v/v_nom)
    v_slope_start_pct = 100.*(v_slope_start/v_nom)
    v_slope_stop_pct = 100.*(v_slope_stop/v_nom)
    if v_pct < v_slope_start_pct:
        p_targ = 100.*pwr
    elif v_pct > v_slope_stop_pct:
        p_targ = 0.
    else:
        p_targ = (pwr - pwr*((v_pct-v_slope_start_pct)/(v_slope_stop_pct-v_slope_start_pct)))*100
    return p_targ


def p_msa_range(v_value, v_msa, p_msa_pct, pwr, v_nom, v_slope_start, v_slope_stop):
    """
    Determine  power target and the min/max p values for pass/fail acceptance based on manufacturer's specified
    accuracies (MSAs).  This assumes that Q1 = 100% and Q2 = 0% Prated

    :param v_value: measured voltage value
    :param v_msa: manufacturer's specified accuracy of voltage
    :param p_msa_pct: manufacturer's specified accuracy of power in percentage
    :param pwr: power level 
    :param v_nom: EUT nominal voltage
    :param v_slope_start: VW V1
    :param v_slope_stop: VW V2
    :return: points for q_target, q_target_min, q_target_max
    """
    p_targ = p_target(v_value, v_nom, v_slope_start, v_slope_stop, pwr)    # target power for the voltage measurement
    p1 = p_target(v_value - v_msa, v_nom, v_slope_start, v_slope_stop, pwr)  # power target from the lower voltage limit
    p2 = p_target(v_value + v_msa, v_nom, v_slope_start, v_slope_stop, pwr)  # power target from the upper voltage limit
    if p1 >= p_targ:
        # if the VW curve has a negative slope
        # add the power MSA to the high side (left point, p1)
        # subtract the power MSA from the low side (right point, p2)
        #
        #                          \ * (v_value - v_msa, p_upper)
        #                           \
        #                            . (v_value - v_msa, p1)
        #                             \
        #                              x (v_value, p_target)
        #                               \
        #                                . (v_value + v_msa, p2)
        #                                 \
        #     (v_value + v_msa, p_lower) * \

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

def v_mean(data, phases):
    """
    Average the EUT voltage from all phases
    :param data: dataset
    :param phases: number of phases in the EUT
    :return: mean EUT voltage
    """
    if phases == 'Single phase':
        volt = data.get('AC_VRMS_1')
    elif phases == 'Split phase':
        volt = (data.get('AC_VRMS_1') + data.get('AC_VRMS_2'))/2
    elif phases == 'Three phase':
        volt = (data.get('AC_VRMS_1') + data.get('AC_VRMS_2') + data.get('AC_VRMS_3'))/3
    else:
        ts.log_error('Inverter phase parameter not set correctly.')
        raise

    return volt

def test_run():

    result = script.RESULT_FAIL
    daq = None
    data = None
    grid = None
    pv = None
    eut = None
    chil = None
    result_summary = None
    v_steps_dic = collections.OrderedDict()

    try:
        """
        Configuration
        """
        # Initiliaze VW EUT specified parameters variables
        vw_mode = ts.param_value('eut_vw.vw_mode')
        phases = ts.param_value('eut_vw.phases')
        v_nom = ts.param_value('eut_vw.v_nom')
        v_min = ts.param_value('eut_vw.v_min')
        v_max = ts.param_value('eut_vw.v_max')
        MSA_V = ts.param_value('eut_vw.MSA_V')
        MSA_P = ts.param_value('eut_vw.MSA_P')
        t_settling = ts.param_value('eut_vw.ts')
        v_start_max = ts.param_value('eut_vw.vstart_max')
        v_start_min = ts.param_value('eut_vw.vstart_min')
        # the k_power_slopes are the slope in %Prated/V
        k_p_slope_max = ts.param_value('eut_vw.k_p_v_max')
        k_p_slope_min = ts.param_value('eut_vw.k_p_v_min')
        hysteresis = ts.param_value('eut_vw.hysteresis')
        MSA_t = ts.param_value('eut_vw.MSA_t')
        v_stop_max = ts.param_value('eut_vw.vstop_max')
        v_stop_min = ts.param_value('eut_vw.vstop_min')
        # t_return is the delay before releasing the hysteresis power level and increasing power at k_power_rate
        t_return_max = ts.param_value('eut_vw.treturn_max')
        t_return_min = ts.param_value('eut_vw.treturn_min')
        # the k_power_rates are the hysteresis response in %P_rated/sec
        k_p_rate_max = ts.param_value('eut_vw.k_p_rate_max')
        k_p_rate_min = ts.param_value('eut_vw.k_p_rate_min')
        # Initiliaze VW test parameters variables
        curves = ts.param_value('vw.curves')
        pwr_lvl = ts.param_value('vw.power_lvl')
        n_iterations = ts.param_value('vw.n_iter')
        n_points = ts.param_value('vw.n_points')
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
        p_rated = ts.param_value('eut_vw.p_rated')
        eff = {
                1.0 : ts.param_value('eut_vw.efficiency_100')/100,
                0.33:ts.param_value('eut_vw.efficiency_33')/100,
        } 
        MSA_P_pct = (MSA_P/p_rated)*100.        
        pv.power_set(p_rated)
        pv.power_on()  # power on at p_rated
        # DAS soft channels
        das_points = {'sc': ('P_TARGET', 'P_TARGET_MIN', 'P_TARGET_MAX', 'P_ACT', 'V_TARGET','V_ACT','event')}

        # initialize data acquisition system
        daq = das.das_init(ts, sc_points=das_points['sc'])
        daq.sc['P_TARGET'] = 100
        daq.sc['P_TARGET_MIN'] = 100
        daq.sc['P_TARGET_MAX'] = 100
        daq.sc['V_TARGET'] = v_nom
        daq.sc['event'] = 'None'
        """
        EUT Configuration
        """ 
        # Configure the EUT communications
        if eut is not None:
            eut.config()
            ts.log_debug(eut.measurements())
            ts.log_debug('L/HVRT and trip parameters set to the widest range : v_min:{0} V, v_max:{1} V'.format(v_min,v_max))
            
            # TODO : Need to update FRT parameters with SunSpec Model reference
            eut.vrt_stay_connected_high(params={'Ena' : True,'ActCrv':0, 'Tms1':3000,'V1' : f_max,'Tms2':0.16,'V2' : v_max})
            eut.vrt_stay_connected_low(params={'Ena' : True,'ActCrv':0, 'Tms1':3000,'V1' : f_min,'Tms2':0.16,'V2' : v_min})
        else:
            ts.log_debug('Set L/HVRT and trip parameters to the widest range of adjustability possible.')
        """
        Test Configuration
        """
        if curves == 'Both':
            vw_curves = [1, 2]
        elif curves == 'Characteristic Curve 1':
            vw_curves = [1]
        else:  # Characteristic Curve 2
            vw_curves = [2]
        if pwr_lvl == 'All':
            pv_powers = [1., 0.33]
        elif pwr_lvl == '100%':
            pv_powers = [1.]
        else:  # Power at 33%
            pv_powers = [0.33]

        if hysteresis == 'Enabled':
            hyst = [True]
        elif hysteresis == 'Disabled':
            hyst = [False]
        else:
            hyst = [True, False]

        # open result summary file
        result_summary_filename = 'result_summary.csv'
        result_summary = open(ts.result_file_path(result_summary_filename), 'a+')
        ts.result_file(result_summary_filename)
        result_summary.write('Result,Test Name,Power Level,Iteration,direction,V_target,V_actual,Power_target,Power_actual,P_min,P_max,Dataset File\n')
        """
        Test start
        """
        for vw_curve in vw_curves:
            if vw_curve == 1:  # characteristic curve 1
                v_start = round(v_start_min,2)
                k_power_volt = k_p_slope_max
                v_stop = round(v_start + 100./k_power_volt,2)
            else:  # characteristic curve 2
                v_start = round(v_start_max,2)
                k_power_volt = k_p_slope_min
                v_stop = round(v_start + 100./k_power_volt,2)
            # Update hysteresis parameters
            for hys in hyst:
                if hys and vw_curve == 1:
                    k_power_rate = k_p_rate_max
                    v_stop = v_stop_min
                    t_return = t_return_max
                    # todo: hysteresis - der.py and der_sunspec.py and pysunspec must be updated
                    # vw_hyst_params = {'v_stop': v_stop, 'k_power_rate': k_power_rate, 't_return': t_return}
                    # eut.volt_watt(params=vw_hyst_params)
                if hys and vw_curve == 2:
                    k_power_rate = k_p_rate_max
                    v_stop = v_stop_min
                    t_return = t_return_max
                    # todo: hysteresis - der.py and der_sunspec.py and pysunspec must be updated
                    # vw_hyst_params = {'v_stop': v_stop, 'k_power_rate': k_power_rate, 't_return': t_return}
                    # eut.volt_watt(params=vw_hyst_params)
                # TODO : Need to update VW parameters with SunSpec Model reference
                if eut is not None:
                    vw_curve_params = {'v': [v_start, v_stop], 'w': [100., 0], 'DeptRef': 'W_MAX_PCT'}
                    vw_params = {'Ena': True, 'ActCrv': 1, 'curve': vw_curve_params}
                    eut.volt_watt(params=vw_params)
                    ts.log_debug('Initial EUT VW settings are %s' % eut.volt_watt())        
                ts.log('Testing VW function with these parameters : v_start = {0} V , v_stop = {1} V'.format(v_start,v_stop))
                for power in pv_powers:
                    # SA15.3.1(d) to (g)
                    v_stop = round(v_start + (power*100.)/k_power_volt,3) 
                    if v_stop < v_max:
                        v_points_up = list(np.linspace(v_min+MSA_V, v_start-1.5*MSA_V, n_points))                    
                        if MSA_V*1.5 >= round((v_stop - v_start)/2,2):
                            # This is to be sure to create a continuous v_step and handle in case msa is too high              
                            v_points_up += list(np.linspace(v_start+(MSA_V/(2*n_points)), v_stop-(MSA_V/(2*n_points)), 2*n_points))                            
                        else:
                            v_points_up += list(np.linspace(v_start+1.5*MSA_V, v_stop-1.5*MSA_V, 2*n_points))
                        v_points_up += list(np.linspace(v_stop+1.5*MSA_V, v_max-MSA_V, n_points))
                        v_steps_dic['up'] = np.around(v_points_up,  decimals=2)
                    else:  # slope extends past EUT trip point - skip final line segment
                        v_points_up = list(np.linspace(v_min+MSA_V, v_start-1.5*MSA_V, n_points)) + \
                                   list(np.linspace(v_start+1.5*MSA_V, v_max-MSA_V, 2*n_points))
                        v_steps_dic['up'] = np.around(v_points_up,  decimals=2)
                    # Include n_points on each of the 2 line segments for going down
                    v_points_down = list(reversed(v_points_up))
                    v_steps_dic['down'] = np.around(v_points_down,  decimals=2) 
                    if hys :
                        v_points_down = list(np.linspace(v_max-MSA_V, v_stop+1.5*MSA_V, 2*n_points)) + \
                                        list(np.linspace(v_stop-1.5*MSA_V, v_min+MSA_V, n_points))
                        v_steps_dic['down'] = np.around(v_points_down,  decimals=2)
                    ts.log('Testing VW function at the following voltage up points %s' % v_steps_dic['up'])
                    ts.log('Testing VW function at the following voltage down points %s' % v_steps_dic['down'])
                    
                    # Set to the power level and apply effiency correction
                    pv_power_setting = (p_rated*power)/eff[power]
                    ts.log('Set PV simulator power to {} with efficiency at {} %'.format(p_rated*power,eff[power]*100.))
                    pv.power_set(pv_power_setting)
                    for n_iter in range(n_iterations):                      
                        dataset_filename = 'VW_curve_%s_pwr_%0.2f_iter_%s.csv' % (vw_curve, power, n_iter + 1)                        
                        daq.data_capture(True)
                        for direction ,v_steps in v_steps_dic.iteritems():
                            for v_step in v_steps:                                
                                # The above seems to imply that t_return is active with and without hysteresis.
                                # Here we assume this is only enabled with hysteresis. Once reaching V_stop, the timer
                                # begins. After timing out, the EUT should ramp at K_power_rate back to rated power assuming
                                # V_stop < V_start
                                ts.log('        Recording power at voltage %0.3f V for 2*t_settling = %0.1f sec.' %
                                       (v_step, 2 * t_settling))
                                daq.sc['V_TARGET'] = v_step
                                daq.sc['event'] = 'v_step_{}'.format(direction)
                                
                                if hys and direction == 'down' and v_step > v_stop-1.5*MSA_V:
                                    p_targ = 0.
                                    p_min_to_pass = -MSA_P
                                    p_max_to_pass = MSA_P
                                    daq.sc['P_TARGET'] = p_targ
                                    daq.sc['P_TARGET_MIN'] = p_min_to_pass
                                    daq.sc['P_TARGET_MAX'] = p_max_to_pass
                                    grid.voltage(v_step)
                                    # If the t_return functionality is initiated when Vstop - MSAV < V < Vstop + MSAV
                                    # if v_step < v_stop-0.5*MSA_V:
                                    #     t_return_start_time = time.time()
                                    #     daq.sc['event'] = 't_return_Started'
                                    ts.sleep(2 * t_settling)
                                    daq.data_sample()
                                    data = daq.data_capture_read()
                                    AC_W_pct = (power_total(data, phases) / p_rated) * 100.
                                    passfail = 'Pass' if p_min_to_pass <= AC_W_pct <= p_max_to_pass else 'Fail'
                                    
                                elif hys and direction == 'down' and v_step == v_stop-1.5*MSA_V :
                                    p_targ = 100.  # it's ramping from 0 to 100
                                    p_min_to_pass = -MSA_P
                                    p_max_to_pass = MSA_P
                                    daq.sc['P_TARGET'] = p_targ
                                    daq.sc['P_TARGET_MIN'] = p_min_to_pass
                                    daq.sc['P_TARGET_MAX'] = p_max_to_pass
                                    grid.voltage(v_step)
                                    # Hysteresis Pass/Fail analysis
                                    hys_start_time = time.time()
                                    passfail = 'Pass'
                                    while (time.time() - hys_start_time) < (2*t_settling + t_return + 100./k_power_rate):
                                        hys_clock = (time.time() - hys_start_time)
                                        daq.data_sample()
                                        data = daq.data_capture_read()
                                        AC_W_pct = (power_total(data, phases) / p_rated) * 100.
                                        remaining_time = 2*t_settling + t_return + 100./k_power_rate - hys_clock
                                        ts.log('Evaluating hysteresis return. Time=%0.3f, P=%0.2f. t_return=%s s, '
                                               'rate=%0.2f. Additional analysis time = %0.2f' %
                                               (hys_clock, AC_W_pct, t_return, k_power_rate, remaining_time))
                                        # check to see that t_return is used correctly: fail if power increasing early
                                        if hys_clock < t_return + MSA_t and AC_W_pct > MSA_P:
                                            passfail = 'Fail'
                                        # check to see that the k_power_rate ramp is followed with appropriate MSA values
                                        if (t_return + MSA_t) < hys_clock < (t_return + MSA_t + 100./k_power_rate):
                                            # Verify the EUT is not under-performing or over-performing
                                            if (k_power_rate*(hys_clock-MSA_t) - MSA_P) < AC_W_pct or \
                                                            AC_W_pct < (k_power_rate*(hys_clock+MSA_t) + MSA_P):
                                                passfail = 'Fail'
                                    # Do final power check after settling time
                                    daq.data_sample()
                                    data = daq.data_capture_read()
                                    AC_W_pct = (power_total(data, phases) / p_rated) * 100.
                                    if AC_W_pct < p_min_to_pass or AC_W_pct > p_max_to_pass:
                                        passfail = 'Fail'

                                else:  # in all cases without hysteresis and when v_step < v_stop-1.5*MSA_V
                                    p_targ, p_min, p_max = p_msa_range(v_value=v_step,
                                                                       v_msa=MSA_V,
                                                                       p_msa_pct=MSA_P_pct,
                                                                       pwr = power,
                                                                       v_nom=v_nom,
                                                                       v_slope_start=v_start,
                                                                       v_slope_stop=v_stop)
                                    daq.sc['P_TARGET'] = p_targ
                                    daq.sc['P_TARGET_MIN'] = p_min
                                    daq.sc['P_TARGET_MAX'] = p_max
                                    ts.log('        Powers targ, min, max: %s, %s, %s' % (p_targ, p_min, p_max))
                                    grid.voltage(v_step)
                                    ts.sleep(2 * t_settling)
                                    daq.data_sample()
                                    data = daq.data_capture_read()                                    
                                    try:
                                        AC_W_pct = (power_total(data, phases) / p_rated) * 100.
                                        V_ACT = v_mean(data, phases)
                                        passfail = 'Pass' if daq.sc['P_TARGET_MIN'] <= AC_W_pct <= daq.sc['P_TARGET_MAX'] else 'Fail'
                                    except:
                                        AC_W_pct = 'No Data'
                                        V_ACT= 'No Data'
                                        passfail = 'Fail'
                                daq.sc['P_ACT'] = AC_W_pct
                                daq.sc['V_ACT'] = V_ACT
                                daq.sc['event'] = 'T_settling_done_{}'.format(direction)
                                daq.data_sample()
                                ts.log('        Powers measured: {} %'.format(AC_W_pct))
                                result_summary.write('%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n' %
                                                     (passfail,
                                                      ts.config_name(),
                                                      power*100.,
                                                      n_iter+1,
                                                      direction,
                                                      v_step,
                                                      V_ACT,
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
                        'plot.y.points': 'V_TARGET,V_ACT',
                        'plot.y.title': 'Voltage (V)',
                        'plot.V_TARGET.point': 'True',
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
                        ts.result_file(dataset_filename)

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
            if v_nom is not None:
                grid.voltage(v_nom)
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

        result = test_run()

        ts.result(result)
        if result == script.RESULT_FAIL:
            rc = 1

    except Exception, e:
        ts.log_error('Test script exception: %s' % traceback.format_exc())
        rc = 1

    sys.exit(rc)

info = script.ScriptInfo(name=os.path.basename(__file__), run=run, version='1.0.0')

# EUT VW parameters
info.param_group('eut_vw', label='VW EUT specified parameters',glob=True)
info.param('eut_vw.phases', label='Phases', default='Single Phase', values=['Single phase', 'Split phase', 'Three phase'])
info.param('eut_vw.p_rated', label='Output Power Rating (W)', default=10000.)
info.param('eut_vw.v_min', label='Min AC voltage range with function enabled (V)', default=108.)
info.param('eut_vw.v_max', label='Max AC voltage range with function enabled (V)', default=132.)
info.param('eut_vw.v_nom', label='Nominal AC voltage (V)', default=120.)
info.param('eut_vw.MSA_V', label='Manufacturer\'s stated AC voltage accuracy (V)', default=0.1)
info.param('eut_vw.MSA_P', label='Manufacturer\'s stated power accuracy (W)', default=10.)
info.param('eut_vw.ts', label='Settling time (s)', default=1.)
info.param('eut_vw.vstart_max', label='Max start of power reduction, v_start (V)', default=135.6)
info.param('eut_vw.vstart_min', label='Min start of power reduction, v_start (V)', default=123.6)
info.param('eut_vw.k_p_v_max', label='Maximum slope of active power reduction (%Prated/V)', default=100.)
info.param('eut_vw.k_p_v_min', label='Minimum slope of active power reduction (%Prated/V)', default=10.)
info.param('eut_vw.efficiency_33', label='CEC Efficiency list for power level = 33% at nominal VDC', default=97.0)
info.param('eut_vw.efficiency_100', label='CEC Efficiency list for power level = 100% at nominal VDC', default=96.9)
# Not including a 'Both' option because UL 1741 SA is either/or
info.param('eut_vw.hysteresis', label='Hysteresis in the Volt-Watt function', default='Disabled',values=['Enabled', 'Disabled'])  
info.param('eut_vw.MSA_t', label='Manufacturer\'s stated time accuracy (s)', default=0.01, active='vw.hysteresis',active_value=['Enabled'])
info.param('eut_vw.vstop_max', label='Max stop of voltage curtailment (V)', default=135.6, active='vw.hysteresis',active_value=['Enabled'])
info.param('eut_vw.vstop_min', label='Min stop of voltage curtailment (V)', default=123.6, active='vw.hysteresis',active_value=['Enabled'])
info.param('eut_vw.treturn_max', label='Maximum adjustment of a delay before return to normal operation (s)', default=0.1,active='vw.hysteresis', active_value=['Enabled'])
info.param('eut_vw.treturn_min', label='Minimum adjustment of a delay before return to normal operation (s)', default=0.1,active='vw.hysteresis', active_value=['Enabled'])
info.param('eut_vw.k_p_rate_max', label='Max active power rate to return to normal operation (%Prated/Sec)', default=0.1,active='vw.hysteresis', active_value=['Enabled'])
info.param('eut_vw.k_p_rate_min', label='Min active power rate to return to normal operation (%Prated/Sec)', default=0.1,active='vw.hysteresis', active_value=['Enabled'])
# VW test parameters
info.param_group('vw', label='Test Parameters')
info.param('vw.curves', label='Curves to Evaluate', default='Both', values=['Characteristic Curve 1','Characteristic Curve 2', 'Both'])
info.param('vw.power_lvl', label='Power Levels', default='All', values=['100%', '33%', 'All'])
info.param('vw.n_iter', label='Number of iteration for each test', default=3)
info.param('vw.n_points', label='Number of points tested above v_start', default=3)

# Other equipment parameters
der.params(info)
gridsim.params(info)
pvsim.params(info)
das.params(info)
hil.params(info)

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
