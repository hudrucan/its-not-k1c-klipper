# Load Cell Probe
#
# Copyright (C) 2025  Gareth Farrington <gareth@waves.ky>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import logging, math
import mcu
from . import probe, sos_filter, load_cell, hx71x, ads1220

np = None  # delay NumPy import until configuration time

# constants for fixed point numbers
Q2_INT_BITS = 2
Q2_FRAC_BITS = (32 - (1 + Q2_INT_BITS))
Q16_INT_BITS = 16
Q16_FRAC_BITS = (32 - (1 + Q16_INT_BITS))


######## Types

class TapClassifierModule(object):
    def classify(self, tap_analysis):
        pass


class NozzleCleanerModule(object):
    def clean_nozzle(self, attempt, retries, probe_pos):
        pass


# Capture and preserve a Trapezoidal Move as a python type
class TrapezoidalMove(object):
    def __init__(self, move):
        # copy c data to python memory
        self.print_time = float(move.print_time)
        self.move_t = float(move.move_t)
        self.start_v = float(move.start_v)
        self.accel = float(move.accel)
        self.start_x = float(move.start_x)
        self.start_y = float(move.start_y)
        self.start_z = float(move.start_z)
        self.x_r = float(move.x_r)
        self.y_r = float(move.y_r)
        self.z_r = float(move.z_r)

    def to_dict(self):
        return {
            'print_time': float(self.print_time), 'move_t': float(self.move_t),
            'start_v': float(self.start_v), 'accel': float(self.accel),
            'start_x': float(self.start_x), 'start_y': float(self.start_y),
            'start_z': float(self.start_z), 'x_r': float(self.x_r),
            'y_r': float(self.y_r), 'z_r': float(self.z_r)
        }


# point on a time/force graph
class ForcePoint(object):
    def __init__(self, time_t, force):
        self.time = float(time_t)
        self.force = float(force)

    def to_dict(self):
        return {'time': self.time, 'force': self.force}


# slope/intercept based line where x is time and y is force
class ForceLine(object):
    def __init__(self, slope, intercept):
        self.slope = float(slope)
        self.intercept = float(intercept)

    # measure angles between lines at the 25g = 0.1s (100ms) scale
    # Note: this is the same scale used by Prusa
    # returns +/- 0-180. Positive values represent clockwise rotation
    def angle(self, line, time_scale=0.1, gram_scale=25):
        scaling_factor = time_scale / gram_scale
        this_slope = self.slope * scaling_factor
        other_slope = line.slope * scaling_factor
        radians = (math.atan2(this_slope, 1) - math.atan2(other_slope, 1))
        return math.degrees(radians)

    def find_force(self, time):
        return self.slope * time + self.intercept

    def find_time(self, force):
        return (force - self.intercept) / self.slope

    def intersection(self, line):
        numerator = -self.intercept + line.intercept
        denominator = self.slope - line.slope
        # lines are parallel, will not intersect
        if denominator == 0.:
            # to get debuggable data we want to return a clearly bad value here
            return ForcePoint(0., 0.)
        intersection_time = numerator / denominator
        intersection_force = self.find_force(intersection_time)
        return ForcePoint(intersection_time, intersection_force)

    def to_dict(self):
        return {'slope': self.slope, 'intercept': self.intercept}


#########################
# Math Support Functions

# helper class for working with a time/force graph
# work with subsections to find elbows and best fit lines
class ForceGraph:
    def __init__(self, time_nd_64, force_nd_64):
        self.time = time_nd_64
        self.force = force_nd_64
        # linear implementation:
        self._cum_x = np.cumsum(self.time)
        self._cum_y = np.cumsum(self.force)
        self._cum_xx = np.cumsum(self.time * self.time)
        self._cum_xy = np.cumsum(self.time * self.force)
        self._cum_yy = np.cumsum(self.force * self.force)

    def _get_segment_sum(self, arr, start_idx, end_idx):
        prior_sum = 0 if start_idx == 0 else arr[start_idx - 1]
        return arr[end_idx] - prior_sum

    def _get_segment_stats(self, start_idx, end_idx):
        """
        Get statistics for segment [start_idx:end_idx] using cumulative sums
        """
        n = end_idx - start_idx
        sum_x = self._get_segment_sum(self._cum_x, start_idx, end_idx)
        sum_y = self._get_segment_sum(self._cum_y, start_idx, end_idx)
        sum_xx = self._get_segment_sum(self._cum_xx, start_idx, end_idx)
        sum_xy = self._get_segment_sum(self._cum_xy, start_idx, end_idx)
        sum_yy = self._get_segment_sum(self._cum_yy, start_idx, end_idx)
        return n, sum_x, sum_y, sum_xx, sum_xy, sum_yy

    def _least_squares(self, start_idx, end_idx):
        """
        Compute slope/intercept and RSS for a segment of the data
        """
        n = (end_idx - start_idx) + 1
        if n < 2:
            raise ValueError("Error: fewer than 2 points used")

        sum_x = self._get_segment_sum(self._cum_x, start_idx, end_idx)
        sum_y = self._get_segment_sum(self._cum_y, start_idx, end_idx)
        sum_xx = self._get_segment_sum(self._cum_xx, start_idx, end_idx)
        sum_xy = self._get_segment_sum(self._cum_xy, start_idx, end_idx)
        sum_yy = self._get_segment_sum(self._cum_yy, start_idx, end_idx)

        denom = n * sum_xx - sum_x * sum_x
        if abs(denom) < 1e-10:
            return None, np.inf

        slope = (n * sum_xy - sum_x * sum_y) / denom
        intercept = (sum_y - slope * sum_x) / n
        rss = (sum_yy
               - 2 * slope * sum_xy
               - 2 * intercept * sum_y
               + slope * slope * sum_xx
               + 2 * slope * intercept * sum_x
               + n * intercept * intercept)
        return [slope, intercept], max(0.0, rss)

    # search exhaustively for the 2 lines that best fit the data
    # return the elbow index
    def _two_lines_best_fit(self, start_idx, end_idx):
        best_error = float('inf')
        best_fit_index = -1
        for i in range(1 + start_idx, end_idx - 1):
            params1, r1 = self._least_squares(start_idx, i)
            params2, r2 = self._least_squares(i + 1, end_idx)
            if params1 is not None and params2 is not None:
                error = r1 + r2
                if error < best_error:
                    best_error = error
                    best_fit_index = i
        # the index returns is the first point in the second line
        return best_fit_index

    def find_elbow(self, start_idx, end_idx):
        return self._two_lines_best_fit(start_idx, end_idx)

    # finds the index nearest to a time
    def index_near(self, instant):
        idx = int(np.searchsorted(self.time, instant))
        return min(idx, len(self.time) - 1)

    # construct a line from 2 points
    def _points_to_line(self, a, b):
        slope = (b.force - a.force) / (b.time - a.time)
        intercept = a.force - (slope * a.time)
        return ForceLine(slope, intercept)

    # construct a line using a subset of the graph
    def line(self, start_idx, end_idx):
        params, rss = self._least_squares(start_idx, end_idx)
        return ForceLine(params[0], params[1])

    # given a line and a range, calculate the standard deviation of the noise
    def noise_std(self, start_idx, end_idx, line):
        f = self.force[start_idx:end_idx]
        t = self.time[start_idx:end_idx]
        noise = []
        for i in range(len(f)):
            noise.append(f[i] - line.find_force(t[i]))
        return np.std(noise, dtype=np.float64)

    # true if the reference force won't be confused for noise in the graph chunk
    # reference force must be more than 3 standard deviations away from the line
    # at the reference index
    def is_clear_signal(self, start_idx, end_idx, line, reference_idx,
            force_idx):
        noise = self.noise_std(start_idx, end_idx, line)
        noise_3_std = noise * 3
        base_force = line.find_force(self.time[reference_idx])
        return abs(base_force - self.force[force_idx]) > noise_3_std

    # return the first index that exceeds the median force between
    # start and end index
    def _split_by_force(self, start_idx, end_idx):
        start_f = self.force[start_idx]
        end_f = self.force[end_idx]
        median_f = (start_f + end_f) / 2.
        scan = range(start_idx, end_idx)
        # if force is ascending, swap the scan direction
        if start_f > end_f:
            scan = reversed(scan)
        for i in scan:
            if self.force[i] > median_f:
                return i
        return None

    # break a tap event down into 6 points and 5 lines:
    #           |*----*\
    #           |       \
    #    *-----*|        \*-----*
    def tap_decompose(self, homing_end_time, pullback_start_time,
            pullback_cruise_time, pullback_cruise_duration):
        homing_end_idx = self.index_near(homing_end_time)
        # use the pullback duration to trim the amount of approach data used
        homing_start_time = homing_end_time - pullback_cruise_duration
        homing_start_idx = self.index_near(homing_start_time)
        # locate the point where the probe made contact with the bed
        contact_elbow_idx = self.find_elbow(homing_start_idx, homing_end_idx)

        pullback_start_idx = self.index_near(pullback_start_time)
        pullback_cruise_idx = self.index_near(pullback_cruise_time)
        # limit use of additional data after the pullback move ends
        pullback_end_time = (
                pullback_cruise_time + (pullback_cruise_duration * 1.5))
        pullback_end_idx = self.index_near(pullback_end_time)

        # l1 is the approach line
        l1 = self.line(homing_start_idx, contact_elbow_idx)
        # sometime after contact_elbow_idx is the peak force and the start of
        # the dwell line
        dwell_end = self.time[contact_elbow_idx] + pullback_cruise_duration
        dwell_end_idx = min(pullback_start_idx, self.index_near(dwell_end))
        dwell_start_idx = self.find_elbow(contact_elbow_idx, dwell_end_idx)
        # l2 is the compression line
        # also +1 the last index in case its sequential [1, 2]
        l2 = self.line(contact_elbow_idx, dwell_start_idx + 1)
        # l3 is the dwell line
        l3 = self.line(dwell_start_idx, pullback_start_idx)

        # find the approximate elbow location
        break_contact_idx = self.find_elbow(pullback_cruise_idx,
            pullback_end_idx)
        # l5 is the line after decompression ends
        l5 = self.line(break_contact_idx, pullback_end_idx)
        # split the points between the elbow and the start of movement by force
        midpoint_idx = self._split_by_force(pullback_cruise_idx,
            break_contact_idx)
        # elbow finding success depends on their being good signal-to-noise
        # this checks if there will be enough clear data to analyze
        use_curve_optimization = False
        if midpoint_idx is not None:
            clear_dwell = self.is_clear_signal(dwell_start_idx, dwell_end_idx,
                l3, dwell_end_idx, midpoint_idx - 1)
            clear_decomp = self.is_clear_signal(break_contact_idx,
                pullback_end_idx, l5, break_contact_idx, midpoint_idx)
            use_curve_optimization = clear_dwell and clear_decomp
        if use_curve_optimization:
            # perform iterative refinement
            l4_start = self.line(pullback_cruise_idx, midpoint_idx)
            # real break contact index
            break_contact_idx = self.find_elbow(midpoint_idx, pullback_end_idx)
            l4_end = self.line(midpoint_idx, break_contact_idx)
            l5 = self.line(break_contact_idx, pullback_end_idx)
            # a synthetic l4 is built from 2 points:
            l4 = self._points_to_line(l4_start.intersection(l3),
                l4_end.intersection(l5))
        else:
            # noise is too high, don't use the curve optimization
            l4 = self.line(pullback_cruise_idx, break_contact_idx)
            # log for user debugging
            logging.info('TapAnalysis: curve optimization not used')

        return [l1, l2, l3, l4, l5], homing_start_idx, pullback_end_idx


class TapValidationError(Exception):
    def __init__(self, error_code, message):
        super(TapValidationError, self).__init__(message)
        self.error_code = error_code
        pass

    def to_dict(self):
        return  {
            'error_code': self.error_code,
            'message': str(self)
        }


# Move index constants. The PROBE_START move may be deleted from the trapq if
# the probe takes longer than 30s. Indexing from the end of the list
# is always consistent:
PROBE_START = -6
PROBE_CRUISE = -5
PROBE_HALT = -4
PULLBACK_START = -3
PULLBACK_CRUISE = -2
PULLBACK_END = -1


class TapAnalysis(object):
    def __init__(self, samples, trigger_force):
        self._is_valid = False
        self._tap_pos = None
        self._tap_points = []
        self._tap_lines = []
        self._tap_angles = []
        self._elapsed = 0.
        self._collection_time = 0.
        self._error = None
        self._home_end_time = None
        self._pullback_start_time = None
        self._pullback_end_time = None
        self._pullback_cruise_time = None
        self._pullback_duration = None
        self._homing_start_index = 0
        self._pullback_end_index = -1
        nd_samples = np.asarray(samples, dtype=np.float64)
        self.time = nd_samples[:, 0]
        self.force = nd_samples[:, 1]

    def get_collection_time(self):
        return self._collection_time

    def set_collection_time(self, collection_time):
        self._collection_time = collection_time

    # convert to dictionary for JSON encoder
    def to_dict(self):
        return {
            'time': self._time.tolist(),
            'force': self._force.tolist(),
            'tap_points': [point.to_dict() for point in self._tap_points],
            'tap_lines': [line.to_dict() for line in self._tap_lines],
            'tap_angles': self.get_tap_angles(),
            'tap_pos': self.get_tap_pos(),
            'moves': [move.to_dict() for move in self._moves],
            'home_end_time': self.get_home_end_time(),
            'pullback_start_time': self.get_pullback_start_time(),
            'pullback_end_time': self.get_pullback_end_time(),
            'collection_time': self.get_collection_time(),
            'elapsed': self.get_elapsed(),
            'is_valid': self.is_valid(),
            'error': None if self._error is None else self._error.to_dict(),
        }


# Orchestrate TapAnalysis and TapClassifier. Handle timing, error capture,
# event broadcast, clients & logging
class TapAnalysisHelper:
    def __init__(self, printer, name, tap_classifier):
        self._printer = printer
        self._tap_classifier = tap_classifier
        # webhooks support
        self._clients = load_cell.ApiClientHelper(printer)
        header = {"header": ["probe_tap_event"]}
        self._clients.add_mux_endpoint("load_cell_probe/dump_taps",
            "load_cell_probe", name, header)

    def analyze(self, samples, trigger_force, collection_time):
        t_start = time.time()
        tap_analysis = TapAnalysis(samples, trigger_force)
        try:
            tap_analysis.analyze(self._printer)
        except TapValidationError as ve:
            tap_analysis.set_is_valid(False)
            tap_analysis.set_validation_error(ve)
        # tap classifier always gets to process the data
        try:
            self._tap_classifier.classify(tap_analysis)
        except TapValidationError as ve:
            tap_analysis.set_is_valid(False)
            tap_analysis.set_validation_error(ve)
        # total elapsed time for all calculations
        tap_analysis.set_elapsed(time.time() - t_start)
        tap_analysis.set_collection_time(collection_time)
        # broadcast tap event data:
        self._clients.send({'tap': tap_analysis.to_dict()})
        self._log_errors(tap_analysis)
        return tap_analysis

    # log errors to event log
    def _log_errors(self, tap_analysis):
        # if the tap is valid, don't log any errors
        if tap_analysis.is_valid():
            return
        # log errors
        ve = tap_analysis.get_validation_error()
        logging.info("Bad tap detected: %s - %s" % (ve.error_code, ve))

    # get internal tap events
    def add_client(self, callback):
        self._clients.add_client(callback)


# Access a parameter from config or GCode command via a consistent interface
# stores name and constraints to keep things DRY
class ParamHelper:
    def __init__(self, config, name, type_name, default=None, minval=None,
            maxval=None, above=None, below=None, max_len=None):
        self._config_section = config.get_name()
        self._config_error = config.error
        self.name = name
        self._type_name = type_name
        self.value = default
        self.minval = minval
        self.maxval = maxval
        self.above = above
        self.below = below
        self.max_len = max_len
        # read from config once
        self.value = self.get(config=config)

    def _get_name(self, gcmd):
        return self.name.upper() if gcmd else self.name

    def _validate_float(self, description, error, value, above, below):
        above = above or self.above
        if above is not None and value <= above:
            raise error("%s must be above %s" % (description, above))
        below = below or self.below
        if below is not None and value >= below:
            raise error("%s must be below %s" % (description, below))

    # support for validating individual options in a list of floats
    def _validate_float_list(self, gcmd, values, above, below):
        if gcmd:
            description = ("Error on '%s': %s" % (
                gcmd.get_commandline(), self._get_name(gcmd)))
            error = gcmd.error
        else:
            description = ("Option '%s' in section '%s'" % (
                self._get_name(gcmd), self._config_section))
            error = self._config_error
        if self.max_len is not None and len(values) > self.max_len:
            raise error(
                "%s has maximum length %s" % (description, self.max_len))
        for value in values:
            self._validate_float(description, error, value, above, below)

    def _get_int(self, config, gcmd, minval, maxval):
        get = gcmd.get_int if gcmd else config.getint
        return get(self._get_name(gcmd), self.value, minval or self.minval,
            maxval or self.maxval)

    def _get_float(self, config, gcmd, minval, maxval, above, below):
        get = gcmd.get_float if gcmd else config.getfloat
        return get(self._get_name(gcmd), self.value, minval or self.minval,
            maxval or self.maxval, above or self.above, below or self.below)

    def _get_float_list(self, config, gcmd, above, below):
        # this code defaults to the empty list, never return None
        default = (self.value or [])
        if gcmd:
            # if the parameter isn't part of the command, return the default
            if not self._get_name(gcmd) in gcmd.get_command_parameters():
                return default
            # parameter exists, always prefer whatever is in the command
            value = gcmd.get(self._get_name(gcmd), default='')
            # Return an empty list for empty value
            if len(value.strip()) == 0:
                return []
            try:
                float_list = [float(p.strip()) for p in value.split(',')]
            except:
                raise gcmd.error("Error on '%s': unable to parse %s" % (
                    gcmd.get_commandline(), value))
        else:
            float_list = config.getfloatlist(self._get_name(gcmd),
                default=default)
        if float_list:
            self._validate_float_list(gcmd, float_list, above, below)
        return float_list

    def get(self, gcmd=None, minval=None, maxval=None, above=None, below=None,
            config=None):
        if config is None and gcmd is None:
            return self.value
        if self._type_name == 'int':
            return self._get_int(config, gcmd, minval, maxval)
        elif self._type_name == 'float':
            return self._get_float(config, gcmd, minval, maxval, above, below)
        else:
            return self._get_float_list(config, gcmd, above, below)


def intParamHelper(config, name, default=None, minval=None, maxval=None):
    return ParamHelper(config, name, 'int', default, minval=minval,
        maxval=maxval)


def floatParamHelper(config, name, default=None, minval=None, maxval=None,
        above=None, below=None):
    return ParamHelper(config, name, 'float', default, minval=minval,
        maxval=maxval, above=above, below=below)


def floatListParamHelper(config, name, default=None, above=None, below=None,
        max_len=None):
    return ParamHelper(config, name, 'float_list', default, above=above,
        below=below, max_len=max_len)


# container for filter parameters
# allows different filter configurations to be compared
class ContinuousTareFilter:
    def __init__(self, sps=None, drift=None, drift_delay=None, buzz=None,
            buzz_delay=None, notches=None, notch_quality=None):
        self.sps = sps
        self.drift = drift
        self.drift_delay = drift_delay
        self.buzz = buzz
        self.buzz_delay = buzz_delay
        self.notches = notches
        self.notch_quality = notch_quality

    def __eq__(self, other):
        if not isinstance(other, ContinuousTareFilter):
            return False
        return (
                self.sps == other.sps and self.drift == other.drift and
                self.drift_delay == other.drift_delay and self.buzz ==
                other.buzz and self.buzz_delay == other.buzz_delay and
                self.notches == other.notches and self.notch_quality ==
                other.notch_quality)

    # create a filter design from the parameters
    def design_filter(self, error_func):
        design = sos_filter.DigitalFilter(self.sps, error_func, self.drift,
            self.drift_delay, self.buzz, self.buzz_delay, self.notches,
            self.notch_quality)
        fixed_filter = sos_filter.FixedPointSosFilter(
            design.get_filter_sections(), design.get_initial_state(),
            Q2_INT_BITS, Q16_INT_BITS)
        return fixed_filter


# Combine ContinuousTareFilter and SosFilter into an easy-to-use class
class ContinuousTareFilterHelper:
    def __init__(self, config, sensor, cmd_queue):
        self._sensor = sensor
        self._sps = self._sensor.get_samples_per_second()
        max_filter_frequency = math.floor(self._sps / 2.)
        # setup filter parameters
        self._drift_param = floatParamHelper(config,
            "drift_filter_cutoff_frequency", default=None, minval=0.1,
            maxval=20.0)
        self._drift_delay_param = intParamHelper(config, "drift_filter_delay",
            default=2, minval=1, maxval=2)
        self._buzz_param = floatParamHelper(config,
            "buzz_filter_cutoff_frequency", default=None,
            above=min(80.0, max_filter_frequency - 1.0),
            below=max_filter_frequency)
        self._buzz_delay_param = intParamHelper(config, "buzz_filter_delay",
            default=2, minval=1, maxval=2)
        self._notches_param = floatListParamHelper(config,
            "notch_filter_frequencies", default=[], above=0.,
            below=max_filter_frequency, max_len=2)
        self._notch_quality_param = floatParamHelper(config,
            "notch_filter_quality", default=2.0, minval=0.5, maxval=6.0)
        # filter design specified in the config file, used for defaults
        self._config_design = ContinuousTareFilter()  # empty filter
        self._config_design = self._build_filter()
        # filter design currently inside the MCU
        self._active_design = self._config_design
        self._sos_filter = self._create_filter(
            self._active_design.design_filter(config.error), cmd_queue)

    def _build_filter(self, gcmd=None):
        drift = self._drift_param.get(gcmd)
        drift_delay = self._drift_delay_param.get(gcmd)
        buzz = self._buzz_param.get(gcmd)
        buzz_delay = self._buzz_delay_param.get(gcmd)
        # notches must be between drift and buzz:
        notches = self._notches_param.get(gcmd, above=drift, below=buzz)
        notch_quality = self._notch_quality_param.get(gcmd)
        return ContinuousTareFilter(self._sps, drift, drift_delay, buzz,
            buzz_delay, notches, notch_quality)

    def _create_filter(self, fixed_filter, cmd_queue):
        return sos_filter.SosFilter(self._sensor.get_mcu(), cmd_queue,
            fixed_filter, 4)

    def update_from_command(self, gcmd, cq=None):
        gcmd_filter = self._build_filter(gcmd)
        # if filters are identical, no change required
        if self._active_design == gcmd_filter:
            return
        # update MCU filter from GCode command
        self._sos_filter.change_filter(
            self._active_design.design_filter(gcmd.error))

    def get_sos_filter(self):
        return self._sos_filter


# check results from the collector for errors and raise an exception is found
def check_sensor_errors(results, printer):
    samples, errors = results
    if errors:
        raise printer.command_error("Load cell sensor reported errors while"
                                    " probing: %i errors, %i overflows" % (
                                        errors[0], errors[1]))
    return samples


class LoadCellProbeConfigHelper:
    def __init__(self, config, load_cell_inst):
        self._printer = config.get_printer()
        self._load_cell = load_cell_inst
        self._sensor = load_cell_inst.get_sensor()
        self._rest_time = 1. / float(self._sensor.get_samples_per_second())
        # Collect 4 x 60hz power cycles of data to average across power noise
        self._tare_time_param = floatParamHelper(config, 'tare_time',
            default=4. / 60., minval=0.01, maxval=1.0)
        # triggering options
        self._trigger_force_param = intParamHelper(config, 'trigger_force',
            default=75, minval=10, maxval=250)
        self._force_safety_limit_param = intParamHelper(config,
            'force_safety_limit', minval=100, maxval=5000, default=2000)

    def get_tare_samples(self, gcmd=None):
        tare_time = self._tare_time_param.get(gcmd)
        sps = self._sensor.get_samples_per_second()
        return max(2, math.ceil(tare_time * sps))

    def get_trigger_force_grams(self, gcmd=None):
        return self._trigger_force_param.get(gcmd)

    def get_safety_limit_grams(self, gcmd=None):
        return self._force_safety_limit_param.get(gcmd)

    def get_rest_time(self):
        return self._rest_time

    def get_safety_range(self, gcmd=None):
        counts_per_gram = self._load_cell.get_counts_per_gram()
        # calculate the safety band
        zero = self._load_cell.get_reference_tare_counts()
        safety_counts = int(counts_per_gram * self.get_safety_limit_grams(gcmd))
        safety_min = int(zero - safety_counts)
        safety_max = int(zero + safety_counts)
        # don't allow a safety range outside the sensor's real range
        sensor_min, sensor_max = self._load_cell.get_sensor().get_range()
        if safety_min <= sensor_min or safety_max >= sensor_max:
            cmd_err = self._printer.command_error
            raise cmd_err("Load cell force_safety_limit exceeds sensor range!")
        return safety_min, safety_max

    # calculate 1/counts_per_gram in Q2 fixed point
    def get_grams_per_count(self):
        counts_per_gram = self._load_cell.get_counts_per_gram()
        # The counts_per_gram could be so large that it becomes 0.0 when
        # converted to Q2 format. This would mean the ADC range only measures a
        # few grams which seems very unlikely. Treat this as an error:
        if counts_per_gram >= 2**Q2_FRAC_BITS:
            raise OverflowError("counts_per_gram value is too large to filter")
        return sos_filter.to_fixed_32((1. / counts_per_gram), Q2_INT_BITS)


# McuLoadCellProbe is the interface to `load_cell_probe` on the MCU
# This also manages the SosFilter so all commands use one command queue
class McuLoadCellProbe:
    WATCHDOG_MAX = 3
    ERROR_SAFETY_RANGE = mcu.MCU_trsync.REASON_COMMS_TIMEOUT + 1
    ERROR_OVERFLOW = mcu.MCU_trsync.REASON_COMMS_TIMEOUT + 2
    ERROR_WATCHDOG = mcu.MCU_trsync.REASON_COMMS_TIMEOUT + 3

    def __init__(self, config, load_cell_inst, sos_filter_inst, config_helper,
            trigger_dispatch):
        self._printer = config.get_printer()
        self._load_cell = load_cell_inst
        self._sos_filter = sos_filter_inst
        self._config_helper = config_helper
        self._sensor = load_cell_inst.get_sensor()
        self._mcu = self._sensor.get_mcu()
        # configure MCU objects
        self._dispatch = trigger_dispatch
        self._cmd_queue = self._dispatch.get_command_queue()
        self._oid = self._mcu.create_oid()
        self._config_commands()
        self._home_cmd = None
        self._query_cmd = None
        self._set_range_cmd = None
        self._mcu.register_config_callback(self._build_config)
        self._printer.register_event_handler("klippy:connect", self._on_connect)

    def _config_commands(self):
        self._sos_filter.create_filter()
        self._mcu.add_config_cmd(
            "config_load_cell_probe oid=%d sos_filter_oid=%d" % (
                self._oid, self._sos_filter.get_oid()))

    def _build_config(self):
        # Lookup commands
        self._query_cmd = self._mcu.lookup_query_command(
            "load_cell_probe_query_state oid=%c",
            "load_cell_probe_state oid=%c is_homing_trigger=%c "
            "trigger_ticks=%u", oid=self._oid, cq=self._cmd_queue)
        self._set_range_cmd = self._mcu.lookup_command(
            "load_cell_probe_set_range"
            " oid=%c safety_counts_min=%i safety_counts_max=%i tare_counts=%i"
            " trigger_grams=%u grams_per_count=%i", cq=self._cmd_queue)
        self._home_cmd = self._mcu.lookup_command(
            "load_cell_probe_home oid=%c trsync_oid=%c trigger_reason=%c"
            " error_reason=%c clock=%u rest_ticks=%u timeout=%u",
            cq=self._cmd_queue)

    # the sensor data stream is connected on the MCU at the ready event
    def _on_connect(self):
        self._sensor.attach_load_cell_probe(self._oid)

    def get_oid(self):
        return self._oid

    def get_mcu(self):
        return self._mcu

    def get_load_cell(self):
        return self._load_cell

    def get_dispatch(self):
        return self._dispatch

    def set_endstop_range(self, tare_counts, gcmd=None):
        # update the load cell so it reflects the new tare value
        self._load_cell.tare(tare_counts)
        # update internal tare value
        safety_min, safety_max = self._config_helper.get_safety_range(gcmd)
        args = [self._oid, safety_min, safety_max, int(tare_counts),
            self._config_helper.get_trigger_force_grams(gcmd),
            self._config_helper.get_grams_per_count()]
        self._set_range_cmd.send(args)
        self._sos_filter.reset_filter()

    def home_start(self, print_time):
        clock = self._mcu.print_time_to_clock(print_time)
        rest_time = self._config_helper.get_rest_time()
        rest_ticks = self._mcu.seconds_to_clock(rest_time)
        self._home_cmd.send([self._oid, self._dispatch.get_oid(),
            mcu.MCU_trsync.REASON_ENDSTOP_HIT, self.ERROR_SAFETY_RANGE, clock,
            rest_ticks, self.WATCHDOG_MAX], reqclock=clock)

    def clear_home(self):
        params = self._query_cmd.send([self._oid])
        # The time of the first sample that triggered is in "trigger_ticks"
        trigger_ticks = self._mcu.clock32_to_clock64(params['trigger_ticks'])
        # clear trsync from load_cell_endstop
        self._home_cmd.send([self._oid, 0, 0, 0, 0, 0, 0, 0])
        return self._mcu.clock_to_print_time(trigger_ticks)


# Execute probing moves using the McuLoadCellProbe
class LoadCellProbingMove:
    ERROR_MAP = {
        mcu.MCU_trsync.REASON_COMMS_TIMEOUT: "Communication timeout during "
                                             "homing",
        McuLoadCellProbe.ERROR_SAFETY_RANGE: "Load Cell Probe Error: load "
                                             "exceeds safety limit",
        McuLoadCellProbe.ERROR_OVERFLOW: "Load Cell Probe Error: fixed point "
                                         "math overflow",
        McuLoadCellProbe.ERROR_WATCHDOG: "Load Cell Probe Error: timed out "
                                         "waiting for sensor data"
    }

    def __init__(self, config, mcu_load_cell_probe, param_helper,
            continuous_tare_filter_helper, config_helper):
        self._printer = config.get_printer()
        self._mcu_load_cell_probe = mcu_load_cell_probe
        self._param_helper = param_helper
        self._continuous_tare_filter_helper = continuous_tare_filter_helper
        self._config_helper = config_helper
        self._mcu = mcu_load_cell_probe.get_mcu()
        self._load_cell = mcu_load_cell_probe.get_load_cell()
        self._z_min_position = probe.lookup_minimum_z(config)
        self._dispatch = mcu_load_cell_probe.get_dispatch()
        probe.LookupZSteppers(config, self._dispatch.add_stepper)
        # internal state tracking
        self._tare_counts = 0
        self._last_trigger_time = 0

    def _start_collector(self):
        toolhead = self._printer.lookup_object('toolhead')
        # homing uses the toolhead last move time which gets special handling
        # to significantly buffer print_time if the move queue has drained
        print_time = toolhead.get_last_move_time()
        collector = self._load_cell.get_collector()
        collector.start_collecting(min_time=print_time)
        return collector

    # pauses for the last move to complete and then
    # sets the endstop tare value and range
    def _pause_and_tare(self, gcmd):
        collector = self._start_collector()
        num_samples = self._config_helper.get_tare_samples(gcmd)
        # use collect_min collected samples are not wasted
        results = collector.collect_min(num_samples)
        tare_samples = check_sensor_errors(results, self._printer)
        tare_counts = np.average(np.array(tare_samples)[:, 2].astype(float))
        # update sos_filter with any gcode parameter changes
        self._continuous_tare_filter_helper.update_from_command(gcmd)
        self._mcu_load_cell_probe.set_endstop_range(tare_counts, gcmd)

    def _home_start(self, print_time):
        # start trsync
        trigger_completion = self._dispatch.start(print_time)
        self._mcu_load_cell_probe.home_start(print_time)
        return trigger_completion

    def home_start(self, print_time, sample_time, sample_count, rest_time,
            triggered=True):
        return self._home_start(print_time)

    def home_wait(self, home_end_time):
        self._dispatch.wait_end(home_end_time)
        # trigger has happened, now to find out why...
        res = self._dispatch.stop()
        # clear the homing state so it stops processing samples
        self._last_trigger_time = self._mcu_load_cell_probe.clear_home()
        if res >= mcu.MCU_trsync.REASON_COMMS_TIMEOUT:
            error = "Load Cell Probe Error: unknown reason code %i" % (res,)
            if res in self.ERROR_MAP:
                error = self.ERROR_MAP[res]
            raise self._printer.command_error(error)
        if res != mcu.MCU_trsync.REASON_ENDSTOP_HIT:
            return 0.
        return self._last_trigger_time

    def get_steppers(self):
        return self._dispatch.get_steppers()

    # Probe towards z_min until the load_cell_probe on the MCU triggers
    def probing_move(self, gcmd):
        # do not permit probing if the load cell is not calibrated
        if not self._load_cell.is_calibrated():
            raise self._printer.command_error("Load Cell not calibrated")
        # tare the sensor just before probing
        self._pause_and_tare(gcmd)
        # get params for the homing move
        toolhead = self._printer.lookup_object('toolhead')
        pos = toolhead.get_position()
        pos[2] = self._z_min_position
        speed = self._param_helper.get_probe_params(gcmd)['probe_speed']
        phoming = self._printer.lookup_object('homing')
        # start collector after tare samples are consumed
        collector = self._start_collector()
        # do homing move
        return phoming.probing_move(self, pos, speed), collector

    # Wait for the MCU to trigger with no movement
    def probing_test(self, gcmd, timeout):
        self._pause_and_tare(gcmd)
        toolhead = self._printer.lookup_object('toolhead')
        print_time = toolhead.get_last_move_time()
        self._home_start(print_time)
        return self.home_wait(print_time + timeout)

    def get_status(self, eventtime):
        return {
            'tare_counts': self._tare_counts,
            'last_trigger_time': self._last_trigger_time,
        }


# Perform a single complete tap
class TappingMove:
    def __init__(self, config, mcu, load_cell_probing_move, tap_analysis_helper,
            config_helper):
        self._printer = config.get_printer()
        self._mcu = mcu
        self._load_cell_probing_move = load_cell_probing_move
        self._config_helper = config_helper
        # track results of the last tap
        self._last_result = None
        self._is_last_result_valid = False
        # webhooks support
        self._clients = load_cell.ApiClientHelper(config.get_printer())
        name = config.get_name()
        header = {"header": ["probe_tap_event"]}
        self._clients.add_mux_endpoint("load_cell_probe/dump_taps",
            "load_cell_probe", name, header)

    # perform a probing move and a pullback move
    def run_tap(self, gcmd):
        # do the descending move
        epos, collector = self._load_cell_probing_move.probing_move(gcmd)
        # collect samples from the tap
        results = collector.collect_until(pullback_end_time)
        # calculate how long we waited to get the data
        t_end = self._printer.get_reactor().monotonic()
        t_end = self._mcu.estimated_print_time(t_end)
        logging.error(f"t_start {pullback_end_time}, t_end: {t_end}")
        collection_time = t_end - pullback_end_time
        # check for data errors
        samples = check_sensor_errors(results, self._printer)
        # Analyze the tap data
        tap_analysis = self._tap_analysis_helper.analyze(samples,
            trigger_force, collection_time)
        self._last_analysis = tap_analysis
        self._is_last_result_valid = tap_analysis.is_valid()
        if self._is_last_result_valid:
            epos[2] = tap_analysis.get_tap_pos()[2]
        self._last_result = epos[2]
        return epos, self._is_last_result_valid

    def get_status(self, eventtime):
        return {
            'last_z_result': self._last_result,
            'is_last_tap_valid': self._is_last_result_valid
        }


# ProbeSession that implements Tap logic
class TapSession:
    def __init__(self, config, tapping_move, probe_params_helper):
        self._printer = config.get_printer()
        self._tapping_move = tapping_move
        self._probe_params_helper = probe_params_helper
        # Session state
        self._results = []

    def start_probe_session(self, gcmd):
        return self

    def end_probe_session(self):
        self._results = []

    # probe until a single good sample is returned or retries are exhausted
    def run_probe(self, gcmd):
        epos, is_good = self._tapping_move.run_tap(gcmd)
        self._results.append(epos)

    def pull_probed_results(self):
        res = self._results
        self._results = []
        return res


class LoadCellProbeCommands:
    def __init__(self, config, load_cell_probing_move):
        self._printer = config.get_printer()
        self._load_cell_probing_move = load_cell_probing_move
        self._register_commands()

    def _register_commands(self):
        # Register commands
        gcode = self._printer.lookup_object('gcode')
        gcode.register_command("LOAD_CELL_TEST_TAP",
            self.cmd_LOAD_CELL_TEST_TAP, desc=self.cmd_LOAD_CELL_TEST_TAP_help)

    cmd_LOAD_CELL_TEST_TAP_help = "Tap the load cell probe to verify operation"

    def cmd_LOAD_CELL_TEST_TAP(self, gcmd):
        taps = gcmd.get_int("TAPS", 3, minval=1, maxval=10)
        timeout = gcmd.get_float("TIMEOUT", 30., minval=1., maxval=120.)
        gcmd.respond_info("Tap the load cell %s times:" % (taps,))
        reactor = self._printer.get_reactor()
        for i in range(0, taps):
            result = self._load_cell_probing_move.probing_test(gcmd, timeout)
            if result == 0.:
                # notify of error, likely due to timeout
                raise gcmd.error("Test timeout out")
            gcmd.respond_info("Tap Detected!")
            # give the user some time for their finger to move away
            reactor.pause(reactor.monotonic() + 0.2)
        gcmd.respond_info("Test complete, %s taps detected" % (taps,))


class LoadCellPrinterProbe:
    def __init__(self, config):
        cfg_error = config.error
        try:
            global np
            import numpy as np
        except:
            raise cfg_error("[load_cell_probe] requires the NumPy module")
        self._printer = config.get_printer()
        # Sensor types supported by load_cell_probe
        sensors = {}
        sensors.update(hx71x.HX71X_SENSOR_TYPES)
        sensors.update(ads1220.ADS1220_SENSOR_TYPE)
        sensor_class = config.getchoice('sensor_type', sensors)
        sensor = sensor_class(config)
        self._load_cell = load_cell.LoadCell(config, sensor)
        # Read all user configuration and build modules
        config_helper = LoadCellProbeConfigHelper(config, self._load_cell)
        self._mcu = self._load_cell.get_sensor().get_mcu()
        trigger_dispatch = mcu.TriggerDispatch(self._mcu)
        continuous_tare_filter_helper = ContinuousTareFilterHelper(config,
            sensor, trigger_dispatch.get_command_queue())
        # Probe Interface
        self._param_helper = probe.ProbeParameterHelper(config)
        self._cmd_helper = probe.ProbeCommandHelper(config, self)
        self._probe_offsets = probe.ProbeOffsetsHelper(config)
        self._mcu_load_cell_probe = McuLoadCellProbe(config, self._load_cell,
            continuous_tare_filter_helper.get_sos_filter(), config_helper,
            trigger_dispatch)
        load_cell_probing_move = LoadCellProbingMove(config,
            self._mcu_load_cell_probe, self._param_helper,
            continuous_tare_filter_helper, config_helper)
        self._tapping_move = TappingMove(config, self._mcu,
            load_cell_probing_move, self._tap_analysis_helper, config_helper)
        tap_session = TapSession(config, self._tapping_move, self._param_helper,
            nozzle_cleaner, config_helper)
        self._probe_session = probe.ProbeSessionHelper(config,
            self._param_helper, tap_session.start_probe_session)
        # printer integration
        LoadCellProbeCommands(config, load_cell_probing_move)
        probe.ProbeVirtualEndstopDeprecation(config)
        self._printer.add_object('probe', self)

    def get_probe_params(self, gcmd=None):
        return self._param_helper.get_probe_params(gcmd)

    def get_offsets(self):
        return self._probe_offsets.get_offsets()

    def start_probe_session(self, gcmd):
        return self._probe_session.start_probe_session(gcmd)

    def get_status(self, eventtime):
        status = self._cmd_helper.get_status(eventtime)
        status.update(self._load_cell.get_status(eventtime))
        status.update(self._tapping_move.get_status(eventtime))
        return status


def load_config(config):
    return LoadCellPrinterProbe(config)
