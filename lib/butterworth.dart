import 'dart:math' as math;

import 'package:complex/complex.dart';
import 'package:iirjdart/src/band_pass_transform.dart';
import 'package:iirjdart/src/band_stop_transform.dart';
import 'package:iirjdart/src/cascade.dart';
import 'package:iirjdart/src/direct_form.dart';
import 'package:iirjdart/src/high_pass_transform.dart';
import 'package:iirjdart/src/layout_base.dart';
import 'package:iirjdart/src/low_pass_transform.dart';
import 'package:iirjdart/src/peak_transform.dart';

class _AnalogLowPass extends LayoutBase {
  int _nPoles;

  _AnalogLowPass(int nPoles)
      : _nPoles = nPoles,
        super(nPoles) {
    _nPoles = nPoles;
    setNormal(0, 1);
  }

  void design() {
    reset();
    double n2 = (2 * _nPoles).toDouble();
    int pairs = _nPoles ~/ 2;
    for (int i = 0; i < pairs; ++i) {
      Complex c =
          Complex.polar(1.0, math.pi / 2.0 + (2 * i + 1) * math.pi / n2);
      addPoleZeroConjugatePairs(c, Complex.infinity);
    }

    if ((_nPoles & 1) == 1) add(Complex(-1), Complex.infinity);
  }
}

/// Create a Butterworth filter.
///
/// Example:
/// ```dart
///   Butterworth butterworth = new Butterworth();
///   butterworth.bandPass(2,250,50,5);
/// ```
class Butterworth extends Cascade {
  /// Butterworth Low-pass filter.
  ///
  /// Params:
  /// * order - The order of the filter
  /// * sampleRate - The sampling rate of the system
  /// * cutoffFrequency - the cutoff frequency
  /// * directFormType - The filter topology. Default direct_form_II
  void lowPass(int order, double sampleRate, double cutoffFrequency,
      [int directFormType = DirectFormAbstract.direct_form_II]) {
    _AnalogLowPass analogProto = _AnalogLowPass(order);
    analogProto.design();
    LayoutBase digitalProto = LayoutBase(order);
    LowPassTransform(cutoffFrequency / sampleRate, digitalProto, analogProto);
    setLayout(digitalProto, directFormType);
  }

  /// Butterworth High-pass filter.
  ///
  /// Params:
  /// * order - Filter order (ideally only even orders)
  /// * sampleRate - Sampling rate of the system
  /// * cutoffFrequency - Cutoff of the system
  /// * directFormType - The filter topology. Default direct_form_II
  void highPass(int order, double sampleRate, double cutoffFrequency,
      [int directFormType = DirectFormAbstract.direct_form_II]) {
    _AnalogLowPass analogProto = _AnalogLowPass(order);
    analogProto.design();
    LayoutBase digitalProto = LayoutBase(order);
    HighPassTransform(cutoffFrequency / sampleRate, digitalProto, analogProto);
    setLayout(digitalProto, directFormType);
  }

  /// Butterworth Bandstop filter.
  ///
  /// Params:
  /// * order - Filter order (actual order is twice)
  /// * sampleRate - Sampling rate of the system
  /// * centerFrequency - Center frequency
  /// * widthFrequency - Width of the notch
  /// * directFormType - The filter topology. Default direct_form_II
  void bandStop(int order, double sampleRate, double centerFrequency,
      double widthFrequency,
      [int directFormType = DirectFormAbstract.direct_form_II]) {
    _AnalogLowPass analogProto = _AnalogLowPass(order);
    analogProto.design();
    LayoutBase digitalProto = LayoutBase(order * 2);
    BandStopTransform(centerFrequency / sampleRate, widthFrequency / sampleRate,
        digitalProto, analogProto);
    setLayout(digitalProto, directFormType);
  }

  /// Butterworth Bandpass filter.
  ///
  /// Params:
  /// * order - Filter order
  /// * sampleRate - Sampling rate
  /// * centerFrequency - Center frequency
  /// * widthFrequency - Width of the notch
  /// * directFormType - The filter topology. Default direct_form_II
  void bandPass(int order, double sampleRate, double centerFrequency,
      double widthFrequency,
      [int directFormType = DirectFormAbstract.direct_form_II]) {
    _AnalogLowPass analogProto = _AnalogLowPass(order);
    analogProto.design();
    LayoutBase digitalProto = LayoutBase(order * 2);
    BandPassTransform(centerFrequency / sampleRate, widthFrequency / sampleRate,
        digitalProto, analogProto);
    setLayout(digitalProto, directFormType);
  }

  /// Butterworth Peak filter.
  ///
  /// Params:
  /// * order - Filter order
  /// * sampleRate - Sampling rate
  /// * centerFrequency - Center frequency
  /// * widthFrequency - Width of the notch
  /// * directFormType - The filter topology. Default direct_form_II
  void peak(int order, double sampleRate, double centerFrequency,
      double widthFrequency,
      [int directFormType = DirectFormAbstract.direct_form_II]) {
    _AnalogLowPass analogProto = _AnalogLowPass(order);
    analogProto.design();
    LayoutBase digitalProto = LayoutBase(order * 2);
    PeakTransform(centerFrequency / sampleRate, widthFrequency / sampleRate,
        digitalProto, analogProto);
    setLayout(digitalProto, directFormType);
  }
  /*
    /*
     * Formulas from https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
     */
user defined parameters:

    Fs (the sampling frequency)

    f0 ("wherever it's happenin', man."  Center Frequency or
        Corner Frequency, or shelf midpoint frequency, depending
        on which filter type.  The "significant frequency".)

    dBgain (used only for peaking and shelving filters)

    Q (the EE kind of definition, except for peakingEQ in which A*Q is
        the classic EE Q.  That adjustment in definition was made so that
        a boost of N dB followed by a cut of N dB for identical Q and
        f0/Fs results in a precisely flat unity gain filter or "wire".)

     _or_ BW, the bandwidth in octaves (between -3 dB frequencies for BPF
        and notch or between midpoint (dBgain/2) gain frequencies for
        peaking EQ)

     _or_ S, a "shelf slope" parameter (for shelving EQ only).  When S = 1,
        the shelf slope is as steep as it can be and remain monotonically
        increasing or decreasing gain with frequency.  The shelf slope, in
        dB/octave, remains proportional to S for all other values for a
        fixed f0/Fs and dBgain.



Then compute a few intermediate variables:

    A  = sqrt( 10^(dBgain/20) )
       =       10^(dBgain/40)     (for peaking and shelving EQ filters only)

    w0 = 2*pi*f0/Fs

    cos(w0)
    sin(w0)

    alpha = sin(w0)/(2*Q)                                       (case: Q)
          = sin(w0)*sinh( ln(2)/2 * BW * w0/sin(w0) )           (case: BW)
          = sin(w0)/2 * sqrt( (A + 1/A)*(1/S - 1) + 2 )         (case: S)

        FYI: The relationship between bandwidth and Q is
             1/Q = 2*sinh(ln(2)/2*BW*w0/sin(w0))     (digital filter w BLT)
        or   1/Q = 2*sinh(ln(2)/2*BW)             (analog filter prototype)

        The relationship between shelf slope and Q is
             1/Q = sqrt((A + 1/A)*(1/S - 1) + 2)

    2*sqrt(A)*alpha  =  sin(w0) * sqrt( (A^2 + 1)*(1/S - 1) + 2*A )
        is a handy intermediate variable for shelving EQ filters.


Finally, compute the coefficients for whichever filter type you want:
   (The analog prototypes, H(s), are shown for each filter
        type for normalized frequency.)

    // H(s) = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)

peakingEQ:  H(s) = (s^2 + s*(A/Q) + 1) / (s^2 + s/(A*Q) + 1)

            b0 =   1 + alpha*A
            b1 =  -2*cos(w0)
            b2 =   1 - alpha*A
            a0 =   1 + alpha/A
            a1 =  -2*cos(w0)
            a2 =   1 - alpha/A


    alpha = sin(w0)/(2*Q)                                       (case: Q)
          = sin(w0)*sinh( ln(2)/2 * BW * w0/sin(w0) )           (case: BW)
          = sin(w0)/2 * sqrt( (A + 1/A)*(1/S - 1) + 2 )         (case: S)


            b0 + b1*z^-1 + b2*z^-2
    H(z) = ------------------------                                  (Eq 1)
            a0 + a1*z^-1 + a2*z^-2


    peak: function (params) {
      var coeffs = initCoeffs()
      var p = preCalcGain(params)
      coeffs.k = 1
      coeffs.a0 = 1 + p.alpha / p.A
      coeffs.a.push(-2 * p.cw / coeffs.a0)
      coeffs.a.push((1 - p.alpha / p.A) / coeffs.a0)
      coeffs.b.push((1 + p.alpha * p.A) / coeffs.a0)
      coeffs.b.push(-2 * p.cw / coeffs.a0)
      coeffs.b.push((1 - p.alpha * p.A) / coeffs.a0)
      return coeffs
    },


    // H(s) = 1 / (s^2 + s/Q + 1)
    lowpass: function (params) {
      var coeffs = initCoeffs()
      if (params.BW) {
        delete params.BW
      }
      var p = preCalc(params, coeffs)
      if (params.preGain) {
        coeffs.k = (1 - p.cw) * 0.5
        coeffs.b.push(1 / (p.a0))
      } else {
        coeffs.k = 1
        coeffs.b.push((1 - p.cw) / (2 * p.a0))
      }
      coeffs.b.push(2 * coeffs.b[0])
      coeffs.b.push(coeffs.b[0])
      return coeffs
    },

    // H(s) = s^2 / (s^2 + s/Q + 1)
    highpass: function (params) {
      var coeffs = initCoeffs()
      if (params.BW) {
        delete params.BW
      }
      var p = preCalc(params, coeffs)
      if (params.preGain) {
        coeffs.k = (1 + p.cw) * 0.5
        coeffs.b.push(1 / (p.a0))
      } else {
        coeffs.k = 1
        coeffs.b.push((1 + p.cw) / (2 * p.a0))
      }
      coeffs.b.push(-2 * coeffs.b[0])
      coeffs.b.push(coeffs.b[0])
      return coeffs
    },

    // H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
    allpass: function (params) {
      var coeffs = initCoeffs()
      if (params.BW) {
        delete params.BW
      }
      var p = preCalc(params, coeffs)
      coeffs.k = 1
      coeffs.b.push((1 - p.alpha) / p.a0)
      coeffs.b.push(-2 * p.cw / p.a0)
      coeffs.b.push((1 + p.alpha) / p.a0)
      return coeffs
    },

    // H(s) = s / (s^2 + s/Q + 1)
    bandpassQ: function (params) {
      var coeffs = initCoeffs()
      var p = preCalc(params, coeffs)
      coeffs.k = 1
      coeffs.b.push(p.alpha * params.Q / p.a0)
      coeffs.b.push(0)
      coeffs.b.push(-coeffs.b[0])
      return coeffs
    },

    // H(s) = (s/Q) / (s^2 + s/Q + 1)
    bandpass: function (params) {
      var coeffs = initCoeffs()
      var p = preCalc(params, coeffs)
      coeffs.k = 1
      coeffs.b.push(p.alpha / p.a0)
      coeffs.b.push(0)
      coeffs.b.push(-coeffs.b[0])
      return coeffs
    },

    // H(s) = (s^2 + 1) / (s^2 + s/Q + 1)
    bandstop: function (params) {
      var coeffs = initCoeffs()
      var p = preCalc(params, coeffs)
      coeffs.k = 1
      coeffs.b.push(1 / p.a0)
      coeffs.b.push(-2 * p.cw / p.a0)
      coeffs.b.push(coeffs.b[0])
      return coeffs
    },


  */
}
