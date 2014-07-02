
"""
This module wraps the process of "linking" to GSL.
Doesn't have to be a module itself; it could just be top matter in
classMarker.py.
"""

import ctypes as c
from ctypes.util import find_library

# Find the full library names in a (platform-independent way).
CBLAS_LIBNAME = find_library('gslcblas')
GSL_LIBNAME   = find_library('gsl')

# Load the BLAS implementation on which GSL depends (but we do not).
# On Windohs the RTLD_GLOBAL flag will be ignored.
c.CDLL(CBLAS_LIBNAME,mode=c.RTLD_GLOBAL)
libgsl = c.CDLL(GSL_LIBNAME)

# Get "pointers" to the desired functions from libgsl
norm_cdf = libgsl.gsl_cdf_ugaussian_P
chi2_cdf = libgsl.gsl_cdf_chisq_P

# Inform Python that both functions take and return doubles
# (ctypes defaults everything to int's unless otherwise informed.)
norm_cdf.argtypes = [ c.c_double, ]
norm_cdf.restype  =   c.c_double
chi2_cdf.argtypes = [ c.c_double, c.c_double ]
chi2_cdf.restype  =   c.c_double

# The C functions are now available to importers of this module as
# gsl.norm_cdf and gsl.chi2_cdf.

"""
Historical record.
The following was used to generate the following data with which to compare
ctypes-mediated calls to GSL with scipy.stats calls.

import sys
from scipy import stats
import random

n = int(sys.argv[1])
while n > 0:
	n -= 1
	chisq = random.uniform(0.0,20.0)
	df    = random.randint(1,4)
	print "stats.chi2.cdf( %f, %d) == %f" % ( chisq, df, stats.chi2.cdf( chisq, df ) )

n = int(sys.argv[1])
while n > 0:
	n -= 1
	rv = random.uniform(-4,4)
	print "stats.norm.cdf( %f) == %f" % ( rv, stats.norm.cdf( rv ) )


stats.chi2.cdf( 3.837189, 4) == 0.571512
stats.chi2.cdf( 10.134084, 2) == 0.993699
stats.chi2.cdf( 8.909714, 1) == 0.997163
stats.chi2.cdf( 9.817288, 3) == 0.979815
stats.chi2.cdf( 11.453361, 1) == 0.999286
stats.chi2.cdf( 1.614816, 2) == 0.553987
stats.chi2.cdf( 2.113167, 4) == 0.285046
stats.chi2.cdf( 15.292713, 3) == 0.998417
stats.chi2.cdf( 19.419882, 4) == 0.999350
stats.chi2.cdf( 4.994303, 4) == 0.712117
stats.chi2.cdf( 11.189038, 3) == 0.989254
stats.chi2.cdf( 5.168579, 2) == 0.924550
stats.chi2.cdf( 3.180021, 3) == 0.635306
stats.chi2.cdf( 5.142184, 2) == 0.923548
stats.chi2.cdf( 4.756880, 3) == 0.809510
stats.chi2.cdf( 3.310248, 2) == 0.808932
stats.chi2.cdf( 2.179981, 2) == 0.663780
stats.chi2.cdf( 14.259638, 4) == 0.993489
stats.chi2.cdf( 10.695632, 1) == 0.998926
stats.chi2.cdf( 13.579863, 1) == 0.999771
stats.chi2.cdf( 16.422434, 3) == 0.999071
stats.chi2.cdf( 13.993411, 4) == 0.992684
stats.chi2.cdf( 16.793328, 4) == 0.997880
stats.chi2.cdf( 7.584739, 3) == 0.944579
stats.chi2.cdf( 15.531982, 2) == 0.999576
stats.chi2.cdf( 14.474858, 4) == 0.994076
stats.chi2.cdf( 13.790116, 4) == 0.992004
stats.chi2.cdf( 10.696579, 4) == 0.969806
stats.chi2.cdf( 13.866394, 2) == 0.999025
stats.chi2.cdf( 9.512379, 2) == 0.991402
stats.chi2.cdf( 9.311826, 3) == 0.974580
stats.chi2.cdf( 12.803786, 3) == 0.994919
stats.chi2.cdf( 11.116299, 3) == 0.988887
stats.chi2.cdf( 2.101988, 1) == 0.852892
stats.chi2.cdf( 2.146142, 3) == 0.457366
stats.chi2.cdf( 6.150674, 1) == 0.986864
stats.chi2.cdf( 6.923475, 2) == 0.968625
stats.chi2.cdf( 15.870639, 1) == 0.999932
stats.chi2.cdf( 16.031727, 3) == 0.998883
stats.chi2.cdf( 15.345634, 2) == 0.999535
stats.chi2.cdf( 17.830873, 4) == 0.998668
stats.chi2.cdf( 0.381699, 3) == 0.056003
stats.chi2.cdf( 12.154899, 1) == 0.999510
stats.chi2.cdf( 11.468149, 2) == 0.996766
stats.chi2.cdf( 11.655731, 2) == 0.997056
stats.chi2.cdf( 11.971961, 1) == 0.999460
stats.chi2.cdf( 6.068403, 1) == 0.986238
stats.chi2.cdf( 4.750680, 4) == 0.686151
stats.chi2.cdf( 18.567197, 2) == 0.999907
stats.chi2.cdf( 3.499957, 4) == 0.522115
stats.chi2.cdf( 18.390078, 4) == 0.998965
stats.chi2.cdf( 14.275222, 3) == 0.997447
stats.chi2.cdf( 5.602101, 3) == 0.867342
stats.chi2.cdf( 1.984884, 1) == 0.841123
stats.chi2.cdf( 2.424721, 4) == 0.341836
stats.chi2.cdf( 12.534791, 3) == 0.994241
stats.chi2.cdf( 15.683913, 2) == 0.999607
stats.chi2.cdf( 4.933558, 2) == 0.915142
stats.chi2.cdf( 2.189354, 4) == 0.299021
stats.chi2.cdf( 1.701873, 2) == 0.572985
stats.chi2.cdf( 3.810791, 4) == 0.567782
stats.chi2.cdf( 0.756727, 3) == 0.140213
stats.chi2.cdf( 1.717149, 4) == 0.212400
stats.chi2.cdf( 7.657541, 4) == 0.895040
stats.chi2.cdf( 15.969858, 2) == 0.999659
stats.chi2.cdf( 16.185774, 4) == 0.997220
stats.chi2.cdf( 19.553636, 1) == 0.999990
stats.chi2.cdf( 0.595575, 1) == 0.559729
stats.chi2.cdf( 7.330577, 4) == 0.880584
stats.chi2.cdf( 16.752028, 1) == 0.999957
stats.chi2.cdf( 13.859777, 2) == 0.999022
stats.chi2.cdf( 19.105605, 1) == 0.999988
stats.chi2.cdf( 1.348244, 3) == 0.282290
stats.chi2.cdf( 9.854555, 4) == 0.957050
stats.chi2.cdf( 1.045307, 2) == 0.407055
stats.chi2.cdf( 11.977312, 3) == 0.992539
stats.chi2.cdf( 4.067660, 3) == 0.745750
stats.chi2.cdf( 4.847933, 1) == 0.972321
stats.chi2.cdf( 15.662935, 1) == 0.999924
stats.chi2.cdf( 0.144975, 1) == 0.296616
stats.chi2.cdf( 15.061462, 2) == 0.999464
stats.chi2.cdf( 3.664815, 3) == 0.699993
stats.chi2.cdf( 17.048486, 1) == 0.999964
stats.chi2.cdf( 18.053538, 4) == 0.998795
stats.chi2.cdf( 4.188982, 4) == 0.618966
stats.chi2.cdf( 19.136888, 3) == 0.999744
stats.chi2.cdf( 16.928155, 4) == 0.998004
stats.chi2.cdf( 6.017325, 1) == 0.985834
stats.chi2.cdf( 11.961238, 2) == 0.997473
stats.chi2.cdf( 17.165805, 1) == 0.999966
stats.chi2.cdf( 1.563068, 3) == 0.332207
stats.chi2.cdf( 13.150381, 3) == 0.995678
stats.chi2.cdf( 10.889144, 3) == 0.987659
stats.chi2.cdf( 19.312613, 4) == 0.999318
stats.chi2.cdf( 3.794773, 1) == 0.948587
stats.chi2.cdf( 10.839998, 3) == 0.987377
stats.chi2.cdf( 1.678256, 2) == 0.567913
stats.chi2.cdf( 6.463182, 2) == 0.960505
stats.chi2.cdf( 9.388858, 2) == 0.990854
stats.chi2.cdf( 15.479695, 2) == 0.999565
stats.norm.cdf( -3.222622) == 0.000635
stats.norm.cdf( 1.574812) == 0.942350
stats.norm.cdf( -0.596397) == 0.275455
stats.norm.cdf( -3.856535) == 0.000058
stats.norm.cdf( -1.077643) == 0.140596
stats.norm.cdf( 2.917213) == 0.998234
stats.norm.cdf( 1.243000) == 0.893066
stats.norm.cdf( -1.876978) == 0.030261
stats.norm.cdf( -2.715503) == 0.003309
stats.norm.cdf( 2.337331) == 0.990289
stats.norm.cdf( -3.236118) == 0.000606
stats.norm.cdf( -2.735808) == 0.003111
stats.norm.cdf( 2.749945) == 0.997020
stats.norm.cdf( 1.144781) == 0.873850
stats.norm.cdf( 2.397000) == 0.991735
stats.norm.cdf( 2.625147) == 0.995669
stats.norm.cdf( -0.559216) == 0.288007
stats.norm.cdf( 2.946350) == 0.998392
stats.norm.cdf( 0.722608) == 0.765039
stats.norm.cdf( 3.460202) == 0.999730
stats.norm.cdf( -2.594993) == 0.004730
stats.norm.cdf( -1.942357) == 0.026047
stats.norm.cdf( 3.623008) == 0.999854
stats.norm.cdf( -0.547400) == 0.292052
stats.norm.cdf( 1.154172) == 0.875785
stats.norm.cdf( 3.187919) == 0.999283
stats.norm.cdf( 0.713406) == 0.762203
stats.norm.cdf( 1.416465) == 0.921680
stats.norm.cdf( -3.113847) == 0.000923
stats.norm.cdf( -3.764906) == 0.000083
stats.norm.cdf( -1.302427) == 0.096385
stats.norm.cdf( -2.765894) == 0.002838
stats.norm.cdf( -3.109307) == 0.000938
stats.norm.cdf( 1.648201) == 0.950344
stats.norm.cdf( 0.489791) == 0.687859
stats.norm.cdf( -3.914136) == 0.000045
stats.norm.cdf( -2.559358) == 0.005243
stats.norm.cdf( 1.537968) == 0.937972
stats.norm.cdf( -1.566247) == 0.058645
stats.norm.cdf( -1.629578) == 0.051595
stats.norm.cdf( 0.462983) == 0.678312
stats.norm.cdf( -3.375428) == 0.000369
stats.norm.cdf( -1.940325) == 0.026170
stats.norm.cdf( -0.005296) == 0.497887
stats.norm.cdf( -3.486196) == 0.000245
stats.norm.cdf( -1.450274) == 0.073491
stats.norm.cdf( 0.003383) == 0.501350
stats.norm.cdf( -2.029284) == 0.021215
stats.norm.cdf( -1.687629) == 0.045741
stats.norm.cdf( 3.730042) == 0.999904
stats.norm.cdf( 2.397392) == 0.991744
stats.norm.cdf( 1.870999) == 0.969327
stats.norm.cdf( 3.918413) == 0.999955
stats.norm.cdf( 0.715808) == 0.762945
stats.norm.cdf( 3.547510) == 0.999806
stats.norm.cdf( -1.878236) == 0.030174
stats.norm.cdf( 1.752621) == 0.960166
stats.norm.cdf( -2.930430) == 0.001692
stats.norm.cdf( 1.451507) == 0.926681
stats.norm.cdf( -2.865447) == 0.002082
stats.norm.cdf( -3.114890) == 0.000920
stats.norm.cdf( 2.760925) == 0.997118
stats.norm.cdf( -1.569998) == 0.058208
stats.norm.cdf( 0.477666) == 0.683556
stats.norm.cdf( 1.669559) == 0.952497
stats.norm.cdf( -0.980465) == 0.163428
stats.norm.cdf( 2.795951) == 0.997413
stats.norm.cdf( 0.455132) == 0.675493
stats.norm.cdf( -2.079290) == 0.018795
stats.norm.cdf( -2.410700) == 0.007961
stats.norm.cdf( 2.318680) == 0.989794
stats.norm.cdf( -1.788332) == 0.036861
stats.norm.cdf( -3.149278) == 0.000818
stats.norm.cdf( -1.894779) == 0.029061
stats.norm.cdf( 0.320028) == 0.625526
stats.norm.cdf( -1.352479) == 0.088111
stats.norm.cdf( -0.749966) == 0.226638
stats.norm.cdf( -1.480541) == 0.069365
stats.norm.cdf( -1.135005) == 0.128187
stats.norm.cdf( 0.487274) == 0.686968
stats.norm.cdf( -3.530569) == 0.000207
stats.norm.cdf( -1.497954) == 0.067073
stats.norm.cdf( -1.741887) == 0.040764
stats.norm.cdf( -0.535194) == 0.296258
stats.norm.cdf( -2.587910) == 0.004828
stats.norm.cdf( -0.605991) == 0.272260
stats.norm.cdf( 3.573184) == 0.999824
stats.norm.cdf( -3.877482) == 0.000053
stats.norm.cdf( 0.817312) == 0.793125
stats.norm.cdf( 3.896229) == 0.999951
stats.norm.cdf( 2.839063) == 0.997738
stats.norm.cdf( -1.148386) == 0.125405
stats.norm.cdf( -0.343865) == 0.365474
stats.norm.cdf( 0.201503) == 0.579847
stats.norm.cdf( 2.175061) == 0.985187
stats.norm.cdf( 3.927386) == 0.999957
stats.norm.cdf( -0.398771) == 0.345031
stats.norm.cdf( -0.490063) == 0.312045
stats.norm.cdf( -2.305063) == 0.010582
tats.norm.cdf( 1.451831) == 0.926726
"""

