This package contains 2 different implementations of the multi-band summary correlogram(MBSC)-based pitch detector for noisy speech.

They are:
1. mbsc.m
2. fast_mbsc_fixedWinlen_tracking.m


- mbsc.m is the code that is used to generate the pitch detection results reported in [1]. 


- fast_mbsc_fixedWinlen_tracking.m is the fast version of mbsc.m that uses a fixed window length instead of multiple window lengths used in the original mbsc.m.  User is advised to tune the window length and voicing detection threshold (vThres) to obtain a desirable performance for the input signals.  The pitch detection performance of this fast version is sub-optimal compared to the original mbsc.m, but it run ~20 times faster.  
The fast_mbsc_fixedWinlen_tracking.m also incorporates an option to activate a pitch continuity tracking algorithm after multi-band summary correlogram computation, to select MBSC peak candidates with F0 value continuity within each voiced segment (initially detected with an SNR-adaptive threshold whose range is dependent on the input vThres).  The boundaries of the initial voice segment can expand or contract depending on the F0 continuity across time, and tracking begins from the frame with the highest MBSC peak amplitude in the voicing segment.  This tracking scheme should be beneficial when the noise present is not harmonic/tonal in nature.  This tracking algorithm is also used in the MBSC pitch detector version in SRI's prosodic speaker ID sub-system submitted for phase II evaluation of the DARPA RATS project.  


[1] L. N. Tan, and A. Alwan, "Multi-Band Summary Correlogram-based Pitch Detection for Noisy Speech", Speech Communication, in press.