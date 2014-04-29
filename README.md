qHTSPipelineGUI
===============

A R Shiny Interface for running Curvep and the relevant R codes

Input
-----

Tab-delimited files with four sets of required columns:  

- pathway unique column: *uniqueID* column (i.e., an ID unique to the row based on pathway)
- concentration columns: *conc[0-9]+* columns (e.g., conc0, conc1, conc2, ... etc.)
  - unit: nM, uM, M, or log10(M) 
- response columns: *resp[0-9]+* columns (e.g., resp0, resp1, resp2, ... etc.)  
  - unit: -100% - 0% - 100%; 0% is the baseline
- pathway column: *pathway* column


Additional Input
----------------

- mask column: *Mask.Flags* or *curvep_mask* (Mask.Flags column has higher priority)
- cytotoxicity data: the required columns plus the *curvep_r[0-9]+* columns 
- unused columns will be kept


Parameters
----------

The capital letters are the parameters used in Curvep command line. 

- primary
  - max.deviation (MXDV, default=5%): the size of spike allowed
  - baseline threshold (THR, default=15%): the amount of baseline noise allowed
  - max.response (RNG, default=-100%): the highest response allowed (also an indication of monotonicity direction)
  
- secondary
  - carryover threshold (CRO, default=80%): the amount of low-concentration response required to be potentially considered as true signals
  - min.#points for u-shape (USHAPE, default=4): the minimum number of points to form a U-Shape (i.e. non-monotonic response); otherwise noise
  - min.#flat points to detect baseline-shift (BSHIFT, default=3): the minimum number of points needed on the flat portion of the curve to call it a baseline shift
  - baseline shift correction mode (BLFX, default=FALSE): if TRUE, the whole curve will vertically move so the lower asymptote reaches the baseline instead only the lower asymptote (FALSE) 
  - favors corrections based on low conc-s (BYHI, default=TRUE): when corrections are needed, the corrections with more points at the lower test concentrations are retained
  - allow extrapolation beyond test conc. boundaries (XPLAX, default=FALSE): particular for the calculation of POD/AC50 in the case of the curve where only the upper asymptote is available (e.g., when the potency exceeds the tested concentration range)  
  - use cytotoxicity data to generate mask (default=FALSE): cytotoxicity curvep data (when available) are used to create mask for the non-monotonic curves (caused by cytotoxicity by assumption) in the activation-type assays.  

Output
----------------

- *curvep_wauc* (weight area under the curve)
- *curvep_pod* (point-of-depature)
- *curvep_remark* (curve indication)
- *curvep_n_corrections* (total number of correntions have been made)
- *curvep_mask* (outlier points generated)
- *curvep_r[0-9]+* (the responses after curvep treatment)
- more (AC50, C% ... under construction but provided in Curvep command line)
