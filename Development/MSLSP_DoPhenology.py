import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import datetime 
from rpy2.robjects import r
import re
from itertools import groupby
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from csaps import csaps
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# # Loading in Rdata into python
# robjects.r['load']('MSLSP_test_input.RData')

# b2 = np.array(robjects.r['b2'])
# b3 = np.array(robjects.r['b3'])
# b4 = np.array(robjects.r['b4'])
# b5 = np.array(robjects.r['b5'])
# b6 = np.array(robjects.r['b6'])
# b7 = np.array(robjects.r['b7'])
# vi = np.array(robjects.r['vi'])
# snowPix = np.array(robjects.r['snowPix'])
# dates = np.array(robjects.r['dates'])
# imgYrs = np.array(robjects.r['imgYrs'])
# phenYrs = np.array(robjects.r['phenYrs'])
# splineStart = np.array(robjects.r['splineStart'])
# splineEnd = np.array(robjects.r['splineEnd'])
# numDaysFit = np.array(robjects.r['numDaysFit'])
# pheno_pars = robjects.r('pheno_pars')
# pheno_pars = dict(zip(pheno_pars.names, list(pheno_pars)))

def CheckSpike_MultiBand(blue, red, vi, dates, pheno_pars):
    ''' Despike time-series
        Uses the Maja approach to get positive time-series spikes associated with clouds
        See Section 2.1.4 in this document: https://www.theia-land.fr/sites/default/files/imce/produits/atbd_maja_071217.pdf
        The maja approach relies on both the blue and red bands (If change is red is 1.5x the change in blue, it is likely a land condition change)
        In addition, a 3 point method is used to look for negative spikes in EVI2. This step is added because the Maja approach won't capture
        spikes associated with cloud shadow
        Vectorized to avoid looping through dates
        Douglas Bolton
    '''
    blue_og = blue # Preserve original vector
    red_og = red # Preserve original vector
    vi_og = vi # Preserve original vector
    dates_og = dates

    good = ~np.isnan(blue_og) & ~np.isnan(red_og) & ~np.isnan(vi_og) 

    x_outs = np.zeros(len(blue_og), dtype=bool) # Create the outlier output vector
    count = 0
    while (count < pheno_pars['maxDespikeIterations'][0]):
        count += 1

        bS = blue_og[good] # Subset to non missing values
        rS = red_og[good]
        eS = vi_og[good]
        dS = dates_og[good] # Subset date vector

        ind1 = np.arange(0, len(dS)-2) # Get indices for first, second, and third images
        ind2 = np.arange(1, len(dS)-1)
        ind3 = np.arange(2, len(dS))

        dDiff1 = dS[ind2] - dS[ind1]
        bDiff1 = bS[ind2] - bS[ind1]
        rDiff1 = rS[ind2] - rS[ind1]
        bTest1 = bDiff1 > (0.03 * (1 + dDiff1/30))
        rTest1 = rDiff1 < (1.5 * bDiff1)

        dDiff2 = dS[ind3] - dS[ind2]
        bDiff2 = bS[ind2] - bS[ind3] # 2 minus 3 because we are investigating 2 as a peak
        rDiff2 = rS[ind2] - rS[ind3]
        bTest2 = bDiff2 > (0.03 * (1 + dDiff2/30))
        rTest2 = rDiff2 < (1.5 * bDiff2)

        majaTest = bTest1 & rTest1 & bTest2 & rTest2

        dayFrac = (dS[ind2] - dS[ind1]) / (dS[ind3] - dS[ind1]) # Calculate time fraction of date 1 to 2 compared to date 1 to 3
        fitVal = eS[ind1] + (eS[ind3] - eS[ind1]) * dayFrac     # Calculate value at point 2 if a straight line is drawn from point 1 to 3

        dev = fitVal - eS[ind2]
        devRatio = dev / (eS[ind3] - eS[ind1])
        dDiff = dS[ind3] - dS[ind1]

        # Look for negative spikes in vi
        eTest = (dev > pheno_pars['minResid'][0]) & (abs(devRatio) > pheno_pars['spikeThresh'][0]) & (dDiff < pheno_pars['maxDistance'][0])

        check = majaTest | eTest
        check = np.r_[False, check, False]
        check[np.isnan(check) | np.isinf(check)] = False
        if (np.sum(check) == 0):
            break # Break if no observations have been despiked

        x_outs[good] = check # Expand to size of original x, accounting for missing values
        good[x_outs] = False # Remove the despiked values from the pixels of interest and try again
    
    return x_outs

def Screen_SnowFills(x, dormVal, snowPix, dates, pheno_pars):
    ''' Check the snow flags to see if the dormant value would be a good fit
        Specifically, if filling with the dormant value creates a negative spike in the time series, we will remove the snow flag
        Outputs an updated vector of snow flags
        Douglas Bolton
    '''
    # Reduce vectors to valid observations and snow pixels
    sORg = ~np.isnan(x) | snowPix == True
    snowSub = snowPix[sORg]
    xSub = x[sORg]
    d = dates[sORg]
    padded = np.r_[0, np.convolve(snowSub, np.ones(3,dtype=int),'valid'), 0]         # Get rolling sum of snow detections
    padded2 = np.r_[3, pd.Series(padded).rolling(3).max().dropna().astype('int'), 3] # Get rolling max of the sums. Setting first and last to 3 keeps these values if they are snow
    
    # If the rolling max is less than 3, means center value isn't part of three consecutive snow values
    # So we will only check values less than 3 to determine if they are spikes
    toCheck = np.where((padded2 < 3) & (snowSub == 1))[0] 
    noNA = np.where(~np.isnan(xSub))[0]
    rmIt = np.full(len(snowSub), False)
    for ind2 in toCheck:
        inds1 = noNA[noNA < ind2]
        inds3 = noNA[noNA > ind2]
        if len(inds1) == 0 or len(inds3) == 0:  # If there are no observations
            rmIt[ind2] = True
            continue
        ind1 = max(inds1)
        ind3 = min(inds3)

        dayFrac = (d[ind2]-d[ind1]) / (d[ind3]-d[ind1])           # Calculate time fraction of date 1 to 2 compared to date 1 to 3
        fitVal = xSub[ind1] + (xSub[ind3] - xSub[ind1]) * dayFrac # Calculate value at point 2 if a straight line is drawn from point 1 to 3
        dev = fitVal - dormVal
        devRatio = dev / (xSub[ind3] - xSub[ind1])
    
        rmIt[ind2] = ((dev > pheno_pars['minResid'][0]) & ((abs(devRatio)) > pheno_pars['spikeThresh']))

    rmFull = np.full(len(snowPix), False)
    rmFull[sORg] = rmIt
    snowPix[rmFull] = False
    return snowPix

# def Smooth_VI(x, dates, pred_dates, weights, pheno_pars, dormant_value):
#     ''' Fit a cubic spline to the HLS time-series
#         Written by Josh Gray and Douglas Bolton
#     '''
#     # Get index of pixels with good values
#     ind = ~np.isnan(x)
#     # Smooth with a spline to get continuous daily series
#     r_smooth_spline = robjects.r['smooth.spline'] #extract R function
#     spl = r_smooth_spline(dates[ind], x[ind], spar = pheno_pars['splineSpar'][0], w=weights[ind])
#     # Weighted version
#     xSmooth = robjects.r['predict'](spl, pred_dates).rx2('y')
#     # Screen and fill values less than the dormant value
#     xSmooth[xSmooth < dormant_value] = dormant_value 
#     return(xSmooth)

# def Smooth_VI2(x, dates, pred_dates, weights, pheno_pars, dormant_value):
#     ''' CSAPS Version (Not using)
#     '''
#     #Get index of pixels with good values
#     ind = ~np.isnan(x)
    
#     # Sort dates, x, and weights to the same indices
#     datesInd = np.argsort(dates[ind])
#     x = x[datesInd]
#     weights = weights[datesInd]

#     # Get unique date values and corresponding indices
#     uniqueDates, indicesList = np.unique(dates[datesInd], return_index=True)
#     # Smooth with a spline to get continuous daily series
#     spl = csaps(uniqueDates, x[indicesList], weights=weights[indicesList], smooth=0.25)
#     xSmooth = spl(pred_dates)
    
#     # Screen and fill values less than the dormant value
#     xSmooth[xSmooth < dormant_value] = dormant_value 
#     return xSmooth

def Smooth_VI3(x, dates, pred_dates, weights, dormant_value):
    ''' Testing smooth spline with scipy (first call)
    '''
    #Get index of pixels with good values
    ind = ~np.isnan(x)
    
    # Sort dates, x, and weights to the same indices
    datesInd = np.argsort(dates[ind])
    x = x[datesInd]
    weights = weights[datesInd]

    # Get unique date values and corresponding indices
    uniqueDates, indicesList = np.unique(dates[datesInd], return_index=True)

    # Smooth with a spline to get continuous daily series
    spl = UnivariateSpline(uniqueDates, x[indicesList], w=weights[indicesList], s=0.0235)
    # Weighted version
    xSmooth = spl(pred_dates)

    # Screen and fill values less than the dormant value
    xSmooth[xSmooth < dormant_value] = dormant_value

    return xSmooth

def Smooth_VI4(x, dates, pred_dates, weights, dormant_value):
    ''' Testing smooth spline with scipy (second call)
    '''
    #Get index of pixels with good values
    ind = ~np.isnan(x)
    
    # Sort dates, x, and weights to the same indices
    datesInd = np.argsort(dates[ind])
    x = x[datesInd]
    weights = weights[datesInd]

    # Get unique date values and corresponding indices
    uniqueDates, indicesList = np.unique(dates[datesInd], return_index=True)

    # Smooth with a spline to get continuous daily series
    spl = UnivariateSpline(uniqueDates, x[indicesList], w=weights[indicesList], s=0.03)
    # Weighted version
    xSmooth = spl(pred_dates)

    # Screen and fill values less than the dormant value
    xSmooth[xSmooth < dormant_value] = dormant_value

    return xSmooth


# def Smooth_Bands(x, dates, pred_dates, weights, pheno_pars):
#     # Get index of pixels with good values
#     ind = ~np.isnan(x)
#     # Smooth with a spline to get continuous daily series
#     r_smooth_spline = robjects.r['smooth.spline'] #extract R function
#     print(np.sort(dates[ind]))
#     # print(x[ind])
#     spl = r_smooth_spline(dates[ind], x[ind], spar = pheno_pars['splineSpar'][0], w=weights[ind])
#     # Weighted version
#     xSmooth = robjects.r['predict'](spl, pred_dates).rx2('y')

#     # plt.plot(dates, x, 'o', pred_dates, xSmooth, '-')
#     # plt.title("Original")
#     # plt.show()

#     # Determine upper and lower bound of good data (set in json file)
#     lowBound = np.nanquantile(x, pheno_pars['bandLimits'][0])
#     upperBound = np.nanquantile(x, pheno_pars['bandLimits'][1])

#     xSmooth[xSmooth < lowBound] = lowBound
#     xSmooth[xSmooth > upperBound] = upperBound
#     return xSmooth



def FindPeaks(x, mag_order = True):
    '''New from Josh - 2018-10-31
       Finds time-series peaks
       Josh Gray
    '''
    # Function to identify peaks in time series x (or troughs if x=-x), supports "flat top" peaks
    # if mag_order is TRUE, peaks are returned in order of increasing magnitude (of x)
    d = np.diff(x)
    d_code = (d > 0) + (2 * (d < 0)) # 0=no change, 1=inc, 2=dec
    peaks = np.array([i.start() for i in re.finditer("12", ''.join(d_code.astype(str)))]) # No match is empty list
    if len(peaks) == 0:
        peaks = np.array([np.nan])
    flat_peaks = np.array([i.start() for i in re.finditer("10+2", ''.join(d_code.astype(str)))]) # No match is empty list
    if len(flat_peaks) == 0:
        flat_peaks = np.array([np.nan])
    d_code_rle = np.array([(k, sum(1 for i in g)) for k,g in groupby(d_code)])
    d_code_rlel = np.array([d_code_rle[i][1] for i in range(len(d_code_rle))])
    flat_peaks = flat_peaks + np.round(d_code_rlel[np.where(np.cumsum(d_code_rlel)==flat_peaks)[0] + 1] / 2)
    pre_peaks = np.array([peaks + 1, flat_peaks + 1], dtype=object) # Extra step to create an array combining peaks and flat_peaks
    peaks = np.sort((pre_peaks[np.array([i.size>0 for i in pre_peaks])]))[0]
    if mag_order:
        return peaks[np.argsort(x[peaks])]
    return peaks


def GetSegs(peaks, x, pars, peak=np.nan):
    '''New from Josh - 2018-10-31
       Determines valid segments
       Josh Gray 
    '''
    # Identifies valid increasing segments in x subject to the parameters in pars
    # Returns a list of segments: c(start, peak, end). Don't call directly w/ peak!=NA
    # NOTE: returned segments will not necessarily be in order, and may not completely partition x

    # Ensure that peaks are in increasing order of x's magnitude
    tmp_peaks = peaks[np.argsort(x[peaks])] # So we only have to sort once if they're in the wrong order
    if ~(tmp_peaks==peaks).all():
        peaks = tmp_peaks

    # If no peak is specified, we start at the beginning
    if np.isnan(peak):
        peak = peaks[0]

    # Get the next largest peak; will be NA if this peak is the highest (last one to do)
    if len(peaks) > 1:
        if np.where(peaks == peak)[0][0] + 1 >= len(peaks):
            next_highest_peak = np.nan
        else:
            next_highest_peak = peaks[np.where(peaks == peak)[0][0] + 1]
    else:
        next_highest_peak = np.nan

    # We could have any combination of rel_amp_frac, rel_peak_frac, and min_seg_amplitude specified
    # Initialize seg_thresh and peak_thresh to zero
    # Determine the "global max/min", if peak_frac is specified, set it, if amp_frac is specified, set it
    # If min_seg_amplitude is set, choose the max of that and amp_frac
    seg_thresh = peak_thresh = 0
    global_max = np.nanmax(x[(int(pars['splineBuffer'][0])+1):(int(pars['splineBuffer'][0])+365)]) #find global min/max within a target year
    global_min = np.nanmin(x[(int(pars['splineBuffer'][0])+1):(int(pars['splineBuffer'][0])+365)])
    if ~np.isnan(pars['rel_amp_frac'][0]):
        seg_thresh = (global_max - global_min) * pars['rel_amp_frac'][0]

    if ~np.isnan(pars['min_seg_amplitude'][0]):
        seg_thresh = max(pars['min_seg_amplitude'][0], seg_thresh)

    # Checks if the period preceding the peak covers enough amplitude
    # Search before the peak up to the maximum of: previous peak, the head of x, or the peak - max_increase_length
    previous_peaks = peaks[peaks - peak < 0]
    previous_peak = np.nan
    if len(previous_peaks) > 0:
        previous_peak = max(previous_peaks)
    search_start = int(np.nanmax([0, peak - pars['max_increase_length'][0], previous_peak]))
    search_end = peak
    # Get the index of the closest minimum value within the search window
    # NOTE: should maybe retrieve the troughs here with FindPeaks(-x) instead
    # In the event of repeated minimum values, we take the closest one here
    inc_min_ind = np.nanmax(np.where(x[search_start:search_end] == np.nanmin(x[search_start:search_end]))[0] + search_start) #? Took out -1 for python
    seg_amp = x[peak] - x[inc_min_ind] # get the increasing segment amplitude

    if (seg_amp >= seg_thresh) & (x[peak] >= peak_thresh):
        # Check for a valid decreasing segment
        next_peaks = peaks[peaks - peak > 0]
        next_peak = np.nan
        if len(next_peaks) > 0:
            next_peak = min(next_peaks)
        # Search after the peak up to the minimum of: next peak, the tail of x, or the max_decrease_length
        search_start = peak
        search_end = int(np.nanmin([len(x), peak + pars['max_decrease_length'][0] + 1, next_peak]))
        # Get the index of the closest minimum value within the search window
        # NOTE: see above note about finding troughs instead
        dec_min_ind = np.nanmin(np.where(x[search_start:search_end] == np.nanmin(x[search_start:search_end]))[0] + search_start) #? Added +1 for python Come here if more error arise
        seg_amp = x[peak] - x[dec_min_ind] # get the decreasing segment amplitude
        if seg_amp >= seg_thresh:
            # We found a valid segment, store it as a list with a single vector: c(start, peak , end)
            tmp_seg = np.array([np.r_[inc_min_ind, peak, dec_min_ind]])
            # If this isn't the last peak, then call CheckSegRec again w/ next highest peak
            if ~np.isnan(next_highest_peak):
                return np.vstack((tmp_seg[0], GetSegs(peaks, x, pars, peak=next_highest_peak)))
            else:
                return tmp_seg
        else:
            # Increase was valid, but decrease was not
            if np.isnan(next_highest_peak):
                return np.array([-9999, -9999, -9999])
            elif np.where(peaks == peak)[0][0] == 0:
                peaks = peaks[1:]
            else:
                peaks = np.array([peaks[-np.where(peaks == peak)[0][0]]]) # Remove peak from peaks list
            # If this isn't the last peak, then call CheckSegRec again w/ next highest peak
            if ~np.isnan(next_highest_peak):
                return GetSegs(peaks, x, pars, peak=next_highest_peak)


    else:
        # Increase segment not valid
        if np.isnan(next_highest_peak):
            return np.array([-9999, -9999, -9999])
        elif np.where(peaks == peak)[0][0] == 0:
            peaks = peaks[1:]
        else:
            peaks = np.array([peaks[-np.where(peaks == peak)[0][0]]]) # Remove peak from peaks list
        if ~np.isnan(next_highest_peak): 
            return GetSegs(peaks, x, pars, peak=next_highest_peak)


def GetPhenoDates(segs, x, dates, pheno_pars):
    '''Get phenology dates from segments. Also pull the peak date
       Josh Gray. Updated by Douglas Bolton to include peak date and cleaned 
    '''
    pheno_dates = []

    #Pull greenup dates
    for gup_thresh in pheno_pars['gup_threshes']:
        pheno_dates = np.r_[pheno_dates, dates[[GetSegThresh(segs[i], x, gup_thresh, gup=True) for i in range(len(segs))]]]

    #Pull peak dates
    pheno_dates = np.r_[pheno_dates, dates[[segs[i][1] for i in range(len(segs))]]]

    #Pull greendown dates
    for gdown_thresh in pheno_pars['gdown_threshes']:
        pheno_dates = np.r_[pheno_dates, dates[[GetSegThresh(segs[i], x, gdown_thresh, gup=False) for i in range(len(segs))]]]

    return pheno_dates


def GetThresh(thresh_value, x, first_greater=True, gup=True):
    '''
    Josh Gray
    '''
    # Returns the index of the first/last value of x that is greater/less than the value of thresh
    # If gup is False (greendown) then it returns the first/last value of x that is less/greater than
    # the value of thresh. first/last and greater/less determined by first_greater
    # NOTE: if thresh is 1 or 0, rounding error can be a problem. Now we round the threshold and each
    # of the evi values to 6 decimal places to compensate

    if gup:
        if first_greater:
            return (np.min(np.where(np.round(x, decimals=6) >= np.round(thresh_value, decimals=6))))
        else:  
            return (np.max(np.where(np.round(x, decimals=6) <= np.round(thresh_value, decimals=6))))
    else:
        if first_greater:
            return (np.min(np.where(np.round(x, decimals=6) <= np.round(thresh_value, decimals=6))))
        else:  
            return (np.max(np.where(np.round(x, decimals=6) >= np.round(thresh_value, decimals=6))))


def GetSegThresh(seg, x, thresh, gup=True):
    '''
    Josh Gray
    '''
    if gup:
        # Check for valid greenup segment
        if ~np.isnan(seg[0]) & ~np.isnan(seg[1]):
            gup_thresh = x[seg[0]] + ((x[seg[1]] - x[seg[0]]) * thresh)
            gup_thresh_index = GetThresh(gup_thresh, x[seg[0]:seg[1]], first_greater=True, gup=True)
            return (gup_thresh_index + seg[0])
        else:
            return np.nan

    else:
        # Check for valid greenup segment
        if ~np.isnan(seg[1]) & ~np.isnan(seg[2]):
            gdown_thresh = x[seg[2]] + ((x[seg[1]] - x[seg[2]]) * thresh)
            gdown_thresh_index = GetThresh(gdown_thresh, x[seg[1]:seg[2]], first_greater=False, gup=False)
            return (gdown_thresh_index + seg[1] + 1)
        else:
            return np.nan


def GetSegMetrics(seg, x_smooth, x_raw, smooth_dates, raw_dates):
    if np.isnan(seg).any():
        return np.nan
    # Get the subset of the smoothed and original time series
    tmp_seg_smooth = x_smooth[seg[0]:seg[2]+1]
    tmp_gup_smooth = x_smooth[seg[0]:seg[1]+1]
    tmp_gdown_smooth = x_smooth[seg[1]:seg[2]+1]

    # Get the full segment minimum/maximum SVI
    seg_min = np.nanmin(tmp_seg_smooth)
    seg_max = np.nanmax(tmp_seg_smooth)
    seg_amp = seg_max - seg_min

    # Get the segment integrated SVI: the sum of values.
    # For MODIS C6, this is the sum of values above the minimum evi.
    seg_int = np.sum(tmp_seg_smooth)


    # Organize greenup segment
    ######################################
    gup_raw_date_inds = np.where((raw_dates >= smooth_dates[seg[0]]) & (raw_dates <= smooth_dates[seg[1]])) # indices in raw data of gup segment
    gup_smooth_date_inds = np.where(np.isin(smooth_dates, raw_dates[gup_raw_date_inds])) # indices of raw dates in smooth dates

    raw_dates_gup = raw_dates[gup_raw_date_inds]
    gup_raw_data = x_raw[gup_raw_date_inds] # Get the raw data associated with the gup segment (this is the pre-filled, despiked version)
    gup_smooth_data = x_smooth[gup_smooth_date_inds] # Get the smoothed values associated with each raw data value

    gup_numObs = np.sum(~np.isnan(gup_raw_data))


    # Organize greendown segment
    ######################################
    gdown_raw_date_inds = np.where((raw_dates >= smooth_dates[seg[1]]) & (raw_dates <= smooth_dates[seg[2]])) # Indices in raw data of gdown segment
    gdown_smooth_date_inds = np.where(np.isin(smooth_dates, raw_dates[gdown_raw_date_inds])) # Indices of raw dates in smooth dates

    raw_dates_gdown = raw_dates[gdown_raw_date_inds]
    gdown_raw_data = x_raw[gdown_raw_date_inds] # Get the raw data associated with the gdown segment (this is the pre-filled, despiked version)
    gdown_smooth_data = x_smooth[gdown_smooth_date_inds] # Get the smoothed values associated with each raw data value

    gdown_numObs = np.sum(~np.isnan(gdown_raw_data))


    if (gup_numObs == 0 | gdown_numObs == 0):
        rep = np.empty(9); rep[:] = np.nan
        return rep


    ###Get the observation density for each period
    #This approach counts snow filled values as good values, since snow images are valuable for pinning down dormant period
    ###Find the biggest gap between images
    gup_seg_rsquared = np.array(1 - (np.nansum((gup_raw_data - gup_smooth_data)*(gup_raw_data - gup_smooth_data)) / np.nansum((gup_raw_data - np.nanmean(gup_raw_data))*(gup_raw_data - np.nanmean(gup_raw_data)))))
    gup_seg_rsquared[np.isinf(gup_seg_rsquared)] = np.nan
    gup_maxgap = max(np.diff(np.r_[smooth_dates[seg[0]],raw_dates_gup[~np.isnan(gup_raw_data)],smooth_dates[seg[1]]]))
    gup_maxgap_frac = gup_maxgap / (seg[1] - seg[0])

    gdown_seg_rsquared = np.array(1 - (np.nansum((gdown_raw_data - gdown_smooth_data)*(gdown_raw_data - gdown_smooth_data)) / np.nansum((gdown_raw_data - np.nanmean(gdown_raw_data))*(gdown_raw_data - np.nanmean(gdown_raw_data)))))
    gdown_seg_rsquared[np.isinf(gdown_seg_rsquared)] = np.nan
    gdown_maxgap = max(np.diff(np.r_[smooth_dates[seg[1]],raw_dates_gdown[~np.isnan(gdown_raw_data)],smooth_dates[seg[2]]]))
    gdown_maxgap_frac = gdown_maxgap / (seg[2] - seg[1])

    return np.r_[seg_amp, seg_max, seg_int, gup_seg_rsquared, gup_numObs, gup_maxgap_frac, gdown_seg_rsquared, gdown_numObs, gdown_maxgap_frac]


def GetSegMetricsLight(seg, smooth_dates, raw_dates):
    if np.isnan(seg).any():
        return (np.nan, np.nan)
    
    # Organize greenup segment
    ######################################
    gup_raw_date_inds = np.where((raw_dates >= smooth_dates[seg[0]]) & (raw_dates <= smooth_dates[seg[1]])) # indices in raw data of gup segment
    raw_dates_gup = raw_dates[gup_raw_date_inds]

    # Organize greendown segment
    ######################################
    gdown_raw_date_inds = np.where((raw_dates >= smooth_dates[seg[1]]) & (raw_dates <= smooth_dates[seg[2]])) # indices in raw data of gdown segment
    raw_dates_gdown = raw_dates[gdown_raw_date_inds]

    ###Get the observation density for each period
    #This approach counts snow filled values as good values, since snow images are valuable for pinning down dormant period
    ###Find the biggest gap between images
    gup_maxgap = max(np.diff(np.r_[smooth_dates[seg[0]],raw_dates_gup,smooth_dates[seg[1]]]))
    gup_maxgap_frac = gup_maxgap / (seg[1] - seg[0])

    gdown_maxgap = max(np.diff(np.r_[smooth_dates[seg[1]],raw_dates_gdown,smooth_dates[seg[2]]]))
    gdown_maxgap_frac = gdown_maxgap / (seg[2] - seg[1])

    return np.r_[gup_maxgap_frac, gdown_maxgap_frac]


def AnnualMetrics(smoothed_vi, pred_dates, filled_vi, baseWeight, yr, pheno_pars):
    '''When a cycle is not detected, return a subset of metrics for the calendar year
       For now, returning evi amplitude and evi maximum
       Douglas Bolton 
    '''
    out = np.empty((pheno_pars['numLyrs'][0])); out[:] = np.nan
    try:
        inyear = np.array([date[0:4] for date in date_to_string2(pred_dates)]) == yr.astype(str) 

        vi_inyear = smoothed_vi[inyear]
        seg_min = np.nanmin(vi_inyear)
        seg_max = np.nanmax(vi_inyear)
    
        #* Subtracting 1 for accurate indexing
        out[pheno_pars['loc_numCycles'][0]-1] = 0 # Zero cycles detected
        out[pheno_pars['loc_max'][0]-1] = seg_max * 10000 # Cycle Maximum
        out[pheno_pars['loc_amp'][0]-1] = (seg_max - seg_min) * 10000 # Cycle amplitude


        # And get calendar year metrics with snow counted as good
        ind = ~np.isnan(filled_vi) & inyear
        out[pheno_pars['loc_numObs_count_snow'][0]-1] = np.sum(ind)
        out[pheno_pars['loc_maxGap_count_snow'][0]-1] = np.max(np.diff(np.r_[np.min(pred_dates[inyear]), pred_dates[ind], np.max(inyear)]))

        # Now get the full segment metrics, not counting snow and not counting gap filled
        ind = (ind) & (baseWeight == 1)   # If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation
        out[pheno_pars['loc_numObs']-1] = np.sum(ind)
        out[pheno_pars['loc_maxGap']-1] = np.max(np.diff(np.r_[np.min(pred_dates[inyear]), pred_dates[ind], np.max(inyear)]))

        return out

    except:
        print("An exception occurred")


def MakeComposites(compositeDates, pred_dates, smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6, smooth_b7):
    '''Make image composites for specific dates from a smoothed time-series
       Does all bands at the same time
       Douglas Bolton 
    '''
    days = np.where(np.isin(pred_dates, compositeDates))

    # Ensure that all days were found. If not, don't return any, as this means something is incorrect
    if len(days[0]) != len(compositeDates):
        return np.nan

    if np.isnan(smooth_b2).all():
        out = np.vstack((np.full([1, len(days[0])], np.nan), np.full([1, len(days[0])], np.nan), np.full([1, len(days[0])], np.nan), np.full([1, len(days[0])], np.nan), np.full([1, len(days[0])], np.nan), np.full([1, len(days[0])], np.nan)))
        return out

    out = np.vstack((smooth_b2[days],smooth_b3[days],smooth_b4[days],smooth_b5[days],smooth_b6[days],smooth_b7[days]))

    return out


def calculateWeights(smoothMat_Masked, numDaysFit, numYrs, pheno_pars):
    outWeights = np.zeros((numYrs, numDaysFit[0], numYrs))
    for y in range(numYrs):
        # Only compare on dates that have splined data in target year
        ind = ~np.isnan(smoothMat_Masked[:,y])
        sub_vi = smoothMat_Masked[ind,]
        numGoodDays = np.sum(~np.isnan(sub_vi), axis=0)

        # What approach to use for weighting
        #####
        eucl = np.array([np.nansum((sub_vi[:,x] - sub_vi[:,y])*(sub_vi[:,x] - sub_vi[:,y]))**0.5 for x in range(np.shape(sub_vi)[1])])

        # Now calculate euclidean distance assuming the average through the year
        # Scale euclidean distances between this value and a perfect fit (0 to 1)
        theAvg = np.full((len(sub_vi[:,y]), numYrs), np.nanmean(sub_vi[:,y]))
        theAvg[np.isnan(sub_vi)] = np.nan   # Only calculate for days that have data

        max_eucl = np.array([np.nansum((theAvg[:,x] - sub_vi[:,y])*(theAvg[:,x] - sub_vi[:,y]))**0.5 for x in range(np.shape(sub_vi)[1])])
        scaled_eucl = 1 - (eucl / max_eucl)
        scaled_eucl[scaled_eucl < 0] = 0

        # Weigh as the scaled euclidean distance (0 = same/worse than assuming average, 1 = perfect fit)
        weight = pheno_pars['maxWeight'][0] * scaled_eucl
        weight[numGoodDays < pheno_pars['minDaysForSplineComparison'][0]] = 0
        weight[np.isnan(weight)] = 0
        weight[np.isinf(weight)] = 0
        outWeights[y,:,:] = np.resize(weight, (numDaysFit[0], numYrs))
    
    return outWeights




def date_to_string(dates):
    ''' Converts floats to string date
    '''
    final_dates = []
    origin = datetime.datetime(1970, 1, 1, 0, 0)
    for date in dates:
        date = datetime.timedelta(date)
        new_date = origin + date
        final_dates.append(new_date.strftime("%Y-%m-%d"))
    
    final_dates = np.array(final_dates)

    return final_dates 

def date_to_string2(dates):
    ''' Converts floats to string date
    '''
    origin = np.datetime64('1970-01-01')
    dates = dates.astype('timedelta64[D]')

    result = origin + dates
    return np.datetime_as_string(result)


def date_to_string3(date):
    ''' Converting pheno pars dates to epoch dates
    '''
    origin = datetime.datetime(1970, 1, 1, 0, 0)
    new_date = date - origin
    return new_date.days


def DoPhenologyHLS(b2, b3, b4, b5, b6, b7, vi, snowPix, dates, imgYrs, phenYrs, splineStart, splineEnd, numDaysFit, pheno_pars, dormStart, dormEnd):
    ''' Calculate pheno metrics for each pixel
        This version using alternate years to gap fill
        Code adapted by Douglas Bolton
        Based on MODIS C6 algorithm developed by Josh Gray
    '''

    # Despike, calculate dormant value, fill snow with dormant value, despike poorly fit snow values
    try:
        spikes = CheckSpike_MultiBand(b2/10000, b4/10000, vi, dates, pheno_pars)

        # Set despiked values to NA
        vi[spikes] = np.nan
        b2[spikes] = np.nan; b3[spikes] = np.nan; b4[spikes] = np.nan; 
        b5[spikes] = np.nan; b6[spikes] = np.nan; b7[spikes] = np.nan; 

        dormIms = (dates >= dormStart) & (dates <= dormEnd)

        vi_dorm = np.nanquantile(vi[dormIms], pheno_pars['dormantQuantile'][0]) # Calc vi dormant value
        snowPix = np.array(Screen_SnowFills(vi, vi_dorm, snowPix, dates, pheno_pars), dtype=bool) # Screen poorly filled snow values

        # Now calculate dormancy values and fill individual bands
        dormObs = dormIms & (vi < vi_dorm)
        b2_dorm = np.median(b2[dormObs]); b2[snowPix] = b2_dorm
        b3_dorm = np.median(b3[dormObs]); b3[snowPix] = b3_dorm
        b4_dorm = np.median(b4[dormObs]); b4[snowPix] = b4_dorm
        b5_dorm = np.median(b5[dormObs]); b5[snowPix] = b5_dorm
        b6_dorm = np.median(b6[dormObs]); b6[snowPix] = b6_dorm
        b7_dorm = np.median(b7[dormObs]); b7[snowPix] = b7_dorm

        vi[snowPix] = vi_dorm # Fill remaining snow values with dormant value for vi



        # Determine gaps that require filling
        gDates = dates[~np.isnan(vi)]
        dDiff = np.diff(gDates) > pheno_pars['gapLengthToFill'][0] # Gaps greater than 30 days will be filled
        dStart = gDates[np.r_[dDiff, False]]
        dEnd = gDates[np.r_[False,dDiff]]


        # Locate gaps in date vector
        all_dates = np.linspace(min(splineStart), max(splineEnd), num=int(max(splineEnd)-min(splineStart)+1))

        fill_locations = np.full(len(all_dates), False)
        for d in range(len(dStart)):
            fill_locations[(all_dates >= dStart[d]) & (all_dates < dEnd[d])] = True

        fill_dates = all_dates[fill_locations]

        yToDo = np.array(range(len(imgYrs)))
        yToDo = yToDo[np.isin(imgYrs, phenYrs)]
        yrsWithGaps = np.array([]) 
        for y in yToDo:
            pred_dates = np.linspace(splineStart[y], splineEnd[y], num=int(splineEnd[y]-splineStart[y]+1))
            if np.sum(np.isin(pred_dates, fill_dates)) > 0:
                yrsWithGaps = np.r_[yrsWithGaps, imgYrs[y]]  


        #If there are gaps to be filled then we will spline all years
        #If not, just spline product years
        if len(yrsWithGaps) > 0:
            yrs = imgYrs
        else:
            yrs = phenYrs
            splineStart = splineStart[np.isin(imgYrs, phenYrs)]
            splineEnd = splineEnd[np.isin(imgYrs, phenYrs)]

        numYrs = len(yrs)
        daysVec = np.array(range(numDaysFit[0])) 





        #First, we will fit splines to each year individually
        #To line up observations from each year, we will create a matrix for vi and each band (numDaysFit x numYears)

        smoothMat = np.empty((numDaysFit[0], numYrs)); smoothMat[:] = np.nan
        maskMat = np.empty((numDaysFit[0], numYrs)); maskMat[:] = 0
        fillMat = np.copy(smoothMat)
        b2Mat = np.copy(smoothMat); b3Mat = np.copy(smoothMat); b4Mat = np.copy(smoothMat)
        b5Mat = np.copy(smoothMat); b6Mat = np.copy(smoothMat); b7Mat = np.copy(smoothMat)
        baseWeights = np.copy(maskMat)

        for y in range(numYrs):
            #Use try statement, because we don't want to stop processing if only an error in one year
            try:
                
                dateRange = (dates >= splineStart[y]) & (dates <= splineEnd[y]) & (~np.isnan(vi))
                dateSub = dates[dateRange]; viSub = vi[dateRange]; snowSub = snowPix[dateRange] 

                # Get weights
                weights = np.full(len(snowSub), 1.0) 
                weights[snowSub == 1] = pheno_pars['snowWeight'][0] 

                pred_dates = np.linspace(splineStart[y], splineEnd[y], num=int(splineEnd[y]-splineStart[y]+1)) 
                # Assign weights and run cubic spline
                smoothed = Smooth_VI3(viSub, dateSub, pred_dates, weights, vi_dorm) 

                # Mask spline in gaps, and before/after first/last image
                maskMat[fill_locations[np.isin(all_dates, pred_dates)],y] = 1    # Mask spline in gaps
                maskMat[pred_dates < dateSub[0],y] = 1                           # Mask spline before first image
                maskMat[pred_dates > dateSub[len(dateSub)-1],y] = 1

                # Mask spline in the buffer years (only interested in comparing splines in target year)

                maskMat[np.array([date[0:4] for date in date_to_string2(pred_dates)]) != np.str_(yrs[y]),y] = 1
                fillDs = np.isin(pred_dates, dateSub)

                smoothMat[:,y] = smoothed
                baseWeights[fillDs, y] = weights
                fillMat[fillDs,y] = viSub
                b2Mat[fillDs, y] = b2[dateRange]; b3Mat[fillDs, y] = b3[dateRange]; b4Mat[fillDs, y] = b4[dateRange]
                b5Mat[fillDs, y] = b5[dateRange]; b6Mat[fillDs, y] = b6[dateRange]; b7Mat[fillDs, y] = b7[dateRange]

            except:
                print("Error!")

        xs = np.resize(daysVec, len(daysVec)*numYrs)
        ys = fillMat.flatten('F') #! vecLength
        ysGood = ~np.isnan(ys)
        baseW = baseWeights.flatten('F')   # Base Weights are 1=clear observation, 0.5=snow-filled

        b2Mat = b2Mat.flatten('F');b3Mat = b3Mat.flatten('F')
        b4Mat = b4Mat.flatten('F');b5Mat = b5Mat.flatten('F')
        b6Mat = b6Mat.flatten('F');b7Mat = b7Mat.flatten('F')

        smoothMat_Masked = np.copy(smoothMat)
        smoothMat_Masked[maskMat.astype(bool)] = np.nan 


        # Loop through years, compare spline to other years, weight each year based on similarity, fit spline, calculate phenology
        # Just product years now
        yToDo = np.array(range(numYrs))
        yToDo = yToDo[np.isin(yrs, phenYrs)]


        weightArray = calculateWeights(smoothMat_Masked, numDaysFit, numYrs, pheno_pars)

        prevYear = daysVec <= pheno_pars['splineBuffer'][0]-1
        inYear = (daysVec > pheno_pars['splineBuffer'][0]-1) & (daysVec <= (pheno_pars['splineBuffer'][0]+364))
        nextYear = (daysVec > (pheno_pars['splineBuffer']+364))
    except:
        #If there is an error despiking or other initial steps, return NAs
        end_result = np.empty((pheno_pars['numLyrs'][0]*len(phenYrs),1)); end_result.fill(np.nan) 
        return end_result.reshape(1, pheno_pars['numLyrs'][0]*len(phenYrs)) # []

    outAll = np.array([])
    for y in yToDo:
        try:
            pred_dates = np.linspace(splineStart[y], splineEnd[y], num=int(splineEnd[y]-splineStart[y]+1))

            if np.isin(yrs[y], yrsWithGaps):

                indPrev = np.array(y-1); indPrev[indPrev<0] = 0 
                indNext = np.array(y+1); indNext[indNext>numYrs-1] = numYrs-1 
                weights = np.concatenate([weightArray[indPrev, prevYear==True], weightArray[y, inYear==True], weightArray[indNext, nextYear==True]])

                # Where are the gaps?
                toFill = fill_locations[np.isin(all_dates, pred_dates)]

                weights[~toFill,] = 0   # Set weight to zero for observations that aren't in a gap
                weights[:,y] = 1        # Set weights in target year to 1
    

                # Now that we have weights, calculate phenology
                #######################
                weights = weights.flatten('F') * baseW  # Multiple weights by base weight (1=good,0.5=snow-filled)
                theInds = ysGood & (weights > 0)
                xs_sub = xs[theInds]; w_sub = weights[theInds]
               
                smoothed_vi = Smooth_VI4(ys[theInds], xs_sub, daysVec, w_sub, vi_dorm) 

            else:

                # Variables needed for next steps if the above gap filling was not done
                theInds = np.full(len(ysGood), False)
                theInds[((y-1+1)*numDaysFit[0]):((y+1)*numDaysFit[0])] = True
                xs_sub = xs[theInds]; w_sub = baseW[theInds]

                smoothed_vi = smoothMat[:,y] # If no gaps to fill, just use existing spline 
            
            # Fit phenology
            peaks = FindPeaks(smoothed_vi)
            if np.isnan(peaks).all():
                outAll = np.r_[outAll, AnnualMetrics(smoothed_vi, pred_dates, fillMat[:,y], baseWeights[:,y], yrs[y], pheno_pars)] #If no peaks, calc annual metrics, and move to next year
                continue

            #Find full segments
            full_segs = GetSegs(peaks, smoothed_vi, pheno_pars)
            # No Null in python so just removing [-9999, -9999, -9999]

            if np.array_equal(full_segs, [-9999, -9999, -9999]):
                outAll = np.r_[outAll, AnnualMetrics(smoothed_vi, pred_dates, fillMat[:,y], baseWeights[:,y], yrs[y], pheno_pars)] #If no valid segments, calc annual metrics, and move to next year
                continue

            if -9999 in full_segs:
                full_segs = np.array(full_segs[0:-1], dtype=int)

            #Only keep segments with peaks within year *****
            full_segs = full_segs[inYear[[full_segs[i][1] for i in range(len(full_segs))]]]

            if len(full_segs) == 0:
                outAll = np.r_[outAll, AnnualMetrics(smoothed_vi, pred_dates, fillMat[:,y], baseWeights[:,y], yrs[y], pheno_pars)] #If no valid peaks, calc annual metrics, and move to next year
                continue

            #Get PhenoDates
            pheno_dates = GetPhenoDates(full_segs, smoothed_vi, pred_dates, pheno_pars)
            phen = pheno_dates #* In R we unlist but no need here because it's already on array
            if np.isnan(phen).all():
                outAll = np.r_[outAll, AnnualMetrics(smoothed_vi, pred_dates, fillMat[:,y], baseWeights[:,y], yrs[y], pheno_pars)] #If no dates detected, calc annual metrics, and move to next year
                continue

            # Get metrics that describe the segments and the year


            # First, get metrics counting gap filled observations as "good" observations
            seg_metricsFill = np.array([GetSegMetricsLight(full_segs[i], daysVec, np.sort(xs_sub)) for i in range(len(full_segs))])
            un = np.concatenate(seg_metricsFill)
            ln = len(un)
            gup_maxgap_frac_filled = un[np.arange(0, ln, 2)] * 100
            gdown_maxgap_frac_filled = un[np.arange(1, ln, 2)] * 100


            # Second, get segment metrics with snow observations counted as "good" observations
            filled_vi = fillMat[:,y]
            seg_metricsFill = np.array([GetSegMetricsLight(full_segs[i], daysVec, daysVec[~np.isnan(filled_vi)]) for i in range(len(full_segs))])
            un = np.concatenate(seg_metricsFill)
            ln = len(un)
            gup_maxgap_frac_count_snow = un[np.arange(0, ln, 2)] * 100
            gdown_maxgap_frac_count_snow = un[np.arange(1, ln, 2)] * 100

            #And get calendar year metrics with snow counted as good
            numObs_count_snow = np.sum(~np.isnan(filled_vi) & inYear)
            maxGap_annual_count_snow = max(np.diff(np.r_[pheno_pars['splineBuffer'], daysVec[~np.isnan(filled_vi) & inYear], 364+pheno_pars['splineBuffer']])) # Got rid of +1 in python


            # Now get the full segment metrics, not counting snow and not counting gap filled
            filled_vi[baseWeights[:,y] < 1] = np.nan # If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. So set to NA.
            numObs = np.sum(~np.isnan(filled_vi) & inYear) # Number of observations in year
            maxGap_annual = max(np.diff(np.r_[pheno_pars['splineBuffer'], daysVec[~np.isnan(filled_vi) & inYear], 365+pheno_pars['splineBuffer']])) # Max gap (in days) during year (Got rid of +1 in python)
            seg_metrics = np.array([GetSegMetrics(full_segs[i], smoothed_vi, filled_vi[~np.isnan(filled_vi)], pred_dates, pred_dates[~np.isnan(filled_vi)]) for i in range(len(full_segs))]) # Full segment metrics


            # Unlist and scale the seg metrics
            un = np.concatenate(seg_metrics)
            ln = len(un)
            seg_amp = un[np.arange(0, ln, 9)] * 10000
            seg_max = un[np.arange(1, ln, 9)] * 10000
            seg_int = un[np.arange(2, ln, 9)] * 100
            gup_rsq = un[np.arange(3, ln, 9)] * 10000
            gup_numObs = un[np.arange(4, ln, 9)]
            gup_maxgap_frac = un[np.arange(5, ln, 9)] * 100
            gdown_rsq = un[np.arange(6, ln, 9)] * 10000
            gdown_numObs = un[np.arange(7, ln, 9)]
            gdown_maxgap_frac = un[np.arange(8, ln, 9)] * 100



            numRecords = len(seg_amp) # How many cycles were recorded

            naCheck = np.isnan(seg_amp)
            numCyc = np.sum(naCheck == 0) #How many cycles have good data (seg metrics has valid observations)

            # If no cycles have good data, record Na output and move to next
            if numCyc == 0:
                outAll = np.r_[outAll, AnnualMetrics(smoothed_vi, pred_dates, fillMat[:,y], baseWeights[:,y], yrs[y], pheno_pars)]
                continue


            # Fit spline to individual bands. Waiting until now so that we don't fit for pixels with no cycle
            # Since this takes the most time, there is an option to skip this and just fill NAs
            # if pheno_pars['doComposites'][0]:
            #     smooth_b2 = Smooth_Bands(b2Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
            #     smooth_b3 = Smooth_Bands(b3Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
            #     smooth_b4 = Smooth_Bands(b4Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
            #     smooth_b5 = Smooth_Bands(b5Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
            #     smooth_b6 = Smooth_Bands(b6Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
            #     smooth_b7 = Smooth_Bands(b7Mat[theInds], xs_sub, daysVec, w_sub, pheno_pars)
            # else: 
                smooth_b2 = np.array(np.nan);smooth_b3 = np.array(np.nan);smooth_b4 = np.array(np.nan)
                smooth_b5 = np.array(np.nan);smooth_b6 = np.array(np.nan);smooth_b7 = np.array(np.nan)

            if numRecords == 1:
                comp = (MakeComposites(phen,pred_dates,smooth_b2,smooth_b3,smooth_b4,smooth_b5,smooth_b6,smooth_b7)).flatten('F')

                # If only one cycle was recorded, report it, and fill NA values for second cycle
                out = np.r_[1, phen, seg_max, seg_amp, seg_int, comp, gup_numObs, gup_maxgap_frac, gup_maxgap_frac_count_snow, gup_maxgap_frac_filled, gup_rsq, gdown_numObs, gdown_maxgap_frac, gdown_maxgap_frac_count_snow, gdown_maxgap_frac_filled, gdown_rsq, np.full([1, pheno_pars['numLyrsCycle2'][0]], np.nan)[0], numObs, numObs_count_snow, maxGap_annual, maxGap_annual_count_snow]

            else:
                # If there are multiple cycles, sort by amplitude and report two highest amplitudes (highest amplitude first)
                theOrd = np.argsort(seg_amp)

                phen1 = phen[np.arange(theOrd[0], len(phen), numRecords)]
                phen2 = phen[np.arange(theOrd[1], len(phen), numRecords)]

                comp1 = (MakeComposites(phen1,pred_dates,smooth_b2,smooth_b3,smooth_b4,smooth_b5,smooth_b6,smooth_b7)).flatten('F')
                comp1 = (MakeComposites(phen2,pred_dates,smooth_b2,smooth_b3,smooth_b4,smooth_b5,smooth_b6,smooth_b7)).flatten('F')

                if naCheck[theOrd[1]]:
                    # If the second cycle did not have enough observations (seg_metrics = NA), only report the first cycle
                    out = np.r_[numCyc, phen1, seg_max[theOrd[0]], seg_amp[theOrd[0]], seg_int[theOrd[0]], comp1, gup_numObs[theOrd[0]], gup_maxgap_frac[theOrd[0]], gup_maxgap_frac_count_snow[theOrd[0]], gup_maxgap_frac_filled[theOrd[0]], gup_rsq[theOrd[0]], gdown_numObs[theOrd[0]], gdown_maxgap_frac[theOrd[0]], gdown_maxgap_frac_count_snow[theOrd[0]], gdown_maxgap_frac_filled[theOrd[0]], gdown_rsq[theOrd[0]],numObs, numObs_count_snow, maxGap_annual, maxGap_annual_count_snow]

                else: 
                    # If the second cycle had enough observations, report both
                    out = np.r_[numCyc, phen1, seg_max[theOrd[0]], seg_amp[theOrd[0]], seg_int[theOrd[0]], comp1, gup_numObs[theOrd[0]], gup_maxgap_frac[theOrd[0]], gup_maxgap_frac_count_snow[theOrd[0]], gup_maxgap_frac_filled[theOrd[0]], gup_rsq[theOrd[0]], gdown_numObs[theOrd[0]], gdown_maxgap_frac[theOrd[0]], gdown_maxgap_frac_count_snow[theOrd[0]], gdown_maxgap_frac_filled[theOrd[0]], gdown_rsq[theOrd[0]], phen2, seg_max[theOrd[1]], seg_amp[theOrd[1]], seg_int[theOrd[1]], comp1, gup_numObs[theOrd[1]], gup_maxgap_frac[theOrd[1]], gup_maxgap_frac_count_snow[theOrd[1]], gup_maxgap_frac_filled[theOrd[1]], gup_rsq[theOrd[1]], gdown_numObs[theOrd[1]], gdown_maxgap_frac[theOrd[1]], gdown_maxgap_frac_count_snow[theOrd[1]], gdown_maxgap_frac_filled[theOrd[1]], gdown_rsq[theOrd[1]], numObs, numObs_count_snow, maxGap_annual, maxGap_annual_count_snow]
            outAll = np.r_[outAll, out]; del(out)

        except:
            outAll = np.r_[outAll, np.full([1, pheno_pars['numLyrs'][0]], np.nan)[0]]
            return outAll.reshape(1, len(outAll))

    return outAll.reshape(1, len(outAll))





