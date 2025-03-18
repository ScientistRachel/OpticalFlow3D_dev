"""
calc_flow contains the functions necessary to calculate optical flow in 2D and
3D, as well as a function for parsing files and their metadata.
"""

import math
import numpy as np
from scipy.ndimage import correlate1d

def calc_flow2D(images,xySig=3,tSig=1,wSig=4):
    """
    The calc_flow2D function calculates optical flow velocities for a single
    z-slice over time. Surrounding images in time are necessary to perform
    the calculations. To peform calculations on an entire timelapse, see
    parse_flow

    This script uses the convention that (0,0) is located in the upper-left
    corner of an image. This is inline with conventions used in other
    programs (e.g., ImageJ/FIJI), but note that it means that positive
    y-velocities point down, which can be non-intuitive in some cases.

    ARGS:
    images = 3D matrix with dimensions N_T, N_Y, N_X
            N_T should be odd as only the central timepoint will be analyzed.
            N_T must be greater than or equal to 3*tSig+1.
    xySig  = sigma value for smoothing in all spatial dimensions. Default 3.
            Larger values remove noise but remove spatial detail.
    tSig   = sigma value for smoothing in the temporal dimension. Default 1.
            Larger values remove noise but remove temporal detail.
    wSig   = sigma value for Lucas-Kanade neighborhood. Default is 4.
            Larger values include a larger neighboorhood in the
            Lucas-Kanade constraint and will smooth over small features.

    RETURNS:
    vx     = Velocity in the x direction, reported as pixels/frame
    vy     = Velocity in the y direction, reported as pixels/frame
    rel    = Reliability, the smallest eigenvalue of (A'wA)
            This is a measure of confidence in the linear algebra solution
            and can be used to mask the velocities for downstream analysis.
    """

    ### Check the function inputs ################################################
    # Check that the images are 2D + time
    if not(len(images.shape)==3): 
        exit('ERROR: Input image must be a 3D matrix with dimensions N_T, N_Y, N_X')
    # Check image size against tSig
    Nt = images.shape[0] 
    if Nt < 3*tSig+1:
        exit('ERROR: Input images will lead to edge effects. N_T must be >= 3*tSig+1')
    # Check for an odd number of frames
    if not(Nt % 2):
        exit('ERROR: Input images must have an odd number of timepoints. Only the central time point is analyzed')
    NtSlice = math.ceil(Nt/2)-1 # -1 because python indexing starts from 0


    ### Set up filters ###########################################################
    # Common terms
    x = np.arange(-math.ceil(3*xySig),math.ceil(3*xySig)+1)
    xySig2 = xySig/4
    y = np.arange(-math.ceil(3*xySig2),math.ceil(3*xySig2)+1)
    fderiv = np.exp(-x*x/2/xySig/xySig)/math.sqrt(2*math.pi)/xySig
    fsmooth = np.exp(-y*y/2/xySig2/xySig2)/math.sqrt(2*math.pi)/xySig2
    gderiv = x/xySig/xySig
    gsmooth = 1

    # Build y-gradient filter kernels (along first spatial dimension)
    yFil1 = (fderiv*gderiv)
    xFil1 = (fsmooth*gsmooth)
    # Build x-gradient filter kernels (along second spatial dimension)
    yFil2 = (fsmooth*gsmooth)
    xFil2 = (fderiv*gderiv)

    # Build t-gradient filter kernels (t = third dimension)
    t = np.arange(-math.ceil(3*tSig),math.ceil(3*tSig)+1)
    fx = np.exp(-x*x/2/xySig/xySig)/math.sqrt(2*math.pi)/xySig
    ft = np.exp(-t*t/2/tSig/tSig)/math.sqrt(2*math.pi)/tSig
    gx = 1
    gt = t/tSig/tSig
    yFil3 = (fx*gx)
    xFil3 = yFil3
    tFil3 = ft*gt

    # Structure tensor -- Lucas Kanade neighborhood filter
    wRange = np.arange(-math.ceil(3*wSig),math.ceil(3*wSig))
    gw = np.exp(-wRange*wRange/2/wSig/wSig)/math.sqrt(2*math.pi)/wSig
    yFil4 = gw
    xFil4 = gw

    # Throughout will use del to keep the memory clear as this processing is memory intensive
    del gderiv, gsmooth, gt, gw, gx, ft, fx, fsmooth, fderiv, x, y, t, wRange


    ### Spatial and Temporal Gradients ###########################################
    # Spatial gradients require at least N_T = 3 to avoid edge effets, while 
    # the temporal gradient requires N_T >= 3*tSig+1.
    # After gradients are calculated, only the central 3 slices are kept
    dyI = correlate1d(correlate1d(images, yFil1, axis=1, mode='reflect'), xFil1, axis=2, mode='reflect')
    dyI = dyI[NtSlice-1:NtSlice+1,:]
    del xFil1, yFil1

    dxI = correlate1d(correlate1d(images, yFil2, axis=1, mode='reflect'), xFil2, axis=2, mode='reflect')
    dxI = dxI[NtSlice-1:NtSlice+1,:]
    del xFil2, yFil2

    dtI = correlate1d(correlate1d(correlate1d(images, yFil3, axis=1, mode='reflect'), xFil3, axis=2, mode='reflect'), tFil3, axis=0, mode='reflect')
    dtI = dtI[NtSlice-1:NtSlice+1,:]
    del xFil3, yFil3, tFil3

    del images


    ### Structure Tensor Inputs ##################################################
    # The following calculations are for the individual elements of the
    # matrices required for the optical flow calculation, incorporating
    # Gaussian weighting into the Lucas-Kanade constraint.
    # Once filtering is done, only the central timepoint is kept to save
    # memory. Because all outputs of the above section are N_X x N_Y x 3,
    # the central slice = 1.

    # Time components
    wdtx = correlate1d(correlate1d(dxI*dtI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdtx = wdtx[1,:]
    wdty = correlate1d(correlate1d(dyI*dtI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdty = wdty[1,:]
    del dtI

    # Spatial Components
    wdxy = correlate1d(correlate1d(dxI*dyI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdxy = wdxy[1,:]
    wdx2 = correlate1d(correlate1d(dxI*dxI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdx2 = wdx2[1,:]
    del dxI
    wdy2 = correlate1d(correlate1d(dyI*dyI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdy2 = wdy2[1,:]
    del dyI
    del xFil4, yFil4


    ### Optical Flow Calculations ################################################
    # Equation is v = (A' w A)^-1 A' w b
    # A = -[dxI dyI]
    # b = [dtI]
    # w multiplication is incorporated in the structure tensor inputs above
    # A' w b = -[wdtx wdty]  (minus sign because of negative sign on A)
    # (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]
    # A^-1 = [a b ; c d]^-1 = (1/det(A))[d - b; -c a]
    determinant = (wdx2*wdy2) - (wdxy*wdxy)
    vx = ((determinant+np.finfo(float).eps)**-1)*((wdy2*-wdtx)+(-wdxy*-wdty))
    vy = ((determinant+np.finfo(float).eps)**-1)*((-wdxy*-wdtx)+(wdx2*-wdty))
    del wdtx, wdty, wdxy


    ### Eigenvalues for Reliability ################################################
    # solve det(A^T w A - lamda I) = 0
    # (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]
    trace = wdx2 + wdy2
    del wdx2, wdy2

    L1 = (trace + np.sqrt(trace**2 - 4*determinant))/2
    L2 = (trace - np.sqrt(trace**2 - 4*determinant))/2
    rel = np.real(np.minimum(L1,L2))
    del L1, L2

    return vx, vy, rel


    """
    The calc_flow2D function calculates optical flow velocities for a single
    z-slice over time. Surrounding images in time are necessary to perform
    the calculations. To peform calculations on an entire timelapse, see
    parse_flow

    This script uses the convention that (0,0) is located in the upper-left
    corner of an image. This is inline with conventions used in other
    programs (e.g., ImageJ/FIJI), but note that it means that positive
    y-velocities point down, which can be non-intuitive in some cases.

    ARGS:
    images = 3D matrix with dimensions N_T, N_Y, N_X
            N_T should be odd as only the central timepoint will be analyzed.
            N_T must be greater than or equal to 3*tSig+1.
    xyzSig  = sigma value for smoothing in all spatial dimensions. Default 3.
            Larger values remove noise but remove spatial detail.
    tSig   = sigma value for smoothing in the temporal dimension. Default 1.
            Larger values remove noise but remove temporal detail.
    wSig   = sigma value for Lucas-Kanade neighborhood. Default is 4.
            Larger values include a larger neighboorhood in the
            Lucas-Kanade constraint and will smooth over small features.

    RETURNS:
    vx     = Velocity in the x direction, reported as pixels/frame
    vy     = Velocity in the y direction, reported as pixels/frame
    rel    = Reliability, the smallest eigenvalue of (A'wA)
            This is a measure of confidence in the linear algebra solution
            and can be used to mask the velocities for downstream analysis.
    """

    ### Check the function inputs ################################################
    # Check that the images are 2D + time
    if not(len(images.shape)==3): 
        exit('ERROR: Input image must be a 3D matrix with dimensions N_T, N_Y, N_X')
    # Check image size against tSig
    Nt = images.shape[0] 
    if Nt < 3*tSig+1:
        exit('ERROR: Input images will lead to edge effects. N_T must be >= 3*tSig+1')
    # Check for an odd number of frames
    if not(Nt % 2):
        exit('ERROR: Input images must have an odd number of timepoints. Only the central time point is analyzed')
    NtSlice = math.ceil(Nt/2)-1 # -1 because python indexing starts from 0


    ### Set up filters ###########################################################
    # Common terms
    x = np.arange(-math.ceil(3*xyzSig),math.ceil(3*xyzSig)+1)
    xyzSig2 = xyzSig/4
    y = np.arange(-math.ceil(3*xyzSig2),math.ceil(3*xyzSig2)+1)
    fderiv = np.exp(-x*x/2/xyzSig/xyzSig)/math.sqrt(2*math.pi)/xyzSig
    fsmooth = np.exp(-y*y/2/xyzSig2/xyzSig2)/math.sqrt(2*math.pi)/xyzSig2
    gderiv = x/xyzSig/xyzSig
    gsmooth = 1

    # Build y-gradient filter kernels (along first spatial dimension)
    yFil1 = (fderiv*gderiv)
    xFil1 = (fsmooth*gsmooth)
    # Build x-gradient filter kernels (along second spatial dimension)
    yFil2 = (fsmooth*gsmooth)
    xFil2 = (fderiv*gderiv)

    # Build t-gradient filter kernels (t = third dimension)
    t = np.arange(-math.ceil(3*tSig),math.ceil(3*tSig)+1)
    fx = np.exp(-x*x/2/xyzSig/xyzSig)/math.sqrt(2*math.pi)/xyzSig
    ft = np.exp(-t*t/2/tSig/tSig)/math.sqrt(2*math.pi)/tSig
    gx = 1
    gt = t/tSig/tSig
    yFil3 = (fx*gx)
    xFil3 = yFil3
    tFil3 = ft*gt

    # Structure tensor -- Lucas Kanade neighborhood filter
    wRange = np.arange(-math.ceil(3*wSig),math.ceil(3*wSig))
    gw = np.exp(-wRange*wRange/2/wSig/wSig)/math.sqrt(2*math.pi)/wSig
    yFil4 = gw
    xFil4 = gw

    # Throughout will use del to keep the memory clear as this processing is memory intensive
    del gderiv, gsmooth, gt, gw, gx, ft, fx, fsmooth, fderiv, x, y, t, wRange


    ### Spatial and Temporal Gradients ###########################################
    # Spatial gradients require at least N_T = 3 to avoid edge effets, while 
    # the temporal gradient requires N_T >= 3*tSig+1.
    # After gradients are calculated, only the central 3 slices are kept
    dyI = correlate1d(correlate1d(images, yFil1, axis=1, mode='reflect'), xFil1, axis=2, mode='reflect')
    dyI = dyI[NtSlice-1:NtSlice+1,:]
    del xFil1, yFil1

    dxI = correlate1d(correlate1d(images, yFil2, axis=1, mode='reflect'), xFil2, axis=2, mode='reflect')
    dxI = dxI[NtSlice-1:NtSlice+1,:]
    del xFil2, yFil2

    dtI = correlate1d(correlate1d(correlate1d(images, yFil3, axis=1, mode='reflect'), xFil3, axis=2, mode='reflect'), tFil3, axis=0, mode='reflect')
    dtI = dtI[NtSlice-1:NtSlice+1,:]
    del xFil3, yFil3, tFil3

    del images


    ### Structure Tensor Inputs ##################################################
    # The following calculations are for the individual elements of the
    # matrices required for the optical flow calculation, incorporating
    # Gaussian weighting into the Lucas-Kanade constraint.
    # Once filtering is done, only the central timepoint is kept to save
    # memory. Because all outputs of the above section are N_X x N_Y x 3,
    # the central slice = 1.

    # Time components
    wdtx = correlate1d(correlate1d(dxI*dtI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdtx = wdtx[1,:]
    wdty = correlate1d(correlate1d(dyI*dtI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdty = wdty[1,:]
    del dtI

    # Spatial Components
    wdxy = correlate1d(correlate1d(dxI*dyI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdxy = wdxy[1,:]
    wdx2 = correlate1d(correlate1d(dxI*dxI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdx2 = wdx2[1,:]
    del dxI
    wdy2 = correlate1d(correlate1d(dyI*dyI, yFil4, axis=1, mode='reflect'), xFil4, axis=2, mode='reflect')
    wdy2 = wdy2[1,:]
    del dyI
    del xFil4, yFil4


    ### Optical Flow Calculations ################################################
    # Equation is v = (A' w A)^-1 A' w b
    # A = -[dxI dyI]
    # b = [dtI]
    # w multiplication is incorporated in the structure tensor inputs above
    # A' w b = -[wdtx wdty]  (minus sign because of negative sign on A)
    # (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]
    # A^-1 = [a b ; c d]^-1 = (1/det(A))[d - b; -c a]
    determinant = (wdx2*wdy2) - (wdxy*wdxy)
    vx = ((determinant+np.finfo(float).eps)**-1)*((wdy2*-wdtx)+(-wdxy*-wdty))
    vy = ((determinant+np.finfo(float).eps)**-1)*((-wdxy*-wdtx)+(wdx2*-wdty))
    del wdtx, wdty, wdxy


    ### Eigenvalues for Reliability ################################################
    # solve det(A^T w A - lamda I) = 0
    # (A' w A) = [a=wdx2 b=wdxy ; c=wdxy d=wdy2]
    trace = wdx2 + wdy2
    del wdx2, wdy2

    L1 = (trace + np.sqrt(trace**2 - 4*determinant))/2
    L2 = (trace - np.sqrt(trace**2 - 4*determinant))/2
    rel = np.real(np.minimum(L1,L2))
    del L1, L2

    return vx, vy, rel