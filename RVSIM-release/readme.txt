RVSIM Software release.

========================================================================

-----------COPYRIGHT NOTICE ENDS WITH THIS LINE------------%

Author  : Guangyi Yang, Deshi Li, Fan Lu, Yue Liao and Wen Yang
Version : 1.0

The authors are with the School of Electronic Information, Wuhan University, 
Wuhan 430072, China.

Kindly report any suggestions or corrections to ygy@whu.edu.cn

========================================================================

This is a demonstration of the Riesz transform and Visual contrast sensitivity-based feature SIMilarity (RVSIM) index.
The algorithm is described in:

Yang G, Li D, Lu F, et al. RVSIM: a feature similarity method for full-reference image quality assessment[J]. EURASIP Journal on Image and Video Processing, 2018, 2018(1): 6.


You can change this program as you like and use it anywhere, but please
refer to its original source (cite our paper and our web page at
https://sites.google.com/site/jacobygy/).

========================================================================

Running on Matlab 


Input : Two test images in grayscale form

Output: An index indicating the similarity between the 2 test images. The index ranges between 0 and 1 (0 represents the lowest similarity, 1 the highest).
  
Usage:
    RVSIM_index = RVSIM(img1, img2);
    
    The function takes 2 grayscale images img1 and img2 as input. Besides
    it requires that many values of visual Contrast Sensitivity
    Function(CSF) to accomplish band weighting. In file
    "RVSIM_csf-8bands-m2.1min3.mat" we have already stored a set of CSF 
    values for a log-gabor filter bank whose arguments are shown below.

    nscale          = 5;
    norient         = 1;
    minWaveLength   = 3;
    mult            = 2.1;
    sigmaOnf        = 0.55;
    dThetaOnSigma   = 1.5;
    
    If any of the arguments(except nscale) should be changed, the band
    distribution of the filter bank might be changed and hence a new set
    of CSF values would be necessary in matching with the new filter bank.

Dependencies: 

Data files: RVSIM_csf-5bands-m2.1min3.mat (provided with release)

Image Files: testimage1.bmp, testimage2.bmp

========================================================================
