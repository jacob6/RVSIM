# RVSIM

This is the official repository of [*RVSIM: a feature similarity method for full-reference image quality assessment*](https://jivp-eurasipjournals.springeropen.com/articles/10.1186/s13640-018-0246-1). The paper is published on the open access [*EURASIP Journal on Image and Video Processing*](https://jivp-eurasipjournals.springeropen.com/) and the full text is available [online](https://jivp-eurasipjournals.springeropen.com/track/pdf/10.1186/s13640-018-0246-1). 

### Citation
If you find our work useful for your research, we kindly suggest you to cite 

> Yang G, Li D, Lu F, et al. RVSIM: a feature similarity method for full-reference image quality assessment[J]. EURASIP Journal on Image and Video Processing, 2018, 2018(1): 6.

```
@Article{Yang2018,
author="Yang, Guangyi
and Li, Deshi
and Lu, Fan
and Liao, Yue
and Yang, Wen",
title="RVSIM: a feature similarity method for full-reference image quality assessment",
journal="EURASIP Journal on Image and Video Processing",
year="2018",
month="Jan",
day="19",
volume="2018",
number="1",
pages="6",
abstract="Image quality assessment is an important topic in the field of digital image processing. In this study, a full-reference image quality assessment method called Riesz transform and Visual contrast sensitivity-based feature SIMilarity index (RVSIM) is proposed. More precisely, a Log-Gabor filter is first used to decompose reference and distorted images, and Riesz transform is performed on the decomposed images on the basis of monogenic signal theory. Then, the monogenic signal similarity matrix is obtained by calculating the similarity of the local amplitude/phase/direction characteristics of monogenic signal. Next, we weight the summation of these characteristics with visual contrast sensitivity. Since the first-order Riesz transform cannot clearly express the corners and intersection points in the image, we calculate the gradient magnitude similarity between the reference and distorted images as a feature, which is combined with monogenic signal similarity to obtain a local quality map. Finally, we conduct the monogenic phase congruency using the Riesz transform feature matrix from the reference image and utilize it as a weighted function to derive the similarity index. Extensive experiments on five benchmark IQA databases, namely, LIVE, CSIQ, TID2008, TID2013, and Waterloo Exploration, indicate that RVSIM is a robust IQA method.",
issn="1687-5281",
doi="10.1186/s13640-018-0246-1",
url="https://doi.org/10.1186/s13640-018-0246-1"
}
```

### Instruction
This is an implementation of *the Riesz transform and Visual contrast sensitivity-based feature SIMilarity* (`RVSIM`) index. Given a pair of original and distorted images in grayscale, this method estimates the image quality of the bad one. The output index ranges between `0` and `1`, with a lower score represents a poorer image quality. All the `MATLAB` code provided is fully tested on `Windows 7 x64` and `MATLAB R2016a`. 

## Author
> Guangyi Yang, Deshi Li, Fan Lu, Yue Liao and Wen Yang  
> from the School of Electronic Information, Wuhan University, 
Wuhan 430072, China.

Kindly report any suggestions or corrections to `ygy@whu.edu.cn`

## Performance
`Fig.7` taken from the [paper](https://jivp-eurasipjournals.springeropen.com/articles/10.1186/s13640-018-0246-1) to demonstrate the performance of the `RVSIM` method. As is seen, the `RVSIM` index is rather in line with `DMOS`, which indicates good consistency with the subjective perception of HVS. 

![fig7](https://github.com/jacob6/RVSIM/blob/master/fig7.png)

## Usage
for detailed usage and required dependencies, please view `RVSIM_realease/readme.txt`

