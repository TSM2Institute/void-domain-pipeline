J/ApJ/785/104    redMaPPer. I. Algorithm applied to SDSS DR8    (Rykoff+, 2014)
================================================================================
redMaPPer. I. Algorithm and SDSS DR8 catalog.
    Rykoff E.S., Rozo E., Busha M.T., Cunha C.E., Finoguenov A., Evrard A.,
    Hao J., Koester B.P., Leauthaud A., Nord B., Pierre M., Reddick R.,
    Sadibekova T., Sheldon E.S., Wechsler R.H.
   <Astrophys. J., 785, 104 (2014)>
   =2014ApJ...785..104R    (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: Clusters, galaxy ; Photometry ; Models ; Redshifts  ; Surveys
Keywords: galaxies: clusters: general

Abstract:
    We describe redMaPPer, a new red sequence cluster finder specifically
    designed to make optimal use of ongoing and near-future large
    photometric surveys. The algorithm has multiple attractive features:
    (1) it can iteratively self-train the red sequence model based on a
    minimal spectroscopic training sample, an important feature for
    high-redshift surveys. (2) It can handle complex masks with varying
    depth. (3) It produces cluster-appropriate random points to enable
    large-scale structure studies. (4) All clusters are assigned a full
    redshift probability distribution P(z). (5) Similarly, clusters can
    have multiple candidate central galaxies, each with corresponding
    centering probabilities. (6) The algorithm is parallel and numerically
    efficient: it can run a Dark Energy Survey-like catalog in ~500 CPU
    hours. (7) The algorithm exhibits excellent photometric redshift
    performance, the richness estimates are tightly correlated with
    external mass proxies, and the completeness and purity of the
    corresponding catalogs are superb. We apply the redMaPPer algorithm to
    ~10000deg^2^ of SDSS DR8 data and present the resulting catalog of
    ~25000 clusters over the redshift range z{isin}[0.08,0.55]. The
    redMaPPer photometric redshifts are nearly Gaussian, with a scatter
    {sigma}_z_~0.006 at z~0.1, increasing to {sigma}_z_~0.02 at z~0.5 due
    to increased photometric noise near the survey limit. The median value
    for |{Delta}z|/(1+z) for the full sample is 0.006. The incidence of
    projection effects is low (<= 5%). Detailed performance comparisons of
    the redMaPPer DR8 cluster catalog to X-ray and Sunyaev-Zel'dovich
    catalogs are presented in a companion paper.

Description:
    The redMaPPer algorithm is designed to handle an arbitrary photometric
    galaxy catalog, with an arbitrary number of photometric bands (>=3).
    Of course, the quality of the output depends on the quality of the
    photometry. As a case study, in this paper we run redMaPPer on SDSS
    DR8 data (Aihara et al. 2011ApJS..193...29A), due to its large area
    and uniform coverage.

    Although our cluster finder uses only photometric data, we require
    spectroscopic data to calibrate the red sequence and to validate our
    photometric redshifts. For this purpose we use the SDSS DR9
    spectroscopic catalog (Ahn et al. 2012, V/139).

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
table1.dat       829    25325   redMaPPer DR8 cluster catalog
table2.dat       153  1318360   redMaPPer DR8 member catalog
--------------------------------------------------------------------------------

See also:
 V/139   : The SDSS Photometric Catalog, Release 9 (Adelman-McCarthy+, 2012)
 II/306  : The SDSS Photometric Catalog, Release 8 (Adelman-McCarthy+, 2011)
 VII/110 : Rich Clusters of Galaxies (Abell+ 1989)
 J/ApJS/224/1    : redMaPPer cluster catalog from DES data (Rykoff+, 2016)
 J/ApJ/761/22    : NIR galaxy cluster candidates in the SPT survey (Song+, 2012)
 J/MNRAS/423/3561 : X-ray clusters from XMM (Clerc+, 2012)
 J/ApJS/199/34   : Clusters of galaxies in SDSS-III (Wen+, 2012)
 J/ApJ/746/178   : The augmented maxBCG cluster catalog (Rykoff+, 2012)
 J/A+A/535/A65   : Galaxy clusters in the 4 CFHTLS Wide fields (Durret+, 2011)
 J/A+A/534/A109  : MCXC Meta-Catalogue X-ray galaxy Clusters (Piffaretti+, 2011)
 J/ApJ/736/21    : Gal. clusters optical cat. from AMF on SDSS DR6 (Szabo+ 2011)
 J/ApJS/191/254  : GMBCG galaxy cluster catalog from SDSS DR7 (Hao+, 2010)
 J/AJ/137/2981   : Northern Optical Cluster Survey. III. (Gal+, 2009)
 J/ApJ/684/905   : Galaxy clusters from IRAC Shallow Survey (Eisenhardt+, 2008)
 J/MNRAS/379/867 : BCG C4 cluster catalog (von der Linden+, 2007)
 J/ApJ/660/239   : MaxBCG catalog of galaxy clusters from SDSS (Koester+, 2007)
 J/AJ/128/1017   : Northern Optical Cluster Survey. IV. (Lopes+, 2005)
 J/ApJS/148/243  : Catalog of clusters of galaxies from SDSS (Bahcall+, 2003)
 J/AJ/123/1807   : SDSS galaxy clusters redshifts (Goto+, 2002)
 http://www.sdss3.org/ : SDSS-III home page
 http://risa.stanford.edu/redmapper/ : redMaPPer home page

Byte-by-byte Description of file: table1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label   Explanations
--------------------------------------------------------------------------------
   1-  6 I6     ---     ID      [1/183668] redMaPPer Cluster identifier
   8- 27 A20    ---     Name    redMaPPer Cluster name (RMJHHMMSS.s+DDMMSS.s)
  29- 39 F11.7  deg     RAdeg   Right Ascension in decimal degrees (J2000)
  41- 51 F11.7  deg     DEdeg   Declination in decimal degrees (J2000)
  53- 58 F6.4   ---     zlambda [0.08/0.55] Cluster photo-z; z_{lambda}_
  60- 65 F6.4   ---   e_zlambda [0.004/0.03] Gaussian error estimate in zlambda
  67- 72 F6.2   ---     lambda  [20/300.2] Richness estimate
  74- 78 F5.2   ---   e_lambda  [1/74] Gaussian error estimate in lambda
  80- 84 F5.3   ---     S       [1/3.4] Richness scale factor; see Equ. 22
  86- 93 F8.5   ---     zspec   [0.05/0.9]?=-1 SDSS spectroscopic redshift (1)
  95-113 I19    ---     ObjID   SDSS DR8 CAS object identifier
 115-120 F6.3   mag     imag    [13.2/20.2] Dereddened i-band cmodel magnitude
                                 (for most likely central galaxy)
 122-126 F5.3   mag   e_imag    [0.002/0.4] Error in imag
 128-133 F6.3   mag     umagm   [13.4/29.8] Dereddened u-band model magnitude
                                 (for most likely central galaxy)
 135-140 F6.3   mag   e_umagm   [0.006/13] Error in umagm
 142-147 F6.3   mag     gmagm   [14.6/29] Dereddened g-band model magnitude
                                 (for most likely central galaxy)
 149-154 F6.3   mag   e_gmagm   [0.002/20] Error in gmagm
 156-161 F6.3   mag     rmagm   [13.6/22.5] Dereddened r-band model magnitude
                                 (for most likely central galaxy)
 163-167 F5.3   mag   e_rmagm   [0.002/6] Error in rmagm
 169-174 F6.3   mag     imagm   [13.2/21.7] Dereddened i-band model magnitude
                                 (for most likely central galaxy)
 176-180 F5.3   mag   e_imagm   [0.002/0.3] Error in imagm
 182-187 F6.3   mag     zmagm   [12.9/27.6] Dereddened z-band model magnitude
                                 (for most likely central galaxy)
 189-194 F6.3   mag   e_zmagm   [0.003/34] Error in zmagm
 196-202 F7.3   Lsun    iLum    [8.8/169.3] Total membership-weighted i-band
                                 luminosity
 204-212 E9.3   ---     PCen0   [0.2/1] Centering probability (2)
 214-224 F11.7  deg     RA0deg  Right Ascension in decimal degrees (J2000) (2)
 226-236 F11.7  deg     DE0deg  Declination in decimal degrees (J2000) (2)
 238-256 I19    ---     ID0     DR8 CAS object identifier (2)
 258-266 E9.3   ---     PCen1   [0/0.5] Centering probability (3)
 268-278 F11.7  deg     RA1deg  Right Ascension in decimal degrees (J2000) (3)
 280-290 F11.7  deg     DE1deg  Declination in decimal degrees (J2000) (3)
 292-310 I19    ---     ID1     DR8 CAS object identifier (3)
 312-320 E9.3   ---     PCen2   [0/0.4] Centering probability (4)
 322-332 F11.7  deg     RA2deg  Right Ascension in decimal degrees (J2000) (4)
 334-344 F11.7  deg     DE2deg  Declination in decimal degrees (J2000) (4)
 346-364 I19    ---     ID2     DR8 CAS object identifier (4)
 366-374 E9.3   ---     PCen3   [0/0.3] Centering probability (5)
 376-386 F11.7  deg     RA3deg  Right Ascension in decimal degrees (J2000) (5)
 388-398 F11.7  deg     DE3deg  Declination in decimal degrees (J2000) (5)
 400-418 I19    ---     ID3     DR8 CAS object identifier (5)
 420-428 E9.3   ---     PCen4   [0/0.2] Centering probability (6)
 430-440 F11.7  deg     RA4deg  Right Ascension in decimal degrees (J2000) (6)
 442-452 F11.7  deg     DE4deg  Declination in decimal degrees (J2000) (6)
 454-472 I19    ---     ID4     DR8 CAS object identifier (6)
 474-479 F6.4   ---     PZbin1  [0.03/0.42] PZ bin number 1
 481-486 F6.4   ---     PZbin2  [0.04/0.44] PZ bin number 2
 488-493 F6.4   ---     PZbin3  [0.04/0.45] PZ bin number 3
 495-500 F6.4   ---     PZbin4  [0.04/0.46] PZ bin number 4
 502-507 F6.4   ---     PZbin5  [0.05/0.47] PZ bin number 5
 509-514 F6.4   ---     PZbin6  [0.05/0.49] PZ bin number 6
 516-521 F6.4   ---     PZbin7  [0.06/0.50] PZ bin number 7
 523-528 F6.4   ---     PZbin8  [0.06/0.51] PZ bin number 8
 530-535 F6.4   ---     PZbin9  [0.07/0.53] PZ bin number 9
 537-542 F6.4   ---     PZbin10 [0.07/0.54] PZ bin number 10
 544-549 F6.4   ---     PZbin11 [0.08/0.55] PZ bin number 11
 551-556 F6.4   ---     PZbin12 [0.08/0.57] PZ bin number 12
 558-563 F6.4   ---     PZbin13 [0.08/0.59] PZ bin number 13
 565-570 F6.4   ---     PZbin14 [0.09/0.61] PZ bin number 14
 572-577 F6.4   ---     PZbin15 [0.09/0.63] PZ bin number 15
 579-584 F6.4   ---     PZbin16 [0.10/0.65] PZ bin number 16
 586-591 F6.4   ---     PZbin17 [0.10/0.67] PZ bin number 17
 593-598 F6.4   ---     PZbin18 [0.10/0.69] PZ bin number 18
 600-605 F6.4   ---     PZbin19 [0.11/0.71] PZ bin number 19
 607-612 F6.4   ---     PZbin20 [0.11/0.73] PZ bin number 20
 614-619 F6.4   ---     PZbin21 [0.12/0.75] PZ bin number 21
 621-629 E9.3   ---     PZ1     [0/0.14] The P(z) evaluated at PZbin1 (7)
 631-639 E9.3   ---     PZ2     [0/0.37] The P(z) evaluated at PZbin2 (7)
 641-649 E9.3   ---     PZ3     [0/0.82] The P(z) evaluated at PZbin3 (7)
 651-659 E9.3   ---     PZ4     [0/1.62] The P(z) evaluated at PZbin4 (7)
 661-669 E9.3   ---     PZ5     [0/3] The P(z) evaluated at PZbin5 (7)
 671-679 E9.3   ---     PZ6     [0.0001/4.7] The P(z) evaluated at PZbin6 (7)
 681-689 E9.3   ---     PZ7     [0.01/8.4] The P(z) evaluated at PZbin7 (7)
 691-699 E9.3   ---     PZ8     [0.4/24.4] The P(z) evaluated at PZbin8 (7)
 701-709 E9.3   ---     PZ9     [3.6/52.7] The P(z) evaluated at PZbin9 (7)
 711-719 E9.3   ---     PZ10    [10.6/82] The P(z) evaluated at PZbin10 (7)
 721-729 E9.3   ---     PZ11    [14/94.2] The P(z) evaluated at PZbin11 (7)
 731-739 E9.3   ---     PZ12    [11.2/82.8] The P(z) evaluated at PZbin12 (7)
 741-749 E9.3   ---     PZ13    [5.6/53.5] The P(z) evaluated at PZbin13 (7)
 751-759 E9.3   ---     PZ14    [1.3/25.8] The P(z) evaluated at PZbin14 (7)
 761-769 E9.3   ---     PZ15    [0.2/9] The P(z) evaluated at PZbin15 (7)
 771-779 E9.3   ---     PZ16    [0.03/2.2] The P(z) evaluated at PZbin16 (7)
 781-789 E9.3   ---     PZ17    [0.005/0.4] The P(z) evaluated at PZbin17 (7)
 791-799 E9.3   ---     PZ18    [0.0007/0.06] The P(z) evaluated at PZbin18 (7)
 801-809 E9.3   ---     PZ19    [0.0001/0.005] The P(z) evaluated at PZbin19 (7)
 811-819 E9.3   ---     PZ20    [0/0.00031] The P(z) evaluated at PZbin20 (7)
 821-829 E9.3   ---     PZ21    The P(z) evaluated at PZbin21 (7)
--------------------------------------------------------------------------------
Note (1): For most likely center.
Note (2): For most likely central.
Note (3): For second most likely central.
Note (4): For third most likely central.
Note (5): For fourth most likely central.
Note (6): For fifth most likely central.
Note (7): Full redshift probability distribution P(z);
           see Sheldon et al. (2012ApJS..201...32S)
--------------------------------------------------------------------------------

Byte-by-byte Description of file: table2.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label  Explanations
--------------------------------------------------------------------------------
   1-  6 I6     ---     ID     [1/183668] redMaPPer Cluster identifier
   8- 18 F11.7  deg     RAdeg  Right Ascension in decimal degrees (J2000)
  20- 30 F11.7  deg     DEdeg  Declination in decimal degrees (J2000)
  32- 36 F5.3   ---     R      [0/1.8] Distance from cluster center; h^-1^Mpc
  38- 42 F5.3   ---     PMem   [0/1] Membership probability
  44- 49 F6.3   mag     imag   [12.9/21] Dereddened i-band cmodel magnitude
  51- 55 F5.3   mag   e_imag   [0.002/7.6] Error in imag
  57- 62 F6.3   mag     umagm  [13.4/31.3] Dereddened u-band model magnitude
  64- 69 F6.3   mag   e_umagm  [0.006/38] Error in umagm
  71- 76 F6.3   mag     gmagm  [14/31.1] Dereddened g-band model magnitude
  78- 83 F6.3   mag   e_gmagm  [0.002/96] Error in gmagm
  85- 90 F6.3   mag     rmagm  [13/26] Dereddened r-band model magnitude
  92- 96 F5.3   mag   e_rmagm  [0.002/8.2] Error in rmagm
  98-103 F6.3   mag     imagm  [12.8/26] Dereddened i-band model magnitude
 105-110 F6.3   mag   e_imagm  [0.002/13] Error in imagm
 112-117 F6.3   mag     zmagm  [12.5/29] Dereddened z-band model magnitude
 119-124 F6.3   mag   e_zmagm  [0.003/99] Error in zmagm
 126-133 F8.5   ---     zspec  [0.02/0.92]?=-1 SDSS spectroscopic redshift
 135-153 I19    ---     ObjID  SDSS DR8 CAS object identifier
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

References:
    Rozo & Rykoff   Paper II.       2014ApJ...783...80R
    Rozo et al.     Paper III.      2015MNRAS.450..592R
    Rozo et al.     Paper IV.       2015MNRAS.453...38R
    Rykoff et al.   from DES data   2016ApJS..224....1R   Cat. J/ApJS/224/1

================================================================================
(End)                 Prepared by [AAS], Emmanuelle Perret [CDS]    16-Nov-2016
