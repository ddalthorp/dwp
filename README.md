# dwp
Estimating the number of birds or bats killed at a wind power facility involves
adjusting observed carcass counts to account for imperfect detection. To account
for carcasses that lie outside the searched area, the count is adjusted by 
dividing by the estimated proportion of carcasses lying within the searched area 
or the *density-weighted proportion* (dwp). The dwp package provides tools for 
estimating dwp for use in mortality estimators such as GenEst (Dalthorp et al. 
2018) and Evidence of Absense (Dalthorp et al. 2017).


## DISCLAIMER
12345678901234567890123456789012345678901234567890123456789012345678901234567890
Version 1.0 of the software has been approved for release by the U.S. Geological Survey (USGS). 
Although the software has been subjected to rigorous review, the USGS reserves 
the right to update the software as needed pursuant to further analysis and 
review. No warranty, expressed or implied, is made by the USGS or the U.S. 
Government as to the functionality of the software and related material nor shall
the fact of release constitute any such warranty. Furthermore, the software is
released on condition that neither the USGS nor the U.S. Government shall be held
liable for any damages resulting from its authorized or unauthorized use.

Minor modifications and bug fixes in subsequent versions have not reviewed
by USGS.

## Revisions:
v1.1 
Method for estimating the area searched in each ring has been updated. The new method does not
rely on `gpclib::area.poly`, which is no longer available on CRAN. It estimates the area searched in 
each ring by averaging the arc lengths of the intersections of the inner and outer ring boundaris
with the search polygon.

## Further Reading
dwp User Guide: Daniel Dalthorp, Manuela Huso, Mark Dalthorp, Jeff Mintz
Accounting for the Fraction of Carcasses outside the Searched Area and the Estimation of Bird and
Bat Fatalities at Wind Energy Facilities
https://arxiv.org/abs/2201.10064

Dalthorp, D, Simonis, J, Madsen, L, Huso, M, Rabie, P, Mintz, J, Wolpert, R, 
Studyvin, J, and Korner-Nievergelt, F. 2018. R package. GenEst: Generalized Mortality Estimator.
https://CRAN.R-project.org/package=GenEst

Dalthorp, D, Huso, M, Dail, D, and Kenyon, J. 2017. Evidence of Absence (v2.0)
Software and user Guide. https://pubs.er.usgs.gov/publication/ds1055