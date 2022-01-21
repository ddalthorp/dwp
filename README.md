# dwp
Estimating the number of birds or bats killed at a wind power facility involves
adjusting observed carcass counts to account for imperfect detection. To account
for carcasses that lie outside the searched area, the count is adjusted by 
dividing by the estimated proportion of carcasses lying within the searched area 
or the *density-weighted proportion* (dwp). The dwp package provides tools for 
estimating dwp for use in mortality estimators such as GenEst (Dalthorp et al. 
2018) and Evidence of Absense (Dalthorp et al. 2017).


## DISCLAIMER

This software has been approved for release by the U.S. Geological Survey (USGS). 
Although the software has been subjected to rigorous review, the USGS reserves 
the right to update the software as needed pursuant to further analysis and 
review. No warranty, expressed or implied, is made by the USGS or the U.S. 
Government as to the functionality of the software and related material nor shall
the fact of release constitute any such warranty. Furthermore, the software is
released on condition that neither the USGS nor the U.S. Government shall be held
liable for any damages resulting from its authorized or unauthorized use.

## Installation
Setup and installation require several steps. Do not skip any steps.

### Updated version of R (>= 4.1.0, released May 2021):
R is free and open source software for statistical computing. If R is not
installed on your computer or if your version of R is <4.1.0, download and 
install the latest version from https://cran.r-project.org/, following the 
instructions provided at the site. In particular, "Download" and then "install 
R for the first time" (if working in Windows), or "Download" and then follow 
the further instructions on the subsequent web page (if working on Mac OS or 
Linux-like OS). If you already have an older copy of R installed on your 
computer, the new version will be installed alongside the old. Unless you know 
a reason why you want to keep both versions, it is usually a good idea to 
uninstall the old version to avoid confusion and clutter. 

NOTE TO EXPERIENCED R USERS: When you install a new version of R, packages that 
you previously installed under an older version may not be immediately 
available to the new R. If not, the easiest way to make them available is to 
copy the package folders in your old "library" folder into the "library" folder 
in your new R installation. Then, enter `update.packages()` in R. If asked 
about a CRAN mirror, choose the nearest location. If you are working in Windows 
OS and are asked whether you want to install packages "from source", choose 
"No".


### Third-party packages: 
Several third-party pacakges are required; all are free and open source and 
available from CRAN. The easiest way to install them is to run the following 
commands in R (with guidance concerning potential dialog boxes given below the 
commands):


```

package_req <- c("boot", "expint", "GenEst", "gpclib", "gtools", "invgamma",
"magrittr", "MASS", "matrixStats", "methods", "mvtnorm", "numDeriv", "plotrix",
"pracma", "sf", "statmod", "VGAM")

package_new <- package_req[!(package_req %in% installed.packages()[,"Package"])] 
if(length(package_new) > 0) install.packages(package_new)

```
-- If asked about a "CRAN mirror", choose the nearest location.

-- If asked whether you want to use a "personal library", choose "Yes"

-- If you are on Windows and are asked whether you want to install packages and 
their dependencies "from source", choose "No" (unless you are ready to go to 
lunch, in which case, you can select "Yes" and the installation may well be
done by the time you get back).

### dwp: 
Click on "Tags" under the "Repository" tab on the left sidebar at 
https://code.usgs.gov/ecosystems/dwp and then click the link for the specific 
release you want. 

-- For Windows, download the compressed folder dwp_1.x.x.zip (do not unzip) and
note where it is stored. You will install from the local .zip folder. 

-- For Mac OS or Unix-like OS, download the compressed file dwp_1.x.x.tar.gz
and note where it is stored. You will install from the local .tar.gz file. 

If you are working directly in R (not R Studio), run the following command:
```
install.packages(file.choose()) # and navigate to the package archive file you just downloaded: dwp_1.x.x.xxx
```

If you are working in R Studio:

Click "Install" in the Packages pane.

Select "Package Archive File (.zip; .tar.gz)" as "Install from:" in the dialog 
box.

Browse to where you saved the zip file, and open it so it appears in the
Package archive" space.

Click the Install button on the dialog box.

## Getting Started

Download the User Guide from <??? we have it but it has to go through the copy editor. With GenEst, we had a copy of the pre-final (i.e., version sent to the copy editor for formatting and style) but approved copy of the user guide posted at code.usgs.gov and then later linked to the T&M docs after they made it through Tacoma.>

In R, type the following commands for an overview of the main commands in the 
package:

```
library(dwp)
?dwp

```

Help files for dwp functions are accessible in the standard R way, for example:

```
?ddFit
```

## Further Reading
dwp User Guide: <in press...can we put a temp copy up at code.usgs.gov?>

Dalthorp, D, Simonis, J, Madsen, L, Huso, M, Rabie, P, Mintz, J, Wolpert, R, 
Studyvin, J, and Korner-Nievergelt, F. 2018. Generalized Mortality Estimator.
R package version 1.4.6.
https://CRAN.R-project.org/package=GenEst

Dalthorp, D, Huso, M, Dail, D, and Kenyon, J. 2017. Evidence of Absence (v2.0)
Software and user Guide. https://pubs.er.usgs.gov/publication/ds1055