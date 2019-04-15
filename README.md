moralizing-gods-reanalysis
============

[![License: CC BY-NC-SA 4.0](https://licensebuttons.net/l/by-nc-sa/4.0/80x15.png)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

The following code was based on the analyses in: 

> Whitehouse, H., François, P., Savage, P. E., Currie, T. E., Feeney, K. C., Cioni, E., Purcell, R., Ross, R. M., Larson, J., Baines, J., ter Haar, B., Covey, A., Turchin, P. (2019). Complex societies precede moralizing gods throughout world history. Nature.

with data file "exportdat.csv" scraped from the Seshat database on 10 Jan 2018. Steps 1 to 8 are the original analysis code with some minor reorganization, while steps 9 to 12 are extensions and modifications.

The code was developed on R v3.5.3 with the following packages for steps 1 to 11:

- maps (3.3.0)
- plotrix (3.7-4)
- dplyr (0.8.0.1)
- plry (1.8.4)
- testthat (2.0.1)
- viridis (0.5.1)
- rstan (2.18.2)
- ggplot (3.1.0)
- rethinking (1.88)

Step 12 was developed with the following packages:

- glmmADMB (0.8.5)
- glmmTMB (0.2.3)
- nlme (3.1-137)
- lme4 (1.1-21)
- stargazer (5.2.2)
- DHARMa (0.2.4)
- bbmle (1.0.20)
- knitr (1.22)
- car (3.0-2)
- bootpredictlme4 (0.1)
- reshape (0.8.8)
- yarrr (0.1.5)

Step 13 was developed with the following packages:

- dplyr (0.8.0.1)
- ggplot2 (3.1.0.9000)

To run the primary analyses:

1. Set the working directory to "moralizing-gods-reanalysis"
2. In R, type `source("run_project.r")`.

This code is derivative modification of [moralizing-gods](https://github.com/pesavage/moralizing-gods) by Patrick Savage and Peter Turchin, used under Creative Commons License [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/), modified and extended by Bret Beheim, Martin Lang and Rachel Spicer under the same license. See LICENSE.md for details.
