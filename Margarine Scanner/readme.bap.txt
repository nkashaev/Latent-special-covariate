Dries F. Benoit, Stefan Van Aelst, and Dirk Van den Poel, "Outlier 
Robust Bayesian Multinomial Choice Modeling", Journal of Applied 
Econometrics, Vol. 31, No. 7, 2016, pp. 1445-1466.

Both real and simulated data sets were used in the paper. The artificial 
data were simulated as described in the paper.

The real dataset used in the paper (margarine scanner dataset) was 
obtained from the R-package "bayesm" that comes with the following book:

  Rossi, P., Allenby, G., McCulloch, R., 2005. Bayesian Statistics and 
  Marketing, John Wiley & Sons, New York.

And can be found at the following url (19/06/2015):

  Rossi, P., 2015. bayesm: Bayesian Inference for 
  Marketing/Micro-econometrics. R package version 3.0-1. 
  http://cran.r-project.org/web/packages/bayesm

The complete dataset is comprised of 9196 purchases of ten brands of 
margarine by 517 households in Springfield, MO. The panel structure of
the data is transformed into a cross-sectional structure by retaining
only the first purchase of the household. This results in a dataset
with 242 households and an equal number of observed choices. The five
most purchased brands were retained.

The resulting dataset is named "margarine_subset.csv". It is an ASCII
file in DOS format and is zipped in the file bap-data.zip. Unix/Linux
users should use "unzip -a".

The dataset contains 8 variables:

"id"       : the unique identifier of the household
"choice"   : the brand that was purchased (with 1=Blue Bonnet, 
             2=Fleischman's, 3=house brand, 4=generic brand, 5=Shed Spread)
"pBB"      : price of Blue Bonnet
"pFM"      : price of Fleischmann's
"pHse"     : price of house brand
"pGen"     : price of generic brand
"pSS"      : price of Shed Spread
"income"   : log of household income

Prof. Dr. Dries Benoit
Assistant Professor of Data Analytics
Faculty of Economics & Business Administration
Tweekerkenstraat 2
9000 Gent
Belgium
