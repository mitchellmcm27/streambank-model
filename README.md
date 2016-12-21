# streambank-model
## Streambank erosion model in Matlab

Paper currently under review.

### Run the model with supplied data
 * Clone the repo (in terminal)
 ````dos
 cd my/fav/directory
 git clone https://github.com/mitchellmcm27/streambank-model.git
 ````
 
 * Start up Matlab, navigate to the newly created ```streambank-model``` folder in your Matlab file tree
 * load "model_data.mat" into matlab, which contains the input data table:
 ````matlab
 load('model_data')
 ````
 
 * Run the model on a single site:
 ````matlab
 train_model_monthly(model_data(44,:), 'animate')
 ````
  * 44 is a particular row (site), used for an example
  * ```animate``` (optional) tells the function to plot animations of simulations
  * ```plt``` (optional) plot major results but no animations
  * Animations are saved as .gifs in the ```gifs``` folder
  
 * Run the model on all sites:
 ````matlab
 output = train_model_monthly(model_data)
 ````
 
 * The ```output``` struct contains 3 objects:
  * ```output.tbl```: data table with rows for each modeled site and columns for input data, erosion rates, and predicted rates
  * ```output.mdl1```: statistical model using K1 erodibility (see paper)
  * ```output.mdl2```: statistical model using K2 erodibility (see paper)
 
 * After fitting a model, you can predict erosion rates at any location:
 ````matlab
 eval_model_monthly(model_data, 44, output.mdl)
 ````
  * The arguments are: input data table, row number, and fitted model
  * Note: The fitted model coefficients/exponents are currently hard-coded into the function and passing ```output.mdl``` doesn't actually do anything.
 
 * Many of the other functions in the repo are outdated or are used by the functions I mentioned here. Be careful deleting functions (sorry).
 
## Literature folder
 
The ```literature``` folder contains the manuscript based on this model (currently under review). The folder also contains a bibliographic database (.bib file), which contains references for all of the literature cited in the paper. You can browse the .bib file in [JabRef](https://www.fosshub.com/JabRef.html) or import it into [Mendeley](https://www.mendeley.com/) or [Zotero](https://www.zotero.org/). The .bib file is a [BiBTeX](http://www.bibtex.org/) database generated by Mendeley Desktop and cleaned up in JabRef.

## Other folders

The other folders contain supporting functions. Some are required for the model calculations, others are only for plotting, colorbars, etc.
