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
 
 * After fitting a model, you can predict erosion rates at any location (note, the fitted model exponents are currently hard-coded into this function and passing ```output.mdl``` doesn't actually do anything):
 ````matlab
 eval_model_monthly(model_data,44,output.mdl)
 ````
 
 * Many of the other functions are outdated or are used by the functions I mentioned here. Be careful deleting functions (sorry).
