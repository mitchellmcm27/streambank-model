# streambank-model
## Streambank erosion model in Matlab

### Run the model with supplied data
 * Clone the repo (in terminal)
 ````dos
 cd my/fav/directory
 git clone https://github.com/mitchellmcm27/streambank-model
 ````
 
 * Start up Matlab, navigate to the newly created ```streambank-model``` folder in your Matlab file tree
 * load "model_data.mat" into matlab, which contains the input data table:
 ````matlab
 load('model_data')
 ````
 
 * Run the model on a single site, animating the simulations:
 ````matlab
 train_model_monthly(model_data(44,:), 'animate')
 ````
 
 * Run the model on all sites without plotting or animating results:
 ````matlab
 output = train_model_monthly(model_data)
 ````
  * Animations are saved as .gifs in the ```gifs``` folder
 
 * The other optional argument is ```plt``` which will plot major simulation data, but no animations
 
 * The ```output``` struct contains 3 objects:
  1. ```output.tbl```: data table with rows for each modeled site and columns for input data, erosion rates, and predicted rates
  2. ```output.mdl1```: statistical model using K1 erodibility (see paper)
  3. ```output.mdl2```: statistical model using K2 erodibility (see paper)
 
 * After fitting a model, you can predict erosion rates at any location (note, the fitted model exponents are currently hard-coded into this function and passing ```output.mdl``` doesn't actually do anything):
 ````matlab
 eval_model_monthly(model_data,44,output.mdl)
 ````
