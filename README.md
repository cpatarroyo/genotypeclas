# Genotype classificator

This is a machine learning-based classifier for the classification of <i>Phytophthora infestans</i> microsatellite genotypes into their corresponding clonal lineages. This classifier uses [AdaBag](https://doi.org/10.18637/JSS.V054.I02) as its machine learning algorithm and it uses the 12 standard microsatellite markers described by [Li et al. 2013](https://doi.org/10.1016/j.mimet.2012.11.021) for the characterization of <i>P. infestans</i> isolates. This model was trained using a 1392 microsatellite genotype dataset from isolates described previously from Colombia [Chaves et al. 2020](https://doi.org/10.1094/PHYTO-05-19-0175-R), Peru [Lindqvist-Kreuze et al. 2020](https://doi.org/10.1111/ppa.13125), United States [Wang et al. 2017](https://doi.org/10.1111/mec.14000), India and Europe [Dey et al. 2018](https://doi.org/10.1038/s41598-018-22192-1). 
This model was trained with a dataset of 1392 microsatellite genotypes in GenAlEx format (Training_DB2.csv), and can be used to classify unknown genotypes as it is, or it can be re-trained with additional genotype data. Instructions below. 

## Usage

This genotype classifier is divided into two functional parts:
1) Model training
2) Genotype classification

The model training is only done when there are new classified genotypes to add to the training dataset or when you want to change the training dataset. Detailed instructions on how to do both below.

### Instalation

To use this classifier you need to have R (and [RStudio](https://posit.co/download/rstudio-desktop/) for easier local use) installed. Once you have it you will need to install several libraries using the following commands:
```
install.packages(shiny)
install.packages(poppr)
install.packages(adabag)
install.packages(caret)
```

After the installation you can clone this GitHub repository to your local machine. This repository can also be cloned to a server for remote access/data processing. 

### Model training

For training the model the [ML_train_script.R] can be executed. This will create/overwrite the [classifier/Trained_model.rds] object with the trained AdaBag model. By default this script uses the [Training_DB.csv] file for training the model. There are two simple options to change the dataset used for model training:
1) Add the new dataset to the local folder where the repository was cloned and rename it to "Training_DB.csv".
2) in the [ML_train_script.R] line 8 change the file route "Training_DB.csv" to the route corresponding to your training set.

### Genotype classification

For local usage, if you are using RStudio open the [classifier/app.R] file, and click the "Run App" button. This will open the webapp UI where you can upload a GenAlEx file (or any poppr supported format) with the <i>P. infestans</i> microsatellite genotypes to be classified. Remember that the used markers must be the within the markers described by [Li et al. 2013](https://doi.org/10.1016/j.mimet.2012.11.021). 

For remote use from an own server you have to access from your web browser the ip address of the server and the port desginated for the Shiny Server service, together with the route to the folder containing the present repository. 