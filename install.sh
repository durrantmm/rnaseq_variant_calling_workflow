#!/usr/bin/env bash

IS_CONDA="$(which conda)"

if [ ${#IS_CONDA} -eq "0" ]
    then
        echo "It seems that you do not have conda installed. Download it from https://www.continuum.io/downloads"
        exit 1
fi

echo "Removing preexisting rnavcw environment if exists"
source deactivate
conda remove --name rnavcw --all --yes

echo "Creating the new rnavcw environment from environment.yml"
conda env create -f environment.yml

if [ $? -ne 0 ]
    then
        echo "Problem creating the conda environment. I don't know what to tell you."
        exit 1
fi

echo "Downloading the test data."
wget https://s3-us-west-1.amazonaws.com/mdurrant/biodb/bundles/rnaseq_variant_calling_workflow/test.tar.gz
tar -zxvf test.tar.gz
rm test.tar.c0gz

if [ $? -ne 0 ]
    then
        echo "Could not download the test data"
        exit 1
fi



echo "------------------------------------------------------------------------------------------------------"
echo "INSTALLATION COMPLETE"
echo "------------------------------------------------------------------------------------------------------"
echo "NOTE: You must specify the required information in the config.yaml file to run the workflow properly."
echo "NOTE: Once configured, enter:"
echo ""
echo "> source activate rnavcw"
echo "> snakemake"
echo ""
echo "In this directory to run the workflow."
echo "------------------------------------------------------------------------------------------------------"