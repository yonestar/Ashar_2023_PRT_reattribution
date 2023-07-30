# Semantic analyses of pain reattribution in pain reprocessing therapy (PRT)

This study aims to understand how people with chronic back pain think about the underlying causes of their pain, to test whether pain reprocessing therapy (PRT) changes those attributions, and whether reattribution may be a mechanism by which PRT leads to pain relief.

## Team Lead
Yoni Ashar

## Associated publication

revise and resubmit at the moment...

## Structure of this github repository 

* `scripts` contains MATLAB and R scripts for analyses
* `data` contains de-identified survey data
* `figures` contains figures, published and unpublished, from this project

The main analyses of mind-brain attribution scores are in `analyze_categories.m`

The analyses of changes in specific words are in `analyses_of_words.m`

The text-scaling analyses are in the `parrot` subfolder. See also the parrot toolbox: XXX

## Dependencies

* CanlabCore (for various plotting and modelling functions): https://github.com/canlab/CanlabCore
* Parrot (for the analyses using that toolbox): https://github.com/wilryh/parrot
* brewermap for plot colors: https://github.com/DrosteEffect/BrewerMap
* Mediation Toolbox for mediation analyses: https://github.com/canlab/MediationToolbox
