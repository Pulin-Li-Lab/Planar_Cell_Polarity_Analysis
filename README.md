Analysis code from Wallach et el. “Synthetic reconstitution of planar polarity initiation reveals collective migration as a symmetry-breaking cue”


Visualization and analysis of live recordings of cell polarity

Analysis work flow is as follows:

Track cells and identify polarization using FIJI Trackmate and Quantify polarity 

Use QuantifyPolarityTrackmateCombiner to intergrate the two outputs

QPTMImportandProcessing sets up the previous output into formats for following scripts as well as offering general visualization

BulkMarkovModels uses preceding inputs to train and implement Markov models on binned data

SingleCellMarkovModels implements the previously trained models on single cell traces.

…

RosePlotGen was used to identify migrating cells and visualize their polarization angle via rose plot
