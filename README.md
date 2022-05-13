# MaMuTLibrary
MaMuT Dataset Manipulation Script Library

These Perl scripts are used to filter, color, annotate, filter, subset, merge, and export track data from MaMuT datasets in order to visualize or quantify cell behaviors.

# Script listing

## MaMuT break divisions
### Break all divisions in a MaMuT dataset, not preserving even the closest daughter links<br>
 usage: `perl MaMuT_dataset_break_divisions.pl mamut_dataset.xml`<br>
 
## MaMuT color spots
### Colors spots and links in a MaMuT dataset, according to instruction provided...<br>
usages: `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml col=-3197141` #(fixed Java color all spots/links)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml col=255,255,128` #(fixed Java color all spots/links)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml rnd` #(random color by track)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml vec_xz=4` #(color track directionality converting angles to hue in single plane xy yz xz, by moving window size in timeframes)<br>`
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml tis` #(color by tissueID i.e. populated by SVF)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml div` #(color division nodes)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml den=4` #(color by density, by number of radii units around each cell to count other cells)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml vel=4` #(color by velocity, by moving window size in timeframes)<br>
### note this script is single-threaded and slow!<br>
 
## MaMuT dataset SVF tissue type annotations to vanilla TGMM/MaMuT dataset
### Add tissue type annotations to a vanilla TGMM->MaMuT dataset, provided an equivalent TGMM->SVF->MaMuT dataset that is already annotated:<br>
 usage: `perl MaMuT_dataset_dataset_add_nearestneighbor_tissue_types_from_SVF.pl mamut_dataset_TGMMplain.xml mamut_dataset_SVF.xml`<br>
 
## MaMuT dataset radius annotations, export table of spots (CSV format: Annotation,Timepoint,Position_T,ID) with changed radii comparing one dataset to another
### Export cell radius modifications between two MaMuT datasets, to a table:<br>
 usage: `perl MaMuT_dataset_dataset_export_radius_annotations.pl mamut_dataset_original.xml mamut_dataset_annotated.xml`<br>
	where mamut_dataset_annotated.xml is the modified MaMuT dataset file, with cells modified by increased radius ('E'), or reduced radius('Q')<br>
	
## MaMuT subtract background movement and print track coordinates in time
### Reconstructs each [track] lineage coordinates from a MaMuT dataset into a .tsv table, after subtracting/correcting for nearby motion provided in the second MaMuT .xml file<br>
###  it prints each [track] lineage in a row of the output table<br>
###  columns of the output table represent each timepoint<br>
###  each cell's XYZ coordinates are printed in the appropriate cell of the table for its lineage/track and timepoint<br>
 usage: `perl MaMuT_dataset_print_track_coordinates_in_time.pl dataset_mamut.xml dataset_background_mamut.xml`<br>
### note this script is single-threaded and slow!<br>
  
## MaMuT dataset TGMM division annotation xml generator
### Input MaMuT .xml file and a matching MaMuT .xml file with changes to the radius of the cells as follows to generate classifierAnnotations_YYYYMMDDTHHMMSS.xml files<br>
###  (you still need MATLAB to generate final classifier matrix files for TGMM use)<br>
 usage: `perl MaMuT_dataset_dataset_TGMM_div_annotate.pl mamut_dataset_original.xml mamut_dataset_annotated.xml Path/to/GMEMtracking3D/XML_finalResult_lht` <br>
   where mamut_dataset_annotated.xml is the saved modified MaMuT dataset file, where division nodes (division parents) have been modified with increased radius ('E') for correct divisions, and reduced radius('Q') for incorrect divisions<br>
   
## MaMuT downsample density
### Downsamples a MaMuT dataset by factor X, going toward uniform cell density starting backwards and working forwards, with tracks allocated randomly among the bins<br>
 usage: `perl MaMuT_dataset_downsample_density_backwards_in_time.pl dataset_mamut.xml X`<br>
	where X is the desired fold reduction<br>
### note this script is single-threaded and slow!<br>
   
   
   
 
