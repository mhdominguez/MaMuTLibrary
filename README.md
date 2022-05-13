# MaMuTLibrary
MaMuT Dataset Manipulation Script Library

These Perl scripts are used to filter, color, annotate, filter, subset, merge, and export track data from MaMuT datasets in order to visualize or quantify cell behaviors.



# Script listing

## MaMuT break divisions
### Break all divisions in a MaMuT dataset, not preserving even the closest daughter links
 usage: `perl MaMuT_dataset_break_divisions.pl mamut_dataset.xml`<br>
 
<br><br>
## MaMuT color spots
### Colors spots and links in a MaMuT dataset, according to instruction provided...
usages: `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml col=-3197141` #(fixed Java color all spots/links)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml col=255,255,128` #(fixed Java color all spots/links)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml rnd` #(random color by track)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml vec_xz=4` #(color track directionality converting angles to hue in single plane xy yz xz, by moving window size in timeframes)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml tis` #(color by tissueID i.e. populated by SVF)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml div` #(color division nodes)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml den=4` #(color by density, by number of radii units around each cell to count other cells)<br>
                `perl MaMuT_dataset_color_spots.pl mamut_dataset.xml vel=4` #(color by velocity, by moving window size in timeframes)<br>
### note this script is single-threaded and slow!<br>
 
<br><br>
## MaMuT dataset SVF tissue type annotations to vanilla TGMM/MaMuT dataset
### Add tissue type annotations to a vanilla TGMM->MaMuT dataset, provided an equivalent TGMM->SVF->MaMuT dataset that is already annotated:
 usage: `perl MaMuT_dataset_dataset_add_nearestneighbor_tissue_types_from_SVF.pl mamut_dataset_TGMMplain.xml mamut_dataset_SVF.xml`<br>

<br><br>
## MaMuT dataset radius annotations, export table of spots (CSV format: Annotation,Timepoint,Position_T,ID) with changed radii comparing one dataset to another
### Export cell radius modifications between two MaMuT datasets, to a table:<br>
 usage: `perl MaMuT_dataset_dataset_export_radius_annotations.pl mamut_dataset_original.xml mamut_dataset_annotated.xml`<br>
	where mamut_dataset_annotated.xml is the modified MaMuT dataset file, with cells modified by increased radius ('E'), or reduced radius('Q')<br>
	
<br><br>
## MaMuT subtract background movement and print track coordinates in time
### Reconstructs each [track] lineage coordinates from a MaMuT dataset into a .tsv table, after subtracting/correcting for nearby motion provided in the second MaMuT .xml file
###   it prints each [track] lineage in a row of the output table<br>
###   columns of the output table represent each timepoint<br>
###   each cell's XYZ coordinates are printed in the appropriate cell of the table for its lineage/track and timepoint<br>
 usage: `perl MaMuT_dataset_print_track_coordinates_in_time.pl dataset_mamut.xml dataset_background_mamut.xml`<br>
### note this script is single-threaded and slow!<br>
  
<br><br>
## MaMuT dataset TGMM division annotation xml generator
### Input MaMuT .xml file and a matching MaMuT .xml file with changes to the radius of the cells as follows to generate classifierAnnotations_YYYYMMDDTHHMMSS.xml files
###   (you still need MATLAB to generate final classifier matrix files for TGMM use)<br>
 usage: `perl MaMuT_dataset_dataset_TGMM_div_annotate.pl mamut_dataset_original.xml mamut_dataset_annotated.xml Path/to/GMEMtracking3D/XML_finalResult_lht` <br>
   where mamut_dataset_annotated.xml is the saved modified MaMuT dataset file, where division nodes (division parents) have been modified with increased radius ('E') for correct divisions, and reduced radius('Q') for incorrect divisions<br>
   
<br><br>
## MaMuT downsample density
### Downsamples a MaMuT dataset by factor X, going toward uniform cell density starting backwards and working forwards, with tracks allocated randomly among the bins
 usage: `perl MaMuT_dataset_downsample_density_backwards_in_time.pl dataset_mamut.xml X`<br>
	where X is the desired fold reduction<br>
### note this script is single-threaded and slow!
 
<br><br>
## MaMuT dataset radius annotations, export list of spots annotated
### List all spots in a dataset according to radius annotation format: Annotation,Timepoint,Position_T,ID,Name
 usage: `perl MaMuT_dataset_export_allspots.pl mamut_dataset.xml`<br>

<br><br>
## MaMuT print track coordinates in time
### Reconstructs a MaMuT dataset by [track] lineage coordinates into a .tsv table
###   it prints each [track] lineage in a row of the output table
###   columns of the output table represent each timepoint
###   each cell's XYZ coordinates are printed in the appropriate cell of the table for its lineage/track and timepoint
 usage: `perl MaMuT_dataset_print_track_coordinates_in_time.pl dataset_mamut.xml`<br>
### note this script is single-threaded and slow!

<br><br>
## MaMuT split dataset by manual color (or tissue type)
### Splits a MaMuT dataset into separate datasets, with tracks allocated by MANUAL_COLOR or TISSUE_TYPE annotations
 usage: `perl MaMuT_dataset_split_manual_color.pl dataset_mamut.xml`<br>
### this is the quick-n-dirty version of the script and does not exhaustively ensure that entire tracks stay together -- therefore it is best used for SVF2MM output where TISSUE_TYPE is forward propagated through each entire track
### note this script uses 8 threads in parallel, and requires package libparallel-forkmanager-perl on Ubuntu

<br><br>
## MaMuT split dataset by manual color (or tissue type)
### Splits a MaMuT dataset into separate datasets, with tracks allocated by MANUAL_COLOR or TISSUE_TYPE annotations
 usage: `perl MaMuT_dataset_split_manual_color.pl dataset_mamut.xml`<br>
### note this script is single-threaded and slow!

<br><br>
## MaMuT split dataset by tracks, randomly into X bins
### Splits a MaMuT dataset into separate datasets, with tracks allocated randomly among the bins
 usage: `perl MaMuT_dataset_split_random_subgroups.pl dataset_mamut.xml X`<br>
	where X is the desired number of bins<br>
### note this script is single-threaded and slow!

<br><br>
## MaMuT split dataset by timepoints
### Splits a MaMuT dataset into two datasets, one containing spots and tracks bounded by start/stop, the other with spots that are out-of-bounds
 usage: `perl MMaMuT_dataset_split_spot_timepoint.pl dataset_mamut.xml start=0 stop=40`<br>
	where 0 and 40 create the timepoint bounds for the resulting dataset<br>
### note this script is relatively fast, since it just splits the spots and doesn't really reconstruct lineages -- user may need to open and re-save the dataset in order for everything to work properly after splitting by timepoint

<br><br>
## MaMuT split dataset by track filter(s)
### Filter a MaMuT dataset by splitting into 2 datasets (.0.xml and .1.xml), with .0.xml containing rejected tracks/cells, and .1.xml containing filtered tracks/cells according to criteria specified...
 usages: `perl MaMuT_dataset_split_track_filter.pl dataset_mamut.xml track_start_min=20 track_stop_max=150 track_duration_min=50` #(filters tracks starting after timepoint 20, ending before timepoint 150, and with minimum duration 50 timeframes<br>
         `perl MaMuT_dataset_split_track_filter.pl dataset_mamut.xml track_displacement_min=5 track_displacement_max=20` #(filters tracks with net displacement between 5 and 20 units)<br>
         `perl MaMuT_dataset_split_track_filter.pl dataset_mamut.xml avg_vel_min=5 avg_vel_max=20` #(filters tracks with average velocity between between 5 and 20 units)<br>
         `perl MaMuT_dataset_split_track_filter.pl dataset_mamut.xml number_splits_max=3 rad_min=4` #(filters tracks with fewer than 3 splits and cell radii greater than 4)<br>
         `perl MaMuT_dataset_split_track_filter.pl dataset_mamut.xml pos_x_min=-40 pos_y_max=200 pos_z_min=-1700 pos_z_max=0` #(filters tracks with cells bounded by max/min coordinates as specified)<br>
         `perl MaMuT_dataset_split_track_filter.pl dataset_mamut.xml ann=spotRadiusAnnotations_20211124T093902_downOnly.txt` #(filters tracks with spots on a radius-annotation CSV list i.e. generated by MaMuT_dataset_dataset_export_radius_annotations.pl)<br>
### note this script is single-threaded and slow!

<br><br>
## MaMuT split dataset by XYZ coordinates into bins, backwards in time
### Splits MaMuT dataset into n bins by dividing it into mean/stdev or percentile bins along a particular axis, using chronologic last cell in the track to assign the bin
 usage: `perl MaMuT_dataset_split_XYZ_coordinate_backward_in_time.pl mamut_dataset.xml bin=2 med_z` #med_z can be replaced with mean_x, mean_y, mean_z, med_x, or med_y<br>

<br><br>
## MaMuT split dataset by XYZ coordinates into bins
### Splits MaMuT dataset into n bins by dividing it into mean/stdev or percentile bins along a particular axis, using chronologic first cell in the track to assign the bin
 usage: `perl MaMuT_dataset_split_XYZ_coordinate.pl mamut_dataset.xml bin=2 med_z` #med_z can be replaced with mean_x, mean_y, mean_z, med_x, or med_y<br>

<br><br>
## MaMuT dataset subtract spots from input list (using Name)
### Subtracts spots from a MaMuT dataset, with spots inputted from file(s) in radius annotation CSV format: Annotation,Timepoint,Position_T,ID,Name
 usage: `perl MaMuT_dataset_subtract_spots_quickndirty.pl dataset_mamut.xml spots1.txt spots2.txt spots3.txt ...`<br>
### this is quick-n-dirty, and does not exhaustively ensure that entire tracks stay together -- therefore it is best used for SVF2MM output where TISSUE_TYPE is forward propagated through each entire track

<br><br>
## MaMuT dataset subtract spots from input list (using ID)
### Subtracts spots from a MaMuT dataset, with spots inputted from file(s) in radius annotation CSV format: Annotation,Timepoint,Position_T,ID
 usage: `perl MaMuT_dataset_subtract_spots_quickndirty.pl dataset_mamut.xml spots1.txt spots2.txt spots3.txt ...`<br>
### this is quick-n-dirty, and does not exhaustively ensure that entire tracks stay together -- therefore it is best used for SVF2MM output where TISSUE_TYPE is forward propagated through each entire track

<br><br>
## MaMuT print radius annotation up/down counts
### Prints up and down counts for radius annotation files (CSV format: Annotation,Timepoint,Position_T,ID), as would be generated by MaMuT_dataset_dataset_export_radius_annotations.pl
 usage: `bash MaMuT_radius_annotations_print_UpDown_counts.sh`<br>

<br><br>
## MaMuT return java / RGB colors
### Interconverts Java colors (used by MaMuT) and RGB
 usages: `perl MaMuT_return_java_color_for_RGB.pl -235908` #(converts this java color to RGB)<br>
         `perl MaMuT_return_java_color_for_RGB.pl 128,255,72` #(converts this RGB color to Java)<br>
         `perl MaMuT_return_java_color_for_RGB.pl 128 255 72` #(converts this RGB color to Java)<br>

<br><br>
## MaMuT track coordinates daughter separation analyze
### Analyzes tracks for separation of daughters (best used with manually annotated tracks; TGMM is not sufficiently accurate with cell divisions to be relied on for rigorous quantification)
###   prints a table of daughter-mother distances, and daughter-daughter distances
 usage: `perl MaMuT_track_coordinates_daughter_separation_analyze.pl dataset_mamut.track_coordinates_in_time.tsv start=20 stop=120`<br>

<br><br>
## MaMuT track coordinates pairwise analyze
### Analyzes track pairs exported from MaMuT_dataset_print_track_coordinates_in_time.pl, printing:
###   a table of features such as Begin Coord, Timepoint start, Begin Diff, End Diff
###   a summary table of single-axis track crossing events, and relative positions of track pair begin or end along that axis
 usage: `perl MaMuT_track_coordinates_pairwise_analyze.pl dataset_mamut.track_coordinates_in_time.tsv start=15 stop=90 admit_stop=45 cell_max_dist=250 xpos` #(analyzes track pairs along x-axis, starting at timepoint 15, stopping at timepoint 90, and without admission of new tracks after timepoint 45, using maximum distance of 250 between the cells as a filter<br>
### this script is very slow!

<br><br>
## MaMuT track coordinates single data export
### Summarizes all tracks exported from MaMuT_dataset_print_track_coordinates_in_time.pl, printing a table of features such as begin/end XYZT, peak displacement, average density, etc
 usage: `perl MaMuT_track_coordinates_single_data_export.pl dataset_mamut.track_coordinates_in_time.tsv`<br>
