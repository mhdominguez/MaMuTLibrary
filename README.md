# MaMuTLibrary
MaMuT Dataset Manipulation Script Library

These Perl scripts are used tofilter, color, annotate, filter, subset, merge, and export track data from MaMuT datasets in order to visualize or quantify cell behaviors.

# Script listing

## MaMuT break divisions
 Break all divisions in a MaMuT dataset, not preserving even the closest daughter links
 usage: perl MaMuT_dataset_break_divisions.pl mamut_dataset.xml
 
## MaMuT color spots
Colors spots and links in a MaMuT dataset, according to instruction provided...
usages: perl MaMuT_dataset_color_spots.pl mamut_dataset.xml col=-3197141 #(fixed Java color all spots/links)
                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml col=255,255,128 #(fixed Java color all spots/links)
                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml rnd #(random color by track)
                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml vec_xz=4 #(color track directionality converting angles to hue in single plane xy yz xz, by moving window size in timeframes)
                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml tis #(color by tissueID i.e. populated by SVF)
                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml div #(color division nodes)
                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml den=4 #(color by density, by number of radii units around each cell to count other cells)
                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml vel=4 #(color by velocity, by moving window size in timeframes)
 note this script is single-threaded and slow!
 
## MaMuT dataset SVF tissue type annotations to vanilla TGMM/MaMuT dataset
 Add tissue type annotations to a vanilla TGMM->MaMuT dataset, provided an equivalent TGMM->SVF->MaMuT dataset that is already annotated:
 usage: perl MaMuT_dataset_dataset_add_nearestneighbor_tissue_types_from_SVF.pl mamut_dataset_TGMMplain.xml mamut_dataset_SVF.xml
 
## MaMuT dataset radius annotations, export table of spots (CSV format: Annotation,Timepoint,Position_T,ID) with changed radii comparing one dataset to another
 Export cell radius modifications between two MaMuT datasets, to a table:
 usage: perl MaMuT_dataset_dataset_export_radius_annotations.pl mamut_dataset_original.xml mamut_dataset_annotated.xml
	where mamut_dataset_annotated.xml is the modified MaMuT dataset file, with cells modified by increased radius ('E'), or reduced radius('Q')
  
  
