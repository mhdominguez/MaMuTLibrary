#!/usr/bin/perl

# MaMuT dataset TGMM division annotation xml generator
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Input MaMuT .xml file and a matching MaMuT .xml file with changes to the radius of the cells as follows to generate classifierAnnotations_YYYYMMDDTHHMMSS.xml files
#  (you still need MATLAB to generate final classifier matrix files for TGMM use)
# usage: perl MaMuT_dataset_dataset_TGMM_div_annotate.pl mamut_dataset_original.xml mamut_dataset_annotated.xml Path/to/GMEMtracking3D/XML_finalResult_lht 
#   where mamut_dataset_annotated.xml is the saved modified MaMuT dataset file, where division nodes (division parents) have been modified with increased radius ('E') for correct divisions, and reduced radius('Q') for incorrect divisions

use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;
#use constant {
#	LOCK_SH => 1,
#	LOCK_NB => 4,
#	LOCK_UN => 8,
#	LOCK_EX => 2,
#};

sub datetoday {

	my ( $year, $month, $day );
	
	my @chains = localtime(time);
	$chains[5] += 1900;

	$year = $chains[5];
	$month = sprintf("%02d", $chains[4]+1 );
	$day = sprintf("%02d", $chains[3] );
	undef @chains;	
	return ( $year, $month, $day );
}
sub timenow {

	my ( $sc, $mn, $hr, undef ) = localtime(time);

	$sc = sprintf("%02d", $sc );
	$mn = sprintf("%02d", $mn );
	$hr = sprintf("%02d", $hr );
	return ( $hr, $mn, $sc );

}

use constant PI    => 4 * atan2(1, 1);
sub vec2_distance_and_angle { #returns the distance between two spots, and the angle of the vector (0deg is horizontal right to left)
	my ( $x_displacement, $y_displacement ) = @_;
	return ( sqrt(($x_displacement**2) + ($y_displacement**2)), atan2($y_displacement,$x_displacement)*180/PI );
}

sub vec3_distance { #returns the distance between two spots, and the angle of the vector (0deg is horizontal right to left)
	my ( $x_displacement, $y_displacement, $z_displacement ) = @_;
	return ( sqrt(($x_displacement**2) + ($y_displacement**2) + ($z_displacement**2)) );
}


	
#==================
#Main subroutine
#==================

sub main {
	#my (  @arguments ) = @_;
	my $dataset_file_pre = $_[0];
	my $dataset_file_ann = $_[1];
	my $GMEMfinalResult_folder = $_[2];
	my @lines_main; #unchanging lines during reconstruction
	my $line_number_spots = -1; #indicates where to put <AllSpots> blocks when reconstructing the file at the end
	my @lines_spots_pre;
	my @lines_spots_pre_ids;
	my @lines_spots_ann;
	my @lines_spots_ann_ids;	
	my @frames_spots_pre;
	my @frames_spots_ann;
	my $line_number_tracks = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my $line_number_filters = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my $going;
	#my $this_spot_id;
	
	unless ( defined($dataset_file_pre) && defined($dataset_file_pre) && defined($GMEMfinalResult_folder) && $dataset_file_pre =~ /\w/ && $dataset_file_ann =~ /\w/ && $GMEMfinalResult_folder =~ /\w/ ) {
		print "This script needs three arguments, in format: mamut_dataset_original.xml mamut_dataset_annotated.xml Path/to/GMEMtracking3D/XML_finalResult_lht \n";
		return 0;
	}

	my $average_velocity = 0, $number_of_edges = 0;
	#---------------
	#parse dataset.mamut files
	#---------------	
	if ( open(FILE, "<$path/$dataset_file_pre" ) ) {
		print "Reading spots from $dataset_file_pre.\n";
		flock(FILE, LOCK_EX);
		while( <FILE> ) {
			chomp;
			if ( $_ =~ /(^|^\s+)<AllSpots/ ) {
				push( @lines_main, $_ );
				$line_number_spots = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<SpotsInFrame\s+frame=\"(\d+)\"/i ) {
						push( @frames_spots_pre, $2 );
						my @this_block = ( $_ );
						my @this_block_ids = ( -1 );
						while( <FILE> ) {
							chomp;
							if ( $_ =~ /<Spot/ && $_ =~ /ID=\"(\d+)\"/ ) {
								push( @this_block_ids, $1 );
							} else {
								push( @this_block_ids, -1 );
							}
							push( @this_block, $_ );
							last if ( $_ =~ /(^|^\s+)<\/SpotsInFrame/ );
						}
						push( @lines_spots_pre, \@this_block );
						push( @lines_spots_pre_ids, \@this_block_ids );
					} elsif ( $_ =~ /(^|^\s+)<\/AllSpots/ ) {
						push( @lines_main, $_ );
						last;
					} else {
						#not inside a SpotsInFrame segment, and yet not finished with AllSpots, so do nothing here put keep reading next line
					}
				}
			} elsif ( $_ =~ /(^|^\s+)<AllTracks/ && $_ !~ /\/>\s*$/ ) {
				push( @lines_main, $_ );
				$line_number_tracks = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<\/AllTracks/ ) {
						push( @lines_main, $_ );
						last;
					} else {
						#not inside a Track segment, and yet not finished with AllTracks, so do nothing here put keep reading next line
					}
				}
			} elsif ( $_ =~ /(^|^\s+)<FilteredTracks/ && $_ !~ /\/>\s*$/ ) {
				push( @lines_main, $_ );
				$line_number_filters = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<\/FilteredTracks/ ) {
						push( @lines_main, $_ );	
						last;
					}
				}
			}
		}
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error opening $path/$dataset_file_pre!\n";
		return;
	}
	
	if ( open(FILE, "<$path/$dataset_file_ann" ) ) {
		print "Reading spots from $dataset_file_ann.\n";
		flock(FILE, LOCK_EX);
		while( <FILE> ) {
			chomp;
			if ( $_ =~ /(^|^\s+)<AllSpots/ ) {
				push( @lines_main, $_ );
				$line_number_spots = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<SpotsInFrame\s+frame=\"(\d+)\"/i ) {
						push( @frames_spots_ann, $2 );
						my @this_block = ( $_ );
						my @this_block_ids = ( -1 );
						while( <FILE> ) {
							chomp;
							if ( $_ =~ /<Spot/ && $_ =~ /ID=\"(\d+)\"/ ) {
								push( @this_block_ids, $1 );
							} else {
								push( @this_block_ids, -1 );
							}
							push( @this_block, $_ );
							last if ( $_ =~ /(^|^\s+)<\/SpotsInFrame/ );
						}
						push( @lines_spots_ann, \@this_block );
						push( @lines_spots_ann_ids, \@this_block_ids );
					} elsif ( $_ =~ /(^|^\s+)<\/AllSpots/ ) {
						push( @lines_main, $_ );
						last;
					} else {
						#not inside a SpotsInFrame segment, and yet not finished with AllSpots, so do nothing here put keep reading next line
					}
				}
			} elsif ( $_ =~ /(^|^\s+)<AllTracks/ && $_ !~ /\/>\s*$/ ) {
				push( @lines_main, $_ );
				$line_number_tracks = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<\/AllTracks/ ) {
						push( @lines_main, $_ );
						last;
					} else {
						#not inside a Track segment, and yet not finished with AllTracks, so do nothing here put keep reading next line
					}
				}
			} elsif ( $_ =~ /(^|^\s+)<FilteredTracks/ && $_ !~ /\/>\s*$/ ) {
				push( @lines_main, $_ );
				$line_number_filters = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<\/FilteredTracks/ ) {
						push( @lines_main, $_ );	
						last;
					}
				}
			}
		}
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error opening $path/$dataset_file_ann!\n";
		return;
	}
	
	my $out_xml == "";
	my $this_m;
	my $this_intensity;
	my $this_sv_idx;
	my $this_name;
	my $tgmm_spot_id;
	#print "Hello, here are your timepoints: " . join(',', @frames_spots_pre ) . "\n";
	#---------------
	#At each timepoint, match spots and determine radii
	#---------------
	for (my $t=0; $t<scalar(@frames_spots_ann); $t++ ) {
		
		#match the timepoints between MaMuT dataset files
		my $pre_frame_pointer = -1;
		for (my $u=0; $u<scalar(@frames_spots_pre); $u++ ) {
			if ( $frames_spots_pre[$u] eq $frames_spots_ann[$t] ) {
				$pre_frame_pointer = $u;
				last;
			}
		}
		#print "A: \$t = $t\n";
		if ( $pre_frame_pointer < 0 ) {
			print "Problem matching frame $frames_spots_ann[$t] within frames provided in original MaMuT dataset file, proceeding without any annotations to this frame!\n";
			next;
		}
		
		#now, open GMEMtracking file for this frame
		$frame_four_digits = sprintf("%04d", $frames_spots_ann[$t]);
		my %frame_spot_hash = ();
		if ( open(FILE, "<$GMEMfinalResult_folder/GMEMfinalResult_frame" . $frame_four_digits . ".xml" ) ) {
			flock(FILE, LOCK_EX);
			chomp;
			$_ = $/;
			$/ = undef;
			my $xml = <FILE>; 
			flock(FILE, LOCK_UN);
			close(FILE);
			$/ = $_;
			$xml =~ s/\n//g; #these files have newlines in between XML blocks
			
			while( $xml =~ /<GaussianMixtureModel id="(\d+)"\s+(.*?)>/ig ) {  #.*?m=\"(.*?)\"
				$frame_spot_hash{$1} = $2;
				#print "  GMEM spot added\n";
			}
			#print "Read " . scalar(keys %frame_spot_hash) . " spots in GMEMfinalResult_frame" . $frame_four_digits . ".xml.\n";
		} else {
			print "Error opening $path/$dataset_file_ann!\n";
			next;
		}
		#print "B: \$t = $t\n";
		#now, match all spots
		for ( my $k=0; $k<scalar(@{$lines_spots_ann_ids[$t]}); $k++ ) {
			next if ( $lines_spots_ann_ids[$t][$k] < 0);
			$going = FALSE;
			for ( my $j=0; $j<scalar(@{$lines_spots_pre_ids[$pre_frame_pointer]}); $j++ ) {
				if ( $lines_spots_pre_ids[$pre_frame_pointer][$j] eq $lines_spots_ann_ids[$t][$k] ) { 
					#matched the spot
					$going = TRUE;
					
					#now, compare radii
					my $radius_pre;
					if ( $lines_spots_pre[$pre_frame_pointer][$j] =~ /RADIUS=\"(\d+\.\d+)\"/ ) {
						$radius_pre = $1;
					}
					
					my $radius_ann;
					if ( $lines_spots_ann[$t][$k] =~ /RADIUS=\"(\d+\.\d+)\"/ ) {
						$radius_ann = $1;
					}
					
					my $class_output = "";
					if ( defined($radius_pre) && $radius_pre ne "" && defined($radius_ann) && $radius_ann ne "" ) {
						next if ( $radius_pre eq $radius_ann ); #they match radii, so not annotated
					
						if ( $radius_pre > $radius_ann + 0.0001 ) {
							$class_output = "cellDivisionWrong";
						} elsif ( $radius_ann > $radius_pre + 0.0001 ) {
							$class_output = "cellDivisionCorrect";
						}
					
					} else {
						print "Could not find radius for spot $lines_spots_pre_ids[$pre_frame_pointer][$j] $lines_spots_ann_ids[$t][$k], so not going to include.\n";
						print "  $lines_spots_ann[$t][$k]\n  $lines_spots_pre[$pre_frame_pointer][$j]\n";
						next;
					}
					
					#don't bother looking it up if there is no annotation
					next if ( $class_output eq "" );
					$this_name = undef;
					#capture value of ID for looking up m (centroid matrix) from GMEMfinalResult_frame XML file
					if ( $lines_spots_ann[$t][$k] =~ /name=(".*?")/ ) {
						$this_name = $1;
						$tgmm_spot_id = undef;
						
						if ( $this_name =~ /\((\d+)\)/ ) {
							$tgmm_spot_id = $1;
							$this_m = undef;
							$this_intensity = undef;
							$this_sv_idx = undef;							
							#print $frame_spot_hash{$tgmm_spot_id} . "\n";

							if ( $frame_spot_hash{$tgmm_spot_id} =~ /nu=(\".*?\s*\")/ ) {
								$this_intensity = $1;
							}
							if ( $frame_spot_hash{$tgmm_spot_id} =~ /m=\"(.*?)\s*\"/ ) { #note quotes not included in capture
								$this_m = $1;
							}							
							if ( $frame_spot_hash{$tgmm_spot_id} =~ /svIdx=(\".*?\s*\")/ ) {
								$this_sv_idx = $1;
							}
							if ( defined($this_m) && $this_m ne "" && defined($this_intensity) && $this_intensity ne "" && defined($this_sv_idx) && $this_sv_idx ) {
								
								#output what we need here -- most of this is just fluff
								$GMEMfinalResult_folder =~ s/\/$//g;
								$out_xml .= qq~<Surface name="Ellipsoid" id="~ . $tgmm_spot_id . qq~" numCoeffs="9" coeffs="0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 ~ . $this_m . qq~ " intensity=~ . $this_intensity . qq~ covarianceMatrixSize="3" svFilename="~ . "$GMEMfinalResult_folder/GMEMfinalResult_frame" . $frame_four_digits . ".svb" . qq~" svIdx=~ . $this_sv_idx . qq~ imFilename="t0~ . $frame_four_digits . qq~_s01.klb" class="~ . $class_output . qq~"></Surface>\n~;
							} else {
								print "Cannot find centroid matrix ($this_m), supervoxel index ($this_sv_idx), or intensity ($this_intensity) within hash table for MaMuT spot $lines_spots_pre_ids[$pre_frame_pointer][$j] / TGMM spot $tgmm_spot_id at frame $frames_spots_ann[$t], not using!\n";
							}
						} else {
							print "Cannot find name () of spot $lines_spots_pre_ids[$pre_frame_pointer][$j], not going to process this one.\n";
						}
					} else {
						print "Cannot find name of spot $lines_spots_pre_ids[$pre_frame_pointer][$j], not going to process this one.\n";
					}
					
					last;
				}
			}
			unless ( $going ) {
				print "Problem matching spot $lines_spots_ann_ids[$t][$k] in frame $frames_spots_ann[$t], proceeding without annotating this spot!\n";			
			}
		}
		#print "C: \$t = $t\n";
		#return;
	}
	
	#output file
	if ( $out_xml ne "" ) {
		mkdir "annForCellDivDiscrWithTempWin";
		if ( open( FILE, ">annForCellDivDiscrWithTempWin/classifierAnnotations_" . join( '', datetoday() ) . "T" . join( '', timenow() ) . ".xml"  ) ) {
			flock( FILE, LOCK_EX );
			print FILE qq~<?xml version="1.0" encoding="UTF-8"?>\n<document>\n~;
			print FILE $out_xml;
			print FILE qq~<\/document>\n~;
			flock( FILE, LOCK_UN );
			close( FILE );
		} else {
			print "Could not open annotation XML file for writing!\n";
		}
	} else {
		print "Nothing to output, most likely this means there were no annotations noted in MaMuT dataset files!\n";
	}
}

#==================
#Array remove duplicates
#==================
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}


#==================
#Comparison subroutines
#==================
sub max {
	my $max = shift;
	$max = $max > $_ ? $max : $_ for @_;
	return $max
}

sub min {
	my $min = shift;
	$min = $min < $_ ? $min : $_ for @_;
	return $min
}

sub set_stats {

	my $avg = 0;
	my $std = 0;
	my $med = 0;
	#my $mdn;
	#use Data::Dumper;	print Dumper( \@_ );
	for ( my $i = 0; $i < scalar(@_); $i++ ) {
		#find the outliers of theta, and the average theta -> the average/outlier differential will determine how excessively to adjust the exponent
		$avg += @_[$i];#abs(@_[$i]);
	}
	$avg /= scalar(@_);
	for ( my $i = 0; $i < scalar(@_); $i++ ) {
		#find the outliers of theta, and the average theta -> the average/outlier differential will determine how excessively to adjust the exponent
		$std += (@_[$i] - $avg) ** 2;
	}
	$std = sqrt( $std / (scalar(@_)-1) );
	
	#TODO: fix this median-determining function!
	@_ = sort { $a <=> $b } @_;#{ $a cmp $b };
	#@_ = sort(@_);#{ $a cmp $b };
	if ( scalar(@_) & 1 ) { #odd number of numbers
		$med = int( scalar(@_)/2);
		$med = @_[$med];
		#$med = "me!";
	} else {
		$med = scalar(@_) / 2;
		$med = ( @_[$med] + @_[$med -1] ) / 2;
	}
	
	return ( $avg, $med, $std );
}


#my $MANUAL_COLOR="-65281";
#print join( ',', hsl_to_rgb((240-(0.5 * 240))/360,1,.5)) . "\n";
main(@ARGV);



