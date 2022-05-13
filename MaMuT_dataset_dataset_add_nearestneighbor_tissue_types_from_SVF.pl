#!/usr/bin/perl

# MaMuT dataset SVF tissue type annotations to vanilla TGMM/MaMuT dataset
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Add tissue type annotations to a vanilla TGMM->MaMuT dataset, provided an equivalent TGMM->SVF->MaMuT dataset that is already annotated:
# usage: perl MaMuT_dataset_dataset_add_nearestneighbor_tissue_types_from_SVF.pl mamut_dataset_TGMMplain.xml mamut_dataset_SVF.xml

use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;


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
	#my $GMEMfinalResult_folder = $_[2];
	my @lines_main; #unchanging lines during reconstruction
	my $line_number_spots = -1; #indicates where to put <AllSpots> blocks when reconstructing the file at the end
	my @lines_spots_pre;
	my @lines_spots_pre_ids;
	my @lines_spots_pre_names;
	my @lines_spots_ann;
	my @lines_spots_ann_ids;
	my @lines_spots_ann_names;
	my @frames_spots_pre;
	my @frames_spots_ann;
	my $line_number_tracks = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my $line_number_filters = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my $going;
	#my $this_spot_id;
	
	my @cloud_x_TGMM;
	my @cloud_x_SVF;
	my @cloud_y_TGMM;
	my @cloud_y_SVF;
	my @cloud_z_TGMM;
	my @cloud_z_SVF;
	my @cloud_tissue_SVF;
	my @cloud_manual_color_SVF;
	my @cloud_t_SVF;
	my @cloud_t_SVF_stop;
	my @cloud_t_SVF_start;
	
	my @TGMM_xml;
	
	unless ( defined($dataset_file_pre) && defined($dataset_file_ann) && $dataset_file_pre =~ /\w/ && $dataset_file_ann =~ /\w/ )  {
		print "This script needs three arguments, in format: mamut_dataset_original.xml mamut_dataset_annotated.xml\n";
		return 0;
	}

	my $average_velocity = 0, $number_of_edges = 0;
	my $image_data = "";
	#---------------
	#parse dataset.mamut files
	#---------------	
	my $this_frame;
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
						$this_frame = $2;
						push( @frames_spots_ann, $this_frame );
						$cloud_t_SVF_start[$this_frame] = scalar(@cloud_x_SVF);
						my @this_block = ( $_ );
						my @this_block_ids = ( -1 );
						my @this_block_names = ( "" );
						while( <FILE> ) {
							chomp;
							if ( $_ =~ /<Spot/ && $_ =~ /ID=\"(\d+)\"/ ) {
								push( @this_block_ids, $1 );
								if ( $_ =~ /NAME=\"([\w\ ]+)?\"/i ) {
									push( @this_block_names, $1 );
									#print $1 . "\n";
								} else {
									push( @this_block_names, "" );
								}
								#capture coordinates
								my @xyz;
								if ($_ =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ ) {
									$xyz[0] = $1; #sprintf("%.3f", $1);
									#push( @cloud_x_SVF, $1 );
								} else {
								
								}
								if ($_ =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ ) {
									$xyz[1] = $1; #sprintf("%.3f", $1);
									#push( @cloud_y_SVF, $1 );
								}
								if ($_ =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
									$xyz[2] = $1; #sprintf("%.3f", $1);
									#push( @cloud_z_SVF, $1 );
								}
								if ($_ =~ /POSITION_T=\"([+-]?\d+(\.?\d*))\"/ ) {
									$xyz[3] = $1; #sprintf("%.3f", $1);
									#push( @cloud_z_SVF, $1 );
								}
								if ( defined($xyz[0]) && defined($xyz[1]) && defined($xyz[2]) && $xyz[0] =~ /\d/ && $xyz[1] =~ /\d/ && $xyz[2] =~ /\d/ && defined($xyz[3]) && $xyz[3] =~ /\d/ ) {
									push( @cloud_x_SVF, $xyz[0] );
									push( @cloud_y_SVF, $xyz[1] );
									push( @cloud_z_SVF, $xyz[2] );
									push( @cloud_t_SVF, $xyz[3] );
									if ( $_ =~ /MANUAL_COLOR=\"([+-]?\d+(\.?\d*))\"/ ) {
										push( @cloud_manual_color_SVF, $1 );
									} else {
										push( @cloud_manual_color_SVF, "" );
									}
									if ( $_ =~ /TISSUE_TYPE=\"([+-]?\d+(\.?\d*))\"/ ) {
										push( @cloud_tissue_SVF, $1 );
									} else {
										push( @cloud_tissue_SVF, "" );
										
										
									}
									#print "here2\n";
								} else {
									#print join(',', @xyz ) . "\n";
								}
								#if ($_ =~ /RADIUS=\"([+-]?\d+(\.?\d*))\"/ ) {
								#	#$xyz[3] = sprintf("%.3f", $1);
								#} else {
								#	#$xyz[3] = "N\/A";
								#}
								

							
							} else {
								push( @this_block_ids, -1 );
								push( @this_block_names, "" );
							}
							push( @this_block, $_ );
							if ( $_ =~ /(^|^\s+)<\/SpotsInFrame/ ) {
								$cloud_t_SVF_stop[$this_frame] = scalar(@cloud_x_SVF);
								last;
							}
						}
						push( @lines_spots_ann, \@this_block );
						push( @lines_spots_ann_ids, \@this_block_ids );
						push( @lines_spots_ann_names, \@this_block_names );
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
			} elsif ( $_ =~ /<ImageData/ ) {
				$image_data = $_;
			}
		}
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error opening $path/$dataset_file_ann!\n";
		return;
	}
	
	my $max_t = max( @cloud_t_SVF );
	my $min_t = min( @cloud_t_SVF );
	

	if ( open(FILE, "<$path/$dataset_file_pre" ) ) {
		print "Reading spots from $dataset_file_pre.\n";
		flock(FILE, LOCK_EX);
		while( <FILE> ) {
			chomp;
			push( @TGMM_xml, $_ );
			if ( $_ =~ /(^|^\s+)<AllSpots/ ) {
				push( @lines_main, $_ );
				$line_number_spots = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					push( @TGMM_xml, $_ );
					if ( $_ =~ /(^|^\s+)<SpotsInFrame\s+frame=\"(\d+)\"/i ) {
						$this_frame = $2;
						push( @frames_spots_pre, $this_frame );
						my @this_block = ( $_ );
						my @this_block_ids = ( -1 );
						my @this_block_names = ( "" );
						while( <FILE> ) {
							chomp;
							push( @TGMM_xml, $_ );
							if ( $_ =~ /<Spot/ && $_ =~ /ID=\"(\d+)\"/ ) {
								push( @this_block_ids, $1 );
								if ( $_ =~ /NAME=\"([\w\ ]+)?\"/i ) {
									push( @this_block_names, $1 );
									#print $1 . "\n";
								} else {
									push( @this_block_names, "" );
								}
								
								#capture coordinates
								if ( $this_frame <= $max_t && $this_frame >= $min_t ) { #don't compare coordinates unless we are within the same timepoint range
								#my @xyz;
								if ($_ =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ ) {
									#$xyz[0] = sprintf("%.3f", $1);
									push( @cloud_x_TGMM, $1 );
								}					
								if ($_ =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ ) {
									#$xyz[1] = sprintf("%.3f", $1);
									push( @cloud_y_TGMM, $1 );
								}
								if ($_ =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
									#$xyz[2] = sprintf("%.3f", $1);
									push( @cloud_z_TGMM, $1 );
								}
								#if ($_ =~ /RADIUS=\"([+-]?\d+(\.?\d*))\"/ ) {
								#	#$xyz[3] = sprintf("%.3f", $1);
								#} else {
								#	#$xyz[3] = "N\/A";
								#}
								}
								
							} else {
								push( @this_block_ids, -1 );
								push( @this_block_names, "" );
							}
							push( @this_block, $_ );
							last if ( $_ =~ /(^|^\s+)<\/SpotsInFrame/ );
						}
						push( @lines_spots_pre, \@this_block );
						push( @lines_spots_pre_ids, \@this_block_ids );
						push( @lines_spots_pre_names, \@this_block_names );
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
					push( @TGMM_xml, $_ );
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
					push( @TGMM_xml, $_ );
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
	my @X_SVF = set_stats_plus( @cloud_x_SVF );
	my @Y_SVF = set_stats_plus( @cloud_y_SVF );
	my @Z_SVF = set_stats_plus( @cloud_z_SVF );
	#print "hello";
	#now, figure out transformations from SVF to TGMM
	my @X_TGMM = set_stats_plus( @cloud_x_TGMM );
	my @Y_TGMM = set_stats_plus( @cloud_y_TGMM );
	my @Z_TGMM = set_stats_plus( @cloud_z_TGMM );
	#print "hello".scalar(@TGMM_xml) . "\n";

	
	#print "Stats X:\n" . join ( ",", @X_TGMM ) . "\n" . join ( ",", @X_SVF ) . "\n\n";
	#print "Stats Y:\n" . join ( ",", @Y_TGMM ) . "\n" . join ( ",", @Y_SVF ) . "\n\n";
	#print "Stats Z:\n" . join ( ",", @Z_TGMM ) . "\n" . join ( ",", @Z_SVF ) . "\n\n";
	
	#use median value for each coordinate to determine midpoint i.e. translation, and use min/max to determine scale factor

	my $scale_X = ( $X_SVF[4] - $X_SVF[3] ) / ( $X_TGMM[4] - $X_TGMM[3] );
	my $scale_Y = ( $Y_SVF[4] - $Y_SVF[3] ) / ( $Y_TGMM[4] - $Y_TGMM[3] );
	my $scale_Z = ( $Z_SVF[4] - $Z_SVF[3] ) / ( $Z_TGMM[4] - $Z_TGMM[3] );
	
	my $scale__X = $X_SVF[2] / $X_TGMM[2];
	my $scale__Y = $Y_SVF[2] / $Y_TGMM[2];
	my $scale__Z = $Z_SVF[2] / $Z_TGMM[2];
	
	$scale_X = ($scale_X+ $scale__X)/2;
	$scale_Y = ($scale_Y+ $scale__Y)/2;
	$scale_Z = ($scale_Z+ $scale__Z)/2;
	my $scale_R = ( $scale_Y + $scale_X + $scale_Z ) / 3;
	my $std_X = $X_SVF[2] / 4;
	my $std_Y = $Y_SVF[2] / 4;
	my $std_Z = $Z_SVF[2] / 4;
	
	my $translate_X = $X_SVF[1] - ( $scale_X * $X_TGMM[1] );
	my $translate_Y = $Y_SVF[1] - ( $scale_Y * $Y_TGMM[1] );
	my $translate_Z = $Z_SVF[1] - ( $scale_Z * $Z_TGMM[1] );	
	#go through original file and scale coordinates, mapping colors and tissues
	my $start_here = 1;
	for (my $i=0; $i<scalar( @TGMM_xml); $i++ ) {
	#map {
		$_ = $TGMM_xml[$i];
							if ( $_ =~ /<Spot/ ) { # && $_ =~ /ID=\"(\d+)\"/ ) {
							
								#capture coordinates
								my $this_t;
								if ($_ =~ /POSITION_T=\"([+-]?\d+(\.?\d*))\"/ ) {
									$this_t = $1;
									next unless ( $1 <= $this_t && $1 >= $this_t ); #don't compare coordinates unless we are within the same timepoint range
								}

								my @xyz;
								if ($_ =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ ) {
									my $new_X = ( $scale_X * $1 ) + $translate_X;
									$_ =~ s/POSITION_X=\"([+-]?\d+(\.?\d*))\"/POSITION_X=\"$new_X\"/;
									$xyz[0] = $new_X;
								}					
								if ($_ =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ ) {
									my $new_Y = ( $scale_Y * $1 ) + $translate_Y;
									$_ =~ s/POSITION_Y=\"([+-]?\d+(\.?\d*))\"/POSITION_Y=\"$new_Y\"/;
									$xyz[1] = $new_Y;
								}
								if ($_ =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
									my $new_Z = ( $scale_Z * $1 ) + $translate_Z;
									$_ =~ s/POSITION_Z=\"([+-]?\d+(\.?\d*))\"/POSITION_Z=\"$new_Z\"/;
									$xyz[2] = $new_Z;
								}
								if ($_ =~ /RADIUS=\"([+-]?\d+(\.?\d*))\"/ ) {
									my $new_R = $scale_R * $1;
									$_ =~ s/RADIUS=\"([+-]?\d+(\.?\d*))\"/RADIUS=\"$new_R\"/;
								}
								
								
								#now find nearest neighbor spot and acquire its color/tissue
								if ( defined($xyz[0]) && defined($xyz[1]) && defined($xyz[2]) && $xyz[0] =~ /\d/ && $xyz[1] =~ /\d/ && $xyz[2] =~ /\d/ ) {
									my $closest_spot_distance = vec3_distance(($cloud_x_SVF[0]-$xyz[0]),($cloud_y_SVF[0]-$xyz[1]),($cloud_z_SVF[0]-$xyz[2]));
									my $closest_spot = 0;
									my $this_distance;
									#my $start_here = 0;
									
									for ( my $k=$cloud_t_SVF_start[$this_t]; $k<$cloud_t_SVF_stop[$this_t]; $k++ ) {
										next unless ($cloud_t_SVF[$k] == $this_t);
										#$start_here = $k if ( $start_here != $k );
										#only consider spots within a fourth standard deviation in all three axes
										next if ( abs($cloud_x_SVF[$k]-$xyz[0]) > $std_X );
										next if ( abs($cloud_y_SVF[$k]-$xyz[1]) > $std_Y );
										next if ( abs($cloud_z_SVF[$k]-$xyz[2]) > $std_Z );
									
										$this_distance = vec3_distance(($cloud_x_SVF[$k]-$xyz[0]),($cloud_y_SVF[$k]-$xyz[1]),($cloud_z_SVF[$k]-$xyz[2]));
										if ( defined($this_distance) && $this_distance > 0 && $this_distance < $closest_spot_distance ) {
											$closest_spot_distance = $this_distance;
											$closest_spot = $k;
										}
									}
									
									#now, we have the closest spot
									$_ =~ s/MANUAL_COLOR=\"([+-]?\d+(\.?\d*))\"//g; #remove prior colors
									$_ =~ s/TISSUE_TYPE=\"([+-]?\d+(\.?\d*))\"//g; #remove prior colors
								
									$_ =~ s/\s\/>/ MANUAL_COLOR=\"$cloud_manual_color_SVF[$closest_spot]\" \/>/ if (defined($cloud_manual_color_SVF[$closest_spot]) && $cloud_manual_color_SVF[$closest_spot] =~ /\d/ );
									$_ =~ s/\s\/>/ TISSUE_TYPE=\"$cloud_tissue_SVF[$closest_spot]\" \/>/ if (defined($cloud_tissue_SVF[$closest_spot]) && $cloud_tissue_SVF[$closest_spot] =~ /\d/ );
								
								}
								
			
								
							}	 elsif ( $_ =~ /<ImageData/ ) {
								$_ = $image_data;
							}

		$TGMM_xml[$i] = $_;
		
		#if ( $i % 10000 == 0 ) {
		if ( rand() < 0.0001 ) {
			print "Progress: " . int(100 * $i / scalar(@TGMM_xml) ) . " % \n";
		
		}
	#} @TGMM_xml;
	}
	
	#output file
	#if ( $out_xml ne "" ) {
		#mkdir "annForCellDivDiscrWithTempWin";
		if ( open( FILE, ">$path/$dataset_file_pre\." . join( '', datetoday() ) . "T" . join( '', timenow() ) . "SVF_scaled_colored.xml"  ) ) {
			flock( FILE, LOCK_EX );
			#print FILE qq~<?xml version="1.0" encoding="UTF-8"?>\n<document>\n~;
			print FILE join( "\n", @TGMM_xml );
			#print FILE qq~<\/document>\n~;
			flock( FILE, LOCK_UN );
			close( FILE );
		} else {
			print "Could not open XML file for writing!\n";
		}
	#} else {
	#	print "Nothing to output, most likely this means there was an error within MaMuT dataset files!\n";
	#}
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

sub set_stats_plus {

	my $avg = 0;
	my $std = 0;
	my $med = 0;
	#my $mdn;
	#use Data::Dumper;	print Dumper( \@_ );
	my $min = @_[0];
	my $max = @_[0];
	for ( my $i = 0; $i < scalar(@_); $i++ ) {
		#find the outliers of theta, and the average theta -> the average/outlier differential will determine how excessively to adjust the exponent
		$avg += @_[$i];#abs(@_[$i]);
		$max = @_[$i] if ( @_[$i] > $max );
		$min = @_[$i] if ( @_[$i] < $min );
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
	
	return ( $avg, $med, $std, $min, $max );
}


#my $MANUAL_COLOR="-65281";
#print join( ',', hsl_to_rgb((240-(0.5 * 240))/360,1,.5)) . "\n";
main(@ARGV);



