#!/usr/bin/perl

# MaMuT subtract background movement and print track coordinates in time
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Reconstructs each [track] lineage coordinates from a MaMuT dataset into a .tsv table, after subtracting/correcting for nearby motion provided in the second MaMuT .xml file
# it prints each [track] lineage in a row of the output table
# columns of the output table represent each timepoint
# each cell's XYZ coordinates are printed in the appropriate cell of the table for its lineage/track and timepoint
# usage: perl MaMuT_dataset_print_track_coordinates_in_time.pl dataset_mamut.xml dataset_background_mamut.xml

# note this script is single-threaded and slow!


use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;



sub vec3_distance { #returns the distance between two spots, and the angle of the vector (0deg is horizontal right to left)
	my ( $x_displacement, $y_displacement, $z_displacement ) = @_;
	return ( sqrt(($x_displacement**2) + ($y_displacement**2) + ($z_displacement**2)) );
}

#==================
#Main subroutine
#==================
sub main {
	#my (  @arguments ) = @_;
	my $dataset_file = $_[0];
	my $dataset_file_ann = $_[1];
	
	#my @lines_main; #unchanging lines during reconstruction
	my $line_number_spots = -1; #indicates where to put <AllSpots> blocks when reconstructing the file at the end
	my @lines_spots;
	my @frames_spots;
	
	my $line_number_tracks = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_tracks;
	#my @track_tracks;
	my @spots_in_tracks;
	my @spots_in_tracks_complete;
	my @track_in_tracks;
	my $line_number_filters = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_filters;
	my @track_in_filters;
	
	my @cloud_x_SVF;
	my @cloud_y_SVF;
	my @cloud_z_SVF;
	my @cloud_t_SVF;
	my @cloud_t_SVF_start;
	my @cloud_t_SVF_stop;
	my @cloud_id_SVF;
	my @cloud_mother_index_list_SVF;
	my %cloud_self_index_list;
	
	if ( open(FILE, "<$path/$dataset_file_ann" ) ) {
		print "Reading spots and track segments from $dataset_file_ann.\n";
		flock(FILE, LOCK_EX);
		while( <FILE> ) {
			chomp;
			if ( $_ =~ /(^|^\s+)<AllSpots/ ) {
				#push( @lines_main, $_ );
				#$line_number_spots = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<SpotsInFrame\s+frame=\"(\d+)\"/i ) {
						$this_frame = $2;
						push( @frames_spots_ann, $this_frame );
						$cloud_t_SVF_start[$this_frame] = scalar(@cloud_x_SVF);
						#my @this_block = ( $_ );
						#my @this_block_ids = ( -1 );
						#my @this_block_names = ( "" );
						while( <FILE> ) {
							chomp;
							if ( $_ =~ /<Spot/ && $_ =~ /ID=\"(\d+)\"/ ) {
								#push( @this_block_ids, $1 );
								my $this_id = $1;
								#if ( $_ =~ /NAME=\"([\w\ ]+)?\"/i ) {
								#	push( @this_block_names, $1 );
								#	#print $1 . "\n";
								#} else {
								#	push( @this_block_names, "" );
								#}
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
									$cloud_self_index_list{$this_id} = scalar(@cloud_x_SVF);
									push( @cloud_x_SVF, $xyz[0] );
									push( @cloud_y_SVF, $xyz[1] );
									push( @cloud_z_SVF, $xyz[2] );
									push( @cloud_t_SVF, $xyz[3] );
									push( @cloud_id_SVF, $this_id );
									push( @cloud_mother_index_list_SVF, -1 ); #default is mother unknown
								} else {
									#print join(',', @xyz ) . "\n";
								}
							} else {
								#push( @this_block_ids, -1 );
								#push( @this_block_names, "" );
							}
							#push( @this_block, $_ );
							if ( $_ =~ /(^|^\s+)<\/SpotsInFrame/ ) {
								$cloud_t_SVF_stop[$this_frame] = scalar(@cloud_x_SVF);
								last;
							}
						}
						#push( @lines_spots_ann, \@this_block );
						#push( @lines_spots_ann_ids, \@this_block_ids );
						#push( @lines_spots_ann_names, \@this_block_names );
					} elsif ( $_ =~ /(^|^\s+)<\/AllSpots/ ) {
						#push( @lines_main, $_ );
						last;
					} else {
						#not inside a SpotsInFrame segment, and yet not finished with AllSpots, so do nothing here put keep reading next line
					}
				}
			} elsif ( $_ =~ /(^|^\s+)<AllTracks/ && $_ !~ /\/>\s*$/ ) {
				#push( @lines_main, $_ );
				#$line_number_tracks = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					#if ( $_ =~ /(^|^\s+)<\/AllTracks/ ) {
					#	#push( @lines_main, $_ );
					#	last;
					#} else {
					#	#not inside a Track segment, and yet not finished with AllTracks, so do nothing here put keep reading next line
					#}
					if ( $_ =~ /(^|^\s+)<Track.*?TRACK_ID=\"(\d+)\"/ ) {
						
						#my $this_track_id = $2;
						#my @this_block = ( $_ );
						#print "inside here: $this_track_id\n";
						while( <FILE> ) {
							chomp;
							#push( @this_block, $_ );
							if ( $_ =~ /(^|^\s+)<\/Track/ ) {
								last;
							} elsif ($_ =~ /(^|^\s+)<Edge/) { #store all spots associated with this track for reconstruction
								if ( $_ =~ /SPOT_SOURCE_ID=\"(\d+)\".*?SPOT_TARGET_ID=\"(\d+)\"/ ) {
									#print "here\n";
									#my $source = 1; my $target = 2;
									#push( @{$spots_in_tracks_complete[scalar(@lines_tracks)]}, $1, $2 );
									$cloud_mother_index_list_SVF[$cloud_self_index_list{$2}] = $cloud_self_index_list{$1};
								}
							}
						}
						#@{$spots_in_tracks[scalar(@lines_tracks)]} = uniq(@{$spots_in_tracks_complete[scalar(@lines_tracks)]});
						#$track_in_tracks[scalar(@lines_tracks)] = $this_track_id;
						#$bins_tracks[scalar(@lines_tracks)] = -1;
						#push( @lines_tracks, \@this_block );
					} elsif ( $_ =~ /(^|^\s+)<\/AllTracks/ ) {
						#push( @lines_tracks, [ ($_) ] );
						#push( @lines_main, $_ );
						last;
					}
				}
			} elsif ( $_ =~ /(^|^\s+)<FilteredTracks/ && $_ !~ /\/>\s*$/ ) {
				#push( @lines_main, $_ );
				#$line_number_filters = scalar(@lines_main);
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<\/FilteredTracks/ ) {
						#push( @lines_main, $_ );	
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
	
	my $max_t = max( @cloud_t_SVF );
	my $min_t = min( @cloud_t_SVF );
	my @X_SVF = set_stats_plus( @cloud_x_SVF );
	my @Y_SVF = set_stats_plus( @cloud_y_SVF );
	my @Z_SVF = set_stats_plus( @cloud_z_SVF );
	my $std_X = $X_SVF[2]/ 2;
	my $std_Y = $Y_SVF[2]/ 2;
	my $std_Z = $Z_SVF[2];# / 2;
	
	#print $std_X . "," . $std_Y . "," . $std_Z . "\n"; return;
	
	#parse dataset.mamut
	if ( open(FILE, "<$path/$dataset_file" ) ) {
		print "Writing track coordinates in time for $dataset_file, subtracting $dataset_file_ann movement, to a new tab-delimited file.\n";
		flock(FILE, LOCK_EX);
		while( <FILE> ) {
			chomp;
			#if ( $_ =~ /^\s*<.+>$/ ) { #one tag per line encountered on this line
			if ( $_ =~ /(^|^\s+)<AllSpots/ ) {
				push( @lines_main, $_ );
				$line_number_spots = scalar(@lines_main);
				#push( @lines_spots, [ ($_) ] );
				
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<SpotsInFrame.*frame=\"(\d+)\"/i ) {
						push( @frames_spots, $2 );
						my @this_block = ( $_ );
						while( <FILE> ) {
							chomp;
							push( @this_block, $_ );
							last if ( $_ =~ /(^|^\s+)<\/SpotsInFrame/ );
						}
						push( @lines_spots, \@this_block );
					} elsif ( $_ =~ /(^|^\s+)<\/AllSpots/ ) {
						push( @lines_main, $_ );
						#push( @lines_spots, [ ($_) ] );
						last;
					} else {
						#not inside a SpotsInFrame segment, and yet not finished with AllSpots, so do nothing here put keep reading next line
					}
				}
			} elsif ( $_ =~ /(^|^\s+)<AllTracks/ && $_ !~ /\/>\s*$/ ) {
				push( @lines_main, $_ );
				$line_number_tracks = scalar(@lines_main);
				#push( @lines_tracks, [ ($_) ] );
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<Track.*?TRACK_ID=\"(\d+)\"/ ) {
						
						my $this_track_id = $2;
						my @this_block = ( $_ );
						#print "inside here: $this_track_id\n";
						while( <FILE> ) {
							chomp;
							push( @this_block, $_ );
							if ( $_ =~ /(^|^\s+)<\/Track/ ) {
								last;
							} elsif ($_ =~ /(^|^\s+)<Edge/) { #store all spots associated with this track for reconstruction
								if ( $_ =~ /SPOT_SOURCE_ID=\"(\d+)\".*?SPOT_TARGET_ID=\"(\d+)\"/ ) {
									#print "here\n";
									push( @{$spots_in_tracks_complete[scalar(@lines_tracks)]}, $1, $2 );
								}
							}
						}
						@{$spots_in_tracks[scalar(@lines_tracks)]} = uniq(@{$spots_in_tracks_complete[scalar(@lines_tracks)]});
						$track_in_tracks[scalar(@lines_tracks)] = $this_track_id;
						$bins_tracks[scalar(@lines_tracks)] = -1;
						push( @lines_tracks, \@this_block );
					} elsif ( $_ =~ /(^|^\s+)<\/AllTracks/ ) {
						#push( @lines_tracks, [ ($_) ] );
						push( @lines_main, $_ );
						last;
					} else {
						#not inside a Track segment, and yet not finished with AllTracks, so do nothing here put keep reading next line
					}
				}
			} elsif ( $_ =~ /(^|^\s+)<FilteredTracks/ && $_ !~ /\/>\s*$/ ) {
				push( @lines_main, $_ );
				$line_number_filters = scalar(@lines_main);
				#push( @lines_filters, [ ($_) ] );
				while( <FILE> ) {
					chomp;
					if ( $_ =~ /(^|^\s+)<\/FilteredTracks/ ) {
						push( @lines_main, $_ );	
						last;
					} elsif ( $_ =~ /TRACK_ID=\"(\d+)\"/ ) {
						#within FilteredTracks segment
						push( @track_in_filters, $1 );
						push( @lines_filters, $_ );
					}
				}	
			} else {
				push( @lines_main, $_ );
			}

		}
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error opening $path/$dataset_file!\n";
		return;
	}
	#return;
	#strip filename of extension, and assume it is XML
	( undef, $dataset_file ) = split( /\./, reverse($dataset_file), 2 );
	$dataset_file = reverse($dataset_file);

	#iterate over spots at t=0, subdivide them, while storing IDs for tree reconstruction
	my @output_table; #array of array, with each array mirroring track_in_tracks, and storing all coordinates for each timepoint 0...n in [0]...[n] positions of the array
	
	for ( my $l=0; $l<scalar(@track_in_tracks); $l++ ) { #iterate over known tracks
		my $xyzprint;
		my $this_spot_id;
		my @lineages = (); #array to store next linked cell for each simultaneous lineage being traced for a given track at a given timepoint
		my @lineages_xyz = (); #array of array to store @xyz for each mother cell, array is manipulated in parallel with @lineages
		my @lineages_adj_xyz = (); #array of array to store adjusted @xyz for each mother cell, array is manipulated in parallel with @lineages
		
		#my @counter = ( 0,0,0,0 ); #DEBUG
		
		for ( my $i=0; $i<scalar(@lines_spots); $i++ ) { #then iterate over timepoints
			$current_timepoint_positions = "";
			$output_table[$l][$frames_spots[$i]] = ""; #clear current timepoint for this track
			
			for ( my $k=scalar(@{$lines_spots[$i]})-1; $k>=0; $k-- ) { #look at all spots at this timepoint
				if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /ID=\"(\d+)\"/ ) {
					#capture ID
					$this_spot_id = $1;
					my @xyz;
					my @mother_xyz; #recovered from @lineages_xyz when the spot's lineage is found
					my @mother_adj_xyz;
					
					#see if this spot belongs in this track
					my $going = 0;
					for ( my $j=scalar(@{$spots_in_tracks[$l]})-1; $j>=0; $j-- ) {	
						if ( $spots_in_tracks[$l][$j] eq $this_spot_id ) {
							$going = 1;
							last;
						}
					}
					next unless ( $going > 0 );
					
					#capture coordinates
					if ($lines_spots[$i][$k] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ ) {
						$xyz[0] = sprintf("%.3f", $1);
					}					
					if ($lines_spots[$i][$k] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ ) {
						$xyz[1] = sprintf("%.3f", $1);
					}
					if ($lines_spots[$i][$k] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
						$xyz[2] = sprintf("%.3f", $1);
					}
					if ($lines_spots[$i][$k] =~ /RADIUS=\"([+-]?\d+(\.?\d*))\"/ ) {
						$xyz[3] = sprintf("%.3f", $1);
					} else {
						$xyz[3] = "N\/A";
					}
					
					next unless ( defined($xyz[0]) && defined($xyz[1]) && defined($xyz[2]) && $xyz[0] =~ /\d/ && $xyz[1] =~ /\d/ && $xyz[2] =~ /\d/ );
					
					#find nearest neighbor in other dataset
					my $closest_spot_distance = vec3_distance(($cloud_x_SVF[0]-$xyz[0]),($cloud_y_SVF[0]-$xyz[1]),($cloud_z_SVF[0]-$xyz[2]));
					my $closest_spot = 0;
					my $closest_spot_mother = -1;
					my $this_distance;
					for ( my $q=$cloud_t_SVF_start[$frames_spots[$i]]; $q<$cloud_t_SVF_stop[$frames_spots[$i]]; $q++ ) {
						next unless ($cloud_t_SVF[$q] == $frames_spots[$i]);
						#$start_here = $q if ( $start_here != $q );
						#only consider spots within a fourth standard deviation in all three axes
						next if ( abs($cloud_x_SVF[$q]-$xyz[0]) > $std_X );
						next if ( abs($cloud_y_SVF[$q]-$xyz[1]) > $std_Y );
						next if ( abs($cloud_z_SVF[$q]-$xyz[2]) > $std_Z );
					
						$this_distance = vec3_distance(($cloud_x_SVF[$q]-$xyz[0]),($cloud_y_SVF[$q]-$xyz[1]),($cloud_z_SVF[$q]-$xyz[2]));
						if ( defined($this_distance) && $this_distance > 0 && $this_distance < $closest_spot_distance ) {
							$closest_spot_distance = $this_distance;
							$closest_spot = $q;
						}
					}
									
					#now, we have the closest spot in the other dataset, find its mother
					$closest_spot_mother = $cloud_mother_index_list_SVF[$closest_spot];
					$closest_spot_mother = -1 unless ( $closest_spot_mother >= 0 && $closest_spot_mother < scalar(@cloud_x_SVF) );
					
					

					#print "spot in track $l: timepoint $i\n";
					#find daughter cells of this
					my @daughters;
					for ( my $q=0; $q<scalar(@{$spots_in_tracks_complete[$l]}); $q+=2 ) {
						if ( $spots_in_tracks_complete[$l][$q] eq $this_spot_id ) {
							push( @daughters, $spots_in_tracks_complete[$l][$q+1] );
						}
					}
					
					#start new lineage if none exists
					my $current_lineage = "";
					if ( scalar(@lineages) == 0 ) {
						push( @lineages, $this_spot_id );
						push( @lineages_xyz, \@xyz );
						push( @lineages_adj_xyz, \@xyz ); #no adjustments made
						$current_lineage = "0";
					}
					
					#find our lineage
					my @adj_xyz = @xyz;
					
					for ( my $q=scalar(@lineages)-1; $q>=0; $q-- ) {
						if ( $lineages[$q] eq $this_spot_id ) {
							$current_lineage = "$q";
							@mother_xyz = @{$lineages_xyz[$q]};
							@mother_adj_xyz = @{$lineages_adj_xyz[$q]};
							if ( defined($mother_xyz[0]) && defined($mother_xyz[1]) && defined($mother_xyz[2]) && $mother_xyz[0] =~ /\d/ && $mother_xyz[1] =~ /\d/ && $mother_xyz[2] =~ /\d/ ) {
								if ( defined($mother_adj_xyz[0]) && defined($mother_adj_xyz[1]) && defined($mother_adj_xyz[2]) && $mother_adj_xyz[0] =~ /\d/ && $mother_adj_xyz[1] =~ /\d/ && $mother_adj_xyz[2] =~ /\d/ ) {
									if ( $closest_spot_mother >= 0 ) {
										#everything is defined, so we can proceed to adjust @adj_xyz based on @mother_adj_xyz, taking movement from (@xyz-@mother_xyz) minus (@nn_xyz-@nn_mother_xyz)
										$adj_xyz[0] = sprintf("%.3f", $mother_adj_xyz[0] + ( ($xyz[0]-$mother_xyz[0]) - ( $cloud_x_SVF[$closest_spot] - $cloud_x_SVF[$closest_spot_mother] ) ) );
										$adj_xyz[1] = sprintf("%.3f", $mother_adj_xyz[1] + ( ($xyz[1]-$mother_xyz[1]) - ( $cloud_y_SVF[$closest_spot] - $cloud_y_SVF[$closest_spot_mother] ) ) );
										$adj_xyz[2] = sprintf("%.3f", $mother_adj_xyz[2] + ( ($xyz[2]-$mother_xyz[2]) - ( $cloud_z_SVF[$closest_spot] - $cloud_z_SVF[$closest_spot_mother] ) ) );
										#$counter[0]++;
									} else {
										#don't have closest spot data, so just propagate motion directly through without adjustment
										$adj_xyz[0] = sprintf("%.3f", $mother_adj_xyz[0] + ( $xyz[0]-$mother_xyz[0]) );
										$adj_xyz[1] = sprintf("%.3f", $mother_adj_xyz[1] + ( $xyz[1]-$mother_xyz[1]) );
										$adj_xyz[2] = sprintf("%.3f", $mother_adj_xyz[2] + ( $xyz[2]-$mother_xyz[2]) );
										
										#print "no closest spot data... xyz: " . join(",",@xyz) . ",$frames_spots[$i] | clostest_spot_xyz: " .  $cloud_x_SVF[$closest_spot] . "," . $cloud_y_SVF[$closest_spot] . "," . $cloud_z_SVF[$closest_spot] . "," . $cloud_t_SVF[$closest_spot] . "\n";
										#$counter[1]++;
									}
								} else {
									#$counter[2]++;
								}
							} else {
								#$counter[3]++;
							}
							
							if ( scalar(@daughters) > 1 ) {
								my $new_lineage_spot = scalar(@lineages);
								$current_lineage .= "->" . $new_lineage_spot . "+" . ($new_lineage_spot+1);
								$lineages[$new_lineage_spot] = $daughters[0];
								$lineages[$new_lineage_spot+1] = $daughters[1];
								$lineages_xyz[$new_lineage_spot] = \@xyz;
								$lineages_xyz[$new_lineage_spot+1] = \@xyz;
								$lineages_adj_xyz[$new_lineage_spot] = \@adj_xyz;
								$lineages_adj_xyz[$new_lineage_spot+1] = \@adj_xyz;
							} elsif (scalar(@daughters) == 1)  {
								$lineages[$q] = $daughters[0];
								$lineages_xyz[$q] = \@xyz;
								$lineages_adj_xyz[$q] = \@adj_xyz;
							}
							last;
						}
					}
					
					
						#shift @xyz before printing to reflect subtraction of movement from background dataset
					
					
					
					
					#} else { #start of new track, disregard subtraction/shifting
						$xyzprint = "{" . $current_lineage . "," . join(',',@adj_xyz) . "}";
						if ( $output_table[$l][$frames_spots[$i]] eq "" ) { #empty so don't try to find the spot
							$output_table[$l][$frames_spots[$i]] = $xyzprint;
						} else { 
							$output_table[$l][$frames_spots[$i]] .= ";" . $xyzprint; #concatenate output
						}
					#}
				}
			}
		}
		
		#print "track $l, counters: " . join(",",@counter) . "\n";
		if ( rand() < 0.01 ) {
			print "Progress: " . 100*($l/scalar(@track_in_tracks)) . " %\n";
		}
	}

	if ( open(FILE, ">$path/$dataset_file\.track_coordinates_in_time_shifted\.tsv" ) ) {
		flock(FILE, LOCK_EX);
		print FILE "Track";
		my $max_timepoints = max(@frames_spots);
		#print $max_timepoints . "\n";
		for ( my $i=min(@frames_spots); $i<$max_timepoints; $i++ ) {
			print FILE "\t Frame " . $i;
		}
		print FILE "\n";
		for ( my $l=0; $l<scalar(@track_in_tracks); $l++ ) {
			print FILE $track_in_tracks[$l];
			for ( my $k=0; $k<scalar(@{$output_table[$l]}); $k++ ) {
				if ( defined($output_table[$l][$k]) && $output_table[$l][$k] =~ /\w/ ) {
					print FILE "\t" . $output_table[$l][$k];
				} else {
					print FILE "\t";
				}
			}
			print FILE "\n";
		}
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error writing to $path/$dataset_file\.track_coordinates_in_time_shifted\.tsv!\n";
		return;
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

main(@ARGV);


