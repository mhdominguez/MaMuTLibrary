#!/usr/bin/perl

# MaMuT split dataset by XYZ coordinates into bins, backwards in time
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Splits MaMuT dataset into n bins by dividing it into mean/stdev or percentile bins along a particular axis, using chronologic last cell in the track to assign the bin
# usage: perl MaMuT_dataset_split_XYZ_coordinate_backward_in_time.pl mamut_dataset.xml bin=2 med_z #med_z can be replaced with mean_x, mean_y, mean_z, med_x, or med_y



use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;




use constant {
	SPLIT_FEATURE_TRACK_ORIGIN_MEAN_Z => 0,
	SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_Z => 1,
	SPLIT_FEATURE_TRACK_ORIGIN_MEAN_Y => 2,
	SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_Y => 3,
	SPLIT_FEATURE_TRACK_ORIGIN_MEAN_X => 4,
	SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_X => 5,
	
	SPLIT_FEATURE_TRACK_MAX => 6,
};



#==================
#Main subroutine
#==================
sub main {
	#my (  @arguments ) = @_;
	my $dataset_file = shift;
	my @params = @_;
	my @lines_main; #unchanging lines during reconstruction
	my $line_number_spots = -1; #indicates where to put <AllSpots> blocks when reconstructing the file at the end
	my @lines_spots;
	my @frames_spots;
	
	my $line_number_tracks = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_tracks;
	#my @track_tracks;
	my @spots_in_tracks;
	my @track_in_tracks;
	my $line_number_filters = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_filters;
	my @track_in_filters;
	
	my $num_bins = 2;
	my @bins_tracks;
	my $overlapping_bins = FALSE;;

	for ( my $k=0; $k<scalar(@params); $k++ ) {
		$params[$k] =~ s/==/=/g;
		@this_param = split( /=/, $params[$k] );

		if ( $this_param[0] =~ /bin/i ) {
			$num_bins = $this_param[1] if (defined($this_param[1]) && $this_param[1] =~ /\d/ );
			next;
		} elsif ( $this_param[0] =~ /overl/i) {
			$overlapping_bins = TRUE;
			next;
		}
		
		$track_filter_index = -1;
		if ( $this_param[0] =~ /av/i || $this_param[0] =~ /mean/i ) {
			if ( $this_param[0] =~ /x/i ) {
				$track_filter_index = SPLIT_FEATURE_TRACK_ORIGIN_MEAN_X;
			} elsif ( $this_param[0] =~ /y/i ) {
				$track_filter_index = SPLIT_FEATURE_TRACK_ORIGIN_MEAN_Y;
			} elsif ( $this_param[0] =~ /z/i ) {
				$track_filter_index = SPLIT_FEATURE_TRACK_ORIGIN_MEAN_Z;
			}
		} elsif ( $this_param[0] =~ /med/i ) {
			if ( $this_param[0] =~ /x/i ) {
				$track_filter_index = SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_X;
			} elsif ( $this_param[0] =~ /y/i ) {
				$track_filter_index = SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_Y;
			} elsif ( $this_param[0] =~ /z/i ) {
				$track_filter_index = SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_Z;
			}
		} 	

	}
	undef @this_param;	
	
	if ( $num_bins < 2 || $track_filter_index < 0 ) {
		print "Invalid number of bins $num_bins, or split feature direction ($track_filter_index, i.e. x_mean, y_median) indicated!  Cannot proceed.\n";
		return;
	} else {
		if ( $overlapping_bins ) {
			print "Going to split dataset by code #$track_filter_index into $num_bins bins, overlapping...\n";
			$num_bins++;
		} else {
			print "Going to split dataset by code #$track_filter_index into $num_bins bins...\n";
		}
	}
	
	#parse dataset.mamut
	if ( open(FILE, "<$path/$dataset_file" ) ) {
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
									push( @{$spots_in_tracks[scalar(@lines_tracks)]}, $1, $2 );
								}
							}
						}
						#print "scalar " . scalar(@{$spots_in_tracks[scalar(@lines_tracks)]}) . "\n";
						@{$spots_in_tracks[scalar(@lines_tracks)]} = uniq(@{$spots_in_tracks[scalar(@lines_tracks)]});
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
	#print join( "\n", @lines_main );
	#print join( "\n", @track_in_tracks );
	#print "Number of main lines: " . scalar(@lines_main) . ", number of spots lines: " . scalar(@lines_spots) . ", number of tracks lines: " . scalar(@lines_tracks) . ", number of filtered tracks lines: " . scalar(@lines_filters) . "\n";
	
	#return;
	
	#iterate over spots over all time points to determine mean Z coordinate
	my @this_line;
	#my $num_spots = 0;
	#my $sum_spots =0;
	my @spot_coordinates;
	my $this_axis = "";
	if ( $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_X || $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEAN_X ) {
		for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
			for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) {
				if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ ) {
					#capture $1 to store position Z in string form, then immediately use in float sum
					push( @spot_coordinates, $1 );
				}
			}
		}
		$this_axis = "X";
	} elsif ( $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_Y || $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEAN_Y ) {
		for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
			for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) {
				if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ ) {
					#capture $1 to store position Z in string form, then immediately use in float sum
					push( @spot_coordinates, $1 );
				}
			}
		}
		$this_axis = "Y";
	} elsif ( $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_Z || $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEAN_Z ) {
		for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
			for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) {
				if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
					#capture $1 to store position Z in string form, then immediately use in float sum
					push( @spot_coordinates, $1 );
				}
			}
		}
		$this_axis = "Z";
	}
	
	my @split_maxes; #stores max/min values for each bin cells
	my @split_mins;
	if ( scalar(@spot_coordinates) < 2 ) {
		print "Not enough spots!\n";
		return;
	}
	if ( $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEAN_X || $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEAN_Y || $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEAN_Z ) {
		#print $track_filter_index."\n"; return;
		#for two-bin splitting, this means mean, otherwise will split by number of standard deviations (above and below); for example, three will be three bins with middle between one standard dev, 
		my $min = min(@spot_coordinates );
		my $max = max(@spot_coordinates );
		my ( $avg, $med, $std ) = set_stats(@spot_coordinates);
		
		#obvious cases of split vectors
		$split_mins[0] = $min;
		$split_maxes[$num_bins-1] = $max;		
		
		if ( $num_bins < 3 ) { #special case == 2, two bins split at average
			$split_maxes[0] = $split_mins[1] = $avg;
		} elsif ( $num_bins % 2 == 0 ) {
			my $standard_devations = $num_bins/2 -1;
			my $middle_plus_one_position = $num_bins/2;
			#even bin number, split middle bins around the average and one standard devation
			
			#set up first and last bins
			$split_maxes[0] = $avg - ( $std * $standard_devations );
			$split_mins[$num_bins-1] = $avg + ( $std * $standard_devations );
			
			#then iterate on all intermediate bins
			for (my $i=0; $i<$standard_devations; $i++ ) {
				#$i maps to positions in @split_mins/@split_maxes
				my $high = $middle_plus_one_position + $i;
				my $low = $standard_devations - $i;

				#now fill high and low slots
				$split_maxes[$high] = $avg + ( $std * ($i+1) );
				$split_mins[$high] = $avg + ( $std * $i );
				$split_maxes[$low] = $avg - ( $std * $i );
				$split_mins[$low] = $avg - ( $std * ($i+1) );	
			}	
			print " ...splitting by average plus or minus $standard_devations from the mean...\n";
		} else { #odd bin number
			my $standard_devations = int($num_bins/2);
			
			#set up first and last bins
			$split_maxes[0] = $avg - ( $std * $standard_devations );
			$split_mins[$num_bins-1] = $avg + ( $std * $standard_devations );
			
			#populate std_dev matrix for binning
			my $this_bin = 1;	
			for ( my $k=1-$standard_devations; $k<$standard_devations; $k++ ) {
				$split_maxes[$this_bin] = $avg + ( $std * $k );
				$split_mins[$this_bin] = $avg - ( $std * $k );
				$this_bin++;
			}
			print " ...splitting by average plus or minus $standard_devations from the mean...\n";
		}
		
		#double check all bin mins/maxes, make sure they make sense
		for ( my $i=0; $i<$num_bins/2+1; $i++ ) {
			if ( $split_mins[$i] < $min ) {
				$split_mins[$i] = $min;
				if ( $split_maxes[$i] < $min ) {
					$split_maxes[$i] = $min;
				}
			}
		}
		for ( my $i=$num_bins-1; $i>=0; $i-- ) {
			if ( $split_maxes[$i] > $max ) {
				$split_maxes[$i] = $max;
				if ( $split_mins[$i] > $max ) {
					$split_mins[$i] = $max;
				}
			}
		}
		for ( my $i=0; $i<$num_bins; $i++ ) { #final fail-safe
			if ( $split_maxes[$i] < $split_mins ) {
				print "   ...odd splitting binning error occurred, bin $i, range $split_mins[$i] to $split_maxes[$i].\n";
				$split_maxes[$i] = $split_mins[$i];
			}
		}
		print "   ...bin boundaries are: [". join(',', @split_mins) . "] to [" . join(',', @split_maxes) . "].\n";
	} elsif ( $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_X || $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_Y || $track_filter_index == SPLIT_FEATURE_TRACK_ORIGIN_MEDIAN_Z ) {
		#for two-bin splitting, this means median, otherwise will use histogram equally-weight bins	
		@spot_coordinates = sort { $a <=> $b } @spot_coordinates;#{ $a cmp $b };
		#my $mean_z;
		my $bin_size = int(scalar(@spot_coordinates)/$num_bins);
		
		#populate split_maxes and split_mins;
		my $position;
		for ( my $i=0; $i<$num_bins; $i++ ) {
			$position = $i*$bin_size;
			$split_mins[$i] = $spot_coordinates[$position];
			$split_maxes[$i] = $spot_coordinates[$position+$bin_size]; 
			if ( $i == $num_bins -1 ) {
				$split_maxes[$i] = $spot_coordinates[$#spot_coordinates];			
			}
		}

		print " ...splitting by percentile rank, ";
		print "bin boundaries are: [". join(',', @split_mins) . "] to [" . join(',', @split_maxes) . "].\n";
	}
	
	#iterate over spots at t=0, subdivide them, while storing IDs for tree reconstruction
	my @lines_spots_divided; #final storage for subdivisions, 3D matrix [timepoint][subdivision][spot], compared with 2D matrix of @lines_spots [timepoint][spot]
	my @tracking_spots_divided; # storage for spot IDs undergoing lineage reconstruction by the splitting algorithm, 2D matrix [subdivision][spotID], each spotID will be removed once the daughter spot is found and added to @lines_spots_divided
	my $this_spot_id;
	my $this_bin;
	my $captured;
	
	#print "t=0, total spots handled: " . scalar(@{$lines_spots[0]}) . "\n";
	my $index = scalar(@{$lines_spots[$#lines_spots]})-1;
	for ( my $k=$index; $k>=0; $k-- ) {
		if ( $lines_spots[$index][$k] =~ /<Spot/ && $lines_spots[$index][$k] =~ /ID=\"(\d+)\"/ ) {
			#capture ID
			$this_spot_id = $1;
			if ($lines_spots[$index][$k] =~ /POSITION_${this_axis}=\"([+-]?\d+(\.?\d*))\"/ ) {
				$captured = $1;
				$this_bin - 1;
				for ( my $ii=0; $ii<$num_bins; $ii++ ) {
					if ( $captured <= $split_maxes[$ii] && $captured >= $split_mins[$ii] ) {
						$this_bin = $ii;
						last;
					}
				}
				if ( $this_bin < 0 ) {
					print "Could not find appropriate bin for spot $this_spot_id on line $k of time 0\n";
				}

				unshift( @{$lines_spots_divided[$index][$this_bin]}, $lines_spots[$index][$k] );
				#print "about to iterate over " . scalar(@spots_in_tracks) . "\n";
				#now find any tracks belonging to this spot, place track in the the divided bin, and add all associated spots to the @tracking_spots_divided bin
				for ( my $l =0; $l<scalar(@spots_in_tracks); $l++ ) {
					#print "  about to iterate over " . scalar(@{$spots_in_tracks[$l]}) . "\n";
					for ( my $j=scalar(@{$spots_in_tracks[$l]})-1; $j>=0; $j-- ) {
						#print "testing $spots_in_tracks[$l][$j] eq $this_spot_id\n";
						if ( $spots_in_tracks[$l][$j] eq $this_spot_id ) {
							splice( @{$spots_in_tracks[$l]}, $j, 1 ); #already dealt with this spot, so remove it
							$bins_tracks[$l] = $this_bin;
							#print "assigning $this_bin to track $track_in_tracks[$l]\n";
							last;
						}
					}
				}

			}
		} else {
			#control blocks get passed to both divisions
			for (my $l=0; $l<$num_bins; $l++ ) {
				unshift( @{$lines_spots_divided[$index][$l]}, $lines_spots[$index][$k] );
			}
		}
	}
	undef @{$lines_spots[$index]}; #start freeing space

	#now, iterate over all spots in time, reconstructing lineages and sorting tracks to the correct bins
	my $this_track; my $this_spot_in_track;
	my $going;
	for ( my $i=scalar(@lines_spots)-2; $i>=0; $i-- ) {
		#print "t=$i, total spots handled: " . scalar(@{$lines_spots[$i]}) . "\n";
		for ( my $k=scalar(@{$lines_spots[$i]})-1; $k>=0; $k-- ) {
			if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /ID=\"(\d+)\"/ ) {
				#capture ID
				$this_spot_id = $1;
				$this_bin = -1;
				#print "  setting up $this_spot_id\n";
				#try to find this ID within tracks for either bin already decided
				$going = FALSE;
				for ( my $l=0; $l<scalar(@spots_in_tracks); $l++ ) {
					#print "   looking in track $track_in_tracks[$l] ($l)...\n";
					next if ($bins_tracks[$l] < 0);
					#print "   looking in track $track_in_tracks[$l] ($l)...\n";
					for ( my $j=scalar(@{$spots_in_tracks[$l]})-1; $j>=0; $j-- ) {
						if ( $spots_in_tracks[$l][$j] eq $this_spot_id ) {			
							splice( @{$spots_in_tracks[$l]}, $j, 1 ); #already dealt with this spot, so remove it
							#print "    found spot $this_spot_id in track #$l, track $track_in_tracks[$l], position $j\n";
							$this_bin = $bins_tracks[$l];
							unshift( @{$lines_spots_divided[$i][$this_bin]}, $lines_spots[$i][$k] );
							$going = TRUE;
							last;
						}
					}
					last if ( $going );
				}
				
				if ( $this_bin < 0 && $lines_spots[$i][$k] =~ /POSITION_${this_axis}=\"([+-]?\d+(\.?\d*))\"/ ) {
					$captured = $1;
					$this_bin - 1;
					for ( my $ii=0; $ii<$num_bins; $ii++ ) {
						if ( $captured <= $split_maxes[$ii] && $captured >= $split_mins[$ii] ) {
							$this_bin = $ii;
							last;
						}
					}
					if ( $this_bin < 0 ) {
						print "Could not find appropriate bin for spot $this_spot_id on line $k of time $i\n";
					}
					unshift( @{$lines_spots_divided[$i][$this_bin]}, $lines_spots[$i][$k] );
					
					#now find any tracks belonging to this spot, place track in the the divided bin, and add all associated spots to the @tracking_spots_divided bin
					for ( my $l =0; $l<scalar(@spots_in_tracks); $l++ ) {
						for ( my $j=scalar(@{$spots_in_tracks[$l]})-1; $j>=0; $j-- ) {
							if ( $spots_in_tracks[$l][$j] eq $this_spot_id ) {
								splice( @{$spots_in_tracks[$l]}, $j, 1 ); #already dealt with this spot, so remove it
								$bins_tracks[$l] = $this_bin;
								last;
							}
						}
					}
				}
			} else {
				#control blocks get passed to both divisions
				for (my $l=0; $l<$num_bins; $l++ ) {
					unshift( @{$lines_spots_divided[$i][$l]}, $lines_spots[$i][$k] );
				}
			}
		}
		undef @{$lines_spots[$i]}; #freeing space
	}
	#return;
	
	#finally, for each subdivision, update nspots and write file
	my $nspots_divided = ( 0 ) x $num_bins;
	for (my $l=0; $l<$num_bins; $l++ ) {

		for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
			$nspots_divided[$l] += scalar(@{$lines_spots_divided[$i][$l]})-2;
		}
		if ( $overlapping_bins && $l == 0 ) {
			next; #we don't do anything with bin 0 with overlapping bins, we start with bin 1 and always include previous bin
		}
		
		if ( $overlapping_bins ) {
			for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) { #look to find and update nspots
				if ( $lines_main[$fl] =~ /(^|^\s+)<AllSpots/ ) {
					my $number_spots = $nspots_divided[$l] + $nspots_divided[$l-1];
					$lines_main[$fl] =~ s/nspots=\"(\d+)\"/nspots=\"$number_spots\"/;
				}
			}
		} else {
			for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) { #look to find and update nspots
				if ( $lines_main[$fl] =~ /(^|^\s+)<AllSpots/ ) {
					$lines_main[$fl] =~ s/nspots=\"(\d+)\"/nspots=\"$nspots_divided[$l]\"/;
				}
			}
		}

		if ( open(FILE, ">$path/$dataset_file\.$l" . "rev" . ".xml" ) ) {
			flock(FILE, LOCK_EX);
			for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) {
				if ( $fl == $line_number_spots ) {
					if ( $overlapping_bins ) {
						for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
							#hack to make this work, assuming first and last lines are control lines
							my $m = $l-1;
							print FILE $lines_spots_divided[$i][$m][0] . "\n";
							for ( my $k=scalar(@{$lines_spots_divided[$i][$m]})-2; $k>0; $k-- ) {
								print FILE $lines_spots_divided[$i][$m][$k] . "\n";
							}
							$m = $l;
							for ( my $k=scalar(@{$lines_spots_divided[$i][$m]})-2; $k>0; $k-- ) {
								print FILE $lines_spots_divided[$i][$m][$k] . "\n";
							}
							print FILE $lines_spots_divided[$i][$m][scalar(@{$lines_spots_divided[$i][$m]})-1] . "\n";
						}
					} else {
						for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
							for ( my $k=0; $k<scalar(@{$lines_spots_divided[$i][$l]}); $k++ ) {
								print FILE $lines_spots_divided[$i][$l][$k] . "\n";
							}
						}
					}
				} elsif ( $fl == $line_number_tracks ) {
					if ( $overlapping_bins ) {
						for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
							if ( $bins_tracks[$i] == $l-1 ) {
								print FILE join( "\n", @{$lines_tracks[$i]} ) . "\n";
							}
						}
					}					
					for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						if ( $bins_tracks[$i] == $l ) {
							print FILE join( "\n", @{$lines_tracks[$i]} ) . "\n";
						}
					}
				} elsif ( $fl == $line_number_filters ) {
					if ( $overlapping_bins ) {
						for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
							if ( $bins_tracks[$i] == $l-1 ) {
								for ( my $k=0; $k<scalar(@lines_filters); $k++ ) {
									if ( $track_in_tracks[$i] eq $track_in_filters[$k] ) {
										print FILE $lines_filters[$k] . "\n";
										last;
									}
								}
							}
						}
					}
					for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						if ( $bins_tracks[$i] == $l ) {
							for ( my $k=0; $k<scalar(@lines_filters); $k++ ) {
								if ( $track_in_tracks[$i] eq $track_in_filters[$k] ) {
									print FILE $lines_filters[$k] . "\n";
									last;
								}
							}
						}
					}	
				}
				print FILE $lines_main[$fl] . "\n";
			}
			flock(FILE, LOCK_UN);
			close(FILE);
		} else {
			print "Error writing to $path/$dataset_file\.$l" . "rev" . ".xml!\n";
			return;
		}		
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

main(@ARGV);


