#!/usr/bin/perl

# MaMuT downsample density
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Downsamples a MaMuT dataset by factor X, going toward uniform cell density, with tracks allocated randomly among the bins
# usage: perl MaMuT_dataset_downsample_density.pl dataset_mamut.xml X
#	where X is the desired fold reduction

# note this script is single-threaded and slow!

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


#==================
#Main subroutine
#==================
sub main {
	#my (  @arguments ) = @_;
	my $dataset_file = $_[0];
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
	my $split_fold = 2;
	$split_fold = $_[1] if ( defined($_[1]) && $_[1] =~ /\d/ );
	$split_fold = 2 if ( $split_fold < 2 );
	my $split_fold_inverse = $split_fold - (1 / $split_fold);
	my @bins_tracks;
	#my $open_block = -1; #will determine actions to take while parsing XML line-by-line
	
	#parse dataset.mamut
	if ( open(FILE, "<$path/$dataset_file" ) ) {
		print "De-densifying dataset $dataset_file track-wise by $split_fold fold.  Can change the fold reduction in parameter that follows dataset filename.\n";
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
	
	##iterate over spots over all time points to determine mean Z coordinate
	#my @this_line;
	##my $num_spots = 0;
	##my $sum_spots =0;
	#my @spot_coordinates;
	#for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
	#	for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) {
	#		if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
	#			#capture $1 to store position Z in string form, then immediately use in float sum
	#			#$sum_spots += $1;
	#			#$num_spots++;
	#			push( @spot_coordinates, $1 );
	#		}
	#	}
	#}
	#
	##write median Z to $mean_z
	#@spot_coordinates = sort { $a <=> $b } @spot_coordinates;#{ $a cmp $b };
	#my $mean_z;
	#if ( scalar(@_) & 1 ) { #odd number of numbers
	#	my $med = int( scalar(@spot_coordinates)/2);
	#	$mean_z = $spot_coordinates[$med];
	#	#$med = "me!";
	#} else {
	#	my $med = scalar(@spot_coordinates) / 2;
	#	$mean_z = ( $spot_coordinates[$med] + $spot_coordinates[$med -1] ) / 2;
	#}	
	##my $mean_z = $sum_spots/$num_spots;
	##print "mean Z: " . $sum_spots/$num_spots . ", " . $num_spots . "\n"; return;
	#unless( defined($mean_z) && $mean_z =~ /\d/ ) {
	#	print "Could not determine median Z from spot coordinates.  Aborting.\n"
	#	return;
	#}
	
	#iterate over spots at t=0, subdivide them, while storing IDs for tree reconstruction
	my @lines_spots_divided; #final storage for subdivisions, 3D matrix [timepoint][subdivision][spot], compared with 2D matrix of @lines_spots [timepoint][spot]
	my @tracking_spots_divided; # storage for spot IDs undergoing lineage reconstruction by the splitting algorithm, 2D matrix [subdivision][spotID], each spotID will be removed once the daughter spot is found and added to @lines_spots_divided
	my $this_spot_id;
	my $this_bin;
	my @cell_positions = (); #array of 3d vectors, one for each cell
	my $distance_numbers;
	my $distance_sum;
		my $this_distance;
		my $this_spot_density;
		my @cell_densities = ();
		my @all_densities = (); #running total for figuring out how to normalize
		#my $all_densities_begin_this_timepoint;
		#my @ranked_densities_this_timepoint;
	#print "t=0, total spots handled: " . scalar(@{$lines_spots[0]}) . "\n";
	my @timepoint_density_cutoffs = ();
	for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
			for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) { #find all coordinates first
				if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ ) {
					my @this_cell_coordinates = ( $1 );
					if ( $lines_spots[$i][$k] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ ) {
						push( @this_cell_coordinates, $1 );
					}
					if ( $lines_spots[$i][$k] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
						push( @this_cell_coordinates, $1 );
					}
					if ( $lines_spots[$i][$k] =~ /RADIUS=\"([+-]?\d+(\.?\d*))\"/ ) {
						push( @this_cell_coordinates, $1 * 10 );
					}					
					if ( scalar(@this_cell_coordinates) >= 3 && $this_cell_coordinates[0] =~ /\d/ && $this_cell_coordinates[1] =~ /\d/ && $this_cell_coordinates[2] =~ /\d/ ) {
						$cell_positions[$i][$k] = \@this_cell_coordinates;
					} else {
						$cell_positions[$i][$k] = undef;
					}
				} else {
					$cell_positions[$i][$k] = undef;
				}
			}
			#$all_densities_begin_this_timepoint = scalar(@all_densities);
			@all_densities = ();
			for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) { #now figure out average distance of other cells, for each cell
				if ( $lines_spots[$i][$k] =~ /<Spot/ && ref($cell_positions[$i][$k]) ) { #got a valid cell
					$distance_sum =0;
					$distance_numbers =0;
					
					for ( my $kk=0; $kk<scalar(@{$lines_spots[$i]}); $kk++ ) {
						next if ( $kk == $k || !(ref($cell_positions[$i][$kk])) ); #only accept valid non-self cells beyond this point
						$this_distance = (($cell_positions[$i][$k][0]-$cell_positions[$i][$kk][0]) ** 2  + ($cell_positions[$i][$k][1]-$cell_positions[$i][$kk][1]) ** 2 + ($cell_positions[$i][$k][2]-$cell_positions[$i][$kk][2]) ** 2);
						
						#if ( $this_distance < 
							$distance_sum += $this_distance;
							$distance_numbers++;
						#}
					}
						
					if ( $distance_numbers > 0 && $distance_sum > 0 ) {
						$this_spot_density = $distance_sum/$distance_numbers;

						push( @all_densities, $this_spot_density );
						$cell_densities[$i][$k] = $this_spot_density;
					} else {
						$cell_densities[$i][$k] = "N\/A";
					}

				} else {
					$cell_densities[$i][$k] = "N\/A";
				}
			}

			if ( $#all_densities > 10) {
				@all_densities = sort{ $a <=> $b } @all_densities;
				$timepoint_density_cutoffs[$i] = $all_densities[$#all_densities-int($#all_densities/$split_fold)]
				
			} else {
				#nothing to do, there are no densities reported for time point $i, just choose some arbitrary high number
				$timepoint_density_cutoffs[$i] = 1000000;
			}
	}
	undef @cell_positions;
	for ( my $k=scalar(@{$lines_spots[0]})-1; $k>=0; $k-- ) {
		if ( $lines_spots[0][$k] =~ /<Spot/ && $lines_spots[0][$k] =~ /ID=\"(\d+)\"/ ) {
			#capture ID
			$this_spot_id = $1;
			#if ($lines_spots[0][$k] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ ) {
			if ($cell_densities[0][$k] =~ /\d/) { #got a valid cell
				$this_bin = 1; #throw away by default
				$this_bin = 0 if ( $cell_densities[0][$k] =~ /\d/ &&  $cell_densities[0][$k] > $timepoint_density_cutoffs[0] * rand() * $split_fold_inverse );
				
				unshift( @{$lines_spots_divided[0][$this_bin]}, $lines_spots[0][$k] );
				#unshift( @{$tracking_spots_divided[$this_bin]}, $this_spot_id );
				#print "about to iterate over " . scalar(@spots_in_tracks) . "\n";
				#now find any tracks belonging to this spot, place track in the the divided bin, and add all associated spots to the @tracking_spots_divided bin
				for ( my $l =0; $l<scalar(@spots_in_tracks); $l++ ) {
					#print "  about to iterate over " . scalar(@{$spots_in_tracks[$l]}) . "\n";
					for ( my $j=scalar(@{$spots_in_tracks[$l]})-1; $j>=0; $j-- ) {
						#print "testing $spots_in_tracks[$l][$j] eq $this_spot_id\n";
						if ( $spots_in_tracks[$l][$j] eq $this_spot_id ) {
							splice( @{$spots_in_tracks[$l]}, $j, 1 ); #already dealt with this spot, so remove it
							#push( @{$tracks_divided[$this_bin]}, splice(@lines_tracks, $l, 1) );
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
				unshift( @{$lines_spots_divided[0][$l]}, $lines_spots[0][$k] );
			}
		}
	}
	undef @{$lines_spots[0]}; #start freeing space
	#print "join " . join("\n", @{$lines_spots[0]}) . "\n";
	#print "num spots in each division: " . scalar(@{$lines_spots_divided[0][0]}) . ", " .  @{$lines_spots_divided[0][1]} . "\n";
	#return;
	#now, iterate over all spots in time, reconstructing lineages and sorting tracks to the correct bins
	my $this_track; my $this_spot_in_track;
	my $going;
	for ( my $i=1; $i<scalar(@lines_spots); $i++ ) {
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
							#push( @{$tracks_divided[$this_bin]}, splice(@lines_tracks, $l, 1) );
							#print "    found spot $this_spot_id in track #$l, track $track_in_tracks[$l], position $j\n";
							$this_bin = $bins_tracks[$l];
							unshift( @{$lines_spots_divided[$i][$this_bin]}, $lines_spots[$i][$k] );
							$going = TRUE;
							last;
						}
					}
					last if ( $going );
				}
				
				if ( $this_bin < 0 ) {
					$this_bin = 1; #throw away by default
					$this_bin = 0 if ( $cell_densities[$i][$k] =~ /\d/ && $cell_densities[$i][$k] > $timepoint_density_cutoffs[$i] * rand() * $split_fold_inverse );
					unshift( @{$lines_spots_divided[$i][$this_bin]}, $lines_spots[$i][$k] );
					#unshift( @{$tracking_spots_divided[$this_bin]}, $this_spot_id );
					
					#now find any tracks belonging to this spot, place track in the the divided bin, and add all associated spots to the @tracking_spots_divided bin
					for ( my $l =0; $l<scalar(@spots_in_tracks); $l++ ) {
						for ( my $j=scalar(@{$spots_in_tracks[$l]})-1; $j>=0; $j-- ) {
							if ( $spots_in_tracks[$l][$j] eq $this_spot_id ) {
								splice( @{$spots_in_tracks[$l]}, $j, 1 ); #already dealt with this spot, so remove it
								#push( @{$tracks_divided[$this_bin]}, splice(@lines_tracks, $l, 1) );
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
	
	#divide the tracks according to bin numbers assigned above -- most likely this code section could be incorporated into file out below, but for now will leave it here
#	my @tracks_divided; #final storage for tracks [subdivision][tracks]
#	my @filters_divided;
#	for ( my $i=scalar(@bins_tracks)-1; $i>=0; $i-- ) {
#		push( @{$tracks_divided[$bins_tracks[$i]]}, splice(@lines_tracks, $i, 1) );
#		for ( my $k=0; $k<scalar(@lines_filters); $k++ ) { #dont worry about working backwards and splicing because this segment takes up very little memory
#			push( @{$filters_divided[$bins_tracks[$i]]}, $lines_filters[$k] );		
#		}
#	}
	
	#finally, for each subdivision, update nspots and write file
	$num_bins = 1; #only print the de-densified dataset
	my $nspots_divided = ( 0 ) x $num_bins;
	for (my $l=0; $l<$num_bins; $l++ ) {
		for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
			$nspots_divided[$l] += scalar(@{$lines_spots_divided[$i][$l]})-2;
		}
		for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) { #look to find and update nspots
			if ( $lines_main[$fl] =~ /(^|^\s+)<AllSpots/ ) {
				$lines_main[$fl] =~ s/nspots=\"(\d+)\"/nspots=\"$nspots_divided[$l]\"/;
			}
		}
		
		if ( open(FILE, ">$path/$dataset_file\.dedensified$split_fold\.xml" ) ) {
			flock(FILE, LOCK_EX);
			for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) {
				if ( $fl == $line_number_spots ) {
					for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
						for ( my $k=0; $k<scalar(@{$lines_spots_divided[$i][$l]}); $k++ ) {
							print FILE $lines_spots_divided[$i][$l][$k] . "\n";
						}
					}
				} elsif ( $fl == $line_number_tracks ) {
					#for ( my $i=0; $i<scalar(@{$tracks_divided[$l]}); $i++ ) {
					#	print FILE $tracks_divided[$l][$i] . "\n";
					#}
					for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						if ( $bins_tracks[$i] == $l ) {
							print FILE join( "\n", @{$lines_tracks[$i]} ) . "\n";
						}
					}
				} elsif ( $fl == $line_number_filters ) {
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
					#for ( my $i=0; $i<scalar(@{$filters_divided[$l]}); $i++ ) {
					#	print FILE $filters_divided[$l][$i] . "\n";
					#}
				} #else {
				
				print FILE $lines_main[$fl] . "\n";
				#}
			}
			flock(FILE, LOCK_UN);
			close(FILE);
		} else {
			print "Error writing to $path/$dataset_file\.dedensified$split_fold\.xml!\n";
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


main(@ARGV);


