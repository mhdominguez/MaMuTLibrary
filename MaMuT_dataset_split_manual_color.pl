#!/usr/bin/perl

# MaMuT split dataset by manual color (or tissue type)
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Splits a MaMuT dataset into separate datasets, with tracks allocated by MANUAL_COLOR or TISSUE_TYPE annotations
# usage: perl MaMuT_dataset_split_manual_color.pl dataset_mamut.xml

# note this script is single-threaded and slow!

use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;


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
	
	my $num_bins = 1;
	my @bin_labels = ("null");
	#$num_bins = $_[1] if ( defined($_[1]) && $_[1] =~ /\d/ );
	my @bins_tracks;
	#my $open_block = -1; #will determine actions to take while parsing XML line-by-line
	
	#parse dataset.mamut
	my %colors_bins; #hash that stores all the important data	
	if ( open(FILE, "<$path/$dataset_file" ) ) {
		print "Splitting dataset $dataset_file by manual color (or tissue type if no color markers present).\n";
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
							if ( $_ =~ /MANUAL_COLOR=\"([+-]?\d+(\.?\d*))\"/ || $_ =~ /TISSUE_TYPE=\"([+-]?\d+(\.?\d*))\"/ ) {
								my $this_color = $1;
								
								unless ( exists($colors_bins{$this_color}) ) {
									$colors_bins{$this_color} = $num_bins;
									$num_bins++;
									push( @bin_labels, int(abs($this_color)) );
								}
							}
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

	#iterate over spots at t=0, subdivide them, while storing IDs for tree reconstruction
	my @lines_spots_divided; #final storage for subdivisions, 3D matrix [timepoint][subdivision][spot], compared with 2D matrix of @lines_spots [timepoint][spot]
	my @tracking_spots_divided; # storage for spot IDs undergoing lineage reconstruction by the splitting algorithm, 2D matrix [subdivision][spotID], each spotID will be removed once the daughter spot is found and added to @lines_spots_divided
	my $this_spot_id;
	my $this_bin;

	#print "t=0, total spots handled: " . scalar(@{$lines_spots[0]}) . "\n";
	for ( my $k=scalar(@{$lines_spots[0]})-1; $k>=0; $k-- ) {
		if ( $lines_spots[0][$k] =~ /<Spot/ && $lines_spots[0][$k] =~ /ID=\"(\d+)\"/ ) {
			#capture ID
			$this_spot_id = $1;
			if ($lines_spots[0][$k] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
				$this_color = "";
				if ($lines_spots[0][$k] =~ /MANUAL_COLOR=\"([+-]?\d+(\.?\d*))\"/ || $lines_spots[0][$k] =~ /TISSUE_TYPE=\"([+-]?\d+(\.?\d*))\"/  ) {
					$this_color = $1;
				}
				if ( $this_color eq "" ) {
					#could not find a color
					$this_bin = 0;
				} elsif ( exists($colors_bins{$this_color}) ) {
					$this_bin = $colors_bins{$this_color};
				} else {
					$colors_bins{$this_color} = $num_bins;
					$this_bin = $num_bins;
					$num_bins++;
					push( @bin_labels, int(abs($this_color)) );
				}

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
				
				if ( $this_bin < 0 && $lines_spots[$i][$k] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ ) {
					$this_color = "";
					if ($lines_spots[$i][$k] =~ /MANUAL_COLOR=\"([+-]?\d+(\.?\d*))\"/ || $lines_spots[$i][$k] =~ /TISSUE_TYPE=\"([+-]?\d+(\.?\d*))\"/  ) {
						$this_color = $1;
					}
					if ( $this_color eq "" ) {
						#could not find a color
						$this_bin = 0;
					} elsif ( exists($colors_bins{$this_color}) ) {
						$this_bin = $colors_bins{$this_color};
					} else {
						$colors_bins{$this_color} = $num_bins;
						$this_bin = $num_bins;
						$num_bins++;
					}
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
	
	#finally, for each subdivision, update nspots and write file
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
		
		if ( open(FILE, ">$path/$dataset_file\.$bin_labels[$l]\.xml" ) ) {
			flock(FILE, LOCK_EX);
			for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) {
				if ( $fl == $line_number_spots ) {
					for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
						for ( my $k=0; $k<scalar(@{$lines_spots_divided[$i][$l]}); $k++ ) {
							#print FILE $i . "," . $k . ": " . $lines_spots_divided[$i][$l][$k] . "\n";
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
			print "Error writing to $path/$dataset_file\.$l\.xml!\n";
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


