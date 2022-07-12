#!/usr/bin/perl

# MaMuT print track coordinates in time
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Reconstructs a MaMuT dataset by [track] lineage coordinates into a .tsv table
# it prints each [track] lineage in a row of the output table
# columns of the output table represent each timepoint
# each cell's XYZ coordinates are printed in the appropriate cell of the table for its lineage/track and timepoint
# usages:
# perl MaMuT_dataset_print_track_coordinates_in_time.pl dataset_mamut.xml
# perl MaMuT_dataset_print_track_coordinates_in_time.pl *.xml parallel=8 #processes 8 files at a time with .xml extension

# note this script is single-threaded and slow!

BEGIN {
	@ARGV = map glob, @ARGV;
}

use Cwd qw( cwd );

use constant TRUE => 1;
use constant FALSE => 0;

#multithreaded actions
my $can_use_threads = eval 'use threads; 1';



#this script takes an input MaMuT .xml file, and reconstructs each [track] lineage, printing each in a new row of the output table; each column of the output table represents each timepoint, and the XYZ coordinates of the cell(s) belonging to that track at that timepoint


#==================
#Main subroutine
#==================
sub main_print_track {
	my ( $dataset_file, $path ) = @_;
	my @lines_main; #unchanging lines during reconstruction
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

	#parse dataset.mamut
	if ( $dataset_file !~ /\/$/ && -e $dataset_file && open(FILE, "<$dataset_file" ) ) {
		print "Examining track coordinates in time for $dataset_file.\n";
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
		print "Error opening $dataset_file!\n";
		return;
	}
	#return;
	#strip filename of extension, and assume it is XML
	( undef, $dataset_file ) = split( /\./, reverse($dataset_file), 2 );
	$dataset_file = reverse($dataset_file);

	#iterate over spots at t=0, subdivide them, while storing IDs for tree reconstruction
	my @output_table; #array of array, with each array mirroring track_in_tracks, and storing all coordinates for each timepoint 0...n in [0]...[n] positions of the array
	
	for ( my $l=0; $l<scalar(@track_in_tracks); $l++ ) { #iterate over known tracks
		my @xyz;
		my $xyzprint;
		my $this_spot_id;
		my @lineages = (); #array of array to store 
		
		for ( my $i=0; $i<scalar(@lines_spots); $i++ ) { #then iterate over timepoints
			$current_timepoint_positions = "";
			$output_table[$l][$frames_spots[$i]] = ""; #clear current timepoint for this track
			
			for ( my $k=scalar(@{$lines_spots[$i]})-1; $k>=0; $k-- ) { #look at all spots at this timepoint
				if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /ID=\"(\d+)\"/ ) {
					#capture ID
					$this_spot_id = $1;
					@xyz = ();
					
					#see if this spot belongs in this track
					my $going = 0;
					for ( my $j=scalar(@{$spots_in_tracks[$l]})-1; $j>=0; $j-- ) {	
						if ( $spots_in_tracks[$l][$j] eq $this_spot_id ) {
							$going = 1;
							last;
						}
					}
					next unless ( $going > 0 );
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
						$current_lineage = "0";
					}
					
					#find our lineage
					for ( my $q=scalar(@lineages)-1; $q>=0; $q-- ) {
						if ( $lineages[$q] eq $this_spot_id ) {
							$current_lineage = "$q";
							if ( scalar(@daughters) > 1 ) {
								my $new_lineage_spot = scalar(@lineages);
								$current_lineage .= "->" . $new_lineage_spot . "+" . ($new_lineage_spot+1);
								$lineages[$new_lineage_spot] = $daughters[0];
								$lineages[$new_lineage_spot+1] = $daughters[1];
							} elsif (scalar(@daughters) == 1)  {
								$lineages[$q] = $daughters[0];
							}
							last;
						}
					}
					 
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
					$xyzprint = "{" . $current_lineage . "," . join(',',@xyz) . "}";
					if ( $output_table[$l][$frames_spots[$i]] eq "" ) { #empty so don't try to find the spot
						$output_table[$l][$frames_spots[$i]] = $xyzprint;
					} else { 
						$output_table[$l][$frames_spots[$i]] .= ";" . $xyzprint; #concatenate output
					}

				}
			}
		}
	}

	unless( scalar(@frames_spots) > 1 && scalar(@track_in_tracks) > 0 ) {
		print "Insufficient track reconstruction for " . $dataset_file . "\n";
		return;
	}

	if ( open(FILE, ">$path/$dataset_file\.track_coordinates_in_time\.tsv" ) ) {
		print "Writing track coordinates in time for $dataset_file.\n";
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
		print "Error writing to $path/$dataset_file\.$l\.xml!\n";
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



#==================
# MAIN
#==================
#main_sub {
	my @files;
	for ( my $i=$#ARGV; $i>=0; $i-- ) {
		if ( $ARGV[$i] !~ /\=/ && $ARGV[$i] !~ /\/$/ && -e $ARGV[$i] ) {
			push( @files, splice( @ARGV, $i, 1 ) );
		}
	}

	#handle batch processing
	if ( scalar(@files) > 1 ) {
		my $path = Cwd::cwd();
		my $num_threads = 4; #default 4 threads
		for ( my $i=$#ARGV; $i>=0; $i-- ) {
			if ( $ARGV[$i] =~ /parallel==?(\d+)$/ ) { #proper command line for number of threads is *
				$num_threads = $1;
				splice( @ARGV, $i, 1 );
			}
		}
		unless ( defined($num_threads) && $num_threads =~ /^\d+$/ && $can_use_threads && $num_threads > 1) {
			$num_threads = 1;
		}

		my @tot_items = @files;

		my $l__ = 0;
		my @run_threads_ = (); #running threads
		my @tot_threads_ = (); #total threads, at the end, scalar of this should equal $l__

		while( $l__ < scalar(@tot_items) ) {
			if ( scalar( @run_threads_ ) < $num_threads) {
				$tot_threads_[$l__] = threads->create( {'context' => 'list'}, sub { return main_print_track( $tot_items[$l__], $path ) });
				$l__++;
				sleep 1; #give at least one second between requests
			}
			@run_threads_ = threads->list(threads::running());
			foreach my $oldth__ (@tot_threads_) {
				if ($oldth__->is_running()) {
				} elsif ($oldth__->is_joinable()) {
					$oldth__->join();
				}
			}
			@run_threads_ = threads->list(threads::running());
		}

		#clean up finished threads
		@run_threads_ = threads->list(threads::running());
		while (scalar @run_threads_ > 0) {
			foreach my $oldth__ (@tot_threads_) {
				if ($oldth__->is_joinable()) {
					$oldth__->join();
				}
			}
			$l__ = scalar @run_threads_;
			@run_threads_ = threads->list(threads::running());
			redo if ( scalar @run_threads_ == 0 && $l__ > 0 );
		}

	} elsif ( scalar(@files) == 1 ) {
		my $path = Cwd::cwd();
		main_print_track($path . qq~/~ . $files[0], $path );
	} else {
		print "Improperly formatted command line, no matching dataset files found.\n";
	}
#}

