#!/usr/bin/perl

# MaMuT split dataset by timepoints
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Splits a MaMuT dataset into two datasets, one containing spots and tracks bounded by start/stop, the other with spots that are out-of-bounds
# usage: perl MMaMuT_dataset_split_spot_timepoint.pl dataset_mamut.xml start=0 stop=40
#	where 0 and 40 create the timepoint bounds for the resulting dataset

# note this script is relatively fast, since it just splits the spots and doesn't really reconstruct lineages -- user may need to open and re-save the dataset in order for everything to work properly after splitting by timepoint


use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;



use constant {
	FEATURE_START => 0,
	FEATURE_STOP => 1,
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
	my @lines_tracks_spots;
	#my @track_tracks;
	#my @spots_in_tracks;
	#my @track_in_tracks;
	my $line_number_filters = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_filters;
	my @track_in_filters;
	
	my $num_bins = 2;
	#my @bins_tracks;
	my $overlapping_bins = FALSE;
	my @feature;
	$feature[FEATURE_START] = -1;
	$feature[FEATURE_STOP] = -1;

	for ( my $k=0; $k<scalar(@params); $k++ ) {
		$params[$k] =~ s/==/=/g;
		my @this_param = split( /=/, $params[$k] );

		if ( $this_param[0] =~ /start/i && defined($this_param[1]) && $this_param[1] =~ /\d+/ ) {
			$feature[FEATURE_START] = $this_param[1];
		} elsif ( $this_param[0] =~ /stop/i && defined($this_param[1]) && $this_param[1] =~ /\d+/ ) {
			$feature[FEATURE_STOP] = $this_param[1];
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
					if ( $_ =~ /(^|^\s+)<SpotsInFrame.*frame=\"(\d+)\"/i && $_ !~ /\/>\s*$/ ) {
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
						my @this_block_spots = ([ "N\/A" ]);
						#print "inside here: $this_track_id\n";
						while( <FILE> ) {
							chomp;
							push( @this_block, $_ );
							push( @this_block_spots, ["N\/A"] );
							if ( $_ =~ /(^|^\s+)<\/Track/ ) {
								last;
							} elsif ($_ =~ /(^|^\s+)<Edge/) { #store all spots associated with this track for reconstruction
								if ( $_ =~ /SPOT_SOURCE_ID=\"(\d+)\".*?SPOT_TARGET_ID=\"(\d+)\"/ ) {
									#print "here\n";
									#push( @{$spots_in_tracks[scalar(@lines_tracks)]}, $1, $2 );
									$this_block_spots[$#this_block_spots] = [$1,$2];
								}
							}
						}
						#print "scalar " . scalar(@{$spots_in_tracks[scalar(@lines_tracks)]}) . "\n";
						#@{$spots_in_tracks[scalar(@lines_tracks)]} = uniq(@{$spots_in_tracks[scalar(@lines_tracks)]});
						#$track_in_tracks[scalar(@lines_tracks)] = $this_track_id;
						#$bins_tracks[scalar(@lines_tracks)] = -1;
						push( @lines_tracks, \@this_block );
						push( @lines_tracks_spots, \@this_block_spots );
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
	
	if ( $feature[FEATURE_STOP] < 0 ) {
		$feature[FEATURE_STOP] = max(@frames_spots);
	}
	if ( $feature[FEATURE_START] < 0 ) {
		$feature[FEATURE_START] = min(@frames_spots);
	}	
	print "Going to split dataset by timepoint, $feature[FEATURE_START] to $feature[FEATURE_STOP]...\n";
	#return;
	#strip filename of extension, and assume it is XML
	( undef, $dataset_file ) = split( /\./, reverse($dataset_file), 2 );
	$dataset_file = reverse($dataset_file);

	#iterate over spots at t=0, subdivide them, while storing IDs for tree reconstruction
	my @lines_spots_divided; #final storage for subdivisions, 3D matrix [timepoint][subdivision][spot], compared with 2D matrix of @lines_spots [timepoint][spot]
	my @bins_spots_ids;
	my @tracking_spots_divided; # storage for spot IDs undergoing lineage reconstruction by the splitting algorithm, 2D matrix [subdivision][spotID], each spotID will be removed once the daughter spot is found and added to @lines_spots_divided
	my $this_spot_id;
	my $this_bin;
	my $captured;

	#now, iterate over all spots in time, reconstructing lineages and sorting tracks to the correct bins
	my $this_track; my $this_spot_in_track;
	my $going;
	#my $count = 0;
	for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
		#print "t=$i, total spots handled: " . scalar(@{$lines_spots[$i]}) . "\n";
		for ( my $k=scalar(@{$lines_spots[$i]})-1; $k>=0; $k-- ) {
			if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /ID=\"(\d+)\"/ ) {
				#capture ID
				$this_spot_id = $1;
				$this_bin = -1;
				
				if ( $frames_spots[$i]<$feature[FEATURE_START] ) {
					$this_bin = 0;
				} else {
					$this_bin = 1;
				}
				
				if ( $this_bin > 0 ) { #not already rejected
					if ( $frames_spots[$i]>$feature[FEATURE_STOP] ) {
						$this_bin = 0;
					}
				}
				
				next unless ( $this_bin >= 0 );
				
				#$count++;
				
				push( @{$lines_spots_divided[$i][$this_bin]}, $lines_spots[$i][$k] );
				push( @{$bins_spots_ids[$this_bin]}, $this_spot_id );

			} else {
				#control blocks get passed to both divisions
				for (my $l=0; $l<$num_bins; $l++ ) {
					push( @{$lines_spots_divided[$i][$l]}, $lines_spots[$i][$k] );
				}
			}
		}
		undef @{$lines_spots[$i]}; #freeing space
	}
	
	#print "count $count\n";
	#return;
	#
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
		my @print_track_ids;
		
		if ( open(FILE, ">$path/$dataset_file\.$l\.xml" ) ) {
			flock(FILE, LOCK_EX);
			for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) {
				#print "$lines_main[$fl]\n";
				if ( $fl == $line_number_spots ) {
						#print $l . "/" . $fl . ": " . scalar(@lines_spots_divided) . "\n";
						for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
							#print "$i: " . scalar(@{$lines_spots_divided[$i][$l]}) . "\n";
							for ( my $k=scalar(@{$lines_spots_divided[$i][$l]})-1; $k>=0; $k-- ) {
								print FILE $lines_spots_divided[$i][$l][$k] . "\n";
							}
						}
				} elsif ( $fl == $line_number_tracks ) {			
					for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						my $to_print = "";
						my $do_print = 0;
						my $this_track_id = "";
						my $curr_z = 0;
						for ( my $q=0; $q<scalar(@{$lines_tracks[$i]}); $q++ ) {
							
							print FILE $lines_tracks[$i][$q] . "\n";
							next;
							
#							#below is incredibly slow and inefficient, even if the resulting file sizes are smaller
#							
#							if ($lines_tracks_spots[$i][$q][0] eq "N\/A" ) { #generally going to be a control block
#								#print FILE join( "\n", @{$lines_tracks[$i]} ) . "\n";
#								$to_print .= $lines_tracks[$i][$q] . "\n";
#								
#								if ( $lines_tracks[$i][$q] =~ /TRACK_ID=\"(\d+)\"/i ) {
#									$this_track_id = $1;
#								}
#								
#							} else {
#								#try to match this line into a bin based on spot bins_spots_ids
#								$going = 0;
#								my $deref0 = $lines_tracks_spots[$i][$q][0];
#								my $deref1 = $lines_tracks_spots[$i][$q][1];
#								#map { $going++ if ( $_ eq $deref0 || $_ eq $deref1 ) } @{$bins_spots_ids[$l]};
#
#								for ( my $z=$curr_z; $z<scalar(@{$bins_spots_ids[$l]}); $z++ ) {
#									if ($deref0 eq $bins_spots_ids[$l][$z]) {
#										$going++;
#										$curr_z = $z;
#										last;
#									}
#								}
#								if ( $going == 0 ) {
#									for ( my $z=0; $z<$curr_z; $z++ ) {
#										if ($deref0 eq $bins_spots_ids[$l][$z]) {
#											$going++;
#											$curr_z = $z;
#											last;
#										}
#									}
#								}
#								
#								$curr_z -= 20;
#								$curr_z = 0 if ( $curr_z < 0 );
#								
#								my $going1 = 0;
#								for ( my $z=$curr_z; $z<scalar(@{$bins_spots_ids[$l]}); $z++ ) {
#									if ($deref1 eq $bins_spots_ids[$l][$z]) {
#										$going1++;
#										$curr_z = $z;
#										last;
#									}
#								}
#								if ( $going1 == 0 ) {
#									for ( my $z=0; $z<$curr_z; $z++ ) {
#										if ($deref1 eq $bins_spots_ids[$l][$z]) {
#											$going1++;
#											$curr_z = $z;
#											last;
#										}
#									}
#								}
#								
#								#for ( my $z=0; $z < scalar(@{$bins_spots_ids[$l]}); $z++ ) {
#								#	if ($deref1 eq $bins_spots_ids[$l][$z]) {
#								#		$going++;
#								#		#$curr_z = $z;
#								#		last;
#								#	}
#								#}
#								#for ( my $z=0; $z < scalar(@{$bins_spots_ids[$l]}); $z++ ) {
#								#	if ($deref0 eq $bins_spots_ids[$l][$z]) {
#								#		$going++;
#								#		#$curr_z = $z;
#								#		last;
#								#	}
#								#}
#								
#								if ( $going + $going1 > 1 ) { #both spots matched!
#									$to_print .= $lines_tracks[$i][$q] . "\n";	
#									$do_print = 1;
#								}
#							
#							}
						}
#						
#						if ( $do_print ) {
#							print FILE $to_print;
#							push( @print_track_ids, $this_track_id );
#						}

					}
				} elsif ( $fl == $line_number_filters ) {
					#for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						#if ( $bins_tracks[$i] == $l ) {
							for ( my $k=0; $k<scalar(@lines_filters); $k++ ) {
#								if ( $lines_filters[$k] =~ /TRACK_ID=\"(\d+)\"/i ) {
#									my $this_track_id = $1;
#									$going = 0;
#									map {$going++ if ($_ eq $this_track_id) } @print_track_ids;
#									
#									print FILE $lines_filters[$k] . "\n" if ( $going > 0 );
#									#last;
#								} else {
									print FILE $lines_filters[$k] . "\n";
#								}
							}
						#}
					#}	
				}
				print FILE $lines_main[$fl] . "\n";
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


