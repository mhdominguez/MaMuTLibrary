#!/usr/bin/perl

# MaMuT track coordinates pairwise analyze
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Analyzes track pairs exported from MaMuT_dataset_print_track_coordinates_in_time.pl, printing:
#   a table of features such as Begin Coord, Timepoint start, Begin Diff, End Diff
#   a summary table of single-axis track crossing events, and relative positions of track pair begin or end along that axis
# usage: perl MaMuT_track_coordinates_pairwise_analyze.pl dataset_mamut.track_coordinates_in_time.tsv start=15 stop=90 admit_stop=45 cell_max_dist=250 xpos #(analyzes track pairs along x-axis, starting at timepoint 15, stopping at timepoint 90, and without admission of new tracks after timepoint 45, using maximum distance of 250 between the cells as a filter

# this script is very slow!


use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;


use constant {
	FEATURE_TRACK_DURATION => 0,
	FEATURE_TRACK_START => 1,
	FEATURE_TRACK_STOP => 2,
	FEATURE_TRACK_ADMIT_STOP => 3,	
	FEATURE_TRACK_DISPLACEMENT => 4,
	FEATURE_MAX_CELL_DISTANCE => 5,
	FEATURE_MIN_CELL_DISTANCE => 6,
	FEATURE_MAX_DISTANCE_BY_XYZ => 7,
	
	FEATURE_TRACK_MAX => 8,
};

use constant {
	FEATURE_CELL_POS_X => 0,
	FEATURE_CELL_POS_Y => 1,
	FEATURE_CELL_POS_Z => 2,
	
	FEATURE_CELL_DISTANCE => 3,
	
	FILTER_DIRECTION_MAX => 4
};


#==================
#Main subroutine
#==================
sub main {
	my $dataset_file = shift;
	my @params = @_;
	#initialize instructions for this analysis
	my @track_feature = (-1) x FEATURE_TRACK_MAX;
	$track_feature[FEATURE_TRACK_START] = 0;
	my $action = FEATURE_CELL_POS_Y; #default is analyze Y movement of paired tracks
	#my $max_distance = -1;
	
	my @this_param = ();
	for ( my $k=0; $k<scalar(@params); $k++ ) {
		$params[$k] =~ s/==/=/g;
		@this_param = split( /=/, $params[$k] );

		if ( $this_param[0] =~ /start/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_START] = $this_param[1]; 
		} elsif ( $this_param[0] =~ /stop/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			if ( $this_param[0] =~ /new/i || $this_param[0] =~ /ad/i ) {
				$track_feature[FEATURE_TRACK_ADMIT_STOP] = $this_param[1]; 
			} else {
				$track_feature[FEATURE_TRACK_STOP] = $this_param[1]; 
			}
		} elsif ( $this_param[0] =~ /track_displacement/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_DISPLACEMENT] = $this_param[1]; 
		} elsif ( $this_param[0] =~ /dist/i && $this_param[0] =~ /cell/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			if ( $this_param[0] =~ /min/i ) {
				$track_feature[FEATURE_MIN_CELL_DISTANCE] = $this_param[1];
			} elsif ( $this_param[0] =~ /max/i ) {
				$track_feature[FEATURE_MAX_CELL_DISTANCE] = $this_param[1];
			}
		} elsif ( $this_param[0] =~ /track_duration/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_DURATION] = $this_param[1]; 
		#} elsif ( $this_param[0] =~ /max/i && $this_param[0] =~ /dist/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
		#	$max_distance = $this_param[1]; 			
		} elsif ( $this_param[0] =~ /pos/i ) {
			if ( $this_param[0] =~ /x/i ) {
				$action = FEATURE_CELL_POS_X;
			} elsif ( $this_param[0] =~ /y/i ) {
				$action = FEATURE_CELL_POS_Y;
			} elsif ( $this_param[0] =~ /z/i ) {
				$action = FEATURE_CELL_POS_Z;
			}
		} elsif ( $this_param[0] =~ /dist/i && $this_param[1] =~ /\w/ ) {
			$action = FEATURE_CELL_DISTANCE;
			#print "here\n";
		} elsif ( $this_param[0] =~ /dist/i && $this_param[0] =~ /xyz/i ) {
			$track_feature[FEATURE_MAX_DISTANCE_BY_XYZ] = 1;
			#print "here2\n";
		}
	}
	undef @this_param;		
	
	
	#read in coordinate data
	unless ( open(FILE, "<$path/$dataset_file" ) ) {
		print( "Cannot open $symb_dir/$inlist: $!\n" );
		return;
	}
	my @lines_infile;
	my $li = 0;
	my @this_line;
	flock( FILE, LOCK_EX );
	while( <FILE> ) {
		chomp;
		@this_line = split( /\t/, $_ );
		#$lines_infile[$li] = [ ($this_line[0], @this_line[3..$#this_line] ];
		$lines_infile[$li] = [ @this_line ];
		$li++;	
	}
	flock( FILE, LOCK_UN );
	close( FILE );
	
	#set up parameters for analysis now that data is available
	my $going;	
	#print $track_feature[FEATURE_TRACK_DURATION] . "\n";
	$track_feature[FEATURE_TRACK_STOP] = scalar(@{$lines_infile[0]}) if ( $track_feature[FEATURE_TRACK_STOP] < 0 );
	if ( $track_feature[FEATURE_TRACK_DURATION] > 1 || $track_feature[FEATURE_TRACK_DISPLACEMENT] > 1 ) { #ignore tracks below a certain minimum duration
		my @final_tracks = ($lines_infile[0]); #where good tracks go
		for ( my $i=1; $i<scalar(@lines_infile); $i++ ) { #look at each track
			$going = 0;
			for ( my $t=1; $t<scalar(@{$lines_infile[0]}); $t++ ) {
				if ( defined($lines_infile[$i][$t]) && $lines_infile[$i][$t] ne "" && $lines_infile[$i][$t] =~ /\{.*?,(.*?),(.*?),(.*?),.*?\}/) {
					my @xyz_begin = ( $1, $2, $3 );
					for ( my $u=$t+1; $u<=scalar(@{$lines_infile[0]}); $u++ ) {
						unless( defined($lines_infile[$i][$u]) && $lines_infile[$i][$u] =~ /\d/ ) { #found end at $u-1
							if ( $lines_infile[$i][$u-1] =~ /\{.*?,(.*?),(.*?),(.*?),.*?\}/ ) {
								my @xyz_end = ( $1, $2, $3 );
								$going = 1; #default is pass
								if ( $track_feature[FEATURE_TRACK_DURATION] > 1 ) {
									if ( $track_feature[FEATURE_TRACK_DURATION] > $u - $t - 1 ) { #doesn't make the cut
										$going = -1;
										last;
									} else {
										#$going = 1;
									}
								}
								if ( $track_feature[FEATURE_TRACK_DURATION] > 1 ) {
									if ( $track_feature[FEATURE_TRACK_DISPLACEMENT] > sqrt((($xyz_end[0]-$xyz_begin[0])**2) + (($xyz_end[1]-$xyz_begin[1])**2) + (($xyz_end[2]-$xyz_begin[2])**2)) ) { #doesn't make the cut
										$going = -1;
										last;
									} else {
										#$going = 1;
									}
								}
								last;
							}
						}
					}
					last;
				}
			}
			if ( $going > 0 ) {
				push (@final_tracks, $lines_infile[$i] );		
			}
		}
		
		@lines_infile = @final_tracks; #replace original lines with new filtered lines
	}
	
	#print scalar(@lines_infile) . "\n";
	my @output_table;
	
	#setup max timepoint for new track admission
	my $track_pickup_stop = $track_feature[FEATURE_TRACK_STOP];
	$track_pickup_stop = $track_feature[FEATURE_TRACK_ADMIT_STOP] if ( $track_feature[FEATURE_TRACK_ADMIT_STOP] > $track_feature[FEATURE_TRACK_START] && $track_feature[FEATURE_TRACK_ADMIT_STOP] < $track_feature[FEATURE_TRACK_STOP] ); 
	
	for ( my $i=1; $i<scalar(@lines_infile); $i++ ) { #look at each track
		$going = 0;
		my @xyz_begin, @xyz_end;
		my @xyz_other_begin, @xyz_other_end;
		
		#pick up start of tracks
		for ( my $t=$track_feature[FEATURE_TRACK_START]+1; $t<$track_pickup_stop; $t++ ) {
			if ( defined($lines_infile[$i][$t]) && $lines_infile[$i][$t] ne "" && $lines_infile[$i][$t] =~ /\{.*?,(.*?),(.*?),(.*?),.*?\}/) {
				$going = $t;
				@xyz_begin = ( $1, $2, $3 );
				last;
			}	
		}
		
		next unless ( $going > 0 );
		#print "HERE\n";
		#now, iterate over non-self tracks from going to end of this track to create triangle matrix so no work is duplicated
		for ( my $j=$i+1; $j<scalar(@lines_infile); $j++ ) {
			my @output_array;
			for ( my $t=$going; $t<$track_feature[FEATURE_TRACK_STOP]; $t++ ) {
				if ( defined($lines_infile[$i][$t]) && $lines_infile[$i][$t] ne "" && $lines_infile[$i][$t] =~ /\{.*?,(.*?),(.*?),(.*?),.*?\}/ ) {
					if ( scalar(@output_array) == 0 ) { #original cell is defined at this timepoint
						@xyz_begin = ( $1, $2, $3 );
					} else {
						@xyz_end = ( $1, $2, $3 );
					}
					if ( defined($lines_infile[$j][$t]) && $lines_infile[$j][$t] ne "" && $lines_infile[$j][$t] =~ /\{.*?,(.*?),(.*?),(.*?),.*?\}/ ) {
						#my @capture = ( $1, $2, $3 );
						#if ( $capture[0] == 0 ) {
						#	print "zero!!";
						#}
						if ( scalar(@output_array) == 0 ) { #new track is defined at this timepoint
							@xyz_other_begin = ( $1, $2, $3 );
							my $cell_distance_min = 0;
							my $cell_distance_max = 0;
								
							if ( $action == FEATURE_CELL_DISTANCE ) {
								$cell_distance_min = sqrt((($xyz_other_begin[0]-$xyz_begin[0])**2) + (($xyz_other_begin[1]-$xyz_begin[1])**2) + (($xyz_other_begin[2]-$xyz_begin[2])**2));
								$cell_distance_max = $cell_distance_min;
							} elsif ( $track_feature[FEATURE_MAX_DISTANCE_BY_XYZ] > 0 ) {
								$cell_distance_min = abs($xyz_other_begin[$action] - $xyz_begin[$action]);
								$cell_distance_max = sqrt((($xyz_other_begin[0]-$xyz_begin[0])**2) + (($xyz_other_begin[1]-$xyz_begin[1])**2) + (($xyz_other_begin[2]-$xyz_begin[2])**2));
							} else {
								$cell_distance_min = abs($xyz_other_begin[$action] - $xyz_begin[$action]);
								$cell_distance_max = $cell_distance_min;
							}
							#if ( $xyz_other_begin[0] == 0 || $xyz_begin[0] == 0 ) {
								#print "zero!";
							#} 
							if ( ($track_feature[FEATURE_MAX_CELL_DISTANCE] <= 0 || $cell_distance_max < $track_feature[FEATURE_MAX_CELL_DISTANCE]) && ($track_feature[FEATURE_MIN_CELL_DISTANCE] <= 0 || $cell_distance_min >= $track_feature[FEATURE_MIN_CELL_DISTANCE]) ) {
								if ( $action < FEATURE_CELL_DISTANCE ) {
									@output_array = ( $lines_infile[$i][0] . "/" . $lines_infile[$j][0] , $xyz_begin[$action] . "/" . $xyz_other_begin[$action], $t, $xyz_other_begin[$action] - $xyz_begin[$action] ); #first array element is the start timepoint where both cells are visible
								} elsif  ( $action == FEATURE_CELL_DISTANCE ) {
									@output_array = ( $lines_infile[$i][0] . "/" . $lines_infile[$j][0] , "{" . join(',',@xyz_begin) . "} / {" . join(',',@xyz_other_begin) . "}", $t, $cell_distance_max ); #first array element is the start timepoint where both cells are visible
								}
							}
						} else {
							@xyz_other_end = ( $1, $2, $3 );
							if ( $action < FEATURE_CELL_DISTANCE ) {
								push( @output_array, $xyz_other_end[$action] - $xyz_end[$action] );
							} elsif  ( $action == FEATURE_CELL_DISTANCE ) {
								push( @output_array, sqrt((($xyz_other_end[0]-$xyz_end[0])**2) + (($xyz_other_end[1]-$xyz_end[1])**2) + (($xyz_other_end[2]-$xyz_end[2])**2)) );
							}
						}
					}
				}
			}
			push( @output_table, \@output_array ) if ( scalar(@output_array) > 2 );
		}
	}
	my $max_columns = 0;
	map { $max_columns = scalar(@{$_}) if ( scalar(@{$_}) > $max_columns ); } @output_table;
	#print $max_columns . "\n";
	my @fisher_table = (0,0,0,0);
	my @begin_diffs;
	my @end_diffs;
	my @num_crosses;
	my $current_diff;
	if ( open(FILE, ">$path/$dataset_file\.pairwise_track_analysis\.tsv" ) ) {
		flock(FILE, LOCK_EX);
		print FILE "Tracks\tBegin Coord\tTimepoint start\tBegin Diff\tEnd Diff\tBE Ratio\tIsCross";
		for ( my $i=0; $i<$max_columns; $i++ ) {
			print FILE "\tFrame " . $i;
		}
		print FILE "\n";
		for ( my $l=0; $l<scalar(@output_table); $l++ ) {
		#	if ( $max_distance < 0 || abs($output_table[$l][3]) < $max_distance ) {
				#print FILE $output_table[$l][0] . "\t";
				my $end_diff = $output_table[$l][scalar(@{$output_table[$l]})-1];
				push( @end_diffs, abs($end_diff));
				my $begin_diff = $output_table[$l][3];
				push( @begin_diffs, abs($begin_diff));
				print FILE join( "\t", (@{$output_table[$l]})[0..2] ) . "\t" . abs($begin_diff);
				#print FILE "\t" . ($output_table[$l][scalar(@{$output_table[$l]})-1] - $output_table[$l][3]);
				if ( $begin_diff >= 0 ) {
					print FILE "\t" . $end_diff;
				} else {
					print FILE "\t" . (0-$end_diff);
				}
				
				if ( $end_diff == 0 ) {
					print FILE "\t";
				} else {
					my $ratio = ($begin_diff/$end_diff);
					print FILE "\t" . ($ratio>1000 ? 1000 : $ratio);
				}
				
				#if ( ($begin_diff > 0 && $end_diff <= 0) || ($begin_diff < 0 && $end_diff >= 0) ) {
				#	print FILE "\t1";
				#} else {
				#	print FILE "\t0";
				#}
				if ( $begin_diff >= 0 ) {
					if ( $end_diff >= 0 ) {
						# + +
						print FILE "\t0";
						$fisher_table[0]++;
					} else {
						# + -
						print FILE "\t1";
						$fisher_table[1]++;
					}
				} else {
					if ( $end_diff >= 0 ) {
						# - +
						print FILE "\t1";
						$fisher_table[2]++;
					} else {
						# - -
						print FILE "\t0";
						$fisher_table[3]++;
					}
				}
				#	print FILE "\t1";
				#} else {
				#	print FILE "\t0";
				#}
				$num_crosses[$l] = 0;
				$current_diff = $begin_diff;
				for ( my $k=4; $k<scalar(@{$output_table[$l]}); $k++ ) {
					$end_diff = $output_table[$l][$k];
					if ( ($current_diff > 0 && $end_diff <= 0) || ($current_diff < 0 && $end_diff >= 0) ) {
						$num_crosses[$l]++;
						$current_diff = $end_diff;
					}
					if ( $begin_diff >= 0 ) {
						print FILE "\t" . $end_diff;
					} else {
						print FILE "\t" . (0-$end_diff);
					}
				}
				print FILE "\n";
			#}
		}
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error writing to $path/$dataset_file\.pairwise_track_analysis\.tsv!\n";
		return;
	}
	
	unless ( -e "$path/pairwise_track_analysis_fisher_table\.tsv" ) {
		if ( open(FILE, ">$path/pairwise_track_analysis_fisher_table\.tsv" ) ) {
			flock(FILE, LOCK_EX);
			print FILE "File\tNoCross\tYesCross\t+,+\t+,-\t-,+\t-,-\tAvgBeginDiff\tAvgEndDiff\tAvgNumCrosses\tAvgNumCrossesIfCross\n";
			flock(FILE, LOCK_UN);
			close(FILE);
		} else {
			print "Error writing to $path/pairwise_track_analysis_fisher_table\.tsv!\n";
			return;
		}
	}
	my @num_crosses_real;
	my @num_crosses_real_greater_than_0;
	map { push( @num_crosses_real, $_ ) if (defined($_) && $_ =~ /\d/ ); } @num_crosses;
	map { push( @num_crosses_real_greater_than_0, $_ ) if ($_>0); } @num_crosses_real;
	#print join ( "\n", @num_crosses_real ) . "\n";
	if ( open(FILE, ">>$path/pairwise_track_analysis_fisher_table\.tsv" ) ) {
		flock(FILE, LOCK_EX);
		#print FILE "\t+\t-\n";
		#print FILE "+\t" . $fisher_table[0] . "\t" . $fisher_table[1] . "\n";
		print FILE $dataset_file . "\t" . ($fisher_table[0]+$fisher_table[3]) . "\t" . ($fisher_table[1]+$fisher_table[2]) . "\t" . $fisher_table[0] . "\t" . $fisher_table[1] . "\t" . $fisher_table[2] . "\t" . $fisher_table[3] . "\t";
		print FILE (set_stats(@begin_diffs))[0] if (scalar(@begin_diffs) > 1 );
		print FILE "\t";
		print FILE (set_stats(@end_diffs))[0] if (scalar(@end_diffs) > 1 );
		print FILE "\t";
		print FILE ((set_stats(@num_crosses_real))[0]) if (scalar(@num_crosses_real) > 1 && scalar(@output_table) > 1);
		print FILE "\t";
		print FILE ((set_stats(@num_crosses_real_greater_than_0))[0]) if (scalar(@num_crosses_real_greater_than_0) > 1 && scalar(@output_table) > 1);
		print FILE "\n";		
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error writing to $path/pairwise_track_analysis_fisher_table\.tsv!\n";
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


