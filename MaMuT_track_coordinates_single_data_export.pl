#!/usr/bin/perl

# MaMuT track coordinates single data export
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Summarizes all tracks exported from MaMuT_dataset_print_track_coordinates_in_time.pl, printing a table of features such as begin/end XYZT, peak displacement, average density, etc
# usages: 
#  perl MaMuT_track_coordinates_single_data_export.pl dataset_mamut.track_coordinates_in_time.tsv
#  perl MaMuT_track_coordinates_single_data_export.pl dataset_mamut.track_coordinates_in_time.tsv start=0 stop=120 velocity_window=20 timepoints_per_hour=10
#    where 30-120 are the timepoint boundaries for consideration, 20 is the velocity moving average window (motility period), 10 timepoints per hour are captured
#  perl MaMuT_track_coordinates_single_data_export.pl dataset_mamut.track_coordinates_in_time.tsv density_radius=5
#    where 5 times cell's radius is used for density calculations 


BEGIN {
	@ARGV = map glob, @ARGV;
}

use Cwd qw( cwd );
use File::Spec;

use constant TRUE => 1;
use constant FALSE => 0;

#multithreaded actions
my $can_use_threads = eval 'use threads; 1';

use constant {
	FEATURE_TRACK_DURATION => 0,
	FEATURE_TRACK_START => 1,
	FEATURE_TRACK_STOP => 2,
	FEATURE_TRACK_DISPLACEMENT => 3,

	FEATURE_TRACK_MAX => 4,
};

my $density_radius = 12;
my $velocity_window = 5;
my $timeframe_per_hour = 10;

#==================
#Main subroutine
#==================
sub main_single_data_export {
	my $dataset_file = shift;
	my $path = shift;
	my @params = @_;
	#initialize instructions for this analysis
	my @track_feature = (-1) x FEATURE_TRACK_MAX;
	$track_feature[FEATURE_TRACK_START] = 0;
	
	my @this_param = ();
	for ( my $k=0; $k<scalar(@params); $k++ ) {
		$params[$k] =~ s/==/=/g;
		@this_param = split( /=/, $params[$k] );

		if ( $this_param[0] =~ /start/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_START] = $this_param[1]; 
		} elsif ( $this_param[0] =~ /stop/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_STOP] = $this_param[1]; 
		} elsif ( $this_param[0] =~ /track_displacement_max/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_DISPLACEMENT] = $this_param[1]; 
		} elsif ( $this_param[0] =~ /track_duration_min/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_DURATION] = $this_param[1]; 
		} elsif ( $this_param[0] =~ /den/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$density_radius = $this_param[1]; 
		} elsif ( $this_param[0] =~ /vel/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$velocity_window = $this_param[1]; 
		} elsif ( ( $this_param[0] =~ /frames_per_hour/i || $this_param[0] =~ /points_per_hour/i ) && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$timeframe_per_hour = $this_param[1]; 
		}
	}
	undef @this_param;			
	
	unless (defined($dataset_file) && $dataset_file =~ /\w/ ) {
		print "Track coordinate input file not provided!\n";
		return;
	}
	
	#read in coordinate data
	unless ( open(FILE, "<" . File::Spec->catfile( $path, $dataset_file ) ) ) {
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
	
	$track_feature[FEATURE_TRACK_STOP] = scalar(@{$lines_infile[0]}) if ( $track_feature[FEATURE_TRACK_STOP] < 0 );
	if ( $track_feature[FEATURE_TRACK_DISPLACEMENT] > 1 ) { #ignore tracks below a certain minimum duration
		my @final_tracks = ($lines_infile[0]); #where good tracks go
		for ( my $i=1; $i<scalar(@lines_infile); $i++ ) { #look at each track
			$going = 0;
			for ( my $t=1; $t<scalar(@{$lines_infile[0]}); $t++ ) {
				if ( defined($lines_infile[$i][$t]) && $lines_infile[$i][$t] ne "" && $lines_infile[$i][$t] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/) {
					my @xyz_begin = ( $2, $3, $4, $5 );
					for ( my $u=$t+1; $u<=scalar(@{$lines_infile[0]}); $u++ ) {
						unless( defined($lines_infile[$i][$u]) && $lines_infile[$i][$u] =~ /\d/ ) { #found end at $u-1
							if ( $lines_infile[$i][$u-1] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/ ) {
								my @xyz_end = ( $2, $3, $4, $5 );
								$going = 1; #default is pass
									if ( $track_feature[FEATURE_TRACK_DISPLACEMENT] > sqrt((($xyz_end[0]-$xyz_begin[0])**2) + (($xyz_end[1]-$xyz_begin[1])**2) + (($xyz_end[2]-$xyz_begin[2])**2)) ) { #doesn't make the cut
										$going = -1;
										last;
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
	
	#for densities, a prior fill table of cell positions for each timepoint
	my @xyz_all;
	for ( my $t=$track_feature[FEATURE_TRACK_START]+1; $t<=$track_feature[FEATURE_TRACK_STOP]; $t++ ) {
		for ( my $i=1; $i<scalar(@lines_infile); $i++ ) { #look at each track
			if ( defined($lines_infile[$i][$t]) && $lines_infile[$i][$t] ne "" && $lines_infile[$i][$t] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/) {
							my @multi_cell = ($lines_infile[$i][$t]);
							if ( $lines_infile[$i][$t] =~ /;/ ) {
								@multi_cell = split( /;/, $lines_infile[$i][$t] );
							}
								
							#search for all cells and place them on the array
							for ( my $q=0; $q<scalar(@multi_cell);$q++ ) { #already did first one above
								if (  $multi_cell[$q] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/ ) {
									push( @{$xyz_all[$t-1]}, [$2,$3,$4] );
								}
							}
			}
		}
	}
	
	#print scalar(@lines_infile) . "\n"; #" $track_feature[FEATURE_TRACK_STOP]\n";
	print "Summarizing $dataset_file, (".scalar(@lines_infile)." lines)" . "\n" . "  from " . $track_feature[FEATURE_TRACK_START] . " to " . $track_feature[FEATURE_TRACK_STOP] . " (" . $timeframe_per_hour . " timepoints per hour) " . "\n" . "  density radius: " . $density_radius . "x radii units, velocity moving average window: " . $velocity_window . " timeframes\n";
	my @output_table;
	my $split;
	my @velocities;
	my @densities;
	for ( my $i=1; $i<scalar(@lines_infile); $i++ ) { #look at each track
		@velocities = ();
		@densities = ();
		$going = 0;
		my @xyzt_begin, @xyzt_multi_end, @xyzt_end, @lineage_splits;
		my @density_numbers;
		my $peak_track_displacement = 0;
		#print "track $i\n";
	
		#pick up start of tracks
		for ( my $t=$track_feature[FEATURE_TRACK_START]+1; $t<$track_feature[FEATURE_TRACK_STOP]-$velocity_window-1; $t++ ) {
			if ( defined($lines_infile[$i][$t]) && $lines_infile[$i][$t] ne "" && $lines_infile[$i][$t] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/) {
				$going = $t;
				$split = $1;
				@xyzt_begin = ( $2, $3, $4, $5, $t-1 );
				@xyzt_end = @xyzt_begin;
				
				if ( $split =~ /(\d+)->(\d+)\+(\d+)/ ) {
					push( @lineage_splits, [($1,$2)],[($1,$3)] );
					$split = $1;
				}
				
				last;
			}	
		}
		
		next unless ( $going > 0 );
		#print "HERE\n";
		#now, iterate over self tracks from going to end of this track to create triangle matrix so no work is duplicated
		my @output_array;
		for ( my $t=$going+1; $t<$track_feature[FEATURE_TRACK_STOP]; $t++ ) {
			if ( defined($lines_infile[$i][$t]) && $lines_infile[$i][$t] ne "" && $lines_infile[$i][$t] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/ ) {
				#original cell's track is defined at this timepoint
				
				#see if more than one cell belongs to this track at this timepoint
				@xyzt_end = ( $2, $3, $4, $5, $t-1 );	 #only one cell case
				$split = $1;
				if ( $split =~ /(\d+)->(\d+)\+(\d+)/ ) {
					push( @lineage_splits, [($1,$2)],[($1,$3)] );
					$split = $1;
				}
				
				if ( $lines_infile[$i][$t] =~ /;/ ) {
					my @multi_cell = split( /;/, $lines_infile[$i][$t] );
					my $max_displacement = 0, $this_displacement = 0;
					
					#search for max displacement
					for ( my $q=0; $q<scalar(@multi_cell);$q++ ) {
						if (  $multi_cell[$q] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/ ) {
							my @xyzt_this = ( $2, $3, $4, $5, $t-1 );
							$split = $1;
							if ( $split =~ /(\d+)->(\d+)\+(\d+)/ ) {
								push( @lineage_splits, [($1,$2)],[($1,$3)] );
								$split = $1;
							}							
							$this_displacement = sqrt((($xyzt_this[0]-$xyzt_begin[0])**2) + (($xyzt_this[1]-$xyzt_begin[1])**2) + (($xyzt_this[2]-$xyzt_begin[2])**2));
							
							if ( $this_displacement > $max_displacement ) {
								@xyzt_end = @xyzt_this;
								$max_displacement = $this_displacement;
								
								$peak_track_displacement = $this_displacement if ( $this_displacement > $peak_track_displacement );
							}
							push(@xyzt_multi_end, [$split,@xyzt_this]);
						}
					}
				} else {
					my $this_displacement = sqrt((($xyzt_end[0]-$xyzt_begin[0])**2) + (($xyzt_end[1]-$xyzt_begin[1])**2) + (($xyzt_end[2]-$xyzt_begin[2])**2));
					$peak_track_displacement = $this_displacement if ( $this_displacement > $peak_track_displacement );
					@xyzt_multi_end = ([$split,@xyzt_end]);
				}
				
				#check velocity and grab number for our average if needed
				for ( my $q=0; $q<scalar(@xyzt_multi_end); $q++ ) {
					if ( $xyzt_multi_end[$q][5] > $velocity_window ) {
						#okay, find cell's ancestor, record displacement
						my $tt = 1+$xyzt_multi_end[$q][5] - $velocity_window; #go back in time
						
						my @lineages_to_find = ( $xyzt_multi_end[$q][0] ); 
						#my $updating = 1;
						#my $yy = 0;
						#while ( $updating > 0 ) { #keep adding lineages until we don't add any more
						#	$updating = 0;
							for ( $y = 0; $y<scalar(@lineages_to_find); $y++ ) {
								for ( my $z=0; $z<scalar(@lineage_splits); $z++ ) {
									if ( $lineages_to_find[$y] eq $lineage_splits[$z][1] ) {
										push( @lineages_to_find, $lineage_splits[$z][0] );
										#$yy = $y+1;
										#$updating = 1;
										last;
									}
								}
							}
						#}
						#print "lineages to find: " . scalar(@lineages_to_find).  "\n";
						if ( defined($lines_infile[$i][$tt]) && $lines_infile[$i][$tt] ne "" && $lines_infile[$i][$tt] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/) {
							my @xyzt_multi_multi_end;
							my @multi_cell = ($lines_infile[$i][$tt]);
							if ( $lines_infile[$i][$tt] =~ /;/ ) {
								@multi_cell = split( /;/, $lines_infile[$i][$tt] );
							}
								
							#search for all cells and place them on the array
							for ( my $q=0; $q<scalar(@multi_cell);$q++ ) { #already did first one above
								if (  $multi_cell[$q] =~ /\{(.*?),(.*?),(.*?),(.*?),(.*?)\}/ ) {
									my @xyzt_multi_this = ( $2, $3, $4, $5, $tt-1 );	 #only one cell case
									$split = $1;
									if ( $split =~ /(\d+)->(\d+)\+(\d+)/ ) {
										#push( @lineage_splits, [($1,$2)],[($1,$3)] );
										$split = $1;
									}
									#unshift( @xyzt_multi_this, $split );
									#push( @xyzt_multi_multi_end, \@xyzt_multi_this );
									
									my $found = "";
									map { $found = $_ if ( $split eq $_ ) } @lineages_to_find;
									
									if ( $found ne "" ) {
										#print "found ancestor\n";
										push( @velocities, sqrt((($xyzt_multi_end[$q][1]-$xyzt_multi_this[0])**2) + (($xyzt_multi_end[$q][2]-$xyzt_multi_this[1])**2) + (($xyzt_multi_end[$q][3]-$xyzt_multi_this[2])**2)) );
										last;
									} else {
										#print "not pushing $split: " . join(',',@lineages_to_find) . "\n";
									}
									
								}
							}
						}
						
					}
						#do density
						my $density_count = 0;
						for ( my $z=0; $z<scalar(@{$xyz_all[$t-1]}); $z++ ) {
							if ( sqrt((($xyzt_multi_end[$q][1]-$xyz_all[$t-1][$z][0])**2) + (($xyzt_multi_end[$q][2]-$xyz_all[$t-1][$z][1])**2) + (($xyzt_multi_end[$q][3]-$xyz_all[$t-1][$z][2])**2)) < $xyzt_multi_end[$q][4] * $density_radius ) {
								$density_count++;# print "hit\n";
							}
						}
						push( @densities, $density_count-1 ); #don't count self
				}

				$going = $t;
			} else {
				#original cell's track not defined at this timepoint, so it is the real end
				last;
			}
		}
		
		#average velocities here
		if ( $going > $xyzt_begin[4] + $track_feature[FEATURE_TRACK_DURATION] ) { #here, track must last at least track_duration
			push( @output_table, [($i-1,@xyzt_begin[0..4],@xyzt_end[0..4],$peak_track_displacement,$timeframe_per_hour*(set_stats(@velocities))[0]/$velocity_window),(set_stats(@densities))[0]] );
			#print "track $i, average vel: " . $output_table[$#output_table][13] . ", " .scalar(@densities)."\n";
		}
	}

	unless( scalar(@output_table) > 1 && scalar(@{$output_table[0]}) > 13 ) {
		print "Not enough track data for " . $dataset_file . "\n";
		return;
	}

	if ( open(FILE, ">" . File::Spec->catfile( $path, $dataset_file ) . "\.single_track_data\.tsv" ) ) {
		flock(FILE, LOCK_EX);
		print FILE "Tracks\tBegin X\tBegin Y\tBegin Z\tBegin R\tBegin T\tEnd X\tEnd Y\tEnd Z\tBegin R\tEnd T\tPeak Displacement\tAvg Sliding Velocity (micron/hr)\tAvg Density (cells per $density_radius radii)\n";
		map { print FILE join("\t",@{$_}) . "\n"; } @output_table;
		flock(FILE, LOCK_UN);
		close(FILE);

		return "$dataset_file\.single_track_data\.tsv";
	} else {
		print "Error writing to " . File::Spec->catfile( $path, $dataset_file ) . "\.single_track_data\.tsv!\n";
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

	return ( 0, 0, 0 ) if ( scalar(@_) < 1 );
	return ( @_[0], @_[0], 0 ) if ( scalar(@_) < 2 );
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
	my $thread_output; #keeps filenames for concatenating data after all threads finished
	my @files_concatenate = ();

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
				$tot_threads_[$l__] = threads->create( {'context' => 'list'}, sub { return main_single_data_export( $tot_items[$l__], $path, @ARGV ) });
				$l__++;
				sleep 1; #give at least one second between requests
			}
			@run_threads_ = threads->list(threads::running());
			foreach my $oldth__ (@tot_threads_) {
				if ($oldth__->is_running()) {
				} elsif ($oldth__->is_joinable()) {
					$thread_output = $oldth__->join();
					push( @files_concatenate, $thread_output ) if ( defined($thread_output) && -e File::Spec->catfile( $path, $thread_output ) );
				}
			}
			@run_threads_ = threads->list(threads::running());
		}

		#clean up finished threads
		@run_threads_ = threads->list(threads::running());
		while (scalar @run_threads_ > 0) {
			foreach my $oldth__ (@tot_threads_) {
				if ($oldth__->is_joinable()) {
					$thread_output = $oldth__->join();
					push( @files_concatenate, $thread_output ) if ( defined($thread_output) && -e File::Spec->catfile( $path, $thread_output ) );
				}
			}
			$l__ = scalar @run_threads_;
			@run_threads_ = threads->list(threads::running());
			redo if ( scalar @run_threads_ == 0 && $l__ > 0 );
		}

		my $out_summary = File::Spec->catfile( $path, "track_data_summary_" . time() . ".tsv" );
		my $header = "";
		my @lines_outfile;
		for ( my $i=0; $i<scalar(@files_concatenate); $i++ ) {

			if ( open(FILE, "<" . File::Spec->catfile( $path, $files_concatenate[$i] ) ) ) {
				flock( FILE, LOCK_EX );
				$_ = <FILE>; #first line is header
				chomp;
				if ( $header eq "" ) {
					$header = $_;
				} elsif ( $header ne $_ ) {
					print "Mismatch between header when concatenating " . File::Spec->catfile( $path, $files_concatenate[$i] ) . ", ignoring this file.\n";
					flock( FILE, LOCK_UN );
					close( FILE );
					next;
				}

				$files_concatenate[$i] =~ s/\.single_track_data\.tsv$//;
				$files_concatenate[$i] =~ s/\.track_coordinates_in_time\.tsv$//;
				while( <FILE> ) {
					chomp;
					push( @lines_outfile, $files_concatenate[$i] . "\t" . $_ );
				}
				flock( FILE, LOCK_UN );
				close( FILE );

			} else {
				print( "Cannot open " . File::Spec->catfile( $path, $files_concatenate[$i] ) . " $!\n" );
				#return;
			}
		}

		if ( $header ne "" && scalar(@lines_outfile) > 0 && open(FILE, ">$out_summary" ) ) {
			flock(FILE, LOCK_EX);
			print FILE "Source\t" . $header;
			map { print FILE $_ . "\n"; } @lines_outfile;
			flock(FILE, LOCK_UN);
			close(FILE);
		} else {
			print "Error writing to $out_summary!\n";
			#return;
		}

	} elsif ( scalar(@files) == 1 ) {
		my $path = Cwd::cwd();
		main_single_data_export( $files[0], $path, @ARGV ); #do not capture returned value since no summary
	} else {
		print "Improperly formatted command line, no matching track coordinates files found.\n";
	}
#}




