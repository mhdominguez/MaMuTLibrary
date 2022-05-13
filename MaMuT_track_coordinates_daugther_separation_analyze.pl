#!/usr/bin/perl

# MaMuT track coordinates daughter separation analyze
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Analyzes tracks for separation of daughters (best used with manually annotated tracks; TGMM is not sufficiently accurate with cell divisions to be relied on for rigorous quantification)
#   prints a table of daughter-mother distances, and daughter-daughter distances
# usage: perl MaMuT_track_coordinates_daughter_separation_analyze.pl dataset_mamut.track_coordinates_in_time.tsv start=20 stop=120

# this script is very slow!

use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;

use constant {
	FEATURE_TRACK_DURATION => 0,
	FEATURE_TRACK_START => 1,
	FEATURE_TRACK_STOP => 2,
	FEATURE_TRACK_DISPLACEMENT => 3,
	FEATURE_CELL_DISTANCE => 4,

	FEATURE_TRACK_MAX => 5,
};

my $density_radius = 12;
my $velocity_window = 5;
my $timeframe_per_hour = 10;

#==================
#Main subroutine
#==================
sub main {
	my $dataset_file = shift;
	
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
		} elsif ( $this_param[0] =~ /track_displacement/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_DISPLACEMENT] = $this_param[1]; 
		} elsif ( $this_param[0] =~ /cell_dist/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_CELL_DISTANCE] = $this_param[1]; 			
		} elsif ( $this_param[0] =~ /track_duration/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$track_feature[FEATURE_TRACK_DURATION] = $this_param[1]; 
		} elsif ( $this_param[0] =~ /den/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$density_radius = $this_param[1]; 
		} elsif ( $this_param[0] =~ /vel/i && $this_param[1] =~ /^[+-]?\d*\.?\d*([eE][-+]?\d+)?$/ ) {
			$velocity_window = $this_param[1]; 
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
								#if ( $track_feature[FEATURE_TRACK_DURATION] > 1 ) {
								#	if ( $track_feature[FEATURE_TRACK_DURATION] > $u - $t - 1 ) { #doesn't make the cut
								#		$going = -1;
								#		last;
								#	} else {
								#		#$going = 1;
								#	}
								#}
								#if ( $track_feature[FEATURE_TRACK_DURATION] > 1 ) {
									if ( $track_feature[FEATURE_TRACK_DISPLACEMENT] > sqrt((($xyz_end[0]-$xyz_begin[0])**2) + (($xyz_end[1]-$xyz_begin[1])**2) + (($xyz_end[2]-$xyz_begin[2])**2)) ) { #doesn't make the cut
										$going = -1;
										last;
									} else {
										#$going = 1;
									}
								#}
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
	
	print scalar(@lines_infile) . "\n"; #" $track_feature[FEATURE_TRACK_STOP]\n";
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
					my @xyzt_start;
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
							#print "here\n";
							if ($q==0) {
								push( @output_table, [($i-1, sqrt((($xyzt_this[0]-$xyzt_begin[0])**2) + (($xyzt_this[1]-$xyzt_begin[1])**2) + (($xyzt_this[2]-$xyzt_begin[2])**2)), "" )] );
								@xyzt_start = @xyzt_this;
							} else {
								push( @output_table, [($i-1, sqrt((($xyzt_this[0]-$xyzt_begin[0])**2) + (($xyzt_this[1]-$xyzt_begin[1])**2) + (($xyzt_this[2]-$xyzt_begin[2])**2)), sqrt((($xyzt_this[0]-$xyzt_start[0])**2) + (($xyzt_this[1]-$xyzt_start[1])**2) + (($xyzt_this[2]-$xyzt_start[2])**2)) )] );
							}
							push(@xyzt_multi_end, [@xyzt_this]);
						}
					}
				} else {
					print "Two daughters not found for track $i at timepoint $t!\n";
				}


				$going = $t;
			}# else {
			#	#original cell's track not defined at this timepoint, so it is the real end
			#	last;
			#}
		}
		
		#average velocities here
		#print scalar(@velocities) . ": " . $xyzt_begin[4] . "->" . "$xyzt_end[4] \n";
		#if (scalar(@velocities)==0) {
		#	print "problem!$velocity_window\n";
		#}
		#if ( $going > $xyzt_begin[4] + $track_feature[FEATURE_TRACK_DURATION] ) { #here, track must last at least track_duration
		#	push( @output_table, [($i-1,@xyzt_begin[0..4],@xyzt_end[0..4],$peak_track_displacement,$timeframe_per_hour*(set_stats(@velocities))[0]/$velocity_window),(set_stats(@densities))[0]] );
		#	#print "track $i, average vel: " . $output_table[$#output_table][13] . ", " .scalar(@densities)."\n";
		#}
	}

	if ( open(FILE, ">$path/$dataset_file\.daughter_separation_data\.tsv" ) ) {
		flock(FILE, LOCK_EX);
		print FILE "Tracks\tMother-Daughter Distance\tDaughter-Daugther Distance\n";
		map { print FILE join("\t",@{$_}) . "\n"; } @output_table;
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error writing to $path/$dataset_file\.daughter_separation_data\.tsv!\n";
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


