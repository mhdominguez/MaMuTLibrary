#!/usr/bin/perl

# MaMuT break divisions
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Break all divisions in a MaMuT dataset, not preserving even the closest daughter links
# usage: perl MaMuT_dataset_break_divisions.pl mamut_dataset.xml

use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;


#==================
#Main subroutine
#==================
sub main {
	#my (  @arguments ) = @_;
	my $dataset_file = shift;
	#my @subtract_files = @_;
	my @lines_main; #unchanging lines during reconstruction
	my $line_number_spots = -1; #indicates where to put <AllSpots> blocks when reconstructing the file at the end
	my @lines_spots;
	#my @lines_spots_bins;
	my @lines_spots_ids;
	my @frames_spots;
	my @bins_number_of_spots;
	
	my $line_number_tracks = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_tracks;
	#my @track_tracks;
	#my @spots_in_tracks;
	my @track_in_tracks;
	my $line_number_filters = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_filters;
	my @track_in_filters;
	
	my $num_bins = 1;
	my @bins_colors = ();
	#$num_bins = $_[1] if ( defined($_[1]) && $_[1] =~ /\d/ );
	my @bins_tracks;
	my $this_bin;
	#my $open_block = -1; #will determine actions to take while parsing XML line-by-line
	

	
	#parse dataset.mamut
	#my $start_time = time();
	my %colors_bins; #hash that stores all the important data	
	if ( open(FILE, "<$path/$dataset_file" ) ) {
		print "Subtracting spot list contained in " .join(',',@subtract_files)." from dataset $dataset_file.\n";
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
						my $this_frame = $2;
						push( @frames_spots, $this_frame );
						my @this_block = ( $_ );
						#my @this_block_bins = ( -1 );
						my @this_block_ids = ( -1 );
						while( <FILE> ) {
							chomp;
							if ( $_ =~ /ID=\"(\d+)\"/ ) { # a spot 
								my $this_id = $1;
								#my $going = TRUE;
							
								#if ( $going ) {
									push( @this_block, $_ );
									$bins_number_of_spots[0]++;
									#push( @this_block_bins, 0 );
									push( @this_block_ids, $this_id );
								#}
							} else {
								#probably a control block
								push( @this_block, $_ );
								#push( @this_block_bins, -1 );
								push( @this_block_ids, -1 );
							}
							last if ( $_ =~ /(^|^\s+)<\/SpotsInFrame/ );
						}
						push( @lines_spots, \@this_block );
						#push( @lines_spots_bins, \@this_block_bins );
						push( @lines_spots_ids, \@this_block_ids );
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
							}
							#} elsif ($_ =~ /(^|^\s+)<Edge/) { #store all spots associated with this track for reconstruction
							#	if ( $_ =~ /SPOT_SOURCE_ID=\"(\d+)\".*?SPOT_TARGET_ID=\"(\d+)\"/ ) {
							#		#print "here\n";
							#		push( @{$spots_in_tracks[scalar(@lines_tracks)]}, $1, $2 );
							#	}
							#}
						}
						#print "scalar " . scalar(@{$spots_in_tracks[scalar(@lines_tracks)]}) . "\n";
						#@{$spots_in_tracks[scalar(@lines_tracks)]} = uniq(@{$spots_in_tracks[scalar(@lines_tracks)]});
						$track_in_tracks[scalar(@lines_tracks)] = $this_track_id;
						$bins_tracks[scalar(@lines_tracks)] = "-1";
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

	#finally, for each subdivision, update nspots and write file
	my $nspots_divided = ( 0 ) x $num_bins;
	my $going;
	my @spots_this_bin;
	my @spots_this_track;
	my $this_xml = "";
	my $print_this;
	#for (my $l=0; $l<$num_bins; $l++ ) {
	#	print "bins spots $l = $bins_number_of_spots[$l]\n";
	#}
	for (my $l=0; $l<$num_bins; $l++ ) {
		#print "bins spots $l = $bins_number_of_spots[$l]\n";
		
		next unless ( defined($bins_number_of_spots[$l]) && $bins_number_of_spots[$l] =~ /\d/ && $bins_number_of_spots[$l] > 0 );
		for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) { #look to find and update nspots
			if ( $lines_main[$fl] =~ /(^|^\s+)<AllSpots/ ) {
				$lines_main[$fl] =~ s/nspots=\"(\d+)\"/nspots=\"$bins_number_of_spots[$l]\"/;
			}
		}
		
		#make a record of all spots in this bin
		@spots_this_bin = ();
		$this_xml = "";
		$this_spot_list = ""; #Annotation,Timepoint,Position_T,ID,Name
		#$print_this;
		#print "output bin $l...\n";

			for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) {
				if ( $fl == $line_number_spots ) {
					for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
						for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) {
							#if ( $lines_spots_bins[$i][$k] == $l) {
								$this_xml .= $lines_spots[$i][$k] . "\n";
								if ( $lines_spots_ids[$i][$k] >= 0 )  {
									push( @{$spots_this_bin[$i]}, $lines_spots_ids[$i][$k] );
									$this_spot_list .= "N\/A,$i,$frames_spots[$i],$lines_spots_ids[$i][$k]\n";
								}
							#}
						}
					}
					
					#print spots for cross-dataset subtraction i.e. slicing and dicing
					#map { $total_spots += scalar(@{$_}) } @spots_this_bin;
					#my $total_spots = 0;
					#map { $total_spots += scalar(@{$_}) } @spots_this_bin;
					#if ( $total_spots < 2 ) {
					#	$print_this = FALSE;
					#	last;
					#}
				} elsif ( $fl == $line_number_tracks ) {
					for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						#print "  bin $l, track: $i, spots to go at frame 0: " .scalar(@{$spots_this_bin[0]}) ." \n";
						#my @track_print_lines = ();
						#$print_this = "";
						my $num_edges = 0;
						@spots_this_track = ();
						my $track_start_search_timepoint_index = 0;
						#my $spot_start_search_spot_index = 0;
						for ( my $k=0; $k<scalar(@{$lines_tracks[$i]}); $k++ ) {
							#$going = TRUE;
							if ($lines_tracks[$i][$k] =~ /(^|^\s+)<Edge/ && $lines_tracks[$i][$k] =~ /SPOT_SOURCE_ID=\"(\d+)\".*?SPOT_TARGET_ID=\"(\d+)\"/) { #store all spots associated with 	
								$spots_this_track[$k] = $1;
								#$num_edges++;
							} else {
								$spots_this_track[$k] = "______control_block______";
							}
							#push( @track_print_lines, $lines_tracks[$i][$k] );
						}
						
						$print_this = "";
						for ( my $k=0; $k<scalar(@{$lines_tracks[$i]}); $k++ ) {
							$going = TRUE;
							if ( $spots_this_track[$k] ne "______control_block______" ) { 
								for ( my $kk=0; $kk<scalar(@{$lines_tracks[$i]}); $kk++ ) { 
									next if ( $k == $kk );
									if ( $spots_this_track[$k] eq $spots_this_track[$kk] ) {
										$going = FALSE;
										last;
									}
								}
							}
							unless ( $going ) {
								$print_this .= $lines_tracks[$i][$k] . "\n";
								$num_edges++;
							}
						}
						
						if ( $num_edges > 0 ) {
							#now, delete the lines with splits!

						
						
							$this_xml .= $print_this;
							#print "  $num_edges edges\n";
						} else {
							#print "  ..track has no edges\n";
						}
					}
				} elsif ( $fl == $line_number_filters ) {
					#for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						#if ( $bins_tracks[$i] == $l || $bins_tracks[$i] = -1 ) {
							for ( my $k=0; $k<scalar(@lines_filters); $k++ ) {
								#if ( $bins_tracks[$i] =~ /(^|,)$l($|,)/ || $bins_tracks[$i] eq "-1" ) {
									$this_xml .= $lines_filters[$k] . "\n";
									#last;
								#}
							}
						#}
					#}	
					#for ( my $i=0; $i<scalar(@{$filters_divided[$l]}); $i++ ) {
					#	$this_xml .= $filters_divided[$l][$i] . "\n";
					#}
				} #else {
				
				$this_xml .= $lines_main[$fl] . "\n";
				#print "  " . $lines_main[$fl] . "\n";
				#}
			}
		my $timestamp = time();
		if ( open(FILE, ">$path/$dataset_file\.broken_div_$timestamp\.xml" ) ) {
			flock(FILE, LOCK_EX);	
			print FILE $this_xml;
			flock(FILE, LOCK_UN);
			close(FILE);
		} else {
			print "Error writing to $path/$dataset_file\.broken_div_$timestamp\.xml!\n";
			return;
		}
		if ( open(FILE, ">$path/$dataset_file\.broken_div_$timestamp.spotList.txt" ) ) {
			flock(FILE, LOCK_EX);	
			print FILE "Annotation,Timepoint,Position_T,ID\n";
			print FILE $this_spot_list;
			flock(FILE, LOCK_UN);
			close(FILE);
		} else {
			print "Error writing to $path/$dataset_file\.broken_div_$timestamp\.spotList.txt!\n";
			return;
		}
		

	}
	
	#$start_time = time();
	#print "Time D: " . abs($end_time - $start_time) . " s\n";

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
