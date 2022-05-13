#!/usr/bin/perl

# MaMuT dataset subtract spots from input list (using Name)
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Subtracts spots from a MaMuT dataset, with spots inputted from file(s) in radius annotation CSV format: Annotation,Timepoint,Position_T,ID,Name
# usage: perl MaMuT_dataset_subtract_spots_quickndirty.pl dataset_mamut.xml spots1.txt spots2.txt spots3.txt ...

# this is quick-n-dirty, and does not exhaustively ensure that entire tracks stay together -- therefore it is best used for SVF2MM output where TISSUE_TYPE is forward propagated through each entire track


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
	my @subtract_files = @_;
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
	
	#inpuut all spots to be subtracted
	my @annotated_spot_list;
	for ( my $i=0; $i<scalar(@subtract_files); $i++ ) {
			if ( open(FILE, "<$path/$subtract_files[$i]" ) ) {
				flock(FILE, LOCK_EX);
				my @file_lines;
				while( <FILE> ) {
					chomp;
					push( @file_lines, $_ );
				}
				flock(FILE, LOCK_UN);
				close(FILE);
				my @line_split;
				for (my $qq=1; $qq<scalar(@file_lines); $qq++ ) { #ignore header line
					@line_split = split( /,/, $file_lines[$qq] );
					push( @{$annotated_spot_list[$line_split[2]]}, $line_split[4] ); #timepoint, spot name
				
				}
			} else {
				print "Cannot open annotation file $path/$subtract_files[$i]!\n";
			}
			
			#my $count = 0;
			#map {$count += scalar(@{$_}); } @annotated_spot_list;
			#print "After adding $path/$subtract_files[$i], $count spots\n";
	}
	
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
							if ( $_ =~ /ID=\"(\d+)\"/i ) { # a spot 
								my $this_id = $1;
								my $going = TRUE;
								
								my $this_name;
								if ( $_ =~ /NAME=\"([\w\ ]+)?\"/i ) { # a spot name
									$this_name = $1;
								}
								
								for ( my $q=scalar(@{$annotated_spot_list[$this_frame]})-1; $q>=0; $q-- ) {
									if ( $annotated_spot_list[$this_frame][$q] eq $this_name ) {
										#spot on subtract list
										#print "found delete spot $this_name\n";
										$going = FALSE;
										splice( @{$annotated_spot_list[$this_frame]}, $q, 1 ); #remove from spots list since now deleted
										last;
									}
								}
								
								if ( $going ) {
									push( @this_block, $_ );
									$bins_number_of_spots[0]++;
									#push( @this_block_bins, 0 );
									push( @this_block_ids, $this_id );
								}
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
	my $print_this = TRUE;
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
		$print_this;
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
						$print_this = "";
						my $num_edges = 0;
						@spots_this_track = ();
						my $track_start_search_timepoint_index = 0;
						#my $spot_start_search_spot_index = 0;
						for ( my $k=0; $k<scalar(@{$lines_tracks[$i]}); $k++ ) {
							$going = TRUE;
							if ($lines_tracks[$i][$k] =~ /(^|^\s+)<Edge/ && $lines_tracks[$i][$k] =~ /SPOT_SOURCE_ID=\"(\d+)\".*?SPOT_TARGET_ID=\"(\d+)\"/) { #store all spots associated with 	
								my @spots_go = ( $1, $2 );
								my $found_spots = 0;
								#my $kill_loop = FALSE;
								#now, try to find each spot and ensure it is in this bin; if not, don't print the line to the file
								
								#for ( my $q=0; $q<scalar(@spots_go); $q++ ) {
								#	$kill_loop = FALSE;
								
								for ( my $q=$track_start_search_timepoint_index; $q<scalar(@spots_this_bin)+$track_start_search_timepoint_index; $q++ ) {
									#unless ( $q == $track_start_search_timepoint_index ) {
									#	$spot_start_search_spot_index = 0;
									#}
								
									#for ( my $qq=$spot_start_search_spot_index; $qq<scalar(@{$spots_this_bin[$q]})+$spot_start_search_spot_index; $qq++ ) {
									for ( my $qq=0; $qq<scalar(@{$spots_this_bin[$q]}); $qq++ ) {
										#print " Cmp $spots_this_bin[$qq] eq $spots_go[0] eq $spots_go[1] \n";
										if ( $spots_this_bin[$q][$qq] eq $spots_go[0] ) {
											#print "success at $track_start_search_timepoint_index: $spots_this_bin[$q][$qq] eq $spots_go[0]\n";
											$found_spots++;
											$bins_tracks[$i] .= "," . $l unless ( $bins_tracks[$i] =~ /(^|,)$l($|,)/ );
											#for ( my $rr=$qq+1; $rr<scalar(@{$spots_this_bin[$q]}); $rr++ ) {
											for ( my $r=$q+1; $r<scalar(@spots_this_bin); $r++ ) {
												for ( my $rr=0; $rr<scalar(@{$spots_this_bin[$r]}); $rr++ ) {
													if ( $spots_this_bin[$r][$rr] eq $spots_go[1] ) {
														$found_spots++;
														last;
													}
												}
												last if ( $found_spots >1);
											}
											
											#$kill_loop = TRUE;
											#print "found!\n";
											#last if ( $found_spots >1);
											#if ( $qq >= 0 ) {
											#	$spot_start_search_spot_index = $qq - scalar(@{$spots_this_bin[$q]});
											#} else {
											#	$spot_start_search_spot_index = $qq;
											#}
											
											if ( $q >= 0 ) {
												$track_start_search_timepoint_index = $q - scalar(@spots_this_bin);
											} else {
												$track_start_search_timepoint_index = $q;
											}
											last;
											
											
										}
									}
									last if ( $found_spots >0);
								}
								#	last if ( $kill_loop );
								#}
								
								if ( $found_spots < 2 ) {
									$going = FALSE;
									$num_edges--;
									
									last if ( $num_edges < -2 ); #if can't find 3 of edges, just declare this track not part of this bin
								} else {
									#print "yay $found_spots\n";
									push ( @spots_this_track, @spots_go );
									$num_edges++;
								}
							}
					
							if ( $going ) {
								#print "  " . $lines_tracks[$i][$k] . "\n";
								$print_this .= $lines_tracks[$i][$k] . "\n";
							}
						}
						
						if ( $num_edges > 0 ) {
							$this_xml .= $print_this;
							#print "  $num_edges edges\n";
						} else {
							#print "  ..track has no edges\n";
						}
						
						#remove seen spots for this track so later tracks go faster
						my $track_start_search_timepoint_index = 0;
						@spots_this_track = uniq ( @spots_this_track );
						for ( my $q=$track_start_search_timepoint_index; $q<scalar(@spots_this_bin)+$track_start_search_timepoint_index; $q++ ) {
						#for ( my $q=0; $q<scalar(@spots_this_bin); $q++ ) {
							#my $index = -1;
							for ( my $qq = scalar(@{$spots_this_bin[$q]}) -1 ; $qq >= 0; $qq-- ) {
								for ( my $p = scalar(@spots_this_track)-1; $p>=0; $p-- ) {
									if ( $spots_this_track[$p] eq $spots_this_bin[$q][$qq] ) {
										#$index = $qq;
										splice( @{$spots_this_bin[$q]}, $qq, 1 );
										splice( @spots_this_track, $p, 1 );
											
										if ( $q >= 0 ) {
											$track_start_search_timepoint_index = $q - scalar(@spots_this_bin);
										} else {
											$track_start_search_timepoint_index = $q;
										}
										last;
									}
								}
								
							}
						}
					}
				} elsif ( $fl == $line_number_filters ) {
					#for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						#if ( $bins_tracks[$i] == $l || $bins_tracks[$i] = -1 ) {
							for ( my $k=0; $k<scalar(@lines_filters); $k++ ) {
								if ( $bins_tracks[$i] =~ /(^|,)$l($|,)/ || $bins_tracks[$i] eq "-1" ) {
									$this_xml .= $lines_filters[$k] . "\n";
									#last;
								}
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
		if ( open(FILE, ">$path/$dataset_file\.subtracted_$timestamp\.xml" ) ) {
			flock(FILE, LOCK_EX);	
			print FILE $this_xml;
			flock(FILE, LOCK_UN);
			close(FILE);
		} else {
			print "Error writing to $path/$dataset_file\.subtracted_$timestamp\.xml!\n";
			return;
		}
		if ( open(FILE, ">$path/$dataset_file\.subtracted_$timestamp.spotList.txt" ) ) {
			flock(FILE, LOCK_EX);	
			print FILE "Annotation,Timepoint,Position_T,ID\n";
			print FILE $this_spot_list;
			flock(FILE, LOCK_UN);
			close(FILE);
		} else {
			print "Error writing to $path/$dataset_file\.subtracted_$timestamp\.spotList.txt!\n";
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
