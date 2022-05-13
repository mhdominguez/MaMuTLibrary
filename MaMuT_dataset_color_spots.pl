#!/usr/bin/perl

# MaMuT color spots
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Colors spots and links in a MaMuT dataset, according to instruction provided...
# usages: perl MaMuT_dataset_color_spots.pl mamut_dataset.xml col=-3197141 #(fixed Java color all spots/links)
#                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml col=255,255,128 #(fixed Java color all spots/links)
#                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml rnd #(random color by track)
#                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml vec_xz=4 #(color track directionality converting angles to hue in single plane xy yz xz, by moving window size in timeframes)
#                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml tis #(color by tissueID i.e. populated by SVF)
#                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml div #(color division nodes)
#                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml den=4 #(color by density, by number of radii units around each cell to count other cells)
#                perl MaMuT_dataset_color_spots.pl mamut_dataset.xml vel=4 #(color by velocity, by moving window size in timeframes)

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
use constant {
	COLOR_ACTION_RANDOM => 0,
	COLOR_ACTION_SINGLE_FILL => 1,
	COLOR_ACTION_DIVISION_NODE => 2,
	COLOR_ACTION_TISSUE_ID => 3, #will hash or otherwise code TISSUE_TYPE parameter within XML to a numerical color
	COLOR_ACTION_XY_VECTOR_FRAMES => 4,
	COLOR_ACTION_XZ_VECTOR_FRAMES => 5,
	COLOR_ACTION_YZ_VECTOR_FRAMES => 6,
	COLOR_ACTION_XYZ_VELOCITY_FRAMES => 7, #does not show vector, just velocity
	COLOR_ACTION_CELL_RADIAL_NESTING => 8, #this is pretty much deprecated, better to use cell density
	COLOR_ACTION_CELL_DENSITY => 9,
	
	COLOR_ACTION_MAX => 10,
};

use constant PALETTE_COLORS => ( 
	[ 226, 196, 33 ],
	[ 133, 234, 245 ],
	[ 247, 215, 162 ],
	[ 0, 151, 174 ],
[ 207, 55, 43 ],
[ 129, 232, 0 ],
[ 236, 53, 254 ],
	[ 255, 135, 50 ],
[ 0, 121, 153 ],
[ 92, 141, 35 ],
[ 197, 183, 171 ],
[ 226, 243, 202 ],
[ 132, 108, 89 ],
[ 0, 194, 249 ],
[ 167, 221, 152 ],
[ 214, 214, 214 ],

);

#==================
#Color value processing subroutines
#==================
sub MaMuTdec_to_hex {
	return sprintf("0x%06X", shift()+16777216 );	
}

sub hex_to_MaMuTdec {
	return hex(shift)-16777216;
}

sub hex_to_RGB {
	my $hex = shift;
	if ( $hex =~ /0x([0-9A-Fa-f][0-9A-Fa-f])([0-9A-Fa-f][0-9A-Fa-f])([0-9A-Fa-f][0-9A-Fa-f])/ ) {
		my ( $rh, $gh, $bh ) = ( $1, $2, $3 );
		#my ( $rd, $gd, $bd ) = 
		return ( hex("0x$rh"), hex("0x$gh"), hex("0x$bh") );		
	} else {
		print "hex_to_RGB: bad number format: $hex!\n";
	}
}

sub RGB_to_hex {
	my ( $rd, $gd, $bd ) = @_;
	if ( $rd >= 0 && $gd >= 0 && $bd >= 0 && $rd < 256 && $gd < 256 && $bd < 256 ) {
		return "0x" . sprintf("%02X", $rd ) . sprintf("%02X", $gd ) . sprintf("%02X", $bd );
	} else {
		print "RGB_to_hex: bad number format: ($rd, $gd, $bd)!\n";
	}
}

sub random_color {
	return(20+int(rand()*236),20+int(rand()*236),20+int(rand()*236));
}

sub hsv_to_rgb {
	my ($h, $s, $v) = @_;
	my ($red, $green, $blue);
 
	if($v == 0) {
		($red, $green, $blue) = (0, 0, 0);
	} elsif($s == 0) {
		($red, $green, $blue) = ($v, $v, $v);
	} else {
		my $hf = $h / 60;
		my $i = int($hf);
		my $f = $hf - $i;
		my $pv = $v * (1 - $s);
		my $qv = $v * (1 - $s * $f);
		my $tv = $v * (1 - $s * (1 - $f));
 
		if($i == 0) {
			$red = $v;
			$green = $tv;
			$blue = $pv;
		} elsif($i == 1) {
			$red = $qv;
			$green = $v;
			$blue = $pv;
		} elsif($i == 2) {
			$red = $pv;
			$green = $v;
			$blue = $tv;
		} elsif($i == 3) {
			$red = $pv;
			$green = $qv;
			$blue = $v;
		} elsif($i == 4) {
			$red = $tv;
			$green = $pv;
			$blue = $v;
		} elsif($i == 5) {
			$red = $v;
			$green = $pv;
			$blue = $qv;
		} elsif($i == 6) {
			$red = $v;
			$blue = $tv;
			$green = $pv;
		} elsif($i == -1) {
			$red = $v;
			$green = $pv;
			$blue = $qv;
		} else {
			die('Invalid HSV -> RGB conversion.')
		}
	}
 
		return (int($red*255+0.5),int($green*255+0.5),int($blue*255+0.5));
}

sub hsv_to_rgb_broken { #input is H:0to360(degrees), S:0to1, V:0to1 -- subroutine is now broken!!

	my ($h, $s, $v) = @_;
	my ($red, $green, $blue);

	if($v == 0) {
		($red, $green, $blue) = (0, 0, 0);
	} elsif($s == 0) {
		($red, $green, $blue) = ($v, $v, $v);
	} else {
		my $hf = ( (int($h)% 6) / 60 ) ;
		my $i = int($hf) ;
		my $f = $hf - $i;
		my $pv = $v * (1 - $s);
		my $qv = $v * (1 - $s * $f);
		my $tv = $v * (1 - $s * (1 - $f));
		
		if($i == 0) {
			$red = $v;
			$green = $tv;
			$blue = $pv;
		} elsif($i == 1) {
			$red = $qv;
			$green = $v;
			$blue = $pv;
		} elsif($i == 2) {
			$red = $pv;
			$green = $v;
			$blue = $tv;
		} elsif($i == 3) {
			$red = $pv;
			$green = $qv;
			$blue = $v;
		} elsif($i == 4 ) {
			$red = $tv;
			$green = $pv;
			$blue = $v;
		} elsif($i == 5 ) {
			$red = $v;
			$green = $pv;
			$blue = $qv;
		#} elsif($i == 6) {
			#$red = $v;
			#$blue = $tv;
			#$green = $pv;
		#} elsif($i == -1) {
		#	$red = $v;
		#	$green = $pv;
		#	$blue = $qv;
		} else {
			die("Invalid HSV -> RGB conversion. $i")
		}
	}

	return (int($red*255+0.5),int($green*255+0.5),int($blue*255+0.5));
}
sub hue_2_rgb($$$){
   my ($v1,$v2,$vh)=@_;

   $vh+=1.0 if $vh < 0.0;
   $vh-=1.0 if $vh > 1.0;

   return $v1 + ($v2 - $v1) * 6.0 * $vh if 6.0 * $vh < 1;
   return $v2                           if 2.0 * $vh < 1;
   return $v1 + ($v2 - $v1) * ((2.0/3.0) - $vh)*6.0 if 3.0 * $vh < 2;
   return $v1;
}

sub hsl_to_rgb {

  my ($h,$s,$l) = @_;
  my ($r,$g,$b)=(0,0,0);

  if($h<0 || $s>1
  || $s<0 || $h>1
  || $l<0 || $l>1){
    die "bad input for hsl_to_rgb $h $s $l";
  }
  
  #$h = 1-$h;

  my ($var_1,$var_2);

  if (!$s){
      $r=$g=$b=$l;
  }else{
      if( $l < 0.5 ){ $var_2 = $l * ( 1 + $s ); }
      else          { $var_2 = ( $l + $s ) - ( $s * $l ) };

      $var_1 = 2.0 * $l - $var_2;

      $r = hue_2_rgb( $var_1, $var_2, $h + ( 1.0 / 3.0 ) ) ;
      $g = hue_2_rgb( $var_1, $var_2, $h );
      $b = hue_2_rgb( $var_1, $var_2, $h - ( 1.0/ 3.0 ) );
  }
  return (int($r*255+0.5),int($g*255+0.5),int($b*255+0.5));
  #return ($r,$g,$b)
}


use constant PI    => 4 * atan2(1, 1);
sub vec2_distance_and_angle { #returns the distance between two spots, and the angle of the vector (0deg is horizontal right to left)
	my ( $x_displacement, $y_displacement ) = @_;
	return ( sqrt(($x_displacement**2) + ($y_displacement**2)), atan2($y_displacement,$x_displacement)*180/PI );
}

sub vec3_distance { #returns the distance between two spots, and the angle of the vector (0deg is horizontal right to left)
	my ( $x_displacement, $y_displacement, $z_displacement ) = @_;
	return ( sqrt(($x_displacement**2) + ($y_displacement**2) + ($z_displacement**2)) );
}

sub rotate_rgb_color_hue {
	my ( $r, $g, $b, $angle ) = @_;
	$r/= 255; $g /= 255; $b /= 255;
	my $vsu = cos($angle * PI / 180 );
	my $vsw = sin($angle * PI / 180 );
print $vsu . " " . $vsw . "\n";
    my $ret_r = (.299 + .701*$vsu + .168*$vsw)*$r
        +   (.587 - .587*$vsu + .330*$vsw)*$g
        +   (.114 - .114*$vsu - .497*$vsw)*$b;
    my $ret_g = (.299 - .299*$vsu - .328*$vsw)*$r
        +   (.587 + .413*$vsu + .035*$vsw)*$g
        +   (.114 - .114*$vsu + .292*$vsw)*$b;
    my $ret_b = (.299 - .300*$vsu + 1.25*$vsw)*$r
        +   (.587 - .588*$vsu - 1.05*$vsw)*$g
        +   (.114 + .886*$vsu - .203*$vsw)*$b;
print join(',', ($r,$g,$b) ) . "\n";
	return ( 255* $ret_r, 255*$ret_g, 255*$ret_b );

}


	
#==================
#Main subroutine
#==================

sub main {
	#my (  @arguments ) = @_;
	my $dataset_file = $_[0];
	my $action = $_[1];
	my @lines_main; #unchanging lines during reconstruction
	my $line_number_spots = -1; #indicates where to put <AllSpots> blocks when reconstructing the file at the end
	my @lines_spots;
	my @frames_spots;
	
	my $line_number_tracks = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_tracks;
	#my @track_tracks;
	my @spots_in_tracks;
	my @spots_in_tracks_timepoints;
	my @spots_in_tracks_line_spots_index;
	my @track_in_tracks;
	my $line_number_filters = -1; #indicates where to put <AllTracks> blocks when reconstructing the file at the end
	my @lines_filters;
	my @track_in_filters;
	
	#my $num_bins = 2;
	my @bins_tracks;
	my $going;
	#my $open_block = -1; #will determine actions to take while parsing XML line-by-line
	
	#global vars for use with setting color correctly
	my $color_action;
	my $color_action_input_value;
	my $color_action_input_value_2;
	my @random_color_in_tracks;
	my $this_color_text;
	
	#command line parameter handler
	if ( defined($action) && $action =~ /\w/ ) {
		if ($action =~ /(ran|rnd)/i ) {
			$color_action = COLOR_ACTION_RANDOM;
		} elsif ( $action =~ /vec/i && $action =~ /(xy|yz|xz)/i ) {
			my $axis = $1;
			if ( $axis =~ /xz/ ) {
				$color_action = COLOR_ACTION_XZ_VECTOR_FRAMES;
			} elsif ( $axis =~ /yz/ ) {
				$color_action = COLOR_ACTION_YZ_VECTOR_FRAMES;
			} else { #default to xy
				$color_action = COLOR_ACTION_XY_VECTOR_FRAMES;
			}
			
			if ( $action =~ /=(\d+)/ ) {
				$color_action_input_value = $1;
			} else {
				$color_action_input_value = 4; #default is 4 frames back to see vector
			}
		} elsif ( $action =~ /vel/i || $action =~ /(xyz)/i ) {
			$color_action = COLOR_ACTION_XYZ_VELOCITY_FRAMES;
			if ( $action =~ /=(\d+)/ ) {
				$color_action_input_value = $1;
			} else {
				$color_action_input_value = 4; #default is 4 frames back to see vector
			}			
		} elsif ( $action =~ /den/i ) {
			$color_action = COLOR_ACTION_CELL_DENSITY;
			if ( $action =~ /=(\d+\.?\d*)/ ) {
				$color_action_input_value = $1;
			} else {
				$color_action_input_value = 12; #default is 12 radii distance for density estimate per cell
			}			
		} elsif ( $action =~ /col/i || $action =~ /div/i ) {
			$color_action = COLOR_ACTION_SINGLE_FILL;
			if ($action =~ /div/i) { $color_action = COLOR_ACTION_DIVISION_NODE; }
			
			if ( $action =~ /=(\d+),(\d+),(\d+)/ ) {
				$color_action_input_value = hex_to_MaMuTdec(RGB_to_hex( $1, $2, $3 ));
				#$color_action_input_value_2 = hex_to_MaMuTdec(RGB_to_hex( int($1/2), int($2/2), int($3/2) ));
				$color_action_input_value_2 = hex_to_MaMuTdec(RGB_to_hex( 255-$1, 255-$2, 255-$3 ));
			} elsif ( $action =~ /=\#(\w\w\w\w\w\w)/ ) {
				my $hexcode = $1;
				$color_action_input_value = hex_to_MaMuTdec( "0x" . $hexcode );
				my ( $r,$g,$b) = hex_to_RGB( "0x" . $hexcode );
				$color_action_input_value_2 = hex_to_MaMuTdec(RGB_to_hex( 255-$r, 255-$g, 255-$b ));
			} elsif ( $action =~ /=(-\d+)/ ) {
				$color_action_input_value = $1;
				my @rev_rgb = hex_to_RGB(MaMuTdec_to_hex($color_action_input_value));
				$color_action_input_value_2 = hex_to_MaMuTdec(RGB_to_hex( 255-$rev_rgb[0], 255-$rev_rgb[1], 255-$rev_rgb[2] ));
			} else {
				print "When using col= or color=, need to specify single solid fill color as RGB col=rrr,ggg,bbb in range 0-255\n";
				return;
			}						
		} elsif ($action =~ /tis/i) {
			$color_action = COLOR_ACTION_TISSUE_ID; 
			$action = "tissue ID";
		} else {
			$action = "random";
			$color_action = COLOR_ACTION_RANDOM;
		}
	} else {
		$action = "random";
		$color_action = COLOR_ACTION_RANDOM;
	}
	
	my $this_spot_id;
	my @lines_spots_ids;
	my @lines_spots_tissues;
	
	my $average_velocity = 0, $number_of_edges = 0;
	#parse dataset.mamut
	if ( open(FILE, "<$path/$dataset_file" ) ) {
		print "Coloring spots in dataset $dataset_file by $action ($color_action_input_value).\n";
		my $cell1; my $cell2; my $cell1_tp; my $cell2_tp; my $cell1_tp_idx; my $cell2_tp_idx;
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
						my @this_block_ids = ( -1 );
						my @this_block_tissues = (-1 );
						while( <FILE> ) {
							chomp;
							$_ =~ s/MANUAL_COLOR=\"([+-]?\d+(\.?\d*))\"//g; #remove prior colors
							if ( $_ =~ /<Spot/ && $_ =~ /ID=\"(\d+)\"/ ) {
								push( @this_block_ids, $1 );
							} else {
								push( @this_block_ids, -1 );
							}
							if ( $color_action == COLOR_ACTION_TISSUE_ID && $_ =~ /TISSUE_TYPE=\"(.*?)\"/ ) {
								push( @this_block_tissues, $1 );
							}
							push( @this_block, $_ );
							last if ( $_ =~ /(^|^\s+)<\/SpotsInFrame/ );
						}
						push( @lines_spots, \@this_block );
						push( @lines_spots_ids, \@this_block_ids );
						push( @lines_spots_tissues, \@this_block_tissues ) if ( $color_action == COLOR_ACTION_TISSUE_ID );
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
						my $this_track_tissue;
						my $this_track_id = $2;
						$_ =~ s/MANUAL_COLOR=\"([+-]?\d+(\.?\d*))\"//g; #remove prior colors
						my @this_block = ( $_ );
						#print "inside here: $this_track_id\n";
						while( <FILE> ) {
							chomp;
							push( @this_block, $_ );
							if ( $_ =~ /(^|^\s+)<\/Track/ ) {
								last;
							} elsif ($_ =~ /(^|^\s+)<Edge/) { #store all spots associated with this track for reconstruction
								if ( $_ =~ /SPOT_SOURCE_ID=\"(\d+)\".*?SPOT_TARGET_ID=\"(\d+)\"/ ) {
									
									$cell1=$1; $cell2=$2;
									$cell1_tp = -1; $cell2_tp=-1; $cell1_tp_idx = -1; $cell2_tp_idx = -1;
									#now find out where these spots are in time
									for ( my $i=0; $i < scalar(@lines_spots_ids); $i++ ) {
										
										for ( my $k=0; $k<scalar(@{$lines_spots_ids[$i]}); $k++ ) {
											if ( $lines_spots_ids[$i][$k] eq $cell1 ) {
												$cell1_tp = $i;
												$cell1_tp_idx = $k;
												if ($color_action == COLOR_ACTION_TISSUE_ID && defined($lines_spots_tissues[$i][$k]) ) {
													if ( $this_block[$#this_block] !~ /TISSUE_TYPE=/ ) {
														$this_block[$#this_block] =~ s/\/>/TISSUE_TYPE="$lines_spots_tissues[$i][$k]" \/>/;
													} elsif ( $this_block[$#this_block] =~ /TISSUE_TYPE=\"0\"/ ) {
														$this_block[$#this_block] =~ s/TISSUE_TYPE=\"(\d+)\"/TISSUE_TYPE="$lines_spots_tissues[$i][$k]"/;
													}
												}
											}
											if ( $lines_spots_ids[$i][$k] eq $cell2 ) {
												$cell2_tp = $i;
												$cell2_tp_idx = $k;
												#if ($color_action == COLOR_ACTION_TISSUE_ID && defined($lines_spots_tissues[$i][$k]) && $this_block[$#this_block] !~ /TISSUE_TYPE=/ ) {
												#	$this_block[$#this_block] =~ s/\/>/TISSUE_TYPE="$lines_spots_tissues[$i][$k]" \/>/;
												#}
												if ($color_action == COLOR_ACTION_TISSUE_ID && defined($lines_spots_tissues[$i][$k]) ) {
													if ( $this_block[$#this_block] !~ /TISSUE_TYPE=/ ) {
														$this_block[$#this_block] =~ s/\/>/TISSUE_TYPE="$lines_spots_tissues[$i][$k]" \/>/;
													} elsif ( $this_block[$#this_block] =~ /TISSUE_TYPE=\"0\"/ ) {
														$this_block[$#this_block] =~ s/TISSUE_TYPE=\"(\d+)\"/TISSUE_TYPE="$lines_spots_tissues[$i][$k]"/;
													}
												}
											}
											last if ( $cell1_tp >=0 && $cell2_tp >=0 );
										}
										last if ( $cell1_tp >=0 && $cell2_tp >=0 );
									}
									
									push( @{$spots_in_tracks[scalar(@lines_tracks)]}, $cell1, $cell2 ); #each pair is important even,odd because we will reconstruct lineages this way
									push( @{$spots_in_tracks_timepoints[scalar(@lines_tracks)]}, $cell1_tp, $cell2_tp ); #each pair of timepoints for each spots_in_tracks cell
									push( @{$spots_in_tracks_lines_spots_index[scalar(@lines_tracks)]}, $cell1_tp_idx, $cell2_tp_idx ); #each pair of timepoints for each spots_in_tracks cell
									#print "Pushing on spots_in_tracks_lines_spots_index: $cell1_tp_idx, $cell2_tp_idx\n";
								}
								if ( $_ =~ /VELOCITY=\"([+-]?\d+(\.?\d*))\"/ ) {
									$average_velocity += $1;
									$number_of_edges++;
								}
							}
						}
						#print "scalar " . scalar(@{$spots_in_tracks[scalar(@lines_tracks)]}) . "\n";
						#@{$spots_in_tracks[scalar(@lines_tracks)]} = uniq(@{$spots_in_tracks[scalar(@lines_tracks)]});
						$track_in_tracks[scalar(@lines_tracks)] = $this_track_id;
						$bins_tracks[scalar(@lines_tracks)] = 0; #all tracks go to bin 1, track-less spots go to bin 0
						if ( $color_action == COLOR_ACTION_RANDOM ) {
							$this_color_text = "MANUAL_COLOR=\"" . hex_to_MaMuTdec(RGB_to_hex(random_color())) . "\"";
							map{ $_ =~ s/\s\/>/ $this_color_text \/>/ if ($_ =~ /(^|^\s+)<Edge/); } @this_block;
							$random_color_in_tracks[scalar(@lines_tracks)] = $this_color_text;
						} elsif ( $color_action == COLOR_ACTION_SINGLE_FILL || $color_action == COLOR_ACTION_DIVISION_NODE ) {
							$this_color_text = "MANUAL_COLOR=\"" . $color_action_input_value . "\"";
							map{ $_ =~ s/\s\/>/ $this_color_text \/>/ if ($_ =~ /(^|^\s+)<Edge/); } @this_block;
							$random_color_in_tracks[scalar(@lines_tracks)] = $this_color_text;
						} elsif ( $color_action == COLOR_ACTION_TISSUE_ID ) {
							#do nothing, have to find-and-replace at the very end
						}
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
				if ( $_ =~ /<Feature feature=\"POSITION_X\"/ ) { #add additional cell features here
					if( $color_action == COLOR_ACTION_CELL_DENSITY ) {
						push( @lines_main, qq~        <Feature feature="SPOT_DENSITY" name="Spot Density ~ . $color_action_input_value . qq~ Radii" shortname="Density~ . $color_action_input_value . qq~" dimension="QUALITY" isint="false" />~ );
					} elsif( $color_action == COLOR_ACTION_CELL_RADIAL_NESTING ) {
						push( @lines_main, qq~        <Feature feature="SPOT_RADIAL_NESTING" name="Average Spot Distance" shortname="Density" dimension="LENGTH" isint="false" />~ );
					} elsif( $color_action == COLOR_ACTION_XYZ_VELOCITY_FRAMES ) {
						push( @lines_main, qq~        <Feature feature="AVG_VELOCITY_~ . $color_action_input_value . qq~" name="Track Velocity (~ . $color_action_input_value . qq~ frames avg)" shortname="Trk Velocity" dimension="LENGTH" isint="false" />~ );
					}
				}
				push( @lines_main, $_ );
			}

		}
		flock(FILE, LOCK_UN);
		close(FILE);
	} else {
		print "Error opening $path/$dataset_file!\n";
		return;
	}
	
	$average_velocity /= $number_of_edges if ( $number_of_edges > 0 );
	$average_velocity *= $color_action_input_value;
	
	my @lines_spots_divided; #final storage for subdivisions, 3D matrix [timepoint][subdivision][spot], compared with 2D matrix of @lines_spots [timepoint][spot]
	
	if ( scalar(@spots_in_tracks) <= 0 && ( $color_action == COLOR_ACTION_SINGLE_FILL || $color_action == COLOR_ACTION_DIVISION_NODE ) ) { #here, don't worry about tracks, need to fill even track-less spots
		for ( my $i=0; $i<scalar(@lines_spots_ids); $i++ ) {
			for ( my $k=0; $k<scalar(@{$lines_spots_ids[$i]}); $k++ ) {
				push( @{$spots_in_tracks[0]}, $lines_spots_ids[$i][$k] );
				push( @{$spots_in_tracks_timepoints[0]}, $i );
				push( @{$spots_in_tracks_lines_spots_index[0]}, $k );
			}
		}
		push( @bins_tracks, 0 );
		$this_color_text = "MANUAL_COLOR=\"" . $color_action_input_value . "\"";
		$random_color_in_tracks[0] = $this_color_text;		
	}

	my %tissue_id_color;
	if ( $color_action == COLOR_ACTION_TISSUE_ID ) {
		#set up hash for tissue id
		my $current_id = 0; #index of PALETTE_COLORS that we are on
		for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
			for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) { #find all coordinates first
				if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /TISSUE_TYPE=\"(.*?)\"/ ) {
					my $this_type = $1;
					unless ( exists($tissue_id_color{$this_type}) ) {
						if ( $current_id >= scalar(PALETTE_COLORS) ) {
							print "Too many tissue types for current palette!\n";
							return;						
						} else {
							$tissue_id_color{$this_type} = hex_to_MaMuTdec(RGB_to_hex( @{(PALETTE_COLORS)[$current_id]} ));
						}
						$current_id++;
					}
					
					$lines_spots[$i][$k] =~ s/\/>/MANUAL_COLOR="$tissue_id_color{$this_type}" \/>/;
				}
			}
		}
		#assume all tissue IDs are represented in the cells themselves, none are only present in just the edges
		for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
			for ( my $k=0; $k<scalar(@{$lines_tracks[$i]}); $k++ ) {
				if ( $lines_tracks[$i][$k] =~ /<Edge/ && $lines_tracks[$i][$k] =~ /TISSUE_TYPE=\"(.*?)\"/ ) {
					my $this_type = $1;
					$lines_tracks[$i][$k] =~ s/\/>/MANUAL_COLOR="$tissue_id_color{$this_type}" \/>/ if exists($tissue_id_color{$this_type});
				}
			}
		}
	} elsif ( $color_action == COLOR_ACTION_CELL_DENSITY || $color_action == COLOR_ACTION_CELL_RADIAL_NESTING ) {
		my $distance_sum;
		my $distance_numbers;
		my $this_distance;
		my $this_spot_density;
		#my $density_text;
		my @cell_positions = (); #array of 3d vectors, one for each cell
		my @cell_densities = ();
		my @all_densities = (); #running total for figuring out how to normalize
		my $all_densities_begin_this_timepoint;
		my ( $avg, $med, $std );
		for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
			#@cell_positions = ();
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
						push( @this_cell_coordinates, $1 );
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
			
			$all_densities_begin_this_timepoint = scalar(@all_densities);
			for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) { #now figure out average distance of other cells, for each cell
				if ( $lines_spots[$i][$k] =~ /<Spot/ && ref($cell_positions[$i][$k]) ) { #got a valid cell
					$distance_sum =0;
					$distance_numbers =0;
					
					if ($color_action == COLOR_ACTION_CELL_DENSITY) {
						for ( my $kk=0; $kk<scalar(@{$lines_spots[$i]}); $kk++ ) {
							next if ( $kk == $k || !(ref($cell_positions[$i][$kk])) ); #only accept valid non-self cells beyond this point
							$this_distance = vec3_distance($cell_positions[$i][$k][0]-$cell_positions[$i][$kk][0],$cell_positions[$i][$k][1]-$cell_positions[$i][$kk][1],$cell_positions[$i][$k][2]-$cell_positions[$i][$kk][2]);
							if ( $this_distance < $cell_positions[$i][$k][3] * $color_action_input_value ) { #in range, count it
								$distance_numbers++;
								$distance_sum += ($cell_positions[$i][$k][3] * 2) ** 3;
							}
						}
					} else { #case COLOR_ACTION_CELL_RADIAL_NESTING
						for ( my $kk=0; $kk<scalar(@{$lines_spots[$i]}); $kk++ ) {
							next if ( $kk == $k || !(ref($cell_positions[$i][$kk])) ); #only accept valid non-self cells beyond this point
							$distance_sum += vec3_distance($cell_positions[$i][$k][0]-$cell_positions[$i][$kk][0],$cell_positions[$i][$k][1]-$cell_positions[$i][$kk][1],$cell_positions[$i][$k][2]-$cell_positions[$i][$kk][2]);
							$distance_numbers++;
						}
					}
						
					if ( $distance_numbers > 0 && $distance_sum > 0 ) {
						if ($color_action == COLOR_ACTION_CELL_DENSITY) { 
							#$this_spot_density = 1000 * ( $distance_numbers ** 2 ) / $distance_sum;
							$this_spot_density = $distance_numbers;
						} else { #case COLOR_ACTION_CELL_RADIAL_NESTING
							$this_spot_density = $distance_sum/$distance_numbers;
						}
						
						push( @all_densities, $this_spot_density );
						#$density_text = "SPOT_DENSITY=\"" . sprintf("%.8f", $this_spot_density) . "\"";
						#$lines_spots[$i][$k] =~ s/\s\/>/ $density_text \/>/;
						$cell_densities[$i][$k] = $this_spot_density;
					} else {
						if ($color_action == COLOR_ACTION_CELL_DENSITY) { 
							$cell_densities[$i][$k] = 0;
						} else {
							$cell_densities[$i][$k] = "N\/A";
						}
					}

				} else {
					$cell_densities[$i][$k] = "N\/A";
				}
			}
			#figure out average (or max) distance for this timeframe -- only for radial nesting algorithm
			next unless ($color_action == COLOR_ACTION_CELL_RADIAL_NESTING);
			
			if ( $#all_densities - $all_densities_begin_this_timepoint > 0) {
				( $avg, $med, $std ) = set_stats( @all_densities[$all_densities_begin_this_timepoint..$#all_densities] );
				$avg = max(@all_densities[$all_densities_begin_this_timepoint..$#all_densities]);
				if ( defined($avg) && defined($std) && $avg > 0 && $std > 0 ) {
					#normalize to the average (could consider normalizing to max but will do average for now
					map{ $_ /= $avg; } @all_densities[$all_densities_begin_this_timepoint..$#all_densities];
					map{ $_ /= $avg if ( $_ =~ /\d/); } @{$cell_densities[$i]};				
				} else {
					#throw out densities for time point $i
					splice( @all_densities, $all_densities_begin_this_timepoint, scalar(@all_densities) - $all_densities_begin_this_timepoint );
				}
			} else {
				#nothing to do, there are no densities reported for time point $i
			
			}
		}
		 
		( $avg, $med, $std ) = set_stats( @all_densities ) if ( scalar(@all_densities) > 0);
		if ( defined($avg) && defined($std) && $avg > 0 && $std > 0 ) {
			print "  ...average cell density number ($dataset_file) $avg (stdev $std, median $med)\n";
			my $max = $avg + ( 4*$std);
			my $min = $avg - ( 4*$std);
			$min = 0 if ( $min < 0 );
			undef @cell_positions;
			$max -= $min;
			if ( defined($max) && defined($min) && $max > 0 && $min >= 0 ) {
				#print "  ...painting spots\n";
				for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
					for ( my $k=0; $k<scalar(@{$lines_spots[$i]}); $k++ ) { #now paint the spots
						if ( $lines_spots[$i][$k] =~ /<Spot/ && defined($cell_densities[$i][$k]) && $cell_densities[$i][$k] =~ /\d/ ) { #got a valid cell
							$this_color_text = "SPOT_DENSITY=\"" . sprintf("%.8f", $cell_densities[$i][$k]) . "\"";
							#$lines_spots[$i][$k] =~ s/\s\/>/ $density_text \/>/;
							#print "     ....$this_color_text\n";
						
							$this_cell_density = ($cell_densities[$i][$k] - $min) / $max;
							$this_cell_density = 1 if ( $this_cell_density > 1);
							$this_cell_density = 0 if ( $this_cell_density < 0);
							
							$this_color_text .= " MANUAL_COLOR=\"" . hex_to_MaMuTdec(RGB_to_hex(hsl_to_rgb((240-($this_cell_density * 240))/360,1,0.5))) . "\"";
							#$this_color_text = hex_to_MaMuTdec(RGB_to_hex(255,0,0));
							$lines_spots[$i][$k] =~ s/\s\/>/ $this_color_text \/>/;
							#$lines_spots[$i][$k] =~ s/SPOT_DENSITY=/MANUAL_COLOR=\"$this_color_text\" SPOT_DENSITY=/;
						 
						}
					}
				}
			}
		}
	} else { #usual flow here, everything that requires lineage reconstruction
		#reconstruct lineages here
		my @lineage_tree; #array of array of array refreshed for each track, each item of this array represents a different time point, and is a pointer to an array of arrays containing each cell 
		my @check_for_linked_spots;
		my $this_timepoint_index;
		my @related_spots;
		my @related_spots_timepoints;
		my @related_spots_timepoints_indices;
		
		for ( my $l=0; $l<scalar(@spots_in_tracks); $l++ ) { #iterate over tracks first, spots in track second; ignore timepoints for now ##TODO: multithread this
			#push( @lines_spots_ids, \@this_block_ids )
			next if ($bins_tracks[$l] < 0 );
			#@finished_spots = ();
			
			for ( my $j=0; $j<scalar(@{$spots_in_tracks[$l]}); $j++ ) {
				#see if we have seen this spot before
				$going =FALSE;
				for ( my $k=0; $k<$j; $k++ ) {
					$going = TRUE if ( $spots_in_tracks[$l][$k] eq $spots_in_tracks[$l][$j] );
				
				}
				next if ( $going);
	
				my @add_to_array = ( $spots_in_tracks[$l][$j], $spots_in_tracks_timepoints[$l][$j], $spots_in_tracks_lines_spots_index[$l][$j] );
				#print "Add to array init: $spots_in_tracks[$l][$j], $spots_in_tracks_timepoints[$l][$j], $spots_in_tracks_lines_spots_index[$l][$j]\n";
				@check_for_linked_spots = ($j);
				@related_spots = ();
				@related_spots_timepoints = ();
				@related_spots_timepoints_indices = ();
				#find all identical spot instances to this one
				for ( my $k=$j+1; $k<scalar(@{$spots_in_tracks[$l]}); $k++ ) {
					if ( $spots_in_tracks[$l][$k] eq $spots_in_tracks[$l][$j] ) {
						#$going = FALSE;
						#map {} @
						push( @check_for_linked_spots, $k );# if ($going);
					}# && $k != $j );
				}
				
				#okay now go through all instances of this spot to recoved linked spots 
				for ( my $k=0; $k<scalar(@check_for_linked_spots); $k++ ) {
					if ( $check_for_linked_spots[$k] % 2 == 0 ) { #is on even item, its pair is $r+1
						push( @related_spots, $spots_in_tracks[$l][$check_for_linked_spots[$k]+1] );
						push( @related_spots_timepoints, $spots_in_tracks_timepoints[$l][$check_for_linked_spots[$k]+1] );
						push( @related_spots_timepoints_indices, $spots_in_tracks_lines_spots_index[$l][$check_for_linked_spots[$k]+1] );
						
						#if ( $spots_in_tracks_timepoints[$l][$check_for_linked_spots[$k]+1] > $spots_in_tracks_timepoints[$l][$j]  ) {
						#	print "Up Adding spot $spots_in_tracks[$l][$check_for_linked_spots[$k]+1], $spots_in_tracks_timepoints[$l][$check_for_linked_spots[$k]+1], $spots_in_tracks_lines_spots_index[$l][$check_for_linked_spots[$k]+1]\n";
						#	print "  $related_spots_timepoints_indices[$#related_spots_timepoints_indices] vs. $related_spots_timepoints_indices[$#related_spots]\n";
						#}
					} else { #is on odd item, its pair is $r-1 -- by the way this spot will not usually get added since it is typically the one that predates the other in the pair
						push( @related_spots, $spots_in_tracks[$l][$check_for_linked_spots[$k]-1] );
						push( @related_spots_timepoints, $spots_in_tracks_timepoints[$l][$check_for_linked_spots[$k]-1] );
						push( @related_spots_timepoints_indices, $spots_in_tracks_lines_spots_index[$l][$check_for_linked_spots[$k]-1] );
						
						#print "Down Adding spot $spots_in_tracks[$l][$check_for_linked_spots[$k]-1], $spots_in_tracks_timepoints[$l][$check_for_linked_spots[$k]-1], $spots_in_tracks_lines_spots_index[$l][$check_for_linked_spots[$k]-1]\n" if ( $spots_in_tracks_timepoints[$l][$check_for_linked_spots[$k]-1] > $spots_in_tracks_timepoints[$l][$j]  );
					}
	
					if ( $related_spots_timepoints[$#related_spots] > $spots_in_tracks_timepoints[$l][$j]  ) { #this_spot -> other_spot
						#add to lineage tree
						#print "this add_to_array adding: $related_spots[$#related_spots], $related_spots_timepoints[$#related_spots], $related_spots_timepoints_indices[$#related_spots]\n";
						push( @add_to_array, $related_spots[$#related_spots], $related_spots_timepoints[$#related_spots], $related_spots_timepoints_indices[$#related_spots] ); #FYI $#related_spots = $k
						#push( @finished_spots, @check_for_linked_spots[$k] );
					} else {	#other_spot -> this_spot
						#ignore; we could consider using this information in reverse but would have to spider all around the dataset trying to make sure all links were preserved
					}				
				}
				
				#add to growing tree
				push( @{$lineage_tree[$spots_in_tracks_timepoints[$l][$j]]}, \@add_to_array );
				#use Data::Dumper; print Dumper(\@lineage_tree); return if ( $j > 50 ); #output lineage tree, for debug only
			}
			
			#use Data::Dumper; print Dumper(\@lineage_tree); return; #output lineage tree, for debug only
		
			#now that lineage tree is complete, can paint spots accordingly
			if ( $color_action == COLOR_ACTION_RANDOM || $color_action == COLOR_ACTION_SINGLE_FILL) { #entire tree is painted one color
				for ( my $i=0; $i<scalar(@lineage_tree); $i++ ) { #go forward in time
					for ( my $j=0; $j<scalar(@{$lineage_tree[$i]}); $j++ ) {
						#for ( my $k=0; $k<scalar(@{$lineage_tree[$i][$j]}); $k+=3 ) {
							$this_timepoint_index = $lineage_tree[$i][$j][2];
							#print "Replacing color for spot $lineage_tree[$i][$j][0] at timepoint $lineage_tree[$i][$j][1] with index $this_timepoint_index\n  $lines_spots[$i][$this_timepoint_index]\n" if ( $i==0 );
							$lines_spots[$i][$this_timepoint_index] =~ s/\s\/>/ $random_color_in_tracks[$l] \/>/;
							#print "  $lines_spots[$i][$this_timepoint_index]\n" if ( $i==0 );
							#unshift( @{$lines_spots_divided[$i][0]}, $lines_spots[$i][$this_timepoint_index] );
						#}
					}
				}
			} elsif ( $color_action == COLOR_ACTION_DIVISION_NODE ) { #colors only parents of divisions in exact same way as COLOR_ACTION_SINGLE_FILL
				for ( my $i=0; $i<scalar(@lineage_tree); $i++ ) { #go forward in time
					for ( my $j=0; $j<scalar(@{$lineage_tree[$i]}); $j++ ) {

						if (scalar(@{$lineage_tree[$i][$j]}) > 8) { #division node when there are two linked cells						
							#$this_timepoint_index = $lineage_tree[$i][$j][2];
							#$lines_spots[$i][$this_timepoint_index] =~ s/\s\/>/ $random_color_in_tracks[$l] \/>/;
							#next;
							
							$this_timepoint_index = $lineage_tree[$i][$j][2];	
							$lines_spots[$i][$this_timepoint_index] =~ s/\s\/>/ MANUAL_COLOR=\"$color_action_input_value\" \/>/;
							#next;
							#if ($i+1 == $lineage_tree[$i][$j][4]) {
							#	$this_timepoint_index = $lineage_tree[$i][$j][5];	
							#	$lines_spots[$i+1][$this_timepoint_index] =~ s/\s\/>/ MANUAL_COLOR=\"$color_action_input_value\" \/>/;
							#}
							
							#if ($i+1 == $lineage_tree[$i][$j][7]) {
							#	$this_timepoint_index = $lineage_tree[$i][$j][8];	
							#	$lines_spots[$i+1][$this_timepoint_index] =~ s/\s\/>/ MANUAL_COLOR=\"$color_action_input_value_2\" \/>/;
							#}
						} else {
							$this_timepoint_index = $lineage_tree[$i][$j][2];	
							$lines_spots[$i][$this_timepoint_index] =~ s/\s\/>/ MANUAL_COLOR=\"$color_action_input_value_2\" \/>/;
						}
					}
				}
			} elsif ( $color_action >= COLOR_ACTION_XY_VECTOR_FRAMES && $color_action <= COLOR_ACTION_XYZ_VELOCITY_FRAMES ) {
				my @gather_children, @add_to_array; #stores IDs of children during each timepoint to reconstruct each brach
				my @parent_coordinates, @child_coordinates;
				$this_color_text = "MANUAL_COLOR=\"" . hex_to_MaMuTdec(RGB_to_hex(0,0,0)) . "\"";
				for ( my $i=0; $i<scalar(@lineage_tree); $i++ ) { #go forward in time, start with everything black and appropriately added to the output bin
					for ( my $j=0; $j<scalar(@{$lineage_tree[$i]}); $j++ ) {
				#		#for ( my $k=0; $k<scalar(@{$lineage_tree[$i][$j]}); $k+=3 ) {
							$this_timepoint_index = $lineage_tree[$i][$j][2];
							#print "  ..$i: loop testing child $lineage_tree[$i][$j][0]: $lines_spots[$i][$this_timepoint_index]\n";
							#print "this lineage_tree[$i][$j]: " . join(",", @{$lineage_tree[$i][$j]} ) . "\n";
							$lines_spots[$i][$this_timepoint_index] =~ s/\s\/>/ $this_color_text \/>/;
				#			#unshift( @{$lines_spots_divided[$i][0]}, $lines_spots[$i][$this_timepoint_index] );
				#		#}
					}
				}
				#next;
				for ( my $i=0; $i<scalar(@lineage_tree)-$color_action_input_value; $i++ ) { #go forward in time to find spots and then find linked spots $color_action_input_value time frames later to subtract position
					for ( my $j=0; $j<scalar(@{$lineage_tree[$i]}); $j++ ) {
						#for ( my $k=0; $k<scalar(@{$lineage_tree[$i][$j]}); $k+=3 ) {
							#print( "$lineage_tree[$i][$j][2]," );
							next unless ( scalar(@{$lineage_tree[$i][$j]}) % 3 == 0 && scalar(@{$lineage_tree[$i][$j]}) > 3 ); #if there are no children, don't even bother
							$this_timepoint_index = $lineage_tree[$i][$j][2]; #parent cell is the first listed one
							@parent_coordinates = ( );
							if ( $color_action == COLOR_ACTION_XY_VECTOR_FRAMES ) {
								$parent_coordinates[0] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ );
								$parent_coordinates[1] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ );
								#$parent_coordinates[2] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ );
							} elsif ( $color_action == COLOR_ACTION_XZ_VECTOR_FRAMES ) {
								$parent_coordinates[0] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ );
								$parent_coordinates[1] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ );
							} elsif ( $color_action == COLOR_ACTION_YZ_VECTOR_FRAMES ) {
								$parent_coordinates[0] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ );
								$parent_coordinates[1] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ );
							} else { #case COLOR_ACTION_XYZ_VELOCITY
								$parent_coordinates[0] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ );
								$parent_coordinates[1] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ );
								$parent_coordinates[2] = $1 if ( $lines_spots[$i][$this_timepoint_index] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ );
							}
							#print "here1\n";					
							if ( defined($parent_coordinates[0]) && defined($parent_coordinates[1]) && $parent_coordinates[0] =~ /\d/ && $parent_coordinates[1] =~ /\d/ ) { #good coordinates for parent, so now find all children and paint based on their coordinates
								#print "Here2 parent " . $lineage_tree[$i][$j][0] . "\n";
								@gather_children = (); #initialize
								for ( my $k=3; $k<scalar(@{$lineage_tree[$i][$j]}); $k+=3 ) { #gather the immediate children, at time $i
									push( @gather_children, $lineage_tree[$i][$j][$k], $lineage_tree[$i][$j][$k+1], $lineage_tree[$i][$j][$k+2] );
								}
	
								for ( my $tp=$i+1; $tp<$i+$color_action_input_value; $tp++ ) {
									#check gather children for cells at each timepoint, delete them as their own children are added
									#print "  ..checking out timepoint $tp with " . scalar(@gather_children)/3 . " gathered children for tree $l starting at timepoint $i\n";
									@add_to_array = ();
									for ( my $k=scalar(@gather_children)-1; $k>0; $k-=3 ) {
										if ( $gather_children[$k-1] == $tp ) { #process only timepoint we are in
											#print "Debug, correct timepoint $gather_children[$k-1] == $tp\!n";
											#scan lineage tree for its own children
											#print "     ...child $gather_children[$k-2]\n";
											for ( my $jj = 0; $jj < scalar(@{$lineage_tree[$tp]}); $jj++ ) {
												if ( $lineage_tree[$tp][$jj][0] eq $gather_children[$k-2] && scalar(@{$lineage_tree[$tp][$jj]}) > 3 ) {
													#print "this lineage_tree[$tp][$jj]: " . join(",", @{$lineage_tree[$tp][$jj]} ) . "\n";
													for ( my $kk=3; $kk<scalar(@{$lineage_tree[$tp][$jj]}); $kk+=3 ) {
														push( @add_to_array, $lineage_tree[$tp][$jj][$kk], $lineage_tree[$tp][$jj][$kk+1], $lineage_tree[$tp][$jj][$kk+2] );
														#print "Adding to add_to_array($tp,$jj,$kk): $lineage_tree[$tp][$jj][$kk], $lineage_tree[$tp][$jj][$kk+1], $lineage_tree[$tp][$jj][$kk+2]:\n$lines_spots[$tp][$lineage_tree[$tp][$jj][$kk+2]]\n";
													}
												}
											}
											splice( @gather_children, $k-2, 3);
											#print( "Splicing..." . join( ',', splice( @gather_children, $k-2, 3) ) . "\n" ); #delete now that own children added
										} #else {
											#print "Debug, mismatched timepoint $gather_children[$k-1] == $tp\n";
										#}
									}
									push( @gather_children, @add_to_array );
									#print "Pushing " . scalar(@add_to_array)/3 . " new children at $tp...\n";
									#for ( my $k=scalar(@add_to_array)-1; $k>0; $k-=3 ) {
									#	print "  ..$tp: loop testing $k child $add_to_array[$k-2]: $lines_spots[$tp][$add_to_array[$k]]\n";
									#}
								}
								#print "Here3, with " . scalar(@gather_children)/3 . " gather children, first one is: " . $gather_children[0] . ".\n";
								next unless ( scalar(@gather_children) >= 3 );
								#now examine final timepoint
								my $tp = $i+$color_action_input_value;
								for ( my $k=scalar(@gather_children)-1; $k>0; $k-=3 ) {
									if ( $gather_children[$k-1] == $tp ) { #okay, coordinates for this child should be found and spot colored according to vector
										@child_coordinates = ( );
										$this_timepoint_index = $gather_children[$k];
										#print "  ..$tp: loop testing $k child $gather_children[$k-2]: $lines_spots[$tp][$this_timepoint_index]\n";
										if ( $color_action == COLOR_ACTION_XY_VECTOR_FRAMES ) {
											$child_coordinates[0] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ );
											$child_coordinates[1] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ );
										} elsif ( $color_action == COLOR_ACTION_XZ_VECTOR_FRAMES ) {
											$child_coordinates[0] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ );
											$child_coordinates[1] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ );
										} elsif ( $color_action == COLOR_ACTION_YZ_VECTOR_FRAMES ) {
											$child_coordinates[0] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ );
											$child_coordinates[1] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ );
										} else { #case COLOR_ACTION_XYZ_VELOCITY
											$child_coordinates[0] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_X=\"([+-]?\d+(\.?\d*))\"/ );
											$child_coordinates[1] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_Y=\"([+-]?\d+(\.?\d*))\"/ );
											$child_coordinates[2] = $1 if ( $lines_spots[$tp][$this_timepoint_index] =~ /POSITION_Z=\"([+-]?\d+(\.?\d*))\"/ );
										}
										#print "  ...with $child_coordinates[0], $child_coordinates[1]\n";
										if ( defined($child_coordinates[0]) && defined($child_coordinates[1]) && $child_coordinates[0] =~ /\d/ && $child_coordinates[1] =~ /\d/ ) { #good coordinates for parent, so now find all children and paint based on their coordinates
											if ( defined($child_coordinates[2]) && $child_coordinates[2] =~ /\d/ ) { #XYZ velocity, no length
												my $length = vec3_distance($child_coordinates[0]-$parent_coordinates[0],$child_coordinates[1]-$parent_coordinates[1],$child_coordinates[2]-$parent_coordinates[2]);

												$this_color_text = "AVG_VELOCITY_" . $color_action_input_value . "=\"" . sprintf("%.8f", $length/$color_action_input_value) . "\"";
												$lines_spots[$tp][$this_timepoint_index] =~ s/\s\/>/ $this_color_text \/>/;
												#print "length = $length, average vel=$average_velocity: $child_coordinates[0]-$parent_coordinates[0],$child_coordinates[1]-$parent_coordinates[1],$child_coordinates[2]-$parent_coordinates[2]\n";
												$length /= $average_velocity;# * $color_action_input_value/3;
												$length = 1 if ( $length > 1 );
												#printf "\x1b[48;2;%d;%d;%dm", rand 256, rand 256, rand 256;
												$this_color_text = hex_to_MaMuTdec(RGB_to_hex(hsl_to_rgb((240-($length * 240))/360,1,0.5)));
												#printf "\x1b[48;2;%d;%d;%dm",hex_to_RGB(MaMuTdec_to_hex(hex_to_MaMuTdec(RGB_to_hex(@rgb))));
												#print "length: $length";
												#print "\x1b[0m\n";
												##print "hsl length = " . (($length * 240)+120)/360 . ": $this_color_text\n";
											} else {
												my ( $length, $angle ) = vec2_distance_and_angle($child_coordinates[0]-$parent_coordinates[0],$child_coordinates[1]-$parent_coordinates[1]);
												#print "  ...to color child based on $angle and $length.\n";
												$angle = ($angle+180)/360;
												$length /= $average_velocity;# * $color_action_input_value /3;
												$length = 1 if ( $length > 1 );
												$this_color_text = hex_to_MaMuTdec(RGB_to_hex(hsl_to_rgb($angle,$length,0.5)));
											}
											$lines_spots[$tp][$this_timepoint_index] =~ s/MANUAL_COLOR=\"([+-]?\d+(\.?\d*))\"/MANUAL_COLOR=\"$this_color_text\"/;
										} else {
											#print "Failed gather children $gather_children[$k-2] / $gather_children[$k-1]: $lines_spots[$tp][$this_timepoint_index]\n";
										}
										
										splice( @gather_children, $k, 3); #delete now that own children added
									} #else {
										#print "  ..TP mismatch $tp vs. $gather_children[$k-1] for $gather_children[$k]\n";
									#}
								}							
							} #else {
							#	print "No parent coordinates\n";
							#}
							
							#$this_timepoint_index = $lineage_tree[$i][$j][$k+2];
							#$lines_spots[$i][$this_timepoint_index] =~ s/\s\/>/ $random_color_in_tracks[$l] \/>/;
							#unshift( @{$lines_spots_divided[$i][1]}, $lines_spots[$i][$this_timepoint_index] ); #in terms of pushing cells 
						#}
					}
				}
			}
		
			undef @lineage_tree;
		}
	}

	#strip filename of extension, and assume it is XML
	( undef, $dataset_file ) = split( /\./, reverse($dataset_file), 2 );
	$dataset_file = reverse($dataset_file);
	
	#iterate over spots at t=0, subdivide them, while storing IDs for tree reconstruction
	
	
	my $this_bin;

	#now, iterate over all spots in time to sort by bin
	#my $this_track; my $this_spot_in_track;
	
	#my $this_color_text;
	for ( my $i=0; $i<scalar(@lines_spots); $i++ ) {
		#print "t=$i, total spots handled: " . scalar(@{$lines_spots[$i]}) . "\n";
		for ( my $k=scalar(@{$lines_spots[$i]})-1; $k>=0; $k-- ) {
			if( $this_block_ids[$i][$k] >= 0 ) {#if ( $lines_spots[$i][$k] =~ /<Spot/ && $lines_spots[$i][$k] =~ /ID=\"(\d+)\"/ ) {
				#if ( $lines_spots[$i][$k] =~ // ) { #if this spot was part of a track, put it in bin 1, otherwise bin 0
				#	unshift( @{$lines_spots_divided[$i][1]}, $lines_spots[$i][$k] );
				#} else {
					unshift( @{$lines_spots_divided[$i][0]}, $lines_spots[$i][$k] );
				#}	
			} else {
				#control blocks get passed to both divisions
				#for (my $l=0; $l<$num_bins; $l++ ) {
					unshift( @{$lines_spots_divided[$i][0]}, $lines_spots[$i][$k] );
				#}
			}
		}
		undef @{$lines_spots[$i]}; #freeing space
	}
	#return;
	
	#finally, for each subdivision, update nspots and write file
	my $nspots_divided = 0;
		for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
			$nspots_divided += scalar(@{$lines_spots_divided[$i][0]}) -2;
		}
		for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) { #look to find and update nspots
			if ( $lines_main[$fl] =~ /(^|^\s+)<AllSpots/ ) {
				$lines_main[$fl] =~ s/nspots=\"(\d+)\"/nspots=\"$nspots_divided\"/;
			}
		}
		
		if ( open(FILE, ">$path/$dataset_file\.colored\.xml" ) ) {
			flock(FILE, LOCK_EX);
			for ( my $fl=0; $fl<scalar(@lines_main); $fl++ ) {
				if ( $fl == $line_number_spots ) {
					for ( my $i=0; $i<scalar(@lines_spots_divided); $i++ ) {
						for ( my $k=0; $k<scalar(@{$lines_spots_divided[$i][0]}); $k++ ) {
							print FILE $lines_spots_divided[$i][0][$k] . "\n";
						}
					}
				} elsif ( $fl == $line_number_tracks ) {
					for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						if ( $bins_tracks[$i] == 0 ) {
							print FILE join( "\n", @{$lines_tracks[$i]} ) . "\n";
						}
					}
				} elsif ( $fl == $line_number_filters ) {
					for ( my $i=0; $i<scalar(@lines_tracks); $i++ ) {
						if ( $bins_tracks[$i] == 0 ) {
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
			print "Error writing to $path/$dataset_file\.colored\.xml!\n";
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


#my $MANUAL_COLOR="-65281";
#print join( ',', hsl_to_rgb((240-(0.5 * 240))/360,1,.5)) . "\n";
#print join( ',', rotate_rgb_color_hue( @ARGV ) ) . "\n";

main(@ARGV);



