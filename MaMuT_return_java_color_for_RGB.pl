#!/usr/bin/perl

# MaMuT return java / RGB colors
# 2019-2021 Martin H. Dominguez
# Gladstone Institutes

# Interconverts Java colors (used by MaMuT) and RGB
# usages: perl MaMuT_return_java_color_for_RGB.pl -235908 #(converts this java color to RGB)
#         perl MaMuT_return_java_color_for_RGB.pl 128,255,72 #(converts this RGB color to Java)
#         perl MaMuT_return_java_color_for_RGB.pl 128 255 72 #(converts this RGB color to Java)

use Cwd qw( cwd );
my $path = Cwd::cwd();

use constant TRUE => 1;
use constant FALSE => 0;

use constant {
	COLOR_ACTION_RANDOM => 0,
	COLOR_ACTION_SINGLE_FILL => 1,
	COLOR_ACTION_XY_VECTOR_FRAMES => 2,
	COLOR_ACTION_XZ_VECTOR_FRAMES => 3,
	COLOR_ACTION_YZ_VECTOR_FRAMES => 4,
	COLOR_ACTION_XYZ_VELOCITY_FRAMES => 5, #does not show vector, just velocity
	COLOR_ACTION_CELL_RADIAL_NESTING => 6, 
	COLOR_ACTION_CELL_DENSITY => 7,
	
	COLOR_ACTION_MAX => 8,
};

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



#==================
#Main subroutine
#==================
if ( scalar(@ARGV) == 3 ) {
	print "RGB -> Java result: " . hex_to_MaMuTdec(RGB_to_hex(@ARGV)) . "\n";
} elsif ( scalar(@ARGV) == 1 && $ARGV[0] =~ /(\d+),(\d+),(\d+)/ ) {
	print "RGB -> Java result: " . hex_to_MaMuTdec(RGB_to_hex($1,$2,$3)) . "\n";
} elsif ( scalar(@ARGV) == 1 ) {
	my $num = shift;
	$num = 0 - $num if ( $num > 0 );
	print "Java -> RGB result: " . join( ',', hex_to_RGB(MaMuTdec_to_hex($num)) ) . "\n";
}




