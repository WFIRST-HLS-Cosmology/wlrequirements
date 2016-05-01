#!/usr/bin/perl

# Settings
#
# "Skip step" settings (0/1; turn on only to speed up testing; should be all 0's for the actual run)
$SKIP_COEFS = 0;              # If set, forces the code to skip computation of the PSF coefficients
$SKIP_GAL = 0;                # If set, skips computation of the source galaxy-weighted PSF coefficients

# Over-write these choices if set from the command line
if (2<=scalar @ARGV) {
  ($SKIP_COEFS, $SKIP_GAL) = @ARGV;
}

# Filter table:
# (edges in microns)
open(IN, '../InputData/Filters.txt') or die "Can't open filters file\n";
while ($line=<IN>) {
  if ($line !~ m/^\#/) {
    # Un-commented lines in filter table
    @data = split ' ', $line;
    $filters{$data[0]} = "$data[1] $data[2]";
  }
}
close IN;

# Test of filter table input -- leave commented for normal runs
# for $f (keys %filters) {print "$f --> $filters{$f}\n";}

# Table of galaxies to use
$galfile = '../InputData/galaxies.dat';

# Galaxy size range
# [based on log10 reff in arcsec]
$log10reffmin = -1.60;
$dlog10reff = 0.2;
$Nlog10reff = 9;

# Redshift binning
$dz = 0.2;
$Nz = 15;

# Reference parameters
$ref_pupil = 1;
$ref_jit = 400; # image motion in mas^2 per axis
$Nlam = 1; # number of wavelength points
$wfe = 90.; # [nm] static wave front error to use in computing coefficients

# Build image ellipticity code:
#
system "rm main.x";
system "gcc main.c -o main.x -lm -Wall -O3";
if (not -e 'main.x') {die "Error: did not create main.x\n";}
#
# Calling format for the image ellipticity code:
# ./main.x pupil lambdamin lambdamax Nlam hexflags reff g1 g2 ZA0..ZA7 jitxx jitxy jityy Dx Dy alpha alphaplus

if ($SKIP_COEFS==0) {for $f (keys %filters) {

  open(OUT, ">coefficients_$f.dat");

  # Get derivatives of galaxy ellipticity with respect to each parameter
  # for galaxies of each size
  for $j (0..$Nlog10reff-1) {
    $reff = 10**($log10reffmin + $dlog10reff*$j);
    print "Processing ... filter = $f, reff = $reff arcsec\n";
    print OUT (sprintf "%8.5f %9.6f ", log($reff)/log(10.), $reff);

    # Get response1 = d(observed image ellipticity)/d(galaxy shear)
    $dg = .05;
    $Z = '0 'x8;
    @data0 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0 0.0 $Z  $ref_jit 0 $ref_jit 0 0 0 0`;
    @data1 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff -$dg 0.0 $Z  $ref_jit 0 $ref_jit 0 0 0 0`;
    @data2 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff $dg 0.0 $Z  $ref_jit 0 $ref_jit 0 0 0 0`;
    $response1 = ($data2[6]-$data1[6])/(2*$dg);
    print OUT (sprintf "%9.6f ", $response1);
    #print "@data1\n@data2\n\n";  # <-- TEST OUTPUT

    # Get d(observed image ellipticity)/d(foc*astig) [nm^-2]
    $Zt = sprintf "%10.7f %10.7f 0 0 0 0 0 0 ", $wfe/sqrt(2.), -$wfe/sqrt(2.);
    @data1 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Zt  $ref_jit 0 $ref_jit 0 0 0 0`;
    $Zt = sprintf "%10.7f %10.7f 0 0 0 0 0 0 ", $wfe/sqrt(2.), $wfe/sqrt(2.);
    @data2 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Zt  $ref_jit 0 $ref_jit 0 0 0 0`;
    #print "@data1\n@data2\n\n";  # <-- TEST OUTPUT
    $ellip_coef_fa = ($data2[6]-$data1[6])/$wfe**2;

    # Get d(observed image ellipticity)/d(coma*coma) [nm^-2]
    $Zt = sprintf "0 0 0 %10.7f 0 0 0 0 ", $wfe;
    @data1 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Zt  $ref_jit 0 $ref_jit 0 0 0 0`;
    $Zt = sprintf "0 0 0 0 %10.7f 0 0 0 ", $wfe;
    @data2 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Zt  $ref_jit 0 $ref_jit 0 0 0 0`;
    #print "@data1\n@data2\n\n";  # <-- TEST OUTPUT
    $ellip_coef_cc = ($data2[6]-$data1[6])/(-2*$wfe**2);

    # Get d(observed image ellipticity)/d(coma*trefoil) [nm^-2]
    $Zt = sprintf "0 0 0 %10.7f 0 0 %10.7f 0 ", $wfe/sqrt(2.), $wfe/sqrt(2.);
    @data1 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Zt  $ref_jit 0 $ref_jit 0 0 0 0`;
    $Zt = sprintf "0 0 0 %10.7f 0 %10.7f 0 0 ", $wfe/sqrt(2.), $wfe/sqrt(2.);
    @data2 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Zt  $ref_jit 0 $ref_jit 0 0 0 0`;
    #print "@data1\n@data2\n\n";  # <-- TEST OUTPUT
    $ellip_coef_ct = ($data2[6]-$data1[6])/$wfe**2;

    # Get d(observed image ellipticity)/d(sph*astig) [nm^-2]
    $Zt = sprintf "0 %10.7f 0 0 0 0 0 %10.7f ", $wfe/sqrt(2.), -$wfe/sqrt(2.);
    @data1 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Zt  $ref_jit 0 $ref_jit 0 0 0 0`;
    $Zt = sprintf "0 %10.7f 0 0 0 0 0 %10.7f ", $wfe/sqrt(2.), $wfe/sqrt(2.);
    @data2 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Zt  $ref_jit 0 $ref_jit 0 0 0 0`;
    #print "@data1\n@data2\n\n";  # <-- TEST OUTPUT
    $ellip_coef_sa = ($data2[6]-$data1[6])/$wfe**2;

    # Get d(observed image ellipticity)/d(jit*jit) [mas^-2]
    $jit1 = 1.5*$ref_jit; $jit2 = 0.5*$ref_jit;
    @data1 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Z  $jit1 0 $jit2 0 0 0 0`;
    @data2 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Z  $jit2 0 $jit1 0 0 0 0`;
    #print "@data1\n@data2\n\n";  # <-- TEST OUTPUT
    $ellip_coef_jj = ($data2[6]-$data1[6])/2./($jit2-$jit1);

    # Get d(observed image ellipticity)/d(alpha plus)
    $ap = .01;
    @data1 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Z  $ref_jit 0 $ref_jit 0 0 0 -$ap`;
    @data2 = split ' ', `./main.x $ref_pupil $filters{$f} $Nlam 0 $reff 0.0 0.0 $Z  $ref_jit 0 $ref_jit 0 0 0 $ap`;
    #print "@data1\n@data2\n\n";  # <-- TEST OUTPUT
    $ellip_coef_ap = ($data2[6]-$data1[6])/2./$ap;

    print OUT (sprintf "%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e ",
      $ellip_coef_fa/$response1, $ellip_coef_cc/$response1, $ellip_coef_ct/$response1, $ellip_coef_sa/$response1,
      $ellip_coef_jj/$response1, $ellip_coef_ap/$response1);

    print OUT "\n";
  }

  close OUT;
}}

# Galaxy tables
if ($SKIP_GAL==0) {for $f (keys %filters) {

  # Read galaxy file
  $script = q{sed -n '/^[^\#]/p'}; # Strip comments 
  @galaxies = split ' ', `$script $galfile`;
  $Ngal = scalar @galaxies;
  $Ngal = int($Ngal/4.+.00001);

  # Read coefficients into the mydata array
  open(IN, "coefficients_$f.dat");
  for $j (0..$Nlog10reff-1) {
    $line=<IN>;
    @data = split ' ', $line;
    $Ncol = scalar @data;
    for $k (0..$Ncol-1) {$mydata[$j][$k] = $data[$k];}
  }
  close IN;

  # Average of coefficients in each column (skip first 2)
  for $iz (0..$Nz-1) {$wtot[$iz] = 0; for $k (2..$Ncol-1) {$sum[$iz][$k] = 0;}}
  for $igal (0..$Ngal-1) {
    $iz = int($galaxies[4*$igal+1]/$dz+.00000001);
    if ($iz>=0 and $iz<$Nz) {
      $w = 1.;   # weight of this galaxy

      # Get size bin j:
      $jc = (log($galaxies[4*$igal+2]+1e-8)/log(10.) - $log10reffmin)/$dlog10reff;
      $j = int($jc+.5);
      if ($j<0) {$j=0;}
      if ($j>=$Nlog10reff) {$j=$Nlog10reff-1;}

      # TEST OUTPUT
      #if ($igal%10000==3500) {@temp = @galaxies[4*$igal..4*$igal+3]; print STDERR "$f: $igal, $iz, $j ($temp[0] $temp[1] $temp[2] $temp[3])\n";}

      $wtot[$iz] += $w;
      for $k (2..$Ncol-1) {$sum[$iz][$k] += $w*$mydata[$j][$k];}
    }
  }
  open(OUT, ">coefficients_gal_$f.dat");
  for $iz (0..$Nz-1) {
    print OUT (sprintf "%2d %6.4f %6.4f %12.5e  ", $iz, $iz*$dz, ($iz+1)*$dz, $wtot[$iz]);
    for $k (2..$Ncol-1) {
      print OUT (sprintf " %12.5e", $sum[$iz][$k]/$wtot[$iz]);
    }
    print OUT "\n";
  }
  close OUT;
}}
