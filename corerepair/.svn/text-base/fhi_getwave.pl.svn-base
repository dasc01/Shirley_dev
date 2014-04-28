#!/usr/bin/perl -w

$usage = " Usage:
  fhi_getwave.pl legend xv.agr
    where legend is a string referring to legend of a wave function in
          xv.agr -- an xmgr output file from fhi98PP
";

if( @ARGV != 2 ) {
  print $usage;
  exit;
}

$legend=shift @ARGV;
$xvfile=shift @ARGV;
open XV, "<$xvfile" or die "unable to open:\n  $xvfile\n";

while($line=<XV>) {
  if( $line=~"^@ legend string" ) {
    chomp($line);
    @data1 = split /"/, $line;
    if( $legend =~ $data1[1] ) {
      @data2 = split ' ', $data1[1];
      $label = $data2[0];
      if( $label=~"s" ) {
        $l=0;
      } elsif( $label=~"p" ) {
        $l=1;
      } elsif( $label=~"d" ) {
        $l=2;
      } elsif( $label=~"f" ) {
        $l=3;
      } else {
        print "there are g-functions or greater in this file\n";
        print "this implementation only recognizes s,p,d,f\n";
        print "possible error... exiting\n";
        exit;
      }
      if( $#data2 > 0 ) {
        @data3 = split /=/, $data2[$#data2];
        $rcut = $data3[$#data3];
      } else {
        $rcut = 0.0;
      }
      push( @file, "$l $rcut\n" );
      $ngrid=0;
      while($line=<XV>) {
        if( $line=~"&" ) {
          last;
        }
        $file[$#file] .= $line;
        $ngrid++;
      }
      $file[$#file] = "$ngrid $file[$#file]";
    }
  }
}
$nwfc=@file;
print "$nwfc\n";
print @file;

exit;
