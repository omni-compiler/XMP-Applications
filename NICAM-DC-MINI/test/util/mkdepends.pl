#!/usr/bin/perl
#===============================================================================
#
#   mkdepends.pl: dependency list generator     
#
#===============================================================================
#
# Description
#   This program is called from Makefile in each custum project
#   Originated from "Make generator (S.Iga)"
#
# history
#    2011.8.8 (H.yashiro): ver.1
#
#===============================================================================

$topprg = trim($ARGV[0]).".f90";
$dirprg = "../".trim($ARGV[0]);

################################################################################
#
################################################################################

@dirs = ();
@pathlist = split(/:/,$ARGV[1]);

unshift( @dirs, $dirprg );
foreach $vpath ( @pathlist ) {
  unshift( @dirs, trim($vpath) );
};

################################################################################
# search all dependency
################################################################################

@alllist = ();
@hlist   = ();
$outhfile = '0';

# search all directory
foreach $dir (@dirs) {
  my $s  = "";
  my $s2 = "";

  # open directory
  opendir(DIR, $dir) or die "directory is not found! :".$dir;

  # read directory and search all files
  hoge:while( defined( $file = readdir(DIR) ) ) {
    if ( $file =~ /\.f90$/ or $file =~ /\.f$/ or $file =~ /\.c$/ ) {
      # open file
      open(FILE, "$dir".'/'."$file") or die "file is not found! : "."$dir".'/'."$file";

      $filename = $file;
      $dirname  = $dir;
      $modname  = '';
      $usedmods = '';

      # read all lines in the file
      while ( defined($line = <FILE>) ) {

        # file contains program?
        if ( $line =~ /^\s*program/i ) {
          $s  = $line;
          $s  =~ s/^\s*program//i;
          $s  =~ s/!.*//;
          $s2 = trim($s);

          $modname = $s2;
        };
        # file contains module?
        if ( $line =~ /^\s*module/i ) {
          $s = $line;
          $s =~ s/^\s*module//i;
          $s =~ s/!.*//;
          $s2 = trim($s);

          $modname = $s2;
        };

        # store all dependency
        if ( $line =~ /^\s*use/i ) {
          $s  = $line;
          $s  =~ s/^\s*use//i;
          $s  =~ s/!.*//;
          $s  =~ s/,.*//;
          $s2 = trim($s);

          $flag = '0';
          @filelist = split(/&/,$usedmods);

          foreach $item (@filelist) { 
            if ( $s2 eq $item ) { $flag='1' }
          };

          if ( $flag == '0' ) {
            $usedmods = $usedmods.$s2.'&';
          };
        };
      };

      # special treatment for f77 file (isccp)
      if ( $modname eq "mod_rd_mstrnx_ar5" or $modname eq "mod_rd_mstrnx_ar5_chaser" ) {
        $usedmods = $usedmods.'sub_icarus&sub_scops&';
      };

      # special treatment for c file (fio)
      if ( $modname eq "mod_fio" ) {
        $usedmods = $usedmods.'fio&fiof&';
      };

      if ( $file =~ /\.f$/ ) {
        $s = $file;
        $s =~ s/\.f$//;

        $modname = $s;
      };

      if ( $file =~ /\.c$/ ) {
        $s = $file;
        $s =~ s/\.c$//;

        $modname  = $s;
      };

      unshift( @alllist, [$filename,$dirname,$modname,$usedmods] );
#     print $dirname {$i},"/",$filename{$i}," : ",$modname {$i},"\n";
#     print $usedmods{$i},"\n";

      close(FILE);
    };

    # header file
    if ( $file =~ /\.h$/ ) {
      unshift( @hlist, [$file,$dir] );
    };
  };

  closedir(DIR)
}

################################################################################
## genertate depend.info
################################################################################

@filelist = ();
@dirlist  = ();

open(MAKEFILE, "> depends.info" );

print MAKEFILE '### Generated by mkdepend.pl'."\n\n";

# print dependency
&printdep( $topprg, $dirprg );

# print MODS list
print MAKEFILE "\nMODS =";
foreach $item (@filelist) { 
  if ( $item ne $topprg and $item ne /\.h$/ ){
    $s = $item;
    $s =~ s/\.f90$/.o/;
    $s =~ s/\.f$/.o/;
    $s =~ s/\.c$/.o/;

    print MAKEFILE "\t\\\n\t",$s;
  }
}
print MAKEFILE "\n\n";

# print MODSRC list
print MAKEFILE "\MODSRC =";

if ( $outhfile == '1' ) {
  $printdup = "";
  foreach $item (@hlist) { 
    $flag = '0';
    @duplist = split(/&/,$printdup);
    foreach $dupitem (@duplist) { 
      if ( @{$item}[0] eq $dupitem ) { $flag='1' }
    };

    if ( $flag == '0' ) {
      print MAKEFILE "\t\\\n\t",@{$item}[1]."/".@{$item}[0];
      $printdup = $printdup.@{$item}[0].'&'
    };
  };
};

foreach $item (@dirlist) { 
  if ( $item ne $dirprg."/".$topprg ){
    print MAKEFILE "\t\\\n\t",$item;
  }
}
print MAKEFILE "\n";

print "dependency generation OK.\n";

# END

#===============================================================================
# subroutines
#===============================================================================

sub printdep{
  my $targetfile = @_[0];
  my $targetdir  = @_[1];

  my $parent_key;
  my @parent_usedmod;
  my $usedmod;
  my $key;
  my $item;
  my @ufilelist;
  my $iufile;
  my $outline;

  # ignore duplicate
  foreach $item (@filelist) {
    if ( $item eq $targetfile ) {
      return;
    };
  };

  push( @filelist, $targetfile );
  push( @dirlist , $targetdir."/".$targetfile );

  foreach $allitem (@alllist) {
    if ( @{$allitem}[0] eq $targetfile ) {
      $s = $targetfile;
      $s =~ s/\.f90$/.o/;
      $s =~ s/\.f$/.o/;
      $s =~ s/\.c$/.o/;

      $outline  = $s."\t: ".$targetfile;
      @ufilelist = ();
      $iufile    = 0;

      # special treatment for c header
      if ( $targetfile eq "fio.c" ) {
        $outline  = $outline." fio.h fio_def.h";
        $outhfile = '1';
      }
      if ( $targetfile eq "fiof.c" ) {
        $outline  = $outline." fiof.h fio.h fio_def.h";
        $outhfile = '1';
      }

      @parent_usedmod = split(/&/,@{$allitem}[3]);

      foreach $usedmod (@parent_usedmod) { 
        foreach $item (@alllist) {
          if ( @{$item}[2] eq $usedmod ){
            $ufilelist[$iufile] = @{$item}[0];
            $iufile             = $iufile + 1;

            &printdep( @{$item}[0], @{$item}[1] );

            last;
          };
        };
      };

      foreach $item (@ufilelist) {
        $s = $item;
        $s =~ s/\.f90$/.o/;
        $s =~ s/\.f$/.o/;
        $s =~ s/\.c$/.o/;

        $outline = $outline." ".$s;
      };
      $outline = $outline."\n";

      print MAKEFILE $outline;

      last;
    };
  };
}

sub trim{
  my ($out) = @_;

  $out =~ s/^\s+//;
  $out =~ s/\s+$//;

  return $out
}

sub repeatstr{
  my $a = $_[0];
  my $b = $_[1];
  my $out = '';

  while ( $b>0 ) {
    $out = $out."$a";
    $b = $b - 1;
  };

  return $out
}

sub shrinkdir{  ## x/y/z x/y/g $out  ==> $out is ../g/
  my $a = $_[0];
  my $b = $_[1];
  my $out = '';
  my $i;
  my $c;
  my $d;
  my @d1;
  my @d2;

  @d1 = split(/\//,$a);
  @d2 = split(/\//,$b);

  while ( (@d1[0] eq @d2[0]) and (@d1) and (@d2) ) {
    $c = shift(@d1);
    $c = shift(@d2);
  }
  $d = @d1;
  $c = repeatstr('../',$d);
  $out = $c.join('/',@d2);

  if (@d2) { $out = $out.'/' };

  return $out
}  
