#!/usr/bin/perl

use strict;

open(F,$ARGV[0]) || die "Cant open Haplotect summary file!";
<F>;
$_ = <F>;
chomp;
my @l = split("\t",$_);
my ($sampleid,$mle,$mle_ci,$totalSites,$meanCov) = ($l[0],$l[1],$l[2],$l[7],$l[8]);
close F;

my $informativeSites = 0;
my $contamCounts = 0;
my $totalReads = 0;

my %loci = ();
open(F,$ARGV[1]) || die "Cant open Haplotect count file!";
while(<F>){
    next if /^#/;
    chomp;
    my @l = split("\t",$_);

    next if defined $loci{"$l[0]:$l[1]:$l[2]"};

    $loci{"$l[0]:$l[1]:$l[2]"} = 1;

    my %haps = ();
    while($l[10] =~ /([ACGT]{2}):(\d+)/g){
	$haps{$1} = $2 if $2 > 1;
    }
    if (scalar keys %haps == 3){
#	print join("\t",@l),"\n";
	$informativeSites++;
	my @c = sort { $a <=> $b } values %haps;
	$contamCounts += $c[0];
	$totalReads += $l[9];
    } elsif (scalar keys %haps == 4){
#	print join("\t",@l),"\n";
	$informativeSites++;
	my @c = sort { $a <=> $b } values %haps;
	$contamCounts += (($c[0]+$c[1]) / 2);
	$totalReads += $l[9];
    }   
}

print join("\t",qw(#sample totalSites informativeSites meanCoverage contaminatingHaplotypeCount totalReads contaminationEstimate snpEstimate snpEstimateCI)),"\n";
print join("\t",$sampleid,$totalSites,$informativeSites,$meanCov,$contamCounts,$totalReads,$totalReads > 0 ? 2 * $contamCounts/$totalReads : 0,$mle,$mle_ci),"\n";
