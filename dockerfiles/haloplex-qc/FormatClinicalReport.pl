#!/usr/bin/perl

use strict;
use JSON;
use File::Basename;

my $REPORTVERSION = basename($0,".pl");

my %aa3to1 = qw(Ala A Arg R Asn N Asp D Asx B Cys C Glu E Gln Q Glx Z Gly G His H Ile I Leu L Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Xxx X Ter *);

my $variants = shift @ARGV;
my $qc = shift @ARGV;
my $reportinfo = shift @ARGV;

open(F,$qc) || die "Cant open QC file $qc";
my $dat = <F>;
chomp $dat;
my $h = decode_json $dat;

print join("\n","PATIENT ID: " . $h->{sample}, "SAMPLE ID: " . $h->{library},"DATE REPORTED: " . $h->{timestamp}),"\n\n";

my @tier12 = ();
my @tier12_full = ();
my @tier3 = ();
my @tier3_full = ();

my %mutgenes = ();
open(F,$variants) || die;
$_ = <F>;
chomp;
s/\s+$//;
s/\S+\.CVAF/CVAF/;

my @headers = split("\t",$_);
while(<F>){
    chomp;
    s/\s+$//;
    my @F = split("\t",$_);

    my %F = map { $headers[$_] => $F[$_] } 0..$#F;

    next if $F{Consequence}=~/synonymous/ || $F{FILTER} ne 'PASS';

    my $HGVSpShort = $F{HGVSp};
    $HGVSpShort =~ s/\S+:p\.(\S+?)/\1/g; 

    while( $HGVSpShort and my ( $find, $replace ) = each %aa3to1 ) {
	eval "\$HGVSpShort =~ s{$find}{$replace}g";
    }
    $F{HGVSp_Short} = $HGVSpShort;
    
    $F{Consequence} = (split("&",$F{Consequence}))[0];
    $F{Consequence}=~/^(\S+?)_/; 
    my $var=uc($1); 
    $var .= " INSERTION" if length($F{ALT})>length($F{REF}); 
    $var .= " DELETION" if length($F{ALT})<length($F{REF});
    $var .= " SITE VARIANT" if $var =~ /SPLICE/;
    $F{EXON}=~s/\/\d+//g;
    $F{INTRON}=~s/\/\d+//g; 

    my $cat = ($F{MYELOSEQ_MDS_AC}>0 || $F{MYELOSEQ_TCGA_AC}>0) ? "TIER 1-2" : "TIER 3";

    $F{CVAF} = sprintf("%d%",$F{CVAF}*100);

    # special case to handle FLT3 ITD
    if ($F{SYMBOL} eq 'FLT3' && $F{EXON} == 14 && length($F{ALT})>length($F{REF}) && length($F{ALT}) > 6){
	$var = "INSERTION VARIANT";
	$F{HGVSp_Short} = "ITD DETECTED (" . (length($F{ALT})-1) . " bp)";
	$F{CVAF} = '';
	$cat = "TIER 1-2" ;
    }
    $mutgenes{$F{SYMBOL}} = 1;

    $F{MAX_AF} = $F{MAX_AF} ? sprintf("%.3f\%",$F{MAX_AF} * 100) : 'none';
    push @tier12, join("\t",$F{SYMBOL},$F{EXON} ne '' ? "EXON " . $F{EXON} : "INTRON " . $F{INTRON},uc($var),$F{HGVSp_Short},$F{CVAF},"YES") if $cat eq 'TIER 1-2';
    push @tier12_full, join("\t",$F{SYMBOL},"chr".$F{CHROM}.":".$F{POS},$F{REF},$F{ALT},$F{HGVSc},$F{MAX_AF}) if $cat eq 'TIER 1-2';

    push @tier3, join("\t",$F{SYMBOL},$F{EXON} ne '' ? "EXON " . $F{EXON} : "INTRON " . $F{INTRON},uc($var),$F{HGVSp_Short},$F{CVAF},"UNKNOWN") if $cat eq 'TIER 3';
    push @tier3_full, join("\t",$F{SYMBOL},"chr".$F{CHROM}.":".$F{POS},$F{REF},$F{ALT},$F{HGVSc},$F{MAX_AF}) if $cat eq 'TIER 3';

}
close F;

print join("\t",qw(GENE REGION VARIANT), "PROTEIN CHANGE", "VAF" ,"CLINICAL SIGNIFICANCE"),"\n";
print join("\n",@tier12,@tier3),"\n\n";

my @wtgenes = grep { !defined($mutgenes{$_}) } sort keys %{$h->{GENECOV}};

print "No mutations were detected in the following sequenced genes and hotspots:\t",join(", ",@wtgenes),"\n\n";

my %g = ();
foreach my $i (@{$h->{PASSED_HOTSPOTQC}}){
  $i =~ /(\S+)_codon_([0-9,-]+)/;
  push @{$g{$1}}, $2;
}
my @out = ();
foreach my $j (sort keys %g){
    my $k = join(",",sort { $a <=> $b } @{$g{$j}});
#    $k =~ s/(?<!\d)(\d+)(?:,((??{$++1}))(?!\d))+/$1-$+/g; #s/(\d+)(?:,((??{$++1})))+/$1-$+/gx;                                                                                             
    push @out, "$j (" . ($k =~ /,|-/ ? "codons" : "codon") . " $k)";
}
@out = qw(NONE) if scalar @out == 0;
print "The following hotspots passed sequencing QC metrics:\t",join(", ",@out),"\n\n";

print "Variant detail:\n";
print join("\t","Gene","Location (hg19)","Reference allele","Variant allele","Transcript:Coding change","Population Frequency (%)"),"\n";
print join("\n",@tier12_full,@tier3_full),"\n\n";

print "The following target genes failed minimum sequencing QC metrics (>90% of positions with >50x coverage):\t",join(", ",map { $_ . " (" . sprintf("%.2f",$h->{FAILEDGENES}{$_}) . "% at 50x)" } sort keys %{$h->{FAILEDGENES}}),"\n\n";

#my %g = ();
#foreach my $i (@{$h->{FAILED_EXONS}}){
#  $i =~ /(\S+)_exon_(\d+)/;
#  push @{$g{$1}}, $2;
#}
#my @out = ();
#foreach my $j (sort keys %g){
#    my $k = join(",",sort { $a <=> $b } @{$g{$j}});
#    $k =~ s/(?<!\d)(\d+)(?:,((??{$++1}))(?!\d))+/$1-$+/g; #s/(\d+)(?:,((??{$++1})))+/$1-$+/gx;  
#    push @out, "$j (" . ($k =~ /,|-/ ? "exons" : "exon") . " $k)";
#}
#print "The following exons failed minimum sequencing QC metrics (>90% of positions with >50x coverage):\t",join(", ",@out),"\n\n";

open(R,$reportinfo) || die "Cant open report info file: $reportinfo";
while(<R>){
    print;
}
close R;

print "\n<", $REPORTVERSION, ">\n";
