#!/usr/bin/perl

use strict;
use Getopt::Long;
  
use File::Basename;
use Statistics::Basic qw(:all);
use JSON;

sub sum {
    my $s = 0;
    map { $s+=$_ } @_;
    $s;
}

sub print_qc {
    my ($id,$dat,$qcDat) = @_;

    my $qcflag = $$qcDat{PASSFLAG};
    my @out = ($id,$dat);
    if (defined($qcDat) && defined($qcDat->{$id})){
	my ($lo,$hi) = @{$qcDat->{$id}};

	if ($hi ne '.' && $lo ne '.'){
	    push @out, '['.join("-",$lo,$hi).']';
	    $qcflag = $qcDat->{FAILFLAG} if $dat > $hi || $dat < $lo;
	} elsif ($hi eq '.' && $lo ne '.'){
	    push @out, "[>$lo]";
	    $qcflag = $qcDat->{FAILFLAG} if $dat < $lo;
	} elsif ($lo eq '.' && $hi ne '.'){
	    push @out, "[<$hi]";
	    $qcflag = $qcDat->{FAILFLAG} if $dat > $hi;
	} else {
	    push @out, "[]";
	}
    } else {
	push @out, "[]";
    }
    return (@out,$qcflag);
}

sub print_gene_qc {
    my ($g,$h,$dat,$qcDat) = @_;

    my $qcflag = $$qcDat{PASSFLAG};
    my @out = ($g);
    foreach my $i (@$h){

	push @out, $dat->{$i};

	if (defined($qcDat) && defined($qcDat->{$i})){
	    my ($lo,$hi) = @{$qcDat->{$i}};
	    if ($hi ne '.' && $lo ne '.'){
		push @out, '['.join("-",$lo,$hi).']';
		$qcflag = $qcDat->{FAILFLAG} if $dat->{$i} > $hi || $dat->{$i} < $lo;
	    } elsif ($hi eq '.' && $lo ne '.'){
		push @out, "[>$lo]";
		$qcflag = $qcDat->{FAILFLAG} if $dat->{$i} < $lo;
	    } elsif ($lo eq '.' && $hi ne '.'){
		push @out, "[<$hi]";
		$qcflag = $qcDat->{FAILFLAG} if $dat->{$i} > $hi;
	    } else {
		push @out, "[]";
	    }
	} else {
	    push @out, "[]";
	} 
    }
    return (@out,$qcflag);
}

my $QCVERSION = basename($0,".pl");

my $qcfail = '***FAILED***';
my $qcpass = '';

my $REFFASTA='';
my $TARGETBED='';
my $AMPLICONBED='';
my $COVERAGEBED='';
my $DEMUXFILE='';
my $ABAM='';
my $CBAM='';
my $VCF='';
my $VARIANTFILE='';
my $HAPLOTECT='';
my $HAPLOTECTSITES='';
my $COVTHRESH = '50,100,200';
my $MINCOV1=50;

my $QCMETRICS = '';
my $INFOFILE = '';

my $BEDTOOLS="/usr/local/bin/bedtools";
my $SAMTOOLS="/usr/local/bin/samtools";

my $minBaseQual = 10;
my $minMapQual = 1;

my $help = '';

my $usage = '';

my $usage = <<END;

  HaloplexQC -r|reference <reference fasta> -q|qcmetrics <qc metrics file> -a|amplicons <amplicon bed file> -t|targets <target bed file> -b|bedfile -d|demux <demux info file> -c|consensusbam <consensus bam> -w|rawbam <raw aligned bam> -v|-variants <variant file> -i|infofile <assay info>

    Optional arguments:

      -t|bedtools [/usr/local/bin/samtools]
      -s|samtools [/usr/local/bin/bedtools]

      -h|help prints help

END

die $usage if $#ARGV < 2;

GetOptions("r|reference=s" => \$REFFASTA,
	   "a|amplicons=s" => \$AMPLICONBED,
           "t|targets=s" => \$TARGETBED,
           "b|bedfile=s" => \$COVERAGEBED,
	   "d|demux=s" => \$DEMUXFILE,
	   "c|consensusbam=s" => \$CBAM,
	   "v|variants=s" => \$VARIANTFILE,
	   "w|rawbam=s" => \$ABAM,
	   "q|qc=s" => \$QCMETRICS,
	   "i|infofile=s" => \$INFOFILE,
	   "m|haplotect=s" => \$HAPLOTECT,
	   "p|haplotectsites=s" => \$HAPLOTECTSITES,
	   "l|bedtools=s" => \$BEDTOOLS,
	   "s|samtools=s" => \$SAMTOOLS,
	   "h|help" => \$help);

die $usage if $help;

die "samtools location not valid: $SAMTOOLS\n" if !-s $SAMTOOLS;
die "bedtools location not valid: $BEDTOOLS\n" if !-s $BEDTOOLS;
die "amplicon bed file location not valid: $AMPLICONBED\n" if !-s $AMPLICONBED;
die "reference location not valid: $REFFASTA\n" if !-s $REFFASTA;
die "target bed file location not valid: $TARGETBED\n" if !-s $TARGETBED;
die "coverage qc bed file location not valid: $COVERAGEBED\n" if !-s $COVERAGEBED;
die "consensus bam not valid: $CBAM\n" if !-s $CBAM;
die "variant file not valid: $VARIANTFILE\n" if !-e $VARIANTFILE;
die "variant file not valid: $HAPLOTECT\n" if !-s $HAPLOTECT;
die "variant file not valid: $HAPLOTECTSITES\n" if !-e $HAPLOTECTSITES;
die "aligned bam not valid: $ABAM\n" if !-s $ABAM;
die "demux file location not valid: $DEMUXFILE\n" if !-s $DEMUXFILE;
die "qc file location not valid: $QCMETRICS\n" if !-s $QCMETRICS;
die "info file location not valid: $INFOFILE\n" if !-s $INFOFILE;

my @covThresholds = split(",",$COVTHRESH);
$MINCOV1 = $covThresholds[0];

# get QC metrics

my %QC = (PASSFLAG => $qcpass, FAILFLAG => $qcfail);
open(Q,$QCMETRICS) || die;
while(<Q>){
    chomp;
    next if /^#/;
    s/\s+$//g;
    my @l = split("\t",$_);
    $QC{$l[0]} = [ $l[1], $l[2] ];
}
close Q;

my $rg = `$SAMTOOLS view -H $CBAM`;

($rg =~ /ID:(\S+)/) || die "No instrument id found";
my $instrumentid = $1;

($rg =~ /PU:(\S+)/) || die "No lane id found";
my $laneid = $1;

($rg =~ /SM:(\S+)/) || die "No sample name found";
my $sample = $1;

($rg =~ /LB:(\S+)/) || die "No library id found";
my $library = $1;

my %out = (sample => $sample,
	   library => $library,
	   instrument_id => $instrumentid,
	   lane_id => $laneid,
	   version => $QCVERSION, 
	   lowcov => $MINCOV1,
	   coverageQCLevels => $COVTHRESH,
	   coverage_bed => $COVERAGEBED,
	   target_bed => $TARGETBED,
	   amplicon_bed => $AMPLICONBED,
	   reference => $REFFASTA);

# read in demux file and get demux stats
my $flowcellreads = 0;
my $unknownreads = 0;
my $nsamples = 0;
open(D,$DEMUXFILE) || die;
while(<D>){
    chomp;
    
    die "demux file $DEMUXFILE not in proper format. Must be:\nname\tindex\tfastq1\tfastq2\tnumreads\n" if $_ !~ /\S+\t\S+\t\S+\t\S+\d+/;
    
    my @l = split("\t",$_);
    $flowcellreads += $l[4];
    if(/unknown/){
	$unknownreads = $l[4];
    } else {
	$nsamples++;
    }
}

$out{FLOWCELL_READS} = $flowcellreads;
$out{FLOWCELL_SAMPLES} = $nsamples;
$out{UNKNOWN_READS} = $unknownreads;

# run samtools to get mapping stats
my $TOTAL_READS=`$SAMTOOLS view -c -F 0x900 $ABAM`; # no secondary or supp alignments
chomp $TOTAL_READS;
my $ALIGNED_READS=`$SAMTOOLS view -c -F 0x904 $ABAM`; # no secondary or supp alignments or unmapped reads
chomp $ALIGNED_READS;
my $ONTARGET_READS=`$SAMTOOLS view -c -F 0x904 -L $TARGETBED $ABAM`; # overlap target bed file
chomp $ONTARGET_READS;
my $CONDENSED_ONTARGET_READS=`$SAMTOOLS view -c -F 0x904 -L $TARGETBED $CBAM`; # overlap target bed file in condensed bam
chomp $CONDENSED_ONTARGET_READS;

$out{TOTAL_READS} = $TOTAL_READS;
$out{ALIGNED_READS} = $ALIGNED_READS;
$out{PERCENT_ALIGNED_READS} = sprintf("%.1f",$ALIGNED_READS / $TOTAL_READS * 100);
$out{ONTARGET_READS} = $ONTARGET_READS;
$out{PERCENT_ONTARGET_READS} = sprintf("%.1f",$ONTARGET_READS / $ALIGNED_READS * 100);
$out{CONDENSED_ONTARGET_READS} = $CONDENSED_ONTARGET_READS;
$out{PERCENT_CONDENSED_ONTARGET_READS} = sprintf("%.1f",$CONDENSED_ONTARGET_READS / $ALIGNED_READS * 100);

# reads per umi
open(BC,"$SAMTOOLS view -F 0x904 -L $TARGETBED $CBAM |") || die;
my @h = ();
while(<BC>){
    $_ =~ /BC:Z:(\d+)/;
    push @h, $1;
}
close BC;
my %bc = ();
map { $bc{$_}++ } @h;

$out{MEAN_READS_PER_UMI} = sprintf("%.2f", mean(@h));
$out{MEDIAN_READS_PER_UMI} = sprintf("%.1f",median(@h));
$out{PERCENT_UMI_COVERAGE_1x} = sprintf("%.1f",$bc{1} / sum(values %bc) * 100);
$out{PERCENT_UMI_COVERAGE_2x} = sprintf("%.1f",$bc{2} / sum(values %bc) * 100);
$out{PERCENT_UMI_COVERAGE_3x} = sprintf("%.1f",$bc{3} / sum(values %bc) * 100);
$out{PERCENT_UMI_COVERAGE_OVER_4x} = sprintf("%.1f",(1 - (($bc{1}+$bc{2}+$bc{3}) / sum(values %bc))) * 100);

my %amps = ();
open(A,$AMPLICONBED) || die "cant open file: $AMPLICONBED";
while(<A>){
    chomp;
    my @l = split("\t",$_);
    $amps{$l[3]} = 0;
}
close A;

# count of proper read pairs per amplicon
open(AMPS,"$SAMTOOLS view -F 0x904 -f 0x42 -L $TARGETBED $CBAM |") || die;
while(<AMPS>){
    $_ =~ /X1:Z:(\S+)/;
    $amps{$1}++;
}
close AMPS;

my %ampcnt = ();
map { $ampcnt{$_}++ } values %amps;

$out{MEAN_READ_PAIRS_PER_AMPLICON} = sprintf("%.2f",mean(values %amps));
$out{MEDIAN_READ_PAIRS_PER_AMPLICON} = sprintf("%.1f",median(values %amps));
$out{PERCENT_AMPLICON_COVERAGE_0x} = sprintf("%.1f",$ampcnt{0} / scalar(keys %amps) * 100);
$out{PERCENT_AMPLICON_COVERAGE_10x} = sprintf("%.1f",(1 - (sum(map { $ampcnt{$_} } 0..9) / scalar(keys %amps))) * 100);
$out{PERCENT_AMPLICON_COVERAGE_50x} = sprintf("%.1f",(1 - (sum(map { $ampcnt{$_} } 0..49) / scalar(keys %amps))) * 100);
$out{PERCENT_AMPLICON_COVERAGE_100x} = sprintf("%.1f",(1 - (sum(map { $ampcnt{$_} } 0..99) / scalar(keys %amps))) * 100);

# get consensus coverage depth for each position in the CoverageQC bed file (that is, the coding portions that will be reported)
my %cov = ();
my %hotspot = ();
my @consensusCoverage = ();
my @rawCoverage = ();
my $totalSize =0;
my %coverageLevels = ();
open(C,"$SAMTOOLS depth -d 100000 -a -q $minBaseQual -Q $minMapQual -b $COVERAGEBED $CBAM $ABAM | /usr/bin/awk -v OFS=\"\t\" '{ print \$1,\$2-1,\$2,\$3,\$4; }' | $BEDTOOLS intersect -a stdin -b $COVERAGEBED -wo |") || die "error running QC script $0";
while(<C>){
    chomp;
    my @l = split("\t",$_);
    if (/HOTSPOTQC/){
	push @{$hotspot{$l[9] . '_codon_' . $l[11]}}, $l[3];
	
    } else {

	$totalSize++;
	push @{$cov{$l[9]}{consensusCoverage}}, $l[3];
	push @{$cov{$l[9]}{rawCoverage}}, $l[4];
	$cov{$l[9]}{geneSize}++;
	map { $cov{$l[9]}{passedCoverage}{$_}++ if $l[3] >= $_ } @covThresholds;

	push @{$cov{$l[9]}{exons}{$l[10]}{consensusCoverage}}, $l[3];
	push @{$cov{$l[9]}{exons}{$l[10]}{rawCoverage}}, $l[4];
	$cov{$l[9]}{exons}{$l[10]}{exonSize} = $l[7] - $l[6];
	map { $cov{$l[9]}{exons}{$l[10]}{passedCoverage}{$_}++ if $l[3] >= $_ } @covThresholds;

    }
}
close C;

# assay wide coverage metrics
$out{MEAN_UNIQUE_COVERAGE} = sprintf("%.1f",mean(map { @{$cov{$_}{consensusCoverage}} } keys %cov));
$out{MEDIAN_UNIQUE_COVERAGE} = sprintf("%.1f",median(map { @{$cov{$_}{consensusCoverage}} } keys %cov));
foreach my $i (@covThresholds){
    $out{"PERCENT_UNIQUE_COVERAGE_" . $i . "x"} = sprintf("%.1f",sum(map { $cov{$_}{passedCoverage}{$i} } keys %cov) / $totalSize * 100);
}
$out{MEAN_RAW_COVERAGE} = sprintf("%.1f",mean(@rawCoverage));
$out{MEDIAN_RAW_COVERAGE} = sprintf("%.1f",mean(@rawCoverage));

my @failedexons = ();

# gene level coverage: mean and median for consensus and raw bam
foreach my $g (sort {$a cmp $b} keys %cov){
    my $con_med = sprintf("%.1f",median(@{$cov{$g}{consensusCoverage}}));
    my $raw_med = sprintf("%.1f",median(@{$cov{$g}{rawCoverage}}));
    my $con_mean = sprintf("%.1f",mean(@{$cov{$g}{consensusCoverage}}));
    my $raw_mean = sprintf("%.1f",mean(@{$cov{$g}{rawCoverage}}));

    $out{GENECOV}{$g} = { MEAN_RAW_COVERAGE => $raw_mean, MEAN_TARGET_COVERAGE => $con_mean,
			  MEDIAN_RAW_COVERAGE => $raw_med, MEDIAN_TARGET_COVERAGE => $con_med, 
			  SIZE => $cov{$g}{geneSize} };
    map { $out{GENECOV}{$g}{"PERCENT_TARGET_COVERAGE_" . $_ . "x"} = sprintf("%.1f",$cov{$g}{passedCoverage}{$_} / $cov{$g}{geneSize} * 100) } @covThresholds;

    $out{FAILEDGENES}{$g} = $cov{$g}{passedCoverage}{$covThresholds[0]} / $cov{$g}{geneSize} * 100 if ($cov{$g}{passedCoverage}{$covThresholds[0]} / $cov{$g}{geneSize} * 100 < $QC{PERCENT_TARGET_COVERAGE_50x}->[0]);
    
    $out{GENECOV}{$g}{FAILED_EXONS} = [];
    foreach my $e (sort keys %{$cov{$g}{exons}}){
	my $con_med = sprintf("%d",median(@{$cov{$g}{exons}{$e}{consensusCoverage}}));
	my $raw_med = sprintf("%d",median(@{$cov{$g}{exons}{$e}{rawCoverage}}));
	my $con_mean = sprintf("%.1f",mean(@{$cov{$g}{exons}{$e}{consensusCoverage}}));
	my $raw_mean = sprintf("%.1f",mean(@{$cov{$g}{exons}{$e}{rawCoverage}}));
	$out{GENECOV}{$g}{exons}{$e} = { MEAN_RAW_COVERAGE => $raw_mean, MEAN_TARGET_COVERAGE => $con_mean,
					 MEDIAN_RAW_COVERAGE => $raw_med, MEDIAN_TARGET_COVERAGE => $con_med,
					 SIZE => $cov{$g}{exons}{$e}{exonSize} };
	map { $out{GENECOV}{$g}{exons}{$e}{"PERCENT_TARGET_COVERAGE_" . $_ . "x"} = sprintf("%.1f",$cov{$g}{exons}{$e}{passedCoverage}{$_} / $cov{$g}{exons}{$e}{exonSize} * 100) } @covThresholds;

	push @{$out{GENECOV}{$g}{FAILED_EXONS}}, $e if $cov{$g}{exons}{$e}{passedCoverage}{$MINCOV1} / $cov{$g}{exons}{$e}{exonSize} * 100 < $QC{PERCENT_TARGET_COVERAGE_50x}->[0];
    }
    $out{GENECOV}{$g}{FAILED_EXON_COUNT} = scalar @{$out{GENECOV}{$g}{FAILED_EXONS}};
}

# hotspot QC
my @passed = ();
my @failed = ();
foreach my $i (sort keys %hotspot){
    if ((sort {$a<=>$b} @{$hotspot{$i}})[0] < $MINCOV1){
	push @failed, $i;
    } else {
	push @passed, $i;
    }
}

$out{FAILED_HOTSPOTQC} = \@failed;
$out{PASSED_HOTSPOTQC} = \@passed;
$out{FAILED_EXONS} = [ map { @{$out{GENECOV}{$_}{FAILED_EXONS}} } sort keys %{$out{GENECOV}} ];

# get variant information 
my @data;

format qc_format =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     @<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<     @<<<<<<<<<<<<<
@data
.

format gene_qc_format_header =
@<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<
@data
.

format gene_qc_format =
@<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<< @<<<<<<<<<<
@data
.

$out{timestamp} = scalar localtime();


# print QC report
open(F,">$library.qc.txt") || die;
select(F);
$~="qc_format";

print F "SAMPLE ID: $out{sample}\tLIBRARY ID: $out{library}\nLANE ID: $out{lane_id}\tINSTRUMENT ID: $out{instrument_id}\tDATE REPORTED: ", $out{timestamp}, "\n\n";

print F "FLOWCELL QC:\n\n";
map { @data = print_qc($_,$out{$_},\%QC); write; } qw(FLOWCELL_READS FLOWCELL_SAMPLES UNKNOWN_READS);
print F "\n";

print F "SEQUENCE DATA QC:\n\n";
map { @data = print_qc($_,$out{$_},\%QC); write; } qw(TOTAL_READS ALIGNED_READS PERCENT_ALIGNED_READS ONTARGET_READS PERCENT_ONTARGET_READS CONDENSED_ONTARGET_READS PERCENT_CONDENSED_ONTARGET_READS);
print F "\n";

print F "LIBRARY QC:\n\n";
map { @data = print_qc($_,$out{$_},\%QC); write; } qw(MEAN_READS_PER_UMI MEDIAN_READS_PER_UMI PERCENT_UMI_COVERAGE_1x PERCENT_UMI_COVERAGE_2x PERCENT_UMI_COVERAGE_3x PERCENT_UMI_COVERAGE_OVER_4x);
print F "\n";

map { @data = print_qc($_,$out{$_},\%QC); write; } qw(MEAN_READ_PAIRS_PER_AMPLICON MEDIAN_READ_PAIRS_PER_AMPLICON PERCENT_AMPLICON_COVERAGE_0x PERCENT_AMPLICON_COVERAGE_10x PERCENT_AMPLICON_COVERAGE_50x PERCENT_AMPLICON_COVERAGE_100x);
print F "\n";

print F "HOTSPOT QC:\n\n";
@data = ("FAILED_HOTSPOTQC",scalar @{$out{FAILED_HOTSPOTQC}},"[0]", scalar @{$out{FAILED_HOTSPOTQC}} > $QC{FAILED_HOTSPOTQC}[0] ? $qcfail : $qcpass);
write;

@data = ("PASSED_HOTSPOTQC",scalar @{$out{PASSED_HOTSPOTQC}},"[".(scalar keys %hotspot). "]",scalar @{$out{PASSED_HOTSPOTQC}} < $QC{PASSED_HOTSPOTQC}[0] ? $qcfail : $qcpass);
write;
print F "\n";

if (scalar @{$out{FAILED_HOTSPOTQC}} > 0){
  print F "**FAILED_HOTSPOTS**\t ". join(",",@{$out{FAILED_HOTSPOTQC}}),"\n\n";
}

print F "ASSAY-WIDE COVERAGE METRICS:\n\n";
map { @data = print_qc($_,$out{$_},\%QC); write; } (qw(MEAN_UNIQUE_COVERAGE MEDIAN_UNIQUE_COVERAGE),map { "PERCENT_UNIQUE_COVERAGE_" . $_ . "x" } @covThresholds);
print F "\n";

print F "GENE/TARGET LEVEL COVERAGE QC:\n\n";

my @headers = (qw(MEAN_TARGET_COVERAGE MEDIAN_TARGET_COVERAGE),(map { "PERCENT_TARGET_COVERAGE_" . $_ . "x" } @covThresholds),"FAILED_EXON_COUNT");
$~="gene_qc_format_header";
@data = ("GENE",@headers);
write;

$~="gene_qc_format";
foreach my $g (sort keys %{$out{GENECOV}}){
    @data = print_gene_qc($g,\@headers,$out{GENECOV}{$g},\%QC);
    write;
}
print F "\n";

print F "FAILED EXONS:\n\n";
@headers = (qw(MEAN_TARGET_COVERAGE MEDIAN_TARGET_COVERAGE),(map { "PERCENT_TARGET_COVERAGE_" . $_ . "x" } @covThresholds));
$~="gene_qc_format_header";
@data = ("EXON",@headers,"");
write;

$~="gene_qc_format";
foreach my $g (sort keys %{$out{GENECOV}}){
    foreach my $e (sort @{$out{GENECOV}{$g}{FAILED_EXONS}}){
	@data = (print_gene_qc($e,\@headers,$out{GENECOV}{$g}{exons}{$e},\%QC),"","","","","");
        write;
    }
}
print F "\n";


my %aa3to1 = qw(Ala A Arg R Asn N Asp D Asx B Cys C Glu E Gln Q Glx Z Gly G His H Ile I Leu L Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Xxx X Ter *);

my %vars = ("TIER 1-2" => [], "TIER 3" => [], "TIER 4" => [], "FILTERED" => []);
open(VF,$VARIANTFILE) || die;
$_ = <VF>;
chomp;
s/\s+$//;
s/\S+\.CVAF/VAF/;
s/\S+\.TAMP/AMPLICONS/;
s/\S+\.SAMP/SUPPORTING_AMPLICONS/;
s/\S+\.NR/COVERAGE/;
s/\S+\.NV/VARIANT_READS/;

@headers = split("\t",$_);

while(<VF>){
    chomp;
    s/\s+$//;
    my @F = split("\t",$_);
    my %F = map { $headers[$_] => $F[$_] } 0..$#F;

    $F{Consequence} = (split("&",$F{Consequence}))[0];
    $F{Consequence} =~ /^(\S+?)_/;
    $F{EXON} =~ s/\/\d+//g;
    $F{INTRON} =~ s/\/\d+//g;

    my $HGVSpShort = $F{HGVSp};
    $HGVSpShort =~ s/\S+:p\.(\S+?)/\1/g;

    while( $HGVSpShort and my ( $find, $replace ) = each %aa3to1 ) {
        eval "\$HGVSpShort =~ s{$find}{$replace}g";
    }
    $F{HGVSp} = $HGVSpShort;

    my $cat = ($F{MYELOSEQ_MDS_AC}>0 || $F{MYELOSEQ_TCGA_AC}>0) ? "TIER 1-2" : "TIER 3";
    $cat = "TIER 4" if $F{Consequence}=~/synonymous/;

    # special case to handle FLT3 ITD
    if ($F{SYMBOL} eq 'FLT3' && $F{EXON} == 14 && length($F{ALT})>length($F{REF}) && length($F{ALT}) > 6){
        $cat = "TIER 1-2" ;
    }

    $cat = "FILTERED" if $F{FILTER} ne 'PASS';

    $F{MAX_AF} = $F{MAX_AF} ? sprintf("%.3f\%",$F{MAX_AF} * 100) : 'none';

    push @{$vars{$cat}}, [$cat,$F{FILTER},$F{set},$F{SYMBOL},"chr".$F{CHROM},$F{POS},$F{REF},$F{ALT},$F{Consequence},$F{HGVSp},$F{HGVSc},$F{EXON} || 'NA',$F{INTRON} || 'NA',$F{MAX_AF},$F{COVERAGE},$F{VARIANT_READS},sprintf("%.2f\%",$F{VAF} * 100),$F{AMPLICONS},$F{SUPPORTING_AMPLICONS}];

}
close VF;

print F "TIER 1-3 Variant detail:\n\n";
if (scalar @{$vars{"TIER 1-2"}} > 0 || scalar @{$vars{"TIER 3"}} > 0){
   print F uc(join("\t", qw(category filter callers gene chrom position ref alt consequence p_syntax c_syntax exon intron pop_freq coverage variant_reads vaf amplicons supporting_amplicons))), "\n";
   print F join("\n", map { join("\t", @{$_}) } (@{$vars{"TIER 1-2"}},@{$vars{"TIER 3"}})), "\n\n";
} else {
 print F "NO VARIANTS\n\n";
}

print F "TIER 4 Variants:\n\n";
if (scalar @{$vars{"TIER 4"}} > 0){
   print F uc(join("\t", qw(category filter callers gene chrom position ref alt consequence p_syntax c_syntax exon intron pop_freq coverage variant_reads vaf amplicons supporting_amplicons))), "\n";
   print F join("\n", map { join("\t", @{$_}) } @{$vars{"TIER 4"}}), "\n\n";
} else {
 print F "NO VARIANTS\n\n";
}

print F "Filtered Variants:\n\n";
if (scalar @{$vars{"FILTERED"}} > 0){
   print F uc(join("\t", qw(category filter callers gene chrom position ref alt consequence p_syntax c_syntax exon intron pop_freq coverage variant_reads vaf amplicons supporting_amplicons))),"\n";
   print F join("\n", map { join("\t", @{$_}) } @{$vars{"FILTERED"}}) . "\n\n";
} else {
 print F "NO VARIANTS\n\n";
}

print F "CONTAMINATION ESTIMATE\n\n";
open(H,$HAPLOTECT) || die "Cant open haplotect file: $HAPLOTECT";
$_ = <H>;
print F $_;
$_ = <H>;
chomp;
close H;
my @h = split("\t",$_);
print F join("\t",(@h,($h[2] >= $QC{HAPLOTECT_SITES} && $h[6] > $QC{HAPLOTECT_CONTAM} ? $qcfail : $qcpass))),"\n";
$out{HAPLOTECT_SITES} = $h[2];
$out{HAPLOTECT_CONTAM} = $h[6];
$out{HAPLOTECT_CONTAMSNP} = $h[7];

print F "\n\n";

print F "CONTAMINATING HAPLOTYPE INFORMATION\n\n";
open(H,$HAPLOTECTSITES) || die "Cant open haplotect sites file: $HAPLOTECTSITES";
while(<H>){
  print F $_;
}
close H;

print F "\n\n";

open(R,$INFOFILE) || die "Cant open report info file: $INFOFILE";
while(<R>){
    print F $_;
}
close R;
print F "\n";

print F "#qc version: $out{version} mincov value: $out{lowcov}\n";
print F "#coverage bed file: $out{coverage_bed}\n#target bed file: $out{target_bed}\n#amplicon bed file: $out{amplicon_bed}\n";

close F;


# print QC json file
open(F,">$library.qc.json") || die;
print F encode_json \%out, "\n";
close F;
