#! /usr/bin/perl

#Copyright (C) 2015 Feiyu Du <fdu@genome.wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;
use File::Spec;
use File::Basename;
use IO::File;

umask 002;

die "Provide batch_name, index_mapping_tsv, instrument_data_id, flow_cell_id, lane_number, csf_raw_fastq_dir" unless @ARGV and @ARGV == 6;
my ($name, $map_tsv, $instrument_id, $flow_cell_id, $lane_num, $fastq_dir) = @ARGV;

unless (-d $fastq_dir) {
    die "csq_fastq_dir $fastq_dir is not valid";
}

unless ($lane_num =~ /^\d+$/) {
    die "lane_number $lane_num is not a number";
}

my $dir = '/gscmnt/gc13016/cle/54f8f7b915cb472aa183c721307369ab_scratch_space/myeloseq/RUN';
my $out_dir = $dir . '/wdl_out';
my $log_dir = $dir . '/log';
my $wdl_path = $dir . '/Haloplex.wdl';
my $template_json = $dir.'/template.json';
my $conf_path = $dir . '/application.conf';

my $run_dir = $out_dir . "/$name";
unless (mkdir $run_dir) {
    die "Failed to make directory $run_dir";
}

my $pu = join '.', $flow_cell_id, $lane_num;

my $sample_index  = $run_dir . '/sample_index';
my $json_file     = $run_dir . '/Haloplex.json';

my $in_fh    = IO::File->new($map_tsv);
my $index_fh = IO::File->new(">$sample_index");

while (my $line = $in_fh->getline) {
    chomp $line;
    my ($index, $lib_name) = split /\t/, $line;
    my ($sample) = $lib_name =~ /^(\S+?)\-lib\d+/;
    my $rg_str = '@RG\tPU:'.$pu.'\tLB:'.$lib_name.'\tID:'.$instrument_id.'\tSM:'.$sample.'\tPL:ILLUMINA\tCN:WUGSC';
    $index_fh->print(join "\t", $index, $lib_name, $rg_str);
    $index_fh->print("\n");
}
$in_fh->close;
$index_fh->close;

my $template_json_fh = IO::File->new($template_json);
my $json_fh  = IO::File->new(">$json_file");

while (my $l = $template_json_fh->getline) {
    if ($l =~ /SampleSheet/) {
        $l =~ s/PATH/$sample_index/;
    }
    elsif ($l =~ /OutputDir/) {
        $l =~ s/PATH/$run_dir/;
    }
    elsif ($l =~ /IlluminaDir/) {
        $l =~ s/PATH/$fastq_dir/;
    }
    $json_fh->print($l);
}
$template_json_fh->close;
$json_fh->close;

my $out_log = File::Spec->join($dir, 'log', $name.'.out');
my $err_log = File::Spec->join($dir, 'log', $name.'.err');

`bsub -oo $out_log -eo $err_log -q compute-mgi-cle -R "select[mem>8000] rusage[mem=8000]" -M 8000000 -a "docker(sleongmgi/cromwell:develop-with-mysql)" /usr/bin/java -Dconfig.file=$conf_path -jar /cromwell/cromwell.jar run $wdl_path $json_file`;
    
print "$name job submitted\n";

