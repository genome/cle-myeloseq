#! /usr/bin/perl

#Copyright (C) 2018 Feiyu Du <fdu@wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;
use File::Spec;
use File::Basename;
use File::Copy qw(copy);

umask 002;

unless (@ARGV >= 1) {
    die "Batch directory to be transferred are required. Sample case dirs that are excluded for transfer are optional"; 
}

my $batch_dir = shift @ARGV;
my @exclude_cases = @ARGV;
my $batch_name = basename $batch_dir;

unless (-d $batch_dir) {
    die "The provided batch dir $batch_dir is not valid";
}

my $dir     = '/gscmnt/gc13016/cle/54f8f7b915cb472aa183c721307369ab_scratch_space/myeloseq/RUN';
my $log_dir = $dir.'/log';
my $tmp_dir = $dir.'/xfer_tmp/'.$batch_name;

if (-d $tmp_dir) {
    die "Batch tmp dir: $tmp_dir already existing. Remove or rename it before retry";
}
else {
    mkdir $tmp_dir;
    unless (-d $tmp_dir) {
        die "Failed to mkdir $tmp_dir";
    }
}

my @files_xfer = qw(
    annotated.vcf.gz  annotated_filtered.vcf.gz 
    consensus.bam     consensus.bam.bai 
    variant_report.txt   igv.xml  qc.txt
);

opendir(my $dir_h, $batch_dir);
for my $case (readdir $dir_h) {
    next if $case =~ /^\./;
    next if grep {$_ eq $case} @exclude_cases;
    next unless -d File::Spec->join($batch_dir, $case);

    my ($real_name) = $case =~ /^(\S+lib\d+)_[ATCG]{8}/;
    #my $xfer_file_str = join ' ', map{File::Spec->join($batch_dir, $case, $real_name.'.'.$_)}@files_xfer;
    my $local_tmp_dir = File::Spec->join($tmp_dir, $case);
    mkdir $local_tmp_dir;

    for my $file_name (@files_xfer) {
        my $file_path = File::Spec->join($batch_dir, $case, $real_name.'.'.$file_name);
        copy $file_path, $local_tmp_dir;
    }
}

closedir $dir_h;

my $dest_dir = $batch_name.'/';
my $err_log  = File::Spec->join($log_dir, $batch_name.'_xfer.err');
my $out_log  = File::Spec->join($log_dir, $batch_name.'_xfer.out');

my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
my $private_key = "/gscuser/$username/.ssh/id_rsa_";
my %info = (
   fdu => {
       suffix => 'pathology',
       user   => 'myloseq_transfer',
   },
   mharriso => { 
       suffix => 'myeloseq',
       user   => 'myloseq_transfer2',
   }
);
$private_key .= $info{$username}->{suffix};

unless (-s $private_key) {
    die "Private key: $private_key is not valid";
}

my $host  = $info{$username}->{user}.'@128.252.17.197';
my $image = 'registry.gsc.wustl.edu/dataxfer/data-transfer-helper';
my $queue = 'compute-mgi-cle';

my $cmd = "bsub -e $err_log -o $out_log -q $queue -a 'docker($image)' scp -i $private_key -r $tmp_dir $host:$dest_dir";
system $cmd;
