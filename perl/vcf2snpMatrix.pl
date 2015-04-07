#!/usr/bin/env perl
#
# 
#
# for pilon vcfs processed to remove passing ref bases or full vcfs (ambiguous sites need to be present)
# note - currently only works for single sequence reference
# includes transfer of vcf annotations
# copied from ~cdesjard/plscripts/snp_tools/pilonVCF2SNPmatrixII.pl

use strict;
use warnings;
use Data::Dumper;

my $usage = "usage: $0 file_of_vcf_filenames binary_a0/binary_a1/binary_a- [snps_only] > out_file\n";

# base = print actual base
# binary_a0 = print 0 ref, 0 amb, 1 alt
# binary_a1 = print 0 ref, 1 amb, 1 alt
# binary_a1 = print 0 ref, - amb, 1 alt

my $infilename = $ARGV[0] or die $usage;
my $type = $ARGV[1] or die $usage;
chomp $type;
my $events = 'all';
if ($ARGV[2]) {
	if ($ARGV[2] eq 'snps_only') {
		$events = 'snps';
	} else {
		die "$ARGV[2] is an invalid argument.\n";
	}
}

if ($type ne "binary_a0" and $type ne "binary_a1" and $type ne "binary_a-") {
	die "Error, $type is not a valid option for ambiguity coding.\n";
}

unless (open (INFILE, $infilename)) {
    die "File $infilename not found.\n";
}

my %ref_bases;
my %alt_bases;
my %passed_snp_positions;
my %all_snp_positions;
my %genome_list;
my %annotations;

foreach my $in (<INFILE>) {
    chomp $in;
    if ($in =~ /(\w+)(\.annotated)*\.vcf$/) {
        my $genome = $1;
        
        unless (open (VCFFILE, $in)) {
            die "File $in not found.\n";
        }
        
        print STDERR "Searching $in...\n";

        $genome_list{$genome} = 1; 
        
        foreach my $vcfline (<VCFFILE>) {
            chomp $vcfline;
            if ($vcfline =~ /^\#/) {
            } else {
                my @x = split (/\t/, $vcfline);
                if (($x[4] eq '.' and $x[6] eq 'PASS') or $x[4] =~ /^<\w+>$/) {
                } elsif (($x[6] eq 'PASS' and $events eq 'all') or ($x[3] =~ /^[ACGT]$/ and $x[4] =~ /^[ACGT]$/ and $x[6] eq 'PASS' and $events eq 'snps')) {
                        $alt_bases{$genome}{$x[1]} = $x[4];
                        $ref_bases{$x[1]} = $x[3];
                        $passed_snp_positions{$x[1]} = 1;
                        if ($x[10]) {
                            $annotations{$x[1]} = $x[10];
                        }
                } else {
                    $all_snp_positions{$genome}{$x[1]} = 1;
                }
            }           
        }       
    }
}
close INFILE;

#print Dumper(\%genome_list);

my @positions = sort {$a <=> $b} keys %passed_snp_positions;
my @genome_list = sort keys %genome_list;

print "#snp_pos";
foreach my $position (@positions) {
    print "\t$position";
}
print "\n";

print "#annotation";
foreach my $position(@positions) {
    if ($annotations{$position}) {
        print "\t$annotations{$position}";
    } else {
        print "\t.";
    }
}
print "\n";

print "reference";
foreach my $position(@positions) {
    if ($type eq 'base') {
        print "\t$ref_bases{$position}";
    } elsif ($type =~ /binary/) {
        print "\t0";
    }
}
print "\n";

foreach my $genome (@genome_list) {
    print "$genome";
    foreach my $position(@positions) {
        if ($alt_bases{$genome}{$position}) {
            print "\t1";
        } elsif ($all_snp_positions{$genome}{$position}) {
            if ($type eq 'binary_a1') {
                print "\t1";
            } elsif ($type eq 'binary_a0') {
                print "\t0";
            } elsif ($type eq 'binary_a-') {
                print "\t-";
            }         
        } else {
            print "\t0";                         
        }
    }
    print "\n";
}



exit;

####


