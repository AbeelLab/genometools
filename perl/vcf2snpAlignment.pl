#!/usr/bin/env perl
#
# 
#
# for pilon vcfs (reduced or not reduced)
# copied from ~cdesjard/plscripts/snp_tools/plscripts/snp_tools/pilonVCF2fastaII.pl

use strict;
use warnings;
use Data::Dumper;

my $usage = "usage: $0 file_of_vcf_filenames [positions_to_exclude] > out_file\n";

# vcf filenames must end in .vcf

my $infilename = $ARGV[0] or die $usage;
chomp $infilename;

my %sites_to_exclude;
if ($ARGV[1]) {
	my $exfile = $ARGV[1];

	unless (open (EXFILE, $exfile)) {
    	die "File $exfile not found.\n";
	}

	foreach my $in (<EXFILE>) {
    	if ($in =~ /^#/) {
    	} else {
        	chomp $in;
        	my @x = split(/\t/, $in);
        	$sites_to_exclude{$x[0]} = 1;
        }
    } 
    close EXFILE;
}

#print Dumper(\%sites_to_exclude);

unless (open (INFILE, $infilename)) {
    die "File $infilename not found.\n";
}

my %ref_bases;
my %alt_bases;
my %passed_snp_positions;
my %all_snp_positions;
my %genome_list;

foreach my $in (<INFILE>) {
    chomp $in;
    if ($in =~ /([\w-]+)(\.annotated)*\.vcf$/) { 
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
                if (($x[4] eq '.' and $x[6] eq 'PASS') or $sites_to_exclude{$x[1]}) {
                	#### print STDERR "excluded site $x[1]\n"; ####
                } elsif ($x[3] =~ /^[ACGT]$/ and $x[4] =~ /^[ACGT]$/ and $x[6] eq 'PASS') {
                    $alt_bases{$genome}{$x[1]} = $x[4];
                    $ref_bases{$x[1]} = $x[3];
                    $passed_snp_positions{$x[1]} = 1;
                } else {
                	$all_snp_positions{$genome}{$x[1]} = 1;
                }
            }           
        }       
    }
}
close INFILE;

#print Dumper(\%all_snp_positions);

my @positions = sort {$a <=> $b} keys %passed_snp_positions;
my @genome_list = sort keys %genome_list;

print ">reference\n";
my $sequence = '';
foreach my $position(@positions) {
    $sequence .= $ref_bases{$position};
}
print_fasta($sequence, 60);

foreach my $genome(@genome_list) {
    print ">$genome\n";
    my $sequence = '';
    foreach my $position(@positions) {
        if ($alt_bases{$genome}{$position}) {
            $sequence .= $alt_bases{$genome}{$position};
        } elsif ($all_snp_positions{$genome}{$position}) {
            $sequence .= 'N';
        } else {
            $sequence .= $ref_bases{$position};
        }
    }
    print_fasta($sequence, 60);
}

exit;

####

sub print_fasta {
    my($seq, $length) = @_;
    for (my $i = 0; $i < length($seq); $i += $length) {
        print substr($seq, $i, $length), "\n";
    }
}

exit;

