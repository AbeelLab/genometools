#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");

use Fasta_reader;
use Gene_obj;
use GFF3_utils;
use Nuc_translator;

use Data::Dumper;

my $debug = 0;


#DB::disable_profile();  # Devel::NYTProf

my $usage = <<__EOUSAGE__;

##############################################################################
#
#  --gff3          gene annotations in gff3 format
#  --genome        genome sequence in fasta format.
#  --vcf           SNP data in vcf format
#  --codon-based   the effect of substitutions will be evaluated in the context of all
#                  substitutions in the codon. 
#                  If this parameter is not given the effect of each substitution 
#                  will be evaluated independently, disregarding other substitutions in
#                  the same codon.
#
#  -X              write the .coding_mutations_described.txt file
#
##############################################################################


__EOUSAGE__

;

my ($gff3_file, $genome_file, $vcf_file );
my $codon_based = 0;

my $WRITE_CODING_MUTANT_DESCRIPTIONS = 0;

&GetOptions ( 'gff3=s' => \$gff3_file,
	'genome=s' => \$genome_file,
	'vcf=s' => \$vcf_file,
	'codon-based' => \$codon_based,
	'X' => \$WRITE_CODING_MUTANT_DESCRIPTIONS,
	);

unless ($gff3_file && $genome_file && $vcf_file) {
	die $usage;
}

foreach my $file ($gff3_file, $genome_file, $vcf_file) {
	unless (-s $file) {
		die "Error, cannot find $file";
	}
}




my %CDS_SEQUENCES;  # storing all CDS sequences here... very hacky, not proud.

my $SNP_LINE_COUNTER = 0;

main: {
	
	
	# get gene annotations.
	my $gene_obj_indexer = {};
	my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer);
	
	my $fasta_reader = new Fasta_reader($genome_file);
	my %genome = $fasta_reader->retrieve_all_seqs_hash();
	
	my $curr_contig = "";
	my $snp_text = "";
	
	my $coding_mutations_explained_file = "$vcf_file.coding_mutations_described.txt";
	my $cod_mut_exp_fh;
	if ($WRITE_CODING_MUTANT_DESCRIPTIONS) {
		open ( $cod_mut_exp_fh, ">$coding_mutations_explained_file") or die "Error, cannot write to $coding_mutations_explained_file";
		
	}
	
	
	#DB::enable_profile();
	
	my $counter = 0;
	print STDERR "-reading vcf file: $vcf_file\n";
	
	
	open (my $fh, $vcf_file) or die "Error, cannot open file $vcf_file";
	while (<$fh>) {
		#print STDERR $_;
		my $line = $_;
		
		
		if (/^\#/) { 
			## retain header for Gustavo's importer.
			print;
			next; 
		}
		unless (/\w/) { next; }
		chomp;
		
                
		my @x = split(/\t/);
		
		
		my $contig = $x[0];
		my $snp_call = $x[4];
		my $snp_quality_tag = $x[6];
		
		#if ($snp_call eq ".") { next; }
		#if ($snp_quality_tag =~ /GATKStandard|LowQual/) { next; } # low quality, ignore.
		
		if ($curr_contig ne $contig) {
			print STDERR "\n-analyzing snps on contig: $curr_contig\n" if $curr_contig;
			&analyze_snps($curr_contig, $snp_text, $gene_obj_indexer, $contig_to_gene_list_href, \%genome, $cod_mut_exp_fh ) if $curr_contig;
			$curr_contig = $contig;
			$snp_text = "";
			$counter = 0;
			print STDERR "\n----\nReading VCF for contig: $contig\n";
		}
		$snp_text .= $line;
		
		$counter++;
		print STDERR "\r[$counter vcf lines read for $contig]     " if ($counter % 1000 == 0);
		
		
	}
	close $fh;
	
	&analyze_snps($curr_contig, $snp_text, $gene_obj_indexer, $contig_to_gene_list_href, \%genome, $cod_mut_exp_fh) if $curr_contig;
	
	print STDERR "\n\nDone.\n\n";
	
	exit(0);
}


####
sub analyze_snps {
	my ($curr_contig, $snp_text, $gene_obj_indexer, $contig_to_gene_list_href, $genome_href, $cod_mut_exp_fh) = @_;
	
	my $contig_seq = $genome_href->{$curr_contig} or die "Error, no sequence for contig: $curr_contig";
	my @contig_chars = split(//, $contig_seq);
	
	## annotate every base.
	my @base_annotations; # populated as genes are encountered below. Maintained as a FIFO
	
	## decorate with annotations.
	
	
	
	my @gene_structs;
	if (exists $contig_to_gene_list_href->{$curr_contig}) {
		
		my @gene_ids = @{$contig_to_gene_list_href->{$curr_contig}};		
		
		
		foreach my $gene_id (@gene_ids) {
			my $gene_obj = $gene_obj_indexer->{$gene_id};
			my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_gene_span();
			push (@gene_structs, { lend => $lend,
					rend => $rend,
					gene => $gene_obj,
			} );
		}
		
		@gene_structs = sort {$a->{lend}<=>$b->{lend}} @gene_structs;
	}
	
	
	my $prev_lend_gene_struct;  # track for intergenics
	
	####
	# annotate the vcf file
	my @snp_lines = split(/\n/, $snp_text);
	
	my $counter = 0;
	print STDERR "-annotating SNPs on contig: $curr_contig\n";
	
	for( my $indSNP = 0; $indSNP < scalar(@snp_lines); $indSNP++ ){  
		my $snp_line = $snp_lines[ $indSNP ];

		
		## define annotations
		my $annotation = "";
		
		# print STDERR "SNPLINE: $snp_line\n";
		$counter++;
		$SNP_LINE_COUNTER++;
		if ($counter % 1000 == 0) { 
			print STDERR "\r[contig_snp_count: $counter, total_snp_count: $SNP_LINE_COUNTER]   ";
		}
		
		my @x = split(/\t/, $snp_line);
		my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $snp_info, @rest) = @x;
		if( $alt =~ /<[\w\W]+>/ ){
			$annotation = '';
			goto SKIP_ANNOTATION;
		}
		
		if (length($ref) == 1 && $contig_chars[$pos-1] ne $ref) {
			print STDERR "\n!! ERROR !! $snp_line\nError, $chrom $pos $ref  != actual base: " . $contig_chars[$pos-1] . "\n";
			die;
		}
		
		## encode base annotations according to gene structure
		while (@gene_structs) {
			if ($gene_structs[0]->{rend} < $pos) {
				# we're past the point of this feature in the contig.
				$prev_lend_gene_struct = shift @gene_structs;
				
			}
			elsif ($gene_structs[0]->{lend} <= $pos && $gene_structs[0]->{rend} >= $pos) {
				# gene overlaps SNP
				$prev_lend_gene_struct = shift @gene_structs;
				&decorate_with_gene_annotations($prev_lend_gene_struct->{gene}, \$contig_seq, \@base_annotations);
			}
			elsif ($gene_structs[0]->{lend} > $pos) {
				# not there yet
				last;
			}
			#print STDERR "Got gene: " . Dumper($gene_structs[0]) . " and position: $pos\n";
		}
		
		
		if (@base_annotations) {
			
			## adjust cursor into base annotations
			
			my $left_index = $base_annotations[0]->{pos};
			unless (defined $left_index) {
				die "Error, no index defined: " . Dumper(\@base_annotations);
			}
			
			while (@base_annotations && ($base_annotations[0]->{pos} < $pos)) {
				my $base_annot = shift @base_annotations;
				unless (defined $base_annot->{pos}) {
					die "Error, no index defined: " . Dumper(\@base_annotations);
				}
			}
			
			if (@base_annotations && ($base_annotations[0]->{pos} < $pos) ) {
				die "Wassup!";
			}
		}
		
		if (@base_annotations) {
			
			my $left_index = $base_annotations[0]->{pos};
			#print Dumper($base_annotations[0]);
			
			if ($left_index == $pos) {
				my $annot_string = $base_annotations[0]->{annot};
				
				my @annotations = split(/\t/, $annot_string);
				foreach my $indiv_annot (@annotations) {
					
					## if a coding mutation, examine the impact on the codon
					if ($indiv_annot =~ /^CDS/) {
						
						## See if two alternative genotypes (REF field containing two values
						## separated by comma, then issue the function with each alternative mutation
						if ($alt =~ /,/ ) {
							my $concat_alt_gen_annot = '';
							foreach my $alt_genotype ( split /,/, $alt ){		
								my $alt_gen_annot;
								#print STDERR "\nALT:$alt_genotype\n";
								#print STDERR $indiv_annot , "\n";
								$alt_gen_annot = &examine_coding_mutation($indiv_annot, $alt_genotype, $cod_mut_exp_fh, \@snp_lines, $indSNP );
								$alt_gen_annot =~ s/[\w\W]+?codon:/codon:/;
								#print STDERR $alt_gen_annot . "\n";
								
								$concat_alt_gen_annot .= $alt_gen_annot . ',';
								#getc();
							}
							$concat_alt_gen_annot =~ s/,$//;
							$indiv_annot =~ s/codon:[\w\W]+$//;
							$indiv_annot .= $concat_alt_gen_annot;
							#print STDERR  $indiv_annot . "\n";
							#getc();
						## see if indel
						}elsif ($snp_info && ((length($alt) > 1 or length($ref) > 1) || $snp_info =~ /INDEL/)) {
							my $indel_len = length($alt) - length($ref);
							my $indel_type = ($indel_len > 0) ? "INSERTION" : "DELETION";
							$indiv_annot .= ",$indel_type\[$indel_len]";
						}
						elsif ($alt ne '.') {
							$indiv_annot = &examine_coding_mutation($indiv_annot, $alt, $cod_mut_exp_fh, \@snp_lines, $indSNP );
						}
					}
				}
				
				$annotation = join("\t", @annotations); ## build back into the single string.
				
				
			}
			else {
				die "Error, should have a base annotation";
			}
		}
		unless ($annotation) {
			# must be intergenic
			$annotation = "intergenic";
			if ($prev_lend_gene_struct) {
				my $prev_gene_obj = $prev_lend_gene_struct->{gene};
				my $prev_gene_id = $prev_gene_obj->{TU_feat_name};
				my $prev_gene_lend = $prev_lend_gene_struct->{lend};
				my $prev_gene_rend = $prev_lend_gene_struct->{rend};
				my $prev_gene_orient = $prev_gene_obj->get_orientation();
				my $delta_pos = $pos - $prev_gene_rend;
				my $prev_gene_name = $prev_gene_obj->{com_name};
				$annotation .= " -- prev_gene($prev_gene_lend-$prev_gene_rend)\[<-($delta_pos)\]: $prev_gene_id $prev_gene_name ($prev_gene_orient) ";
			}
			if (my $next_gene_struct = $gene_structs[0]) {
				my $next_gene_obj = $next_gene_struct->{gene};
				my $next_gene_id = $next_gene_obj->{TU_feat_name};
				my $next_gene_lend = $next_gene_struct->{lend};
				my $next_gene_rend = $next_gene_struct->{rend};
				my $next_gene_orient = $next_gene_obj->get_orientation();
				my $delta_pos = $next_gene_lend - $pos;
				my $next_gene_name = $next_gene_obj->{com_name};
				$annotation .= " -- next_gene($next_gene_lend-$next_gene_rend)\[($delta_pos)-->]: $next_gene_id $next_gene_name ($next_gene_orient) ";
			}
		}

		SKIP_ANNOTATION:
		print join("\t", @x, $annotation) . "\n";
	}
	
	print STDERR "\n"; 
}


####
sub decorate_with_gene_annotations {
	my ($gene_obj, $contig_seq_sref, $base_annotations_aref) = @_;
	
	my ($lend, $rend);
	{
		my @coords;
		foreach my $isoform ($gene_obj, $gene_obj->get_additional_isoforms()) {
			my ($iso_lend, $iso_rend) = $isoform->get_transcript_span();
			push (@coords, $iso_lend, $iso_rend);
		}
		
		@coords = sort {$a<=>$b} @coords;
		$lend = shift @coords;
		$rend = pop @coords;
	}
    	
	
	#@$base_annotations_aref = (); # clear it out   ## NO!! bad idea...  retain it.
	
	my $base_annots_start_add_pos = $lend;
	if (@$base_annotations_aref) {
		$base_annots_start_add_pos = $base_annotations_aref->[$#$base_annotations_aref]->{pos} + 1;
	}
        
	for (my $i = $base_annots_start_add_pos; $i <= $rend; $i++) {
		push (@$base_annotations_aref, { pos => $i,
				annot => "",
		} );
	}
	
	
	my $left_index = $base_annotations_aref->[0]->{pos};
	
	$gene_obj->create_all_sequence_types($contig_seq_sref);
	
	my $gene_id = $gene_obj->{TU_feat_name};
	my $orient = $gene_obj->get_orientation();
	
	foreach my $isoform ($gene_obj, $gene_obj->get_additional_isoforms()) {
		
		my $model_id = $isoform->{Model_feat_name};
		my $com_name = $isoform->{com_name};
		
		my $protein = $isoform->get_protein_sequence();
		my $cds_seq = $isoform->get_CDS_sequence();
		$CDS_SEQUENCES{$model_id} = $cds_seq; # for studying coding mutations later, if they exist
		
		my $complete_protein = ($protein && $protein =~ /^M/ && $protein =~ /\*$/) ? 1 : 0;
		
		## process UTR information.
		if ($isoform->has_UTRs()) {
			
			my @coordsets5p = $isoform->get_5prime_UTR_coords();
			my @coordsets3p = $isoform->get_3prime_UTR_coords();
			
			#print "Got UTRs:" . Dumper(\@coordsets5p) . Dumper(\@coordsets3p) . $isoform->toString();
			
			
			
			
			foreach my $coordset ($isoform->get_5prime_UTR_coords()) {
				my ($lend, $rend) = sort {$a<=>$b} @$coordset;
				#print STDERR "\tprocessing 5pUTR $lend-$rend.\n";
				
				for (my $i = $lend; $i <= $rend; $i++) {
					my $annot_text = join(",", "p5UTR", $gene_id, $model_id, $com_name, $orient);
					&add_annotation($base_annotations_aref, $i, $annot_text);
				}
			}
			
			foreach my $coordset ($isoform->get_3prime_UTR_coords()) {
				my ($lend, $rend) = sort {$a<=>$b} @$coordset;
				#print STDERR "\tprocessing 3pUTR $lend-$rend.\n";
				
				for (my $i = $lend; $i <= $rend; $i++) {
					my $annot_text = join(",", "p3UTR", $gene_id, $model_id, $com_name, $orient);
					&add_annotation($base_annotations_aref, $i, $annot_text);
				}
			}
		}
		
		## annotate introns:
		if (my @introns = $isoform->get_intron_coordinates()) {
			
			foreach my $intron (@introns) {
				
				my ($lend, $rend) = sort {$a<=>$b} @$intron;

				# Dinucleotide signal on the left
				my $left_signal_start = $lend;
				my $left_signal_end   = ($rend == $lend)?$rend:$lend+1;

				# Dinucleotide signal on the right
				my $right_signal_start   = ($rend == $lend)?$lend:$rend-1;
				my $right_signal_end     = $rend;
				

				# print STDERR "\tprocessing intron $lend-$rend\n";
				for (my $i = $lend; $i <= $rend; $i++) {
					
					my $qualifier = 'intron';
					
					if( $i >= $left_signal_start && $i <= $left_signal_end ){
						my $signal_type = ( $orient eq '+' )?'donor':'acceptor';
						$qualifier .= "/splice_$signal_type";
					}
					
					if( $i >= $right_signal_start && $i <= $right_signal_end ){
						my $signal_type = ( $orient eq '+' )?'acceptor':'donor';
						$qualifier .= "/splice_$signal_type";
					}
					
					
					my $annot_text = join(",", $qualifier, $gene_id, $model_id, $com_name, $orient);
					&add_annotation($base_annotations_aref, $i, $annot_text);
					
					
				}
				
				
			}
			
		}
		
		## Process coding regions of exons.
		
		
		my @exons = sort {$a->{end5}<=>$b->{end5}} $isoform->get_exons();
		
		if ($orient eq '-') {
			@exons = reverse @exons;
		}
		
		
		
		my $sum_cds_length = 0;
		foreach my $exon (@exons) {
			
			if (my $cds_obj = $exon->get_CDS_obj()) {
				
				my $end5 = $cds_obj->{end5};
				my $end3 = $cds_obj->{end3};
				# print STDERR "\tprocessing exon $end5-$end3\n";
				
				if ($orient eq '+') {
					for (my $i = $end5; $i <= $end3; $i++) {
						
						$sum_cds_length++;
						my $frame = (($sum_cds_length -1) % 3 ) + 1;
						my $codon = &get_codon(\$cds_seq, $sum_cds_length);
						my $annot = join(",", "CDS", $gene_id, $model_id, $com_name, "trans_orient:$orient", "loc_in_cds:$sum_cds_length", "codon_pos:$frame", "codon:$codon");
						# print "$annot\n";
						&add_annotation($base_annotations_aref, $i, $annot);
					}
				}
				else {
					# minus orientation
					for (my $i = $end5; $i >= $end3; $i--) {
						
						$sum_cds_length++;
						
						my $frame = (($sum_cds_length -1) % 3 ) + 1;
						my $codon = &get_codon(\$cds_seq, $sum_cds_length);
						my $annot = join(",", "CDS", $gene_id, $model_id, $com_name, "trans_orient:$orient", "loc_in_cds:$sum_cds_length", "codon_pos:$frame", "codon:$codon");
						&add_annotation($base_annotations_aref, $i, $annot);
						
					}
				}
				
			}
		}
		
		
	}
	
	
	
	return;
}


####
sub add_annotation {
	my ($base_annotations_aref, $position, $annotation) = @_;
	
	my $left_pos = $base_annotations_aref->[0]->{pos};
	my $right_pos = $base_annotations_aref->[$#$base_annotations_aref]->{pos};
	
	if (! defined ($left_pos) || ! defined ($right_pos)) {
		croak "Error, existing base_annotations_aref lacks boundary positions defined: " . Dumper(\$base_annotations_aref);
	}
	
	if ($position > $right_pos) {
		croak "Error, position: $position is out of range: $left_pos-$right_pos" . Dumper($base_annotations_aref);
	}
	
	if ($position < $left_pos) {
		croak "Error, position: $position is out of range: $left_pos-$right_pos" . Dumper($base_annotations_aref);
	}
	
	my $delta = $position - $left_pos;
	if ($base_annotations_aref->[$delta]->{pos} != $position) {
		croak "Error, position: $position is not at $delta in the existing list: " . Dumper($base_annotations_aref->[$delta]);
		
	}
	
	## if got this far, all should be OK AFAICT
	
	if ($base_annotations_aref->[$delta]->{annot}) {
		$base_annotations_aref->[$delta]->{annot} .= "\t";
	}
	
	$base_annotations_aref->[$delta]->{annot} .= $annotation;
	
	return;
}

####
sub get_codon {
	my ($cds_seq_sref, $pos) = @_;
	
	
	my $last_pos = length($$cds_seq_sref) - 1;
	
	my $codon_pos = int(($pos-1)/3) * 3;
	
	my $codon_string = "";
	for (my $i = $codon_pos; $i <= $codon_pos + 2; $i++) {
		if ($i > $last_pos) { 
			last;
		}
		$codon_string .= substr($$cds_seq_sref, $i, 1);
	}
	
	return($codon_string);
}

####
sub examine_coding_mutation {
	my ($annotation, $alt_base, $cod_mut_exp_fh, $refArrSNPs, $indSNP ) = @_;
	
	my @x = split(/,/, $annotation);
	
	my $model_id = $x[2];
	my $cds_seq = $CDS_SEQUENCES{$model_id} or die "Error, no CDS sequence for $model_id, [annot: $annotation]";
	
	
	my $codon = pop @x;
	$codon =~ s/codon://;
	
	
	
	if (length($codon) != 3) {
		return($annotation); # unchanged
	}
	
	$annotation = join(",", @x);  #stripped of codon. We'll add it back shortly, and decorate it.
	
	
	my $codon_pos = pop @x;
	$codon_pos =~ s/^\S+://;
	
	my $cds_base_pos = pop @x;
	$cds_base_pos =~ s/^\S+://;
	
	my $aa_pos = int(($cds_base_pos-1)/3)+1;
	
	my $transcript_orientation = pop @x;
	$transcript_orientation =~ s/^\S+://;
	
	if ($transcript_orientation eq '-') {
		$alt_base = &reverse_complement($alt_base);
	}
	
	
	my $original_aa = &translate_sequence($codon, 1);
	
	## mutate the codon	
	my $mutated_codon;
	
	if( $codon_based ){
		$mutated_codon = mutate_codon( $alt_base, $transcript_orientation, $codon_pos, $codon, $refArrSNPs, $indSNP );	
	}else{
		my @codon_chars = split(//, $codon);
		$codon_chars[$codon_pos-1] = $alt_base;
		$mutated_codon = join("", @codon_chars);
	}
	
	
	my $mutated_aa = &translate_sequence($mutated_codon, 1);
	
	{
		## make codons pretty and case-informative
		my @codon_chars = split(//, lc $codon);
		$codon_chars[$codon_pos-1] = uc $codon_chars[$codon_pos-1];
		$codon = join("", @codon_chars);
		
		my @mutated_codon_chars = split(//, lc $mutated_codon);
		$mutated_codon_chars[$codon_pos-1] = uc $mutated_codon_chars[$codon_pos-1];
		$mutated_codon = join("", @mutated_codon_chars);
	}
	
	
	my $mutation_type = "";
	if ($mutated_aa eq $original_aa) {
		$mutation_type = "SYN";
	}
	elsif ($mutated_aa eq '*' && $original_aa ne '*') {
		$mutation_type = "NON";
	}
	elsif ($original_aa eq '*' && $mutated_aa ne '*') {
		$mutation_type = "RTH";
	}
	else {
		$mutation_type = "NSY";
	}
	
	$annotation .= ",codon:$codon-$mutated_codon,pep:$original_aa\->$mutated_aa," 
        . &get_amino_acid_triplet_label($original_aa) . "-" . $aa_pos . "-" 
        . &get_amino_acid_triplet_label($mutated_aa) . ",($mutation_type)";
        
        
        
        if ($cod_mut_exp_fh) {
        	
        	my $coding_mutation_explained_text = 
        	&explain_coding_mutation_pedantically($cds_seq, $cds_base_pos, $codon, $mutated_codon, $annotation);
        	
        	print $cod_mut_exp_fh ">$model_id $annotation\n$coding_mutation_explained_text\n";
        }
        
        
        return($annotation);
}


####
sub explain_coding_mutation_pedantically {
	my ($cds_seq, $cds_base_pos, $codon, $mutated_codon, $annotation) = @_;
	
	my $ret_text = "";
	
	my @codons;
	while ($cds_seq =~ /(\S{3})/g) {
		push (@codons, $1);
	}
	
	my $mutated_codon_index = int(($cds_base_pos-1)/3);
	
	my $codon_counter = 0;
	
	my $number_line = "";
	my $codon_line = "";
	my $translation_line = "";
	my $mutation_codon_line = "";
	my $mutation_translation_line = "";
	
	my $codons_per_line = 15;
	
	for (my $i = 0; $i < $#codons; $i++) {
		if ($i > 0 && $i % $codons_per_line == 0) {
			foreach my $line ($number_line, $codon_line, $translation_line, $mutation_codon_line, $mutation_translation_line) {
				$ret_text .= "$line\n" if ($line =~ /\w/);
			}
			$ret_text .= "\n";
			
			## reinit
			$number_line = "";
			$codon_line = "";
			$translation_line = "";
			$mutation_codon_line = "";
			$mutation_translation_line = "";
			
			
		}
		
		$number_line .= sprintf("%5d", $i+1); # use codon numbering starting from 1 rather than zero
		$codon_line .= sprintf("%5s", $codons[$i]);
		$translation_line .= sprintf("%4s ", &translate_sequence($codons[$i], 1));
		
		if ($i == $mutated_codon_index) {
			$mutation_codon_line .= sprintf("%5s", $mutated_codon);
			$mutation_translation_line .= sprintf("%4s ", &translate_sequence($mutated_codon, 1)) . "\t$annotation";
		}
		else {
			$mutation_codon_line .= " " x 5;
			$mutation_translation_line .= " " x 5;
		}
	}
	
	if ($number_line) {
		foreach my $line ($number_line, $codon_line, $translation_line, $mutation_codon_line, $mutation_translation_line) {
			$ret_text .= "$line\n" if ($line =~ /\w/);
		}
		$ret_text .= "\n";
	}
	
	return($ret_text);
}

####
sub get_amino_acid_triplet_label {
	my ($aa_letter) = @_;
	
	my %triplets = ( 'A' => "Ala",
		'C' => "Cys",
		'D' => 'Asp',
		'E' => 'Glu',
		'F' => 'Phe',
		'G' => 'Gly',
		'H' => 'His',
		'I' => 'Ile',
		'K' => 'Lys',
		'L' => 'Leu',
		'M' => 'Met',
		'N' => 'Asn',
		'P' => 'Pro',
		'Q' => 'Gln',
		'R' => 'Arg',
		'S' => 'Ser',
		'T' => 'Thr',
		'V' => 'Val',
		'W' => 'Trp',
		'Y' => 'Tyr',
		'X' => 'XXX',
		'*' => 'STP',
		);
	
	my $triplet = $triplets{ uc $aa_letter };
	
	unless ($triplet) {
		die "Error, residue($aa_letter) is unknown";
	}
	
	return($triplet);
}

sub mutate_codon{
	my ( $alt_base, $transcript_orientation, $codon_pos, $codon, $refArrSNPs, $indSNP ) = @_;
	
	# Set initially the mutate codon as the original codon
	my @alt_codon = split(//, $codon);
	
	$alt_codon[  $codon_pos - 1 ] =  $alt_base;
	
	# If alternate base is the first position in the codon
	my $alt = '';
	if( $codon_pos == 1 ){
		# get substitution in the next position
		print STDERR "get substitution in the next position\n" if $debug;;
		$alt = find_next_substitution( 1, $transcript_orientation, $refArrSNPs, $indSNP );
		$alt_codon[ 1 ] = $alt if $alt ne '';
		
		# get substitution 2 positions apart from it 
		$alt = find_next_substitution( 2, $transcript_orientation, $refArrSNPs, $indSNP );
		$alt_codon[ 2 ] = $alt if $alt ne '';
	}
	
	if( $codon_pos == 2 ){
		# get substitution in the previous position
		print STDERR "get substitution in the previous position\n" if $debug;;
		$alt = find_next_substitution( - 1, $transcript_orientation, $refArrSNPs, $indSNP );
		$alt_codon[ 0 ] = $alt if $alt ne '';
		
		# get substitution in the next position
		$alt = find_next_substitution( 1  , $transcript_orientation, $refArrSNPs, $indSNP );
		$alt_codon[ 2 ] = $alt if $alt ne '';
	}
	
	if( $codon_pos == 3 ){
		# get substitution two positions apart (upstream)
		print STDERR "get substitution two positions apart (upstream)\n" if $debug;;
		$alt = find_next_substitution( - 2, $transcript_orientation, $refArrSNPs, $indSNP );
		$alt_codon[ 0 ] = $alt if $alt ne '';
		
		# get substitution in the previous position
		$alt = find_next_substitution( - 1, $transcript_orientation, $refArrSNPs, $indSNP );
		$alt_codon[ 1 ] = $alt if $alt ne '';
	}
	
	
	my $mutated_codon = join("", @alt_codon);
	return $mutated_codon;
}


sub find_next_substitution{
	my ( $genomic_dist, $transcript_orientation, $refArrSNPs, $indSNP ) = @_;
	my $snp_line = $refArrSNPs->[ $indSNP ];
	
	my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $snp_info, @rest) = split(/\t/, $snp_line);;
	
	$genomic_dist = $genomic_dist * -1 if( $transcript_orientation eq '-' );
	
	my $sought_pos = $pos + $genomic_dist;

	print STDERR "Current position: $pos  sought position: $sought_pos  Genomic dist: $genomic_dist \n" if $debug;

	my $seek_ind = $indSNP;
	
	my $next_subst = '';
	
	# Seeking forward
	if( $genomic_dist > 0 ){
		
		# Set position to minus so it can always pass the first
		# iteration of the while
		$pos = -1;
		$seek_ind++;
		while( $pos != $sought_pos && $seek_ind < scalar(@{$refArrSNPs} ) ){
			my $line = $refArrSNPs->[ $seek_ind ];	
			($chrom, $pos, $id, $ref, $alt, $qual, $filter, $snp_info, @rest) = split(/\t/, $line);
	
			print STDERR "Seeking forward, current position: $pos\n" if $debug;;			
			
			# Return nothing if pos passed the position that we are looking for
			# but no substitution was found
			if( $pos > $sought_pos ){ 
				print STDERR "Sought position, NOT FOUND!!!\n" if $debug;;				
				return '';
			}
			$seek_ind++;
		}
		
		# Seeking backward 
	}else{
		
		# Set position to minus so it can always pass the first
		# iteration of the while
		$pos = -1;
		$seek_ind--;
		while( $pos != $sought_pos && $seek_ind >= 0 ){
			my $line = $refArrSNPs->[ $seek_ind ];	
			($chrom, $pos, $id, $ref, $alt, $qual, $filter, $snp_info, @rest) = split(/\t/, $line);

			print STDERR "Seeking backwards, current position: $pos\n" if $debug;;			
			
			# Return nothing if pos passed the position that we are looking for
			# but no substitution was found
			if( $pos < $sought_pos ){
				print STDERR "Sought position, NOT FOUND!!!\n" if $debug;				
				return '';
			}
			$seek_ind--;
		};		
		
	}
	 
	
	print STDERR "Found sought position. Curr pos: $pos, Sought position: $sought_pos alt: $alt!!!\n" if( $pos == $sought_pos && $debug );
	
	# If pos found and if it is a substitution
	if( $pos == $sought_pos && $alt ne '.' && length($alt) == 1 ){
		
		if ($transcript_orientation eq '-') {
			$alt = &reverse_complement($alt);
		}
		
		print STDERR "Returning alt: $alt\n" if $debug;
		return $alt;
	}
	
	return '';
	
}


