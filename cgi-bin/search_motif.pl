#!/usr/bin/perl -w
=pod

=head1 NAME

RTScan.pl 

=head1 SYNOPSIS

perl RTScan.pl --help

=head1 REQUIRES

Perl5
Bioperl
NCBI-Blast
EMBOSS

=head1 DESCRIPTION

=cut

use lib '/home/gdroc/lib';
use strict;
use FindBin;
use File::Temp qw/ tempfile tempdir /;
use CGI;
use JSON; 
use lib "$FindBin::Bin/lib";
use warnings; 
use Bio::SeqIO;
use File::Basename;
use Bio::Index::Fasta;
use Bio::SearchIO;
use FileUpload;
use Spreadsheet::WriteExcel;
use Data::Dumper;
our $WIDTH      = 50;
our $BLASTALL   = "/home/gdroc/bin/blastall";
our $FORMATDB   = "/home/gdroc/bin/formatdb";
our $WORDMATCH  = "/home/gdroc/bin/bin/wordmatch";
our $WORDFINDER = "/home/gdroc/bin/bin/wordfinder";
our $CDHIT      = "/home/gdroc/bin/cdhit";
our $DEBUG      = 0;

my $param = new CGI;
my $program             = "blastn"; 
my $format              = 0;
my $thread              = 8;
#my $extend_hit_region   = 150;
my $interval            = 3;
my $threshold_hit_polyA = 20;
my $threshold_tsd       = 10;

my $global 		      = $param->param("global");
my $query_name 		  = $param->param("query_name");
my $hit_name   		  = $param->param("hit_name");
my $abbrev     		  = $param->param("abbrev");
my $hsp_start_hit  	  = $param->param("hsp_start_hit");
my $species	     	  = $param->param("species");
my $hsp_end_hit    	  = $param->param("hsp_end_hit");
my $hsp_start_query   = $param->param("hsp_start_query");
my $percent_identity  = $param->param("percent_identity");
my $hsp_end_query     = $param->param("hsp_end_query");
my $hsp_length   	  = $param->param("hsp_length");
my $gaps              = $param->param("gaps");
my $strand     		  = $param->param("strand");
my $cpt_hsp    		  = $param->param("id");
my $width      		  = $param->param("width");
my $database   		  = $param->param("database"); 
my $output_dir 		  = $param->param("output_dir");
my $result_dir 		  = $param->param("result_dir");
my $query_length 	  = $param->param("query_length");
my $tsd_length 		  = $param->param("tsd_length");
my $repeat_max_length = $param->param("repeat_max_length");

my $file_txt		  = $param->param("file_txt");
my $file_tsd_txt  	  = $param->param("file_tsd_txt");
my $file_fasta		  = $param->param("file_fasta");
my $file_xls		  = $param->param("file_xls");
my $file_url_txt	  = $param->param("file_url_txt");
my $file_url_fasta	  = $param->param("file_url_fasta");
my $file_url_xls	  = $param->param("file_url_xls");
my $file_url_tsd_txt  = $param->param("file_url_tsd_txt"); 


	


our $repeat = "/var/www/bank/repbase/repbase.fna";
 
open( TXT, ">>$file_txt" ); 
print $param->header('application/json');
my $cpt_tsd           = 0;
my $overlap           = 8; 
my $extend_sequence   = $width * 2; 
my $extend_hit_region = 150;
my $seqobj = &get_index( $hit_name, $database );
my $hit_length = $seqobj->length;

my @sequence;
$strand = $strand == 1 ? "+" : "-";
my $id = join( "_", $query_name, $hit_name, $hsp_start_hit, $hsp_end_hit );
my $query_id = $species eq "" ? join( "_", $abbrev, sprintf( "%03d",$cpt_hsp)) : join( "_", $abbrev, sprintf( "%03d",$cpt_hsp));
my $id_dir = $output_dir . "/" . $query_id;
system("mkdir $id_dir"); 

my $start_extend5 = $strand eq "+" ? $hsp_start_hit - $extend_hit_region : $hsp_end_hit + 1 - $overlap;
my $end_extend5 = $strand eq "+" ? $hsp_start_hit - 1 + $overlap : $hsp_end_hit + $extend_hit_region;
my $start_extend3 = $strand eq "+" ? $hsp_end_hit + 1 - $overlap : $hsp_start_hit - $extend_hit_region;
my $end_extend3 = $strand eq "+"  ? $hsp_end_hit + $extend_hit_region : $hsp_start_hit - 1 + $overlap;
 	
$start_extend5 = 1           if $start_extend5 < 1;
$start_extend3 = 1           if $start_extend3 < 1;
$end_extend5   = $hit_length if $end_extend5 > $hit_length;
$end_extend3   = $hit_length if $end_extend3 > $hit_length;
my ( $file_five, $sequence_five_prime )   = &get_sequence( $id_dir, $query_id . "_5", $start_extend5, $end_extend5, $seqobj, $strand, $width );
my ( $file_three, $sequence_three_prime ) = &get_sequence( $id_dir, $query_id . "_3", $start_extend3, $end_extend3, $seqobj, $strand, $width );
my ( $file_hit, $sequence_hit )           = &get_sequence( $id_dir, $query_id . "_hit",$hsp_start_hit, $hsp_end_hit, $seqobj, $strand, $width );
my $hit_seq = substr( $sequence_hit, 0, 5 ) . "..". substr( $sequence_hit, -5 );
system("rm $file_hit") if -e $file_hit;
my ( $current_start5_tsd, $current_end5_tsd, $current_start3_tsd, $current_end3_tsd, $current_length ) = &wordmatch( $file_five, $file_three, $query_id ); 
my $start5_tsd = 0;
my $end5_tsd   = 0;
my $start3_tsd = 0;
my $end3_tsd   = 0;  
my $mismatch   = 0; 
my $start_cleavage_before;
my $end_cleavage_before;
my $start_cleavage_after;
my $end_cleavage_after;
my $repeat_size = "";
my $tsd_overlap5 = 0;
my $tsd_overlap3 = 0;
my $found_tsd    = 0;
my $found_polyA  = 0;	
my $found_repeat = 0;
my @repeat;
my $cpt_repeat = 0;
my %polyA; 
my $start_polyA  = 0;
my $end_polyA    = 0;	
my $subseq_polyA = "";
my $length_polyA = 0;  
my $seq_polyA = "";
my $algorithm = "none";
$algorithm    = "wordmatch" if $current_start5_tsd > 0;
my ( $file_cleavage_before_reverse, $sequence_cleavage_before_reverse, $reverse_strand, $file_cleavage_after_reverse, $sequence_cleavage_after_reverse);
my $polyA_ref = "";
my $file_five_tsd = "";
my $sequence_five = "";
my $cleavage = "";
my $file_three_tsd = ""; 
my $sequence_three = "";
unless ( $current_start5_tsd > 0 ) {
	$start_extend3 = $strand eq "+" ? $hsp_end_hit + 1 - $overlap : $hsp_start_hit - $extend_hit_region - $repeat_max_length;
	$end_extend3   = $strand eq "+" ? $hsp_end_hit + $extend_hit_region + $repeat_max_length : $hsp_start_hit - 1 + $overlap;
	$start_extend5 = 1           if $start_extend5 < 1;
	$start_extend3 = 1           if $start_extend3 < 1;
	$end_extend5   = $hit_length if $end_extend5 > $hit_length;
	$end_extend3   = $hit_length if $end_extend3 > $hit_length;
	( $file_three, $sequence_three_prime ) = &get_sequence( $id_dir, $query_id . "_3_1", $start_extend3, $end_extend3, $seqobj, $strand, $width );
	( $current_start5_tsd, $current_end5_tsd, $current_start3_tsd, $current_end3_tsd, $current_length ) = &wordmatch( $file_five, $file_three, $query_id );
	$algorithm     = "wordmatch" if $current_start5_tsd > 0;
	unless ( $current_start5_tsd > 0 ) {
		$start_extend5 = $strand eq "+" ? $hsp_start_hit - 30 : $hsp_end_hit + 1;
		$end_extend5   = $strand eq "+" ? $hsp_start_hit - 1 : $hsp_end_hit + 30;
		$start_extend5 = 1 if $start_extend5 < 1;
		$end_extend5 = $hit_length if $end_extend5 > $hit_length;
		my ( $file_five, $sequence_five_prime ) = &get_sequence( $id_dir, $query_id . "_5", $start_extend5, $end_extend5, $seqobj, $strand, $width );
		( $current_start5_tsd, $current_end5_tsd, $current_start3_tsd, $current_end3_tsd, $current_length, $mismatch ) = &wordfinder( $file_five, $file_three, $query_id );
		if ( $current_start5_tsd > 0 ) {
			$algorithm = "wordfinder";
		}
	}
}
 
if ( $current_start5_tsd > 0 ) { 
	$found_tsd  = 1;
	$start5_tsd = $strand eq "+" ? ( $start_extend5 + $current_start5_tsd - 1 ) : ( $end_extend5 - $current_end5_tsd + 1 );
	$end5_tsd   = $strand eq "+" ? ( $start_extend5 + $current_end5_tsd - 1 ) : ( $end_extend5 - $current_start5_tsd + 1 );
	$start3_tsd = $strand eq "+" ? ( $start_extend3 + $current_start3_tsd - 1 ) : ( $end_extend3 - $current_end3_tsd + 1 );
	$end3_tsd = $strand eq "+" ? ( $start_extend3 + $current_end3_tsd - 1 ) : ( $end_extend3 - $current_start3_tsd + 1 ); 
	$start_cleavage_before = $strand eq "+" ? $start5_tsd : $end5_tsd - 4;
	$end_cleavage_before =  $strand eq "+" ? $start5_tsd + 4 : $end5_tsd;
	$start_cleavage_after = $strand eq "+" ? $start5_tsd - 3 : $end5_tsd + 1;
	$end_cleavage_after = $strand eq "+" ? $start5_tsd - 1 : $end5_tsd + 3;
	$reverse_strand = $strand eq "+" ? "-" : "+";
	( $file_cleavage_before_reverse, $sequence_cleavage_before_reverse ) = &get_sequence( $id_dir, $query_id . "_cleavage_before", $start_cleavage_before, $end_cleavage_before, $seqobj, $reverse_strand, $width );
	( $file_cleavage_after_reverse, $sequence_cleavage_after_reverse )  = &get_sequence( $id_dir, $query_id . "_cleavage_after", $start_cleavage_after, $end_cleavage_after, $seqobj, $reverse_strand, $width );
	$cleavage = join( "-", $sequence_cleavage_before_reverse, $sequence_cleavage_after_reverse );
 
	( $file_five_tsd, $sequence_five ) =  &get_sequence( $id_dir, $query_id . "_5tsd", $start5_tsd, $end5_tsd, $seqobj, $strand, $width );
	( $file_three_tsd, $sequence_three ) = &get_sequence( $id_dir, $query_id . "_3tsd", $start3_tsd, $end3_tsd, $seqobj, $strand, $width ); 
	$start_extend5 = $start5_tsd;
	$end_extend5   = $end5_tsd;
	$start_extend3 = $start3_tsd;
	$end_extend3   = $end3_tsd;
	$tsd_overlap5  = $strand eq "+" ? $end5_tsd - $hsp_start_hit + 1 : $hsp_end_hit - $start5_tsd + 1;
	$tsd_overlap3 = $strand eq "+" ? $hsp_end_hit - $start3_tsd + 1 : $end3_tsd - $hsp_start_hit + 1;
	$tsd_overlap5 = 0 if $tsd_overlap5 < 0;
	$tsd_overlap3 = 0 if $tsd_overlap3 < 0;
	$repeat_size  = $strand eq "+" ? $start3_tsd - $hsp_end_hit + 1 : $hsp_start_hit - $end3_tsd + 1; 
}
else {
	$end_extend3 = $strand eq "+" ? $hsp_end_hit + 1 - $overlap : $hsp_start_hit - $extend_hit_region - $repeat_max_length;
	$start_extend3 = $strand eq "+" ? $hsp_end_hit + $extend_hit_region + $repeat_max_length : $hsp_start_hit - 1 + $overlap;
	$end_extend3 = 1 if $end_extend3 < 1;
	$start_extend3 = $hit_length if $start_extend3 > $hit_length;  
}
my $sequence_repeat_polyA = "";
my $start_hit_tsd3 = $strand eq "+" ? $hsp_end_hit + 1 : $end_extend3 + 1;
my $end_hit_tsd3 = $strand eq "+" ? $start_extend3 - 1 : $hsp_start_hit - 1;
my $length           = $end_hit_tsd3 - $start_hit_tsd3;
my $start_hit_repeat = $start_hit_tsd3;
my $end_hit_repeat   = $end_hit_tsd3;
my $repeat_length;
my $file_repeat;
my %repeat_family;
if ( $length > 0 ) {
	my ( $file_repeat, $sequence_repeat ) = &get_sequence( $id_dir, $query_id . "_repeat_blast", $start_hit_tsd3, $end_hit_tsd3, $seqobj, $strand, $width );
	my $file_repeat_out =  $output_dir . "/blast/" . $query_id . "_repeat.out";
	my $max_end_repeat   = 0;
	my $min_start_repeat = 0;
	$file_repeat_out =  &run_blast( $file_repeat, $repeat, $file_repeat_out, "1e-10" , " -v 50 -b 50 -K 50");
	my @hsp = sort { $a->{hsp}->start('query') <=> $b->{hsp}->start('query') } @{ &parse_blast_result($file_repeat_out) };
	system("rm $file_repeat") if -e $file_repeat;
	if ( scalar(@hsp) > 0 ) {
    	my $start_repeat_init = $strand eq "+" ? $start_hit_tsd3 + $hsp[0]->{hsp}->start('query') - 1 : $start_hit_tsd3 + $hsp[-1]->{hsp}->end('query');
		my $size = $strand eq "+" ? $start_repeat_init - $start_hit_tsd3 : $end_hit_tsd3 - $start_repeat_init;
		if ( $size < 30 ) {
			foreach my $repeat_hsp (@hsp) {
				$cpt_repeat++;
				my $start_repeat = $start_hit_tsd3 + $repeat_hsp->{hsp}->start('query') - 1;
				my $end_repeat = $start_hit_tsd3 + $repeat_hsp->{hsp}->end('query');
				my $family =  ( split( /\#/, $repeat_hsp->{hit} ) )[1];
				$repeat_family{$query_id}{$family}++;
				push @repeat,
				{
				    start_repeat => $start_repeat,
				    hit_name     => $repeat_hsp->{hit},
				    end_repeat   => $end_repeat,
				    hsp_length   => $repeat_hsp->{hsp}->hsp_length(),
				    family => $family
				    };
				if ($min_start_repeat) {
				    $max_end_repeat = $end_repeat + 1
					if $end_repeat > $max_end_repeat;
				    $min_start_repeat = $start_repeat - 1
					if $start_repeat < $min_start_repeat;    
				}
				else {
				    $max_end_repeat   = $end_repeat + 1;
				    $min_start_repeat = $start_repeat - 1;
				}
			}
		}
	}
	@repeat = sort { $b->{start_repeat} <=> $a->{start_repeat} } @repeat if $strand eq "-";
	$found_repeat = 1 if @repeat;
	if ( $max_end_repeat > 0 && $strand eq "+" ) {
		$start_hit_repeat = $max_end_repeat;
	}
	if ( $min_start_repeat > 0 && $strand eq "-" ) {
		$end_hit_repeat = $min_start_repeat;
	}
	if ( $end_hit_repeat - $start_hit_repeat > 1 ) { 
	    ( $file_repeat, $sequence_repeat_polyA ) = &get_sequence( $id_dir, $query_id . "_polyA", $start_hit_repeat, $end_hit_repeat, $seqobj, $strand, $width );
		( $polyA_ref, $found_polyA )      = &polyA_search($sequence_repeat_polyA);
		%polyA = %$polyA_ref;  
	} 
	
}
if ( $found_polyA == 1 ) {
	my $cpt_polyA = 0;
	foreach my $index ( sort { $a <=> $b } keys %polyA ) {
		$cpt_polyA++;
		$subseq_polyA = $polyA{$index}{subseq};
		$length_polyA = $polyA{$index}{length};
		$start_polyA  = $strand eq "+" ? $start_hit_repeat + $index : $end_hit_repeat - $index - $polyA{$index}{length} + 1;
		$end_polyA    = $strand eq "+" ? $start_polyA + $polyA{$index}{length} - 1 : $end_hit_repeat - $index;
		last if $found_tsd == 0; 
	}
} 
my $start_extend       = $hsp_start_hit - $extend_sequence;
my $end_extend         = $hsp_end_hit + $extend_sequence;
my $interval_hit_polyA = 0;
		
if ( $found_tsd == 1 ) {
	$start_extend = $start5_tsd - $extend_sequence if $strand eq "+";
	$end_extend = $end5_tsd + $extend_sequence if $strand eq "-";
}
if (@repeat) {
	for ( my $i = 0 ; $i <= $#repeat ; $i++ ) {
		$repeat_length = $repeat[$i]->{hsp_length};
		$end_extend    = $repeat[$i]->{end_repeat} + $extend_sequence if $strand eq "+";
		$start_extend = $repeat[$i]->{start_repeat} - $extend_sequence if $strand eq "-";
	}
}
if ( $start_polyA > 0 ) {
	$end_extend = $end_polyA + $extend_sequence  if $strand eq "+";
	$start_extend = $start_polyA - $extend_sequence if $strand eq "-";
	$interval_hit_polyA = $strand eq "+" ? $start_polyA - $hsp_end_hit : $hsp_start_hit - $end_polyA;
}
if ( $found_tsd == 1 ) {
	$end_extend   = $end3_tsd + $extend_sequence if $strand eq "+";
	$start_extend = $start3_tsd - $extend_sequence if $strand eq "-";
}
$end_extend   = $hit_length if $end_extend > $hit_length;
$start_extend = 1 if $start_extend < 1;

my ( $file_complete, $sequence_complete ) =  &get_sequence( $result_dir, $query_id, $start_extend, $end_extend, $seqobj, $strand, $width );
system("rm -f $file_complete") if -e $file_complete;
 ### Get Sequence
my ($file_sequence_before,$sequence_before);
my ($file_sequence_after,$sequence_after);
my ($file_between_hit_polyA);
my $sequence_between_hit_polyA = "";
my $new_start_hit = $hsp_start_hit;
my $new_end_hit   = $hsp_end_hit;
my $file_between_hit_tsd;
my $sequence_between_hit_tsd = "";
#$sequence_three .= " (TSD3 : ".$start3_tsd."..".$end3_tsd.")"  if $found_tsd==1;
#$hit_seq .= " (Hit : ".$hsp_start_hit."..".$hsp_end_hit.")" ;
my $file_between_polyA_tsd ;
my $sequence_between_polyA_tsd = "";
if ($found_tsd == 1){
	my $start_extend_before = $strand eq "+" ? $start5_tsd - $extend_sequence : $end5_tsd + 1;
	my $end_extend_before   = $strand eq "+" ? $start5_tsd - 1 : $end5_tsd + $extend_sequence;
	($file_sequence_before,$sequence_before) = &get_sequence( $id_dir, $query_id . "_hit",$start_extend_before, $end_extend_before , $seqobj, $strand, $width );
	system("rm -f $file_sequence_before") if -e $file_sequence_before;
	#$sequence_five .= " (Extend5 : ".$start_extend_before."..".$end_extend_before.") + " ;
	my $start_extend_after = $strand eq "+" ? $end3_tsd + 1 : $start3_tsd - $extend_sequence;
	my $end_extend_after   = $strand eq "+" ? $end3_tsd + $extend_sequence : $start3_tsd - 1;
	($file_sequence_after,$sequence_after) = &get_sequence( $id_dir, $query_id . "_hit",$start_extend_after, $end_extend_after , $seqobj, $strand, $width );
	system("rm -f $file_sequence_after") if -e $file_sequence_after;
	#$sequence_three .= " + (Extend3 : ".$start_extend_after."..".$end_extend_after.")" ;
	my $interval_tsd_hit = $strand eq "+" ? $hsp_start_hit - $end5_tsd : $start5_tsd - $hsp_end_hit;
	if ($tsd_overlap5) {
		$new_start_hit =  $strand eq "+" ? $hsp_start_hit + $tsd_overlap5 : $hsp_start_hit;
		$new_end_hit   = $strand eq "+" ? $hsp_end_hit : $hsp_end_hit - $tsd_overlap5;
		( $file_hit, $sequence_hit )           = &get_sequence( $id_dir, $query_id . "_hit",$new_start_hit, $new_end_hit, $seqobj, $strand, $width ) ;
		system("rm -f $file_hit") if -e $file_hit;
		#$hit_seq .= " + (Overlap : ".$new_start_hit."..".$new_end_hit.")" ;
	} 
	if ($tsd_overlap3) {
		$new_start_hit =  $strand eq "+" ? $hsp_start_hit : $hsp_start_hit + $tsd_overlap3;
		$new_end_hit   = $strand eq "+" ? $hsp_end_hit - $tsd_overlap3 : $hsp_end_hit;
		( $file_hit, $sequence_hit )           = &get_sequence( $id_dir, $query_id . "_hit",$new_start_hit, $new_end_hit, $seqobj, $strand, $width ) ; 
		system("rm -f $file_hit") if -e $file_hit;
		#$hit_seq .= " + (Overlap : ".$new_start_hit."..".$new_end_hit.")" ;
	}
	if ($interval_tsd_hit > 1) { 
		my $new_start_hit1 =  $strand eq "+" ? $end5_tsd + 1 : $new_start_hit;
		my $new_end_hit1   = $strand eq "+" ? $new_end_hit : $start5_tsd - 1;
		( $file_hit, $sequence_hit )           = &get_sequence( $id_dir, $query_id . "_hit",$new_start_hit1, $new_end_hit1, $seqobj, $strand, $width ) ;
		system("rm -f $file_hit") if -e $file_hit;
		#$hit_seq .= " + (Intervale TSD5 - Hit : ".$new_start_hit1 ."..".$new_end_hit1.")" ;
	}
	if ($start_polyA > 0 ) {
		my $start_between_hit_polyA  = $strand eq "+" ? $new_end_hit + 1 : $end_polyA + 1;
		my $end_between_hit_polyA    = $strand eq "+" ? $start_polyA  - 1 : $new_start_hit - 1;
		if ($end_between_hit_polyA - $start_between_hit_polyA > 0){
			( $file_between_hit_polyA, $sequence_between_hit_polyA )           = &get_sequence( $id_dir, $query_id . "_hit",$start_between_hit_polyA, $end_between_hit_polyA, $seqobj, $strand, $width );		
			system("rm -f $file_between_hit_polyA") if -e $file_between_hit_polyA;
			#$hit_seq .= " + (Hit - PolyA ".$start_between_hit_polyA."..".$end_between_hit_polyA.")";
		}
		my $start_between_polyA_tsd  = $strand eq "+" ? $end_polyA + 1 : $end3_tsd + 1;
		my $end_between_polyA_tsd    = $strand eq "+" ? $start3_tsd - 1 : $start_polyA -1;		
		if ($end_between_polyA_tsd - $start_between_polyA_tsd > 0) {
			#$hit_seq .= " + (PolyA - TSD : ".$start_between_polyA_tsd ."..".$end_between_polyA_tsd .")";
			( $file_between_polyA_tsd ,$sequence_between_polyA_tsd )           = &get_sequence( $id_dir, $query_id . "_hit",$start_between_polyA_tsd, $end_between_polyA_tsd, $seqobj, $strand, $width ) ; 
			system("rm -f $file_between_polyA_tsd") if -e $file_between_polyA_tsd;
		}
	}
	else {
		my $start_between_hit_tsd  = $strand eq "+" ? $new_end_hit + 1 : $end3_tsd + 1;
		my $end_between_hit_tsd    = $strand eq "+" ? $start3_tsd  - 1 : $new_start_hit - 1;
		if ($end_between_hit_tsd - $start_between_hit_tsd > 0){
			( $file_between_hit_tsd, $sequence_between_hit_tsd )           = &get_sequence( $id_dir, $query_id . "_hit",$start_between_hit_tsd, $end_between_hit_tsd, $seqobj, $strand, $width );		
			system("rm -f $file_between_hit_tsd") if -e $file_between_hit_tsd;
			#$hit_seq .= " + (Hit - TSD ".$start_between_hit_tsd."..".$end_between_hit_tsd.")";
		}
	}	
}
else {
	my $start_extend_before = $strand eq "+" ? $hsp_start_hit - $extend_sequence : $hsp_end_hit + 1;
	my $end_extend_before   = $strand eq "+" ? $hsp_start_hit - 1 : $hsp_end_hit + $extend_sequence;
	($file_sequence_before,$sequence_before) = &get_sequence( $id_dir, $query_id . "_hit",$start_extend_before, $end_extend_before , $seqobj, $strand, $width );
	system("rm -f $file_sequence_before") if -e $file_sequence_before;
	#$hit_seq .= " (Extend5 :".$start_extend_before."..".$end_extend_before.")" ;
	
	if ($start_polyA > 0 ) {
		my $start_extend_polyA = $strand eq "+" ? $end_polyA + 1 : $start_polyA - $extend_sequence;
		my $end_extend_polyA   = $strand eq "+" ? $end_polyA + $extend_sequence : $start_polyA - 1;
		($file_sequence_after,$sequence_after) = &get_sequence( $id_dir, $query_id . "_hit",$start_extend_polyA, $end_extend_polyA , $seqobj, $strand, $width );
		system("rm -f $file_sequence_after") if -e $file_sequence_after;
		my $start_between_hit_polyA  = $strand eq "+" ? $new_end_hit + 1 : $end_polyA + 1;
		my $end_between_hit_polyA    = $strand eq "+" ? $start_polyA  - 1 : $new_start_hit - 1;
		if ($end_between_hit_polyA - $start_between_hit_polyA > 0) {
		
			( $file_between_hit_polyA, $sequence_between_hit_polyA )           = &get_sequence( $id_dir, $query_id . "_hit",$start_between_hit_polyA, $end_between_hit_polyA, $seqobj, $strand, $width );		
			system("rm -f $file_between_hit_polyA") if -e $file_between_hit_polyA;	
		#	$hit_seq .= " + (Hit - PolyA : ".$start_between_hit_polyA."..".$end_between_hit_polyA.") ";
		}
		#$hit_seq .= " + (Extend3 : ".$start_extend_polyA."..".$end_extend_polyA.")" ;
	
	} 
	else {
		my $start_extend_polyA = $strand eq "+" ? $new_end_hit + 1 : $new_start_hit - $extend_sequence;
		my $end_extend_polyA   = $strand eq "+" ? $new_end_hit + $extend_sequence : $new_start_hit - 1;
		($file_sequence_after,$sequence_after) = &get_sequence( $id_dir, $query_id . "_hit",$start_extend_polyA, $end_extend_polyA , $seqobj, $strand, $width );
		system("rm -f $file_sequence_after") if -e $file_sequence_after;
		#$hit_seq .= " + (Extend3 : ".$start_extend_polyA."..".$end_extend_polyA.")" ;
	}
} 
my $class = "";
my $cptN = &count($sequence_complete,"N"); 
if ( $found_tsd == 0 ) { 
	if ( $hsp_end_query + $threshold_tsd > $query_length ) {
		if ( $hsp_start_query == 1 ) {
			if ( $found_repeat == 1 ) { 
				$class = $cptN > 5 ? "Gaps" : "Repeat";
			}
			else { 
				if ( $found_polyA == 1 ) {
				    if ( $interval_hit_polyA > $threshold_hit_polyA ) {
						$class = $cptN > 5 ? "Gaps" : "Alone";
				    }
				    else {
						$class = $cptN > 5 ? "Gaps" : "PolyA";
				    }
				}
				else {
				    $class = $cptN > 5 ? "Gaps" :"Alone";
				}
			}
		}
		else {
			$class = $cptN > 5 ? "Gaps" : "To Check";
		}
	} 
	else  {
		$class = $cptN > 5 ? "Gaps" : "To Check";
	}
}
else {
	if ( $hsp_end_query < $query_length - $threshold_tsd ) {
		if ( $hsp_start_query > 2 ) { 
			$class = $cptN > 5 ? "Gaps" :"To Check";
		}
		else {
			if ( $found_repeat == 1 ) { 
				$class = $cptN > 5 ? "Gaps" : "Repeat";	
			}
			else {
				if ( $found_polyA == 1 ) {
				    if ( $interval_hit_polyA > $threshold_hit_polyA ) {
						$class = $cptN > 5 ? "Gaps" : "To Check";
				    }
				    else {
						$class = $cptN > 5 ? "Gaps" : "PolyA";
					}
				}
				else {
					$class = $cptN > 5 ? "Gaps" : "3' truncated";
				}
			}
		}
	}
	else {
		if ( $found_repeat == 1 ) { 
			$class = $cptN > 5 ? "Gaps" : "Repeat";	
		}
		else {
			if ( $found_polyA == 1 ) {
				if ( $interval_hit_polyA > $threshold_hit_polyA ) {
					$class = $cptN > 5 ? "Gaps" : "Repeat";
				}
				else {
					$class = $cptN > 5 ? "Gaps" : "PolyA";
				}
			}
			else {
				$class = $cptN > 5 ? "Gaps" : "Alone";
			}
		}
	}
}  

system("rm -Rf $id_dir") if -d $id_dir;  
 
my $tsd_status = $found_tsd == 1 ? "TSD found" : "TSD not found";

my $json;

#$sequence_five .= " (TSD5 : ".$start5_tsd."..".$end5_tsd.")" if $found_tsd==1;
#$subseq_polyA.= " (PolyA : ".$start_polyA."..".$end_polyA.")" if $start_polyA > 1;

print TXT join( "\t" , $hit_name, $hsp_start_hit, $hsp_end_hit,$percent_identity,$gaps,$hit_seq,$strand,$class,    $query_id,  $query_name,  $hsp_start_query,  $hsp_end_query, $query_length, $hsp_length,  $algorithm, $found_tsd,  $start5_tsd, $end5_tsd,  $start3_tsd,  $end3_tsd,  $sequence_five,  $sequence_three,$repeat_size,  $tsd_overlap5, $tsd_overlap3, $cleavage,  $found_repeat, $found_polyA, $start_polyA,$end_polyA, $length_polyA,  $subseq_polyA,  $file_url_fasta,$sequence_complete,$sequence_hit, $sequence_repeat_polyA,$sequence_before,$sequence_after , $sequence_between_hit_polyA,$sequence_between_polyA_tsd,$sequence_between_hit_tsd, join(",",@repeat)),"\n";
close TXT;
$json->{"id"} = $cpt_hsp;

$json->{"file_url_txt"}   	 = $file_url_txt;
$json->{"file_url_fasta"}    = $file_url_fasta;
$json->{"file_url_xls"}   	 = $file_url_xls;
$json->{"file_url_tsd_txt"}  = $file_url_tsd_txt; 

$json->{"result"}= {
					hit_name     => $hit_name,
					hit_start    => $hsp_start_hit,
					hit_end      => $hsp_end_hit,
					hit_seq      => $hit_seq,
					sequence_hit => $sequence_hit,
					strand       => $strand,
					class        => $class,
					query_id     => $query_id,
					query_name   => $query_name,
					query_start  => $hsp_start_query,
					query_end    => $hsp_end_query,
					query_length => $query_length,
					hsp_length   => $hsp_length,
					algorithm    => $algorithm,
					tsd          => $tsd_status,
					tsd5_start   => $start5_tsd,
					tsd5_end     => $end5_tsd,
					tsd3_start   => $start3_tsd,
					tsd3_end     => $end3_tsd,
					tsd5_seq     => $sequence_five,
					tsd3_seq     => $sequence_three,
					repeat_size  => $repeat_size,
					tsd_overlap5 => $tsd_overlap5,
					tsd_overlap3 => $tsd_overlap3,
					cleavage     => $cleavage, 
					repeat       => $found_repeat,
					polyA        => $found_polyA,
					start_polyA  => $start_polyA,
					end_polyA    => $end_polyA,
					length_polyA => $length_polyA,
					seq_polyA    => $subseq_polyA, 
					seq_complete => $sequence_complete,
					repeat_hsp   => \@repeat,
					sequence_repeat=>$sequence_repeat_polyA,
					sequence_before => $sequence_before,
					sequence_between_hit_polyA => $sequence_between_hit_polyA,
					sequence_after => $sequence_after,
					sequence_between_polyA_tsd=>$sequence_between_polyA_tsd,
					sequence_between_hit_tsd=>$sequence_between_hit_tsd
					};

my $json_val = to_json($json );
print $json_val;
sub get_index {
    my ( $id, $database ) = @_;
    $ENV{BIOPERL_INDEX_TYPE} = "SDBM_File";
    $ENV{BIOPERL_INDEX}      = "index";
    my $file_index = $database . ".inx"; 
    my $index = Bio::Index::Fasta->new( -filename => $file_index  );
    my $seqobj = $index->fetch($id);
    return $seqobj;
}
 
sub get_sequence {
    my ( $dir, $id, $start, $end, $seqobj, $strand, $width ) = @_;
    my $seqout = $seqobj->trunc( $start, $end );
    if ( $strand eq "-" ) {
		$seqout = $seqout->revcom();
    }
    my $file = $dir . "/" . $id . ".fna";
    my $out  = new Bio::SeqIO(
			      -format => 'fasta',
			      -file   => ">$file",
			      -width  => $width
			      
			      );
    $seqout->display_id($id);
    $seqout->desc("");
    $out->write_seq($seqout);
    $out->close;
    return ( $file, $seqout->seq() );
}
 
sub run_blast {
    my ( $query, $database, $output, $evalue , $other ) = @_;
    #$evalue = 10 unless $evalue;
    if ( $BLASTALL =~ /blastall/ ) {   
   		system("$BLASTALL -i $query -d $database -a $thread -p $program -o $output -m 8 -e $evalue $other 2>/dev/null 1>/dev/null");
     }
    else {
		system("$BLASTALL -query $query -db $database -task $program -num_threads $thread -out $output -outfmt 6 -evalue $evalue") unless -e $output;
    }
    return $output;
}



sub wordfinder {
    my ( $file_5, $file_3, $id ) = @_;
    my $wordfinder_dir = $output_dir . "/wordfinder";
    my $output         = $wordfinder_dir . "/" . $id . ".wordfinder";
    my $file_error     = $wordfinder_dir . "/" . $id . ".error";
    system( "$WORDFINDER -asequence $file_5 -bsequence $file_3 -outfile $output -gapopen 30 -gapextend 2 -width $tsd_length -errorfile $file_error  2>/dev/null" );
  
    open( IN, $output );
    my @line;
    my $cpt = 0;
    my $start_five_tsd  = 0;
    my $end_five_tsd    = 0;
    my $start_three_tsd = 0; 
    my $end_three_tsd   = 0;
    my $mismatch = 0;
    my $length;
   	local $/ = "\n";
   	local $_; 
    while ( my $line = <IN> ) {
		chomp($line);
		next if $line =~ /^#/; 
		my @data = ( split( /\s+/, $line ) );
		if ( scalar(@data) == 4 ) {
	    	if ( $cpt == 0 ) {
				$cpt++;
				if ($start_five_tsd) {
		    		$start_five_tsd = $data[1] if $data[1] < $start_five_tsd;
		    		$end_five_tsd   = $data[3] if $data[3] > $end_five_tsd;
				}
				else {
		    		$start_five_tsd = $data[1];
		    		$end_five_tsd   = $data[3];
				}
	    	}
	    	else {
				$cpt = 0;
				if ($start_three_tsd) {
		    		$start_three_tsd = $data[1] if $data[1] < $start_three_tsd;
		    		$end_three_tsd   = $data[3] if $data[3] > $end_three_tsd;
				}
				else {
		   			$start_three_tsd = $data[1];
		    		$end_three_tsd   = $data[3];
				}
	    	}
		}
		else {
	    	if ( $line ne "" ) {
				$line =~ s/ //g;
				my @data = ( split( //, $line ) );
				for ( my $i = 0 ; $i <= $#data ; $i++ ) {
		    		if ( $data[$i] eq "." ) {
						$mismatch++;
		    		}
				}
	    	}
		}
    } 
    close IN;
    if ( 30 - $end_five_tsd < 3 && $mismatch < 3 ) {
		return ( $start_five_tsd, $end_five_tsd, $start_three_tsd, $end_three_tsd,  $length, $mismatch );
    }
    else {
		return -1;
    }
}

sub wordmatch {
    my ( $file_5, $file_3, $id ) = @_;
   	local $/ = "\n";
   	local $_;
    my $wordmatch_dir = $output_dir . "/wordmatch";
    my $output        = $wordmatch_dir . "/" . $id . ".wordmatch";
    my $file_5_gff    = $wordmatch_dir . "/" . $id . "_5.gff";
    my $file_3_gff    = $wordmatch_dir . "/" . $id . "_3.gff";
    my $file_log      = $wordmatch_dir . "/" . $id . ".log";
    my $min_location  = $extend_hit_region - $interval;  
    system("$WORDMATCH -asequence $file_5 -bsequence $file_3 -wordsize $tsd_length -outfile $output -aoutfeat $file_5_gff -boutfeat $file_3_gff  -logfile $file_log 2>/dev/null" ) ;
      system("rm -f $file_5_gff") if -e $file_5_gff;
    system("rm -f $file_3_gff") if -e $file_3_gff;
    system("rm -f $file_log") if -e $file_log;
    open( IN, $output );
    my $current_start_five_prime_tsd  = 0;
    my $current_end_five_prime_tsd    = 0;
    my $current_start_three_prime_tsd = 0;
    my $current_end_three_prime_tsd   = 0;
    my $current_length                = 0; 
    while ( my $line = <IN> ) {
		chomp($line);
		next if $line =~ /^#/;  
		if ( $line =~ /^\s+(\d*.*)/ ) { 
	    	my ( $length, $id, $five_prime_tsd, $three_prime_tsd ) = ( split( /\s+/, $1 ) )[ 0, 1, 3, 6 ];
	    	my ( $start_five_prime_tsd, $end_five_prime_tsd )   = $five_prime_tsd =~ /(\d*)\.\.(\d*)/;
	    	my ( $start_three_prime_tsd, $end_three_prime_tsd ) = $three_prime_tsd =~ /(\d*)\.\.(\d*)/; 
		    if ($current_end_five_prime_tsd) {
				if ( $end_five_prime_tsd > $current_end_five_prime_tsd ) {
		    		$current_start_five_prime_tsd  = $start_five_prime_tsd;
		    		$current_end_five_prime_tsd    = $end_five_prime_tsd;
		    		$current_start_three_prime_tsd = $start_three_prime_tsd;
		    		$current_end_three_prime_tsd   = $end_three_prime_tsd;
		    		$current_length                = $length;
				}
	    	}
	    	else {
				$current_start_five_prime_tsd  = $start_five_prime_tsd;
				$current_end_five_prime_tsd    = $end_five_prime_tsd;
				$current_start_three_prime_tsd = $start_three_prime_tsd;
				$current_end_three_prime_tsd   = $end_three_prime_tsd;
				$current_length                = $length;
		
	    	}
		}
    }
    #exit;
    close IN; 
    if ( $current_end_five_prime_tsd >= $min_location ) { 
		return ($current_start_five_prime_tsd,  $current_end_five_prime_tsd,$current_start_three_prime_tsd, $current_end_three_prime_tsd,$current_length);
    }
    else { 
		return -1;
    }
}



sub polyA_search {
    my $sequence_polyA = shift;
    my $offset = 0;
    my $char = "A";
    my $num = 10;
    my $cutoff = 6;
    my %polyA;
    my $found_polyA = 0;
    my $nbA = &count($sequence_polyA,"A");
    my $start_polyA = 0;
    if ($nbA == length($sequence_polyA)) {
		$polyA{$start_polyA} = {
	    	length => length($sequence_polyA),
	    	subseq => $sequence_polyA
	    };
		$found_polyA = 1;
    }
    else { 
		my $cpt = 0;
		my $result = index($sequence_polyA, $char, $offset);  
		if ($result >= 0) {
	    	my @seq = (split(//,$sequence_polyA));
	    	my $start_polyA = $result;
	    	my $nucl = 0;
	    	for (my $i=$result;$i<=$#seq;$i++) {
				my $subseq = substr($sequence_polyA,$i,$num); 
				my $nbA = &count($subseq,$char);
				my $percentA = $nbA / length($subseq); 
				if ($nbA > $cutoff) {
		    		$start_polyA = $i if $i < $start_polyA;
		    		$nucl++;
		    		if ($i == $#seq){
						my $end = $nucl + $num if $nucl > 0;
						my $subseq = substr($sequence_polyA,$start_polyA,$end);	 
						if ($nucl > 0) {
			    			$polyA{$start_polyA} = {
								length => length($subseq),
								subseq => $subseq
							};
			    			$found_polyA = 1;
						}
		    		}
				}
				else {
		    		$cpt++;
		    		if ($nucl > 0) {
						my $end = $nucl + $num;
						my $subseq = substr($sequence_polyA,$start_polyA,$end);	 
						$polyA{$start_polyA} = {
			    			length => length($subseq),
			    			subseq => $subseq
			    		};
						$found_polyA = 1;
		    		}
		    		$nucl = 0;
		    		$start_polyA = $i;
				}
			}
		}
    }
    return ( \%polyA, $found_polyA );
}

sub count {
    my ($subseq,$nucl) = @_;
    my @nucl   = ( split( //, $subseq ) );
    my $cpt   = 0;
    for ( my $i = 0 ; $i <= $#nucl ; $i++ ) {
		if ( $nucl[$i] eq $nucl ) {
	    	$cpt++;
		}
    }
    return $cpt;
}


 
sub parse_blast_result {
    my $file = shift;
    my $in   = new Bio::SearchIO(
				 -format => 'blasttable',
				 -file   => $file
				 );
    my $result = $in->next_result();
    my @hsp;
    my @result;
    if ($result) {
		my ( $start_query, $end_query, $start_hit, $end_hit );
		while ( my $hit = $result->next_hit ) {
	    	while ( my $hsp = $hit->next_hsp ) {
				if ($start_query) {
		    		if ( $hsp->end('query') < $result[0]->{start} ) {
						push @result, {
			    			start => $hsp->start('query'),
			    			end   => $hsp->end('query')
			   			};
						@result = sort { $a->{start} <=> $b->{start} } @result;
						push @hsp, {
			   				hit => $hit->name(),
			    			hsp => $hsp
			    		};
		    		}
		    		if ( $hsp->start('query') > $result[-1]->{end} ) {
						push @result, {
			    			start => $hsp->start('query'),
			    			end   => $hsp->end('query')
			    		};
						@result = sort { $a->{start} <=> $b->{start} } @result;
						push @hsp, {
			    			hit => $hit->name(),
			    			hsp => $hsp
			    		};
		   			}
		    		else {
						for ( my $i = 0 ; $i < $#result ; $i++ ) {
			    			my $tmp_start = $result[$i]->{end};
			    			my $tmp_end   = $result[ $i + 1 ]->{start};
			    			if ( $hsp->start('query') > $tmp_start && $hsp->end('query') < $tmp_end ) {
								push @result,{
				    				start => $hsp->start('query'),
				    				end   => $hsp->end('query')
				    			};
								@result = sort { $a->{start} <=> $b->{start} } @result;
								push @hsp,  {
				    				hit => $hit->name(),
				    				hsp => $hsp
				    			};
			    			}
			    			else { next; }
						}
		    		}
				}
				else {
		    		$start_query = $hsp->start('query');
		    		$end_query   = $hsp->end('query');
		    		push @result, {
						start => $hsp->start('query'),
						end   => $hsp->end('query')
					};
		    		push @hsp, {
						hit => $hit->name(),
						hsp => $hsp
					};
				}
	    	}
		}
    }
    return \@hsp;
}


1;
