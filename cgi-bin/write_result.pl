#!/usr/bin/perl 
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
	    
use SVG;
use Tie::IxHash;

tie my %chr, "Tie::IxHash"
    or die "could not tie %chr";

our $CDHIT      = "/home/gdroc/bin/cdhit";
our $DEBUG      = 0;

my $param = new CGI; 
 

my $file_txt		  = $param->param("file_txt");
my $file_tsd_txt  	  = $param->param("file_tsd_txt");
my $file_fasta		  = $param->param("file_fasta");
my $file_xls		  = $param->param("file_xls");
my $file_url_txt	  = $param->param("file_url_txt");
my $file_url_fasta	  = $param->param("file_url_fasta");
my $file_url_xls	  = $param->param("file_url_xls");
my $file_url_tsd_txt  = $param->param("file_url_tsd_txt");
my $species           = $param->param("species");
my $release           = $param->param("release");
my $width_seq             = $param->param("width");
my $num_hsps          = $param->param("num_hsps");
my $min_identity      = $param->param("min_identity");
my $min_hsp_length    = $param->param("min_hsp_length"); 
my $file_svg    	  = $param->param("file_svg");
my $file_url_svg      = $param->param("file_url_svg");
my $file_chr_length   = $param->param("file_chr_length");
my $num_hits          = $param->param("num_hits");
 
my $cp = $file_txt .".bk";

open( TXT, $file_txt ); 
system("cp $file_txt $cp");

my $file_complete = $file_fasta .".bk";
my $out = new Bio::SeqIO(
	-file => ">$file_complete",
	-format=> "fasta",
	-width => $width_seq
	);
open(FASTA,">$file_fasta");
print $param->header('application/json');

my $cpt_tsd = 0; 
my %count;
my %count_tsd;
my $count_only_tsd = 0;
my %rt_result;
my %seq_id;
my $line = 0;
my $file_before =  $file_fasta	.".verif";
open(OUT,">$file_before");
while(<TXT>) {
	chomp; 
	$line++;
	my ($hit_name, $hsp_start_hit, $hsp_end_hit,$percent_identity, $gaps,$hit_seq,$strand,$class,   $query_id,  $query_name,  $hsp_start_query,  $hsp_end_query, $query_length, $hsp_length,  $algorithm, $found_tsd,  $start5_tsd, $end5_tsd,  $start3_tsd,  $end3_tsd,  $sequence_five,  $sequence_three,$repeat_size,  $tsd_overlap5, $tsd_overlap3, $cleavage,  $found_repeat, $found_polyA, $start_polyA,$end_polyA, $length_polyA,  $subseq_polyA,  $file_url_fasta, $sequence_complete,$sequence_hit,$sequence_repeat,$sequence_before,$sequence_after , $sequence_between_hit_polyA,$sequence_between_polyA_tsd,$sequence_between_hit_tsd,$repeat) = (split(/\t/,$_));
	$cpt_tsd += $found_tsd;
	$rt_result{$query_id} = {
					hit_name     => $hit_name,
					hit_start    => $hsp_start_hit,
					hit_end      => $hsp_end_hit,
					hit_seq      => $hit_seq,
					strand       => $strand,
					class        => $class, 
					query_id     => $query_id,
					query_name   => $query_name,
					query_start  => $hsp_start_query,
					query_end    => $hsp_end_query,
					query_length => $query_length,
					hsp_length   => $hsp_length,
					percent_identity=> $percent_identity,
					algorithm    => $algorithm,
					tsd          => $found_tsd,
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
					gaps         => $gaps,
					repeat       => $found_repeat,
					polyA        => $found_polyA,
					start_polyA  => $start_polyA,
					end_polyA    => $end_polyA,
					length_polyA => $length_polyA,
					seq_polyA    => $subseq_polyA,
					seq_repeat    =>$sequence_repeat,
					sequence_hit => $sequence_hit,
					sequence_complete => $sequence_complete,
					repeat_hsp   => $repeat,
					sequence_before => $sequence_before,
					sequence_between_hit_polyA => $sequence_between_hit_polyA,
					sequence_after => $sequence_after,
					sequence_between_polyA_tsd=>$sequence_between_polyA_tsd,
					sequence_between_hit_tsd=>$sequence_between_hit_tsd
				}; 
	$count{$class}++;
	$count_tsd{$class}{$found_tsd}++;
	$count_only_tsd++ if $found_tsd == 1;
	print OUT ">$query_id\n$sequence_before\n";
	push @{$seq_id{$class}},$query_id;
	
}
close TXT;
close OUT;
  
open( TXT, ">$file_txt" ); 
open( TXT_TSD, ">$file_tsd_txt" ); 
my $duplication = &run_cdhit($file_before);
my %duplication1 = %$duplication;
my %remove_class; 
foreach my $dupli ( keys %duplication1 ) {
	$remove_class{ $rt_result{$dupli}{class} }{$rt_result{$dupli}{tsd}}++;  
	foreach my $similar ( keys %{$duplication1{$dupli}}) { 
		$duplication{$dupli} = $similar;
	}
}    
my @species_to_remove = ("microcebus_murinus");

my %species_to_remove = map { $_ => 1 } @species_to_remove;
 
open(IN,$file_chr_length);
my $height = 0; 
my $count_chr = 0;
while(<IN>){
    chomp;
    $count_chr++;
    my ($label,$length) = (split(/\t/,$_));
    $chr{$label} = {
    					length => $length,
    					cpt    => $count_chr
    					};
    $height = $length if $length > $height;
}
close IN; 
$height = int($height / 100000) + 50;
$height =~ s/,/./g;
my $nb_chr = keys(%chr); 
open(OUT,">$file_svg");

my $width = ($nb_chr + 1) * 50;
my $svg = SVG->new(
		   width  => $width,
		   height => $height
		   );
my $defs = $svg->defs(id => 'defs');
my $gradient1 = $defs->gradient(
				-type => "linear",
				id    => "gradient1"
				);
my $stop1 = $gradient1->stop(
			     id => 'stop1',
			     offset=>0,
			     style => {
				 'stop-color'=>'#a1b0c8',
				 'stop-opacity'=> 1
				 });
my $stop2 = $gradient1->stop(
			     id     => 'stop2',
			     offset => 1,
			     style  => {
				 'stop-color'=>'#a1b0c8',
				 'stop-opacity'=> 0
				 });
my $gradient2 = $defs->gradient(
				-type => 'linear',
				id    => 'gradient2'
				);
my $stop3 = $gradient2->stop(
			     id => 'stop3',
			     offset=>0,
			     style => {
				 'stop-color'=>'#9ac3e9',
				 'stop-opacity'=> 1
				 });
my $gradient3 = $defs->gradient(  
				  -type => 'linear',
				  id    => 'gradient3'
				  );
my $stop4 = $gradient3->stop(
			     id => 'stop4',
			     offset=>0,
			     style => {
				 'stop-color'=>'#2696ff',
				 'stop-opacity'=> 1
				 });
    
my $cpt = 0;
foreach my $label (keys %chr) {
    $cpt++;
    my $x = $cpt * 50;
    my $x2 = $x + 15;
    my $length = $chr{$label}{length}; 
    my $height = $length / 100000;
    
    $height =~ s/,/./g;
    my $gradient_id = join("_","gradient",$label);
    my $text_id = join("_","text",$label);
    my $gradient = $defs->gradient(
                                    -type => 'linear',
                                    id    => $gradient_id,
                                    x1 => $x,
                                    y1 => '158.5618',
                                    x2 => $x2,
                                    y2 => '158.5618',
                                    'xlink:href' => '#gradient1',
                                    'gradientUnits' => 'userSpaceOnUse'
                                    );
    
    my $chr = $svg->rectangle(
			      x     => $x, 
			      y      => 20,
			      width => 15,
			      height => $height,
			      rx    => 5, 
			      ry     => 5,
			      id    => $label,
			      style => {
				  'opacity' => 0.6,
				  'fill' => 'url(#'.$gradient_id.')',
				  'fill-opacity' => 1,
				  'fill-rule' => 'evenodd',
				  'stroke' => '#000000',
				  'stroke-width' => '0.80000001',
				  'stroke-linecap' => 'round',
				  'stroke-linejoin' => 'round',
				  'stroke-dasharray' => 'none'
				  }
			      );
    my $text = $svg->text(
                           id=> $text_id, 
                           x=> $x - 5,
                           y=> '15',   
                           style=> {
                               'font-size' => '10px',
                               'font-style' => 'normal',
                               'font-variant' => 'normal',
                               'font-weight' => 'normal',
                               'font-stretch' => 'normal',
                               'text-align' => 'start',
                               'line-height' => '125%',
                               'writing-mode' => 'lr-tb',
                               'text-anchor' => 'start',
                               'fill' => '#000000',
                               'font-family' => 'Lucida Sans',
                               '-inkscape-font-specification' => 'Lucida Sans'
                               }
			  )->cdata($label);
}

 
my $workbook = Spreadsheet::WriteExcel->new( $file_xls ); 
my $worksheet1 = $workbook->add_worksheet("result");
my $worksheet2 = $workbook->add_worksheet("cleavage");
my $worksheet3 = $workbook->add_worksheet("report");
my $bold = $workbook->add_format(); 
$bold->set_bold();  
my $red    = $workbook->add_format();   # T
my $green  = $workbook->add_format();   # A
my $orange = $workbook->add_format();   # G
my $blue   = $workbook->add_format();   # C
$red->set_bg_color( 'red' );
$green->set_bg_color( 'green' );
$orange->set_bg_color( 'orange' );
$blue->set_bg_color( 'blue' );
my $row  = 0;
my $row1 = 0;
$worksheet2->write( 0,0, 'Query id', $bold );
$worksheet2->write( 0,1, 'Class', $bold ); 
my $num_hsp         = scalar( keys %rt_result );
my $num_duplication = scalar( keys %duplication );
my $num_gaps        = $count{'Gaps'} ? $count{'Gaps'} : 0;
my $num_total       = $num_hits - $num_duplication - $num_gaps;
    
print TXT "Species : " . $species, "\n";
print TXT "Release : " . $release, "\n";
print TXT "Number of hsps total = $num_hsps\n";
print TXT "Number of hsps with at least ", $min_identity, "% ID and ", $min_hsp_length, " bp : ",$num_hits, "\n";
print TXT "Number of hsps with at least ", $min_identity, "% ID and TSD found ", $count_only_tsd, "\n\n" if $count_only_tsd;
if ($count_only_tsd) {
	print TXT_TSD "Species : " . $species, "\n";
	print TXT_TSD "Release : " . $release, "\n";
	print TXT_TSD "Number of hsps total = $num_hsps\n";
	print TXT_TSD "Number of hsps with at least ", $min_identity, "% ID and ", $min_hsp_length, " bp : ", $num_hits, "\n";
	print TXT_TSD "Number of hsps with at least ", $min_identity, "% ID and TSD found ", $count_only_tsd, "\n\n";
}
$worksheet1->write( $row1,0, 'Species', $bold );
$worksheet1->write( $row1,1, $species ); 
$row1++;
$worksheet1->write( $row1,0, 'Release', $bold );
$worksheet1->write( $row1,1, $release ); 
$row1++;
$worksheet1->write( $row1,0, 'Number of hsps total', $bold );
$worksheet1->write( $row1,1, $num_hsps ); 
$row1++;
$worksheet1->write( $row1,0, 'Number of hsps with at least '. $min_identity. ' % ID and '. $min_hsp_length. ' bp', $bold );
$worksheet1->write( $row1,1, $num_hits ); 
if ( $count_only_tsd > 0) {
	$row1++;
	$worksheet1->write( $row1,0, 'Number of hsps with at least '. $min_identity . '% ID and TSD found', $bold );
	$worksheet1->write( $row1,1, $count_only_tsd );
} 
$row1+=2;
my $total_hit_with_tsd = 0;
foreach my $class ( "Alone","Repeat","PolyA","3' truncated"  ) { 
	my $nb_class = $count{$class} ? $count{$class} : 0; 
	if ( defined $remove_class{$class} ) { 
	    if ($remove_class{$class}{1}) {
			$nb_class          =  $nb_class - $remove_class{$class}{1};
	    }
	    if ($remove_class{$class}{0}) {
			$nb_class          =  $nb_class - $remove_class{$class}{0};
	    }
	    if ($remove_class{$class}{1}) {
			$nb_class_with_tsd = $count_tsd{$class}{1} ? $count_tsd{$class}{1} - $remove_class{$class}{1} : 0;
	    }
	    else {
			$nb_class_with_tsd = $count_tsd{$class}{1} ? $count_tsd{$class}{1} : 0;
	    }
	    if ($remove_class{$class}{0}) {
			$nb_class_without_tsd = $count_tsd{$class}{0} ? $count_tsd{$class}{0} - $remove_class{$class}{0} : 0;			
	    }
	    else {
			$nb_class_without_tsd = $count_tsd{$class}{0} ? $count_tsd{$class}{0}  : 0;
	    }
	}
	else { 
	    $nb_class_with_tsd = $count_tsd{$class}{1} ? $count_tsd{$class}{1} : 0;
	    $nb_class_without_tsd = $count_tsd{$class}{0} ? $count_tsd{$class}{0} : 0;
	}
	$count{$class}{"with"} = $nb_class;
	$count{$class}{"without"} = $nb_class_with_tsd;
	$total_hit  += $nb_class;
	$total_hit_with_tsd += $nb_class_with_tsd;
	
}
 

$worksheet1->write( $row1,0, 'Species', $bold );
$worksheet1->write( $row1,1, 'Total hits', $bold );
$worksheet1->write( $row1,2, 'Total hits > %identity',$bold );
$worksheet1->write( $row1,3, 'Total hits analyzed',$bold );
$worksheet1->write( $row1,4, 'Total hits analyzed with tsd', $bold );
$worksheet1->write( $row1,5, 'Alone', $bold );
$worksheet1->write( $row1,7, 'Repeat', $bold );
$worksheet1->write( $row1,9, 'PolyA', $bold );
$worksheet1->write( $row1,11, "3' trunc", $bold ); 
$row1++;
$worksheet1->write( $row1,0, $species );
$worksheet1->write( $row1,1, $num_hsps );
$worksheet1->write( $row1,2, $line ); #135
$worksheet1->write( $row1,3, $total_hit);
$worksheet1->write( $row1,4, $total_hit_with_tsd );

my $row3 = 0; 
$worksheet3->write( $row3,0, 'Species', $bold );
$worksheet3->write( $row3,1, 'Total hits', $bold );
$worksheet3->write( $row3,2, 'Total hits > %identity',$bold );
$worksheet3->write( $row3,3, 'Total hits analyzed',$bold );
$worksheet3->write( $row3,4, 'Total hits analyzed with tsd', $bold );
$worksheet3->write( $row3,5, 'Alone', $bold );
$worksheet3->write( $row3,7, 'Repeat', $bold );
$worksheet3->write( $row3,9, 'PolyA', $bold );
$worksheet3->write( $row3,11, "3' trunc", $bold ); 
$row3++;
$worksheet3->write( $row3,0, $species );
$worksheet3->write( $row3,1, $num_hsps );
$worksheet3->write( $row3,2, $line );
$worksheet3->write( $row3,3, $total_hit);
$worksheet3->write( $row3,4, $total_hit_with_tsd );



my $cpt_class = 4;
foreach my $class ( "Alone","Repeat","PolyA","3' truncated"  ) {
	$cpt_class++;
	$worksheet1->write( $row1,$cpt_class,$count{$class}{"with"} );
	$worksheet3->write( $row3,$cpt_class,$count{$class}{"with"}); 
	$cpt_class++; 
	$worksheet1->write( $row1,$cpt_class ,$count{$class}{"without"} );
	$worksheet3->write( $row3,$cpt_class , $count{$class}{"without"} );
}
$row1+=2;
$worksheet1->write( $row1,0, 'Class', $bold );
$worksheet1->write( $row1,1, 'Number', $bold );
$worksheet1->write( $row1,2, '%', $bold ); 
$row1++;  
my $json;
foreach my $class ( "Alone","Repeat","PolyA","3' truncated","To Check" ) {
	next if $class eq "Gaps"; 
	my $nb_class_with_tsd    = 0;
	my $nb_class_without_tsd = 0;
	my $nb_class  			 = $count{$class} ? $count{$class} : 0; 
	if ( defined $remove_class{$class} ) { 
	    if ($remove_class{$class}{1}) {
			$nb_class          =  $nb_class - $remove_class{$class}{1};
	    }
	    if ($remove_class{$class}{0}) {
			$nb_class          =  $nb_class - $remove_class{$class}{0};
	    }
	    if ($remove_class{$class}{1}) {
			$nb_class_with_tsd = $count_tsd{$class}{1} ? $count_tsd{$class}{1} - $remove_class{$class}{1} : 0;
	    }
	    else {
			$nb_class_with_tsd = $count_tsd{$class}{1} ? $count_tsd{$class}{1} : 0;
	    }
	    if ($remove_class{$class}{0}) {
			$nb_class_without_tsd = $count_tsd{$class}{0} ? $count_tsd{$class}{0} - $remove_class{$class}{0} : 0;			
	    }
	    else {
			$nb_class_without_tsd = $count_tsd{$class}{0} ? $count_tsd{$class}{0}  : 0;
	    }
	}
	else { 
	    $nb_class_with_tsd = $count_tsd{$class}{1} ? $count_tsd{$class}{1} : 0;
	    $nb_class_without_tsd = $count_tsd{$class}{0} ? $count_tsd{$class}{0} : 0;
	}
	my $percent = sprintf( "%.2f", ( $nb_class / $line ) * 100 );
	$worksheet1->write( $row1,0,$class );
	$worksheet1->write( $row1,1, $nb_class );
	$worksheet1->write( $row1,2, $percent);
	$row1++;
	$worksheet1->write( $row1,0,"");
	$worksheet1->write( $row1,1, "w tsd" );
	$worksheet1->write( $row1,2, $nb_class_with_tsd);
	$row1++;
	$worksheet1->write( $row1,0,"");
	$worksheet1->write( $row1,1, "w/o tsd");
	$worksheet1->write( $row1,2, $nb_class_without_tsd);
	$row1+=2;  
}
    
$worksheet1->write( $row1,0,"Total",$bold);
$worksheet1->write( $row1,1,$num_total);
$row1++;
$worksheet1->write( $row1,0,"Gaps",$bold);
$worksheet1->write( $row1,1,$num_gaps);
$row1++;
$worksheet1->write( $row1,0,"Duplication",$bold);
$worksheet1->write( $row1,1,$num_duplication);
$row1+=2;	
$worksheet1->write($row1,0,"Query id",$bold);
$worksheet1->write($row1,1,"Class",$bold);
$worksheet1->write($row1,2,"% Identity",$bold);
$worksheet1->write($row1,3,"Query start",$bold);
$worksheet1->write($row1,4,"Query end",$bold);
$worksheet1->write($row1,5,"TSD Overlap 5'",$bold);
$worksheet1->write($row1,6,"TSD Overlap 3'",$bold);
$worksheet1->write($row1,7,"TSD Status",$bold);
$worksheet1->write($row1,8,"TSD length",$bold);
$worksheet1->write($row1,9,"Repeat Length",$bold);
$worksheet1->write($row1,10,"PolyA",$bold);
$row1++;  
my $cpt_hsp = 0;  
my $cpt_hsp_tsd = 0; 
my %sequence_by_class;
foreach my $query_id ( sort { $a cmp $b } keys %rt_result ) {
	if ($query_id ) {
		$cpt_hsp++; 
		print TXT "\n***********\n";
		print TXT $rt_result{$query_id}{query_id} . "\n"; 
		my $space = " ";
    	$space x= 50;
    	my $sequence_before = $rt_result{$query_id}{sequence_before};
   	 	my $sequence_after = $rt_result{$query_id}{sequence_after}; 
		my @seq;
		push @seq , $sequence_before;
		if ($rt_result{$query_id}{tsd5_seq}){
			push @seq, $rt_result{$query_id}{tsd5_seq};
		}
		if ($rt_result{$query_id}{sequence_hit}){
			push @seq, $rt_result{$query_id}{sequence_hit};
		}
		if ($rt_result{$query_id}{sequence_between_hit_polyA}){
			push @seq, $rt_result{$query_id}{sequence_between_hit_polyA};
		}	
		if ($rt_result{$query_id}{seq_polyA}) {
			push @seq, $rt_result{$query_id}{seq_polyA};
		}
		if ($rt_result{$query_id}{sequence_between_polyA_tsd}){
			push @seq, $rt_result{$query_id}{sequence_between_polyA_tsd};
		}
		if ($rt_result{$query_id}{sequence_between_hit_tsd}){
			push @seq, $rt_result{$query_id}{sequence_between_hit_tsd};
		}
		if ($rt_result{$query_id}{tsd3_seq}){
			push @seq, $rt_result{$query_id}{tsd3_seq};
		}
		push @seq,   $sequence_after;
		my $sequence  = join($space,@seq);
		$sequence_complete = &format_sequence($sequence ,$width_seq); 
		print FASTA  ">$query_id\n$sequence_complete\n" unless ( $duplication{$query_id}); 
		
		my $class1 = $rt_result{$query_id}{class};
		$class1 =~s/\s/\_/g;
		$class1 =~s/\'//g;
		$class1 =~s/3\_//;
		push @{$sequence_by_class{$class1}},  ">$query_id\n$sequence_complete\n" unless ( $duplication{$query_id}); 
		my $seqobj = new Bio::PrimarySeq(
			-id => $query_id,
			-seq => $rt_result{$query_id}{sequence_complete}
		);
		$out->write_seq($seqobj) unless ( $duplication{$query_id});
		my $length_tsd = $rt_result{$query_id}{tsd} == 1 ? length( $rt_result{$query_id}{tsd5_seq} ) : "na";
		print TXT "Query Name            = ", $rt_result{$query_id}{query_name}, " from ", $rt_result{$query_id}{query_start}, " to ", $rt_result{$query_id}{query_end}, "\n";
		print TXT "Hit Name              = ", $rt_result{$query_id}{hit_name}, " from ", $rt_result{$query_id}{hit_start}, " to ",$rt_result{$query_id}{hit_end}, "\n";
		if (  $duplication{$query_id} ) {
	    	my $dupli_query_id = $duplication{$query_id}  ;
	    	print TXT "Duplication with : ", $dupli_query_id, "\n";
		}
		my $class = $rt_result{$query_id}{class};
		$rt_result{$query_id}{class} = "Repeat" if $class eq "L1";
		my $tsd_status = $rt_result{$query_id}{tsd} == 1 ? "TSD found" : "TSD not found"; 
		print TXT "Strand                = ", $rt_result{$query_id}{strand}, "\n";
		print TXT "Class                 = ", $rt_result{$query_id}{class},  "\n";
		print TXT "HSP length            = ", $rt_result{$query_id}{hsp_length}, "\n";
		print TXT "$tsd_status\n";
		print TXT "Algorithm used        = ", $rt_result{$query_id}{algorithm}, "\n" if $rt_result{$query_id}{algorithm};
		print TXT "Match found in 5' TSD = ", $rt_result{$query_id}{tsd_overlap5}, "\n" if $rt_result{$query_id}{tsd_overlap5} > 0;
		print TXT "Match found in 3' TSD = ", $rt_result{$query_id}{tsd_overlap3}, "\n" if $rt_result{$query_id}{tsd_overlap3} > 0;
		my $location = $rt_result{$query_id}{hit_start} / 100000;
    	$location += 20;
		$location =~ s/,/./g;
    	my $x1 = ($chr{$rt_result{$query_id}{hit_name}}{cpt} * 50);
    	my $x2 = ($chr{$rt_result{$query_id}{hit_name}}{cpt} * 50) + 15;
    	my $line = $svg->anchor(
                            id      => $query_id
                            )->line(
                                    x1     => $x1, 
                                    y1     => $location,
                                    x2     => $x2,
                                    y2     => $location,
                                    style=> {
                                            'stroke' => 'rgb(255,0,0)',
                                            'stroke-width' => 2
                                            }
                                    );
  		my $desc = $line->desc(id => 'document'.$query_id)->cdata($query_id);
    	my $title = $line->title(  id => 'document1'.$query_id)->cdata($query_id);
		if ($rt_result{$query_id}{tsd} == 1) {
	    	$cpt_hsp_tsd++;
	    	print TXT_TSD "\n***********\n";
	    	print TXT_TSD $rt_result{$query_id}{query_id} . "\n";
	    	my $length_tsd = $rt_result{$query_id}{tsd} == 1 ? length( $rt_result{$query_id}{tsd5_seq} ) : "na";
	    	print TXT_TSD "Query Name            = ", $rt_result{$query_id}{query_name}, " from ", $rt_result{$query_id}{query_start}, " to ", $rt_result{$query_id}{query_end}, "\n";
	    	print TXT_TSD "Hit Name              = ", $rt_result{$query_id}{hit_name}, " from ", $rt_result{$query_id}{hit_start}, " to ",$rt_result{$query_id}{hit_end}, "\n";
	    	if (  $duplication{$query_id} ) {
	   			my $dupli_query_id = $duplication{$query_id}  ;
	    		print TXT_TSD "Duplication with : ", $dupli_query_id, "\n";
			}
	    	print TXT_TSD "Strand                = ", $rt_result{$query_id}{strand}, "\n";
	    	print TXT_TSD "Class                 = ", $rt_result{$query_id}{class},  "\n";
	    	print TXT_TSD "HSP length            = ", $rt_result{$query_id}{hsp_length}, "\n";
	    	print TXT_TSD "$tsd_status\n";
	    	print TXT_TSD "Algorithm used        = ", $rt_result{$query_id}{algorithm}, "\n" if $rt_result{$query_id}{algorithm};
	    	print TXT_TSD "Match found in 5' TSD = ", $rt_result{$query_id}{tsd_overlap5}, "\n" if $rt_result{$query_id}{tsd_overlap5} > 0;
	    	print TXT_TSD "Match found in 3' TSD = ", $rt_result{$query_id}{tsd_overlap3}, "\n" if $rt_result{$query_id}{tsd_overlap3} > 0;
		}
		my $tsd_overlap5 = $rt_result{$query_id}{tsd_overlap5} > 0 ? $rt_result{$query_id}{tsd_overlap5} : "";
		my $tsd_overlap3 = $rt_result{$query_id}{tsd_overlap3} > 0 ? $rt_result{$query_id}{tsd_overlap3} : "";
		if ( $rt_result{$query_id}{cleavage} ) { 
	    	print TXT "Cleavage              = ", $rt_result{$query_id}{cleavage}, "\n";
	    	print TXT_TSD "Cleavage              = ", $rt_result{$query_id}{cleavage}, "\n" if ($rt_result{$query_id}{tsd} == 1); 
	    	my ($five,$three) = (split(/\-/, $rt_result{$query_id}{cleavage}));
	    	$row++;
	    	$worksheet2->write( $row,0, $query_id );
	    	$worksheet2->write( $row,1, $rt_result{$query_id}{class} );
	    	my @five  = (split(//,$five));
	    	push @five ,"-";
	    	push @five, (split(//,$three));
	    	for (my $i=0;$i<=$#five;$i++) {
			my $color = "" ;
			my $nucl = $five[$i];
			if ($nucl eq "T") {
		    	$color = $red;
			}
			if ($nucl eq "A") {
		    	$color = $green;
			}
			if ($nucl eq "C") {
		    	$color = $blue;
			}
			if ($nucl eq "G") {
		    	$color = $orange;
			}
			$worksheet2->write( $row,$i+2, $nucl ,$color);	
	    }
	}
	my $repeat_size = $rt_result{$query_id}{repeat_size}  ? $rt_result{$query_id}{repeat_size} : 0;
	print TXT "Gaps                  = ", $rt_result{$query_id}{gaps}, "\n" if $rt_result{$query_id}{gaps};
	print TXT "Repeat Size                  = ", $repeat_size, "\n" if $repeat_size; 
	if ($rt_result{$query_id}{tsd} == 1) {
	    print TXT_TSD "Gaps                  = ", $rt_result{$query_id}{gaps}, "\n" if $rt_result{$query_id}{gaps};
	    print TXT_TSD "Repeat Size                  = ", $repeat_size, "\n" if $repeat_size; 
	}
	if ( $rt_result{$query_id}{tsd} == 1 ) {
	    print TXT join( "\t", "5TSD", $rt_result{$query_id}{tsd5_start}, $rt_result{$query_id}{tsd5_end}, length( $rt_result{$query_id}{tsd5_seq} ), $rt_result{$query_id}{tsd5_seq} ), "\n";
	    print TXT_TSD join( "\t", "5TSD", $rt_result{$query_id}{tsd5_start}, $rt_result{$query_id}{tsd5_end}, length( $rt_result{$query_id}{tsd5_seq} ), $rt_result{$query_id}{tsd5_seq} ), "\n" if ($rt_result{$query_id}{tsd} == 1);
	}
	print TXT join( "\t" ,"Hit", $rt_result{$query_id}{hit_start}, $rt_result{$query_id}{hit_end}, $rt_result{$query_id}{hit_seq} ), "\n";
	print TXT_TSD join( "\t" ,"Hit", $rt_result{$query_id}{hit_start}, $rt_result{$query_id}{hit_end}, $rt_result{$query_id}{hit_seq} ), "\n" if ($rt_result{$query_id}{tsd} == 1);
	my $repeat_length = 0;
	if ( $rt_result{$query_id}{repeat} == 1 ) {
	    my @repeat = @{ $rt_result{$query_id}{repeat_hsp} };
	    for ( my $i = 0 ; $i <= $#repeat ; $i++ ) {
			print TXT join( "\t", "Repeat", $repeat[$i]->{start_repeat}, $repeat[$i]->{end_repeat}, $repeat[$i]->{hit_name}, $repeat[$i]->{hsp_length} ), "\n";
			print TXT_TSD join( "\t", "Repeat", $repeat[$i]->{start_repeat}, $repeat[$i]->{end_repeat}, $repeat[$i]->{hit_name}, $repeat[$i]->{hsp_length} ), "\n" if ($rt_result{$query_id}{tsd} == 1);
			$repeat_length += $repeat[$i]->{hsp_length};
	    }
	}
	$repeat_length = "" if $repeat_length == 0;
	my $polyA_def = "na";
	if ( $rt_result{$query_id}{polyA} == 1 ) {
	    my $length_polyA = $rt_result{$query_id}{length_polyA};
	    my @seq_polyA    = ( split( //, $rt_result{$query_id}{seq_polyA} ) ); 
	    my $cptA         = 0;
	    for ( my $i = 0 ; $i <= $#seq_polyA ; $i++ ) {
			if ( $seq_polyA[$i] eq "A" ) {
		    	$cptA++;
			}
	    } 
	    $polyA_def = $length_polyA . " bp / " . $cptA . " As";
	    print TXT join( "\t", "PolyA", $rt_result{$query_id}{start_polyA}, $rt_result{$query_id}{end_polyA}, $rt_result{$query_id}{length_polyA}, $rt_result{$query_id}{seq_polyA} ), "\n";
	    print TXT_TSD join( "\t", "PolyA", $rt_result{$query_id}{start_polyA}, $rt_result{$query_id}{end_polyA}, $rt_result{$query_id}{length_polyA}, $rt_result{$query_id}{seq_polyA} ), "\n" if ($rt_result{$query_id}{tsd} == 1);
	}
	if ( $rt_result{$query_id}{tsd} == 1 ) {
	    print TXT join( "\t", "3TSD", $rt_result{$query_id}{tsd3_start}, $rt_result{$query_id}{tsd3_end}, length( $rt_result{$query_id}{tsd3_seq} ), $rt_result{$query_id}{tsd3_seq} ),  "\n";
	    print TXT_TSD join( "\t", "3TSD", $rt_result{$query_id}{tsd3_start}, $rt_result{$query_id}{tsd3_end}, length( $rt_result{$query_id}{tsd3_seq} ), $rt_result{$query_id}{tsd3_seq} ),  "\n" if ($rt_result{$query_id}{tsd} == 1);
	}
	$tsd_status = $tsd_status eq "TSD found" ? "y" : "n";
	$worksheet1->write($row1,0,$query_id);
	$worksheet1->write($row1,1,$class);
	$worksheet1->write($row1,2,$rt_result{$query_id}{percent_identity});
	$worksheet1->write($row1,3,$rt_result{$query_id}{query_start});
	$worksheet1->write($row1,4,$rt_result{$query_id}{query_end});
	$worksheet1->write($row1,5,$tsd_overlap5);
	$worksheet1->write($row1,6,$tsd_overlap3);
	$worksheet1->write($row1,7,$tsd_status);
	$worksheet1->write($row1,8,$length_tsd);
	$worksheet1->write($row1,9,$repeat_length);
	$worksheet1->write($row1,10,$polyA_def);
	$row1++; 
	}
} 
$workbook->close;

 
foreach my $class ( keys %sequence_by_class ) { 
	my $file_class = $file_fasta;
	$file_class =~ s/\.fna/\_$class\.fna/;
	$file_url_class = $file_url_fasta;
	$file_url_class =~ s/\.fna/\_$class\.fna/;
	$json->{$class}    = $file_url_class;
	open(FASTA,">$file_class"); 
	print FASTA join("\n",@{$sequence_by_class{$class}}),"\n";
	close FASTA;
}
close TXT;
close TXT_TSD;  
print OUT $svg->xmlify; 
#my $json;
$json->{"id"} = $cpt_hsp; 
$json->{"file_url_txt"}   	 = $file_url_txt;
$json->{"file_url_svg"}   	 = $file_url_svg;
$json->{"file_url_fasta"}    = $file_url_fasta;
$json->{"file_url_xls"}   	 = $file_url_xls;
$json->{"file_url_tsd_txt"}  = $file_url_tsd_txt;
$json->{"file_url_svg"}  = $file_url_svg;
$json->{"count_tsd"}  		 = $cpt_tsd;
my $json_val = to_json($json );
print $json_val;


sub format_sequence {
	my ($sequence, $lenght_line) = @_;
	my $length_seq  = length($sequence);
    my $nb          = int $length_seq / $lenght_line;
    my @formated_sequence;
    if ( $nb < 1 ) {
    	push @formated_sequence, substr( $sequence, 0, $lenght_line );
    }
    else {
        for ( my $j = 0 ; $j <= $nb ; $j++ ) {
            push @formated_sequence,
            substr( $sequence, $lenght_line * $j, $lenght_line );
        }
    }
    return ( join( "\n", @formated_sequence ) . "\n" );
}

sub run_cdhit {
    my $file           = shift;
    my $file_cdhit     = $file . ".cdhit";
    my $file_bak_clstr = $file_cdhit . ".bak.clstr";
    my $file_clstr     = $file_cdhit . ".clstr";
    system("$CDHIT -i $file -o $file_cdhit -T 1 -c 0.95 2>/dev/null 1>/dev/null"); 
    
    if (-e $file_bak_clstr) {
		open( IN, $file_bak_clstr );
    	my %exist_id; 
    
    	while (<IN>) {
			chomp;
			my ( $cpt, $value ) = ( split( /\t/, $_ ) );
			my $id;
			if ( $value =~ /\d+.*\>(.*)\.{3}\s.*/ ) {
	   			$id = $1; 
			}
			if ( $exist_id{$cpt} ) {
	    		$duplication{$id}{$exist_id{$cpt}}=1; 
			}
			else {
	    		$exist_id{$cpt} = $id;
			} 
    	}
    	close IN;
    	#system("mv $file_cdhit $file"); 
    	#system("rm $file_clstr $file_bak_clstr");
	} 
	return \%duplication;
}


1;
