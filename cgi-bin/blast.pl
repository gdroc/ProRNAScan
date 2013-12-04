#!/usr/bin/perl -w 


use lib '/home/gdroc/lib';
#use strict;
use FindBin;
use File::Temp qw/ tempfile tempdir /;
use CGI qw(:standard);use CGI::Carp qw(fatalsToBrowser);
use JSON;
use JSON::Any;
use lib "$FindBin::Bin/lib";
use warnings; 
use Bio::SeqIO;
use File::Basename;
use Bio::Index::Fasta;
use Bio::SearchIO;
use FileUpload; 
use Data::Dumper; 
our $BLASTALL   = "/home/gdroc/bin/blastall"; 
our $DEBUG      = 0;

my $param = new CGI;
my $program             = "blastn"; 
my $format              = 0;
my $thread              = 8;
my $extend_hit_region   = 150;
my $interval            = 3;
my $threshold_hit_polyA = 20;
my $threshold_tsd       = 10;
 
my $clipboard           = $param->param("clipboard");
#my $file                = $param->param("file");
my $species             = $param->param("species");
my $min_identity        = $param->param("min_identity");
my $min_hsp_length      = $param->param("hsp_length");
my $tsd_length          = $param->param("tsd_length");
my $evalue              = $param->param("evalue"); 
our $width              = $param->param("width");
our $release;

my $current_list  = "/var/www/species.json";
 
local $/;
open( my $fh, '<', $current_list );
my  $json_text   = <$fh>;
my  $perl_scalar = decode_json( $json_text );

my $common_name;
my $optgroup;
foreach my $category (keys %$perl_scalar) { 
	for my $item( @{$perl_scalar->{$category}} ){
		if ($item->{directory} eq $species) {
			$common_name = $item->{common_name};
			$release     = $item->{release};
			$optgroup    = $category; 
		} 
	}
}

my $temp_path = "/var/www/tmp/tmp$$";
mkdir($temp_path,0777);  
my $directory = $temp_path;  
my $base_path = "http://endorphine.igh.cnrs.fr/tmp/tmp$$";
my $input = $temp_path ."/temp_upload$$";
#if ($file ne "") {
#    $input = &upload_file($file,$input);
#}
#if ($clipboard ne "") {
    $input = &clipboard($clipboard,$input);
#}
system("chmod 777 $input");
our $abbrev;
our %rt_result;
our %count;
our %count_tsd; 
our $count_only_tsd;
our $repeat = "/var/www/bank/repbase/repbase.fna";
our $database = "/var/www/bank/".$species."/".$species."_".$release.".fna";
our $file_chr_length = "/var/www/bank/".$species."/".$species."_".$release.".txt";
our $species_rename = $species;
$species_rename =~ s/_/ /;
$species_rename = ucfirst($species_rename);      
 
our ($repeat_id, $repeat_max_length) = &get_length($repeat);
our ($query_id , $query_length) = &get_length($input); 

our  ( $result_dir, $output_dir , $file_txt, $file_tsd_txt, $file_fasta , $file_xls , $blast_output,$file_url_txt,$file_url_fasta,$file_url_xls,$file_url_tsd_txt,$file_svg,$file_url_svg) = &create_directory($directory,$release,$species); 
 
&run_blast( $input, $database, $blast_output );


print $param->header('application/json');

my ($json_ref,$num_hsp_ok,$num_hsp_total) = &get_result($blast_output);

my $json;
$json->{"num_result"} 		 = $num_hsp_ok;
$json->{"num_result_total"}  = $num_hsp_total;
$json->{"result_dir"} 		 = $result_dir;
$json->{"output_dir"}   	 = $output_dir;
$json->{"file_txt"}   	 	 = $file_txt;
$json->{"file_svg"}   	 	 = $file_svg;
$json->{"file_chr_length"}   = $file_chr_length;
$json->{"file_tsd_txt"}   	 = $file_tsd_txt;
$json->{"file_fasta"}   	 = $file_fasta;
$json->{"file_xls"}	         = $file_xls; 
$json->{"file_url_txt"}   	 = $file_url_txt;
$json->{"file_url_svg"}   	 = $file_url_svg;
$json->{"file_url_fasta"}    = $file_url_fasta;
$json->{"file_url_xls"}   	 = $file_url_xls;
$json->{"file_url_tsd_txt"}  = $file_url_tsd_txt;
$json->{"species"}       	 = $species;
$json->{"repeat_max_length"} = $repeat_max_length;
$json->{"query_length"}      = $query_length;
$json->{"query_name"}        = $query_id;
$json->{"tsd_length"}        = $tsd_length;
$json->{"species_rename"}    = $species_rename;
$json->{"abbrev"}    		 = $abbrev;
$json->{"database"}   		 = $database;
$json->{"repeat"}        	 = $repeat;
$json->{"width"}      		 = $width;
$json->{"result"}     		 = $json_ref;
$json->{"release"}     		 = $release;
$json->{"hsp_length"}     	 = $min_hsp_length;
$json->{"min_identity"}      = $min_identity; 

my $json_val = to_json($json );
print $json_val;

 
 
sub upload_file {
    my ($file,$file_temp) = @_;
    my $upload = new FileUpload($file); 
    if ($upload) {
        my @chararr = $upload->allow_char('\\\/'); 
        $upload->save_as("$file_temp");
    }
    return $file_temp;
}


sub clipboard {
    my ($clipboard,$file_temp) = @_;
    my @donnees = (split(/\n/,$clipboard));
    my (@sequence,$name);
    open (OUT,">$file_temp") or die "Cannot write in $file_temp!!!<br />";
    my $cpt = 0;
    my $seq = "";
    for (my $i=0;$i<=$#donnees;$i++) {
		if ($i==0 && substr($donnees[$i],0,1) ne ">") {
	    	$cpt++;
	    	$name = ">Query_".$cpt ."\n";
		}
		if (substr($donnees[$i],0,1) eq ">") {
	    	$name = $donnees[$i];
		}
		else {
	    	$seq .= $donnees[$i];
		}
		print OUT $name . $seq ."\n";          
		$name = "";
		$seq = "";
    }
    close OUT;
    return $file_temp;
}


sub create_directory {
    my ($directory,$release,$species) = @_; 
    my $result_dir =  $directory . "/result_" . $release;
    my $output_dir =  $directory . "/output_" . $release;
    my $result_url_dir =  $base_path . "/result_" . $release;
    my $output_url_dir =  $base_path . "/output_" . $release;
    
    my @dir = ($directory, $output_dir, $result_dir,$output_dir . "/wordfinder",$output_dir . "/wordmatch",$output_dir . "/blast"); 
    foreach my $dir (@dir) {
		system("mkdir $dir ;chmod -Rf 777 $dir") unless -e $dir;
    }
    my $blast_output;
    my ($file_txt,$file_fasta,$file_xls,$file_tsd_txt,$file_svg );
    my ($file_url_txt,$file_url_fasta,$file_url_xls,$file_url_tsd_txt,$file_url_svg );
    $species =~ s/\s/\_/g;
    my ($genus,$species1) = (split(/\_/,$species)); 
    $abbrev       = substr($genus,0,1)."_".substr($species1,0,3)."_".$release ."_" .$query_id;
    $file_txt     = $result_dir . "/" . $species . "_" . $release ."_" .$query_id. ".txt";
    $file_svg     = $result_dir . "/" . $species . "_" . $release ."_" .$query_id. ".svg";
    $file_tsd_txt = $result_dir . "/" . $species . "_" . $release  ."_" .$query_id. "_tsd.txt";
    $file_fasta   = $result_dir . "/" . $species . "_" . $release  ."_" .$query_id. ".fna";
    $file_xls     = $result_dir . "/" . $species . "_" . $release  ."_" .$query_id. ".xls"; 
    $blast_output = $output_dir . "/blast/" . join("_",$species, $release ,$query_id . ".out");  

    $file_url_txt     = $result_url_dir. "/" . $species . "_" . $release ."_" .$query_id. ".txt";
    $file_url_svg     = $result_url_dir. "/" . $species . "_" . $release ."_" .$query_id. ".svg";
    $file_url_tsd_txt = $result_url_dir. "/" . $species . "_" . $release  ."_" .$query_id. "_tsd.txt";
    $file_url_fasta   = $result_url_dir . "/" . $species . "_" . $release  ."_" .$query_id. ".fna";
    $file_url_xls     = $result_url_dir . "/" . $species . "_" . $release  ."_" .$query_id. ".xls"; 
    
    return ( $result_dir, $output_dir , $file_txt, $file_tsd_txt, $file_fasta , $file_xls , $blast_output,$file_url_txt,$file_url_fasta,$file_url_xls,$file_url_tsd_txt,$file_svg,$file_url_svg);
}

sub get_length {
    my $file   = shift;
    my $length = 0;
    my $name;
    my $in     = new Bio::SeqIO(
				-format => 'fasta',
				-file   => $file
				);
    while ( my $seqobj = $in->next_seq ) {
		if ( $seqobj->length > $length ) {
	    	$length = $seqobj->length();
	   		$name   = $seqobj->display_id();
		}
    }
    return ($name,$length);
} 

sub run_blast {
    my ( $query, $database, $output, $evalue ) = @_;
    $evalue = 10 unless $evalue;
    if ( $BLASTALL =~ /blastall/ ) {   
 		system("$BLASTALL -i $query -d $database -a $thread -p $program -o $output -m 8 -e $evalue -F none 2>/dev/null 1>/dev/null");
    }
    else {
		system("$BLASTALL -query $query -db $database -task $program -num_threads $thread -out $output -outfmt 6 -evalue $evalue") unless -e $output;
    }
    return $output;
}

sub get_result {
    my $blast_result = shift;
    my $parser       = new Bio::SearchIO(
					 -file   => $blast_result,
					 -format => 'blasttable'
					 ); 
    my $cpt_hsp_result = 0;
    my $cpt_hsp_ok     = 0;
    my $num_hsps; 
    my @json; 
  
    while ( my $result = $parser->next_result ) {
		foreach my $hit ( sort { $a <=> $b } $result->hits ) {
	    	foreach my $hsp ( sort { $a->start('hit') <=> $b->start('hit') }  $hit->hsps() ) {
				$cpt_hsp_result++;
				if (   $hsp->percent_identity >= $min_identity  && $hsp->hsp_length >= $min_hsp_length ) {
		    		$cpt_hsp_ok++; 
                	push @json , {
                		id          	=> $cpt_hsp_ok,
                		query_name  	=> $result->query_name,
                		hit_name		=> $hit->name,
                		hsp_start_hit   => $hsp->start('hit'),
                		hsp_end_hit     => $hsp->end('hit'),
                		hsp_start_query => $hsp->start('query'),
                		hsp_end_query   => $hsp->end('query'),
                		hsp_length      => $hsp->hsp_length(),
                		strand      	=> $hsp->strand('hit') ,
                		percent_identity=> $hsp->percent_identity,
                		gaps			=> $hsp->gaps 
                	};    
				}
	    	}
		}
    } 
    return (\@json,$cpt_hsp_ok,$cpt_hsp_result); 
}
    

1;
