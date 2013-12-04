$(document).ready(function()
 {
 	var selectoptions_osj = {
    		"blastn": {
    	         "key" : 'blastn',
                 "defaultvalue" : 'TIGR6_pseudochromosome',
    	         "values" : {
                     "Gene sequences. Includes UTRs and introns.": 'TIGR6_seq_20090130',
                     "CDS sequences. Does not include intron sequences or UTRs.": 'TIGR6_cds_20090130',
                     "MSU Pseudochromosomes": 'TIGR6_pseudochromosome'
                     }
              	},
    	    	"blastp": {
                 "key" : 'blastp',
                 "defaultvalue" : 'TIGR6_pep_20090130',
                 "values" : {
                     "Translated protein sequences": 'TIGR6_pep_20090130' 
                     }
              	},
            	"blastx": {
                 "key" : 'blastx',
                 "defaultvalue" : 'TIGR6_pep_20090130',
                 "values" : {
                     "Translated protein sequences": 'TIGR6_pep_20090130' 
                     }
              	},
    		"tblastn": {
    	         "key" : 'tblastn',
                 "defaultvalue" : 'TIGR6_pseudochromosome',
    	         "values" : {
                     "Gene sequences. Includes UTRs and introns.": 'TIGR6_seq_20090130',
                     "CDS sequences. Does not include intron sequences or UTRs.": 'TIGR6_cds_20090130',
                     "MSU Pseudochromosomes": 'TIGR6_pseudochromosome'
                     }
              	},
    		"tblastx": {
    	         "key" : 'tblastx',
                 "defaultvalue" : 'TIGR6_pseudochromosome',
    	         "values" : {
                     "Gene sequences. Includes UTRs and introns.": 'TIGR6_seq_20090130',
                     "CDS sequences. Does not include intron sequences or UTRs.": 'TIGR6_cds_20090130',
                     "MSU Pseudochromosomes": 'TIGR6_pseudochromosome'
                     }
              	}
    	};
	var selectoptions_osi = {
    		"blastn": {
    	         "key" : 'blastn',
                 "defaultvalue" : 'BGI2_pseudochromosome',
    	         "values" : {
                     "CDS sequences": '9311_glean_gene_cds.fa',
                     "BGI Pseudochromosomes.": 'BGI2_pseudochromosome'
                     }
              	},
            	"blastp": {
                 "key" : 'blastp',
                 "defaultvalue" : '9311_glean_gene_pep.fa',
                 "values" : {
                     "Translated protein sequences": '9311_glean_gene_pep.fa' 
                     }
              	},
            	"blastx": {
                 "key" : 'blastx',
                 "defaultvalue" : '9311_glean_gene_pep.fa',
                 "values" : {
                     "Translated protein sequences": '9311_glean_gene_pep.fa' 
                     }
              	},
              	"tblastn": {
    	         "key" : 'tblastn',
                 "defaultvalue" : 'BGI2_pseudochromosome',
    	         "values" : {
                     "CDS sequences": '9311_glean_gene_cds.fa',
                     "BGI Pseudochromosomes.": 'BGI2_pseudochromosome'
                     }
        	},
    		"tblastx": {
    	         "key" : 'tblastx',
                 "defaultvalue" : 'BGI2_pseudochromosome',
    	         "values" : {
                     "CDS sequences": '9311_glean_gene_cds.fa',
                     "BGI Pseudochromosomes.": 'BGI2_pseudochromosome'
                     }
             	}
    	}; 
   	var selectoptions_at = {
    		"blastn": {
    	         "key" : 'blastn',
                 "defaultvalue" : 'TAIR9_pseudochromosome',
    	         "values" : {
                     "Gene sequences. Includes UTRs and introns.": 'TAIR10_seq_20110103_representative_gene_model',
                     "CDS sequences. Does not include intron sequences or UTRs.": 'TAIR10_cds_20110103_representative_gene_model',
                     "TAIR Pseudochromosomes.": 'TAIR9_pseudochromosome'
                     }
              	},
            	"blastp": {
                 "key" : 'blastp',
                 "defaultvalue" : 'TAIR10_pep_20110103_representative_gene_model',
                 "values" : {
                     "Translated protein sequences": 'TAIR10_pep_20110103_representative_gene_model' 
                     }
            	},
            	"blastx": {
                 "key" : 'blastx',
                 "defaultvalue" : 'TAIR10_pep_20110103_representative_gene_model',
                 "values" : {
                     "Translated protein sequences": 'TAIR10_pep_20110103_representative_gene_model' 
                     }
              	},
              	"tblastn": {
    	         "key" : 'tblastn',
                 "defaultvalue" : 'TAIR9_pseudochromosome',
    	         "values" : {
                     "Gene sequences. Includes UTRs and introns.": 'TAIR10_seq_20110103_representative_gene_model',
                     "CDS sequences. Does not include intron sequences or UTRs.": 'TAIR10_cds_20110103_representative_gene_model',
                     "TAIR Pseudochromosomes.": 'TAIR9_pseudochromosome'
                     }
              	},
    		"tblastx": {
    	         "key" : 'tblastx',
                 "defaultvalue" : 'TAIR9_pseudochromosome',
    	         "values" : {
                     "Gene sequences. Includes UTRs and introns.": 'TAIR10_seq_20110103_representative_gene_model',
                     "CDS sequences. Does not include intron sequences or UTRs.": 'TAIR10_cds_20110103_representative_gene_model',
                     "TAIR Pseudochromosomes.": 'TAIR9_pseudochromosome'
                     }
              	}
    	};
    	var selectoptions_sb = {
    		"blastn": {
    	         "key" : 'blastn',
                 "defaultvalue" : 'JGI1_pseudochromosome',
    	         "values" : {
                     "Genes sequences": 'Sbicolor_79_transcript.fa',
                     "CDS sequences": 'Sbicolor_79_cds.fa',
                     "JGI Pseudochromosomes.": 'JGI1_pseudochromosome'
                     }
              	},
            	"blastp": {
                 "key" : 'blastp',
                 "defaultvalue" : 'Sbicolor_79_pep.fa',
                 "values" : {
                     "Translated protein sequences": 'Sbicolor_79_pep.fa' 
                     }
              	},
            	"blastx": {
                 "key" : 'blastx',
                 "defaultvalue" : 'Sbicolor_79_pep.fa',
                 "values" : {
                     "Translated protein sequences": 'Sbicolor_79_pep.fa' 
                     }
              	},
              	"tblastn": {
    	         "key" : 'tblastn',
                 "defaultvalue" : 'JGI1_pseudochromosome',
    	         "values" : {
                     "Genes sequences": 'Sbicolor_79_transcript.fa',
                     "CDS sequences": 'Sbicolor_79_cds.fa',
                     "JGI Pseudochromosomes.": 'JGI1_pseudochromosome'
                     }
              	},
    		"tblastx": {
    	         "key" : 'tblastx',
                 "defaultvalue" : 'JGI1_pseudochromosome',
    	         "values" : {
                     "Genes sequences": 'Sbicolor_79_transcript.fa',
                     "CDS sequences": 'Sbicolor_79_cds.fa',
                     "JGI Pseudochromosomes.": 'JGI1_pseudochromosome'
                     }
              	}
    	};
	$('#program_osj').doubleSelect('db_osj', selectoptions_osj);
	$('#program_osi').doubleSelect('db_osi', selectoptions_osi);
	$('#program_at').doubleSelect('db_at', selectoptions_at);
	$('#program_sb').doubleSelect('db_sb', selectoptions_sb);
 }); 
