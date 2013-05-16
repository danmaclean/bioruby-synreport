$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))


  
  require 'bio-synreport'
  require 'pp'
  
  #this is how you use it... no really! 
  
  db = Bio::Util::SynReport.new(:gff => ARGV[0], :fasta => ARGV[1], :verbose => true)
  chr, pos, ref,alt = 'Chr2', 15973794, 'C', 'T'
  pp db.mutation_info(chr,pos,alt) 
  
  chr, pos, ref,alt = 'Chr3',11934337, 'g', 'a'
  pp db.mutation_info(chr,pos,alt)
  
    chr, pos, ref,alt = 'Chr5',22909367, 'g', 'a'
  pp db.mutation_info(chr,pos,alt) 
  
   chr, pos, ref,alt = 'Chr5',24163458, 'g', 'a'
  pp db.mutation_info(chr,pos,alt)   
  