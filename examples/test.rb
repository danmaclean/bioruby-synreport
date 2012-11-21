$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))


  
  require 'bio-synreport'
  require 'pp'
  
  #this is how you use it...
  
  db = Bio::Util::SynReport.new(:gff => ARGV[0], :fasta => ARGV[1], :verbose => true)
  chr, pos, ref,alt = 'Chr2',7634495, 'a', 't'
  pp db.mutation_info(chr,pos,alt) 
  
  chr, pos, ref,alt = 'Chr3',123456, 'a', 't'
  pp db.mutation_info(chr,pos,alt)
  
    chr, pos, ref,alt = 'Chr2',7626518, 'a', 't'
  pp db.mutation_info(chr,pos,alt) 