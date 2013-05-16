require 'rubygems'
require 'pp'
require 'bio'

module Bio
  class Util
    
    class MrnaModel < Bio::GFF::GFF3::Record
      attr_accessor :seq, :cds
      def initialize(gff_line)
        super gff_line
        @cds = []
      end
      
      def includes?(seq, point)
        return true if self.seqname == seq and point.to_i >= self.cds_start and point.to_i <= self.cds_end
        false
      end
      
      def cds_start
        @cds.flatten.min
      end
  
      def cds_end
        @cds.flatten.max
      end
      
      def get_nt_number_in_cds(point)
        to_count = @cds.sort.select {|a| a.first <= point}
        in_block = to_count.pop
        distance_in = to_count.inject(1) {|tot, b| tot + ( b.last - b.first) + 1 }
        overhang =  point - in_block.first 
        left_section = distance_in + overhang
        
        if self.strand == '-'
            length = @cds.sort.inject(0) {|tot, b| tot + ( b.last - b.first) + 1 }
            return length - left_section + 1
        end
        
        return left_section
      end
      
      def codon_index(dist)
        (dist - 1) / 3
      end
      
      def codon_position(dist)
        (dist - 1) % 3
      end
      
      def codon_array
        codon_array = []; Bio::Sequence::NA.new(self.seq).window_search(3,3) {|b| codon_array << b}
        codon_array
      end
      
      def nt
      end
      
      ##returns codon and position of nucleotide 
      def codon_and_index(point)
         distance_into_cds = get_nt_number_in_cds point
         codon_idx = codon_index distance_into_cds
         codon_list = codon_array
         codon = codon_list[codon_idx]
         pos = codon_position(distance_into_cds)
         [codon,pos]
      end
      
      def substitution_info(point,alt)
            codon, position = codon_and_index(point)
            new_codon = codon.dup
            new_codon[position] = alt.downcase
            
            a = Bio::Sequence::NA.new(codon).translate.codes.first
            b =  Bio::Sequence::NA.new(new_codon).translate.codes.first
            sub_type = a == b ? "SYN" : "NON_SYN"
            return {#:id => self.gffid, 
                    :chr => self.seqname, 
                    :strand => self.strand, 
                    :position => point,
                    :original_codon => codon, 
                    :original_residue => a || 'stop', 
                    :mutant_codon => new_codon, 
                    :mutant_residue =>b || 'stop', 
                    :position_in_codon => position + 1, 
                    :substitution_type => sub_type
                    }

      end
      
    end#class end
    
    
    class SynReport
      #attr_accessor :cdshash, :cds_list, :mRNAhash, :seqhash
      
      def initialize(opts)
        cdses = []
        mrna_list = []
        seqs = Hash.new
        
        Bio::FastaFormat.open(opts[:fasta]).each { |seq| seqs[seq.entry_id] = seq.to_seq }
        $stderr.puts "Loaded Seq..." if opts[:verbose]
        
        
        @mrnas = Hash.new  {|h,k| h[k] = Hash.new}
        File.open(opts[:gff], "r").each do |gffline|
          record = Bio::GFF::GFF3::Record.new(gffline)
          if record.feature_type == 'mRNA'
            mrna_list << Bio::Util::MrnaModel.new(gffline)
          elsif record.feature_type =='CDS'
            cdses << record
          end
        end
        
        mrna_list.each do |mrna|
          mrna_id = mrna.get_attributes("ID")
          $stderr.puts "No ID for #{cds}" if mrna_id.empty?
          mrna_id = mrna_id.first
          @mrnas[mrna.seqname][mrna_id] = mrna
          @mrnas[mrna.seqname][mrna_id].seq = seqs[mrna_id].seq
        end
        
        cdses.each do |cds|
          cds_parent = cds.get_attributes("Parent")
          $stderr.puts "No Parent for #{cds}" if cds_parent.empty?
          cds_parent = cds_parent.first
          @mrnas[cds.seqname][cds_parent].cds << [cds.start,cds.end]
        end
        $stderr.puts "Loaded GFF..." if opts[:verbose]


      end#init end
        
        def is_in_cds?(chr,point)
          self.mutation_info(chr,point,"a") ? true : false
        end
        
        #returns mutation info if point in CDS, if not in CDS returns false
        def mutation_info(chr,pos,alt)
          pos = pos.to_i
          #cant do indels ...
          return nil if alt.length > 1
          begin
            @mrnas[chr].each_pair do |mrna_id, mrna|
              if mrna.includes?(chr,pos)
                return mrna.substitution_info(pos,alt)   
              end
           end
            false
          rescue
            #somthing unpredicatable went wrong and we couldnt do the conversion ...
            return nil
          end
        end
        
        
    end#class end
  end#class util end
end# module end