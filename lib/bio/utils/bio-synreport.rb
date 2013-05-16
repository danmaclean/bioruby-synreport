require 'rubygems'
require 'pp'
require 'bio'

module Bio
  class Util
    
<<<<<<< HEAD
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
=======
    class MrnaModel
      attr_accessor :seqname, :gff_id, :strand, :cds, :sequences
  
      def initialize(chr, id, strand, cds_arr, seq_arr)
        @seqname, @gff_id, @strand, @cds, @sequences = chr, id, strand, cds_arr, seq_arr
      end
  
      def includes?(seq, point)
        @cds.each {|start, stop| return true if @seqname == seq and point.to_i >= start and point.to_i <= stop}
        false
      end
  
      def seq
        @sequences.join
      end
  
      def substitution_info(chr,point,alt)
        cds_start = @cds.first.first
        running_total = 0
        @cds.each do |start,stop|
          if point.to_i >= start and point.to_i <= stop
            offset = case @strand
            when "+"
              #offset = 
              (point.to_i - start) + running_total
            when "-"
              (stop - point.to_i) + running_total 
            end  #offset = how far into cds SNP is
            codon_number = offset / 3
            position_in_codon = offset % 3
            #pp [offset, codon_number, position_in_codon] 
            codon_array = []; Bio::Sequence::NA.new(self.seq).window_search(3,3) {|b| codon_array << b}
            codon = codon_array[codon_number]
            nt = codon[position_in_codon]
            new_codon = codon.dup
            new_codon[position_in_codon] = alt.downcase
            #pp [codon, position_in_codon, nt, new_codon]
            a = Bio::Sequence::NA.new(codon).translate.codes.first
            b =  Bio::Sequence::NA.new(new_codon).translate.codes.first
            sub_type = a == b ? "SYN" : "NON_SYN"
            return {:id => @gff_id, 
                    :chr => @seqname, 
                    :strand => @strand, 
>>>>>>> 188a1a611ad6334046551c7bba186dc1c7ae85af
                    :position => point,
                    :original_codon => codon, 
                    :original_residue => a || 'stop', 
                    :mutant_codon => new_codon, 
                    :mutant_residue =>b || 'stop', 
<<<<<<< HEAD
                    :position_in_codon => position + 1, 
                    :substitution_type => sub_type
                    }
=======
                    :position_in_codon => position_in_codon + 1, 
                    :substitution_type => sub_type
                    }
          end
          running_total += (stop - start)
          running_total += 1 if @strand == '-' #how far we are into the cds
        end
>>>>>>> 188a1a611ad6334046551c7bba186dc1c7ae85af
      end
      
    end#class end
    
    
    class SynReport
      #attr_accessor :cdshash, :cds_list, :mRNAhash, :seqhash
      
      def initialize(opts)
<<<<<<< HEAD
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
=======
        @gene_array = []
        @cdshash = Hash.new {|h,k| h[k] = Hash.new {|a,b| a[b] = [] } }
        @mRNAhash = Hash.new {|h,k| h[k] = Hash.new {|a,b| a[b] = [] } }
        File.open(opts[:gff], "r").each do |gffline|
          record=Bio::GFF::GFF3::Record.new(gffline)
          if(record.feature_type == 'gene')
              @gene_array << [record.seqname, record.id]
          elsif(record.feature_type == 'CDS' or record.feature_type == 'mRNA')
            parents = record.get_attributes('Parent')
            parents.each do |parent|  
              if record.feature_type == 'CDS'
                @cdshash[record.seqname][parent] << record
              else
                @mRNAhash[record.seqname][parent] << record
              end
            end
          end
        end
        $stderr.puts "Loaded GFF..." if opts[:verbose]
        @seqhash = {}
        Bio::FastaFormat.open(opts[:fasta]).each { |seq| @seqhash[seq.entry_id] = seq.to_seq }
        $stderr.puts "Loaded Seq..." if opts[:verbose]
        
        @models = Hash.new {|h,k| h[k] =  []  }
        $stderr.puts "Building models..." if opts[:verbose] 
        @gene_array.each do |gene|

          mRNAs=@mRNAhash[gene.first][gene.last]
          mRNAs.each do |mRNA|
            next if @seqhash[gene.first].nil?
            cdsa = []
            seqs = []
            cdsary=@cdshash[gene.first][mRNA.id]
            cdsary.each {|c| cdsa << [c.start, c.end]} 
            cdsa.sort!
            cdsa.reverse! if mRNA.strand == '-'
            
            cdsa.each do |cds|

              #cdsa << [cds.start, cds.end]
              if mRNA.strand == '+'
                seqs << Bio::Sequence::NA.new(@seqhash[mRNA.seqname].splicing("#{cds.first}..#{cds.last}") )
              elsif mRNA.strand == "-"
                seqs << Bio::Sequence::NA.new(@seqhash[mRNA.seqname].splicing("#{cds.first}..#{cds.last}") ).complement
              end
            end
            @models[mRNA.seqname] << Bio::Util::MrnaModel.new(mRNA.seqname, mRNA.id, mRNA.strand, cdsa, seqs )
            #pp @models[mRNA.seqname][-1].cds if mRNA.id == 'AT2G17530.1' or mRNA.id == 'AT2G17550.1'
          end
        end
        $stderr.puts "Models built..." if opts[:verbose]
      end#init end
        
        def is_in_cds?(chr,point)
          @self.mutation_info(chr,point) ? true : false
>>>>>>> 188a1a611ad6334046551c7bba186dc1c7ae85af
        end
        
        #returns mutation info if point in CDS, if not in CDS returns false
        def mutation_info(chr,pos,alt)
<<<<<<< HEAD
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
=======

          @models[chr].each do |m|
             if m.includes?(chr,pos)
               return m.substitution_info(chr,pos,alt)   
             end
          end
          false
>>>>>>> 188a1a611ad6334046551c7bba186dc1c7ae85af
        end
        
        
    end#class end
  end#class util end
end# module end