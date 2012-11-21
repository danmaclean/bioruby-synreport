require 'rubygems'
require 'pp'
require 'bio'

module Bio
  class Util
    
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
                    :position => point,
                    :original_codon => codon, 
                    :original_residue => a || 'stop', 
                    :mutant_codon => new_codon, 
                    :mutant_residue =>b || 'stop', 
                    :position_in_codon => position_in_codon + 1, 
                    :substitution_type => sub_type
                    }
          end
          running_total += (stop - start)
          running_total += 1 if @strand == '-' #how far we are into the cds
        end
      end
      
    end#class end
    
    
    class SynReport
      #attr_accessor :cdshash, :cds_list, :mRNAhash, :seqhash
      
      def initialize(opts)
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
        end
        
        #returns mutation info if point in CDS, if not in CDS returns false
        def mutation_info(chr,pos,alt)

          @models[chr].each do |m|
             if m.includes?(chr,pos)
               return m.substitution_info(chr,pos,alt)   
             end
          end
          false
        end
        
        
    end#class end
  end#class util end
end# module end