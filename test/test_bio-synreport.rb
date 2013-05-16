require 'helper'


class TestMrnaModel < Test::Unit::TestCase
  
  def setup
    @AT2G17530_1 = Bio::Util::MrnaModel.new('Chr2	TAIR10	mRNA	7626359	7629149	.	+	.	ID=AT2G17530.1;Parent=AT2G17530;Name=AT2G17530.1;Index=1')
    @AT2G17530_1.cds = [ [7626518,7626700],[7627310,7628101],[7628198,7628245], [7628325,7628624] ]
    @AT2G17530_1.seq = "ATGTCGTGTTCATCCTCATCTGGATCAGAGGAAGACGATGAGGGTTTCGATGCTTACCGT
AAAGGTGGATATCACGCCGTTAGAATCGGAGACCAGTTTGCCGGTGGCCGTTACATTGCT
CAGAGAAAGCTTGGTTGGGGCCAATTCTCCACCGTTTGGCTTGCCTATGATACTCGCACT
TCTAATTATGTTGCTTTGAAGATTCAGAAGAGCGCCTTACAATTTGCTCAAGCTGCACTT
CATGAAATCGAACTTCTTCAAGCTGCTGCTGATGGGGATCCTGAAAATACCAAGTGTGTT
ATTCGTCTTATTGATGACTTCAAGCACGCAGGTCCCAACGGGCAGCATTTATGCATGGTG
CTCGAGTTTCTTGGCGATAGCTTGCTGCGTTTGATTAAATATAACCGTTATAAAGGGATG
GAGTTAAGTAAAGTGCGGGAGATATGCAAATGTATACTGACTGGTCTAGATTATTTGCAC
CGTGAACTCGGTATGATTCACTCCGACTTAAAACCCGAAAACATTCTTCTTTGTTCCACC
ATTGACCCTGCCAAGGATCCTATCAGATCCGGACTAACACCGATACTAGAAAAGCCCGAG
GGGAACCAAAACGGTACATCAACAATGAATCTGATTGAGAAGAAGTTGAAGAGGAGAGCA
AAAAAAGCGGCTGCTAAAATATCAGGAAGAAGAGTTTCGATAGTAGGTTTAAGTGAAACA
CCGAAAAAGAACAAGAGAAACTTGGATGGGATTGATATGAGATGCAAAGTTGTCGACTTC
GGGAACGGGTGTTGGGCTGATAACAAATTTGCAGAAGAAATACAAACAAGACAGTACAGA
GCTCCTGAAGTAATACTTCAGTCAGGTTACTCTTACTCTGTTGATATGTGGTCTTTCGCT
TGTACTGCTTTTGAGCTTGCTACAGGCGATATGCTTTTCGCTCCAAAAGAGGGAAATGGT
TACGGAGAAGACGAGGACCACCTTGCTCTTATGATGGAACTCTTAGGAAAAATGCCTCGA
AAGATTGCCATTGGAGGTGCGAGATCAAAGGATTACTTTGACAGACACGGCGACTTGAAG
AGGATCCGGAGATTAAAATACTGGCCACTCGACCGTTTACTGATTGATAAATACAAGCTT
CCAGAAGCAGAAGCACGAGAATTTGCGGATTTTCTCTGCCCGATAATGGATTTTGCACCT
GAGAAACGACCAACTGCACAACAATGTCTGCAACATCCATGGTTGAATCTAAGGACACAG
AACAATGAAGATGATATAGAAGGTCAGATGAGTAACATGCAGATCAAAGGTTCATGTTCT
TGA".gsub(/\n/,"")
    
    @AT2G38130_1 = Bio::Util::MrnaModel.new('Chr2	TAIR10	mRNA	15978512	15980749	.	-	.	ID=AT2G38130.1;Parent=AT2G38130;Name=AT2G38130.1;Index=1')
    @AT2G38130_1.cds = [[15978639, 15978854],[15979551, 15979606],[15979746, 15979866],[15979966, 15980145]]
    @AT2G38130_1.seq = "ATGGAGAAAGAGATGGAAGATAAAGAAGAATTCGATGAGGGTGAGATTGAGTACACGAGT
TATGCTGGTGAGCATCATCTGCCATTGATTATGTCTCTTGTTGACCAAGAACTTAGTGAA
CCTTACTCCATCTTTACTTACCGGTACTTCGTCTACCTCTGGCCGCAGCTATGCTTCCTG
GCCTTTCACAAAGGTAAATGCGTAGGAACCATAGTCTGTAAGATGGGGGATCATCGACAG
ACTTTCAGAGGGTACATCGCTATGTTGGTTGTGATTAAACCATATCGTGGCCGAGGCATA
GCCTCAGAGCTTGTCACAAGAGCGATAAAAGCGATGATGGAATCAGGCTGTGAAGAGGTA
ACTCTGGAGGCAGAAGTGAGTAACAAAGGAGCATTAGCACTATATGGGCGACTCGGGTTT
ATAAGAGCCAAACGGCTATACCACTATTACTTGAATGGGATGGATGCTTTTCGCCTGAAG
CTCTTGTTCCCTAAGCCTCGTGTACCTCAAATACCTTCTCAAGTTCAAACCCAACAAGAG
TATGAGACCTTTCCTAGGCCTCGTGTACCTTAA".gsub(/\n/,"")
  end
  
  def test_points_in_cds
    [7626518,7626528,7627320,7628208,7628335,7628624].each do |point|
        assert @AT2G17530_1.includes?('Chr2', point), "#{point} should be reported in cds"  
    end
  end
  
  def test_gets_end
    assert_equal 7628624, @AT2G17530_1.cds_end, "cds has wrong end"
  end
  
  def test_gets_start
    assert_equal 7626518, @AT2G17530_1.cds_start, "cds has wrong start"
  end
  
  def test_get_nt_number_in_cds
    #at start of first cds segment, distance into cds is 1 
    assert_equal 1, @AT2G17530_1.get_nt_number_in_cds(7626518), "Offset is wrong"
    #at end of first cds segment distance is length of cds segment (183) 
    assert_equal 183, @AT2G17530_1.get_nt_number_in_cds(7626700), "Offset is wrong"
    #at start of second cds segment distances is length of first + distance into second (1) = 184
    assert_equal 184,  @AT2G17530_1.get_nt_number_in_cds(7627310), "Offset is wrong"
    #ten into second cds segment distance is length of first (184) + distance into second (10) = 194
    assert_equal 194,  @AT2G17530_1.get_nt_number_in_cds(7627320), "Offset is wrong"
    #last position is length of all cds segments
    assert_equal 1323, @AT2G17530_1.get_nt_number_in_cds(7628624), "offset is offset"
    
    #now negative strand gene
    #at end of last cds segment, distance into cds is 1 
    assert_equal 1, @AT2G38130_1.get_nt_number_in_cds(15980145), "Offset is wrong"
    #at start of last cds segment distance is length of cds segment (180) 
    assert_equal 180, @AT2G38130_1.get_nt_number_in_cds(15979966), "Offset is wrong"
    #at end of second cds segment distances is length of first + distance into second (1) = 181
    assert_equal 181,  @AT2G38130_1.get_nt_number_in_cds(15979866), "Offset is wrong"
    #ten from end of second cds segment distance is length of first (181) + distance into second (10) = 191
    assert_equal 191,  @AT2G38130_1.get_nt_number_in_cds(15979856), "Offset is wrong"
    #last position is length of all cds segments 
    assert_equal 573, @AT2G38130_1.get_nt_number_in_cds(15978639)
  end
  
  def test_substitution_info
    
    ##first residue, + strand
    result = @AT2G17530_1.substitution_info(7626518, 'a')
    pp result
    assert_equal('atg', result[:original_codon])
    assert_equal('Met', result[:original_residue])
    assert_equal('atg',result[:mutant_codon])
    assert_equal('Met', result[:mutant_residue])
    assert_equal(1, result[:position_in_codon])
    
    result = @AT2G17530_1.substitution_info(7626519, 'a')
    pp result
    assert_equal('atg', result[:original_codon])
    assert_equal('Met', result[:original_residue])
    assert_equal('aag',result[:mutant_codon])
    assert_equal('Lys', result[:mutant_residue])
    assert_equal(2, result[:position_in_codon])
    
    result = @AT2G17530_1.substitution_info(7626520, 'a')
    pp result
    assert_equal('atg', result[:original_codon])
    assert_equal('Met', result[:original_residue])
    assert_equal('ata',result[:mutant_codon])
    assert_equal('Ile', result[:mutant_residue])
    assert_equal(3, result[:position_in_codon])
    
    ##first residue, - strand
    result = @AT2G38130_1.substitution_info(15980145, 'a')
    pp result
    assert_equal('atg', result[:original_codon])
    assert_equal('Met', result[:original_residue])
    assert_equal('atg',result[:mutant_codon])
    assert_equal('Met', result[:mutant_residue])
    assert_equal(1, result[:position_in_codon])

    result = @AT2G38130_1.substitution_info(15980144, 'a')
    pp result
    assert_equal('atg', result[:original_codon])
    assert_equal('Met', result[:original_residue])
    assert_equal('aag',result[:mutant_codon])
    assert_equal('Lys', result[:mutant_residue])
    assert_equal(2, result[:position_in_codon])
    
    result = @AT2G38130_1.substitution_info(15980143, 'a')
    pp result
    assert_equal('atg', result[:original_codon])
    assert_equal('Met', result[:original_residue])
    assert_equal('ata',result[:mutant_codon])
    assert_equal('Ile', result[:mutant_residue])
    assert_equal(3, result[:position_in_codon])    
    
    ##third residue second cds segment, + strand -> start pos = 7627317
    result = @AT2G17530_1.substitution_info(7627317, 'a')
    pp result
    assert_equal('gtt', result[:original_codon])
    assert_equal('Val', result[:original_residue])
    assert_equal('gat',result[:mutant_codon])
    assert_equal('Asp', result[:mutant_residue])
    assert_equal(2, result[:position_in_codon])

    ##third residue second cds segment, - strand -> start pos = 15979753
    result = @AT2G38130_1.substitution_info(15979753, 'a')
    pp result
    assert_equal('cga', result[:original_codon])
    assert_equal('Arg', result[:original_residue])
    assert_equal('cga',result[:mutant_codon])
    assert_equal('Arg', result[:mutant_residue])
    assert_equal(3, result[:position_in_codon])
    
    ##last residue + strand
    result = @AT2G17530_1.substitution_info(7628624, 'a')
    pp result
    assert_equal('tga',result[:original_codon])
    assert_equal('stop', result[:original_residue])
    assert_equal('tga',result[:mutant_codon])
    assert_equal('stop', result[:mutant_residue])  
    assert_equal(3, result[:position_in_codon]) 
    
    result = @AT2G17530_1.substitution_info(7628622, 'a')
    pp result
    assert_equal('tga',result[:original_codon])
    assert_equal('stop', result[:original_residue])
    assert_equal('aga',result[:mutant_codon])
    assert_equal('Arg', result[:mutant_residue])  
    assert_equal(1, result[:position_in_codon]) 
    
    ##last residue - strand 
    result = @AT2G38130_1.substitution_info(15978639, 'a')
    pp result 
    assert_equal('taa',result[:original_codon])
    assert_equal('stop', result[:original_residue])
    assert_equal('taa',result[:mutant_codon])
    assert_equal('stop', result[:mutant_residue])  
    assert_equal(3, result[:position_in_codon])
    
    result = @AT2G38130_1.substitution_info(15978641, 'a')
    pp result 
    assert_equal('taa',result[:original_codon])
    assert_equal('stop', result[:original_residue])
    assert_equal('aaa',result[:mutant_codon])
    assert_equal('Lys', result[:mutant_residue])  
    assert_equal(1, result[:position_in_codon])  
       
  end
  
end
=begin
                    :original_codon => codon, 
                    :original_residue => a || 'stop', 
                    :mutant_codon => new_codon, 
                    :mutant_residue =>b || 'stop', 
                    :position_in_codon => position_in_codon + 1, 
                    :substitution_type => sub_type


class TestBioSynreport < Test::Unit::TestCase

  def setup
    @plus_strand_first_cds = ''
    @plus_strand_second_cds = ''
    @minus_strand_first_cds = ''
    @minus_strand_second_cds = ''
  end

  def test_nothing
    assert true
  end

end
=end
