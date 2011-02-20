#!/usr/bin/env ruby

# A script to take a FASTA file, remove sequences that will fail, and automatically submit it to the cNLS server at http://nls-mapper.iab.keio.ac.jp/cgi-bin/NLS_Mapper_form.cgi
# Unfortunately, the fasta upload seems to fail. 
# and format it so that it can be uploaded to the cNLS mapper (classical(?) nuclear localisation signal mapper).
# The fasta output file can be uploaded to
# http://nls-mapper.iab.keio.ac.jp/cgi-bin/NLS_Mapper_form.cgi

require 'bio'

module Bio
  class CNLS
    class Result
      attr_accessor :signals
      
      def initialize
        @signals = []
      end
      
      class NLS
        attr_accessor :position, :sequence, :score
        
        # sort by score descending
        def <=>(another)
          -(@score<=>another.score)
        end
      end
      class MonopartiteNLS<NLS; end
      class BipartiteNLS<NLS; end
      
      # Is this result a positive prediction or negative prediction?
      def predicted?
        !signals.nil? and !signals.empty? 
      end
    end
    
    # A class used to automatically submit results to the cNLS webserver and parse the HTML results. 
    class Screenscraper
      require 'uri'
      require 'net/http'
      
      ACCEPTABLE_CUTOFFS = %w(2.0 3.0 4.0 5.0 6.0)
      
      # Contact the cNLS prediction server and submit the amino acid sequence for prediction. Return a Bio::CNLS::Result object. Pause after each round for pause milliseconds, so as not to overload the server.
      def self.submit(amino_acid_sequence, cut_off='3.0', seconds_pause=1)
        # contact webserver and sleep
        html = get_raw_html_result(amino_acid_sequence, cut_off, seconds_pause)
        
        # Return the parsed HTML as a CNLS::Result object
        return parse_html_result(html)
      end
      
      def self.get_raw_html_result(amino_acid_sequence, cut_off='3.0', seconds_pause=1)
        unless ACCEPTABLE_CUTOFFS.include?(cut_off)
          raise Exception, "Specified cutoff `#{cut_off}' for the cNLS screenscraper is invalid. Valid cutoffs are #{ACCEPTABLE_CUTOFFS.join(', ')}. They are strings, not floating point values."
        end
        
        # retrieve the webpage
        res = Net::HTTP.post_form(URI.parse('http://nls-mapper.iab.keio.ac.jp/cgi-bin/NLS_Mapper_y.cgi'),
        {'cut_off' => cut_off, 'typedseq' => amino_acid_sequence})
        
        # if there is an error, raise it
        unless res.kind_of?(Net::HTTPOK)
          raise Exception, "Failed to retrieve cNLS, internet connectivity problem? Using cutoff/sequence #{cutoff}/#{amino_acid_sequence}"
        end
        
        # pause the specified number of seconds
        sleep seconds_pause
        
        return res.body
      end
      
      # Given HTML corresponding to a result, return a parse object that is more programmatically palatable.
      def self.parse_html_result(html)
        result = Result.new
        
        # The mono and bi-partite regular expressions are equivalent except for the Predicted X NLS bit at the beginning, thanksfully. However, they sometimes appear to be slightly different, which is rather odd.
        monopartite_regex = /Predicted monopartite NLS<\/th>\s+<\/TR>\s*<TR bgcolor="#d0d0d0">\s*<th>Pos.<\/th>\s*<th>Sequence<\/th>\s*<th>Score<\/th>\s*<\/TR>\s*<TR><td><strong><big><code>(.*?)<\/code><\/big><\/strong><br.{0,2}><strong><big><code.{2,8}><\/big><\/strong><\/td><td><strong><big><code>(.*?)<\/code><\/big><\/strong><br.{0,2}><strong><big><code.{2,8}><\/big><\/strong><\/td><td align="center"><strong><big><code>(.*?)<\/code><\/big><\/strong><br.{0,2}><strong><big><code.{2,8}><\/big><\/strong><\/td><\/TR>/i
        bipartite_regex =     /Predicted bipartite NLS<\/th>\s+<\/TR>\s*<TR bgcolor="#d0d0d0">\s*<th>Pos.<\/th>\s*<th>Sequence<\/th>\s*<th>Score<\/th>\s*<\/TR>\s*<TR><td><strong><big><code>(.*?)<\/code><\/big><\/strong><br.{0,2}><strong><big><code.{2,8}><\/big><\/strong><\/td><td><strong><big><code>(.*?)<\/code><\/big><\/strong><br.{0,2}><strong><big><code.{2,8}><\/big><\/strong><\/td><td align="center"><strong><big><code>(.*?)<\/code><\/big><\/strong><br.{0,2}><strong><big><code.{2,8}><\/big><\/strong><\/td><\/TR>/i
        
        monopartite_no_hits = /Predicted monopartite NLS<\/th>\s*<\/tr>\s*<tr bgcolor="#d0d0d0">\s*<th>Pos.<\/th>\s*<th>Sequence<\/th>\s*<th>Score<\/th>\s*<\/tr>\s*<tr><td><strong><big><code><\/code><\/big><\/strong><\/td><td><strong><big><code><\/code><\/big><\/strong><\/td><td align="center"><strong><big><code><\/code><\/big><\/strong><\/td><\/tr>/i
        bipartite_no_hits =     /Predicted bipartite NLS<\/th>\s*<\/tr>\s*<tr bgcolor="#d0d0d0">\s*<th>Pos.<\/th>\s*<th>Sequence<\/th>\s*<th>Score<\/th>\s*<\/tr>\s*<tr><td><strong><big><code><\/code><\/big><\/strong><\/td><td><strong><big><code><\/code><\/big><\/strong><\/td><td align="center"><strong><big><code><\/code><\/big><\/strong><\/td><\/tr>/i
        monopartite_no_hits2 = /Predicted monopartite NLS<\/th>\s*<\/TR>\s*<TR bgcolor="#d0d0d0">\s*<th>Pos.<\/th>\s*<th>Sequence<\/th>\s*<th>Score<\/th>\s*<\/TR>\s*<TR><td><strong><big><code \/><\/big><\/strong><\/td><td><strong><big><code \/><\/big><\/strong><\/td><td align="center"><strong><big><code \/><\/big><\/strong><\/td><\/TR>/i
        bipartite_no_hits2 =     /Predicted bipartite NLS<\/th>\s*<\/TR>\s*<TR bgcolor="#d0d0d0">\s*<th>Pos.<\/th>\s*<th>Sequence<\/th>\s*<th>Score<\/th>\s*<\/TR>\s*<TR><td><strong><big><code \/><\/big><\/strong><\/td><td><strong><big><code \/><\/big><\/strong><\/td><td align="center"><strong><big><code \/><\/big><\/strong><\/td><\/TR>/i
        
        split_regex = /<\/code><\/big><\/strong><br.{0,2}><strong><big><code>/
        
        
        # parse out monopartite signals
        if matches = html.match(monopartite_regex)
          positions = matches[1].split(split_regex)
          seqs = matches[2].split(split_regex)
          scores = matches[3].split(split_regex)
          
          positions.each_with_index do |pos, i|
            nls = Result::MonopartiteNLS.new
            nls.position = pos.to_i
            nls.sequence = seqs[i]
            nls.score = scores[i].to_f
            result.signals.push nls
          end
        elsif html.match(monopartite_no_hits) or html.match(monopartite_no_hits2)
          # do nothing, except for not raising a parsing exception
        else
          raise Exception, "Could not parse HTML output returned from cNLS prediction server. In particular, looking for monopartite signals, but the whole document is likely problematic.\n#{html}"
        end
        
        
        # parse out the bipartite signals
        if matches = html.match(bipartite_regex)
          positions = matches[1].split(split_regex)
          seqs = matches[2].split(split_regex)
          scores = matches[3].split(split_regex)
          
          positions.each_with_index do |pos, i|
            nls = Result::BipartiteNLS.new
            nls.position = pos.to_i
            nls.sequence = seqs[i]
            nls.score = scores[i].to_f
            result.signals.push nls
          end
        elsif html.match(bipartite_no_hits) or html.match(bipartite_no_hits2)
          # do nothing, except for not raising a parsing exception
        else
          raise Exception, "Could not parse HTML output returned from cNLS prediction server. In particular, looking for bipartite signals, monopartite signals seemed to be parsed OK.\n#{html}"
        end
        
        return result
      end
    end
  end
end  



if __FILE__ == $0
  require 'optparse'
  
  # When entering sequences less than this number of amino acids as a query
  # it fails (if less than 10 it tells you, if less than 19 then it silently fails)
  QUERY_LENGTH_MINIMUM = 19
  ACCEPTABLE_AMINO_ACID_CHARACTERS = Bio::AminoAcid::Data::WEIGHT.keys.push('*')
  
  options = {
  :verbose => true,
  :cache_html => false
  }
  o = OptionParser.new do |opts|
    opts.banner = [
    'Usage: bio-cnls_formatter.rb [-qh] [fasta_filename]',
    '\tfasta file can also be piped in on STDIN.'
    ]
    opts.on('-q','--quiet','Opposite of verbose. Default is not quiet (verbose is on)') do
      options[:verbose] = false
    end
    opts.on('-h','--html','Cache HTML results in the current directory instead of parsing them') do
      options[:cache_html] = true
    end
  end
  o.parse!
  
  Bio::FlatFile.foreach(ARGF) do |entry|
    # Sequences are automatically disqualified if they contain characters that are neither amino acids or stop codons
    fails = entry.seq.gsub(/[#{ACCEPTABLE_AMINO_ACID_CHARACTERS.join('')}]/,'')
    if fails.length > 0
      if options[:verbose]
        $stderr.puts "Found unacceptable characters in #{entry.definition}: #{fails}"
      end
      next
      
      # Sequence length must be greater than the minimum, excluding
      # stop codons
    elsif entry.seq.gsub(/\*/,'').length < QUERY_LENGTH_MINIMUM
      if options[:verbose]
        $stderr.puts "Query sequence too short (less than #{QUERY_LENGTH_MINIMUM} residues excluding stop codons): #{entry.definition}"
      end
    else
      # This sequence passes, run the prediction on it
      if options[:cache_html]
        res = Bio::CNLS::Screenscraper.get_raw_html_result(entry.seq)
        File.open("#{entry.definition}.html",'w') do |f|
          f.puts res
        end
        $stderr.print '.' if options[:verbose]
      else
        res = Bio::CNLS::Screenscraper.submit(entry.seq)
        puts [
        entry.definition,
        res.predicted?
        ].join("\t")
      end
    end
  end
end