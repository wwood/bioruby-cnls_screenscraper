require 'helper'
require 'bio-cnls_screenscraper'

class TestBioCnlsScreenscraper < Test::Unit::TestCase
  @@data_dir = File.join(File.dirname(__FILE__),['data'])
 
  should "correctly parse hit results with no hits" do
    html = File.open(File.join(@@data_dir,'nohits.html')).read
    result = Bio::CNLS::Screenscraper.parse_html_result(html)
    assert_equal [], result.signals
  end
  
  should "correctly parse bipartite-only signals page" do
    html = File.open(File.join(@@data_dir,'bipartiteHitOnly.html')).read
    result = Bio::CNLS::Screenscraper.parse_html_result(html)
    assert_equal 2, result.signals.length
    assert_equal 'KKKRRRAAAAAAAAAAAAAAAAAARKKKRRR', result.signals.sort[1].sequence
    assert_equal 5.0, result.signals.sort[1].score
    assert_equal 1, result.signals.sort[1].position
  end
  
  should "correctly parse results with monopartite signals only" do
    html = File.open(File.join(@@data_dir,'monopartiteHitOnly.html')).read
    result = Bio::CNLS::Screenscraper.parse_html_result(html)
    assert_equal 1, result.signals.length
    assert_equal 'KKKKRRRAA', result.signals.sort[0].sequence
    assert_equal 10.0, result.signals.sort[0].score
    assert_equal 1, result.signals.sort[0].position
  end
  
  should "apply the correct monopartite cutoff" do
    nls = Bio::CNLS::Result::MonopartiteNLS.new
    nls.score = 8.0
    result = Bio::CNLS::Result.new
    result.signals.push nls
    assert_equal true, result.monopartite_predicted?(7.0)
    assert_equal true, result.monopartite_predicted?(8.0)
    assert_equal false, result.monopartite_predicted?(9.0)
  end
end
