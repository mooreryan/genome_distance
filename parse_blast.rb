class BlastRecord
  attr_accessor(:qseqid, :sseqid, :pident, :length, :mismatch,
                :gapopen, :qstart, :qend, :sstart, :send, :evalue,
                :bitscore)

  def initialize info
    @qseqid   = info[0]
    @sseqid   = info[1]
    @pident   = info[2].to_f
    @length   = info[3].to_i
    @mismatch = info[4].to_i
    @gapopen  = info[5].to_i
    @qstart   = info[6].to_i
    @qend     = info[7].to_i
    @sstart   = info[8].to_i
    @send     = info[9].to_i
    @evalue   = info[10].to_f
    @bitscore = info[11].to_f
  end
end

class Btab < File
  def each_record
    self.each_line do |line|
      record = line.chomp.split "\t"
      yield BlastRecord.new record
    end
  end
end

len_f = ARGV[0]
btab_f = ARGV[1]
lens = {}
num_hits = {}
total_bases = {}
shared_len = {}
base_cov = {}

File.open(len_f).each_line do |line|
  unless line.start_with? "name"
    contig, len = line.chomp.split " "

    if lens.has_key? contig
      abort "#{contig} repeated in #{len_f}"
    else
      lens[contig] = len.to_i
    end
  end
end

lens.keys.sort.combination(2).each do |a, b|
  key = [a, b].join("_-|-_")
  num_hits[key] = 0
  shared_len[key] = 0
  total_bases[key] = lens[a] + lens[b]
  base_cov[key] = { a => [],
                    b => [] }
end

Btab.open(btab_f).each_record do |record|
  unless record.qseqid == record.sseqid
    if record.pident > 75 && record.evalue < 1e-10

      # put them in order string alpha order
      if record.qseqid <= record.sseqid
        key = "#{record.qseqid}_-|-_#{record.sseqid}"
      else
        key = "#{record.sseqid}_-|-_#{record.qseqid}"
      end

      unless num_hits.has_key? key
        abort "problem with the key: #{key}"
      end

      num_hits[key] += 1

      qlen = (record.qstart - record.qend).abs + 1
      hlen = (record.sstart - record.send).abs + 1
      tlen = qlen + hlen

      unless shared_len.has_key? key
        abort "problem with the key: #{key}"
      end

      shared_len[key] += tlen

      unless base_cov.has_key? key
        abort "problem with the key: #{key}"
      end

      (record.qstart-1..record.qend-1).each do |n|
        base_cov[key][record.qseqid] << n
      end

      (record.sstart-1..record.send-1).each do |n|
        base_cov[key][record.sseqid] << n
      end
    end
  end
end

# num_hits.sort_by { |k, v| k }.each do |key, hits|
#   a, b = key.split("_-|-_")

#   puts [a, b, hits].join "\t"
# end
# puts
# shared_len.sort_by { |k, v| k }.each do |key, len|
#   a, b = key.split("_-|-_")

#   puts [a, b, len / total_bases[key].to_f].join "\t"
# end
# puts


# the distance is based on 1 - ((total.hit.len + total.query.len) / (query.genome.len + hit.genome.len))
base_cov.sort_by { |k, v| k }.each do |key, val|
  a, b = key.split("_-|-_")

  a_cov, b_cov = val[a].uniq, val[b].uniq

  shared_cov = a_cov.count + b_cov.count

  puts [a, b, 1 - (shared_cov / shared_len[key].to_f)].join "\t"
end
